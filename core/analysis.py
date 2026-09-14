"""
Lightweight, dependency-free analysis functions used by the live
simulation worker. Deliberately avoids MDAnalysis to keep the per-frame
overhead low enough for real-time reporting.
"""
import numpy as np


def kabsch_align(mobile: np.ndarray, ref: np.ndarray) -> np.ndarray:
    """Superpose `mobile` onto `ref` (both Nx3, already centered or not).
    Returns the aligned mobile coordinates.

    BUG FIX: the final translation previously used `ref_c.mean(axis=0)` —
    but ref_c is ref AFTER subtracting its own mean, so ref_c.mean(axis=0)
    is ~0 by construction. That silently discarded ref's actual position
    in space and re-centered the aligned structure near the origin instead
    of onto ref. Since `rmsd()` below then computes `coords - ref` using
    the ORIGINAL (uncentered) ref, every RMSD was dominated by a large,
    essentially constant offset (ref's distance from the origin) rather
    than the real conformational difference — which is exactly why RMSD
    plots looked flat while Rg (which never uses this alignment) showed
    real variation. Fixed by using ref.mean(axis=0), ref's actual centroid.
    """
    mobile_c = mobile - mobile.mean(axis=0)
    ref_c = ref - ref.mean(axis=0)

    H = mobile_c.T @ ref_c
    U, _S, Vt = np.linalg.svd(H)
    d = np.sign(np.linalg.det(Vt.T @ U.T))
    D = np.diag([1, 1, d])
    R = Vt.T @ D @ U.T

    aligned = (R @ mobile_c.T).T + ref.mean(axis=0)
    return aligned


def kabsch_transform(mobile_sub: np.ndarray, ref_sub: np.ndarray):
    """Return (R, mobile_centroid, ref_centroid) from a SUBSET of atoms.

    Separated out from kabsch_align so a rotation derived from one group of
    atoms can be applied to a different group — which is the whole basis of
    a "fit on receptor, measure binder" RMSD.
    """
    mc = mobile_sub.mean(axis=0)
    rc = ref_sub.mean(axis=0)
    H = (mobile_sub - mc).T @ (ref_sub - rc)
    U, _S, Vt = np.linalg.svd(H)
    d = np.sign(np.linalg.det(Vt.T @ U.T))
    R = Vt.T @ np.diag([1, 1, d]) @ U.T
    return R, mc, rc


def rmsd_fit_on(coords: np.ndarray, ref: np.ndarray,
                fit_idx: np.ndarray, measure_idx: np.ndarray) -> float:
    """Superpose on `fit_idx` atoms, then measure RMSD over `measure_idx`.

    With fit_idx == measure_idx this is an ordinary aligned RMSD. The useful
    case is fit_idx = receptor, measure_idx = binder: the receptor is held
    still and the binder's movement RELATIVE TO THE BINDING SITE is what gets
    measured. That is the quantity that tells you whether a pose is holding,
    and it is invisible to both a whole-complex RMSD (which averages the
    binder's drift away across all the receptor atoms that did not move) and
    to a centroid interface distance (which cannot see rotation or sliding
    that preserves the centroid separation).
    """
    R, mc, rc = kabsch_transform(coords[fit_idx], ref[fit_idx])
    moved = (R @ (coords[measure_idx] - mc).T).T + rc
    diff = moved - ref[measure_idx]
    return float(np.sqrt(np.mean(np.sum(diff * diff, axis=1))))


def rmsd(coords: np.ndarray, ref: np.ndarray, align: bool = True) -> float:
    """RMSD in the same length units as input (nm if fed OpenMM positions)."""
    if align:
        coords = kabsch_align(coords, ref)
    diff = coords - ref
    return float(np.sqrt(np.mean(np.sum(diff * diff, axis=1))))


def radius_of_gyration(coords: np.ndarray, masses: np.ndarray) -> float:
    total_mass = masses.sum()
    com = (coords * masses[:, None]).sum(axis=0) / total_mass
    diff = coords - com
    rg_sq = (masses * np.sum(diff * diff, axis=1)).sum() / total_mass
    return float(np.sqrt(rg_sq))


def interface_distance(coords: np.ndarray, receptor_idx: np.ndarray,
                        binder_idx: np.ndarray) -> float:
    """Centroid-centroid distance between two atom-index groups.
    Cheap proxy for a dissociation coordinate while the sim is running."""
    rc = coords[receptor_idx].mean(axis=0)
    bc = coords[binder_idx].mean(axis=0)
    return float(np.linalg.norm(rc - bc))


class RunningRMSF:
    """Per-residue RMSF accumulated incrementally, one frame at a time.

    RMSF is normally computed after the fact from a whole trajectory, which
    is not an option here — the point of this app is that the numbers appear
    while the run is happening. Storing every frame's coordinates would also
    grow without bound on a long run.

    Instead this keeps only running sums, so memory is O(n_residues) no
    matter how many frames arrive:

        RMSF_i = sqrt( <|r_i|^2> - |<r_i>|^2 )

    computed per residue from the sum and sum-of-squares of its CA position.
    That identity is exact, not an approximation, and is what lets the plot
    update live.

    Every frame is Kabsch-aligned to the reference first, otherwise global
    tumbling and drift of the whole molecule would swamp the local
    fluctuation that RMSF is meant to measure.
    """

    def __init__(self, ref_ca: np.ndarray):
        self.ref = np.asarray(ref_ca, dtype=float)
        n = len(self.ref)
        self._sum = np.zeros((n, 3))
        self._sum_sq = np.zeros(n)
        self._n = 0

    def add_frame(self, ca_coords: np.ndarray) -> None:
        aligned = kabsch_align(np.asarray(ca_coords, dtype=float), self.ref)
        self._sum += aligned
        self._sum_sq += np.sum(aligned * aligned, axis=1)
        self._n += 1

    @property
    def n_frames(self) -> int:
        return self._n

    def values(self) -> np.ndarray:
        """Per-residue RMSF, same units as input (nm for OpenMM positions).

        Returns zeros until at least two frames have been added — a single
        frame has no fluctuation to measure, and reporting a number there
        would be meaningless rather than merely imprecise.
        """
        if self._n < 2:
            return np.zeros(len(self.ref))
        mean = self._sum / self._n
        mean_sq = self._sum_sq / self._n
        var = mean_sq - np.sum(mean * mean, axis=1)
        # Tiny negatives are possible from floating-point cancellation when a
        # residue barely moves; clamp rather than emit NaN from sqrt.
        return np.sqrt(np.clip(var, 0.0, None))
