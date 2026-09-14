import csv
import json
import os

import numpy as np
import pyqtgraph as pg
import pyqtgraph.exporters
from PySide6.QtCore import Qt
from PySide6.QtGui import QColor, QIcon
from PySide6.QtWebEngineCore import QWebEnginePage
from PySide6.QtWebEngineWidgets import QWebEngineView
from PySide6.QtWidgets import (
    QCheckBox,
    QColorDialog,
    QComboBox,
    QDoubleSpinBox,
    QFileDialog,
    QFormLayout,
    QGroupBox,
    QHBoxLayout,
    QLabel,
    QLineEdit,
    QMainWindow,
    QPushButton,
    QSpinBox,
    QSplitter,
    QTextEdit,
    QVBoxLayout,
    QWidget,
)

from core.forcefields import (
    DEFAULT_FORCE_FIELD,
    DEFAULT_IMPLICIT,
    FORCE_FIELDS,
    IMPLICIT_MODELS,
    force_field_names,
    water_names,
)
from core.simulation_worker import SimulationWorker

HERE = os.path.dirname(os.path.abspath(__file__))


class DiagnosticWebEnginePage(QWebEnginePage):
    """Overrides javaScriptConsoleMessage to catch EVERYTHING the JS side
    logs or throws — including uncaught errors from async/internal 3Dmol
    code paths that would never reach our own try/catch + document.title
    convention. This was added after a real bug slipped through
    completely silently: a structure-load call appeared to just vanish
    with no success AND no error message anywhere in the log, and the
    existing error-reporting mechanism (document.title polling) had no
    visibility into why. This is a strictly more complete safety net.
    """
    console_message = None  # set by MainWindow to a callback(level, msg, line, source)

    def javaScriptConsoleMessage(self, level, message, lineNumber, sourceID):
        if self.console_message is not None:
            self.console_message(level, message, lineNumber, sourceID)

# ---------------------------------------------------------------------------
# Dark theme — applied as a stylesheet on the whole window rather than a
# QPalette, since QPalette alone doesn't reliably restyle QGroupBox borders,
# QComboBox popups, etc. across platforms.
# ---------------------------------------------------------------------------
DARK_STYLESHEET = """
QMainWindow, QWidget {
    background-color: #16171c;
    color: #e6e6e6;
    font-size: 13px;
}
QGroupBox {
    border: 1px solid #2c2e36;
    border-radius: 6px;
    margin-top: 10px;
    padding-top: 10px;
    font-weight: 600;
    color: #cfd2da;
}
QGroupBox::title {
    subcontrol-origin: margin;
    left: 10px;
    padding: 0 4px;
}
QLineEdit, QSpinBox, QDoubleSpinBox, QComboBox, QTextEdit {
    background-color: #1e2028;
    border: 1px solid #33353f;
    border-radius: 4px;
    padding: 4px 6px;
    color: #e6e6e6;
    selection-background-color: #3a6ea5;
}
QComboBox QAbstractItemView {
    background-color: #1e2028;
    color: #e6e6e6;
    selection-background-color: #3a6ea5;
}
QPushButton {
    background-color: #262933;
    border: 1px solid #383b46;
    border-radius: 5px;
    padding: 6px 10px;
    color: #e6e6e6;
}
QPushButton:hover { background-color: #2f3340; }
QPushButton:pressed { background-color: #1c1e26; }
QPushButton:disabled { color: #6b6e78; background-color: #1c1e24; }
QLabel { color: #cfd2da; }
QCheckBox { color: #e6e6e6; spacing: 6px; }
QCheckBox::indicator {
    width: 14px; height: 14px;
    border: 1px solid #4a4d58;
    border-radius: 3px;
    background: #1e2028;
}
QCheckBox::indicator:checked {
    background: #3a6ea5;
    border: 1px solid #4a86c4;
}
QSplitter::handle { background-color: #22242c; }
QScrollBar:vertical { background: #1a1c22; width: 10px; }
QScrollBar::handle:vertical { background: #383b46; border-radius: 5px; }
"""

DEFAULT_RECEPTOR_COLOR = "#7fc8f8"   # light blue
DEFAULT_BINDER_COLOR = "#f5a742"     # orange


class ColorPickerButton(QPushButton):
    """Small button showing a color swatch; click opens QColorDialog."""

    def __init__(self, initial_hex, parent=None):
        super().__init__(parent)
        self.color_hex = initial_hex
        self.setFixedWidth(48)
        self._refresh()
        self.clicked.connect(self._pick)

    def _refresh(self):
        self.setStyleSheet(
            f"background-color: {self.color_hex}; border: 1px solid #444; border-radius: 4px;"
        )
        self.setText("")

    def _pick(self):
        col = QColorDialog.getColor(QColor(self.color_hex), self, "Pick a color")
        if col.isValid():
            self.color_hex = col.name()
            self._refresh()


class MainWindow(QMainWindow):
    def __init__(self):
        super().__init__()
        self.setWindowTitle("Jiggle — live molecular dynamics viewer")
        self.resize(1450, 880)
        self.setStyleSheet(DARK_STYLESHEET)

        icon_path = os.path.join(HERE, "assets", "icon.png")
        if os.path.exists(icon_path):
            self.setWindowIcon(QIcon(icon_path))

        self.worker = None
        self._steps, self._rmsd, self._rg, self._iface = [], [], [], []
        self._pe, self._ke, self._times = [], [], []
        self._rmsd_rec, self._rmsd_bnd_self, self._rmsd_bnd_on_rec = [], [], []
        self._rmsf, self._rmsf_resids, self._rmsf_labels = [], [], []
        self._pending_chains = ([], [])
        # True from the moment a structure-load call is dispatched until
        # its completion callback fires. Coordinate updates are held back
        # while this is set — see _load_structure_js / _on_frame. This
        # exists because of a real, confirmed-from-an-actual-run failure
        # mode: a run's log showed the 642-atom structure-load call being
        # dispatched AND its completion confirmed, yet a coordinate update
        # still hit a stale 327-atom model — meaning at least one
        # updateCoords() call reached the renderer before the structure
        # load had taken effect there, despite Python having already
        # "sent" the load call first. Rather than rely on assumptions
        # about QtWebEngine's internal script-ordering guarantees (which
        # this evidence contradicts in practice), Python now enforces the
        # ordering itself: no coordinates are sent until the load is
        # confirmed complete.
        self._structure_load_in_flight = False

        root = QWidget()
        layout = QHBoxLayout(root)
        self.setCentralWidget(root)

        splitter = QSplitter(Qt.Horizontal)
        layout.addWidget(splitter)

        splitter.addWidget(self._build_input_panel())
        splitter.addWidget(self._build_viewer_panel())
        splitter.addWidget(self._build_plots_panel())
        splitter.setSizes([340, 700, 410])

        # Populate the force-field-dependent widgets and the steps<->ns pair
        # once everything exists.
        self._on_forcefield_changed(self.forcefield_combo.currentText())
        self._on_total_steps_changed(self.total_steps_spin.value())

    # ---------------- Left: inputs ----------------
    def _build_input_panel(self):
        box = QGroupBox("Setup")
        form = QFormLayout(box)

        self.pdb_path_edit = QLineEdit()
        self.pdb_path_edit.setPlaceholderText("Path to a .pdb file — single protein or complex")
        self.pdb_path_edit.textChanged.connect(self._on_pdb_path_changed)
        browse_btn = QPushButton("Browse...")
        browse_btn.clicked.connect(self._browse_pdb)
        pdb_row = QHBoxLayout()
        pdb_row.addWidget(self.pdb_path_edit)
        pdb_row.addWidget(browse_btn)
        pdb_row_w = QWidget()
        pdb_row_w.setLayout(pdb_row)
        form.addRow("PDB file:", pdb_row_w)

        self.detected_chains_label = QLabel("Detected chains: —")
        self.detected_chains_label.setStyleSheet("color: #8a8f9c; font-style: italic;")
        form.addRow(self.detected_chains_label)

        self.receptor_chains_edit = QLineEdit()
        self.receptor_chains_edit.setPlaceholderText("auto-filled on file load, e.g. A")
        form.addRow("Receptor chain(s):", self.receptor_chains_edit)

        self.binder_chains_edit = QLineEdit()
        self.binder_chains_edit.setPlaceholderText("leave empty for a single protein / monomer")
        form.addRow("Binder chain(s):", self.binder_chains_edit)

        # --- color pickers ---
        self.receptor_color_btn = ColorPickerButton(DEFAULT_RECEPTOR_COLOR)
        self.binder_color_btn = ColorPickerButton(DEFAULT_BINDER_COLOR)
        color_row = QHBoxLayout()
        color_row.addWidget(QLabel("Receptor"))
        color_row.addWidget(self.receptor_color_btn)
        color_row.addSpacing(12)
        color_row.addWidget(QLabel("Binder"))
        color_row.addWidget(self.binder_color_btn)
        color_row.addStretch()
        color_row_w = QWidget()
        color_row_w.setLayout(color_row)
        form.addRow("Colors:", color_row_w)
        # Re-render immediately on color change if a structure is already loaded.
        self.receptor_color_btn.clicked.connect(self._recolor_if_loaded)
        self.binder_color_btn.clicked.connect(self._recolor_if_loaded)

        self.forcefield_combo = QComboBox()
        self.forcefield_combo.addItems(force_field_names())
        self.forcefield_combo.setCurrentText(DEFAULT_FORCE_FIELD)
        self.forcefield_combo.currentTextChanged.connect(self._on_forcefield_changed)
        form.addRow("Force field:", self.forcefield_combo)

        self.forcefield_note = QLabel("")
        self.forcefield_note.setWordWrap(True)
        self.forcefield_note.setStyleSheet("color: #8a8f9c; font-style: italic;")
        form.addRow(self.forcefield_note)

        self.solvent_combo = QComboBox()
        self.solvent_combo.addItems(["implicit", "explicit"])
        self.solvent_combo.currentTextChanged.connect(self._on_solvent_changed)
        form.addRow("Solvent:", self.solvent_combo)

        # One combo serving two roles: the GB model in implicit mode, the water
        # model in explicit mode. Its contents are rebuilt from the catalogue
        # whenever the force field or solvent changes, so the UI can never
        # offer a pairing the worker is unable to build.
        self.solvent_model_label = QLabel("Implicit model:")
        self.solvent_model_combo = QComboBox()
        form.addRow(self.solvent_model_label, self.solvent_model_combo)

        self.strip_heterogens_check = QCheckBox(
            "Strip crystallization additives (glycerol, PEG, sulfate, "
            "crystallographic water, etc.)"
        )
        self.strip_heterogens_check.setChecked(True)
        form.addRow(self.strip_heterogens_check)

        self.temp_spin = QDoubleSpinBox()
        self.temp_spin.setRange(200, 400)
        self.temp_spin.setValue(300)
        form.addRow("Temperature (K):", self.temp_spin)

        self.minimize_spin = QSpinBox()
        self.minimize_spin.setRange(0, 1_000_000)
        self.minimize_spin.setValue(0)
        self.minimize_spin.setSpecialValueText("0 — run to convergence")
        self.minimize_spin.setToolTip(
            "Maximum L-BFGS iterations for the pre-run energy minimization.\n"
            "0 lets OpenMM run until its own convergence tolerance is met, "
            "which is the right choice for generated or docked structures "
            "with local clashes. Cap it only if minimization is taking longer "
            "than you are willing to wait."
        )
        form.addRow("Minimization steps:", self.minimize_spin)

        self.timestep_spin = QDoubleSpinBox()
        self.timestep_spin.setRange(0.5, 4.0)
        self.timestep_spin.setValue(2.0)
        self.timestep_spin.setSingleStep(0.5)
        self.timestep_spin.valueChanged.connect(self._on_timestep_changed)
        form.addRow("Timestep (fs):", self.timestep_spin)

        # Length can be entered either way round. Steps stays the authoritative
        # value that goes to the worker; nanoseconds is a view onto it. Both
        # directions are guarded by _syncing so the two signals can't ping-pong.
        self._syncing = False

        self.duration_ns_spin = QDoubleSpinBox()
        self.duration_ns_spin.setRange(0.0001, 10000.0)
        self.duration_ns_spin.setDecimals(4)
        self.duration_ns_spin.setSingleStep(0.1)
        self.duration_ns_spin.setSuffix(" ns")
        self.duration_ns_spin.valueChanged.connect(self._on_duration_changed)
        form.addRow("Simulation length:", self.duration_ns_spin)

        self.total_steps_spin = QSpinBox()
        self.total_steps_spin.setRange(1000, 100_000_000)
        self.total_steps_spin.setValue(500_000)
        self.total_steps_spin.setGroupSeparatorShown(True)
        self.total_steps_spin.valueChanged.connect(self._on_total_steps_changed)
        form.addRow("Total steps:", self.total_steps_spin)

        self.report_interval_spin = QSpinBox()
        self.report_interval_spin.setRange(1, 10000)
        self.report_interval_spin.setValue(50)
        self.report_interval_spin.valueChanged.connect(self._refresh_length_summary)
        form.addRow("Steps / frame:", self.report_interval_spin)

        self.length_summary = QLabel("")
        self.length_summary.setWordWrap(True)
        self.length_summary.setStyleSheet("color: #8a8f9c; font-style: italic;")
        form.addRow(self.length_summary)

        self.output_dir_edit = QLineEdit("./md_output")
        form.addRow("Output dir:", self.output_dir_edit)

        self.platform_combo = QComboBox()
        self.platform_combo.addItems(["auto", "CUDA", "OpenCL", "CPU", "Reference"])
        form.addRow("Platform:", self.platform_combo)

        self.preview_btn = QPushButton("Preview Structure")
        self.preview_btn.clicked.connect(self._preview_structure)
        form.addRow(self.preview_btn)

        self.start_btn = QPushButton("Start Simulation")
        self.start_btn.clicked.connect(self._start_simulation)
        form.addRow(self.start_btn)

        self.pause_btn = QPushButton("Pause")
        self.pause_btn.setEnabled(False)
        self.pause_btn.clicked.connect(self._toggle_pause)
        self.stop_btn = QPushButton("Stop")
        self.stop_btn.setEnabled(False)
        self.stop_btn.clicked.connect(self._stop_simulation)
        run_row = QHBoxLayout()
        run_row.addWidget(self.pause_btn)
        run_row.addWidget(self.stop_btn)
        run_row_w = QWidget()
        run_row_w.setLayout(run_row)
        form.addRow(run_row_w)

        self.save_frame_btn = QPushButton("Save current frame as PDB")
        self.save_frame_btn.setEnabled(False)
        self.save_frame_btn.clicked.connect(self._save_current_frame)
        form.addRow(self.save_frame_btn)

        self.log_box = QTextEdit()
        self.log_box.setReadOnly(True)
        form.addRow(QLabel("Log:"))
        form.addRow(self.log_box)

        return box

    # ---------------- Force field / solvent model ----------------
    def _on_forcefield_changed(self, name):
        self.forcefield_note.setText(FORCE_FIELDS[name]["note"])
        self._refresh_solvent_models()

    def _on_solvent_changed(self, _solvent):
        self._refresh_solvent_models()

    def _refresh_solvent_models(self):
        """Rebuild the second combo from the catalogue for the current
        (force field, solvent) pair. Preserves the user's choice by name when
        the new list still contains it."""
        previous = self.solvent_model_combo.currentText()
        explicit = self.solvent_combo.currentText() == "explicit"
        if explicit:
            self.solvent_model_label.setText("Water model:")
            options = water_names(self.forcefield_combo.currentText())
        else:
            self.solvent_model_label.setText("Implicit model:")
            options = list(IMPLICIT_MODELS.keys())

        self.solvent_model_combo.blockSignals(True)
        self.solvent_model_combo.clear()
        self.solvent_model_combo.addItems(options)
        if previous in options:
            self.solvent_model_combo.setCurrentText(previous)
        self.solvent_model_combo.blockSignals(False)

    # ---------------- Length: steps <-> nanoseconds ----------------
    # total_steps is what the worker actually receives; nanoseconds is a
    # derived view. Editing either updates the other. _syncing prevents the
    # two valueChanged signals from bouncing off each other, and blockSignals
    # is used as well because QSpinBox rounds its input, which would otherwise
    # feed a slightly different value straight back.
    def _steps_to_ns(self, steps):
        return steps * self.timestep_spin.value() / 1_000_000.0

    def _ns_to_steps(self, ns):
        return round(ns * 1_000_000.0 / self.timestep_spin.value())

    def _on_total_steps_changed(self, steps):
        if self._syncing:
            return
        self._syncing = True
        self.duration_ns_spin.blockSignals(True)
        self.duration_ns_spin.setValue(self._steps_to_ns(steps))
        self.duration_ns_spin.blockSignals(False)
        self._syncing = False
        self._refresh_length_summary()

    def _on_duration_changed(self, ns):
        if self._syncing:
            return
        self._syncing = True
        steps = max(self.total_steps_spin.minimum(),
                    min(self.total_steps_spin.maximum(), self._ns_to_steps(ns)))
        self.total_steps_spin.blockSignals(True)
        self.total_steps_spin.setValue(steps)
        self.total_steps_spin.blockSignals(False)
        self._syncing = False
        self._refresh_length_summary()

    def _on_timestep_changed(self, _fs):
        # Steps stay fixed; the wall of simulated time they represent changes.
        self._on_total_steps_changed(self.total_steps_spin.value())

    def _refresh_length_summary(self, *_):
        steps = self.total_steps_spin.value()
        dt_fs = self.timestep_spin.value()
        interval = self.report_interval_spin.value()
        n_frames = steps // interval
        ps_per_frame = interval * dt_fs / 1000.0
        self.length_summary.setText(
            f"{steps:,} steps x {dt_fs:g} fs = {self._steps_to_ns(steps):.4g} ns  ·  "
            f"{n_frames:,} frames, {ps_per_frame:.3g} ps apart.  "
            f"Heating adds 1500 steps ({1500 * dt_fs / 1000.0:.3g} ps) on top, "
            f"not counted here."
        )

    def _browse_pdb(self):
        path, _ = QFileDialog.getOpenFileName(self, "Select a PDB file", "", "PDB files (*.pdb)")
        if path:
            self.pdb_path_edit.setText(path)

    def _set_preview_enabled(self, enabled):
        self.preview_btn.setEnabled(enabled)

    @staticmethod
    def _detect_chains(pdb_path):
        """Fast, PDBFixer-free scan of ATOM/HETATM records for distinct
        chain IDs, in first-appearance order. Uses the fixed-column PDB
        format (chain ID is column 22, 0-indexed 21) rather than OpenMM,
        so this is instant even before deciding whether to fix/minimize
        anything — good enough just to tell the user what's in the file.
        """
        chains = []
        seen = set()
        try:
            with open(pdb_path, errors="ignore") as f:
                for line in f:
                    if line.startswith(("ATOM", "HETATM")) and len(line) > 21:
                        ch = line[21:22].strip()
                        if ch and ch not in seen:
                            seen.add(ch)
                            chains.append(ch)
        except OSError:
            pass
        return chains

    def _on_pdb_path_changed(self, path):
        self._set_preview_enabled(True)
        path = path.strip()
        if not path or not os.path.exists(path):
            self.detected_chains_label.setText("Detected chains: —")
            return
        if path == getattr(self, "_last_detected_path", None):
            return  # already handled this exact path — don't fight manual edits
        self._last_detected_path = path

        chains = self._detect_chains(path)
        if not chains:
            self.detected_chains_label.setText("Detected chains: none found (not a valid PDB?)")
            return

        self.detected_chains_label.setText(f"Detected chains: {', '.join(chains)}")

        # Auto-populate receptor/binder fields as a starting guess. This only
        # runs once per newly-selected file (guarded above), so it won't
        # clobber chains the user then edits by hand.
        if len(chains) == 1:
            self.receptor_chains_edit.setText(chains[0])
            self.binder_chains_edit.setText("")
            self.log_box.append(
                f"Detected 1 chain ({chains[0]}) — treating as a single "
                f"protein/monomer. No interface-distance metric will be "
                f"computed; leave 'Binder chain(s)' empty."
            )
        else:
            self.receptor_chains_edit.setText(chains[0])
            self.binder_chains_edit.setText(",".join(chains[1:]))
            self.log_box.append(
                f"Detected {len(chains)} chains ({', '.join(chains)}) — "
                f"guessed receptor={chains[0]}, binder={','.join(chains[1:])}. "
                f"Edit the chain fields above if that split is wrong."
            )

    # ---------------- Preview (pre-simulation, raw structure) ----------------
    def _preview_structure(self):
        """Load the raw, as-uploaded PDB into the 3D view so the user can
        rotate/zoom and sanity-check chain IDs before committing to a run.
        Deliberately skips PDBFixer/OpenMM — this is just a look, not a sim."""
        pdb_path = self.pdb_path_edit.text().strip()
        if not pdb_path or not os.path.exists(pdb_path):
            self.log_box.append("Pick a valid complex PDB first.")
            return
        with open(pdb_path) as f:
            pdb_text = f.read()

        receptor_chains = [c.strip() for c in self.receptor_chains_edit.text().split(",") if c.strip()]
        binder_chains = [c.strip() for c in self.binder_chains_edit.text().split(",") if c.strip()]
        self._pending_chains = (receptor_chains, binder_chains)
        self._load_structure_js(pdb_text)

    def _recolor_if_loaded(self):
        # Re-run the same load call with new colors, if something is on screen.
        # No-op silently if nothing has been loaded/started yet.
        if getattr(self, "_last_pdb_text", None):
            self._load_structure_js(self._last_pdb_text)

    def _load_structure_js(self, pdb_text):
        self._last_pdb_text = pdb_text
        receptor_chains, binder_chains = self._pending_chains
        # Count HETATM as well: this number is compared against the viewer's
        # parsed atom count in the log, and any difference between the two is
        # the signature of the parser silently discarding atoms.
        n_atom_lines = sum(
            1 for ln in pdb_text.split("\n") if ln.startswith(("ATOM  ", "HETATM"))
        )
        js = (
            f"loadInitialStructure({json.dumps(pdb_text)}, "
            f"{json.dumps(receptor_chains)}, {json.dumps(binder_chains)}, "
            f"{json.dumps(self.receptor_color_btn.color_hex)}, "
            f"{json.dumps(self.binder_color_btn.color_hex)});"
        )
        if self._viewer_ready:
            self.log_box.append(
                f"[dispatch] Sending structure to 3D viewer now "
                f"({n_atom_lines} ATOM lines in the PDB text being sent)."
            )
            # Block coordinate updates until this completes — see the
            # _structure_load_in_flight docstring in __init__ for why this
            # is necessary rather than assumed-safe.
            self._structure_load_in_flight = True

            def _on_load_confirmed(_result):
                self._structure_load_in_flight = False
                self.log_box.append(
                    "[dispatch] 3D viewer confirmed it executed the structure-load call."
                )
                # If real frames arrived while the load was in flight, catch
                # up immediately with the freshest one now that it's safe.
                pending = getattr(self, "_last_frame_payload", None)
                if pending is not None:
                    self._render_frame_now(pending)

            self.web_view.page().runJavaScript(js, _on_load_confirmed)
        else:
            self.log_box.append(
                f"[dispatch] 3D viewer page not ready yet — queuing structure "
                f"load ({n_atom_lines} ATOM lines) for when it finishes loading."
            )
            self._pending_structure_js = js

    def _render_frame_now(self, payload):
        """Send one frame's coordinates to the viewer immediately, with no
        throttling/gating check — callers are responsible for having
        already confirmed it's safe to do so."""
        coords_nm = payload["positions_nm"]
        coords_angstrom_flat = (coords_nm * 10.0).flatten().tolist()
        self.web_view.page().runJavaScript(f"updateCoords({json.dumps(coords_angstrom_flat)});")

    # ---------------- Center: 3D viewer ----------------
    def _build_viewer_panel(self):
        box = QGroupBox("Live 3D View")
        v = QVBoxLayout(box)
        self.web_view = QWebEngineView()

        # Custom page catches EVERY console message/uncaught error from the
        # JS side, not just the ones that flow through our own
        # document.title convention — see DiagnosticWebEnginePage docstring.
        diag_page = DiagnosticWebEnginePage(self.web_view)
        diag_page.console_message = self._on_js_console_message
        self.web_view.setPage(diag_page)

        self._viewer_ready = False
        self._pending_structure_js = None

        template_path = os.path.join(HERE, "viewer_template.html")
        vendor_js_path = os.path.join(HERE, "vendor", "3Dmol-min.js")

        with open(template_path) as f:
            html = f.read()
        with open(vendor_js_path) as f:
            vendor_js = f.read()

        # Inline the vendored 3Dmol.js so the view never depends on network
        # access or QWebEngineView's local-file loading rules.
        html = html.replace(
            "<!-- 3DMOL_JS_INJECT: main_window.py replaces this comment with the full\n"
            "     contents of vendor/3Dmol-min.js, inlined, so there is zero external\n"
            "     or local-file resource loading for QWebEngineView to block. -->",
            f"<script>\n{vendor_js}\n</script>",
        )

        # Any JS-side error sets document.title to "JSERROR::...". Catch
        # that here and surface it in the log panel instead of failing
        # silently behind a black viewport.
        self.web_view.titleChanged.connect(self._on_viewer_title_changed)
        self.web_view.loadFinished.connect(self._on_viewer_load_finished)

        self.web_view.setHtml(html)
        v.addWidget(self.web_view)
        return box

    def _on_viewer_load_finished(self, ok):
        self._viewer_ready = ok
        if not ok:
            self.log_box.append("[3D viewer error] page failed to load.")
        elif self._pending_structure_js is not None:
            self.web_view.page().runJavaScript(self._pending_structure_js)
            self._pending_structure_js = None

    def _on_viewer_title_changed(self, title):
        if title.startswith("JSERROR::"):
            self.log_box.append(f"[3D viewer error] {title[len('JSERROR::'):]}")
        elif title.startswith("READY::"):
            detail = title[len("READY::"):]
            self.log_box.append(f"[3D viewer] structure loaded OK. {detail}")

    def _on_js_console_message(self, level, message, line_number, source_id):
        # `level` is a QWebEnginePage.JavaScriptConsoleMessageLevel enum,
        # not a plain int — comparing it directly against 0/1/2 raises
        # TypeError (caught by testing before this ever reached a user).
        # Compare against the enum's .value instead.
        level_value = getattr(level, "value", level)
        if level_value >= 2:
            self.log_box.append(f"[3D viewer JS console] ERROR (line {line_number}): {message}")
        elif level_value == 1:
            self.log_box.append(f"[3D viewer JS console] warning (line {line_number}): {message}")
        # Info-level (0) messages are extremely noisy in practice (3Dmol
        # logs a lot at that level) — deliberately not surfaced.

    # ---------------- Right: live plots ----------------
    def _build_plots_panel(self):
        box = QGroupBox("Live Metrics")
        v = QVBoxLayout(box)

        pg.setConfigOptions(antialias=True)
        pg.setConfigOption("background", "#0e0f13")
        pg.setConfigOption("foreground", "#cfd2da")

        self.rmsd_plot = self._make_styled_plot("RMSD (nm)")
        self.rg_plot = self._make_styled_plot("Radius of Gyration (nm)")
        self.energy_plot = self._make_styled_plot("Potential energy (kJ/mol)")
        self.rmsf_plot = self._make_styled_plot("Per-residue RMSF (nm)")
        self.iface_plot = self._make_styled_plot("Interface distance (nm)")

        self.rmsd_plot.addLegend(offset=(-10, 10), labelTextColor="#cfd2da")
        self.rmsd_curve = self._add_glow_curve(self.rmsd_plot, "#39d6e0", name="complex")
        # Only shown once a binder chain is actually present. A single RMSD
        # trace for a two-body system averages the binder's drift across every
        # receptor atom that did not move, so a pose coming loose barely
        # registers — see rmsd_fit_on() in core/analysis.py.
        self.rmsd_rec_curve = self.rmsd_plot.plot(
            [], [], pen=pg.mkPen("#6fe08a", width=2, style=Qt.DashLine), name="receptor")
        self.rmsd_bnd_self_curve = self.rmsd_plot.plot(
            [], [], pen=pg.mkPen("#e0a44a", width=2, style=Qt.DashLine), name="binder (internal)")
        self.rmsd_bnd_on_rec_curve = self.rmsd_plot.plot(
            [], [], pen=pg.mkPen("#e0479f", width=2), name="binder (fit on receptor)")
        self.rg_curve = self._add_glow_curve(self.rg_plot, "#e0c93f")         # amber
        self.energy_curve = self._add_glow_curve(self.energy_plot, "#7ee081") # green
        self.rmsf_curve = self._add_glow_curve(self.rmsf_plot, "#b48ce0")     # violet
        self.iface_curve = self._add_glow_curve(self.iface_plot, "#e0479f")   # magenta

        # RMSF is per residue, not per step — it is the one plot here whose
        # x axis is not time.
        self.rmsf_plot.setLabel("bottom", "Residue")
        for p in (self.rmsd_plot, self.rg_plot, self.energy_plot, self.iface_plot):
            p.setLabel("bottom", "Step")

        self._time_plots = [
            ("rmsd", self.rmsd_plot), ("radius_of_gyration", self.rg_plot),
            ("potential_energy", self.energy_plot),
            ("interface_distance", self.iface_plot),
        ]
        self._all_plots = [*self._time_plots, ("rmsf", self.rmsf_plot)]

        for _, p in self._all_plots:
            v.addWidget(p)

        self.export_btn = QPushButton("Save plots + data")
        self.export_btn.setEnabled(False)
        self.export_btn.clicked.connect(lambda: self._export_results(auto=False))
        v.addWidget(self.export_btn)

        return box

    def _make_styled_plot(self, title):
        p = pg.PlotWidget(title=f"<span style='color:#e6e6e6;font-size:11pt'>{title}</span>")
        p.showGrid(x=True, y=True, alpha=0.15)
        p.getAxis("left").setPen(pg.mkPen("#4a4d58"))
        p.getAxis("bottom").setPen(pg.mkPen("#4a4d58"))
        p.getAxis("left").setTextPen(pg.mkPen("#9ea2ad"))
        p.getAxis("bottom").setTextPen(pg.mkPen("#9ea2ad"))
        p.setBackground("#101116")
        return p

    def _add_glow_curve(self, plot_widget, hex_color, name=None):
        """A soft 'glow' line: a wide, translucent pen underneath a thin,
        bright one on top, plus a subtle gradient fill under the curve."""
        color = pg.mkColor(hex_color)

        fill_color = pg.mkColor(hex_color)
        fill_color.setAlpha(35)
        fill_curve = plot_widget.plot([], [], pen=None, fillLevel=0, brush=fill_color)

        glow_color = pg.mkColor(hex_color)
        glow_color.setAlpha(70)
        glow_curve = plot_widget.plot([], [], pen=pg.mkPen(glow_color, width=7))

        # The legend entry has to be registered at creation. Calling
        # setData(name=...) afterwards does NOT add the curve to an existing
        # legend — it just sets an option nothing reads, which is why the
        # complex trace had no legend entry.
        main_curve = plot_widget.plot([], [], pen=pg.mkPen(color, width=2), name=name)
        # Keep explicit references (not relying on introspecting the plot's
        # item list) so _update_curve can update all three together.
        main_curve._glow_siblings = [fill_curve, glow_curve]
        main_curve._plot_widget = plot_widget

        # BUG FIX: fillLevel=0 makes the fill curve's shape extend down to
        # y=0, and pyqtgraph's default auto-range fits ALL plotted items'
        # bounding boxes — including that fill's extent to zero. For a
        # metric like Rg that sits in a tight band far from zero (e.g.
        # 0.85-0.95 nm), this forces the y-axis to always span from 0,
        # compressing real, meaningful variation into a thin sliver at the
        # top of the plot ("can't see small changes"). Disabling auto-range
        # on the ViewBox and instead setting the y-range explicitly from
        # the real curve's data (see _update_curve) fixes this — the fill
        # still renders as a nice gradient, it just no longer influences
        # how the axis scales.
        plot_widget.enableAutoRange(axis='y', enable=False)
        return main_curve

    def _update_curve(self, curve, x, y):
        curve.setData(x, y)
        for sibling in getattr(curve, "_glow_siblings", []):
            sibling.setData(x, y)

        # Explicit y-range from the REAL data only (not the fill's implied
        # extent to zero — see _add_glow_curve for why that matters). A
        # small padding keeps the trace off the plot edges; a minimum span
        # avoids a division-by-near-zero-range look when the data is
        # nearly flat (e.g. right at the very start of a run).
        # X range explicitly as well. Relying on pyqtgraph's auto-range for x
        # means the axis is only fitted when the ViewBox next decides to
        # update, and any plot that ends up with auto-range disabled on x
        # silently keeps its default 0-1 view while the real data sits far
        # off screen — the trace is present and correct but invisible.
        # Setting both axes from the data removes that whole failure mode.
        if len(x) >= 2:
            x_arr = np.asarray(x, dtype=float)
            x_min, x_max = float(np.min(x_arr)), float(np.max(x_arr))
            if x_max > x_min:
                pw = getattr(curve, "_plot_widget", None)
                if pw is not None:
                    pad_x = (x_max - x_min) * 0.02
                    pw.getPlotItem().setXRange(x_min - pad_x, x_max + pad_x, padding=0)

        if len(y) >= 2:
            y_arr = np.asarray(y, dtype=float)
            y_min, y_max = float(np.min(y_arr)), float(np.max(y_arr))
            span = y_max - y_min
            if span < 1e-6:
                span = max(abs(y_max), 1e-3) * 0.1 or 1e-3
            pad = span * 0.15
            plot_widget = getattr(curve, "_plot_widget", None)
            if plot_widget is not None:
                # NOTE: plot_widget.setYRange(...) directly is unsafe —
                # PlotWidget inherits a DIFFERENT setYRange(rect, padding)
                # from GraphicsView that takes a QRectF, not (min, max).
                # Calling it with two floats would raise at runtime.
                # getPlotItem().setYRange(...) forwards correctly to the
                # ViewBox's real (min, max, padding=...) signature —
                # verified directly against the installed pyqtgraph build.
                plot_widget.getPlotItem().setYRange(y_min - pad, y_max + pad, padding=0)

    # ---------------- Simulation control ----------------
    def _start_simulation(self):
        pdb_path = self.pdb_path_edit.text().strip()
        if not pdb_path or not os.path.exists(pdb_path):
            self.log_box.append("Pick a valid complex PDB first.")
            return

        os.makedirs(self.output_dir_edit.text().strip() or "./md_output", exist_ok=True)

        receptor_chains = [c.strip() for c in self.receptor_chains_edit.text().split(",") if c.strip()]
        binder_chains = [c.strip() for c in self.binder_chains_edit.text().split(",") if c.strip()]

        self._steps, self._rmsd, self._rg, self._iface = [], [], [], []
        self._pe, self._ke, self._times = [], [], []
        self._rmsd_rec, self._rmsd_bnd_self, self._rmsd_bnd_on_rec = [], [], []
        self._rmsf, self._rmsf_resids, self._rmsf_labels = [], [], []
        self._rmsf_breaks_drawn = False
        for c in (self.rmsd_rec_curve, self.rmsd_bnd_self_curve,
                  self.rmsd_bnd_on_rec_curve):
            c.setData([], [])
        self._update_curve(self.energy_curve, [], [])
        self._update_curve(self.rmsf_curve, [], [])
        self._render_frame_counter = 0
        self._structure_load_in_flight = False
        self._pending_chains = (receptor_chains, binder_chains)

        self.worker = SimulationWorker(
            pdb_path=pdb_path,
            receptor_chains=receptor_chains,
            binder_chains=binder_chains,
            temperature_k=self.temp_spin.value(),
            total_steps=self.total_steps_spin.value(),
            report_interval=self.report_interval_spin.value(),
            timestep_fs=self.timestep_spin.value(),
            solvent=self.solvent_combo.currentText(),
            strip_heterogens=self.strip_heterogens_check.isChecked(),
            output_dir=self.output_dir_edit.text().strip(),
            platform_name=self.platform_combo.currentText(),
            force_field=self.forcefield_combo.currentText(),
            implicit_model=(self.solvent_model_combo.currentText()
                            if self.solvent_combo.currentText() == "implicit"
                            else DEFAULT_IMPLICIT),
            water_model=(self.solvent_model_combo.currentText()
                         if self.solvent_combo.currentText() == "explicit"
                         else None),
            minimization_steps=self.minimize_spin.value(),
        )
        self.worker.structure_ready.connect(self._on_structure_ready)
        self.worker.frame_ready.connect(self._on_frame)
        self.worker.log.connect(self.log_box.append)
        self.worker.finished_ok.connect(self._on_finished)
        self.worker.failed.connect(self._on_failed)
        self.worker.start()

        self.start_btn.setEnabled(False)
        self.stop_btn.setEnabled(True)
        self.pause_btn.setEnabled(True)
        self.pause_btn.setText("Pause")
        self.save_frame_btn.setEnabled(False)

    def _stop_simulation(self):
        if self.worker:
            self.worker.request_stop()

    def _toggle_pause(self):
        if not self.worker:
            return
        if self.worker.is_paused:
            self.worker.request_resume()
            self.pause_btn.setText("Pause")
        else:
            self.worker.request_pause()
            self.pause_btn.setText("Resume")

    def _save_current_frame(self):
        """Write the currently displayed coordinates out as a PDB.

        The view structure emitted by the worker and the per-frame coordinate
        array are both in protein_idx order, so the frame can be written by
        substituting the coordinate columns of the stored PDB text rather than
        rebuilding a topology. Columns 31-54 (1-indexed) are x/y/z in
        %8.3f — the fixed-column PDB format, so slicing is exact.
        """
        payload = getattr(self, "_last_frame_payload", None)
        pdb_text = getattr(self, "_last_pdb_text", None)
        if payload is None or not pdb_text:
            self.log_box.append("Nothing to save yet — no frame has been received.")
            return

        coords = np.asarray(payload["positions_nm"], dtype=float) * 10.0  # nm -> A
        lines = pdb_text.split("\n")
        atom_rows = [i for i, ln in enumerate(lines) if ln.startswith(("ATOM  ", "HETATM"))]
        if len(atom_rows) != len(coords):
            self.log_box.append(
                f"Refusing to save: view structure has {len(atom_rows)} atoms but "
                f"the frame has {len(coords)}. These must match."
            )
            return

        default = os.path.join(
            self.output_dir_edit.text().strip() or ".",
            f"frame_{payload['step']}.pdb",
        )
        path, _ = QFileDialog.getSaveFileName(
            self, "Save current frame", default, "PDB files (*.pdb)"
        )
        if not path:
            return

        for row, (x, y, z) in zip(atom_rows, coords, strict=True):
            ln = lines[row].ljust(80)
            lines[row] = f"{ln[:30]}{x:8.3f}{y:8.3f}{z:8.3f}{ln[54:]}".rstrip()

        header = (
            f"REMARK   1 Jiggle — frame at step {payload['step']} "
            f"({payload['time_ps']:.3f} ps)\n"
            f"REMARK   1 RMSD {payload['rmsd_nm']:.4f} nm, "
            f"Rg {payload['rg_nm']:.4f} nm"
        )
        try:
            with open(path, "w") as f:
                f.write(header + "\n" + "\n".join(lines).rstrip() + "\n")
        except OSError as exc:
            self.log_box.append(f"Could not write {path}: {exc}")
            return
        self.log_box.append(f"Saved frame at step {payload['step']} to {path}")

    # ---------------- Saving results ----------------
    def _export_results(self, auto=False):
        """Write every metric series and every plot to the output directory.

        Until now the only thing a run left behind was topology.pdb and the
        trajectory — the numbers on screen vanished when the window closed.
        CSV is written first and separately from the images, so a failure in
        the image exporter can never cost you the data.
        """
        if not self._steps:
            self.log_box.append("Nothing to export yet — no frames received.")
            return

        out_dir = self.output_dir_edit.text().strip() or "./md_output"
        try:
            os.makedirs(out_dir, exist_ok=True)
        except OSError as exc:
            self.log_box.append(f"Could not create {out_dir}: {exc}")
            return

        written = []

        # --- per-frame metrics ---
        metrics_path = os.path.join(out_dir, "metrics.csv")
        cols = [("step", self._steps), ("time_ps", self._times),
                ("rmsd_nm", self._rmsd), ("radius_of_gyration_nm", self._rg),
                ("potential_energy_kjmol", self._pe),
                ("kinetic_energy_kjmol", self._ke),
                ("interface_distance_nm", self._iface)]
        if self._rmsd_rec:
            cols += [("rmsd_receptor_nm", self._rmsd_rec),
                     ("rmsd_binder_internal_nm", self._rmsd_bnd_self),
                     ("rmsd_binder_on_receptor_nm", self._rmsd_bnd_on_rec)]
        try:
            with open(metrics_path, "w", newline="") as f:
                w = csv.writer(f)
                w.writerow([c[0] for c in cols])
                for row in zip(*[c[1] for c in cols], strict=False):
                    w.writerow(["" if v is None else v for v in row])
            written.append(metrics_path)
        except OSError as exc:
            self.log_box.append(f"Could not write {metrics_path}: {exc}")

        # --- per-residue RMSF ---
        if self._rmsf:
            rmsf_path = os.path.join(out_dir, "rmsf.csv")
            try:
                with open(rmsf_path, "w", newline="") as f:
                    w = csv.writer(f)
                    w.writerow(["residue", "label", "rmsf_nm"])
                    labels = self._rmsf_labels or [""] * len(self._rmsf)
                    for r, lab, val in zip(self._rmsf_resids, labels, self._rmsf, strict=False):
                        w.writerow([r, lab, f"{val:.6f}"])
                written.append(rmsf_path)
            except OSError as exc:
                self.log_box.append(f"Could not write {rmsf_path}: {exc}")

        # --- plot images ---
        # Skipping the interface plot for a monomer: an empty axis saved to
        # disk is worse than no file, because it looks like a measurement.
        has_binder = bool(self._rmsd_rec)
        for name, plot in self._all_plots:
            if name == "interface_distance" and not has_binder:
                continue
            if name == "rmsf" and not self._rmsf:
                continue
            png = os.path.join(out_dir, f"{name}.png")
            try:
                ex = pg.exporters.ImageExporter(plot.getPlotItem())
                ex.parameters()["width"] = 1400
                ex.export(png)
                written.append(png)
            except Exception as exc:
                self.log_box.append(f"Could not export {png}: {exc}")

        prefix = "Run finished — saved" if auto else "Saved"
        self.log_box.append(
            f"{prefix} {len(written)} file(s) to {out_dir}: "
            + ", ".join(os.path.basename(p) for p in written)
        )

    def _on_structure_ready(self, pdb_text):
        self._load_structure_js(pdb_text)

    def _on_frame(self, payload):
        self._steps.append(payload["step"])
        self._rmsd.append(payload["rmsd_nm"])
        self._rg.append(payload["rg_nm"])
        self._iface.append(payload["interface_distance_nm"] or 0.0)

        self._pe.append(payload.get("potential_energy_kjmol"))
        self._ke.append(payload.get("kinetic_energy_kjmol"))
        self._times.append(payload.get("time_ps"))

        self._update_curve(self.rmsd_curve, self._steps, self._rmsd)
        self._update_curve(self.rg_curve, self._steps, self._rg)
        self._update_curve(self.iface_curve, self._steps, self._iface)
        if self._pe[-1] is not None:
            self._update_curve(self.energy_curve, self._steps, self._pe)

        if payload.get("rmsd_binder_on_receptor_nm") is not None:
            self._rmsd_rec.append(payload["rmsd_receptor_nm"])
            self._rmsd_bnd_self.append(payload["rmsd_binder_self_nm"])
            self._rmsd_bnd_on_rec.append(payload["rmsd_binder_on_receptor_nm"])
            self.rmsd_rec_curve.setData(self._steps, self._rmsd_rec)
            self.rmsd_bnd_self_curve.setData(self._steps, self._rmsd_bnd_self)
            self.rmsd_bnd_on_rec_curve.setData(self._steps, self._rmsd_bnd_on_rec)
            # Autoscale over every trace, not just the complex one.
            allv = self._rmsd + self._rmsd_rec + self._rmsd_bnd_self + self._rmsd_bnd_on_rec
            lo, hi = min(allv), max(allv)
            pad = (hi - lo) * 0.12 or 0.01
            item = self.rmsd_plot.getPlotItem()
            item.setYRange(lo - pad, hi + pad, padding=0)
            if len(self._steps) >= 2:
                sx = (self._steps[-1] - self._steps[0]) * 0.02
                item.setXRange(self._steps[0] - sx, self._steps[-1] + sx, padding=0)

        rmsf = payload.get("rmsf_nm")
        if rmsf is not None and len(rmsf):
            self._rmsf = list(rmsf)
            self._rmsf_resids = payload.get("rmsf_resids") or list(range(1, len(rmsf) + 1))
            self._rmsf_labels = payload.get("rmsf_labels") or []
            # Flat zeros until frame 2 — plotting them would imply a measured
            # value of zero rather than "not enough data yet".
            if payload.get("rmsf_n_frames", 0) >= 2:
                self._update_curve(self.rmsf_curve, self._rmsf_resids, self._rmsf)
                breaks = payload.get("rmsf_chain_breaks") or []
                if breaks and not getattr(self, "_rmsf_breaks_drawn", False):
                    for xb in breaks:
                        self.rmsf_plot.addItem(pg.InfiniteLine(
                            pos=xb - 0.5, angle=90,
                            pen=pg.mkPen("#5a5f6d", width=1, style=Qt.DotLine)))
                    self._rmsf_breaks_drawn = True

        if not self.export_btn.isEnabled():
            self.export_btn.setEnabled(True)

        self._last_frame_payload = payload
        if not self.save_frame_btn.isEnabled():
            self.save_frame_btn.setEnabled(True)

        # Metrics/plots update every reporting frame (cheap), but the 3D
        # cartoon rebuild is decoupled to a lower rate — a full cartoon
        # regeneration (regen:true) on every single frame is real,
        # avoidable load on the most fragile part of this pipeline for no
        # visual benefit at typical reporting intervals. Rendering every
        # 5th frame still looks smooth and cuts rebuild count by 80%.
        #
        # Also gated on _structure_load_in_flight: a real run showed a
        # coordinate update reaching the viewer before the structure-load
        # it depended on had taken effect, corrupting the atom count. Any
        # frame that arrives while a load is still in flight is simply
        # skipped here — _last_frame_payload still gets updated above, so
        # _on_load_confirmed's catch-up render (in _load_structure_js)
        # sends the freshest coordinates as soon as it's safe to.
        self._render_frame_counter = getattr(self, "_render_frame_counter", 0) + 1
        if self._render_frame_counter % 5 == 0 and not self._structure_load_in_flight:
            self._render_frame_now(payload)

    def _on_finished(self):
        # Force a final render of the true last frame — throttling above
        # means the last-rendered frame could otherwise be up to 4
        # reporting intervals stale relative to what the metrics/DCD
        # actually ended on. Skipped if a load is still in flight for the
        # same reason as above; _on_load_confirmed's catch-up will cover it.
        if getattr(self, "_last_frame_payload", None) is not None and not self._structure_load_in_flight:
            self._render_frame_now(self._last_frame_payload)
        self.start_btn.setEnabled(True)
        self.stop_btn.setEnabled(False)
        self.pause_btn.setEnabled(False)
        self.pause_btn.setText("Pause")
        self.log_box.append("Simulation finished.")
        self._export_results(auto=True)

    def _on_failed(self, msg):
        self.start_btn.setEnabled(True)
        self.stop_btn.setEnabled(False)
        self.pause_btn.setEnabled(False)
        self.pause_btn.setText("Pause")
        self.log_box.append(f"ERROR: {msg}")
