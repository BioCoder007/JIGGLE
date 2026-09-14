"""
Force-field catalogue.

Every entry here was verified against the installed OpenMM by actually
building a System from it — not by reading the docs. Two things that are easy
to get wrong and are handled below:

  * 4-site waters (TIP4P-Ew) fail with a cryptic "No template found for
    residue N (HOH)" unless the water model NAME is passed to
    Modeller.addSolvent(model=...). The default is 'tip3p', which produces
    3-site waters that then have no template in a 4-site force field.
  * Modeller.addSolvent only knows these model names:
        tip3p, spce, tip4pew, tip5p, swm4ndp
    TIP3P-FB is a 3-site water, so it is built with model='tip3p' and
    parameterised by amber14/tip3pfb.xml. OPC has no addSolvent support in
    this OpenMM build at all, so it is deliberately absent.

The GUI reads this module to populate its dropdowns, so it can never offer a
combination the worker cannot build.
"""

# Implicit-solvent (Generalised Born) models. All of these were confirmed to
# build against every protein force field listed below, CHARMM36 included.
IMPLICIT_MODELS = {
    "OBC2 (GBSA, default)": "implicit/obc2.xml",
    "OBC1":                 "implicit/obc1.xml",
    "GBn2":                 "implicit/gbn2.xml",
    "GBn":                  "implicit/gbn.xml",
    "HCT":                  "implicit/hct.xml",
}

# Water models, keyed by display name -> (xml, addSolvent model name).
# The xml path differs per protein force field, so it is resolved through
# each entry's own "water" table below.
_AMBER14_WATERS = {
    "TIP3P-FB (default)": ("amber14/tip3pfb.xml", "tip3p"),
    "TIP3P":              ("amber14/tip3p.xml",   "tip3p"),
    "SPC/E":              ("amber14/spce.xml",    "spce"),
    "TIP4P-Ew":           ("amber14/tip4pew.xml", "tip4pew"),
}

_CLASSIC_AMBER_WATERS = {
    "TIP3P (default)": ("tip3p.xml",   "tip3p"),
    "SPC/E":           ("spce.xml",    "spce"),
    "TIP4P-Ew":        ("tip4pew.xml", "tip4pew"),
}

_CHARMM_WATERS = {
    "CHARMM TIP3P (default)": ("charmm36/water.xml",    "tip3p"),
    "SPC/E":                  ("charmm36/spce.xml",     "spce"),
    "TIP4P-Ew":               ("charmm36/tip4pew.xml",  "tip4pew"),
}

FORCE_FIELDS = {
    "AMBER14 (ff14SB) — default": {
        "protein": "amber14-all.xml",
        "water": _AMBER14_WATERS,
        "note": "Well-tested general-purpose protein force field. Good default.",
    },
    "AMBER19 (ff19SB)": {
        "protein": "amber19-all.xml",
        "water": _AMBER14_WATERS,
        "note": "Newer backbone parameters. Intended for use with OPC water, "
                "which this OpenMM build cannot solvate with — TIP3P-FB is the "
                "closest available substitute.",
    },
    "AMBER99SB-ILDN": {
        "protein": "amber99sbildn.xml",
        "water": _CLASSIC_AMBER_WATERS,
        "note": "Long-standing benchmark force field; widely used in the "
                "literature, so useful for reproducing published work.",
    },
    "AMBER99SB": {
        "protein": "amber99sb.xml",
        "water": _CLASSIC_AMBER_WATERS,
        "note": "Predecessor to ILDN.",
    },
    "AMBER10": {
        "protein": "amber10.xml",
        "water": _CLASSIC_AMBER_WATERS,
        "note": "Legacy.",
    },
    "AMBER03": {
        "protein": "amber03.xml",
        "water": _CLASSIC_AMBER_WATERS,
        "note": "Legacy.",
    },
    "CHARMM36": {
        "protein": "charmm36.xml",
        "water": _CHARMM_WATERS,
        "note": "Independent parameter lineage from AMBER — useful as a "
                "cross-check when a result looks force-field-dependent. Note "
                "that charmm36.xml is a large file and takes noticeably longer "
                "to parse at startup than the AMBER sets.",
    },
}

DEFAULT_FORCE_FIELD = "AMBER14 (ff14SB) — default"
DEFAULT_IMPLICIT = "OBC2 (GBSA, default)"


def force_field_names():
    return list(FORCE_FIELDS.keys())


def water_names(force_field_name):
    return list(FORCE_FIELDS[force_field_name]["water"].keys())


def default_water_name(force_field_name):
    return water_names(force_field_name)[0]


def resolve(force_field_name, solvent, implicit_name=None, water_name=None):
    """Return (xml_files, water_model_name).

    xml_files is the list passed to app.ForceField(...).
    water_model_name is passed to Modeller.addSolvent(model=...); it is None
    for implicit solvent.
    """
    entry = FORCE_FIELDS[force_field_name]
    if solvent == "implicit":
        implicit_xml = IMPLICIT_MODELS[implicit_name or DEFAULT_IMPLICIT]
        return [entry["protein"], implicit_xml], None

    waters = entry["water"]
    key = water_name if water_name in waters else default_water_name(force_field_name)
    water_xml, model = waters[key]
    return [entry["protein"], water_xml], model
