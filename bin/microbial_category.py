#!/usr/bin/env python3
"""Microbial category vocabulary shared by the TaxTriage scripts.

Every organism carries one of five categories: Primary, Opportunistic,
Potential, Commensal, or Uncategorized. "Uncategorized" means the organism
(and no ancestor in its lineage) is listed in the pathogen sheet.

Anything that is not one of the four curated categories -- a blank value, a
typo in the sheet, or the legacy label "Unknown" written by runs before
2026-11 -- is reported as "Uncategorized".
"""

UNCATEGORIZED = "Uncategorized"
MICROBIAL_CATEGORIES = ("Primary", "Opportunistic", "Potential", "Commensal", UNCATEGORIZED)

# Fields that hold a category: match_paths JSON (microbial_category,
# is_pathogen) and report rows (Microbial Category).
CATEGORY_KEYS = ("microbial_category", "is_pathogen", "Microbial Category")

_ALIASES = {
    "primary": "Primary",
    "primary (sterile)": "Primary",
    "opportunistic": "Opportunistic",
    "potential": "Potential",
    "commensal": "Commensal",
    "uncategorized": UNCATEGORIZED,
}


def normalize_category(label):
    """Return the standard category for ``label``.

    Matching ignores case and surrounding whitespace. Blank, unrecognised or
    legacy ("Unknown") labels return "Uncategorized".
    """
    if label is None or isinstance(label, bool):
        return UNCATEGORIZED
    return _ALIASES.get(str(label).strip().lower(), UNCATEGORIZED)


def canonicalize_categories(obj, keys=CATEGORY_KEYS):
    """Normalise every category field in a JSON-like tree, in place; returns ``obj``.

    Only string or null values of the given keys are touched, so a boolean or
    numeric field that happens to share a name is left alone.
    """
    if isinstance(obj, dict):
        for k, v in obj.items():
            if k in keys and (v is None or isinstance(v, str)):
                obj[k] = normalize_category(v)
            elif isinstance(v, (dict, list)):
                canonicalize_categories(v, keys)
    elif isinstance(obj, list):
        for v in obj:
            canonicalize_categories(v, keys)
    return obj
