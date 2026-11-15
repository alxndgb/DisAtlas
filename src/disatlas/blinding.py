"""Stable blinding of licence-restricted identifiers.

DisAtlas ships Open Targets evidence with real target identifiers. Evidence from
licence-restricted sources (METABASE, CTD and NLP-mined literature, pathology-text and
knowledge-graph derived channels) ships with blinded target identifiers: a
stable salted-hash token ``T_<12 hex>``. The same raw identifier always maps to
the same token in every file, so shared-target structure (and therefore the
disease-disease graph built from it) is preserved exactly, while the raw
identifier cannot be recovered without the salt.

The salt is NOT part of the repository. Shipped data under ``data/raw`` are
already blinded, so no user needs the salt to reproduce any result. The salt is
only needed by the maintainers when re-exporting raw channels; it is read from
the environment variable ``DISATLAS_BLIND_SALT`` or from the file named by
``DISATLAS_BLIND_SALT_FILE``.

Rules (see data/README.md):
  * ``blind_target(raw_id)`` -> ``'T_' + sha256(salt + '|' + raw_id)[:12]``.
  * ``is_restricted(source)`` decides per evidence source whether targets are
    blinded; Open Targets sources (``OT_*``) are not restricted.
  * ``source_group(source)`` maps an internal source label to the public
    channel/group label that is shipped.
"""
from __future__ import annotations

import hashlib
import os
import re
from typing import Iterable

TOKEN_PREFIX = "T_"
TOKEN_HEX = 12
TOKEN_RE = re.compile(r"^T_[0-9a-f]{12}$")

# Public group label for each internal evidence-source label. Anything not
# listed and not starting with ``OT_`` is treated as restricted and grouped as
# ``restricted_other``.
_GROUPS = {
    "METABASE": "curated_mechanism",
    "CTD_marker_mechanism": "curated_toxicogenomic",
    "LIT_SCI": "literature",
    "LIT_PREPRINT": "literature",
    "LIT_SCI+LIT_PREPRINT": "literature",
    "LIT_MECH": "literature",
    "PATH_NLP_ANATOMY": "nlp_pathology",
    "PATH_NLP_PHENOTYPE": "nlp_pathology",
    "PATH_NLP_CELLTYPE": "nlp_pathology",
    "AUTOAB_NLP_GENE": "nlp_autoantibody",
}


# Public evidence sources: identifiers ship raw.
PUBLIC_SOURCES = {"OPENTARGETS", "ORPHANET", "OMIM", "OMIM DIRECT", "GENCC", "OT", "REACTOME", "HPO"}

# Evidence-family labels used in the specific-channel module-support table.
FAMILY_GROUPS = {
    "OT": "open_targets",
    "LIT": "literature",
    "METABASE": "curated_restricted",
    "CTD": "curated_restricted",
    "AUTOAB": "literature_autoantibody",
    "PATH_ANAT": "pathology_text",
    "PATH_CELL": "pathology_text",
    "PATH_PHEN": "pathology_text",
    "MECHNODE": "mechanism_node",
}

MODULE_PREFIX = "M_"
MODULE_RE = re.compile(r"^M_[0-9a-f]{12}$")
# Module namespaces whose identifiers ship raw (public vocabularies / author-curated).
PUBLIC_MODULE_NS = ("REACTOME:", "NAMED:", "MECHNODE:")


def is_restricted(source: str) -> bool:
    """True if targets from this evidence source must be blinded."""
    s = str(source).upper()
    if s.startswith(("OT_", "CURATED_OT_")):
        return False
    return s not in PUBLIC_SOURCES


def source_group(source: str) -> str:
    """Public channel/group label for an internal evidence-source label."""
    s = str(source)
    if s.startswith("OT_"):
        return s  # Open Targets datatype names are public
    if s.startswith("CURATED_OT_"):
        return "OT_" + s[len("CURATED_OT_"):]
    if s.upper() in PUBLIC_SOURCES:
        return s.split()[0].lower()
    if s.startswith("BIKG_NLP_"):
        return "literature"
    if s.startswith("CURATED_CTD"):
        return "curated_restricted"
    return _GROUPS.get(s, "restricted_other")


def family_groups(families: str, sep: str = "|") -> str:
    """Map an evidence-family string (e.g. 'CTD|LIT|OT') to sorted unique public groups."""
    if families is None or families != families or families == "":
        return ""
    return sep.join(sorted({FAMILY_GROUPS.get(f, "restricted_other") for f in str(families).split(sep)}))


def blind_module(module_id: str, salt: str | None = None) -> str:
    """Blind a literature/pathology-derived module identifier ('M_<12hex>').
    REACTOME, NAMED and MECHNODE modules stay raw."""
    m = str(module_id)
    if m.startswith(PUBLIC_MODULE_NS) or MODULE_RE.match(m):
        return m
    s = _salt() if salt is None else salt
    return MODULE_PREFIX + hashlib.sha256(f"{s}|module|{m}".encode()).hexdigest()[:TOKEN_HEX]


def _salt() -> str:
    s = os.environ.get("DISATLAS_BLIND_SALT")
    if not s:
        f = os.environ.get("DISATLAS_BLIND_SALT_FILE")
        if f and os.path.exists(f):
            with open(f) as fh:
                s = fh.read().strip()
    if not s:
        raise RuntimeError(
            "blinding salt not available: shipped data are already blinded; "
            "set DISATLAS_BLIND_SALT or DISATLAS_BLIND_SALT_FILE only to re-export raw channels")
    return s


def blind_target(raw_id: str, salt: str | None = None) -> str:
    """Stable token for one raw target identifier (idempotent on tokens)."""
    r = str(raw_id)
    if TOKEN_RE.match(r):
        return r
    s = _salt() if salt is None else salt
    return TOKEN_PREFIX + hashlib.sha256(f"{s}|{r}".encode()).hexdigest()[:TOKEN_HEX]


def blind_series(ids: Iterable[str], salt: str | None = None) -> list[str]:
    """Vectorised blind_target with a per-call cache (unique ids hashed once)."""
    s = _salt() if salt is None else salt
    cache: dict[str, str] = {}
    out = []
    for i in ids:
        t = cache.get(i)
        if t is None:
            t = cache[i] = blind_target(i, s)
        out.append(t)
    return out


def blind_frame(df, target_col: str = "target", source_col: str = "source", salt: str | None = None):
    """Return a copy of df with restricted-source targets blinded and the source
    column replaced by its public group label (column ``source_group``)."""
    out = df.copy()
    restricted = out[source_col].map(is_restricted).to_numpy()
    if restricted.any():
        out.loc[restricted, target_col] = blind_series(out.loc[restricted, target_col].astype(str), salt)
    out["source_group"] = out[source_col].map(source_group)
    return out.drop(columns=[source_col])


def is_blinded(x: str) -> bool:
    return bool(TOKEN_RE.match(str(x)) or MODULE_RE.match(str(x)))
