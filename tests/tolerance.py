"""Registry-wide author tolerance for decimal values (AUROC, percentile, share, ratio, delta).

A computed row with op `dtol` passes when the printed decimal differs from the regenerated value by at most
DECIMAL_TOL (absolute, in the printed unit).  It is not used for counts, p-values, whole-percent figures or
odds ratios; rows that use it carry a `tol_note:` in their note giving both values.  See docs/ERRATA.md.
"""
DECIMAL_TOL = 0.003


def dtol_ok(v, t, token):
    assert "." in token, f"dtol is only for decimal tokens: {token}"
    return abs(v - t) <= DECIMAL_TOL + 1e-12
