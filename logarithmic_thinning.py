"""Public API for the logarithmic-thinning package.

This module exists so that the import path matches the distribution name:

    from logarithmic_thinning import logarithmic_thinning

It is a thin re-export of :mod:`thin_sorted_pvalues`, which keeps its own name because
``run_pipeline.sh`` invokes it by path and the test suite imports it directly. Both import
paths are supported and refer to the same objects.

The thinning rule
-----------------
Given ``total_rows`` p-values already sorted most-significant-first, return the 1-based
ranks to keep: every rank above ``factor / (factor - 1)`` is kept unthinned, and below that
the ranks are sampled geometrically by repeatedly dividing by ``factor``. A Manhattan or QQ
plot built from the kept ranks is visually indistinguishable from one built from all rows,
at a fraction of the size, because nothing near the top of the distribution is discarded.

    >>> from logarithmic_thinning import logarithmic_thinning
    >>> ranks = logarithmic_thinning(7_834_624, 1.0005)
    >>> len(ranks)
    17705

This function is pure and uses only the standard library, so its output is deterministic
across Python versions, platforms and numerical-library versions -- there is nothing here
whose version could change a result.
"""

from thin_sorted_pvalues import (  # noqa: F401
    UNTHINNED_COLUMN,
    determine_threshold,
    logarithmic_thinning,
)

__all__ = [
    "logarithmic_thinning",
    "determine_threshold",
    "UNTHINNED_COLUMN",
]

__version__ = "1.0.0"
