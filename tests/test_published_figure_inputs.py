"""Regression test: the SELECTION RULE behind a published figure.

The Manhattan panel of Figure 3a in the retinal-vasculometry GWAS paper is drawn from three
thinned SNP tables -- ``{mtifs,dtifs,lvs}_thinned_snps.csv`` -- produced by this package.

SCOPE -- read this before citing the test as evidence of reproducibility.

What it verifies: *which rows* thinning keeps. That is a pure function of two numbers, the
row count before thinning and the factor; no p-values, positions or genotypes enter into it.
So the exact published rank selection regenerates in milliseconds, on any machine, with no
access to UK Biobank data.

What it does NOT verify: the *contents* of those rows -- ``chr``, ``pos``, ``pp``,
``original_index``. Producing those means running the pipeline end to end (read the GWAS
outputs, threshold, sort into chunks, k-way merge into one global ordering, then select),
which takes minutes for the TIF sets and up to hours for the 1024 latent variables. Nothing
here exercises the sort or the merge; ``test_large_pipeline_consistency`` covers those on
synthetic data.

So: this test proves the algorithm has not drifted. It does not prove the pipeline still
reproduces the published files. Those are different claims and only the first is made here.

Each case is pinned to the SHA-256 of the comma-joined kept ranks; if the rule ever changes,
these fail and name the figure that would silently change with it.

Cross-checked against the real files on 2026-08-13: the regenerated rank list is identical,
element for element, to the ranks recorded in all three published CSVs -- so the published
tables were selected by exactly this rule at factor 1.0003. See
``examples/reproduce_figure3a_ranks.py`` to repeat that check where the files are available.
"""

import hashlib

from logarithmic_thinning import determine_threshold, logarithmic_thinning

# (label, rows before thinning, factor, rows kept, sha256 of the kept ranks)
PUBLISHED_FIGURE3A_INPUTS = [
    ("mtifs", 2_832_855, 1.0003, 24_412,
     "030586d9e2027358cdaeeb4b2d7ff904ce7af0ccf84d44923bb8f2d15e0760a4"),
    ("dtifs", 2_914_788, 1.0003, 24_507,
     "74d73ff3b8214cad55eda92f0aedc15c0c5ccc4cb188709324d37e11a409c325"),
    ("lvs", 171_151_566, 1.0003, 38_083,
     "4e39ceb539ebed64a2d0431d67aa52a80ce6913653b9a4e70a107ed798448b38"),
]


def _kept_ranks(total_rows, factor):
    return [r for r in logarithmic_thinning(total_rows, factor) if 1 <= r <= total_rows]


def test_reproduces_published_figure3a_thinning():
    for label, total_rows, factor, expected_kept, expected_sha in PUBLISHED_FIGURE3A_INPUTS:
        ranks = _kept_ranks(total_rows, factor)

        assert len(ranks) == expected_kept, (
            f"{label}: kept {len(ranks)} rows, published figure used {expected_kept}"
        )
        digest = hashlib.sha256(",".join(map(str, ranks)).encode()).hexdigest()
        assert digest == expected_sha, (
            f"{label}: the kept-rank set changed. Figure 3a would no longer reproduce."
        )


def test_thinning_keeps_the_strongest_signals_unthinned():
    """The top of the distribution must survive intact -- that is the point of the method."""
    for label, total_rows, factor, _, _ in PUBLISHED_FIGURE3A_INPUTS:
        ranks = _kept_ranks(total_rows, factor)
        head = determine_threshold(factor)  # 3334 for factor=1.0003

        assert ranks[: head - 1] == list(range(1, head)), (
            f"{label}: the {head - 1} most significant rows are no longer kept consecutively"
        )
        assert ranks[0] == 1, f"{label}: the most significant row was dropped"
        assert ranks[-1] == total_rows, f"{label}: the least significant row was dropped"


def test_thinning_is_deterministic():
    """Same inputs, same output -- no RNG, no dict ordering, no platform dependence."""
    for _, total_rows, factor, _, _ in PUBLISHED_FIGURE3A_INPUTS:
        assert _kept_ranks(total_rows, factor) == _kept_ranks(total_rows, factor)


if __name__ == "__main__":
    test_reproduces_published_figure3a_thinning()
    test_thinning_keeps_the_strongest_signals_unthinned()
    test_thinning_is_deterministic()
    print("Published Figure 3a thinning reproduced exactly.")
