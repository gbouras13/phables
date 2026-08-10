#!/usr/bin/env python3

"""Regression test for workflow/scripts/format_koverage_results.py's parsing of
Koverage's sample_coverage.tsv.

Context
-------
postprocess.smk's `koverage_genomes` rule runs `koverage run --reads ... --ref ...`
(no `coverm` subcommand), which is Koverage's *native* "map" mode. That mode's
sample_coverage.tsv header is:

    Sample  Contig  Count  RPM  RPKM  RPK  TPM  Mean  Median  Hitrate  Variance
      0       1       2     3    4     5    6    7      8        9        10

This differs from Koverage's `coverm` subcommand mode header:

    Sample  Contig  Count  RPKM  TPM  Mean  Covered_fraction  Variance
      0       1       2     3    4     5           6             7

It is easy to mix these two up (both are valid Koverage outputs, just from
different subcommands) and silently mislabel RPKM/TPM or Mean/Variance in
phables' sample_genome_rpkm.tsv / sample_genome_mean_coverage.tsv reports.
This test locks in the native "map" mode column indices that
format_koverage_results.py actually relies on, using a header + data row
copied verbatim from Koverage's own coverage.smk / sampleCoverage.py.

Run directly with:
    python3 tests/test_format_koverage_results.py
"""

import importlib.util
import os
import sys

REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SCRIPT_PATH = os.path.join(
    REPO_ROOT, "workflow", "scripts", "format_koverage_results.py"
)

# Import format_koverage_results.py as a module without executing its
# `if __name__ == "__main__":` block (which requires a snakemake object).
spec = importlib.util.spec_from_file_location(
    "format_koverage_results", SCRIPT_PATH
)
format_koverage_results = importlib.util.module_from_spec(spec)
spec.loader.exec_module(format_koverage_results)


# A native "map" mode header, verbatim from Koverage's
# workflow/rules/coverage.smk `all_sample_coverage` rule.
NATIVE_MAP_MODE_HEADER = (
    "Sample\tContig\tCount\tRPM\tRPKM\tRPK\tTPM\tMean\tMedian\tHitrate\tVariance"
)

# A corresponding data row with distinguishable values for every column, in
# the same "{:.4g}"-formatted style sampleCoverage.py writes.
NATIVE_MAP_MODE_ROW = (
    "sample1\tgenome_A\t100\t11.11\t22.22\t33.33\t44.44\t55.55\t66.66\t0.9\t77.77"
)


def test_parse_koverage_row_uses_native_map_mode_indices():
    strings = NATIVE_MAP_MODE_ROW.split("\t")
    sample, contig, count, rpkm_val, mean_val = (
        format_koverage_results.parse_koverage_row(strings)
    )

    assert sample == "sample1", f"expected sample='sample1', got {sample!r}"
    assert contig == "genome_A", f"expected contig='genome_A', got {contig!r}"
    assert count == 100, f"expected Count=100, got {count!r}"
    # Column 4 in native "map" mode is RPKM (not TPM, which is column 6).
    assert rpkm_val == 22.22, f"expected RPKM=22.22, got {rpkm_val!r}"
    # Column 7 in native "map" mode is Mean (not Variance, which is column 10).
    assert mean_val == 55.55, f"expected Mean=55.55, got {mean_val!r}"


def test_header_matches_expected_native_map_mode_layout():
    header_fields = NATIVE_MAP_MODE_HEADER.split("\t")
    assert header_fields[format_koverage_results.IDX_SAMPLE] == "Sample"
    assert header_fields[format_koverage_results.IDX_CONTIG] == "Contig"
    assert header_fields[format_koverage_results.IDX_COUNT] == "Count"
    assert header_fields[format_koverage_results.IDX_RPKM] == "RPKM"
    assert header_fields[format_koverage_results.IDX_MEAN] == "Mean"


def main():
    tests = [
        test_parse_koverage_row_uses_native_map_mode_indices,
        test_header_matches_expected_native_map_mode_layout,
    ]
    failures = 0
    for test in tests:
        try:
            test()
            print(f"PASS: {test.__name__}")
        except AssertionError as e:
            failures += 1
            print(f"FAIL: {test.__name__}: {e}")

    if failures:
        print(f"\n{failures} test(s) failed")
        sys.exit(1)
    else:
        print("\nAll tests passed")


if __name__ == "__main__":
    main()
