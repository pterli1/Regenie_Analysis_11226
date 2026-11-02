"""Tests for `engine.custom_blocks_array`, which groups variant-density blocks into step 2 blocks.

Expected outputs are worked out by hand, covering block grouping, filling gaps,
chromosome filtering, `no_pgs_chroms` and regions mode.

Run with::

    pytest regenie_workflows/federated_analysis/fedcat/tests/test_custom_blocks_array.py -v
"""

from __future__ import annotations

import pytest

# parents: test_custom_blocks_array.py → tests/ → fedcat/ → regenie_workflows/
from fedcat.core.engine import custom_blocks_array


# 100k blocks on chr1: [0,99999], [100000,199999], [200000,299999], [300000,399999]
def _chr1_blocks(n: int) -> list[tuple[str, int, int]]:
    """Return `n` consecutive 100 kb blocks on chromosome 1."""
    return [("1", i * 100_000, i * 100_000 + 99_999) for i in range(n)]


class TestNonRegions:
    """Whole-chromosome blocks."""

    def test_one_block_per_group_no_fill(self):
        var_blocks = _chr1_blocks(3)
        out = custom_blocks_array(
            var_blocks,
            target_block_length=100_000,  # blocks_per_group = 1
            prefix="pgs_",
            no_pgs_chroms=[],
            chroms_to_run=["1"],
            fill_gaps=False,
        )
        # starts[0] is decremented by 1 (→0); each group is one 100k block.
        assert out == [
            ["1", "0", "99999", "pgs_1"],
            ["1", "100001", "199999", "pgs_1"],
            ["1", "200001", "299999", "pgs_1"],
        ]

    def test_fill_gaps_clamps_first_and_last_and_closes_gaps(self):
        var_blocks = _chr1_blocks(3)
        out = custom_blocks_array(
            var_blocks,
            target_block_length=100_000,
            prefix="pgs_",
            no_pgs_chroms=[],
            chroms_to_run=["1"],
            fill_gaps=True,
        )
        # first START→1, gaps closed to prev STOP+1, last STOP→300_000_000
        assert out == [
            ["1", "1", "99999", "pgs_1"],
            ["1", "100000", "199999", "pgs_1"],
            ["1", "200000", "300000000", "pgs_1"],
        ]

    def test_blocks_per_group_grouping(self):
        # target 200k → 2 blocks per group; 4 blocks → 2 groups
        var_blocks = _chr1_blocks(4)
        out = custom_blocks_array(
            var_blocks,
            target_block_length=200_000,
            prefix="pgs_",
            no_pgs_chroms=[],
            chroms_to_run=["1"],
            fill_gaps=False,
        )
        assert out == [
            ["1", "0", "199999", "pgs_1"],
            ["1", "200001", "399999", "pgs_1"],
        ]

    def test_chrom_filter_and_no_pgs(self):
        var_blocks = [("1", 0, 99_999), ("2", 0, 99_999)]
        out = custom_blocks_array(
            var_blocks,
            target_block_length=100_000,
            prefix="pgs_",
            no_pgs_chroms=["1"],
            chroms_to_run=["1"],  # chr2 filtered out
            fill_gaps=False,
        )
        # only chr1; no_pgs → EXCLCOVARS "None"
        assert out == [["1", "0", "99999", "None"]]

    def test_multi_chrom_fill_gaps_independent_per_chrom(self):
        var_blocks = _chr1_blocks(2) + [("2", 0, 99_999), ("2", 100_000, 199_999)]
        out = custom_blocks_array(
            var_blocks,
            target_block_length=100_000,
            prefix="pgs_",
            no_pgs_chroms=[],
            chroms_to_run=["1", "2"],
            fill_gaps=True,
        )
        # each chrom independently: first START→1, last STOP→300_000_000
        assert out == [
            ["1", "1", "99999", "pgs_1"],
            ["1", "100000", "300000000", "pgs_1"],
            ["2", "1", "99999", "pgs_2"],
            ["2", "100000", "300000000", "pgs_2"],
        ]


class TestRegionsMode:
    """Blocks restricted to `regions`."""

    def test_single_region_fill_gaps(self):
        var_blocks = _chr1_blocks(3)
        out = custom_blocks_array(
            var_blocks,
            target_block_length=100_000,
            prefix="pgs_",
            no_pgs_chroms=[],
            chroms_to_run=["1"],
            fill_gaps=True,
            regions=[("1", 50_000, 250_000)],
        )
        # region clamps first START→region_start, last STOP→region_end; gaps closed
        assert out == [
            ["1", "50000", "99999", "pgs_1"],
            ["1", "100000", "199999", "pgs_1"],
            ["1", "200000", "250000", "pgs_1"],
        ]

    def test_region_chrom_not_in_run_is_skipped(self):
        var_blocks = _chr1_blocks(2)
        out = custom_blocks_array(
            var_blocks,
            target_block_length=100_000,
            prefix="pgs_",
            no_pgs_chroms=[],
            chroms_to_run=["1"],
            fill_gaps=True,
            regions=[("2", 0, 100_000)],  # chr2 not in chroms_to_run
        )
        assert out == []

    def test_region_with_no_overlapping_blocks_skipped(self):
        var_blocks = _chr1_blocks(2)  # spans 0..199999
        out = custom_blocks_array(
            var_blocks,
            target_block_length=100_000,
            prefix="pgs_",
            no_pgs_chroms=[],
            chroms_to_run=["1"],
            fill_gaps=True,
            regions=[("1", 500_000, 600_000)],  # beyond all blocks
        )
        assert out == []


class TestExcludeCovars:
    """``exclude_covars`` / ``exclude_covars_dict`` stack additively with the PGS column."""

    @pytest.mark.parametrize(
        "kwargs, expected",
        [
            ({"prefix": "pgs_", "no_pgs_chroms": ["X"]}, ["pgs_1", "pgs_1", "None"]),
            ({"prefix": None, "no_pgs_chroms": []}, ["None", "None", "None"]),
            ({"prefix": None, "no_pgs_chroms": [], "exclude_covars": ["batch", "centre"]}, ["batch,centre"] * 3),
            (
                {"prefix": "pgs_", "no_pgs_chroms": [], "exclude_covars": ["batch"]},
                ["pgs_1,batch", "pgs_1,batch", "pgs_X,batch"],
            ),
            # per-chromosome exclusions apply only on their chromosome
            (
                {"prefix": None, "no_pgs_chroms": [], "exclude_covars_dict": {"X": ["sex"], "2": ["ignored"]}},
                ["None", "None", "sex"],
            ),
            (
                {
                    "prefix": "pgs_",
                    "no_pgs_chroms": [],
                    "exclude_covars": ["batch"],
                    "exclude_covars_dict": {"1": ["hla"], "X": ["sex"]},
                },
                ["pgs_1,batch,hla", "pgs_1,batch,hla", "pgs_X,batch,sex"],
            ),
            # a no-PGS chromosome still gets the other exclusions
            (
                {
                    "prefix": "pgs_",
                    "no_pgs_chroms": ["X"],
                    "exclude_covars": ["batch"],
                    "exclude_covars_dict": {"X": ["sex"]},
                },
                ["pgs_1,batch", "pgs_1,batch", "batch,sex"],
            ),
            # duplicates collapse in order: pgs_1 is the chr1 PGS column but a distinct name on chrX
            (
                {
                    "prefix": "pgs_",
                    "no_pgs_chroms": [],
                    "exclude_covars": ["batch", "pgs_1"],
                    "exclude_covars_dict": {"1": ["batch"]},
                },
                ["pgs_1,batch", "pgs_1,batch", "pgs_X,batch,pgs_1"],
            ),
        ],
    )
    def test_exclusions_per_block(self, kwargs, expected):
        """Exclusion strings for two chr1 blocks and one chrX block."""
        var_blocks = _chr1_blocks(2) + [("X", 0, 99_999)]
        out = custom_blocks_array(
            var_blocks, target_block_length=100_000, chroms_to_run=["1", "X"], fill_gaps=False, **kwargs
        )
        assert [row[3] for row in out] == expected

    def test_regions_mode_applies_the_same_exclusions(self):
        out = custom_blocks_array(
            _chr1_blocks(2),
            target_block_length=100_000,
            prefix="pgs_",
            no_pgs_chroms=[],
            chroms_to_run=["1"],
            fill_gaps=False,
            regions=[("1", 1, 199_999)],
            exclude_covars=["batch"],
            exclude_covars_dict={"1": ["hla"]},
        )
        assert [row[3] for row in out] == ["pgs_1,batch,hla", "pgs_1,batch,hla"]
