"""Tests that jobs are built by each biobank's workflow interface, not its identity.

Step 1, step 2 and conditioning jobs and their inputs follow
`BIOBANK_PROFILES[bb].workflow_interface`, so a new biobank on an existing platform
and interface needs only a `BioBank` member, a profile row and config.
"""

from __future__ import annotations

import pytest

from fedcat.core.engine import _step2_block_length


def test_step2_block_length_keys_on_interface():
    size = 1000

    # Singles: dnanexus uses a larger numerator than gcs → larger blocks.
    wdl_singles = _step2_block_length("wdl_scatter_v1", True, False, size)
    dx_singles = _step2_block_length("dx_applet_v1", True, False, size)
    assert wdl_singles == int(max(1.5e12 / size, 6_000_000))
    assert dx_singles == int(max(3e12 / size, 6_000_000))
    assert dx_singles > wdl_singles

    # Aggregates: coding_only pins a fixed size (differs by interface).
    assert _step2_block_length("wdl_scatter_v1", False, True, size) == 20_000_000
    assert _step2_block_length("dx_applet_v1", False, True, size) == 30_000_000

    # Aggregates, not coding-only: bounded window.
    assert _step2_block_length("wdl_scatter_v1", False, False, size) == int(
        max(min(1.5e11 / size, 20_000_000), 600_000)
    )
    assert _step2_block_length("dx_applet_v1", False, False, size) == int(
        max(min(3e11 / size, 30_000_000), 600_000)
    )


def test_step2_block_length_rejects_unknown_interface():
    # An unhandled interface must fail loudly (case _), not silently take a branch.
    with pytest.raises(ValueError, match="unsupported workflow_interface"):
        _step2_block_length("some_future_interface", True, False, 1000)
