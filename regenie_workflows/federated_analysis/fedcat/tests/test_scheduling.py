"""Tests for the flexible ``--delay-start`` scheduling parser in ``fedcat.__main__``.

``_parse_when`` delegates to dateparser (with a fixed ``RELATIVE_BASE`` so the
future-preference roll is deterministic), covering both relative delays and
absolute times behind ``discover --delay-start``.

Run with::

    pytest regenie_workflows/federated_analysis/fedcat/tests/test_scheduling.py -v
"""

from __future__ import annotations

from datetime import datetime

import pytest

from fedcat.__main__ import _parse_when

NOW = datetime(2026, 7, 21, 10, 0, 0)


@pytest.mark.parametrize(
    "text, expected",
    [
        # Relative delays (anchored to NOW).
        ("in 2 hours", datetime(2026, 7, 21, 12, 0, 0)),
        ("2h", datetime(2026, 7, 21, 12, 0, 0)),
        ("1h30m", datetime(2026, 7, 21, 11, 30, 0)),
        ("90m", datetime(2026, 7, 21, 11, 30, 0)),
        ("30 min", datetime(2026, 7, 21, 10, 30, 0)),
        # Absolute times today.
        ("14:30", datetime(2026, 7, 21, 14, 30, 0)),
        ("2:30pm", datetime(2026, 7, 21, 14, 30, 0)),
        ("tomorrow 9am", datetime(2026, 7, 22, 9, 0, 0)),
        ("2026-07-22 09:00", datetime(2026, 7, 22, 9, 0, 0)),
    ],
)
def test_parse_when_valid(text: str, expected: datetime) -> None:
    assert _parse_when([text], NOW) == expected


def test_parse_when_time_of_day_rolls_to_tomorrow() -> None:
    """A bare clock time already past today prefers the future (tomorrow)."""
    assert _parse_when(["09:00"], NOW) == datetime(2026, 7, 22, 9, 0, 0)


def test_parse_when_joins_unquoted_tokens() -> None:
    """nargs='+' words are joined, so the user need not quote the expression."""
    assert _parse_when(["in", "2", "hours"], NOW) == datetime(2026, 7, 21, 12, 0, 0)


@pytest.mark.parametrize("tokens", [[], [""], ["banana"], ["not", "a", "time"]])
def test_parse_when_invalid(tokens: list[str]) -> None:
    with pytest.raises(ValueError):
        _parse_when(tokens, NOW)
