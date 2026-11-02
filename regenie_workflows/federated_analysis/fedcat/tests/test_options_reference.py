"""The committed CLI reference, `docs/options_reference.md`, is up to date with the argparse help."""

import importlib.util
from pathlib import Path

PACKAGE_ROOT = Path(__file__).resolve().parents[2]


def test_options_reference_up_to_date():
    spec = importlib.util.spec_from_file_location("gen_docs", PACKAGE_ROOT / "tools" / "docs" / "gen_docs.py")
    gen_docs = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(gen_docs)

    committed = (PACKAGE_ROOT / "docs" / "options_reference.md").read_text()
    assert committed == gen_docs.render_reference() + "\n", (
        "docs/options_reference.md is out of date; regenerate it with `python tools/docs/gen_docs.py`"
    )
