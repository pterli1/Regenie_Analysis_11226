# Docstring conventions

fedcat's API reference is to be generated from its docstrings with MkDocs Material and mkdocstrings, so docstrings are Google style with Markdown markup.

## Python

- Every module, class and function has a docstring describing the current code. Leave out how it came to be, unless that explains something surprising.
- The summary is one line ending in a full stop: imperative for functions (`Return ...`, `Build ...`), a noun phrase for classes, properties and attributes.
- Public functions and methods (no leading `_`) have `Args:`, `Returns:` and `Raises:` sections as they apply. Leave out `self`/`cls`, and `Returns:` for `None`.
- Private helpers get a summary that names any argument whose meaning is not obvious, with sections only where they help.
- A method implementing an interface (such as `Platform`) documents only what is specific to it; the arguments are documented on the interface.
- The `check_*` functions in `fedcat.core.checks` share one signature, documented in the module, so they have a summary only.
- Class fields, including config fields, are documented by a string on the line after the field, not an `Attributes:` section:

  ```python
  cond_max_steps: int = -1
  """Most conditioning steps; -1 for no limit."""
  ```

- Markup is Markdown: single backticks for code, and cross-references as ``[`DiscoveryRun`][fedcat.core.engine.DiscoveryRun]``. Use cross-references sparingly.
- Comments explain why, or what is not evident from the code; don't restate the code.
- Tests have a module docstring and a one-line docstring on helpers, fixtures and fake classes. A test function needs one only when its name does not say what it checks.

## Rust

- Each file starts with `//!` module docs, and every function, struct and field has `///` docs.
- Python-facing functions use `# Arguments`, `# Returns`, `# Errors` and `# Panics` sections as they apply.
- The Python API of `fedcat_rs` is documented in `fedcat_rs/fedcat_rs.pyi`, in the Python style above. Keep it in step with the `#[pyfunction]` signatures.

## Command line

Options are documented in their argparse `help`. `docs/options_reference.md` and the options on the documentation site are generated from them:

```bash
python tools/docs/gen_docs.py
```
