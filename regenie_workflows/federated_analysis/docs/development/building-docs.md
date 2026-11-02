# Building these docs

The site is built with [Zensical](https://zensical.org/), from the Material for MkDocs team, with [mkdocstrings](https://mkdocstrings.github.io/) for the API reference. The sources are `mkdocs.yml` and `docs/` in `regenie_workflows/federated_analysis`.

## Building locally

Install FedCAT (with `fedcat_rs`), then the `docs` extra. From `regenie_workflows/federated_analysis`:

```bash
pip install -e ".[docs]"
python tools/docs/gen_docs.py   # generate the reference content
zensical serve                  # live preview at http://localhost:8000
zensical build                  # static site in site/
```

Rerun `gen_docs.py` after changing a command-line option, a config class or the changelog. It imports `fedcat`, so needs its dependencies, including `fedcat_rs`.

The configuration also builds with MkDocs 1.x (`pip install "mkdocs<2" mkdocs-material mkdocs-api-autonav mkdocs-exclude`, then `mkdocs serve`).

## Publishing

The site is published at <https://exetergenetics.github.io/fedcat/> by the `fedcat-docs` GitHub Actions workflow (`.github/workflows/fedcat-docs.yml`). On every push to `main` that touches fedcat or its docs, it builds the site with `zensical build --strict` and pushes it to the `gh-pages` branch of [`ExeterGenetics/fedcat`](https://github.com/ExeterGenetics/fedcat), which serves it. Pull requests build the site without publishing it, to catch broken pages and an out-of-date `docs/options_reference.md`.

The workflow authenticates to `ExeterGenetics/fedcat` with a deploy key, stored in this repository as the `FEDCAT_DOCS_DEPLOY_KEY` secret.

## What is generated

`tools/docs/gen_docs.py` writes Markdown snippets to `build/docs_snippets/` (not committed), which pages include with `--8<-- "<snippet>"`:

| Snippet | Included by | Generated from |
|---|---|---|
| `options/<subcommand>.md` | Each [command page](../commands/index.md) | The subcommand's argparse help |
| `commands.md`, `state-codes.md` | [Commands](../commands/index.md), [discover](../commands/discover.md) | The subcommands, and the `JobState`/`StepState` enums |
| `config-reference.md` | [Config file reference](../reference/config.md) | The attribute docstrings of the config classes in `fedcat.core.structures` |
| `changelog.md` | [Changelog](../reference/changelog.md) | `regenie_workflows/CHANGELOG.md`, with its links pointed at GitHub |

It also stages the `fedcat_rs` type stub (`fedcat_rs/fedcat_rs.pyi`) in `build/docs_stubs/` for the [fedcat_rs page](../reference/api/fedcat_rs.md), and rewrites [`docs/options_reference.md`](https://github.com/ExeterGenetics/biobank_workflows/blob/main/regenie_workflows/federated_analysis/docs/options_reference.md), all options on one page for GitHub. That file is committed; `fedcat/tests/test_options_reference.py` fails when it is out of date.

The `api-autonav` plugin makes an API page per `fedcat` module, so new modules need no docs change. Everything else in `docs/` is written by hand.

## Writing docstrings

See [Docstring conventions](../docstring_conventions.md): Google style with Markdown, and cross-references like ``[`DiscoveryRun`][fedcat.core.engine.DiscoveryRun]``. The build warns about unresolved cross-references.
