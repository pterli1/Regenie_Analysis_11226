#!/usr/bin/env python3
"""Generate the documentation's reference content from fedcat's code.

Run before building or serving the site, from `regenie_workflows/federated_analysis`:

```bash
python tools/docs/gen_docs.py
zensical serve
```

It writes, from the argparse help and the config classes' attribute docstrings:

- `docs/options_reference.md`, every subcommand's options on one page, which is
  committed for reading on GitHub;
- Markdown snippets in `build/docs_snippets/`, which the site's pages include with
  `--8<--`: each subcommand's options, the subcommand table, the state codes, the
  config file reference and the changelog;
- the `fedcat_rs` type stub as a stubs-only package in `build/docs_stubs/`, where
  mkdocstrings finds it before the `fedcat_rs/` source directory.
"""

from __future__ import annotations

import argparse
import enum
import json
import re
import shutil
import sys
from enum import IntEnum
from pathlib import Path

import attrs

PACKAGE_ROOT = Path(__file__).resolve().parents[2]
DOCS_DIR = PACKAGE_ROOT / "docs"
SNIPPETS_DIR = PACKAGE_ROOT / "build" / "docs_snippets"
STUBS_DIR = PACKAGE_ROOT / "build" / "docs_stubs"
CHANGELOG = PACKAGE_ROOT.parent / "CHANGELOG.md"
REPO_URL = "https://github.com/ExeterGenetics/biobank_workflows"

# Import fedcat from this checkout rather than any installed version
if str(PACKAGE_ROOT) not in sys.path:
    sys.path.insert(0, str(PACKAGE_ROOT))


# ── Command-line options ──────────────────────────────────────────────────────


def _format_type(action: argparse.Action) -> str:
    """Return an option's type as shown in the table, e.g. `flag`, `int` or `str[]`."""
    if isinstance(action, argparse._VersionAction):
        return "flag"
    if isinstance(action, (argparse._StoreTrueAction, argparse._StoreFalseAction)):
        return "flag"
    if isinstance(action, argparse._CountAction):
        return "int"
    if action.nargs in ("+", "*"):
        base = getattr(action.type, "__name__", "str") if action.type else "str"
        return f"{base}[]"
    if action.type is not None:
        return getattr(action.type, "__name__", "str")
    return "str"


def _format_default(action: argparse.Action) -> str:
    """Return an option's default as shown in the table, or `""` for none."""
    if isinstance(action, (argparse._StoreTrueAction, argparse._StoreFalseAction)):
        return ""
    if isinstance(action, argparse._VersionAction):
        return ""
    if action.default is None or action.default is argparse.SUPPRESS:
        return ""
    if isinstance(action.default, list):
        if not action.default:
            return "[]"
        return ", ".join(str(x) for x in action.default)
    return str(action.default)


def _format_choices(action: argparse.Action) -> str:
    """Return an option's choices as backticked, comma-separated values, or `""`."""
    if action.choices:
        return ", ".join(f"`{c}`" for c in action.choices)
    return ""


def _render_action(action: argparse.Action) -> dict[str, str]:
    """Return an option's table columns: flags, type, required, default, choices and help."""
    flags = ", ".join(f"`{o}`" for o in action.option_strings) if action.option_strings else f"`{action.dest}`"
    type_str = _format_type(action)
    required = "**yes**" if action.required else ""
    default = _format_default(action)
    choices = _format_choices(action)
    help_text = action.help or ""
    # Expand argparse's %(default)s and %(prog)s
    help_text = help_text.replace("%(default)s", default).replace("%(prog)s", "fedcat")
    # Keep placeholders such as <tab> visible, and table cells intact
    help_text = help_text.replace("<", "&lt;").replace(">", "&gt;").replace("|", "\\|")
    return {
        "flags": flags,
        "type": type_str,
        "required": required,
        "default": default,
        "choices": choices,
        "help": help_text,
    }


def _render_list(rows: list[dict[str, str]]) -> list[str]:
    """Render options as a Markdown definition list, one entry per option."""
    lines: list[str] = []
    for r in rows:
        term = [r["flags"], f"*{r['type']}*"]
        if r["required"]:
            term.append("**required**")
        if r["default"]:
            term.append(f"default `{r['default']}`")
        lines.append(" · ".join(term))
        lines.append(f":   {r['help']}")
        if r["choices"]:
            lines.append("")
            lines.append(f"    Choices: {r['choices']}.")
        lines.append("")
    return lines


def _render_group(
    title: str | None,
    description: str | None,
    actions: list[argparse.Action],
    title_format: str = "**{}**",
    as_list: bool = False,
) -> str:
    """Render an argument group as Markdown, headed by its title and description.

    `title_format` formats the title, e.g. as bold text or a heading. Options are a
    table, or a definition list with `as_list`.
    """
    lines: list[str] = []
    if title and title not in ("positional arguments", "options"):
        lines.append(title_format.format(title))
        if description:
            lines.append("")
            lines.append(f"_{description}_")
        lines.append("")

    rows = []
    for action in actions:
        if isinstance(action, argparse._HelpAction):
            continue
        rows.append(_render_action(action))

    if not rows:
        return ""

    if as_list:
        return "\n".join(lines + _render_list(rows)).rstrip()

    lines.append("| Option | Type | Default | Description |")
    lines.append("|---|---|---|---|")
    for r in rows:
        default_str = r["default"]
        if default_str:
            default_str = f"`{default_str}`"
        help_parts = []
        if r["required"]:
            help_parts.append("**(required)**")
        help_parts.append(r["help"])
        if r["choices"]:
            help_parts.append(f"Choices: {r['choices']}.")
        help_str = " ".join(help_parts)
        lines.append(f"| {r['flags']} | {r['type']} | {default_str} | {help_str} |")

    return "\n".join(lines)


def render_parser(
    parser: argparse.ArgumentParser,
    description: bool = True,
    title_format: str = "**{}**",
    as_list: bool = False,
) -> str:
    """Render a subcommand's argument groups as Markdown.

    Args:
        parser: The subcommand's parser.
        description: Whether to start with the subcommand's description.
        title_format: Format of argument group titles, e.g. `"### {}"` for headings.
        as_list: Render options as a definition list rather than a table.

    Returns:
        The Markdown.
    """
    lines: list[str] = []

    if description and parser.description:
        lines.append(parser.description)
        lines.append("")

    for group in parser._action_groups:
        actions = [
            a
            for a in group._group_actions
            if not isinstance(a, (argparse._HelpAction, argparse._SubParsersAction))
        ]
        if not actions:
            continue
        title = group.title
        rendered = _render_group(title, group.description, actions, title_format, as_list)
        if rendered:
            lines.append(rendered)
            lines.append("")

    return "\n".join(lines)


def _render_enum(cls: type[IntEnum]) -> str:
    """Render an IntEnum as a class definition with one member per line."""
    members = "\n".join(f"    {m.name} = {m.value}" for m in cls)
    return f"class {cls.__name__}(IntEnum):\n{members}"


def render_state_codes() -> str:
    """Render the `JobState` and `StepState` enums as a Python code block."""
    from fedcat.core.structures import JobState, StepState

    return "\n".join(["```python", _render_enum(JobState), "", "", _render_enum(StepState), "```"])


def subcommands() -> dict[str, tuple[str, argparse.ArgumentParser]]:
    """Return each `fedcat` subcommand's one-line help and parser, in CLI order."""
    import fedcat.__main__ as main_mod

    parser = main_mod._build_parser()
    sp_action: argparse._SubParsersAction = next(
        a for a in parser._actions if isinstance(a, argparse._SubParsersAction)
    )
    subcommand_help: dict[str, str] = {a.dest: (a.help or "") for a in sp_action._choices_actions}
    return {name: (subcommand_help.get(name, ""), sub) for name, sub in sp_action.choices.items()}


def render_reference() -> str:
    """Render the whole CLI reference, every subcommand on one page."""
    commands = subcommands()

    out: list[str] = []
    out.append("# Options Reference")
    out.append("")
    out.append("<!-- AUTO-GENERATED by tools/docs/gen_docs.py — do not edit by hand -->")
    out.append("")
    out.append("```")
    out.append("fedcat [--version] <subcommand> [options]")
    out.append("```")
    out.append("")
    out.append("## Subcommands")
    out.append("")
    out.append("| Command | Description |")
    out.append("|---|---|")
    for name, (help_text, _) in commands.items():
        out.append(f"| [`{name}`](#{name}) | {help_text} |")
    out.append("")

    for name, (_, subparser) in commands.items():
        out.append("---")
        out.append("")
        out.append(f"### `{name}`")
        out.append("")
        out.append("```")
        out.append(f"fedcat {name} [options]")
        out.append("```")
        out.append("")
        out.append(render_parser(subparser))

    # State codes, from the enums
    out.append("---")
    out.append("")
    out.append("## State codes")
    out.append("")
    out.append("Integer state codes used by `--reset-jobs` and `--cond-reset-step`:")
    out.append("")
    out.append(render_state_codes())
    out.append("")

    return "\n".join(out)


def render_commands_table() -> str:
    """Render the table of subcommands, linking to their pages."""
    rows = ["| Command | Description |", "|---|---|"]
    rows += [f"| [`{name}`]({name}.md) | {help_text} |" for name, (help_text, _) in subcommands().items()]
    return "\n".join(rows)


def render_subcommand_options(parser: argparse.ArgumentParser, name: str) -> str:
    """Render a subcommand's usage and options for its page, as definition lists."""
    usage = f"```\nfedcat {name} [options]\n```\n\n"
    return usage + render_parser(parser, description=False, title_format="### {}", as_list=True)


# ── Config file reference ─────────────────────────────────────────────────────

# Each config class, the JSON key path it is found at, and a short introduction
CONFIG_SECTIONS: list[tuple[str, str, str]] = [
    ("DiscoveryConfig", "the top level", ""),
    ("CohortCfg", "`cohorts.<biobank>`", "One per biobank to run, keyed `aou` or `ukb`."),
    ("CohortWorkflows", "`cohorts.<biobank>.workflows`", ""),
    (
        "JobResources",
        "`cohorts.<biobank>.workflows.resources.<job>`",
        (
            "Keyed by the workflow field names (`regenie_step1`, ..., `conditioning_aggregate`), "
            "plus `meta_workflow` and `subset_workflow`."
        ),
    ),
    ("AnnotationEntry", "items of `regenie_annotations`", ""),
    ("MetaAnalysisConfig", "`meta_analysis_config`", ""),
    ("CovarStratification", "`stratify_by_covar`", ""),
]

# Fields fedcat sets at runtime, which are not read from the config file
RUNTIME_FIELDS = {("DiscoveryConfig", "covar_stratum"), ("CohortCfg", "output_location")}


def _one_line(text: str) -> str:
    """Join a docstring's lines into one line."""
    return " ".join(line.strip() for line in text.splitlines())


def _json_default(value) -> str:
    """Render a field default as JSON."""
    if isinstance(value, attrs.Factory):
        value = value.factory()
    if isinstance(value, enum.Enum):
        value = value.value
    return json.dumps(value)


def _enum_choices(annotation: str, namespace: dict) -> str:
    """Return the allowed values of an enum named in a type annotation, or `""`."""
    for name in re.findall(r"\b[A-Z]\w+\b", annotation):
        obj = namespace.get(name)
        if isinstance(obj, type) and issubclass(obj, enum.Enum) and name not in ("BioBank", "JobKind"):
            return " One of " + ", ".join(f"`{m.value}`" for m in obj) + "."
    return ""


def render_config_reference() -> str:
    """Render the fields of each config class, for the config file reference page."""
    import griffe

    from fedcat.core import structures

    module = griffe.load("fedcat.core.structures", search_paths=[str(PACKAGE_ROOT)], allow_inspection=False)
    namespace = vars(structures)
    anchors = {name: name.lower() for name, _, _ in CONFIG_SECTIONS}

    out: list[str] = []
    for name, location, intro in CONFIG_SECTIONS:
        cls = getattr(structures, name)
        docs = module.members[name]
        out += [f"## {name} {{#{anchors[name]}}}", "", f"Found at: {location}.", ""]
        if docs.docstring:
            out += [docs.docstring.value.split("\n\n")[0], ""]
        if intro:
            out += [intro, ""]
        for field in attrs.fields(cls):
            if (name, field.name) in RUNTIME_FIELDS:
                continue
            member = docs.members.get(field.name)
            description = _one_line(member.docstring.value) if member is not None and member.docstring else ""
            annotation = str(field.type)
            description += _enum_choices(annotation, namespace)
            for other, anchor in anchors.items():
                if re.search(rf"\b{other}\b", annotation):
                    description += f" See [{other}](#{anchor})."
            default = (
                "**required**" if field.default is attrs.NOTHING else f"default `{_json_default(field.default)}`"
            )
            out += [f"`{field.name}` · *{annotation}* · {default}", f":   {description}", ""]

    return "\n".join(out)


# ── Changelog ─────────────────────────────────────────────────────────────────


_LIST_ITEM = re.compile(r"^( *)[-*] ")


def _python_markdown_lists(text: str) -> str:
    """Rewrite GitHub-style lists for Python-Markdown.

    Python-Markdown needs a blank line before a list and four spaces per nesting level,
    where GitHub accepts a list straight after a paragraph and any deeper indent for
    nesting (2, 4 or 6 spaces in the changelog). Each item's level is worked out from its
    indent relative to the items above it.
    """
    out: list[str] = []
    indents: list[int] = []  # indent of each open nesting level
    for line in text.splitlines():
        item = _LIST_ITEM.match(line)
        if not item:
            if line.strip():
                indents = []
            out.append(line)
            continue
        if out and out[-1].strip() and not _LIST_ITEM.match(out[-1]):
            out.append("")
        indent = len(item[1])
        while indents and indent < indents[-1]:
            indents.pop()
        if not indents or indent > indents[-1]:
            indents.append(indent)
        out.append(" " * (4 * (len(indents) - 1)) + line.lstrip(" "))
    return "\n".join(out)


def render_changelog() -> str:
    """Return the repository changelog for the site, with its links pointing at GitHub."""
    text = re.sub(
        r"\]\((?!https?://|#)(?:\./)?([^)]+)\)",
        rf"]({REPO_URL}/blob/main/regenie_workflows/\1)",
        CHANGELOG.read_text(),
    )
    return _python_markdown_lists(text)


# ── Main ──────────────────────────────────────────────────────────────────────


def _write(path: Path, text: str) -> None:
    """Write a file, creating its directory, with a final newline."""
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(text.rstrip("\n") + "\n")


def main() -> None:
    """Write the options reference, the site's snippets and the staged `fedcat_rs` stub."""
    # print() of the whole reference used to add the final blank line; keep the file identical
    (DOCS_DIR / "options_reference.md").write_text(render_reference() + "\n")

    if SNIPPETS_DIR.exists():
        shutil.rmtree(SNIPPETS_DIR)
    for name, (_, parser) in subcommands().items():
        _write(SNIPPETS_DIR / "options" / f"{name}.md", render_subcommand_options(parser, name))
    _write(SNIPPETS_DIR / "commands.md", render_commands_table())
    _write(SNIPPETS_DIR / "state-codes.md", render_state_codes())
    _write(SNIPPETS_DIR / "config-reference.md", render_config_reference())
    _write(SNIPPETS_DIR / "changelog.md", render_changelog())

    stub = STUBS_DIR / "fedcat_rs" / "__init__.pyi"
    stub.parent.mkdir(parents=True, exist_ok=True)
    shutil.copyfile(PACKAGE_ROOT / "fedcat_rs" / "fedcat_rs.pyi", stub)
    print(
        f"Wrote docs/options_reference.md, {SNIPPETS_DIR.relative_to(PACKAGE_ROOT)}/ and {stub.relative_to(PACKAGE_ROOT)}"
    )


if __name__ == "__main__":
    main()
