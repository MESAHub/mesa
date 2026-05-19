#!/usr/bin/env python3
"""Summarize documentation and tagging coverage for the MESA test suite."""

from __future__ import annotations

import argparse
import json
import re
from dataclasses import dataclass
from pathlib import Path


MODULES = ("star", "binary", "astero")
DO_ONE_RE = re.compile(r"^\s*do_one\s+([^\s#]+)")


@dataclass(frozen=True)
class TestCase:
    module: str
    name: str

    @property
    def label(self) -> str:
        return f"{self.module}/{self.name}"


def repo_root() -> Path:
    return Path(__file__).resolve().parents[1]


def sort_key(case: TestCase) -> tuple[str, str]:
    return (case.module, case.name.lower())


def read_active_tests(root: Path) -> list[TestCase]:
    tests: list[TestCase] = []
    for module in MODULES:
        source = root / module / "test_suite" / "do1_test_source"
        for line in source.read_text(encoding="utf-8").splitlines():
            match = DO_ONE_RE.match(line)
            if match:
                tests.append(TestCase(module, match.group(1)))
    return sorted(tests, key=sort_key)


def read_test_dirs(root: Path) -> list[TestCase]:
    tests: list[TestCase] = []
    for module in MODULES:
        suite = root / module / "test_suite"
        for child in suite.iterdir():
            if child.is_dir():
                tests.append(TestCase(module, child.name))
    return sorted(tests, key=sort_key)


def read_docs_pages(root: Path) -> dict[str, Path]:
    docs_dir = root / "docs" / "source" / "test_suite"
    return {path.stem: path for path in sorted(docs_dir.glob("*.rst"))}


def read_tags(readme: Path) -> list[str]:
    if not readme.exists():
        return []

    tags: list[str] = []
    lines = readme.read_text(encoding="utf-8").splitlines()
    for index, line in enumerate(lines):
        stripped = line.strip()
        if not stripped.startswith(".. tags::"):
            continue

        inline_tags = stripped.split("::", 1)[1]
        tags.extend(split_tags(inline_tags))

        for continuation in lines[index + 1 :]:
            if not continuation.strip():
                break
            if continuation.startswith((" ", "\t")):
                tags.extend(split_tags(continuation.strip()))
                continue
            break
        break

    return sorted(set(tags), key=str.lower)


def split_tags(text: str) -> list[str]:
    return [tag.strip() for tag in text.split(",") if tag.strip()]


def gather_inventory(root: Path) -> dict[str, object]:
    active_tests = read_active_tests(root)
    test_dirs = read_test_dirs(root)
    docs_pages = read_docs_pages(root)

    active_names = {case.name for case in active_tests}
    dir_names = {case.name for case in test_dirs}
    doc_names = set(docs_pages)

    readmes = {
        case: root / case.module / "test_suite" / case.name / "README.rst"
        for case in test_dirs
    }
    tags_by_case: dict[str, list[str]] = {}
    for case, readme in readmes.items():
        if readme.exists():
            tags = read_tags(readme)
            if tags:
                tags_by_case[case.label] = tags

    active_missing_dirs = [case for case in active_tests if case.name not in dir_names]
    active_missing_readme = [
        case
        for case in active_tests
        if not (root / case.module / "test_suite" / case.name / "README.rst").exists()
    ]
    active_missing_docs = [case for case in active_tests if case.name not in doc_names]
    active_missing_tags = [
        case
        for case in active_tests
        if case.label not in tags_by_case
    ]
    inactive_dirs = [case for case in test_dirs if case.name not in active_names]
    documented_inactive = sorted(doc_names - active_names, key=str.lower)
    readmes_without_docs = [
        case for case in test_dirs if readmes[case].exists() and case.name not in doc_names
    ]
    broken_docs = sorted(
        [path.name for path in docs_pages.values() if path.is_symlink() and not path.exists()],
        key=str.lower,
    )

    module_counts = {
        module: {
            "active": sum(1 for case in active_tests if case.module == module),
            "directories": sum(1 for case in test_dirs if case.module == module),
        }
        for module in MODULES
    }

    return {
        "root": str(root),
        "module_counts": module_counts,
        "active_tests": [case.label for case in active_tests],
        "test_directories": [case.label for case in test_dirs],
        "docs_pages": sorted(doc_names, key=str.lower),
        "tagged_readmes": tags_by_case,
        "active_missing_dirs": [case.label for case in active_missing_dirs],
        "active_missing_readme": [case.label for case in active_missing_readme],
        "active_missing_docs": [case.label for case in active_missing_docs],
        "active_missing_tags": [case.label for case in active_missing_tags],
        "inactive_dirs": [case.label for case in inactive_dirs],
        "documented_inactive": documented_inactive,
        "readmes_without_docs": [case.label for case in readmes_without_docs],
        "broken_docs": broken_docs,
    }


def limited(items: list[str], limit: int) -> tuple[list[str], int]:
    if limit <= 0 or len(items) <= limit:
        return items, 0
    return items[:limit], len(items) - limit


def bullet_section(title: str, items: list[str], limit: int) -> list[str]:
    lines = [f"## {title}"]
    if not items:
        lines.append("- None")
        return lines

    shown, remaining = limited(items, limit)
    lines.extend(f"- {item}" for item in shown)
    if remaining:
        lines.append(f"- ... {remaining} more (rerun with --limit 0 to show all)")
    return lines


def render_markdown(inventory: dict[str, object], limit: int) -> str:
    module_counts = inventory["module_counts"]
    assert isinstance(module_counts, dict)

    active_tests = inventory["active_tests"]
    test_dirs = inventory["test_directories"]
    docs_pages = inventory["docs_pages"]
    tagged_readmes = inventory["tagged_readmes"]
    assert isinstance(active_tests, list)
    assert isinstance(test_dirs, list)
    assert isinstance(docs_pages, list)
    assert isinstance(tagged_readmes, dict)

    lines = [
        "# MESA test-suite inventory",
        "",
        f"Repository: `{inventory['root']}`",
        "",
        "## Summary",
        f"- Active `do_one` tests: {len(active_tests)}",
        f"- Test directories: {len(test_dirs)}",
        f"- Test-suite docs pages: {len(docs_pages)}",
        f"- README files with tags: {len(tagged_readmes)}",
    ]
    for module in MODULES:
        counts = module_counts[module]
        lines.append(
            f"- {module}: {counts['active']} active, {counts['directories']} directories"
        )

    sections = [
        ("Active tests missing directories", inventory["active_missing_dirs"]),
        ("Active tests missing README.rst", inventory["active_missing_readme"]),
        ("Active tests missing docs pages", inventory["active_missing_docs"]),
        ("Active tests without tags", inventory["active_missing_tags"]),
        ("Inactive test directories", inventory["inactive_dirs"]),
        ("README files without docs pages", inventory["readmes_without_docs"]),
        ("Docs pages not listed in do1_test_source", inventory["documented_inactive"]),
        ("Broken docs symlinks", inventory["broken_docs"]),
    ]
    for title, items in sections:
        assert isinstance(items, list)
        lines.extend(["", *bullet_section(title, items, limit)])

    return "\n".join(lines) + "\n"


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Summarize MESA test-suite docs and tag coverage."
    )
    parser.add_argument(
        "--format",
        choices=("markdown", "json"),
        default="markdown",
        help="Output format.",
    )
    parser.add_argument(
        "--limit",
        type=int,
        default=40,
        help="Maximum items per Markdown section; use 0 for no limit.",
    )
    args = parser.parse_args()

    inventory = gather_inventory(repo_root())
    if args.format == "json":
        print(json.dumps(inventory, indent=2, sort_keys=True))
    else:
        print(render_markdown(inventory, args.limit), end="")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
