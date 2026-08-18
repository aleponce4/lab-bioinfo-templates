#!/usr/bin/env python3
"""Check that every third-party module the Python templates import is declared in
environment.yml.

This exists because `matplotlib-venn` was imported by template 11, installed in the
render workflow, and missing from environment.yml — so anyone who followed the README's
conda instructions silently lost a figure. Standard-library imports are ignored.

Run from the repository root. Exits non-zero on an undeclared import.
"""

from __future__ import annotations

import ast
import json
import pathlib
import re
import sys

# Import name -> distribution name, where they differ.
DISTRIBUTION_NAMES = {
    "matplotlib_venn": "matplotlib-venn",
    "sklearn": "scikit-learn",
    "yaml": "pyyaml",
    "PIL": "pillow",
    "cv2": "opencv-python",
    "Bio": "biopython",
    "umap": "umap-learn",
}

# Provided by the Quarto jupyter engine rather than imported by template code.
ALWAYS_ALLOWED = {"ipykernel", "nbformat", "nbconvert", "jupyter", "jupyterlab"}

ROOT = pathlib.Path(".")


def top_level_imports(source: str) -> set[str]:
    try:
        tree = ast.parse(source)
    except SyntaxError:
        return set()

    found: set[str] = set()
    for node in ast.walk(tree):
        if isinstance(node, ast.Import):
            found.update(alias.name.split(".")[0] for alias in node.names)
        elif isinstance(node, ast.ImportFrom):
            # `from . import x` has no module; skip relative imports.
            if node.level == 0 and node.module:
                found.add(node.module.split(".")[0])
    return found


def python_chunks(qmd_text: str) -> str:
    """Concatenate the ```{python} chunks of a .qmd."""
    return "\n".join(
        m.group(1) for m in re.finditer(r"^```+\s*\{python[^}]*\}\s*$(.*?)^```+\s*$",
                                        qmd_text, re.MULTILINE | re.DOTALL)
    )


def notebook_source(nb_text: str) -> str:
    try:
        nb = json.loads(nb_text)
    except json.JSONDecodeError:
        return ""
    return "\n".join(
        "".join(cell.get("source", ""))
        for cell in nb.get("cells", [])
        if cell.get("cell_type") == "code"
    )


def declared_packages(env_path: pathlib.Path) -> set[str]:
    """Parse environment.yml without requiring PyYAML (CI has no pip installs)."""
    names: set[str] = set()
    in_dependencies = False
    for line in env_path.read_text(encoding="utf-8").splitlines():
        # Track the top-level section so `channels:` entries aren't mistaken for packages.
        if line and not line[0].isspace():
            in_dependencies = line.split(":", 1)[0].strip() == "dependencies"
            continue
        if not in_dependencies:
            continue

        stripped = line.strip()
        if not stripped.startswith("- ") or stripped.endswith(":"):
            continue
        # "- pandas>=2.0" -> "pandas";  "- python=3.11" -> "python"
        spec = stripped[2:].strip()
        name = re.split(r"[=<>!~\[ ]", spec, maxsplit=1)[0]
        if name:
            names.add(name.lower().replace("_", "-"))
    return names


def pip_requirements(env_path: pathlib.Path) -> set[str]:
    """Return the entries under the `- pip:` block, as pip requirement strings."""
    reqs: set[str] = set()
    in_pip = False
    pip_indent = 0
    for line in env_path.read_text(encoding="utf-8").splitlines():
        if not line.strip():
            continue
        indent = len(line) - len(line.lstrip())
        if line.strip() in ("- pip:", "-pip:"):
            in_pip = True
            pip_indent = indent
            continue
        if in_pip:
            if indent <= pip_indent:
                in_pip = False
                continue
            if line.strip().startswith("- "):
                reqs.add(line.strip()[2:].strip())
    return reqs


def main() -> int:
    env_path = ROOT / "environment.yml"
    if not env_path.exists():
        print("environment.yml not found; run from the repository root.", file=sys.stderr)
        return 2

    declared = declared_packages(env_path)

    # `--pip-list` prints the pip requirements so CI can install exactly what
    # environment.yml declares, instead of maintaining a second hardcoded list.
    if "--pip-list" in sys.argv:
        print(" ".join(sorted(pip_requirements(env_path))))
        return 0

    sources: dict[str, set[str]] = {}
    for path in sorted(ROOT.glob("templates/*/template.qmd")):
        sources[str(path)] = top_level_imports(python_chunks(path.read_text(encoding="utf-8")))
    for path in sorted(ROOT.glob("templates/*/template.ipynb")):
        sources[str(path)] = top_level_imports(notebook_source(path.read_text(encoding="utf-8")))
    for path in sorted(ROOT.glob("templates/*/data/simulate_data.py")):
        sources[str(path)] = top_level_imports(path.read_text(encoding="utf-8"))

    stdlib = sys.stdlib_module_names
    missing: dict[str, set[str]] = {}
    used_distributions: set[str] = set()

    for path, modules in sources.items():
        for module in sorted(modules):
            if module in stdlib or module in ALWAYS_ALLOWED:
                continue
            dist = DISTRIBUTION_NAMES.get(module, module).lower().replace("_", "-")
            used_distributions.add(dist)
            if dist not in declared:
                missing.setdefault(dist, set()).add(path)

    for dist, paths in sorted(missing.items()):
        for path in sorted(paths):
            print(f"MISSING: {dist!r} is imported by {path} but not declared in environment.yml")

    unused = declared - used_distributions - ALWAYS_ALLOWED - {"python", "pip"}
    for dist in sorted(unused):
        print(f"note: {dist!r} is declared in environment.yml but never imported")

    if missing:
        print(f"\n{len(missing)} undeclared dependency/dependencies.", file=sys.stderr)
        return 1

    print(f"All {len(used_distributions)} third-party import(s) are declared in environment.yml.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
