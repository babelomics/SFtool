# Contributing to SFtool

Thank you for your interest in contributing to **SFtool**.
This document describes the contribution workflow, coding standards, and best practices for maintaining a clean and consistent codebase.

If you have any questions, feel free to open a discussion or an issue.

---

## 1. Before You Start

✔ Check whether an issue already exists
✔ Choose the correct issue type
✔ When unsure, open a Discussion instead of an issue

---

## 2. Branching Strategy

**Default branch:** `geneBePharmCAT`

Use branch prefixes:

- feature/
- fix/
- refactor/
- docs/
- perf/
- chore/

Example:

```
git checkout geneBePharmCAT
git pull
git checkout -b feature/<issue-number>-short-description
```

---

## 3. Commit Message Style (Conventional Commits)

Format:

```
<type>: <short description> (#issue)
```

Types: feat, fix, docs, refactor, perf, test, chore

---

## 4. Pull Request Process

Use PR template and reference issues:

```
Closes #87
```

Keep PRs focused. Use Squash & Merge.

---

## 5. Issue Workflow

Use templates (Bug, Feature, Refactor, Documentation, Performance).

Labels assigned automatically.

---

## 6. Code Style

Follow PEP 8, use type hints, ensure tests pass.

---

## 7. Development Environment

```
git clone https://github.com/babelomics/SFtool.git
cd SFtool
conda env create -f environment.yml
conda activate sftool
python SFtool.py --help
```

---

## 8. Documentation Contributions

Use docs/ branches and docs: commits.

---

## 9. Adding Templates

Use a chore branch for metadata changes.

---

## 10. Code of Conduct

Maintain a respectful and constructive environment.

---

## 11. Questions?

Open a Discussion or contact maintainers.

---

Thank you for contributing to SFtool!