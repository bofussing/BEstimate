# Changelog
All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Added
- Python type hints throughout codebase
- Comprehensive Sphinx-style docstrings for all public classes and functions
- Initial integration tests for CLI commands and end-to-end testing framework
- Multi-stage Docker builds with production and development containers
- Complete GitLab CI/CD pipeline with stages to build Docker images, run tests,
  publish images to GitLab Container Registry, and deploy to Python packages to
  the GitLab Package Registry
- End-to-end test introduced with `tests/e2e_test.sh`
- Initial `pre-commit` hooks with Black formatting and code quality checks
    - Many flake8 warnings are commented out for now, to be addressed in future
      refactor work
- Checksum validation for data files to ensure integrity

### Changed
- Backwards compatibility breaking changes:
    - Minimal support for Python 3.8 dropped. Now requires Python 3.12 or higher.
- Main CLI entry point refactored to work with package structure instead of as a script
    - Specifically when it is now possible to run `python3 -m BEstimate ...`
      instead of `cd BEstimate && python3 BEstimate.py ...`
- Converted Python package to use `pyproject.toml` in compliance with
  [python.org](https://packaging.python.org/en/latest/specifications/pyproject-toml/)
  recommendations.
    - The `pyproject.toml` was defined with `poetry 2.*` to bring [PEP 621](https://peps.python.org/pep-0621/)
      metadata support and [compliance with modern Python packaging standards](https://python-poetry.org/blog/announcing-poetry-2.0.0/).
- **Ongoing refactor work**. As part of package streamlining small refactors were done to
  improve code quality and maintainability, with larger dependent on improved test code coverage.
    - Converted all tabs to 4 spaces throughout codebase
    - Global variables in `Bestimate.py` e.g. `path`, `ot_path` and `args` are
      no longer set in the module scope but within the `main()` function. These
      globals are set once and renamed to `OUTPUT_PATH`, `OT_PATH` and `ARGS`
      respectively to indicate their constant nature.
    - Data files like `H_sapiens_interfaces.txt` are included within the package structure and referenced via
      `importlib.resources` to ensure compatibility.
    - Refactored path handling to use `OUTPUT_PATH` global consistently

### Fixed
- Corrected import paths for `crispr_analyser`, `x_genome`, and
  `x_crispranalyser` to work with package structure.

### Infrastructure
- Complete package modernization from script to proper Python package
- CI uses security configurations including IDS-related tasks

## [1.1.0] - 2025-05-29

## [1.0.0] - 2025-05-08

## [0.9.0] - 2025-05-05
