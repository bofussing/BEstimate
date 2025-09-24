# Changelog
All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

## [1.1.1] - 2025-09-24
### Added
- Python type hints throughout codebase
- [README.md](README.md) includes:
  - clearer quick start instructions
  - updated CLI usage examples
  - complete CLI argument reference
  - instructions for developer to set up development environment, run tests, work with Docker
  - explanation of how create a release
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
- Secondary CLI entry points `x_genome` and `x_crispranalyser` added
  replacing `python3 x_genome.py ...` and `python3 crispr_analyser.py ...`
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
- Updated base Python Docker image used to build the public and development Docker images - python:3.12-slim bullseye to trixie

### Fixed
- Corrected import paths for `crispr_analyser`, `x_genome`, and
  `x_crispranalyser` to work with package structure.

### Infrastructure
- Complete package modernization from script to proper Python package
- CI uses security configurations including IDS-related tasks

## [1.1.0] - 2025-05-29
_AI generated summary_

### Major Focus: Off-target reporting improvements and efficiency optimizations

### Added
- Timing output for get_off_targets function
- Array of tuples implementation for improved performance

### Changed
- Large refactor of off-target logic replacing Pandas DataFrame with array of tuples for efficiency
- Updated run_offtargets function to work as module
- README updates and formatting improvements

### Fixed
- Off-targets reporting functionality
- Various formatting issues suggested by flake8 and black

### Removed
- Multiprocessing implementation (was slower than serial processing)

## [1.0.0] - 2025-05-08
_AI generated summary_

### Major Focus: CRISPR-Analyser integration and core logic improvements

### Added
- Replacement logic in BEstimate/crispr_analyser
- Paper reproducibility features (Jupyter notebooks and environment)
- Genome retrieval functionality improvements
- Requirements file for dependencies

### Changed
- Updated Conda environment file (bestimate.yml)
- README updates for x_genome.py and general documentation

### Fixed
- Critical bugs in check_index_file() parameters
- Various module import and path issues

### Removed
- Dependency on external CRISPR-Analyser
- Running of x_index.py

## [0.9.0] - 2025-05-05
_AI generated summary_

### Major Focus: Core algorithm development and biological annotations

### Added
- Initial protein-protein interaction (PPI) analysis
- VEP (Variant Effect Predictor) integration
- Uniprot integration for protein change positions
- Post-translational modification (PTM) analysis
- Multiple transcript and variant effect handling
- CDS regions and UTR analysis with MANE transcripts
- Interface disruption annotation
- Ensembl gene ID handling with chromosome control
- Poly-T information support
- GNU GPL v3 Licence

### Changed
- Enhanced protein analysis with percentage calculations
- Improved interaction site extraction
- Result file organization (CRISPR, Edit, Annotation)
