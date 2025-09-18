# BEstimate

BEstimate, a Python package that systematically identifies guide RNA (gRNA)
targetable sites across given sequences for given Base Editors, functional and
clinical effects of the potential edits on the resulting proteins and off target
consequence of the found sequences.

It has the ability to provide in silico analysis of the sequences to identify
positions that can be editable by Base Editors, and their features before
starting experiments.

|                          Main                          |                         Develop                          |
| :----------------------------------------------------: | :------------------------------------------------------: |
| [![pipeline status][main-pipe-badge]][main-branch] | [![pipeline status][develop-pipe-badge]][develop-branch] |

[main-pipe-badge]: https://gitlab.internal.sanger.ac.uk/sci/BEstimate/badges/main/pipeline.svg
[main-branch]: https://gitlab.internal.sanger.ac.uk/sci/BEstimate/-/commits/main
[develop-pipe-badge]: https://gitlab.internal.sanger.ac.uk/sci/BEstimate/badges/develop/pipeline.svg
[develop-branch]: https://gitlab.internal.sanger.ac.uk/sci/BEstimate/-/commits/develop


## Table of Contents
- [Quick start installation](#quick-start-installation)
- [Run BEstimate](#run-bestimate)
    - [Examples with BEstimate](#examples-with-bestimate)
    - [Off-Target examples](#off-targets-examples)
    - [Command line usage and options](#command-line-usage-and-options)
- [Contact](#contact)
- [License](#license)
- [Development](#development)
    - [Software Requirements](#software-requirements)
    - [One time setup](#one-time-setup)
    - [Adding new dependencies](#adding-new-dependencies)
    - [Formatting and pre-commit hooks](#formatting-and-pre-commit-hooks)
    - [Testing](#testing)
    - [CICD (Gitlab CI)](#cicd-gitlab-ci)
- [Git and Tagged releases](#git-and-tagged-releases)

## Quick start installation

To install BEstimate, you require Python 3.12 or higher. It is recommended to use a virtual environment.

```bash
python3 -m venv .venv
source .venv/bin/activate
pip install git+https://gitlab.internal.sanger.ac.uk/sci/BEstimate.git
```

## Run BEstimate

### Examples with BEstimate

For example, if you would like to run for the *SRY* gene with NGG PAM sequence, with CBE (C to T editing) and without VEP and protein analysis:

```bash
BEstimate -gene SRY -assembly GRCh38 -pamseq NGG -pamwin 21-23 -actwin 4-8 -protolen 20 -edit C -edit_to T -o ../output/ -ofile SRY_CBE_NGG
```

The user also run the same analysis for different PAM only changing -pamseq NGN.

*Warning: Be careful to write the PAM sequence to be in concordance with the length of the -pamwin. Here, NGN is in concordance with 21-23 (3 nucleotides). Otherwise, the user need to write NG -pamseq with 21-22 -pamwin.*

If you would like to run for a specific transcript and run the protein analysis:

```bash
BEstimate -gene SRY -assembly GRCh38 -transcript ENST00000383070 -edit C -edit_to T -vep -o ../output/ -ofile SRY_CBE_NGG
```

If you would like to run with a specific point mutation, with NGN PAM and with VEP and protein analysis:
Prepare a `PIK3CA_mutation_file.txt` for example with 3:g.179218303G>A

```bash
BEstimate -gene PIK3CA -assembly GRCh38 -pamseq NGN -pamwin 21-23 -actwin 4-8 -protolen 20 -mutation_file PIK3CA_mutation_file.txt -edit A -edit_to G -vep -ofile PIK3CA_NGN_ABE_mE545K -o ../output/
```

### Off-Target examples

To run the off-target analysis, first you need to have the [Ensembl](https://www.ensembl.org/) Genome indexed for the interested PAM sequence.

The `x_genome` program will download the required files and index the genome for CRISPRs as follows.
- Download the specified FASTA genome assembly files from the Ensembl project,
- Gather CRISPRs from the FASTA files into CSV files detailing chromosome, position in chromosome, as well as PAM position,
- Generate a binary list of gRNA signatures (accounting for PAM position),
- Insert the CRISPRs into a SQLite database for cross-referencing the gRNAs found in the binary list.

to run the **x_genome** program see [command line options](#command-line-options) section below.

For example:

```bash
x_genome --pamseq NGN --assembly GRCh38 --ensembl_version 113
```

The gathering of CRISPRs from the genome assembly takes a while and requires a fair amount of disk storage. For example, using the GRCh38 genome assembly:

| Pam Sequence | Space (GB) | Run Time  |
| ------------ | ---------- | --------- |
| NGG          | 38         | ~3 Hours  |
| NGN          | 140        | ~9 Hours  |

Then, you can run the off-target analysis, see below for the *BRAF* gene:

```bash
BEstimate -gene BRAF -assembly GRCh38 -pamseq NGN -edit A -edit_to G -vep -ot -o ../output -ot_path ../offtargets -ofile BRAF_ABE_NGN
```

### Command line usage and options

There are three programs when the package is installed available from the command line:
- `BEstimate` - the main program to find and analyse Base Editor sites
- `x_genome` - the program to download and index a genome for off-target analysis
- `x_crispranalyser` - the program to run off-target analysis on guides

<details>
<summary>Expand to see <strong>BEstimate</strong> command line options</summary>

```bash
BEstimate --help
usage: BEstimate [inputs]

********************************** Find and Analyse Base Editor sites **********************************

Mandatory Inputs:
  -h, --help            show this help message and exit
  --version             Show program's version number and exit.
  -gene GENE            The hugo symbol of the interested gene!
  -assembly ASSEMBLY    The genome assembly that will be used!
  -transcript TRANSCRIPT
                        The interested ensembl transcript id
  -uniprot UNIPROT      The interested Uniprot id
  -pamseq PAMSEQ        The PAM sequence in which features used for searching activity window and editable nucleotide.
  -pamwin PAMWINDOW     The index of the PAM sequence when starting from the first index of protospacer as 1.
  -actwin ACTWINDOW     The index of the activity window when starting from the first index of protospacer as 1.
  -protolen PROTOLEN    The total protospacer and PAM length.
  -vep                  The boolean option if user wants to analyse the edits through VEP and Uniprot.
  -mutation_file MUTATION_FILE
                        If you have more than one mutations, a file for the mutations on the interested gene that you need to integrate into guide and/or annotation analysis
  -flank                The boolean option if the user wants to add flanking sequences of the gRNAs
  -flank3 FLAN_3        The number of nucleotides in the 3' flanking region
  -flank5 FLAN_5        The number of nucleotides in the 5' flanking region
  -edit {A,T,G,C}       The nucleotide which will be edited.
  -edit_to {A,T,G,C}    The nucleotide after edition.
  -o OUTPUT_PATH        The path for output. If not specified the current directory will be used!
  -ofile OUTPUT_PATH    The output file name, if not specified "position" will be used!
  -ot                   Whether off targets will be computed or not
  -genome GENOME        (If -ot provided) name of the genome file
  -v_ensembl VERSION    The ensembl version in which genome will be retrieved (if the assembly is GRCh37 then please use <=75)
  -ot_path OT_PATH
```

</details>

<details>
<summary>Expand to see <strong>x_genome</strong> command line options</summary>

```bash
usage: x_genome [inputs]

Script for indexing CRISPRs for finding off-targets

options:
  -h, --help            show this help message and exit
  --version             Show the version number and exit
  --pamseq PAMSEQ, -p PAMSEQ
                        The PAM sequence in which features used for searching activity window and editable nucleotide.
  --assembly {GRCh38,GRCh37}, -a {GRCh38,GRCh37}
                        The genome assembly that will be used!
  --output_path OUTPUT_PATH, -o OUTPUT_PATH
                        The path for output. If not specified the current directory will be used!
  --ensembl_version ENSEMBL_VERSION, -e ENSEMBL_VERSION
                        The ensembl version in which genome will be retrieved (if the assembly is GRCh37 then please use <=75)
  --offtargets_path OFFTARGETS_PATH, -ot OFFTARGETS_PATH
                        The path to the root offtargets output directory
```

</details>

<details>
<summary>Expand to see <strong>x_crispranalyser</strong> command line options</summary>

```bash
usage: x_crispranalyser [inputs]

Script for finding off-targets

options:
  -h, --help            show this help message and exit
  --version             Show the version number and exit
  --input_csv INPUT_CSV, -i INPUT_CSV
                        The input CSV file to be analysed
  --binary_index BINARY_INDEX, -b BINARY_INDEX
                        The CRISPR binary index file generated by x_genome.py
  --output_csv OUTPUT_CSV, -o OUTPUT_CSV
                        The output CSV generated
  --db_file DB_FILE, -d DB_FILE
                        The CRISPR DB file generated by x_genome.py
```

</details>

## Contact

BEstimate is the product of Cansu Dinçer, Matthew Coelho and Mathew Garnett from Garnett Group at the Wellcome Sanger Institute. Off-target analysis has been adapted by Bo Fussing from the Cellular Informatics team within the Wellcome Sanger Institute.

If you have any problems or feedback regarding BEstimate, please contact [here](mailto:cd7@sanger.ac.uk).

## License

GNU AFFERO GENERAL PUBLIC LICENSE

BEstimate: A Python module to design and annotate base editor gRNAs

Copyright (C) 2025 Genome Research Limited

Authors: Cansu Dinçer (cd7@sanger.ac.uk), Bo Fussing (bf15@sanger.ac.uk)

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU Affero General Public License as published
by the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU Affero General Public License for more details.

## Further Disclaimer
This tool is for research purposes and not for clinical use.
For policies regarding the underlying data, please also refer to:
- [Ensembl terms and conditions](https://www.ensembl.org/info/about/legal/code_licence.html#:~:text=Subject%20to%20the%20terms%20and,the%20Work%20and%20such%20Derivative)
- [Uniprot terms and conditions](https://www.uniprot.org/help/license#:~:text=We%20make%20no%20warranties%20regarding,by%20patents%20or%20other%20rights.)
- [Interactome Insider terms and conditions](http://interactomeinsider.yulab.org/)


## Development

### Software Requirements
Development requires are quite minimal with two approachs - with a *Virtual Environment* or with a *VSCode devcontainer*.

<details>
<summary><strong>Virtual Environment Requirements</strong></summary>

- **Git Hubflow** tools:
    - On linux see [hubflow installation instructions](https://github.com/dockstore/hubflow)
    - On MacOS with Homebrew: `brew install hubflow`
- **Python 3.12 or higher**. It is up to the developer how to install a specific
  Python if the system default is not suitable.
    - It's recommended to use [pyenv](https://github.com/pyenv/pyenv/blob/master/README.md#installation) as it offers flexibility in managing multiple Python versions. Works on Linux and MacOS, and installing via the command line is straightforward:
        ```bash
        curl https://pyenv.run | bash
        pyenv install 3.12.4
        cd <project-repository-directory>
        pyenv local 3.12.4
        ```
    - Via poetry itself. Since poetry 2.1.0 running `poetry python install 3.12` will [install a standalone Python](https://python-poetry.org/docs/cli/#python-install)
    - On Ubuntu Linux via `apt-get` typically the default python3 is old so we need to add a PPA for newer versions. Famously is the [Deadsnakes PPA](https://launchpad.net/~deadsnakes/+archive/ubuntu/ppa):
        ```bash
        sudo apt update
        sudo apt install -y software-properties-common

        # Add the Deadsnakes PPA and refresh package lists
        sudo add-apt-repository -y ppa:deadsnakes/ppa
        sudo apt update

        # Install Python 3.12, its dev headers, and venv support
        sudo apt install -y python3.12 python3.12-dev python3.12-venv python3.12-distutils

        # Bootstrap pip for this interpreter, then upgrade basics
        python3.12 --version  # sanity check
        python3.12 -m ensurepip --upgrade
        python3.12 -m pip install --upgrade pip setuptools wheel
        ```
    - On MacOS via [Homebrew](https://brew.sh/)
        ```bash
        brew update
        brew install python@3.12
        python3.12 --version # sanity check
        ```
- **Poetry 2.2 or higher**. Poetry is used for dependency management and packaging, there are [multiple ways to install it](https://python-poetry.org/docs/#installation). The official docs recommend will have up-to-date instructions, but in summary:
    - On Linux one typically `sudo apt-get installs pipx` and then `pipx install poetry` to ensure Poetry is up to date and isolated from system Python packages. *Do not* directly `apt-get install poetry` as that version is out of date and not compatible with this project.
    - On MacOS, if you have Homebrew installed, you can use it to install Poetry `brew update && brew install poetry`.

</details>

<details>
<summary><strong>VSCode devcontainer Requirements</strong></summary>

Working with a devcontainer is the easiest way to get started with development
as all dependencies and tools are pre-installed. Working with VSCode is optional
as other IDEs/editors can attach to a running container or you can run commands
directly in the container with `docker run`.

You will need:
- **Docker**
- (Optional) **VSCode** with the [Remote - Containers extension](https://marketplace.visualstudio.com/items?itemName=ms-vscode-remote.remote-containers)

Up-to-date instructions for installing Docker and VSCode can be found on their
respective websites, but the [VSCode instructions summarise all the main steps](https://code.visualstudio.com/docs/devcontainers/containers).

An important note about SSH and working with GitLab/GitHub from within a
devcontainer: forward your keys via your ssh-agent:

`~/.ssh/config` on your Laptop/Workstation (not OpenStack/Farm or other remote host):
```bash
Host *
  AddKeysToAgent yes
  UseKeychain yes
  IdentityFile ~/.ssh/id_rsa
  TCPKeepAlive yes
  ServerAliveInterval 120
```

</details>

### One time setup

<details>
<summary><strong>Virtual Environment Setup</strong></summary>

The process to setup the development environment is as follows:

```bash
git clone https://gitlab.internal.sanger.ac.uk/sci/BEstimate.git
cd BEstimate
python3.12 -m venv .venv
source .venv/bin/activate
poetry install
pre-commit install --install-hooks
```

</details>

<details>
<summary><strong>DevContainer Setup</strong></summary>

The process to setup the development environment is as follows:

```bash
git clone https://gitlab.internal.sanger.ac.uk/sci/BEstimate.git
cd BEstimate
docker build -t bestimate-dev:local -f .devcontainer/Dockerfile .
docker run -it --rm -v $(pwd):/opt/repo  bestimate-dev:local bash

# Inside the container
pwd # should be /opt/repo
BEstimate --version # should show version
exit # exit the container
```

Next, if using VSCode, open the project folder within VSCode and
open the command palette (`Shift+CMD+P` or `Shift+CTRL+P`) and select
`"Remote-Containers: Rebuild and Reopen in Container"`
</details>

### Adding new dependencies

To add a new dependency, use Poetry to add it to the `pyproject.toml` and
`poetry.lock` files:
```bash
poetry add <package-name>
# or for a development dependency:
poetry add --group dev <package-name>
```

To constrain the version of a package, see the [Poetry versioning docs](https://python-poetry.org/docs/dependency-specification/#version-constraints).

**What does the lockfile `poetry.lock` do ?**

<details>
<summary>Expand for lockfile summary</summary>

The `poetry.lock` file ensures that everyone working on the project uses the
same versions of dependencies, which helps to avoid "it works on my machine"
problems.

It separates dependency-of-dependencies from human-specified dependencies in
`pyproject.toml`, and pins them to specific versions.

However the lockfile is a 'disposable' file in that it can be regenerated from
the `pyproject.toml` file if needed. The lockfile should always be committed to
version control. If there are merge conflicts in the lockfile, discard it and
regenerate it with `poetry lock`.

Finally, downstream users of the `BEstimate` package do not benefit from the
lockfile, and install the dependencies as specified in `pyproject.toml` (but not
dev dependencies).

</details>

**What does the `requirements.txt` do ?**

<details>
<summary>Expand for requirements.txt summary</summary>

The `requirements.txt` file is generated by Poetry and pre-commit as an artifact
and to allow developers to install dependencies with `pip install -r requirements.txt`
in environments where Poetry is not available. It should not be manually edited.

</details>

### Formatting and pre-commit hooks

This project uses [pre-commit](https://pre-commit.com/) to manage elements of
code formatting and linting. See the [One time setup](#one-time-setup) section
for installation and setup.

Pre-commit will run automatically on `git commit`. Generally pre-commit will
modify and correct files, these need to be staged again before the commit can
complete.

To run the pre-commit hooks manually:
```bash
pre-commit run -a
```

To push a commit while bypassing pre-commit (there are reasons to do this):
```bash
git commit --no-verify -m "My commit message"
```

The pre-commit configuration is in `.pre-commit-config.yaml` and includes:
- built-in hooks for checking for end-of-file newlines, trailing whitespace and
  ensuring valid JSON, YAML and TOML files
- [black](https://github.com/psf/black) - Python code formatter
- [flake8](https://flake8.pycqa.org/en/latest/) - Python code linter
- [poetry-check](https://python-poetry.org/docs/pre-commit-hooks/)
- [poetry-export](https://github.com/python-poetry/poetry-plugin-export) - generates `requirements.txt` from `poetry.lock`

### Testing

Tests are in the `tests/` directory and use [pytest](https://docs.pytest.org/).

To run the tests:

```bash
pytest tests/
```

Or with the development Docker image:

```bash
docker build -f Dockerfile-dev -t 'bestimate-dev:local' .
docker run -it --rm bestimate-dev:local pytest tests/
```

Or with the public image:

```bash
docker build -t bestimate:local -f Dockerfile .
# The tests don't exist in the image and pytest is not installed
docker run -it --rm \
    -v ./tests/:/opt/repo/tests \
    -w /opt/repo \
    bestimate:local \
    bash -c  'pip install pytest && python -m pytest tests/'
```

### CICD (Gitlab CI)

The CI in `.gitlab-ci.yml` uses the [CICD template repository](https://gitlab.internal.sanger.ac.uk/team113sanger/common/cicd-template) and includes the following **stages** that:
- **build** two Docker images, from `Dockerfile` and `Dockerfile.dev`
- **tests** runs `e2e`, `pytest` and `pre-commit` against the built images
- **publish**
    - if on a tag e.g `1.2.3` publishes the image to GitLab Container Registry as `<image>:<tag>` and publishes the package to GitLab PyPI as `<package>:<tag>`
    - if on `main` branch, publishes the image to Docker Hub as `<image>:latest`
    - if on `develop` branch, publishes the image to Docker Hub as `<image>:develop-latest`

Certain CI variables are maintained in [this repository's CICD settings](https://gitlab.internal.sanger.ac.uk/sci/BEstimate/-/settings/ci_cd#js-cicd-variables-settings). Of note is:
- `GITLAB_DEPLOY_TOKEN_RW` and `GITLAB_DEPLOY_USERNAME_RW` - used to authenticate with the GitLab container registry and PyPI (the `GITLAB_CI_TOKEN` doesn't have API write permissions)
- `DOCKER_HUB_USER` and `DOCKER_HUB_ACCESS_TOKEN` - used to authenticate with Docker Hub to allow pull images without interfering with the Sanger/DockerHub rate limits

## Git and Tagged releases
This repo the GitFlow branching model and uses [hubflow](https://datasift.github.io/gitflow/TheHubFlowTools.html) as a tool to enable this from the CLI.

```bash
# Switch to the develop branch
git checkout develop

# Start a new release branch e.g. 0.1.0 not v0.1.0
git hf release start <project_version>
```

Now, do the following things:
* `CHANGELOG.md`: Under the heading of the newest release version, describe what was changed, fixed, added.
* `pyproject.toml`: Increment the project version to the current release version
* Commit these changes
* Run `pre-commit run -a` to ensure no formatting issues

Finally

```bash
git hf release finish <project_version>
```
