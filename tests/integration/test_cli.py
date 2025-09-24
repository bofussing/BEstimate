import subprocess
import shlex
import shutil

import pytest

import BEstimate
from BEstimate import constants

MODULE_NAME = BEstimate.__name__
PRIMARY_PROGRAM_NAME = constants.PROGRAM_NAME
ALL_PROGRAM_NAMES = [
    PRIMARY_PROGRAM_NAME,
    constants.SECONDARY_PROGRAM_NAME_X_GENOME,
    constants.SECONDARY_PROGRAM_NAME_X_CRISPRANALYZER,
]

# HELPERS


def get_subprocess_message(subproces_result: subprocess.CompletedProcess) -> str:
    indent = " " * 2

    msg = (
        f"Error running CLI command. "
        f"{indent}Command: {subproces_result.args}\n"
        f"{indent}Return code: {subproces_result.returncode}\n"
        f"{indent}Stdout: {subproces_result.stdout!r}\n"
        f"{indent}Stderr: {subproces_result.stderr!r}"
    )
    return msg


# TESTS


def test_python_dash_m__version():
    # Precondition
    # We assume the system has python installed but occasionally the binary may
    # be named python3 with no python binary.
    python_exec = "python"
    if shutil.which("python") is None:
        assert (
            shutil.which("python3") is not None
        ), "Python is not installed or not in PATH"
        python_exec = "python3"

    # Given
    cmd = f"{python_exec} -m {MODULE_NAME} --version"
    expected_version = BEstimate.__version__

    # When
    subproces_result = subprocess.run(shlex.split(cmd), capture_output=True, text=True)

    # Then
    errmsg = get_subprocess_message(subproces_result)
    assert subproces_result.returncode == 0, errmsg
    assert PRIMARY_PROGRAM_NAME in subproces_result.stdout
    assert expected_version in subproces_result.stdout


@pytest.mark.parametrize(
    "program_name", [pytest.param(name, id=name) for name in ALL_PROGRAM_NAMES]
)
def test_cli_on_path(program_name: str):
    # When
    result = shutil.which(program_name)

    assert result is not None, f"{program_name} is not in PATH, has the name changed?"


@pytest.mark.parametrize(
    "program_name", [pytest.param(name, id=name) for name in ALL_PROGRAM_NAMES]
)
def test_cli__version(program_name: str):
    # Given
    cmd = f"{program_name} --version"
    expected_version = BEstimate.__version__

    # When
    subproces_result = subprocess.run(shlex.split(cmd), capture_output=True, text=True)

    # Then
    errmsg = get_subprocess_message(subproces_result)
    assert subproces_result.returncode == 0, errmsg
    assert program_name in subproces_result.stdout
    assert expected_version in subproces_result.stdout


@pytest.mark.parametrize(
    "program_name", [pytest.param(name, id=name) for name in ALL_PROGRAM_NAMES]
)
def test_cli__help(program_name: str):
    # Given
    cmd = f"{program_name} --help"

    # When
    subproces_result = subprocess.run(shlex.split(cmd), capture_output=True, text=True)

    # Then
    errmsg = get_subprocess_message(subproces_result)
    assert subproces_result.returncode == 0, errmsg
    assert program_name in subproces_result.stdout
