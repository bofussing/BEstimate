import pytest
import re

# https://gitlab.internal.sanger.ac.uk/help/user/packages/pypi_repository/index.md#ensure-your-version-string-is-valid
SEMVER_REGEX = "".join(
    [
        r"(?:",
        r"(?:([0-9]+)!)?",  # epoch
        r"([0-9]+(?:\.[0-9]+)*)",  # release segment
        r"([-_\.]?((a|b|c|rc|alpha|beta|pre|preview))[-_\.]?([0-9]+)?)?",  # pre-release
        r"((?:-([0-9]+))|(?:[-_\.]?(post|rev|r)[-_\.]?([0-9]+)?))?",  # post release
        r"([-_\.]?(dev)[-_\.]?([0-9]+)?)?",  # dev release
        r"(?:\+([a-z0-9]+(?:[-_\.][a-z0-9]+)*))?",  # local version
        r")",
    ]
)
MY_PACKAGE_NAME = "BEstimate"


@pytest.mark.filterwarnings("ignore")
def test_package_structure_by_importing_from_src():
    from BEstimate.crispr_analyser import utils

    assert bool(utils)


def test_package_has_package_name():
    # When
    my_package = pytest.importorskip(MY_PACKAGE_NAME, reason="Package not installed")

    # Then
    err_msg = "Package should have __package_name__ attribute"
    assert hasattr(my_package, "__package_name__"), err_msg
    assert my_package.__package_name__ == MY_PACKAGE_NAME


def test_package_has_compatible_version():
    # Given
    rgx_pattern = re.compile(SEMVER_REGEX, re.VERBOSE)
    my_package = pytest.importorskip(MY_PACKAGE_NAME, reason="Package not installed")

    # When
    has_version = hasattr(my_package, "__version__")

    # Then
    assert has_version, "Package should have __version__ attribute"

    # Finally
    version = my_package.__version__
    msg1 = f"Version {version} is not compatible with PEP440"
    msg2 = "Version is 0.0.0. Did you forget to set the version in pyproject.toml"
    assert rgx_pattern.match(version), msg1
    assert version != "0.0.0", msg2
