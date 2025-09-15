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


@pytest.mark.filterwarnings("ignore")
def test_package_structure_by_importing_from_src():
    from BEstimate.crispr_analyser import utils

    assert bool(utils)

def test_version():
    try:
        import BEstimate
    except ImportError:
        should_raise = True
        version = None
    else:
        should_raise = False
        version = BEstimate.__version__

    if should_raise:
        assert False, "Could not import package. Can't test version."

    err_msg = "Version is 0.0.0. Did you forget to set the version in pyproject.toml"
    assert version != "0.0.0" and version is not None, err_msg


def test_package_name():
    try:
        import BEstimate
    except ImportError:
        should_raise = True
        package_name = None
    else:
        should_raise = False
        package_name = BEstimate.__package_name__

    if should_raise:
        assert False, "Could not import package. Can't test package name."

    expected_package_name = "BEstimate"

    assert package_name == expected_package_name


def test_version_is_compatible():
    import BEstimate

    # Given
    rgx_pattern = re.compile(SEMVER_REGEX, re.VERBOSE)

    # When
    version = BEstimate.__version__

    # Then
    msg = f"Version {version} is not compatible with PEP440"
    assert rgx_pattern.match(version), msg
