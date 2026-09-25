import pytest


def test_cli_help():
    """Test that the CLI help command works"""
    result = pytest.main(["scripts/cli.py", "--help"])
    assert result == 0