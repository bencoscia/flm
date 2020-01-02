"""
Unit and regression test for the flm package.
"""

# Import package, test suite, and other packages as needed
import flm
import pytest
import sys

def test_flm_imported():
    """Sample test, will always pass so long as import statement worked"""
    assert "flm" in sys.modules
