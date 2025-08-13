import pytest

# List of dependencies
dependencies = [
    "GEOS5FP",
    "numpy",
    "pandas",
    "PTJPL",
    "rasters",
]

# Generate individual test functions for each dependency
@pytest.mark.parametrize("dependency", dependencies)
def test_dependency_import(dependency):
    __import__(dependency)
