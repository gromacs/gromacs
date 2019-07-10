"""Configuration and fixtures for pytest."""

import pytest
import tempfile
import os

@pytest.fixture()
def cleandir():
    """Provide a clean working directory for tests.

    Example usage:

        import os
        import pytest

        @pytest.mark.usefixtures("cleandir")
        class TestDirectoryInit(object):
            def test_cwd_starts_empty(self):
                assert os.listdir(os.getcwd()) == []
                with open("myfile", "w") as f:
                    f.write("hello")

            def test_cwd_again_starts_empty(self):
                assert os.listdir(os.getcwd()) == []

    Ref: https://docs.pytest.org/en/latest/fixture.html#using-fixtures-from-classes-modules-or-projects
    """
    newpath = tempfile.mkdtemp()
    os.chdir(newpath)


try:
    from mpi4py import MPI
    withmpi_only = pytest.mark.skipif(not MPI.Is_initialized() or MPI.COMM_WORLD.Get_size() < 2,
                                      reason="Test requires at least 2 MPI ranks, but MPI is not initialized or too small.")
except ImportError:
    withmpi_only = pytest.mark.skip(reason="Test requires at least 2 MPI ranks, but mpi4py is not available.")

@pytest.fixture()
def tpr_filename():
    """Provide a sample TPR file by filename."""
    current_dir = os.path.dirname(__file__)
    file_path = os.path.join(current_dir, 'topol.tpr')
    return os.path.abspath(file_path)
