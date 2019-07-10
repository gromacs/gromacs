# The myplugin module must be locatable by Python.
# If you configured CMake in the build directory ``/path/to/repo/build`` then,
# assuming you are in ``/path/to/repo``, run the tests with something like
#     PYTHONPATH=./cmake-build-debug/src/pythonmodule mpiexec -n 2 python -m mpi4py -m pytest tests/

# This test is not currently run automatically in any way. Build the module, point your PYTHONPATH at it,
# and run pytest in the tests directory.

import logging
import os

import gmxapi as gmx
from gmxapi.simulation.context import ParallelArrayContext
from gmxapi.simulation.workflow import WorkElement, from_tpr
from gmxapi import version as gmx_version
import pytest
from tests.conftest import withmpi_only

logging.getLogger().setLevel(logging.DEBUG)
# create console handler
ch = logging.StreamHandler()
ch.setLevel(logging.DEBUG)
# create formatter and add it to the handler
formatter = logging.Formatter('%(asctime)s:%(name)s:%(levelname)s: %(message)s')
ch.setFormatter(formatter)
# add the handlers to the logger
logging.getLogger().addHandler(ch)

logger = logging.getLogger()


def test_import():
    import myplugin
    assert myplugin


@pytest.mark.usefixtures("cleandir")
def test_harmonic_potential(tpr_filename):
    print("Testing plugin potential with input file {}".format(os.path.abspath(tpr_filename)))

    md = from_tpr(tpr_filename, append_output=False)

    context = ParallelArrayContext(md)
    with context as session:
        session.run()

    # Create a WorkElement for the potential
    params = {'sites': [1, 4],
              'R0': 2.0,
              'k': 10000.0}
    potential_element = WorkElement(namespace="myplugin",
                                    operation="create_restraint",
                                    params=params)
    # Note that we could flexibly capture accessor methods as workflow elements, too. Maybe we can
    # hide the extra Python bindings by letting myplugin.HarmonicRestraint automatically convert
    # to a WorkElement when add_dependency is called on it.
    potential_element.name = "harmonic_restraint"
    before = md.workspec.elements[md.name]
    md.add_dependency(potential_element)
    assert potential_element.name in md.workspec.elements
    assert potential_element.workspec is md.workspec
    after = md.workspec.elements[md.name]
    assert before is not after

    # Context will need to do these in __enter__
    # potential = myplugin.HarmonicRestraint()
    # potential.set_params(1, 4, 2.0, 10000.0)

    context = ParallelArrayContext(md)
    with context as session:
        session.run()


@pytest.mark.usefixtures("cleandir")
def test_ensemble_potential_nompi(tpr_filename):
    """Test ensemble potential without an ensemble.
    """
    print("Testing plugin potential with input file {}".format(os.path.abspath(tpr_filename)))

    assert gmx.version.api_is_at_least(0, 0, 5)
    md = from_tpr([tpr_filename], append_output=False)

    # Create a WorkElement for the potential
    params = {'sites': [1, 4],
              'nbins': 10,
              'binWidth': 0.1,
              'min_dist': 0.,
              'max_dist': 10.,
              'experimental': [1.] * 10,
              'nsamples': 1,
              'sample_period': 0.001,
              'nwindows': 4,
              'k': 10000.,
              'sigma': 1.}
    potential = WorkElement(namespace="myplugin",
                            operation="ensemble_restraint",
                            params=params)
    # Note that we could flexibly capture accessor methods as workflow elements, too. Maybe we can
    # hide the extra Python bindings by letting myplugin.HarmonicRestraint automatically convert
    # to a WorkElement when add_dependency is called on it.
    potential.name = "ensemble_restraint"
    md.add_dependency(potential)

    context = ParallelArrayContext(md)

    with context as session:
        session.run()


@withmpi_only
@pytest.mark.usefixtures("cleandir")
def test_ensemble_potential_withmpi(tpr_filename):
    import os

    from mpi4py import MPI
    rank = MPI.COMM_WORLD.Get_rank()

    rank_dir = os.path.join(os.getcwd(), str(rank))
    os.mkdir(rank_dir)

    logger.info("Testing plugin potential with input file {}".format(os.path.abspath(tpr_filename)))

    assert gmx_version.api_is_at_least(0, 0, 5)
    md = from_tpr([tpr_filename, tpr_filename], append_output=False)

    # Create a WorkElement for the potential
    params = {'sites': [1, 4],
              'nbins': 10,
              'binWidth': 0.1,
              'min_dist': 0.,
              'max_dist': 10.,
              'experimental': [0.5] * 10,
              'nsamples': 1,
              'sample_period': 0.001,
              'nwindows': 4,
              'k': 10000.,
              'sigma': 1.}

    potential = WorkElement(namespace="myplugin",
                            operation="ensemble_restraint",
                            params=params)
    # Note that we could flexibly capture accessor methods as workflow elements, too. Maybe we can
    # hide the extra Python bindings by letting myplugin.HarmonicRestraint automatically convert
    # to a WorkElement when add_dependency is called on it.
    potential.name = "ensemble_restraint"
    md.add_dependency(potential)

    context = ParallelArrayContext(md)
    with context as session:
        session.run()
