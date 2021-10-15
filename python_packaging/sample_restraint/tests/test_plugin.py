# The myplugin module must be locatable by Python.
# If you configured CMake in the build directory ``/path/to/repo/build`` then,
# assuming you are in ``/path/to/repo``, run the tests with something like
#     PYTHONPATH=./cmake-build-debug/src/pythonmodule mpiexec -n 2 python -m mpi4py -m pytest tests/

# This test is not currently run automatically in any way. Build the module, point your PYTHONPATH at it,
# and run pytest in the tests directory.

import logging
import os

try:
    import mpi4py.MPI as _MPI
except (ImportError, ModuleNotFoundError):
    _MPI = None

import gmxapi as gmx
from gmxapi.simulation.context import Context
from gmxapi.simulation.workflow import WorkElement, from_tpr
from gmxapi import version as gmx_version
import pytest

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
    # Suppress inspection warning outside of testing context.
    # noinspection PyUnresolvedReferences
    import myplugin
    assert myplugin


@pytest.mark.usefixtures("cleandir")
def test_binding_protocol(spc_water_box, mdrun_kwargs):
    """Test that gmxapi successfully attaches MD plugins."""
    import myplugin

    if _MPI is not None:
        _size = _MPI.COMM_WORLD.Get_size()
        _rank = _MPI.COMM_WORLD.Get_rank()
    else:
        _size = 1
        _rank = 0

    tpr_filename = spc_water_box
    logger.info("Testing plugin potential with input file {}".format(os.path.abspath(tpr_filename)))

    assert gmx.version.api_is_at_least(0, 2, 1)
    md = from_tpr([tpr_filename] * _size, append_output=False, **mdrun_kwargs)

    potential = WorkElement(namespace="myplugin",
                            operation="null_restraint",
                            params={'sites': [1, 4]})
    potential.name = "null restraint"
    md.add_dependency(potential)

    context = Context(md)

    with context as session:
        session.run()

    # See also #3038, #3145, #4079
    assert isinstance(context.potentials, list)
    assert len(context.potentials) > 0
    for restraint in context.potentials:
        if isinstance(restraint, myplugin.NullRestraint):
            assert gmx.version.api_is_at_least(0, 2, 1)
            assert restraint.count() > 1


@pytest.mark.usefixtures("cleandir")
def test_ensemble_potential_nompi(spc_water_box, mdrun_kwargs):
    """Test ensemble potential without an ensemble.
    """
    tpr_filename = spc_water_box
    logger.info("Testing plugin potential with input file {}".format(os.path.abspath(tpr_filename)))

    assert gmx.version.api_is_at_least(0, 0, 5)
    md = from_tpr([tpr_filename], append_output=False, **mdrun_kwargs)

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

    context = Context(md)

    with context as session:
        session.run()


@pytest.mark.withmpi_only
@pytest.mark.usefixtures("cleandir")
def test_ensemble_potential_withmpi(spc_water_box, mdrun_kwargs):
    tpr_filename = spc_water_box

    logger.info("Testing plugin potential with input file {}".format(os.path.abspath(tpr_filename)))

    assert gmx_version.api_is_at_least(0, 0, 5)
    md = from_tpr([tpr_filename, tpr_filename], append_output=False, **mdrun_kwargs)

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

    context = Context(md)
    with context as session:
        session.run()
