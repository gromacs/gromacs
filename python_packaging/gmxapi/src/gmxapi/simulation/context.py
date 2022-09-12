#
# This file is part of the GROMACS molecular simulation package.
#
# Copyright 2019- The GROMACS Authors
# and the project initiators Erik Lindahl, Berk Hess and David van der Spoel.
# Consult the AUTHORS/COPYING files and https://www.gromacs.org for details.
#
# GROMACS is free software; you can redistribute it and/or
# modify it under the terms of the GNU Lesser General Public License
# as published by the Free Software Foundation; either version 2.1
# of the License, or (at your option) any later version.
#
# GROMACS is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public
# License along with GROMACS; if not, see
# https://www.gnu.org/licenses, or write to the Free Software Foundation,
# Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA.
#
# If you want to redistribute modifications to GROMACS, please
# consider that scientific software is very special. Version
# control is crucial - bugs must be traceable. We will be happy to
# consider code for inclusion in the official distribution, but
# derived work must not be called official GROMACS. Details are found
# in the README & COPYING files - if they are missing, get the
# official version at https://www.gromacs.org.
#
# To help us fund GROMACS development, we humbly ask that you cite
# the research papers on the package. Check out https://www.gromacs.org.

# This file is based on the Kasson Lab gmxapi project release 0.0.7.4.
# https://github.com/kassonlab/gmxapi/blob/v0.0.7.4/src/gmx/context.py
# https://github.com/kassonlab/gmxapi/blob/v0.0.7.4/LICENSE

"""
Execution Context
=================
"""

__all__ = ['Context']

import importlib
import os
import warnings
import tempfile

from gmxapi import exceptions
from gmxapi import logger as root_logger
import gmxapi._gmxapi as _gmxapi

# Module-level logger
from gmxapi.exceptions import ApiError
from gmxapi.exceptions import DataShapeError

logger = root_logger.getChild('simulation.context')
logger.info('Importing {}'.format(__name__))


def _load_tpr(self, element):
    """Implement the gromacs.load_tpr operation.

    Updates the minimum width of the workflow parallelism. Does not add any API object to the graph.

    Arguments:
        self: The Context in which this operation is being loaded.
        element: WorkElement specifying the operation.

    Returns:
        A Director that the Context can use in launching the Session.
    """
    class Builder(object):
        def __init__(self, tpr_list):
            logger.debug("Loading tpr builder for tpr_list {}".format(tpr_list))
            self.tpr_list = tpr_list
            self.subscribers = []
            self.width = len(tpr_list)

        def add_subscriber(self, builder):
            builder.infile = self.tpr_list
            self.subscribers.append(builder)

        def build(self, dag):
            width = len(self.tpr_list)
            for builder in self.subscribers:
                builder.width = width
            if 'width' in dag.graph:
                width = max(width, dag.graph['width'])
            dag.graph['width'] = width

    return Builder(element.params['input'])

def _md(context, element):
    """Implement the gmxapi.md operation by returning a builder that can populate a data flow graph for the element.

    Inspects dependencies to set up the simulation runner.

    The graph node created will have `launch` and `run` attributes with function references, and a `width`
    attribute declaring the workflow parallelism requirement.

    Arguments:
        context: The Context in which this operation is being loaded.
        element: WorkElement specifying the operation.

    Returns:
        A Director that the Context can use in launching the Session.
    """
    class Builder(object):
        """Translate md work element to a node in the session's DAG."""
        def __init__(self, element):
            try:
                self.name = element.name
                # Note that currently the calling code is in charge of subscribing this builder to its dependencies.
                # A list of tpr files will be set when the calling code subscribes this builder to a tpr provider.
                self.infile = None
                # Other dependencies in the element may register potentials when subscribed to.
                self.potential = []
                self.input_nodes = []
                self.runtime_params = element.params
            except AttributeError:
                raise exceptions.ValueError("object provided does not seem to be a WorkElement.")
        def add_subscriber(self, builder):
            """The md operation does not yet have any subscribeable facilities."""
            pass
        def build(self, dag):
            """Add a node to the graph that, when launched, will construct a simulation runner.

            Complete the definition of appropriate graph edges for dependencies.

            The launch() method of the added node creates the runner from the tpr file for the current rank and adds
            modules from the incoming edges.
            """
            if not (hasattr(dag, 'add_node')
                    and hasattr(dag, 'add_edge')
                    and hasattr(dag, 'graph')
                    and hasattr(dag, 'nodes')):
                raise exceptions.ProtocolError("Argument 'dag' must provide the DiGraph interface.")
            name = self.name
            dag.add_node(name)
            for neighbor in self.input_nodes:
                dag.add_edge(neighbor, name)
            infile = self.infile
            assert infile is not None
            potential_list = self.potential

            # We may have to widen the graph according to runtime_params.
            # We also ensure that the ensemble data is represented as a list of dict.
            if isinstance(self.runtime_params, list):
                runtime_params = self.runtime_params
                params_width = len(runtime_params)
            else:
                assert isinstance(self.runtime_params, dict)
                # First, determine the final width
                params_width = 1
                for key, value in self.runtime_params.items():
                    if isinstance(value, list):
                        value_length = len(value)
                        if params_width == 1:
                            params_width = value_length
                        elif (value_length > 1) and (params_width != value_length):
                            raise DataShapeError(
                                f'runtime parameter {key}={value} not compatible with ensemble width '
                                f'{params_width}'
                            )
                # Get the list of flattened params dicts
                runtime_params = [dict() for _ in range(params_width)]
                for param, value in self.runtime_params.items():
                    if not isinstance(value, list):
                        for member in range(params_width):
                            runtime_params[member][param] = value
                    elif len(value) == 1:
                        value = value[0]
                        for member in range(params_width):
                            runtime_params[member][param] = value
                    else:
                        assert isinstance(value, list)
                        value_length = len(value)
                        assert value_length > 1
                        assert value_length == params_width
                        for member in range(params_width):
                            runtime_params[member][param] = value[member]
            assert isinstance(runtime_params, list)
            assert params_width == len(runtime_params)

            input_width = len(infile)

            dag_width = dag.graph['width']
            required_width = max((dag_width, input_width, params_width))

            if required_width > 1:
                if (dag_width != 1) and (dag_width != required_width):
                    raise DataShapeError(
                        f'Inputs {infile} and params {runtime_params} cannot be fit to current data flow '
                        f'graph of width {dag_width}.'
                    )
                if input_width == 1:
                    infile = [infile[0]] * required_width
                elif input_width != required_width:
                    raise DataShapeError(
                        f'Inputs {infile} incompatible with required ensemble width {required_width}.'
                    )
                if params_width == 1:
                    _params: dict = runtime_params[0]
                    runtime_params = [_params.copy() for _ in range(required_width)]
                elif params_width != required_width:
                    raise DataShapeError(
                        f'Params {runtime_params} incompatible with required ensemble width {required_width}.'
                    )

            self.infile = infile
            self.runtime_params = runtime_params
            dag.graph['width'] = required_width

            # Provide closure with which to execute tasks for this node.
            # TODO(#4079): launch is explicitly a collaboration with the Context
            #  to build a Session, and should be decoupled from element.workspec._context
            def launch(rank=None):
                assert rank is not None
                # Note that gmxapi 0.0.7 does not support any consumers of mdrun, so the graph cannot be
                # further widened, and launch can assume that *rank* <= dag width.
                assert dag.graph['width'] == required_width
                assert len(self.infile) == required_width
                assert len(self.runtime_params) == required_width
                if rank > required_width:
                    raise ApiError(
                        f'Trying to launch MD on rank {rank}, but ensemble has width {required_width}.'
                    )
                infile = self.infile[rank]
                runtime_params = self.runtime_params[rank]

                # Copy and update, if required by `end_time` parameter.
                temp_filename = None
                if 'end_time' in runtime_params:
                    # Note that mkstemp returns a file descriptor as the first part of the tuple.
                    # We can make this cleaner in 0.0.7 with a separate node that manages the
                    # altered input.
                    _, temp_filename = tempfile.mkstemp(suffix='.tpr')
                    logger.debug('Updating input. Using temp file {}'.format(temp_filename))
                    _gmxapi.rewrite_tprfile(source=infile,
                                         destination=temp_filename,
                                         end_time=runtime_params['end_time'])
                    tpr_file = temp_filename
                else:
                    tpr_file = infile

                logger.info('Loading TPR file: {}'.format(tpr_file))
                system = _gmxapi.from_tpr(tpr_file)
                dag.nodes[name]['system'] = system
                mdargs = _gmxapi.MDArgs()
                mdargs.set(runtime_params)
                # Workaround to give access to plugin potentials used in a context.
                # See https://gitlab.com/gromacs/gromacs/-/issues/3145
                work_context: Context = element.workspec._context
                work_context.potentials = potential_list
                pycontext: _gmxapi.Context = work_context._api_context
                pycontext.setMDArgs(mdargs)
                for potential in potential_list:
                    pycontext.add_mdmodule(potential)
                dag.nodes[name]['session'] = system.launch(pycontext)
                dag.nodes[name]['close'] = dag.nodes[name]['session'].close

                if 'end_time' in runtime_params:
                    def special_close():
                        dag.nodes[name]['session'].close()
                        logger.debug("Unlinking temporary TPR file {}.".format(temp_filename))
                        os.unlink(temp_filename)
                    dag.nodes[name]['close'] = special_close
                else:
                    dag.nodes[name]['close'] = dag.nodes[name]['session'].close

                def runner():
                    """Currently we only support a single call to run."""
                    def done():
                        raise StopIteration()
                    # Replace the runner with a stop condition for subsequent passes.
                    dag.nodes[name]['run'] = done
                    return dag.nodes[name]['session'].run()
                dag.nodes[name]['run'] = runner
                return dag.nodes[name]['run']

            dag.nodes[name]['launch'] = launch

    return Builder(element)


def _get_mpi_ensemble_communicator(session_communicator, ensemble_size):
    """Get an ensemble communicator from an MPI communicator.

    An ensemble communicator is an object that implements mpi4py.MPI.Comm methods
    as described elsewhere in this documentation.

    :param session_communicator: MPI communicator with the interface described by mpi4py.MPI.Comm
    :param ensemble_size: number of ensemble members
    :return: communicator of described size on participating ranks and null communicator on any others.

    Must be called exactly once in every process in `communicator`. It is the
    responsibility of the calling code to refrain from running ensemble operations
    if not part of the ensemble. The calling code determines this by comparing its
    session_communicator.Get_rank() to ensemble_size. This is not a good solution
    because it assumes there is a single ensemble communicator and that ensemble
    work is allocated to ranks serially from session_rank 0. Future work might
    use process groups associated with specific operations in the work graph so
    that processes can check for membership in a group to decide whether to use
    their ensemble communicator. Another possibility would be to return None
    rather than a null communicator in processes that aren't participating in
    a given ensemble.
    """
    try:
        from mpi4py import MPI
    except ImportError:
        raise exceptions.FeatureNotAvailableError('MPI ensemble communicator requires mpi4py package.')

    session_size = session_communicator.Get_size()
    session_rank = session_communicator.Get_rank()

    # Check the ensemble "width" against the available parallelism
    if ensemble_size > session_size:
        msg = 'ParallelArrayContext requires a work array that fits in the MPI communicator: '
        msg += 'array width {} > size {}.'
        msg = msg.format(ensemble_size, session_size)
        raise exceptions.UsageError(msg)
    if ensemble_size < session_size:
        msg = 'MPI context is wider than necessary to run this work:  array width {} vs. size {}.'
        warnings.warn(msg.format(ensemble_size, session_size))

    # Create an appropriate sub-communicator for the present work. Extra ranks will be in a
    # sub-communicator with no work.
    if session_rank < ensemble_size:
        # The session launcher should maintain an inventory of the ensembles and
        # provide an appropriate tag, but right now we just have a sort of
        # Boolean: ensemble or not.
        color = 0
    else:
        color = MPI.UNDEFINED

    ensemble_communicator = session_communicator.Split(color, session_rank)
    try:
        ensemble_communicator_size = ensemble_communicator.Get_size()
        ensemble_communicator_rank = ensemble_communicator.Get_rank()
    except:
        warnings.warn("Possible API programming error: ensemble_communicator does not provide required methods...")
        ensemble_communicator_size = 0
        ensemble_communicator_rank = None
    logger.info("Session rank {} assigned to rank {} of subcommunicator {} of size {}".format(
        session_rank,
        ensemble_communicator_rank,
        ensemble_communicator,
        ensemble_communicator_size
    ))

    # There isn't a good reason to worry about special handling for a null communicator,
    # which we have to explicitly avoid "free"ing, so let's just get rid of it.
    # To do: don't even get the null communicator in the first place. Use a group and create instead of split.
    if ensemble_communicator == MPI.COMM_NULL:
        ensemble_communicator = _DummyCommunicator()

    return ensemble_communicator


class _DummyCommunicator(object):
    """Placeholder class for trivial communication resources.

    Simplifies logic in this module until communications infrastructure is more robust.
    """

    def __init__(self):
        import numpy
        self._numpy = numpy

    def Dup(self):
        return self

    def Free(self):
        return

    def Allreduce(self, send, recv):
        logger.debug("Faking an Allreduce for ensemble of size 1.")
        send_buffer = self._numpy.array(send, copy=False)
        recv_buffer = self._numpy.array(recv, copy=False)
        recv_buffer[:] = send_buffer[:]

    def Get_size(self):
        return 1

    def Get_rank(self):
        return 0

    def __str__(self):
        return '_DummyCommunicator'

    def __repr__(self):
        return '_DummyCommunicator()'


def _acquire_communicator(communicator=None):
    """Get a workflow level communicator for the session.

    This function is intended to be called by the __enter__ method that creates
    a session get a communicator instance. The `Free` method of the returned
    instance must be called exactly once. This should be performed by the
    corresponding __exit__ method.

    Arguments:
        communicator : a communicator to duplicate (optional)

    Returns:
        A communicator that must be explicitly freed by the caller.

    Currently only supports MPI multi-simulation parallelism dependent on
    mpi4py. The mpi4py package should be installed and built with compilers
    that are compatible with the gmxapi installation.

    If provided, `communicator` must provide the mpi4py.MPI.Comm interface.
    Returns either a duplicate of `communicator` or of MPI_COMM_WORLD if mpi4py
    is available. Otherwise, returns a mock communicator that can only manage
    sessions and ensembles of size 0 or 1.

    gmx behavior is undefined if launched with mpiexec and without mpi4py
    """

    if communicator is None:
        try:
            import mpi4py.MPI as _MPI
        except ImportError:
            _MPI = None
        if _MPI is None:
            logger.info("mpi4py is not available for default session communication.")
            communicator = _DummyCommunicator()
        else:
            communicator = _MPI.COMM_WORLD
    else:
        communicator = communicator

    try:
        new_communicator = communicator.Dup()
    except Exception as e:
        message = "Exception when duplicating communicator: {}".format(e)
        raise exceptions.ApiError(message)

    return new_communicator


def _get_ensemble_communicator(communicator, ensemble_size):
    """Provide ensemble_communicator feature in active_context, if possible.

    Must be called on all ranks in `communicator`. The communicator returned
    must be freed by a call to its `Free()` instance method. This function is
    best used in a context manager's `__enter__()` method so that the
    corresponding `context.Free()` can be called in the `__exit__` method.

    Arguments:
        communicator : session communicator for the session with the ensemble.
        ensemble_size : ensemble size of the requested ensemble communicator

    The ensemble_communicator feature should be present if the Context can
    provide communication between all ensemble members. The Context should
    determine this at the launch of the session and set the
    ``_session_ensemble_communicator`` attribute to provide an object that
    implements the same interface as an mpi4py.MPI.Comm object. Actually, this is
    a temporary shim, so the only methods that need to be available are `Get_size`,
    `Get_rank` and something that can be called as
    Allreduce(send, recv) where send and recv are objects providing the Python
    buffer interface.

    Currently, only one ensemble can be managed in a session.
    """
    # For trivial cases, don't bother trying to use MPI
    # Note: all ranks in communicator must agree on the size of the work!
    # Note: If running with a dummy session communicator in an MPI session (user error)
    # every rank will think it is the only rank and will try to perform the
    # same work.
    if communicator.Get_size() <= 1 or ensemble_size <= 1:
        message = "Getting trivial ensemble communicator for ensemble of size {}".format((ensemble_size))
        message += " for session rank {} in session communicator of size {}".format(
            communicator.Get_rank(),
            communicator.Get_size())
        logger.debug(message)
        ensemble_communicator = _DummyCommunicator()
    else:
        message = "Getting an MPI subcommunicator for ensemble of size {}".format(ensemble_size)
        message += " for session rank {} in session communicator of size {}".format(
            communicator.Get_rank(),
            communicator.Get_size())
        logger.debug(message)
        ensemble_communicator = _get_mpi_ensemble_communicator(communicator, ensemble_size)

    return ensemble_communicator

def _get_ensemble_update(context):
    """Set up a simple ensemble resource.

    The context should call this function once per session to get an `ensemble_update`
    function object.

    This is a draft of a Context feature that may not be available in all
    Context implementations. This factory function can be wrapped as a
    ``ensemble_update`` "property" in a Context instance method to produce a Python function
    with the signature ``update(context, send, recv, tag=None)``.

    This feature requires that the Context is capabable of providing the
    ensemble_communicator feature and the numpy feature.
    If both are available, the function object provided by
    ``ensemble_update`` provides
    the ensemble reduce operation used by the restraint potential plugin in the
    gmxapi sample_restraint repository. Otherwise, the provided function object
    will log an error and then raise an exception.

    gmxapi 0.0.5 and 0.0.6 MD plugin clients look for a member function named
    ``ensemble_update`` in the Context that launched them. In the future,
    clients will use session resources to access ensemble reduce operations.
    In the mean time, a transitional implementation can involve defining a
    ``ensemble_update`` property in the Context object that acts as a factory
    function to produce the reducing operation, if possible with the given
    resources.
    """
    try:
        import numpy
    except ImportError:
        message = "ensemble_update requires numpy, but numpy is not available."
        logger.error(message)
        raise exceptions.FeatureNotAvailableError(message)

    def _ensemble_update(active_context, send, recv, tag=None):
        assert not tag is None
        assert str(tag) != ''
        if not tag in active_context.part:
            active_context.part[tag] = 0
        logger.debug("Performing ensemble update.")
        active_context._session_ensemble_communicator.Allreduce(send, recv)
        buffer = numpy.array(recv, copy=False)
        buffer /= active_context.work_width
        suffix = '_{}.npz'.format(tag)
        # These will end up in the working directory and each ensemble member will have one
        filename = str("rank{}part{:04d}{}".format(active_context.rank, int(active_context.part[tag]), suffix))
        numpy.savez(filename, recv=recv)
        active_context.part[tag] += 1

    def _no_ensemble_update(active_context, send, recv, tag=None):
        message = "Attempt to call ensemble_update() in a Context that does not provide the operation."
        # If we confirm effective exception handling, remove the extraneous log.
        logger.error(message)
        raise exceptions.FeatureNotAvailableError(message)

    if context._session_ensemble_communicator is not None:
        functor = _ensemble_update
    else:
        functor = _no_ensemble_update
    return functor


class Context(object):
    """Manage an array of simulation work executing in parallel.

    .. deprecated:: 0.0.7

    This execution context implementation is imported from the gmxapi 0.0.7
    package for GROMACS 2019 and does not conform to the protocols of gmxapi 0.1+.
    It is used internally to support legacy code.

    The following features are subject to be changed or removed without further
    notice.

    Attributes:
        work :obj:`gmx.workflow.WorkSpec`: specification of work to be performed when a session is launched.
        rank : numerical index of the current worker in a running session (None if not running)
        work_width : Detected number of co-scheduled work elements.
        elements : dictionary of references to elements of the workflow.

    `rank`, `work_width`, and `elements` are empty or None until the work is processed, as during session launch.

    To run multiple simulations at a time, whether ensembles or independent
    simulations, Python should be invoked in an MPI environment with `mpi4py`.
    The Context can then use one MPI rank per simulation. It is recommended that
    GROMACS is compiled without MPI in this case. (Thread-MPI is recommended.)
    MPI domain decomposition is not explicitly supported with this Context
    implementation.

    Example:

        >>> from gmxapi.simulation.context import get_context
        >>> from gmxapi.simulation.workflow import from_tpr
        >>> work = from_tpr([tpr_filename, tpr_filename])
        >>> context = get_context(work)
        >>> with context as session:
        ...    session.run()
        ...    # The session is one abstraction too low to know what rank it is. It lets the spawning context manage
        ...    # such things.
        ...    # rank = session.rank
        ...    # The local context object knows where it fits in the global array.
        ...    rank = context.rank
        ...    output_path = os.path.join(context.workdir_list[rank], 'traj.trr')
        ...    assert(os.path.exists(output_path))
        ...    print('Worker {} produced {}'.format(rank, output_path))

    Warning:
        Do not run a gmxapi script in an MPI environment without `mpi4py`. gmxapi
        will not be able to determine that other processes are running the same
        script.

    Implementation notes:

    To produce a running session, the Context __enter__() method is called, according to the Python context manager
    protocol. At this time, the attached WorkSpec must be feasible on the available resources. To turn the specified
    work into an executable directed acyclic graph (DAG), handle objects for the elements in the work spec are sequenced
    in dependency-compatible order and the context creates a "builder" for each according to the element's operation.
    Each builder is subscribed to the builders of its dependency elements. The DAG is then assembled by calling each
    builder in sequence. A builder can add zero, one, or more nodes and edges to the DAG.

    The Session is then launched from the DAG. What happens next is implementation-dependent, and it may take a while for
    us to decide whether and how to standardize interfaces for the DAG nodes and edges and/or execution protocols. I
    expect each node will at least have a `launch()` method, but will also probably have typed input and output ports as well as some signalling.
    A sophisticated and abstract Session implementation could schedule work only to satisfy data dependencies of requested
    output upon request. Our immediate implementation will use the following protocol.

    Each node has a `launch()` method. When the session is entered, the `launch()` method is called for each node in
    dependency order. The launch method returns either a callable (`run()` function) or None, raising an exception in
    case of an error. The sequence of callables is stored by the Session. When Session.run() is called, the
    sequence of callables is called in order. If StopIteration is raised by the callable, it is removed from the sequence.
    The sequence is processed repeatedly until there are no more callables.

    Note that this does not rigorously handle races or deadlocks, or flexibility in automatically chasing dependencies. A more
    thorough implementation could recursively call launch on dependencies (launch could be idempotent or involve some
    signalling to dependents when complete), run calls could be entirely event driven, and/or nodes could "publish"
    output (including just a completion message), blocking for acknowledgement before looking for the next set of subscribed inputs.
    """

    def __init__(self, work=None, workdir_list=None, communicator=None):
        """Create manager for computing resources.

        Does not initialize resources because Python objects by themselves do
        not have a good way to deinitialize resources. Instead, resources are
        initialized using the Python context manager protocol when sessions are
        entered and exited.

        Appropriate computing resources need to be knowable when the Context is created.

        Keyword Arguments:
            work : work specification with which to initialize this context
            workdir_list : deprecated
            communicator : non-owning reference to a multiprocessing communicator
                from which ensemble and simulation session resources will be derived

        If provided, communicator must implement the mpi4py.MPI.Comm interface. The
        Context will use this communicator as the parent for subcommunicators
        used when launching sessions. If provided, communicator is owned by the
        caller, and must be freed by the caller after any sessions are closed.
        By default, the Context will get a reference to MPI_COMM_WORLD, which
        will be freed when the Python process ends and cleans up its resources.
        The communicator stored by the Context instance will not be used directly,
        but will be duplicated when launching sessions using ``with``.
        """

        # self.__context_array = list([Context(work_element) for work_element in work])
        from .workflow import WorkSpec

        # Until better Session abstraction exists at the Python level, a
        # _session_communicator attribute will be added to and removed from the
        # context at session entry and exit. If necessary, a _session_ensemble_communicator
        # will be split from _session_communicator for simulation ensembles
        # present in the specified work.
        self.__communicator = communicator

        self.__work = WorkSpec()
        self.__workdir_list = workdir_list

        self._session = None
        # Note: this attribute is a detail of MPI-based Contexts. Client access
        # is subject to change.
        self.rank = None

        # `work_width` notes the required width of an array of synchronous tasks to perform the specified work.
        # As work elements are processed, self.work_width will be increased as appropriate.
        self.work_width = None

        # initialize the operations map. May be extended during the lifetime of a Context.
        # Note that there may be a difference between built-in operations provided by this module and
        # additional operations registered at run time.
        self.__operations = dict()
        # The map contains a builder for each operation. The builder is created by passing the element to the function
        # in the map. The object returned must have the following methods:
        #
        #   * add_subscriber(another_builder) : allow other builders to subscribe to this one.
        #   * build(dag) : Fulfill the builder responsibilities by adding an arbitrary number of nodes and edges to a Graph.
        #
        # The gmxapi namespace of operations should be consistent with a specified universal set of functionalities
        self.__operations['gmxapi'] = {'md': lambda element : _md(self, element),
                                      }
        # Even if TPR file loading were to become a common and stable enough operation to be specified in
        # an API, it is unlikely to be implemented by any code outside of GROMACS, so let's not clutter
        # a potentially more universal namespace.
        self.__operations['gromacs'] = {'load_tpr': lambda element : _load_tpr(self, element),
                                       }

        # Right now we are treating workspec elements and work DAG nodes as equivalent, but they are
        # emphatically not intended to be tightly coupled. The work specification is intended to be
        # simple, user-friendly, general, and easy-to-implement. The DAG is an implementation detail
        # and may differ across context types. It is likely to have stronger typing of nodes and/or
        # edges. It is not yet specified whether we should translate the work into a graph before, after,
        # or during processing of the elements, so it is not yet known whether we will need facilities
        # to allow cross-referencing between the two graph-type structures. If we instantiate API objects
        # as we process work elements, and the DAG in a context deviates from the work specification
        # topology, we would need to use named dependencies to look up objects to bind to. Such facilities
        # could be hidden in the WorkElement class(es), too, to delegate code away from the Context as a
        # container class growing without bounds...
        # In actuality, we will have to process the entire workspec to some degree to make sure we can
        # run it on the available resources.
        self.elements = None

        # This setter must be called after the operations map has been populated.
        self.work = work

        try:
            self._api_context = _gmxapi.Context()
        except Exception as e:
            logger.error('Got exception when trying to create default library context: ' + str(e))
            raise exceptions.ApiError('Uncaught exception in API object creation.') from e

    @property
    def work(self):
        return self.__work

    @work.setter
    def work(self, work):
        """Set `work` attribute.

        Raises:
            gmx.exceptions.ApiError: work is not compatible with schema or
                known operations.
            gmx.exceptions.UsageError: Context can not access operations in
                the name space given for an Element
            gmx.exceptions.ValueError: assignment operation cannot be performed
                for the provided object (rhs)

        For discussion on error handling, see https://github.com/kassonlab/gmxapi/issues/125
        """
        from .workflow import WorkSpec, WorkElement
        if work is None:
            return

        if isinstance(work, WorkSpec):
            workspec = work
        elif hasattr(work, 'workspec') and isinstance(work.workspec, WorkSpec):
            workspec = work.workspec
        else:
            raise exceptions.ValueError('work argument must provide a gmx.workflow.WorkSpec.')
        workspec._context = self

        # Make sure this context knows how to run the specified work.
        for e in workspec.elements:
            element = WorkElement.deserialize(workspec.elements[e])

            # Note: Non-built-in namespaces (non-native) are treated as modules to import.
            # Native namespaces may not be completely implemented in a particular version of a particular Context.
            if element.namespace in {'gmxapi', 'gromacs'}:
                assert element.namespace in self.__operations
                if not element.operation in self.__operations[element.namespace]:
                    # The requested element is a built-in operation but not available in this Context.
                    # element.namespace should be mapped, but not all operations are necessarily implemented.
                    logger.error("Operation {} not found in map {}".format(element.operation,
                                                                           str(self.__operations)))
                    # This check should be performed when deciding if the context is appropriate for the work.
                    # If we are just going to use a try/catch block for this test, then we should differentiate
                    # this exception from those raised due to incorrect usage.
                    # The exception thrown here may evolve with https://github.com/kassonlab/gmxapi/issues/125
                    raise exceptions.FeatureNotAvailableError(
                        'Specified work cannot be performed due to unimplemented operation {}.{}.'.format(
                            element.namespace,
                            element.operation))

            else:
                assert element.namespace not in {'gmxapi', 'gromacs'}

                # Don't leave an empty nested dictionary if we couldn't map the operation.
                if element.namespace in self.__operations:
                    namespace_map = self.__operations[element.namespace]
                else:
                    namespace_map = dict()

                # Set or update namespace map iff we have something new.
                if element.operation not in namespace_map:
                    try:
                        element_module = importlib.import_module(element.namespace)
                        element_operation = getattr(element_module, element.operation)
                    except ImportError as e:
                        raise exceptions.UsageError(
                            'Cannot find implementation for namespace {}. ImportError: {}'.format(
                                element.namespace,
                                str(e)))
                    except AttributeError:
                        raise exceptions.UsageError(
                            'Cannot find factory for operation {}.{}'.format(
                                element.namespace,
                                element.operation
                            )
                        )
                    namespace_map[element.operation] = element_operation
                    self.__operations[element.namespace] = namespace_map

        self.__work = workspec

    def add_operation(self, namespace, operation, get_builder):
        # noinspection PyUnresolvedReferences
        """Add a builder factory to the operation map.

        Extends the known operations of the Context by mapping an operation in a namespace to a function
        that returns a builder to process a work element referencing the operation. Must be called before
        the work specification is added, since the spec is inspected to confirm that the Context can run it.

        It may be more appropriate to add this functionality to the Context constructor or as auxiliary
        information in the workspec, or to remove it entirely; it is straight-forward to just add snippets
        of code to additional files in the working directory and to make them available as modules for the
        Context to import.

        Example:

            >>> # Import some custom extension code.
            >>> import myplugin
            >>> myelement = myplugin.new_element()
            >>> workspec = gmx.workflow.WorkSpec()
            >>> workspec.add_element(myelement)
            >>> context = gmx.context.ParallelArrayContext()
            >>> context.add_operation(myelement.namespace, myelement.operation, myplugin.element_translator)
            >>> context.work = workspec
            >>> with get_context() as session:
            ...    session.run()

        """
        if namespace not in self.__operations:
            if namespace in {'gmxapi', 'gromacs'}:
                raise exceptions.UsageError("Cannot add operations to built-in namespaces.")
            else:
                self.__operations[namespace] = dict()
        else:
            assert namespace in self.__operations

        if operation in self.__operations[namespace]:
            raise exceptions.UsageError("Operation {}.{} already defined in this context.".format(namespace, operation))
        else:
            self.__operations[namespace][operation] = get_builder

    # Set up a simple ensemble resource
    # This should be implemented for Session, not Context, and use an appropriate subcommunicator
    # that is created and freed as the Session launches and exits.
    def ensemble_update(self, send, recv, tag=None):
        """Implement the ensemble_update member function that gmxapi through 0.0.6 expects.

        """
        # gmxapi through 0.0.6 expects to bind to this member function during "build".
        # This behavior needs to be deprecated (bind during launch, instead), but this
        # dispatching function should be an effective placeholder.
        if tag is None or str(tag) == '':
            raise exceptions.ApiError("ensemble_update must be called with a name tag.")
        # __ensemble_update is an attribute, not an instance function, so we need to explicitly pass 'self'
        return self.__ensemble_update(self, send, recv, tag)

    def __enter__(self):
        """Implement Python context manager protocol, producing a Session.

        The work specified in this Context is inspected to build a directed
        acyclic dependency graph (DAG). A Session is launched after reconciling
        the configured work with available computing resources. Each time the
        `run()` method of the Session is called, the graph is executed, and any
        nodes that have no more work to do are pruned from the graph.

        In practice, there are not currently any types of work implemented that
        do not complete on the first pass.

        Returns:
            Session object that can be run and/or inspected.

        Additional API operations are possible while the Session is active. When used as a Python context manager,
        the Context will close the Session at the end of the `with` block by calling `__exit__`.

        Note: this is probably where we will have to process the work specification to determine whether we
        have appropriate resources (such as sufficiently wide parallelism). Until we have a better Session
        abstraction, this means the clean approach should take two passes to first build a DAG and then
        instantiate objects to perform the work. In the first implementation, we kind of muddle things into
        a single pass.

        TODO(#3145, #4079): Use an explicit function call for context manager behavior.
        By implementing __enter__ and __exit__ directly on the class, we prevent
        ourselves from cleanly nesting context managers for managed resources,
        and we limit our ability to decouple the code involved in setting up and
        tearing down a Session.
        """
        # Cache the working directory from which we were launched
        # so that __exit__() can give us proper context management behavior.
        self.__initial_cwd = os.getcwd()

        try:
            import networkx as nx
            from networkx import DiGraph as _Graph
        except ImportError:
            raise exceptions.FeatureNotAvailableError("gmx requires the networkx package to execute work graphs.")

        logger.debug("Launching session from {}".format(self.__initial_cwd))

        if self._session is not None:
            raise exceptions.Error('Already running.')
        if self.work is None:
            raise exceptions.UsageError('No work to perform!')

        # Set up the global and local context.
        # Check the global MPI configuration
        # Since the Context doesn't have a destructor, if we use an MPI communicator at this scope then
        # it has to be owned and managed outside of Context.
        self._session_communicator = _acquire_communicator(self.__communicator)
        context_comm_size = self._session_communicator.Get_size()
        context_rank = self._session_communicator.Get_rank()
        self.rank = context_rank
        # self._communicator = communicator
        logger.debug("Context rank {} in context {} of size {}".format(context_rank,
                                                                       self._session_communicator,
                                                                       context_comm_size))

        ###
        # Process the work specification.
        ###
        logger.debug("Processing workspec:\n{}".format(str(self.work)))

        # Get a builder for DAG components for each element
        builders = {}
        builder_sequence = []
        for element in self.work:
            # dispatch builders for operation implementations
            try:
                new_builder = self.__operations[element.namespace][element.operation](element)
                assert hasattr(new_builder, 'add_subscriber')
                assert hasattr(new_builder, 'build')

                logger.info("Collected builder for {}".format(element.name))
            except LookupError as e:
                request = '.'.join([element.namespace, element.operation])
                message = 'Could not find an implementation for the specified operation: {}. '.format(request)
                message += str(e)
                raise exceptions.ApiError(message)
            # Subscribing builders is the Context's responsibility because otherwise the builders
            # don't know about each other. Builders should not depend on the Context unless they
            # are a facility provided by the Context, in which case they may be member functions
            # of the Context. We will probably need to pass at least some
            # of the Session to the `launch()` method, though...
            dependencies = element.depends
            for dependency in dependencies:
                # If a dependency is a list, assume it is an "ensemble" of dependencies
                # and pick the element for corresponding to the local rank.
                if isinstance(dependency, (list, tuple)):
                    assert len(dependency) > context_rank
                    name = str(dependency[context_rank])
                else:
                    name = dependency
                logger.info("Subscribing {} to {}.".format(element.name, name))
                builders[name].add_subscriber(new_builder)
            builders[element.name] = new_builder
            builder_sequence.append(element.name)
            logger.debug("Appended {} to builder sequence.".format(element.name))

        # Call the builders in dependency order
        # Note: session_communicator is available, but ensemble_communicator has not been created yet.
        graph = _Graph(width=1)
        logger.info("Building sequence {}".format(builder_sequence))
        for name in builder_sequence:
            builder = builders[name]
            logger.info("Building {}".format(builder))
            logger.debug("Has build attribute {}.".format(builder.build))
            builder.build(graph)
            logger.debug("Built.")
        self.work_width = graph.graph['width']

        # Prepare working directories. This should probably be moved to some aspect of the Session and either
        # removed from here or made more explicit to the user.
        workdir_list = self.__workdir_list
        if workdir_list is None:
            workdir_list = [os.path.join('.', str(i)) for i in range(self.work_width)]
        self.__workdir_list = list([os.path.abspath(dir) for dir in workdir_list])

        # For gmxapi 0.0.6, all ranks have a session_ensemble_communicator
        self._session_ensemble_communicator = _get_ensemble_communicator(self._session_communicator, self.work_width)
        self.__ensemble_update = _get_ensemble_update(self)
        self.part = {}

        # launch() is currently a method of gmx.core.MDSystem and returns a gmxapi::Session.
        # MDSystem objects are obtained from gmx.core.from_tpr(). They also provide add_potential().
        # gmxapi::Session objects are exposed as gmx.core.MDSession and provide run() and close() methods.
        #
        # Here, I want to find the input appropriate for this rank and get an MDSession for it.
        # E.g. Make a pass that allows meta-objects to bind (setting md_proxy._input_tpr and md_proxy._plugins,
        # and then call a routine implemented by each object to run whatever protocol it needs, such
        # as `system = gmx.core.from_tpr(md._input_tpr); system.add_potential(md._plugins)
        # For future design plans, reference https://github.com/kassonlab/gmxapi/issues/65
        #
        # This `if` condition is currently the thing that ultimately determines whether the
        # rank attempts to do work.
        if context_rank < self.work_width:
            # print(graph)
            logger.debug(("Launching graph {}.".format(graph.graph)))
            logger.debug("Graph nodes: {}".format(str(list(graph.nodes))))
            logger.debug("Graph edges: {}".format(str(list(graph.edges))))

            logger.info("Launching work on context rank {}, subcommunicator rank {}.".format(
                self.rank,
                self._session_ensemble_communicator.Get_rank()))

            # Launch the work for this rank
            self.workdir = self.__workdir_list[self.rank]
            if os.path.exists(self.workdir):
                if not os.path.isdir(self.workdir):
                    raise exceptions.Error('{} is not a valid working directory.'.format(self.workdir))
            else:
                os.mkdir(self.workdir)
            os.chdir(self.workdir)
            logger.info('rank {} changed directory to {}'.format(self.rank, self.workdir))
            sorted_nodes = nx.topological_sort(graph)
            runners = []
            closers = []
            for name in sorted_nodes:
                launcher = graph.nodes[name]['launch']
                # TODO(#4079): launch is explicitly a collaboration with the Context
                #  to build a Session, which should provide SessionResources here.
                runner = launcher(self.rank)
                if runner is not None:
                    runners.append(runner)
                    closers.append(graph.nodes[name]['close'])

            # Get a session object to return. It must simply provide a `run()` function.
            context = self # Part of workaround for bug gmxapi-214
            class Session(object):
                def __init__(self, runners, closers):
                    self.runners = list(runners)
                    self.closers = list(closers)

                def run(self):
                    # Note we are not following the documented protocol of running repeatedly yet.
                    to_be_deleted = []
                    for i, runner in enumerate(self.runners):
                        try:
                            runner()
                        except StopIteration:
                            to_be_deleted.insert(0, i)
                    for i in to_be_deleted:
                        del self.runners[i]
                    return True

                def close(self):
                    for close in self.closers:
                        logger.debug("Closing node: {}".format(close))
                        close()
                    # Workaround for bug gmxapi-214
                    if not _gmxapi.has_feature('0.0.7-bugfix-https://github.com/kassonlab/gmxapi/issues/214'):
                        context._api_object = _gmxapi.Context()

            self._session = Session(runners, closers)
        else:
            logger.info("Context rank {} has no work to do".format(self.rank))

            context = self # Part of workaround for bug gmxapi-214
            class NullSession(object):
                def run(self):
                    logger.info("Running null session on rank {}.".format(self.rank))
                    return True
                def close(self):
                    logger.info("Closing null session.")
                    # Workaround for bug gmxapi-214
                    if not _gmxapi.has_feature('0.0.7-bugfix-https://github.com/kassonlab/gmxapi/issues/214'):
                        context._api_object = _gmxapi.Context()
                    return

            self._session = NullSession()
            self._session.rank = self.rank

        # Make sure session has started on all ranks before continuing?

        self._session.graph = graph
        return self._session

    def __exit__(self, exception_type, exception_value, traceback):
        """Implement Python context manager protocol."""
        logger.info("Exiting session on context rank {}.".format(self.rank))
        if self._session is not None:
            logger.info("Calling session.close().")
            self._session.close()
            self._session = None
        else:
            # Note: we should not have a None session but rather an API-compliant Session that just has no work.
            # Reference: https://github.com/kassonlab/gmxapi/issues/41
            logger.info("No _session known to context or session already closed.")
        if hasattr(self, '_session_ensemble_communicator'):
            if self._session_communicator is not None:
                logger.info("Freeing sub-communicator {} on rank {}".format(
                    self._session_ensemble_communicator,
                    self.rank))
                self._session_ensemble_communicator.Free()
            else:
                logger.debug('"None" ensemble communicator does not need to be "Free"d.')
            del self._session_ensemble_communicator
        else:
            logger.debug("No ensemble subcommunicator on context rank {}.".format(self.rank))

        logger.debug("Freeing session communicator.")
        self._session_communicator.Free()
        logger.debug("Deleting session communicator reference.")
        del self._session_communicator

        os.chdir(self.__initial_cwd)
        logger.info("Session closed on context rank {}.".format(self.rank))
        # Note: Since sessions running in different processes can have different
        # work, sessions have not necessarily ended on all ranks. As a result,
        # starting another session on the same resources could block until the
        # resources are available.

        # Python context managers return False when there were no exceptions to handle.
        return False


def get_context(work=None):
    """Get a concrete Context object.

    Args:
        work (gmx.workflow.WorkSpec): runnable work as a valid gmx.workflow.WorkSpec object

    Returns:
        An object implementing the :py:class:`gmx.context.Context` interface, if possible.

    Raises:
        gmx.exceptions.ValueError if an appropriate context for ``work`` could not be loaded.

    If work is provided, return a Context object capable of running the provided work or produce an error.

    The semantics for finding Context implementations needs more consideration, and a more informative exception
    is likely possible.

    A Context can run the provided work if

      * the Context supports can resolve all operations specified in the elements
      * the Context supports DAG topologies implied by the network of dependencies
      * the Context supports features required by the elements with the specified parameters,
        such as synchronous array jobs.

    """
    # We need to define an interface for WorkSpec objects so that we don't need
    # to rely on typing and inter-module dependencies.
    from .workflow import WorkSpec
    workspec = None
    if work is not None:
        if isinstance(work, WorkSpec):
            workspec = work
        elif hasattr(work, 'workspec') and isinstance(work.workspec,
                                                      WorkSpec):
            workspec = work.workspec
        else:
            raise exceptions.ValueError('work argument must provide a gmx.workflow.WorkSpec.')
    if workspec is not None and \
            hasattr(workspec, '_context') and \
            workspec._context is not None:
        context = workspec._context
    else:
        context = Context(work=workspec)

    return context
