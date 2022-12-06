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

__all__ = ["Context"]

import importlib
import os
import warnings
import tempfile
from typing import Optional, TYPE_CHECKING, Any

import numpy as np

import gmxapi._gmxapi as _gmxapi
import gmxapi.utility
from gmxapi import exceptions
from gmxapi import logger as root_logger
from gmxapi.exceptions import ApiError
from gmxapi.exceptions import DataShapeError
from gmxapi.exceptions import ProtocolError

if TYPE_CHECKING:
    from gmxapi.runtime import (
        ResourceAllocation,
        ResourceAssignment,
    )
    import mpi4py.MPI

# Module-level logger
logger = root_logger.getChild("simulation.context")
logger.info("Importing {}".format(__name__))


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
            if "width" in dag.graph:
                width = max(width, dag.graph["width"])
            dag.graph["width"] = width

    return Builder(element.params["input"])


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
                raise exceptions.ValueError(
                    "object provided does not seem to be a WorkElement."
                )

        def add_subscriber(self, builder):
            """The md operation does not yet have any subscribeable facilities."""
            pass

        def build(self, dag):
            """Add a node to the graph that, when launched, will construct a simulation runner.

            Complete the definition of appropriate graph edges for dependencies.

            The launch() method of the added node creates the runner from the tpr file for the current rank and adds
            modules from the incoming edges.
            """
            if not (
                hasattr(dag, "add_node")
                and hasattr(dag, "add_edge")
                and hasattr(dag, "graph")
                and hasattr(dag, "nodes")
            ):
                raise exceptions.ProtocolError(
                    "Argument 'dag' must provide the DiGraph interface."
                )
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
                                f"runtime parameter {key}={value} not compatible with ensemble width "
                                f"{params_width}"
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

            dag_width = dag.graph["width"]
            required_width = max((dag_width, input_width, params_width))

            if required_width > 1:
                if (dag_width != 1) and (dag_width != required_width):
                    raise DataShapeError(
                        f"Inputs {infile} and params {runtime_params} cannot be fit to current data flow "
                        f"graph of width {dag_width}."
                    )
                if input_width == 1:
                    infile = [infile[0]] * required_width
                elif input_width != required_width:
                    raise DataShapeError(
                        f"Inputs {infile} incompatible with required ensemble width {required_width}."
                    )
                if params_width == 1:
                    _params: dict = runtime_params[0]
                    runtime_params = [_params.copy() for _ in range(required_width)]
                elif params_width != required_width:
                    raise DataShapeError(
                        f"Params {runtime_params} incompatible with required ensemble width {required_width}."
                    )

            self.infile = infile
            self.runtime_params = runtime_params
            dag.graph["width"] = required_width

            # Provide closure with which to execute tasks for this node.
            # TODO(#4079): launch is explicitly a collaboration with the Context
            #  to build a Session, and should be decoupled from element.workspec._context
            def launch(ensemble_member=None):
                assert ensemble_member is not None
                # Note that gmxapi 0.0.7 does not support any consumers of mdrun, so the graph cannot be
                # further widened, and launch can assume that *rank* <= dag width.
                assert dag.graph["width"] == required_width
                assert len(self.infile) == required_width
                assert len(self.runtime_params) == required_width
                if ensemble_member > required_width:
                    raise ApiError(
                        f"Trying to launch MD for simulator {ensemble_member}, "
                        f"but ensemble has width {required_width}."
                    )
                infile = self.infile[ensemble_member]
                runtime_params = self.runtime_params[ensemble_member]

                # Copy and update, if required by `end_time` parameter.
                temp_filename = None
                if "end_time" in runtime_params:
                    # Note that mkstemp returns a file descriptor as the first part of the tuple.
                    # We can make this cleaner in 0.0.7 with a separate node that manages the
                    # altered input.
                    _, temp_filename = tempfile.mkstemp(suffix=".tpr")
                    logger.debug(
                        "Updating input. Using temp file {}".format(temp_filename)
                    )
                    _gmxapi.rewrite_tprfile(
                        source=infile,
                        destination=temp_filename,
                        end_time=runtime_params["end_time"],
                    )
                    tpr_file = temp_filename
                else:
                    tpr_file = infile

                logger.info("Loading TPR file: {}".format(tpr_file))
                system = _gmxapi.from_tpr(tpr_file)
                dag.nodes[name]["system"] = system
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
                dag.nodes[name]["session"] = system.launch(pycontext)
                dag.nodes[name]["close"] = dag.nodes[name]["session"].close

                if "end_time" in runtime_params:

                    def special_close():
                        dag.nodes[name]["session"].close()
                        logger.debug(
                            "Unlinking temporary TPR file {}.".format(temp_filename)
                        )
                        os.unlink(temp_filename)

                    dag.nodes[name]["close"] = special_close
                else:
                    dag.nodes[name]["close"] = dag.nodes[name]["session"].close

                def runner():
                    """Currently we only support a single call to run."""

                    def done():
                        raise StopIteration()

                    # Replace the runner with a stop condition for subsequent passes.
                    dag.nodes[name]["run"] = done
                    return dag.nodes[name]["session"].run()

                dag.nodes[name]["run"] = runner
                return dag.nodes[name]["run"]

            dag.nodes[name]["launch"] = launch

    return Builder(element)


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
        ensemble_rank : numerical index of the simulator task in an ensemble (for participating ranks).
        simulator_rank : numerical index of the current worker in a running simulator (None if not running).
        work_width : Detected number of co-scheduled work elements.
        elements : dictionary of references to elements of the workflow.

    *rank*, *work_width*, and *elements* are empty or None until the work is
    processed, as during session launch.

    To run multiple simulations at a time, whether ensembles or independent
    simulations, Python should be invoked in an MPI environment with `mpi4py`.

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
        ...    member = context.ensemble_rank
        ...    output_path = os.path.join(context.workdir_list[rank], 'traj.trr')
        ...    assert(os.path.exists(output_path))
        ...    print('Worker {} produced {}'.format(rank, output_path))

    Warning:
        Do not run a gmxapi script in an MPI environment without `mpi4py`. gmxapi
        will not be able to determine that other processes are running the same
        script.

    Implementation notes:

    To produce a running session, the Context __enter__() method is called,
    according to the Python context manager protocol. At this time, the attached
    WorkSpec must be feasible  on the available resources. To turn the specified
    work into an executable directed acyclic graph (DAG), handle objects for the
    elements in the work spec are sequenced in dependency-compatible order and
    the context creates a "builder" for each according to the element's operation.
    Each builder is subscribed to the builders of its dependency elements. The
    DAG is then assembled by calling each builder in sequence. A builder can add
    zero, one, or more nodes and edges to the DAG.

    The Session is then launched from the DAG. What happens next is
    implementation-dependent, and it may take a while for us to decide whether
    and how to standardize interfaces for the DAG nodes and edges and/or execution
    protocols. I expect each node will at least have a `launch()` method, but will
    also probably have typed input and output ports as well as some signalling.
    A sophisticated and abstract Session implementation could schedule work only
    to satisfy data dependencies of requested output upon request. Our immediate
    implementation will use the following protocol.

    Each node has a `launch()` method. When the session is entered, the `launch()`
    method is called for each node in dependency order. The launch method returns
    either a callable (`run()` function) or None, raising an exception in
    case of an error. The sequence of callables is stored by the Session. When
    Session.run() is called, the sequence of callables is called in order. If
    StopIteration is raised by the callable, it is removed from the sequence.
    The sequence is processed repeatedly until there are no more callables.

    Note that this does not rigorously handle races or deadlocks, or flexibility
    in automatically chasing dependencies. A more thorough implementation could
    recursively call launch on dependencies (launch could be idempotent or involve
    some signalling to dependents when complete), run calls could be entirely
    event driven, and/or nodes could "publish" output (including just a completion
    message), blocking for acknowledgement before looking for the next set of
    subscribed inputs.
    """

    _api_context: Optional[_gmxapi.Context] = None
    """gmxapi library Context object for managing GROMACS resources.
    
    .. versionchanged:: 0.4
    
        Only processes (MPI ranks) participating in the simulation work will
        acquire an API Context. MPI-enabled GROMACS checks the validity of the
        communicator provided to the Context (and will not issue an unusable
        Context object). It is easier to handle a ``None`` object than to
        develop and use some sort of mock or dummy context.

    """

    __workdir_list: Optional[list] = None
    """Caller-provided working directories (if any).
    
    Support the (deprecated) *workdir_list* key word argument.
    """

    __session_allocation: "Optional[ResourceAllocation[mpi4py.MPI.Comm, Any]]" = None
    """Base resources (if any) from which to assign task resources.
    
    To be removed when we no longer accept *communicator* instead of *resources*.
    """

    _session_resources: "Optional[ResourceAssignment[mpi4py.MPI.Comm, Any]]" = None
    """Session scoped communication resources (where applicable).
    
    When the Context is used as a Python context manager, the state of the
    Context during the `with` block is represented by a Session.
            
    Until better Session abstraction exists at the Python level,
    _session_resources may be added either at Context initialization or during
    session entry. _session_resources will be closed and removed at exit.
    """

    _session_ensemble_communicator: "Optional[mpi4py.MPI.Comm]" = None
    """Cached Session resource for communication between ensemble members.

    Retrieved from _session_resources before a Session for lower overhead
    in `ensemble_update()`.
    """

    _base_rank: Optional[int] = None
    """Cached rank in the parent allocation.

    We release our reference to the resource assignment when exiting the Session,
    but it is convenient to retain the rank that we were on, for debugging purposes.

    This attribute is not necessarily available before the Session is being entered,
    so it is initially None. 
    """

    def __init__(
        self,
        work=None,
        *args,
        resources: "Optional[ResourceAssignment[mpi4py.MPI.Comm, Any]]" = None,
        **kwargs,
    ):
        """Create manager for computing resources.

        Does not initialize resources because Python objects by themselves do
        not have a good way to deinitialize resources. Instead, resources are
        initialized using the Python context manager protocol when sessions are
        entered and exited.

        Appropriate computing resources need to be knowable when the Context is created.

        Keyword Arguments:
            work : work specification with which to initialize this context
            workdir_list (optional): deprecated
            communicator (optional): non-owning reference to a multiprocessing communicator
                from which ensemble and simulation session resources will be derived (deprecated)
            resources (optional): non-owning reference to resources (such as an MPI communicator)
                from which ensemble and simulator resources will be drawn. Caller is responsible
                for providing resources that are compatible with the work load.
                (Use instead of *communicator*.)

        .. versionchanged:: 0.4

            Keyword arguments after *work* are now keyword-only. Positional use
            is deprecated and will be unsupported in a future release for easier
            maintenance of the function signature.

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
        from ..runtime import BaseContext

        # Handle positional arguments until we can require key-word-only.
        # TODO(release-2024): Simplify this logic after a suitable deprecation period to disallow *args.
        if not set(kwargs.keys()).issubset({"workdir_list", "communicator"}):
            raise TypeError("Unsupported key word arguments.")
        if (
            len(args) > 2
            or (len(args) == 2 and len(kwargs) > 0)
            or (len(args) == 1 and "workdir_list" in kwargs)
        ):
            raise TypeError("Unsupported positional arguments.")
        workdir_list = kwargs.get("workdir_list", None)
        communicator: Optional["mpi4py.MPI.Comm"] = kwargs.get("communicator", None)
        if len(args) >= 1:
            warnings.warn(
                "Deprecated use of positional argument. Use explicit key word.",
                DeprecationWarning,
            )
            workdir_list = args[0]
        if len(args) == 2:
            communicator = args[1]
        if communicator is not None:
            warnings.warn(
                "Provide *resources* instead of *communicator*.", DeprecationWarning
            )
        # End: handling for deprecated positional argument handling.

        if resources is not None and communicator is not None:
            raise exceptions.UsageError(
                "Do not provide both *communicator* and *resources*."
            )

        if resources is None:
            if communicator is not None:
                if BaseContext.communicator() is None:
                    BaseContext.set_base_communicator(communicator)
                elif BaseContext.communicator() != communicator:
                    raise exceptions.UsageError(
                        "Default base communicator already configured. "
                        "For non-trivial resource management, "
                        "provide Context with *resources* instead of *communicator*."
                    )
            self.__session_allocation = BaseContext.instance()
        else:
            # We do not intend to take ownership of the resources, but until we can
            # re-write the context manager as a function call (such as a `launch()`
            # method), we should not impose on the caller to maintain a reference to
            # the *resources* if the caller has no need for it and no responsibility
            # to clean it up. Hold a reference at least until the context manager exits.
            self._session_resources = resources

        self.__work = WorkSpec()
        self.__workdir_list = workdir_list

        self._session = None
        self.ensemble_rank = None
        self.simulator_rank = None

        # `work_width` notes the required width of an array of synchronous tasks
        # to perform the specified work.
        # As work elements are processed, self.work_width will be increased as appropriate.
        self.work_width = None

        # initialize the operations map. May be extended during the lifetime of a Context.
        # Note that there may be a difference between built-in operations provided by this module and
        # additional operations registered at run time.
        self.__operations = dict()
        # The map contains a builder for each operation. The builder is created
        # by passing the element to the function in the map. The object returned
        # must have the following methods:
        #
        #   * add_subscriber(another_builder) : allow other builders to subscribe
        #       to this one.
        #   * build(dag) : Fulfill the builder responsibilities by adding an
        #       arbitrary number of nodes and edges to a Graph.
        #
        # The gmxapi namespace of operations should be consistent with a
        # specified universal set of functionalities.
        self.__operations["gmxapi"] = {
            "md": lambda element: _md(self, element),
        }
        # Even if TPR file loading were to become a common and stable enough
        # operation to be specified in an API, it is unlikely to be implemented
        # by any code outside GROMACS, so let's not clutter a potentially more
        # universal namespace.
        self.__operations["gromacs"] = {
            "load_tpr": lambda element: _load_tpr(self, element),
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
        elif hasattr(work, "workspec") and isinstance(work.workspec, WorkSpec):
            workspec = work.workspec
        else:
            raise exceptions.ValueError(
                "work argument must provide a gmx.workflow.WorkSpec."
            )
        workspec._context = self

        # Make sure this context knows how to run the specified work.
        for e in workspec.elements:
            element = WorkElement.deserialize(workspec.elements[e])

            # Note: Non-built-in namespaces (non-native) are treated as modules to import.
            # Native namespaces may not be completely implemented in a particular version of a particular Context.
            if element.namespace in {"gmxapi", "gromacs"}:
                assert element.namespace in self.__operations
                if not element.operation in self.__operations[element.namespace]:
                    # The requested element is a built-in operation but not available in this Context.
                    # element.namespace should be mapped, but not all operations are necessarily implemented.
                    logger.error(
                        "Operation {} not found in map {}".format(
                            element.operation, str(self.__operations)
                        )
                    )
                    # This check should be performed when deciding if the context is appropriate for the work.
                    # If we are just going to use a try/catch block for this test, then we should differentiate
                    # this exception from those raised due to incorrect usage.
                    # The exception thrown here may evolve with https://github.com/kassonlab/gmxapi/issues/125
                    raise exceptions.FeatureNotAvailableError(
                        "Specified work cannot be performed due to unimplemented operation {}.{}.".format(
                            element.namespace, element.operation
                        )
                    )

            else:
                assert element.namespace not in {"gmxapi", "gromacs"}

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
                            "Cannot find implementation for namespace {}. ImportError: {}".format(
                                element.namespace, str(e)
                            )
                        )
                    except AttributeError:
                        raise exceptions.UsageError(
                            "Cannot find factory for operation {}.{}".format(
                                element.namespace, element.operation
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
            if namespace in {"gmxapi", "gromacs"}:
                raise exceptions.UsageError(
                    "Cannot add operations to built-in namespaces."
                )
            else:
                self.__operations[namespace] = dict()
        else:
            assert namespace in self.__operations

        if operation in self.__operations[namespace]:
            raise exceptions.UsageError(
                "Operation {}.{} already defined in this context.".format(
                    namespace, operation
                )
            )
        else:
            self.__operations[namespace][operation] = get_builder

    # Set up a simple ensemble resource
    # This should be implemented for Session, not Context, and use an appropriate subcommunicator
    # that is created and freed as the Session launches and exits.
    def ensemble_update(self, send, recv, tag=None):
        """Implement the ensemble_update member function that gmxapi through 0.0.6 expects.

        This is a draft of a Context feature that may not be available in all
        Context implementations.
        This feature requires that the Context is able to provide the
        ensemble_communicator feature and the numpy feature.
        """
        # gmxapi through 0.0.6 expects to bind to this member function during "build".
        # This behavior needs to be deprecated (bind during launch, instead), but this
        # dispatching function should be an effective placeholder.

        if tag is None or str(tag) == "":
            raise exceptions.ApiError("ensemble_update must be called with a name tag.")

        ensemble_communicator = self._session_ensemble_communicator

        assert not tag is None
        assert str(tag) != ""
        if tag not in self.part:
            self.part[tag] = 0
        logger.debug("Performing ensemble update.")
        ensemble_communicator.Allreduce(send, recv)
        buffer = np.array(recv, copy=False)
        buffer /= self.work_width
        suffix = "_{}.npz".format(tag)
        # These will end up in the working directory and each ensemble member will have one
        filename = str(
            "rank{}part{:04d}{}".format(
                ensemble_communicator.Get_rank(), int(self.part[tag]), suffix
            )
        )
        np.savez(filename, recv=recv)
        self.part[tag] += 1

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
        from ..runtime import ResourceAssignment, ResourceRequirements

        try:
            import networkx as nx
            from networkx import DiGraph as _Graph
        except ImportError:
            raise exceptions.FeatureNotAvailableError(
                "gmx requires the networkx package to execute work graphs."
            )

        logger.debug("Launching session from {}".format(self.__initial_cwd))

        if self._session is not None:
            raise exceptions.Error("Already running.")
        if self.work is None:
            raise exceptions.UsageError("No work to perform!")

        # Set up the global and local context.
        # Check the global MPI configuration
        logger.debug(
            f'Starting session for GROMACS MPI type "{gmxapi.utility.config()["gmx_mpi_type"]}"'
        )

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
                new_builder = self.__operations[element.namespace][element.operation](
                    element
                )
                assert hasattr(new_builder, "add_subscriber")
                assert hasattr(new_builder, "build")

                logger.info("Collected builder for {}".format(element.name))
            except LookupError as e:
                request = ".".join([element.namespace, element.operation])
                message = "Could not find an implementation for the specified operation: {}. ".format(
                    request
                )
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
                # and pick the element corresponding to the local rank.
                if isinstance(dependency, (list, tuple)):
                    if self._session_resources is None:
                        raise exceptions.FeatureNotAvailableError(
                            "Cannot infer dependency graph topology. Provide a ResourceAssignment instead "
                            "of a raw Communicator for complex dependencies."
                        )
                    _member_id = self._session_resources.subtask_id()
                    if _member_id is not None:
                        assert len(dependency) > _member_id
                        name = str(dependency[_member_id])
                        del _member_id
                    else:
                        name = None
                else:
                    name = dependency
                if name is not None:
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
        self.work_width = graph.graph["width"]

        # Support the legacy invocation (received *communicator* instead of *resources*)
        if self.__session_allocation is not None:
            assert self._session_resources is None
            # We will need to own this resource assignment, being sure to close
            # and release in __exit__.
            _ensemble_label = 0
            requirements = ResourceRequirements(
                communication=tuple(
                    {"ensemble_comm": _ensemble_label, "subtask_comm": True}
                    for _ in range(self.work_width)
                )
            )
            self._session_resources = ResourceAssignment.assign(
                allocation=self.__session_allocation, requirements=requirements
            )

        self.ensemble_rank = self._session_resources.subtask_id()
        self.simulator_rank = self._session_resources.subtask_rank()
        self._base_rank = self._session_resources.base_rank()
        member_id = self.ensemble_rank

        assert self._session_resources is not None

        ensemble_comm = self._session_resources.ensemble_communicator()
        if ensemble_comm:
            _ensemble_size = ensemble_comm.Get_size()
            if self.work_width != _ensemble_size:
                raise ProtocolError(
                    f"Resource mismatch: {self.work_width} subtasks received ensemble resources for "
                    f"{_ensemble_size}."
                )
            assert ensemble_comm.Get_rank() == self.ensemble_rank
            self._session_ensemble_communicator = ensemble_comm

        # Provide a Session-scoped gmxapi._gmxapi.Context.
        # Provide an alternative to MPI_COMM_WORLD for the SimulatorCommunicator.
        # libgromacs does not accept a null communicator at run time, but our calling
        # convention is to make the same API calls on all ranks at the WORLD level
        # (or at least at the Context.__communicator level).
        # self._api_context needs to be in place before the launch() functor from
        # the _md Builder is called.
        simulation_communicator = self._session_resources.communicator()

        if simulation_communicator:
            # Prepare working directories.
            # This should probably be moved to some aspect of the Session and either
            # removed from here or made more explicit to the user.
            workdir_list = self.__workdir_list
            if workdir_list is None:
                workdir_list = [
                    os.path.join(".", str(i)) for i in range(self.work_width)
                ]
            workdir_list = [os.path.abspath(directory) for directory in workdir_list]

            self.workdir = simulation_communicator.bcast(workdir_list[member_id])
            del workdir_list

            if self._session_resources.is_subtask_root() and not os.path.exists(
                self.workdir
            ):
                os.mkdir(self.workdir)

            try:
                self._api_context = _gmxapi.create_context(simulation_communicator)
            except Exception as e:
                logger.error(
                    "Got exception when trying to create default library context.",
                    exc_info=e,
                )
                # Warning: Be careful about suppressing exceptions.
                #
                # Avoid collective calls in `finally` blocks or context manager
                # `__exit__` functions!
                #
                # Suppressing exceptions could cause hangs on errors.
                # There is potential divergence in collective MPI calls
                # if all ranks do not proceed as expected.
                # If another rank may be in (or approaching) a collective call,
                # we must exit with an uncaught exception (or explicitly
                # MPI_ABORT). If an exception may cause divergence between ranks,
                # we must resolve it (or allow MPI_ABORT) before entering another
                # collective call.
                #
                # Here, the exception raised is likely more
                # informative than the AssertionError that would occur below,
                # so re-raise.
                #
                # For the present context manager, the collective calls in __exit__
                # are okay because __exit__ will not be called if __enter__ does
                # not complete successfully.
                raise e

        if simulation_communicator:
            assert self._api_context is not None
        else:
            assert self._api_context is None

        # In the first implementation, we only provide the ensemble
        # communicator to the simulation master ranks (as we do for tMPI).
        # In the long term, it may be appropriate for all simulator ranks to
        # have access to the ensemble communicator, but that would be tricky
        # to get right. This requires that the work load manager (the current
        # code block) needs to know (or be able to defer negotiation of)
        # something about the library-internal resource allocation.
        # (See also https://gitlab.com/gromacs/gromacs/-/issues/3718)
        #
        # Also note: The GROMACS core does not have a way to receive a
        # communicator that is consistent with the notion of the "ensemble communicator".
        # The `_session_ensemble_communicator` is currently only a gmxapi-specific
        # Session resource usable by MD extension code based on sample_restraint.
        # We should reconsider this in the context of current and future use of
        # the multisim communicator.
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
        if simulation_communicator:
            # print(graph)
            logger.debug(("Launching graph {}.".format(graph.graph)))
            logger.debug("Graph nodes: {}".format(str(list(graph.nodes))))
            logger.debug("Graph edges: {}".format(str(list(graph.edges))))

            simulator_rank = simulation_communicator.Get_rank()
            logger.info(
                "Launching work for member {}, subcommunicator rank {}.".format(
                    member_id, simulator_rank
                )
            )

            # Launch the work for this rank
            simulation_communicator.barrier()
            if not os.path.isdir(self.workdir):
                raise exceptions.Error(
                    "{} is not a valid working directory.".format(self.workdir)
                )
            os.chdir(self.workdir)
            logger.info(
                "rank {} changed directory to {}".format(member_id, self.workdir)
            )
            sorted_nodes = nx.topological_sort(graph)
            runners = []
            closers = []
            for name in sorted_nodes:
                launcher = graph.nodes[name]["launch"]
                # TODO(#4079): launch is explicitly a collaboration with the Context
                #  to build a Session, which should provide SessionResources here.
                runner = launcher(member_id)
                if runner is not None:
                    runners.append(runner)
                    closers.append(graph.nodes[name]["close"])

            # Get a session object to return. It must simply provide a `run()` function.
            class Session(object):
                _context: _gmxapi.Context = None

                def __init__(self, runners, closers, library_context: _gmxapi.Context):
                    self.runners = list(runners)
                    self.closers = list(closers)
                    if library_context is None:
                        raise ProtocolError(
                            "Simulation Session requires a Context from the C++ extension module."
                        )
                    self._context = library_context

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
                    del self._context

            self._session = Session(runners, closers, self._api_context)
        else:
            logger.warning("Process {} has no work to do".format(os.getpid()))

            class NullSession:
                def __init__(self, base_rank: int):
                    self.rank = base_rank

                def run(self):
                    logger.info(f"Running null session on rank {self.rank}.")
                    return True

                def close(self):
                    logger.info("Closing null session.")

            assert self._api_context is None
            self._session = NullSession(self._session_resources.base_rank())

        # Make sure session has started on all ranks before continuing?

        self._session.graph = graph
        return self._session

    def __exit__(self, exception_type, exception_value, traceback):
        """Implement Python context manager protocol."""
        try:
            logger.info(
                "Exiting session on rank {}.".format(
                    self._session_resources.base_rank()
                )
            )
            if self._session is not None:
                logger.info("Calling session.close().")
                self._session.close()
            else:
                # Note: we should not have a None session but rather an API-compliant Session that just has no work.
                # Reference: https://github.com/kassonlab/gmxapi/issues/41
                logger.info("No _session known to context or session already closed.")

            if self._api_context is not None:
                # It really seems like we should expect an explicit close() method in the future.
                if hasattr(self._api_context, "close"):
                    try:
                        self._api_context.close()
                    except _gmxapi.Exception:
                        logger.exception("Error while closing GROMACS API context.")
                logger.debug(f"De-referencing library context {self._api_context}")
                self._api_context = None

            # We only need to free resources if we allocated them here because the caller
            # did not provide session resources.
            if self.__session_allocation is not None:
                # Support the legacy invocation (*communicator* instead of *resources*)
                if self._session_resources is not None:
                    logger.debug("Freeing session resources.")
                    self._session_resources.close()
        finally:
            self._session_resources = None
            if self.__session_allocation is not None:
                logger.debug("Releasing session allocation.")
                self.__session_allocation = None
            os.chdir(self.__initial_cwd)

        self._session = None
        self._session_ensemble_communicator = None
        self.simulator_rank = None
        # Leave `ensemble_rank` set to help clients find the right results afterwards.
        # I.e. do not set `self.ensemble_rank = None` here.

        logger.info("Session closed on context rank {}.".format(self._base_rank))
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
        elif hasattr(work, "workspec") and isinstance(work.workspec, WorkSpec):
            workspec = work.workspec
        else:
            raise exceptions.ValueError(
                "work argument must provide a gmx.workflow.WorkSpec."
            )
    if (
        workspec is not None
        and hasattr(workspec, "_context")
        and workspec._context is not None
    ):
        context = workspec._context
    else:
        context = Context(work=workspec)

    return context
