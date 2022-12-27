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

"""Manage computing resources and runtime context.

Utilities, context managers, and singletons for handling resource allocation
and lifetime management.

For the purposes of this module, the term "assignment" refers to resources that
have been reserved (usually for exclusive use) within a well-defined scope
(such as a specific function call or phase of program execution).

"Allocation" can be taken to mean "resources that are available to be assigned,"
and may or may not be tightly scoped.

Resource assignment may be nested, in which case an Assignment with a broader
scope in the program is used as the Allocation from which a more local
Assignment is made.

Assignment and Allocation are *roles* that can be served by Context objects.
Abstractly, a Context is a notion of some aspect of program or resource *state*.
A Context may be represented by a concrete class when an object manages the
details needed to participate in a stateful protocol.

Note that resource assignment is usually a collective operation within the
scope of an allocation. For instance, an ``MPI_Comm_split`` call must be made
(synchronously) on all members of the parent MPI communicator. Similar care must
be taken with the Allocation and Assignment protocols in this module.

.. versionadded:: 0.4.0

    This module generalizes some of the resource management previously in
    :py:mod:`gmxapi.simulation.mdrun`, decoupling MPI communicator handling
    from :py:mod:`gmxapi.simulation.context`.
    See also https://gitlab.com/gromacs/gromacs/-/issues/3718

"""

__all__ = (
    "scoped_resources",
    "ResourceAllocation",
    "ResourceAssignment",
    "ResourceRequirements",
)

import dataclasses
import threading
import typing
import warnings
import weakref
from _weakref import ReferenceType
from contextlib import contextmanager

import gmxapi.utility
from gmxapi import exceptions
from gmxapi import logger as root_logger

try:
    from typing import Protocol
except ImportError:
    # typing.Protocol provides helpful static type checking support for
    # interfaces, but is not availbale until Python 3.8.
    # noinspection PyTypeChecker
    Protocol = typing.Generic

# Initialize module-level logger
logger = root_logger.getChild(__name__)
logger.info("Importing {}".format(__name__))

# gmxapi requires mpi4py, but it is useful to allow module-level imports to fail
# so that documentation builds can succeed without package dependencies. It is
# also helpful to retain the accompanying logic in case we want to decouple from
# mpi4py at some point in the future.
try:
    import mpi4py.MPI as _MPI
    from mpi4py.MPI import Comm as Communicator
    from mpi4py.MPI import COMM_NULL as MPI_COMM_NULL
    from mpi4py.MPI import COMM_WORLD as MPI_COMM_WORLD
    from mpi4py.MPI import UNDEFINED as MPI_UNDEFINED
except (ImportError, ModuleNotFoundError):
    _MPI = None
    Communicator = None

_context_lock = threading.Lock()
"""Use to avoid race conditions when changing singleton state."""

CommunicatorT = typing.TypeVar("CommunicatorT")
ProvidesTokenT = typing.TypeVar("ProvidesTokenT")
ReceivedTokenT = typing.TypeVar("ReceivedTokenT")


class ResourceAllocation(Protocol[CommunicatorT, ProvidesTokenT]):
    """The interface we expect for a ResourceAllocation role."""

    def communicator(self) -> typing.Optional[CommunicatorT]:
        """Get a reference to the communicator allocated for the current task scope, if any."""
        raise NotImplementedError

    def subscribe(self, subscriber) -> ProvidesTokenT:
        """Alert the ResourceAllocation that its resources are in use.

        For initial implementation, we cannot support locked resources with
        an ``acquire()`` / ``release()`` protocol. But we can still support
        some bookkeeping and issue warnings for developers.
        """
        raise NotImplementedError

    def unsubscribe(self, token: ProvidesTokenT):
        """Notify the ResourceAllocation that a subscriber is done using resources."""
        raise NotImplementedError


def _finalize_communicator(comm):
    if _MPI is not None:
        if comm == MPI_COMM_NULL:
            return
    # noinspection PyBroadException
    try:
        comm.Free()
    except Exception:
        logger.exception(f"Could not Free {comm}.")


class ResourceAssignment(
    ResourceAllocation[CommunicatorT, int],
    typing.Generic[CommunicatorT, ReceivedTokenT],
):
    """Container for assigned resources.

    If an MPI communicator is provided when initializing the ResourceAssignment,
    ownership of the communicator is **not** assumed, and `comm.Free()` will
    **not** be called automatically when the ResourceAssignment is closed or
    finalized.

    To transfer ownership of the communicator in *task* or *ensemble*, add
    a finalizer using either the *finalizer* key word argument or the
    `add_finalizer` method.
    """

    parent: "ReferenceType[ResourceAllocation[CommunicatorT, ReceivedTokenT]]"
    """Resource allocation from which the current assignment is derived.
    
    Notes that this is a weak reference to the provided ResourceAllocation.
    See :py:mod:`weakref`.
    """

    _resource_id: typing.Optional = None
    """Global identifier for this assignment.
    
    If multiple threads or processes share the assignment, the identifier
    is shared across the threads and processes.

    Across threads or processes, multiple instances of ResourceAssignment
    could refer to the same (multi-core) assignment. This identifier can be
    used to identify callers sharing a resource assignment for a particular
    scope of execution. Callers may be grouped by task (within an ensemble)
    using *_group_id*, and within a task using *_local_id* (such as to
    identify the "root" participant in a task).
    
    This attribute is experimental. If used, it is the callers responsibility
    to set the value uniformly.
    """

    _group_ids: typing.Optional[typing.Tuple[int]] = None
    """The resource subgroup within *_resource_id*.

    This identifier is not generally unique when *_resource_id* varies.
    The identifier _may_ be set even for non-participating callers to indicate
    that resources owned by the caller are being used for the task associated
    with the group.
    
    This attribute is experimental. If used, it is the callers responsibility
    to set the value appropriately.    
    """

    _roots: typing.Optional[typing.Tuple[typing.Union[int, None]]] = None
    """Ensemble rank identifiers for subtask roots.
    
    Generally equal to the _group_ids for participants that are subtask roots
    and None otherwise.
    
    This information could be stored in a more concise structure, but a little
    redundancy is useful for error checking and it is convenient to have the
    same look-up scheme of ``attr[_base_comm_rank]``.
    
    This attribute is experimental. If used, it is the callers responsibility
    to set the value appropriately.    
    """

    _subtask_ranks: typing.Optional[typing.Tuple[typing.Union[int, None]]] = None
    """Instance identifier for this assignment.
    
    Generally, serves as the "participant ID" for a resource assignment
    comprising multiple threads or processes.
    
    This identifier is not generally unique when *_resource_id* varies or
    when *_group_id* varies. The identifier should be `None` for
    non-participating callers.

    This attribute is experimental. If used, it is the callers responsibility
    to set the value with appropriate uniqueness.
    """

    _token: ReceivedTokenT
    _subscribers: list

    __base_comm_rank: typing.Optional[int] = None
    """Rank of the caller in the parent allocation."""

    __communicator: typing.Optional[CommunicatorT]
    __ensemble: typing.Optional[CommunicatorT]
    __finalizers: list
    __next_monotonic_integer = 0

    def __init__(
        self,
        *,
        subtask: CommunicatorT,
        parent: ResourceAllocation[CommunicatorT, ReceivedTokenT],
        ensemble: CommunicatorT = None,
        finalizers: typing.Sequence[weakref.finalize] = (),
    ):
        try:
            self.__base_comm_rank = parent.communicator().Get_rank()
            self._base_comm_size = parent.communicator().Get_size()
            self.parent = weakref.ref(parent)
        except (AttributeError, TypeError) as e:
            raise ValueError(f"Unusable resources: {parent}.") from e
        self.__communicator = subtask
        self.__ensemble = ensemble
        try:
            self._token = parent.subscribe(self)
            logger.debug(f"{self} subscribed to {parent}. Token: {self._token}")
        except AttributeError:
            self._token = None
        self._subscribers = list()
        self.__finalize_dirty = weakref.finalize(
            self,
            warnings.warn,
            message=f"{repr(self)} finalized without an explicit call to close().",
            category=ResourceWarning,
        )
        self.__finalizers = [self.__finalize_dirty]
        if self._token is not None:
            self.__finalizers.append(
                weakref.finalize(self, parent.unsubscribe, self._token)
            )
        self.__finalizers.append(
            weakref.finalize(
                self, _remaining_subscribers_warning, repr(self), self._subscribers
            )
        )
        self.__finalizers.extend(finalizers)

    def __repr__(self):
        return (
            f"<ResourceRequirements:{self.__base_comm_rank} _resource_id:{self._resource_id} _group_ids:"
            f"{self._group_ids} _roots:{self._roots} _subtask_ranks:{self._subtask_ranks} parent:"
            f"{self.parent}:{self._token} communicator:{self.__communicator} ensemble:{self.__ensemble}>"
        )

    def _next_monotonic_integer(self):
        i = self.__next_monotonic_integer
        self.__next_monotonic_integer = i + 1
        return i

    def subscribe(self, subscriber) -> int:
        """Receive a subscription notification.

        Implements ResourceAllocation.
        """
        _id = id(subscriber)
        _repr = repr(subscriber)
        token = self._next_monotonic_integer()
        self._subscribers[token] = (_id, _repr)
        return token

    def unsubscribe(self, token: int):
        """Receive an unsubscription notification.

        Implements ResourceAllocation.
        """
        del self._subscribers[token]

    def communicator(self) -> CommunicatorT:
        """Get the communicator for the current process.

        In an MPI context, other Python interpreter processes participating in
        the current work may have different communicators or null communicators.
        """
        return self.__communicator

    def parent_communicator(self) -> typing.Optional[CommunicatorT]:
        """Get a reference to the communicator from which these resources are derived."""
        parent = self.parent()
        if parent:
            return parent.communicator()
        else:
            return None

    def ensemble_communicator(self) -> typing.Optional[CommunicatorT]:
        """Get a reference to the communicator that connects the members of an ensemble task."""
        return self.__ensemble

    def close(self):
        if self.closed:
            logger.debug(f"{self}.close() called on an already-closed resource.")
        else:
            self.__finalize_dirty.detach()
            for finalizer in reversed(self.__finalizers):
                finalizer()
            del self.__finalizers
            logger.debug(f"Closed {self}.")

    def add_finalizer(self, func, *args, **kwargs):
        if self.closed:
            raise exceptions.UsageError(
                "Cannot add finalizer to a resource that is already closed."
            )
        self.__finalizers.append(weakref.finalize(self, func, *args, **kwargs))

    @property
    def closed(self):
        return not self.__finalize_dirty.alive

    @classmethod
    def assign(
        cls,
        allocation: ResourceAllocation[CommunicatorT, ReceivedTokenT],
        requirements: "ResourceRequirements",
    ):
        """Assign resources to members of a new ensemble task.

        If *allocation* includes an MPI communicator, the resources yielded will
        contain a sub-communicator of *allocation.communicator()*.

        Args:
            allocation: available resources, or *None* to determine automatically.
            requirements: parameters for the resource assignment, such as MPI communicators.

        *allocation* is required if *requirements* is not *None*. *requirements*
        may be validated in terms of the provided *allocation*.

        When GROMACS is not built with an MPI library, the number of ranks used
        is equal to *len(requirements.communication)*. Participating ranks are
        not necessarily consecutive.

        Specifically, ``ResourceAllocation.assign(...).communicator`` will be
        :py:class:`~mpi4py.MPI.COMM_NULL` in some processes
        (consistent with :py:func:`mpi4py.MPI.Comm.Split()`) when
        - ``mpi4py`` is in use (``gmxapi.version.has_feature('mpi_bindings') == True``
          or *allocation.communicator()* is an :py:class:`mpi4py.MPI.Comm`),
        - ``gmx.utility.config()['gmx_mpi_type'] != 'library'``,
        - *requirements* is not ``None``, and
        - ``len(requirements.communication) < allocation.communicator().Get_size()``.

        Callers should inspect the yielded resources to decide whether to
        use the new communicator.

        For MPI-enabled GROMACS simulators
        (``gmx.utility.config()['gmx_mpi_type'] == 'library'``), all available
        ranks in *allocation.communicator()* are distributed as evenly as possible
        to the new subcommunicator(s) in contiguous blocks.
        """
        if requirements is not None:
            if allocation is None:
                raise exceptions.UsageError(
                    "Provide allocated resources from which to request a sub-allocation."
                )

        if allocation is None:
            # We need to think more about race conditions and other structural
            # details before allowing implicit automatic detection of resources.
            raise exceptions.FeatureNotAvailableError(
                "Implicit allocation is not currently available."
            )
        base_communicator = allocation.communicator()
        assert isinstance(base_communicator, Communicator)
        base_comm_size = base_communicator.Get_size()
        base_comm_rank = base_communicator.Get_rank()

        n_subtasks = len(requirements.communication)
        subtask_comm = None
        ensemble_comm = None
        _resource_id = None
        group_ids = None
        _local_id = None
        local_ids = None
        roots = None
        if n_subtasks > 0:
            if n_subtasks > base_comm_size:
                error = f"Cannot satisfy {requirements} from {allocation}."
                logger.error(
                    error
                    + " "
                    + "To get a scoped communicator of a specific size, provide a sufficient "
                    + "communicator supporting the mpi4py.MPI.Comm interface."
                )
                raise exceptions.UsageError(error)

            group_ids = tuple(
                _rank * n_subtasks // base_comm_size for _rank in range(base_comm_size)
            )
            """Label a number of contiguous blocks of processes in an allocation.
    
            Divide the ranks of *allocation* into *count* consecutive blocks of
            (roughly) equal size, labeled from 0 to (*count* - 1). Get a list
            corresponding to the label for each rank.
        
            Warning:
                Blocks may span node boundaries.
                When allocation spans multiple nodes, but blocks are smaller than one node,
                we should try to align blocks at node boundaries for more uniform
                communication latency. This is not currently possible (but behavior
                should be kept consistent with `roots`).
            """

            assert set(group_ids) == set(range(n_subtasks))
            _group_id = group_ids[base_comm_rank]

            # Get the index of the first occurrence of each group ID.
            root_indices = set(group_ids.index(i) for i in range(n_subtasks))
            roots = tuple(
                group_ids[_rank] if _rank in root_indices else None
                for _rank in range(base_comm_size)
            )
            """Sparse assignment for each rank in allocation.
            
            Integers from 0 to (*count* - 1) will be distributed in a sequence with
            as many elements as the communicator size in *allocation*. Only one rank
            from the communicator is associated with each ID. Other ranks are
            labeled `None`.
        
            Warning:
                Blocks may span node boundaries.
                When allocation spans multiple nodes, we should try to align blocks at
                node boundaries so that we can assume that a rank assigned None could
                lend its CPU core to the closest non-None member. This is not currently
                possible. If pursued in the future, it could be appropriate to create
                an additional resource layer comprising a parent MPI.Group for the
                processes "owning" each subset of cores.
            """

            assert len(set(roots) - {None}) == n_subtasks
            assert set(roots) - {None} == set(range(n_subtasks))

            _current_offset = 0
            subtask_ranks = [0] * base_comm_size
            for _rank in range(base_comm_size):
                if roots[_rank] is not None:
                    _current_offset = _rank
                subtask_ranks[_rank] = _rank - _current_offset

            communication_requirements = tuple(
                key for task in requirements.communication for key in task
            )
            for requirement in communication_requirements:
                if requirement == "subtask_comm":
                    # We don't support any options, but the participants must agree.
                    assert all(
                        subtask.get("subtask_comm", False)
                        for subtask in requirements.communication
                    )

                    if gmxapi.utility.config()["gmx_mpi_type"] == "library":
                        color = tuple(subtask_id for subtask_id in group_ids)
                    else:
                        color = tuple(
                            MPI_UNDEFINED if roots[i] is None else subtask_id
                            for i, subtask_id in enumerate(group_ids)
                        )
                    subtask_comm = base_communicator.Split(
                        color[base_comm_rank], key=subtask_ranks[base_comm_rank]
                    )
                    if subtask_comm:
                        _local_id = subtask_comm.Get_rank()
                        # Check expected behavior for current resource allocation scheme.
                        expected_comm_size = color.count(color[base_comm_rank])
                        actual_comm_size = subtask_comm.Get_size()
                        if gmxapi.utility.config()["gmx_mpi_type"] != "library":
                            assert expected_comm_size == 1
                        assert actual_comm_size == expected_comm_size
                        assert _local_id == subtask_ranks[base_comm_rank]
                elif requirement == "ensemble_comm":
                    # We are currently assuming that all task requirements include
                    # a specifier for the same ensemble.
                    if not all(
                        "ensemble_comm" in task for task in requirements.communication
                    ):
                        raise exceptions.ProtocolError(
                            "All task requirements are assumed to specify the same ensemble."
                        )
                    # We allow the caller to provide arbitrary ensemble labels, but we
                    # do not yet track them for uniqueness. For now, just try to use the
                    # label to produce a reasonably unique MPI_INT.
                    ensemble_labels = set(
                        task["ensemble_comm"] for task in requirements.communication
                    )
                    assert len(ensemble_labels) == 1
                    _resource_id = ensemble_labels.pop()
                    logger.debug(f"Resource ID: {_resource_id}")
                    color = tuple(
                        MPI_UNDEFINED if _id is None else _resource_id for _id in roots
                    )

                    ensemble_comm = base_communicator.Split(color[base_comm_rank])
                else:
                    raise ValueError(
                        f"Unrecognized communication requirement: {requirement}"
                    )

            # One collective MPI call serves as a useful synchronization and error
            # checking point, as well as a convenient way to get the following data.
            local_ids = base_communicator.allgather(_local_id)
            assert len(local_ids) == base_comm_size
            assert local_ids[base_comm_rank] == _local_id

        resources = cls(subtask=subtask_comm, parent=allocation, ensemble=ensemble_comm)
        assert base_comm_rank == resources.base_rank()
        # Note: mpi4py.MPI.Comm.__bool__() compares to MPI_COMM_NULL.
        if subtask_comm:
            resources.add_finalizer(
                weakref.finalize(resources, _finalize_communicator, subtask_comm)
            )
        if ensemble_comm:
            resources.add_finalizer(
                weakref.finalize(resources, _finalize_communicator, ensemble_comm)
            )
        if local_ids is None:
            resources._subtask_ranks = tuple([None] * base_comm_size)
        else:
            resources._subtask_ranks = tuple(local_ids)
        resources._group_ids = group_ids
        resources._resource_id = _resource_id
        resources._roots = roots

        return resources

    def subtask_rank(self):
        if self.__communicator:
            return self._subtask_ranks[self.__base_comm_rank]

    def subtask_id(self):
        """Identifier for the (group of) resource participants.

        Within the scope of a `resource_id`, callers participating in the same subtask can be identified by sharing the same subtask ID.

        Returns None for callers that are not participating.
        """
        if self.__communicator:
            return self._group_ids[self.__base_comm_rank]

    def is_subtask_root(self):
        """Whether the participant is the "root" member of the subtask.

        If True, then `subtask_id` is also the rank in the `ensemble_communicator`, if any.
        """
        if self._roots is not None:
            return self._roots[self.__base_comm_rank] is not None

    def base_rank(self):
        """Identify the caller within the parent allocation."""
        return self.__base_comm_rank


def _remaining_subscribers_warning(provider: str, subscribers: list):
    for token, subscriber in subscribers:
        _id, _repr = subscriber
        warnings.warn(
            f"Closing {provider} while {_repr} (id: {_id}, token: {token}) is still subscribed."
        )


class BaseContext:
    """Resource allocation and options.

    Includes default run time execution parameters that may be configured only
    until resources are assigned. That is, class data is only mutable while
    no instance exists.
    """

    __instance: typing.ClassVar[typing.Optional[weakref.ReferenceType]] = None

    __mpi_communicator: typing.ClassVar[typing.Optional[Communicator]] = None

    def __new__(cls, *args, **kwargs):
        # Check the weakref for a singleton instance, if any.
        if cls.__instance is None:
            instance = None
        else:
            instance = cls.__instance()
        # Set a new singleton instance, if necessary.
        if instance is None:
            instance = super(BaseContext, cls).__new__(cls, *args, **kwargs)
            cls.__instance = weakref.ref(instance)
        return instance

    @classmethod
    def communicator(cls) -> Communicator:
        return cls.__mpi_communicator

    @classmethod
    def set_base_communicator(cls, communicator):
        if cls.__instance is not None:
            raise exceptions.UsageError(
                "Cannot edit resource defaults after initialization."
            )
        cls.__mpi_communicator = communicator

    @classmethod
    def instance(cls):
        """Get the singleton instance. (Finalize and instantiate if needed.)"""
        if _MPI is not None and cls.__mpi_communicator is None:
            cls.set_base_communicator(MPI_COMM_WORLD)
        return cls()

    def subscribe(self, subscriber):
        """Accept a subscription notification.

        Supports ResourceAllocation protocol, but has no effect.
        """
        return None

    def unsubscribe(self, token):
        """Accept an un-subscription notification.

        Supports ResourceAllocation protocol, but has no effect.
        """
        return None


# TypedDict is not available until Python 3.8, or by adding `typing-extensions` as a project dependency.
# class _CommunicationRequirements(typing.TypedDict, total=False):
#     ensemble_comm: int
#     subtask_comm: bool
_CommunicationRequirements = dict


@dataclasses.dataclass
class ResourceRequirements:
    """Named requirements for resources to be assigned for a task."""

    communication: typing.Tuple[_CommunicationRequirements, ...] = dataclasses.field(
        default_factory=tuple
    )
    """Requested communicator details.
    
    This field is a container with one element for each block of resources.
    Further details of the container or element type are not yet specified.
    """

    def __post_init__(self):
        for requirement in self.communication:
            if requirement.get("ensemble_comm", 0) // int(2**16) > 0:
                raise ValueError(f"Ensemble communicator labels must fit into uint_16.")


@contextmanager
def scoped_resources(
    allocation: ResourceAllocation[Communicator, typing.Any],
    requirements: ResourceRequirements,
):
    """Manage computing resources with a lifecycle.

    Get a Python context manager for a new resource assignment that is released
    on exiting the `with` block. (Wraps `assign_ensemble()`).

    Args:
        allocation: available resources, or *None* to determine automatically.
        requirements: parameters for the resource assignment, such as MPI communicators.

    *allocation* is required if *requirements* is not *None*. *requirements*
    may be validated in terms of the provided *allocation*.

    .. versionchanged:: 0.4

        Previous versions transformed :py:class:`~mpi4py.MPI.COMM_NULL` into a mock
        communicator before returning for consistency in trivial cases with or
        without :py:mod:`mpi4py`. This translation adds unnecessary indirection,
        and has been removed.

        Similarly, for *original_comm.Get_size()* or *requested_size* of 1,
        an unspecified "fake" communicator could be produced to avoid an
        ``MPI_Comm_split()``. While more efficient, the code complexity was
        problematic and has been removed.

    TODO(#4079): Acquire Task computing resources in a more limited scope than data resources.

        This context manager is called "scoped_resources" instead of "scoped_assignment"
        to avoid implying there is any locking on the usage of the returned resources.
        gmxapi.operation.ResourceManager has a hard-coded protocol for preparing task
        resources that first recursively resolves inputs, then sets up the task runner.
        gmxapi.simulation.mdrun handles communication resources in its input adapter,
        LegacyImplementationSubscription, which ends up being reentrant for chained
        simulations, calling ResourceManager.local_input() _after_ sub-allocating
        Resources. This is a bad design, and it means that there is the appearance
        of resource oversubscription because multiple handles exist to sub-allocations
        from the same base resources. However, there are legitimate reasons for
        components to use multiple communicators on the same ranks for different
        groups of communicating participants, and it is not clear how best to distinguish
        between this sort of usage, either.

    """
    resources = ResourceAssignment.assign(
        allocation=allocation, requirements=requirements
    )

    try:
        yield resources
    finally:
        resources.close()
