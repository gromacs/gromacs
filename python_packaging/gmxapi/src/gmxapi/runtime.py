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
import dataclasses
import threading
import typing
import weakref
from contextlib import contextmanager

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
logger.info('Importing {}'.format(__name__))


try:
    import mpi4py.MPI as _MPI
    from mpi4py.MPI import Comm as Communicator
except (ImportError, ModuleNotFoundError):
    _MPI = None
    Communicator = None


_context_lock = threading.Lock()
"""Use to avoid race conditions when changing singleton state."""


CommunicatorT = typing.TypeVar('CommunicatorT')


class ResourceAllocation(Protocol[CommunicatorT]):
    """The interface we expect for a ResourceAllocation role."""
    def communicator(self) -> CommunicatorT:
        ...


def _finalize_communicator(comm):
    if _MPI is not None:
        if comm == _MPI.COMM_NULL:
            return
    # noinspection PyBroadException
    try:
        comm.Free()
    except Exception:
        logger.exception(f'Could not Free {comm}.')


class ResourceAssignment(ResourceAllocation[CommunicatorT]):
    """Container for assigned resources.

    If an MPI communicator is provided when initializing the ResourceAssignment,
    ownership of the communicator is assumed, and `comm.Free()` will be called
    when the ResourceAssignment is closed or finalized.
    """
    __communicator: typing.Optional[CommunicatorT]
    __finalizers: list

    def __init__(self, _mpi_comm: CommunicatorT = None):
        self.__communicator = _mpi_comm
        self.__finalizers = [weakref.finalize(self, logger.debug, f'Closing {self}.')]
        if hasattr(_mpi_comm, 'Free'):
            self.__finalizers.append(weakref.finalize(self, _finalize_communicator, _mpi_comm))

    def communicator(self) -> CommunicatorT:
        return self.__communicator

    def borrow_communicator(self, communicator: Communicator):
        """Borrow a communicator from the caller.

        Caller is responsible for freeing the provided communicator.

        To transfer ownership of the communicator and let it be freed by the
        ResourceAssignment, provide the communicator when initializing the object
        (or add a finalizer to the ResourceAssignment to free the communicator
        on `close()`).
        """
        if self.__communicator is not None:
            raise exceptions.UsageError('A communicator has already been assigned.')
        self.__communicator = communicator

    def close(self):
        if self.closed:
            logger.debug(f'{self}.close() called on an already-closed resource.')
        else:
            for finalizer in reversed(self.__finalizers):
                finalizer()
            self.__finalizers = []

    def add_finalizer(self, func, *args, **kwargs):
        if self.closed:
            raise exceptions.UsageError('Cannot add finalizer to a resource that is already closed.')
        self.__finalizers.append(weakref.finalize(self, func, *args, **kwargs))

    @property
    def closed(self):
        return len(self.__finalizers) == 0 or all(not finalizer.alive for finalizer in self.__finalizers)


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

    def communicator(self):
        return self.__mpi_communicator

    @classmethod
    def set_base_communicator(cls, communicator):
        if cls.__instance is not None:
            raise exceptions.UsageError(
                'Cannot edit resource defaults after initialization.')
        cls.__mpi_communicator = communicator

    @classmethod
    def instance(cls):
        if _MPI is not None and cls.__mpi_communicator is None:
            cls.set_base_communicator(_MPI.COMM_WORLD)
        return cls()


@dataclasses.dataclass
class ResourceRequirements:
    """Named requirements for resources to be assigned for a task."""

    comm_size: typing.Optional[int] = dataclasses.field(default=None)
    """Requested communicator size."""


@contextmanager
def scoped_resources(allocation: BaseContext, requirements: ResourceRequirements = None):
    """Manage computing resources with a lifecycle.

    Get a Python context manager for a new resource assignment that is released
    on exiting the `with` block.

    For resources that include an MPI communicator, the new communicator is a
    sub-communicator of the communicator in *allocation*, if provided,
    *if requirements are specified*. (Note: for initial implementation without
    resource locking, an allocation without a requirements specification is
    simply borrowed for the context; no new communicator is split or duped and
    freed.)

    Args:
        allocation: available resources, or *None* to determine automatically.
        requirements: parameters for the resource assignment, such as MPI communicator size.

    *allocation* is required if *requirements* is not *None*. *requirements*
    may be validated in terms of the provided *allocation*.

    .. tbd wrt #4423 and #4422
        If ``gmxapi.version.has_feature('mpi_bindings') == True``, we assume that
        ``mpi4py`` can be imported, and the returned object is guaranteed to be an
        instance of :py:class:`mpi4py.MPI.Comm`. If the package was built without
        ``mpi4py`` integration, and *original_comm* is not a
        :py:class:`mpi4py.MPI.Comm` instance, the object produced will be a
        (special module object) :py:class:`_NullCommunicator`.

        Previous versions transformed :py:class:`~mpi4py.MPI.COMM_NULL` into a mock
        communicator before returning for consistency in trivial cases with or
        without :py:mod:`mpi4py`. This translation adds unnecessary indirection,
        and has been removed.

    When ``mpi4py`` is in use (``gmxapi.version.has_feature('mpi_bindings') == True``
    or *allocation.communicator()* is an :py:class:`mpi4py.MPI.Comm`), and
    *requested_size* is not ``None``, processes with a *original_comm* rank
    beyond the *requested_size* will receive :py:class:`~mpi4py.MPI.COMM_NULL`,
    consistent with :py:func:`mpi4py.MPI.Comm.Split()`.

    Callers should use ``allocation.communicator().Get_rank()`` to decide
    whether to use the new communicator.

    ..
        Previously, for *original_comm.Get_size()* or *requested_size* of 1,
        an unspecified "fake" communicator could be produced to avoid an
        ``MPI_Comm_split()``. While more efficient, the code complexity was
        problematic.

    TODO(#4079):
        This context manager is called "scoped_resources" instead of "scoped_assignment"
        to avoid implying there is any locking on the usage of the returned resources.
        gmxapi.operation.ResourceManager has a hard-coded protocol for preparing task
        resources that first recursively resolves inputs, then sets up the task runner.
        gmxapi.simulation.mdrun handles communication resources in its input adapter,
        LegacyImplementationSubscription, which ends up being reentrant for chained
        simulations, calling ResourceManager.local_input() _after_ suballocating
        Resources. This is a bad design.

    """
    if requirements is not None:
        if allocation is None:
            raise exceptions.UsageError(
                'Provide allocated resources from which to request a sub-allocation.')

    if allocation is None:
        # We need to think more about race conditions and other structural
        # details before allowing implicit automatic detection of resources.
        raise exceptions.FeatureNotAvailableError(
            'Implicit allocation is not currently available.'
        )

    base_communicator = allocation.communicator()
    # NOTE: We should decouple ensemble size from comm size in #4422 and allow
    # multiple ranks per ensemble member.
    requested_size = requirements.comm_size
    if requested_size is None:
        # TODO: Encapsulate getting the resources from the allocation.
        resources = ResourceAssignment()
        if base_communicator is not None:
            resources.borrow_communicator(base_communicator)
    else:
        try:
            max_size = base_communicator.Get_size()
        except AttributeError:
            max_size = 1
        if requested_size > max_size:
            error = f'Cannot get Comm size {requested_size} from {repr(allocation)}.'
            logger.error(
                error + ' '
                + 'To get a scoped communicator of a specific size, provide a sufficient '
                + 'communicator supporting the mpi4py.MPI.Comm interface.'
            )
            raise exceptions.UsageError(error)

        try:
            rank = base_communicator.Get_rank()
        except AttributeError:
            rank = None
        # TODO: How do we decide if we participate if we don't have a real communicator?
        if not hasattr(base_communicator, 'Split'):
            # TODO: Encapsulate getting the resources from the allocation.
            resources = ResourceAssignment()
            if base_communicator is not None:
                resources.borrow_communicator(base_communicator)
        else:
            if rank is not None and rank < requested_size:
                color = 0
            else:
                color = _MPI.UNDEFINED
            communicator = base_communicator.Split(color=color)
            resources = ResourceAssignment(communicator)

    try:
        yield resources
    finally:
        resources.close()


def assign_ensemble(allocation: ResourceAllocation, membership: typing.Iterable[int]) -> ResourceAssignment:
    """Assign resources to members of a new ensemble.

    Args:
        allocation: resource allocation from which to assign resources.
        membership: ensemble member ID for each rank in *allocation*

    *membership* IDs are in the *new* ensemble, sequenced by rank in *allocation*.
    *membership* values begin at zero and should be monotonic. I.e. the number
    of ensemble members is ``max(membership) + 1``. If the size of the original
    allocation is larger than the required assignment, let *membership* be a
    shorter sequence than ``allocation.communicator.Get_size()``

    Each ensemble member gets its own sub-communicator, split from the allocation.

    Ranks that are not part of any ensemble member get MPI_COMM_NULL.
    (I.e. when ``len(membership) > allocation.communicator().Get_size()``)

    This is a collective call. All ranks in *allocation* must call this function
    with the same arguments.
    """
    membership = list(membership)
    ensemble_size = max(membership) + 1
    if set(membership) != set(range(ensemble_size)):
        raise exceptions.UsageError(f'Invalid ensemble membership list: {membership}')

    required_comm_size = len(membership)

    base_communicator = allocation.communicator()
    # If requested size is 1, we should not need to require mpi4py.
    #
    # It is not safe to proceed with tMPI if multiple processes are
    # active and cannot be detected, but this is an issue without gmxapi,
    # as well, and we consider it a user error.
    #
    # For MPI-GROMACS and a requested size of 1, but no mpi4py, it is safe
    # to perform exactly one simulation, after which libgromacs will call
    # MPI_Finalize.
    #
    # With resolution of #4423, we will require mpi4py for MPI-GROMACS.
    #
    # As part of resolution of #4422, we will remove this requirement and
    # rely on internal (C++ level) MPI management to detect the MPI context
    # and call Init and Finalize the correct number of times.
    # TODO(#4423): Check for MPI bindings.
    if required_comm_size > 1:
        # Ensembles require mpi4py and a valid base communicator.
        if not hasattr(base_communicator, 'Get_size'):
            if _MPI is None:
                error = f'Ensemble work requires the mpi4py Python package.'
            else:
                error = f'Cannot get Comm size {required_comm_size} from {repr(allocation)}.'
            logger.error(
                error + ' '
                + 'To get a scoped communicator of a specific size, provide a parent '
                + 'communicator supporting the mpi4py.MPI.Comm interface.'
            )
            raise exceptions.UsageError(error)

    if _MPI is None:
        assert required_comm_size == 1
        assert base_communicator is None
        resources = ResourceAssignment()
    else:
        assert base_communicator is not None
        original_comm_size = base_communicator.Get_size()
        padding = original_comm_size - required_comm_size
        if padding < 0:
            raise exceptions.FeatureNotAvailableError(
                f'Cannot produce a sub-communicator of size {required_comm_size}'
                f' from a communicator of size {original_comm_size}.')
        split_color = membership + ([_MPI.UNDEFINED] * padding)
        assert len(split_color) == original_comm_size
        current_rank = base_communicator.Get_rank()
        communicator = base_communicator.Split(color=split_color[current_rank])
        resources = ResourceAssignment(communicator)
    return resources
