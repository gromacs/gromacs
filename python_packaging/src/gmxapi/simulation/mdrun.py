#
# This file is part of the GROMACS molecular simulation package.
#
# Copyright (c) 2019, by the GROMACS development team, led by
# Mark Abraham, David van der Spoel, Berk Hess, and Erik Lindahl,
# and including many others, as listed in the AUTHORS file in the
# top-level source directory and at http://www.gromacs.org.
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
# http://www.gnu.org/licenses, or write to the Free Software Foundation,
# Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA.
#
# If you want to redistribute modifications to GROMACS, please
# consider that scientific software is very special. Version
# control is crucial - bugs must be traceable. We will be happy to
# consider code for inclusion in the official distribution, but
# derived work must not be called official GROMACS. Details are found
# in the README & COPYING files - if they are missing, get the
# official version at http://www.gromacs.org.
#
# To help us fund GROMACS development, we humbly ask that you cite
# the research papers on the package. Check out http://www.gromacs.org.

"""mdrun operation module

The mdrun operation (in its first draft) conforms to the user-level API, but does
not use the Python Context resource manager. It uses either the legacy 0.0.7
Context or its own Context, also implemented in this module.
"""

__all__ = ['mdrun']

import inspect
import os
import typing
from contextlib import contextmanager

import gmxapi
import gmxapi.abc
import gmxapi.operation as _op
from gmxapi import exceptions

# The following imports are not marked as public API.
from .abc import ModuleObject

from . import fileio
from . import workflow

# Initialize module-level logger
from gmxapi import logger as root_logger

logger = root_logger.getChild('mdrun')
logger.info('Importing {}'.format(__name__))


# Output in the gmxapi.operation Context.
# TODO: Consider using a single base class for the DataProxy, but have distinct
#  data descriptor behavior (or different descriptor implementations in different
#  subclasses) so that static code inspection can more easily determine the
#  attributes of the data proxies.
_output_descriptors = (
    _op.OutputDataDescriptor('_work_dir', str),
    _op.OutputDataDescriptor('trajectory', str),
    _op.OutputDataDescriptor('parameters', dict)
)
_publishing_descriptors = {desc._name: gmxapi.operation.Publisher(desc._name, desc._dtype) for desc in
                           _output_descriptors}
_output = _op.OutputCollectionDescription(**{descriptor._name: descriptor._dtype for descriptor in
                                             _output_descriptors})


class OutputDataProxy(_op.DataProxyBase,
                      descriptors=_output_descriptors):
    """Implement the 'output' attribute of MDRun operations."""


class PublishingDataProxy(_op.DataProxyBase,
                          descriptors=_publishing_descriptors):
    """Manage output resource updates for MDRun operation."""


_output_factory = _op.OutputFactory(output_proxy=OutputDataProxy,
                                    output_description=_output,
                                    publishing_data_proxy=PublishingDataProxy)


# Input in the gmxapi.operation Context for the dispatching runner.
# The default empty dictionary for parameters just means that there are no overrides
# to the parameters already provided in _simulation_input.
_input = _op.InputCollectionDescription(
    [('_simulation_input', inspect.Parameter('_simulation_input',
                                             inspect.Parameter.POSITIONAL_OR_KEYWORD,
                                             annotation=str)),
     ('parameters', inspect.Parameter('parameters',
                                      inspect.Parameter.POSITIONAL_OR_KEYWORD,
                                      annotation=dict,
                                      default=dict()))
     ])


def _standard_node_resource_factory(*args, **kwargs):
    """Translate Python UI input to the gmxapi.operation node builder inputs."""
    source_collection = _input.bind(*args, **kwargs)
    logger.info('mdrun input bound as source collection {}'.format(source_collection))
    return source_collection


@contextmanager
def scoped_communicator(original_comm, requested_size: int = None):
    from gmxapi.simulation.context import _acquire_communicator, _get_ensemble_communicator

    if requested_size is None:
        communicator = _acquire_communicator(communicator=original_comm)

    else:
        if original_comm is None or not hasattr(original_comm, 'Get_size'):
            raise exceptions.UsageError('A valid communicator must be provided when requesting a specific size.')
        original_comm_size = original_comm.Get_size()
        if original_comm_size < requested_size:
            raise exceptions.FeatureNotAvailableError(
                'Cannot produce a subcommunicator of size {} from a communicator of size {}.'.format(
                    requested_size,
                    original_comm_size
                ))
        assert original_comm_size >= requested_size
        communicator = _get_ensemble_communicator(original_comm, requested_size)

    try:
        yield communicator
    finally:
        communicator.Free()


class LegacyImplementationSubscription(object):
    """Input type representing a subscription to 0.0.7 implementation in gmxapi.operation.

    This input resource is a subscription to work that is dispatched to a sub-context.
    The resource can be created from the standard data of the simulation module.
    """

    def __init__(self, resource_manager: _op.ResourceManager):
        from .context import Context as LegacyContext
        import gmxapi._gmxapi as _gmxapi
        self._gmxapi = _gmxapi

        assert isinstance(resource_manager, _op.ResourceManager)
        # We use several details of the gmxapi.operation.Context.ResourceManager.
        # These dependencies can evolve into the subscription protocol between Contexts.

        # Configure and run a gmxapi 0.0.7 session.
        # 0. Determine ensemble width.
        # 1. Choose, create/check working directories.
        # 2. Create source TPR.
        # 3. Create workspec.
        # 3.5 Add plugin potentials, if any.
        # 4. Run.
        # 5. Collect outputs from context (including plugin outputs) and be ready to publish them.

        # Determine ensemble width
        ensemble_width = resource_manager.ensemble_width

        # Choose working directories
        # TODO: operation working directory naming scheme should be centrally well-defined.
        # Note that workflow.WorkSpec.uid is currently dependent on the input file parameter,
        # so we cannot create the input file path in the working directory based on WorkSpec.uid.
        workdir_list = ['{node}_{member}'.format(node=resource_manager.operation_id,
                                                 member=member)
                        for member in range(ensemble_width)]
        parameters_dict_list = [{}] * ensemble_width

        # This is a reasonable place to start using MPI ensemble implementation details.
        # We will want better abstraction in the future, but it is best if related filesystem
        # accesses occur from the same processes, consistently. Note that we already
        # handle some optimization and dependency reduction when the ensemble size is 1.
        # TODO: multithread and multiprocess alternatives to MPI ensemble management.

        # TODO: Allow user to provide communicator instead of implicitly getting COMM_WORLD
        with scoped_communicator(None) as context_comm:
            context_rank = context_comm.Get_rank()
            with scoped_communicator(context_comm, ensemble_width) as ensemble_comm:
                # Note that in the current implementation, extra ranks have nothing to do,
                # but they may have a dummy communicator, so be sure to skip those members
                # of the context_comm.
                if context_rank < ensemble_width:
                    assert ensemble_comm.Get_size() == ensemble_width
                    ensemble_rank = ensemble_comm.Get_rank()
                    # TODO: This should be a richer object that includes at least host information
                    #  and preferably the full gmxapi Future interface.
                    self.workdir = os.path.abspath(workdir_list[ensemble_rank])

                    with resource_manager.local_input(member=ensemble_rank) as input_pack:
                        source_file = input_pack.kwargs['_simulation_input']
                        parameters = input_pack.kwargs['parameters']
                        # If there are any other key word arguments to process from the gmxapi.mdrun
                        # factory call, do it here.

                    # TODO: We should really name this file with a useful input-dependent tag.
                    tprfile = os.path.join(self.workdir, 'topol.tpr')

                    expected_working_files = [tprfile]

                    if os.path.exists(self.workdir):
                        if os.path.isdir(self.workdir):
                            # Confirm that this is a restarted simulation.
                            # It is unspecified by the API, but at least through gmxapi 0.1,
                            # all simulations are initialized with a checkpoint file named state.cpt
                            # (see src/api/cpp/context.cpp)
                            checkpoint_file = os.path.join(self.workdir, 'state.cpp')
                            expected_working_files.append(checkpoint_file)

                            for file in expected_working_files:
                                if not os.path.exists(file):
                                    raise exceptions.ApiError(
                                        'Cannot determine working directory state: {}'.format(self.workdir))
                        else:
                            raise exceptions.ApiError(
                                'Chosen working directory path exists but is not a directory: {}'.format(self.workdir))
                    else:
                        # Build the working directory and input files.
                        os.mkdir(self.workdir)
                        sim_input = fileio.read_tpr(source_file)
                        for key, value in parameters.items():
                            try:
                                sim_input.parameters.set(key=key, value=value)
                            except _gmxapi.Exception as e:
                                raise exceptions.ApiError(
                                    'Bug encountered. Unknown error when trying to set simulation '
                                    'parameter {} to {}'.format(key, value)
                                ) from e

                        fileio.write_tpr_file(output=tprfile, input=sim_input)
                    logger.info('Created {} on rank {}'.format(tprfile, context_rank))

                    # Gather the actual outputs from the ensemble members.
                    if hasattr(ensemble_comm, 'allgather'):
                        # We should not assume that abspath expands the same on different MPI ranks.
                        workdir_list = ensemble_comm.allgather(self.workdir)
                        tpr_filenames = ensemble_comm.allgather(tprfile)
                        parameters = fileio.read_tpr(tprfile).parameters.extract()
                        parameters_dict_list = ensemble_comm.allgather(parameters)
                    else:
                        workdir_list = [os.path.abspath(workdir) for workdir in workdir_list]
                        # TODO: If we use better input file names, they need to be updated in multiple places.
                        tpr_filenames = [os.path.join(workdir, 'topol.tpr') for workdir in workdir_list]
                        parameters_dict_list = [fileio.read_tpr(tprfile).parameters.extract() for tprfile in tpr_filenames]

                    logger.debug('Context rank {} acknowledges working directories {}'.format(context_rank,
                                                                                             workdir_list))
                    logger.debug('Operation {}:{} will use {}'.format(resource_manager.operation_id,
                                                                      ensemble_rank,
                                                                      self.workdir
                                                                      ))
                    # TODO: We have not exposed the ability to pass any run time parameters to mdrun.
                    work = workflow.from_tpr(tpr_filenames)
                    self.workspec = work.workspec
                    context = LegacyContext(work=self.workspec, workdir_list=workdir_list, communicator=ensemble_comm)
                    self.simulation_module_context = context
                    # Go ahead and execute immediately. No need for lazy initialization in this basic case.
                    with self.simulation_module_context as session:
                        session.run()
                        # TODO: There may be some additional results that we need to extract...
                    # end: if context_rank < ensemble_width

                # end scoped_communicator: ensemble_comm

            if context_comm.Get_size() > 1:
                context_comm.bcast(workdir_list, root=0)
            # end scoped_communicator: context_comm

        self.workdir = workdir_list
        self.parameters = parameters_dict_list


class SubscriptionSessionResources(object):
    """Input and output run-time resources for a MDRun subscription.

    A generalization of this class is the probably the main hook for customizing the resources
    provided to the operation at run time.

    .. todo:: Better factoring of SessionResources, ResourceFactory, Director.resource_factory.
    """

    def __init__(self, input: LegacyImplementationSubscription, output: PublishingDataProxy):
        assert isinstance(input, LegacyImplementationSubscription)
        assert isinstance(output, PublishingDataProxy)
        self.output = output
        member_id = self.output._client_identifier
        # Before enabling the following, be sure we understand what is happening.
        # if member_id is None:
        #     member_id = 0
        self.workdir = input.workdir[member_id]
        self.parameters = input.parameters[member_id]


class SubscriptionPublishingRunner(object):
    """Handle execution in the gmxapi.operation context as a subscription to the gmxapi.simulation.context."""
    input_resource_factory = LegacyImplementationSubscription

    def __init__(self, resources: SubscriptionSessionResources):
        assert isinstance(resources, SubscriptionSessionResources)
        # Note that the resources contain a reference to a simulation ensemble that has already run.
        self.resources = resources

    def run(self):
        """Operation implementation in the gmxapi.operation module context."""
        publisher = self.resources.output
        publisher._work_dir = self.resources.workdir
        publisher.parameters = self.resources.parameters
        # TODO: Make the return value a trajectory handle rather than a file path.
        # TODO: Decide how to handle append vs. noappend output.
        # TODO: More rigorous handling of the trajectory file(s)
        # We have no way to query the name of the trajectory produced, and we
        # have avoided exposing the ability to specify it, so we have to assume
        # GROMACS default behavior.
        publisher.trajectory = os.path.join(self.resources.workdir, 'traj.trr')


_next_uid = 0


def _make_uid(input) -> str:
    # TODO: Use input fingerprint for more useful identification.
    salt = hash(input)
    global _next_uid
    new_uid = 'mdrun_{}_{}'.format(_next_uid, salt)
    _next_uid += 1
    return new_uid


#
# Implementation
#


class ResourceManager(gmxapi.operation.ResourceManager):
    """Manage resources for the MDRun operation in the gmxapi.operation contexts.

    Extends gmxapi.operation.ResourceManager to tolerate non-standard data payloads.
    Futures managed by this resource manager may contain additional attributes.
    """

    def future(self, name: str, description: _op.ResultDescription):
        tpr_future = super().future(name=name, description=description)
        return tpr_future

    def data(self) -> OutputDataProxy:
        return OutputDataProxy(self)

    def update_output(self):
        """Override gmxapi.operation.ResourceManager.update_output because we handle paralellism as 0.0.7."""
        # For the moment, this is copy-pasted from gmxapi.operation.ResourceManager,
        # but the only part we need to override is the ensemble handling at `for i in range(self.ensemble_width)`
        # TODO: Reimplement as the resource factory and director for the operation target context.
        if not self.done():
            self.__operation_entrance_counter += 1
            if self.__operation_entrance_counter > 1:
                raise exceptions.ProtocolError('Bug detected: resource manager tried to execute operation twice.')
            if not self.done():
                # TODO: rewrite with the pattern that this block is directing and then resolving an operation in the
                #  operation's library/implementation context.

                ###
                # Note: this is the resource translation from gmxapi.operation context
                # to the dispatching runner director. It uses details of the gmxapi.operation.Context
                # and of the operation.

                # TODO: Dispatch/discover this resource factory from a canonical place.
                assert hasattr(self._runner_director, 'input_resource_factory')
                # Create on all ranks.
                input = self._runner_director.input_resource_factory(self)
                # End of action of the InputResourceDirector[Context, MdRunSubscription].
                ###

                # We are giving the director a resource that contains the subscription
                # to the dispatched work.

                publishing_resources = self.publishing_resources()
                for member in range(self.ensemble_width):
                    with publishing_resources(ensemble_member=member) as output:
                        resources = self._resource_factory(input=input, output=output)
                        runner = self._runner_director(resources)
                        runner.run()


class StandardInputDescription(_op.InputDescription):
    """Provide the ReadTpr input description in gmxapi.operation Contexts."""

    def signature(self) -> _op.InputCollectionDescription:
        return _input

    def make_uid(self, input: _op.DataEdge) -> str:
        return _make_uid(input)


class RegisteredOperation(_op.OperationImplementation, metaclass=_op.OperationMeta):
    """Provide the gmxapi compatible ReadTpr implementation."""

    # This is a class method to allow the class object to be used in gmxapi.operation._make_registry_key
    @classmethod
    def name(self) -> str:
        """Canonical name for the operation."""
        return 'mdrun'

    @classmethod
    def namespace(self) -> str:
        """read_tpr is importable from the gmxapi module."""
        return 'gmxapi'

    @classmethod
    def director(cls, context: gmxapi.abc.Context) -> _op.OperationDirector:
        if isinstance(context, _op.Context):
            return StandardDirector(context)


class StandardOperationHandle(_op.AbstractOperation):
    """Handle used in Python UI or gmxapi.operation Contexts."""

    def __init__(self, resource_manager: ResourceManager):
        self.__resource_manager = resource_manager

    def run(self):
        self.__resource_manager.update_output()

    @property
    def output(self) -> OutputDataProxy:
        return self.__resource_manager.data()


class StandardDirector(gmxapi.abc.OperationDirector):
    """Direct the instantiation of a read_tpr node in a gmxapi.operation Context.

    .. todo:: Compose this behavior in a more generic class.

    .. todo:: Describe where instances live.
    """

    def __init__(self, context: _op.Context):
        if not isinstance(context, _op.Context):
            raise gmxapi.exceptions.ValueError('StandardDirector requires a gmxapi.operation Context.')
        self.context = context

    def __call__(self, resources: _op.DataSourceCollection, label: str = None) -> StandardOperationHandle:
        builder = self.context.node_builder(operation=RegisteredOperation, label=label)

        builder.set_resource_factory(SubscriptionSessionResources)
        builder.set_input_description(StandardInputDescription())
        builder.set_handle(StandardOperationHandle)
        # The runner in the gmxapi.operation context is the servicer for the legacy context.
        builder.set_runner_director(SubscriptionPublishingRunner)
        builder.set_output_factory(_output_factory)
        builder.set_resource_manager(ResourceManager)

        # Note: we have not yet done any dispatching based on *resources*. We should
        # translate the resources provided into the form that the Context can subscribe to
        # using the dispatching resource_factory. In the second draft, this operation
        # is dispatched to a SimulationModuleContext, which can be subscribed directly
        # to sources that are either literal filenames or gmxapi.simulation sources,
        # while standard Futures can be resolved in the standard context.
        #
        assert isinstance(resources, _op.DataSourceCollection)
        for name, source in resources.items():
            builder.add_input(name, source)

        handle = builder.build()
        assert isinstance(handle, StandardOperationHandle)
        return handle

    def handle_type(self, context: gmxapi.abc.Context):
        return StandardOperationHandle

    # Developer note: the Director instance is a convenient place to get a dispatching
    # factory. The Director may become generic or more universal, but the resource_factory
    # would likely not be typed on the generic parameters of the Director class.
    # Instead, it is likely a generic function with its own TypeVar parameters.
    def resource_factory(self,
                         source: typing.Union[gmxapi.abc.Context, ModuleObject, None],
                         target: gmxapi.abc.Context = None):
        # Distinguish between the UIContext, in which input is in the form
        # of function call arguments, and the StandardContext, implemented in
        # gmxapi.operation. UIContext is probably a virtual context that is
        # asserted by callers in order to get a factory that normalizes UI input
        # for the StandardContext.
        #
        if target is None:
            target = self.context
        if source is None:
            if isinstance(target, _op.Context):
                # Return a factory that can bind to function call arguments to produce a DataSourceCollection.
                return _standard_node_resource_factory
        if isinstance(source, _op.Context):
            return SubscriptionSessionResources
        if isinstance(source, ModuleObject):
            if isinstance(target, _op.Context):
                # We are creating a node in gmxapi.operation.Context from another gmxapi.simulation operation.
                # This means that we want to subscribe to the subcontext instead of the gmxapi.operation.Context.
                # In the first draft, though, we just access a special payload.
                # Return a factory that will consume *_simulation_input* and *parameters*
                # members of a received object.
                logger.info('Building mdrun operation from source {}'.format(source))

                def simulation_input_workaround(input):
                    source = input
                    if hasattr(source, 'output'):
                        source = input.output
                    assert hasattr(source, '_simulation_input')
                    assert hasattr(source, 'parameters')
                    logger.info('mdrun receiving input {}: {}'.format(source._simulation_input.name,
                                                                      source._simulation_input.description))
                    source_collection = _input.bind(_simulation_input=source._simulation_input,
                                                    parameters=source.parameters)
                    logger.info('mdrun input bound as source collection {}'.format(source_collection))
                    return source_collection

                return simulation_input_workaround

        raise gmxapi.exceptions.ValueError('No dispatching from {} context to {}'.format(source, target))


def mdrun(input, label: str = None, context=None):
    """MD simulation operation.

    Arguments:
        input : valid simulation input

    Returns:
        runnable operation to perform the specified simulation

    The *output* attribute of the returned operation handle contains dynamically
    determined outputs from the operation.

    `input` may be a TPR file name or a an object providing the SimulationInput interface.

    Note:
        New function names will be appearing to handle tasks that are separate

        "simulate" is plausibly a dispatcher or base class for various tasks
        dispatched by mdrun. Specific work factories are likely "minimize,"
        "test_particle_insertion," "legacy_simulation" (do_md), or "simulation"
        composition (which may be leap-frog, vv, and other algorithms)
    """
    handle_context = context
    if handle_context is not None:
        raise gmxapi.exceptions.NotImplementedError(
            'context must be None. This factory is only for the Python UI right now.')

    target_context = _op.current_context()
    assert isinstance(target_context, _op.Context)
    # Get a director that will create a node in the standard context.
    node_director = _op._get_operation_director(RegisteredOperation, context=target_context)
    assert isinstance(node_director, StandardDirector)
    # TODO: refine this protocol
    assert handle_context is None
    # Examine the input.
    if isinstance(input, ModuleObject):
        # The input is from read_tpr or modify_input.
        source_context = input
    else:
        # Allow automatic dispatching
        source_context = None

    resource_factory = node_director.resource_factory(source=source_context, target=target_context)
    resources = resource_factory(input)
    handle = node_director(resources=resources, label=label)
    # Note: One effect of the assertions above is to help the type checker infer
    # the return type of the handle. It is hard to convince the type checker that
    # the return value of the node builder is up-cast. We might be able to resolve
    # this by making both get_operation_director and ReadTprImplementation
    # generics of the handle type, or using the handle type as the key for
    # get_operation_director.
    return handle
