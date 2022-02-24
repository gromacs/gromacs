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

"""modify_input operation module

Provides implementation classes and user interface for gmxapi.modify_input.

The *modify_input* operation accepts simulation input, such as from *read_tpr*,
and allows aspects of the simulation input to be overridden.
"""

__all__ = ['modify_input']

import collections.abc
import inspect
import typing

import gmxapi
import gmxapi.abc
import gmxapi.exceptions
import gmxapi.operation as _op
import gmxapi.typing
from gmxapi.operation import InputCollectionDescription, ResourceManager
from .abc import ModuleObject

# Initialize module-level logger
from gmxapi import logger as root_logger

logger = root_logger.getChild('modify_input')
logger.info('Importing {}'.format(__name__))

_output_descriptors = (
    _op.OutputDataDescriptor('_simulation_input', str),
    _op.OutputDataDescriptor('parameters', dict)
)
_publishing_descriptors = {desc._name: gmxapi.operation.Publisher(desc._name, desc._dtype) for desc in
                           _output_descriptors}
_output = _op.OutputCollectionDescription(**{descriptor._name: descriptor._dtype for descriptor in
                                             _output_descriptors})


class OutputDataProxy(ModuleObject,
                      _op.DataProxyBase,
                      descriptors=_output_descriptors):
    """Implement the 'output' attribute of ReadTpr operations."""

    def __init__(self, *args, **kwargs):
        _op.DataProxyBase.__init__(self, *args, **kwargs)


class PublishingDataProxy(_op.DataProxyBase,
                          descriptors=_publishing_descriptors
                          ):
    """Manage output resource updates for ReadTpr operation."""


_output_factory = _op.OutputFactory(output_proxy=OutputDataProxy,
                                    output_description=_output,
                                    publishing_data_proxy=PublishingDataProxy)


class SessionResources(object):
    """Input and output run-time resources for a ModifyTPR operation."""

    def __init__(self, _simulation_input: str, parameters: dict, publisher: PublishingDataProxy):
        """Initialize resources for a gmxapi.operation context.

        Arguments:
            parameters: Source of new simulation parameters data.
            publisher: Output publishing resource provided by the Context.
            _simulation_input: Unspecified. Reserved for implementation details.

        """
        # For the initial implementation, _simulation_input is just a string that
        # will be interpreted as a TPR file path. The typing will change in the
        # near term as we build the data flow model. Before 1.0, unusual members
        # like this will be removed or hidden in the interface of the other
        # input objects.
        assert isinstance(_simulation_input, str)
        self._simulation_input = _simulation_input
        assert isinstance(parameters, collections.abc.Mapping)
        self.parameters = parameters
        self.output = publisher


#
# Helpers
#

# The '_simulation_input' input is a workaround until the data model is improved and the
# simulation module is more integrated.
# Note that when a read_tpr handle is passed with the "input" key word argument,
# its outputs will be mapped directly to modify_input inputs. The presence of
# a user-supplied 'parameters' argument overrides the 'parameters' value from 'input'.
# The lack of a default on 'parameters' makes it required.
_input = _op.InputCollectionDescription(
    [('_simulation_input', inspect.Parameter('_simulation_input',
                                             inspect.Parameter.POSITIONAL_OR_KEYWORD,
                                             annotation=str)),
     ('parameters', inspect.Parameter('parameters',
                                      inspect.Parameter.POSITIONAL_OR_KEYWORD,
                                      annotation=dict))
     ])


def _session_resource_factory(input: _op.InputPack, output: PublishingDataProxy, **kwargs
                              ) -> SessionResources:
    """Translate resources from the gmxapi.operation Context to the ReadTpr implementation."""
    # TODO: Either get rid of **kwargs or clarify the roadmap and timeline for doing so.
    filename = input.kwargs['_simulation_input']
    parameters = input.kwargs['parameters']
    return SessionResources(_simulation_input=filename, parameters=parameters, publisher=output)


def _standard_node_resource_factory(*args, **kwargs):
    """Translate Python UI input to the gmxapi.operation node builder inputs."""
    return _input.bind(*args, **kwargs)


def _run(resources: SessionResources):
    """Operation implementation in the gmxapi.operation module context."""
    # We combine the input from *input*, and key word arguments into the new output.
    for named_data in resources.output.keys():
        assert isinstance(named_data, str)
        assert hasattr(resources, named_data)
        source_value = getattr(resources, named_data)
        logger.debug('modify_input publishing {} to {}'.format(source_value, named_data))
        setattr(resources.output, named_data, source_value)


# Note: this is a class because we attach the input_description functionality,
# but all of its functionality can be composed from dispatching functions.
# TODO: Consider replacing this unique class with a pattern of composed instances of a common class.
class ResourceFactory(gmxapi.abc.ResourceFactory):
    """ModifyInput resource factory."""

    def __init__(self, target_context, source_context):
        self.source_context = source_context  # Determine input form of *create* method.
        self.target_context = target_context  # Context in which resources are consumed.

    # TODO: clean up the dispatching. What to do when input comes from multiple sources?
    # Use a typing overload or a single-dispatch functor for clearer input/result typing.
    def __call__(self, *args, **kwargs):
        # context is used for dispatching and annotates the Context of the other arguments.
        # context == None implies inputs are from Python UI.
        if self.source_context is None:
            # TODO: De-duplicate re: StandardDirector.resource_factory()
            if isinstance(self.target_context, _op.Context):
                return _standard_node_resource_factory(*args, **kwargs)
        if isinstance(self.source_context, _op.Context):
            # TODO: Check whether the consumer is a Context.NodeBuilder or an operation runner.
            # We don't yet use this dispatcher for building nodes, so assume we are launching a session.
            assert 'input' in kwargs
            assert 'output' in kwargs
            return _session_resource_factory(input=kwargs['input'], output=kwargs['output'])
        raise gmxapi.exceptions.MissingImplementationError(
            'No translation from {} context to {}'.format(self.source_context, self.target_context))

    @typing.overload
    def input_description(self, context: _op.Context) -> _op.InputDescription:
        ...

    def input_description(self, context: gmxapi.abc.Context):
        if isinstance(context, _op.Context):
            return StandardInputDescription()
        raise gmxapi.exceptions.MissingImplementationError('No input description available for {} context'.format(context))


class StandardInputDescription(_op.InputDescription):
    """Provide the ModifyInput input description in gmxapi.operation Contexts."""

    # TODO: Improve fingerprinting.
    # If _make_uid can't make sufficiently unique IDs, use additional "salt".
    # Without fingerprinting, we cannot consistently hash *input* across processes,
    # but we can consistently generate integers in the same sequence the first time
    # we see each distinct input.
    _next_uid: typing.ClassVar[int] = 0
    _uids: typing.ClassVar[typing.MutableMapping[int, int]] = {}

    @classmethod
    def _make_uid(cls, input) -> str:
        # TODO: Use input fingerprint for more useful identification.
        # WARNING: The built-in hash will use memory locations, and so will not be consistent across
        # process ranks, even if the input should be the same.
        salt = hash(input)
        if salt not in cls._uids:
            cls._uids[salt] = cls._next_uid
            cls._next_uid += 1
        else:
            logger.debug(
                f'Reissuing uid for modify_input({input}): {cls._uids[salt]}'
            )
        new_uid = 'modify_input_{}'.format(cls._uids[salt])
        return new_uid

    def signature(self) -> InputCollectionDescription:
        return _input

    def make_uid(self, input: _op.DataEdge) -> str:
        assert isinstance(input, _op.DataEdge)
        return self._make_uid(input)


class RegisteredOperation(_op.OperationImplementation, metaclass=_op.OperationMeta):
    """Provide the gmxapi compatible ModifyInput implementation."""

    # This is a class method to allow the class object to be used in gmxapi.operation._make_registry_key
    @classmethod
    def name(self) -> str:
        """Canonical name for the operation."""
        return 'modify_input'

    @classmethod
    def namespace(self) -> str:
        """modify_input is importable from the gmxapi module."""
        return 'gmxapi'

    @classmethod
    def director(cls, context: gmxapi.abc.Context):
        # Currently, we only have a Directory for the gmxapi.operation.Context
        if isinstance(context, _op.Context):
            return StandardDirector(context)


class StandardOperationHandle(_op.AbstractOperation, ModuleObject):
    """Handle used in Python UI or gmxapi.operation Contexts."""

    def __init__(self, resource_manager: ResourceManager):
        self.__resource_manager = resource_manager

    def run(self):
        self.__resource_manager.update_output()

    @property
    def output(self) -> OutputDataProxy:
        return self.__resource_manager.data()


class StandardDirector(gmxapi.abc.OperationDirector):
    """Direct the instantiation of a modify_input node in a gmxapi.operation Context.

    .. todo:: Compose this behavior in a more generic class.

    .. todo:: Describe where instances live.
    """

    def __init__(self, context: _op.Context):
        if not isinstance(context, _op.Context):
            raise gmxapi.exceptions.ValueError('StandardDirector requires a gmxapi.operation Context.')
        self.context = context

    def __call__(self, resources: _op.DataSourceCollection, label: str = None) -> StandardOperationHandle:
        builder = self.context.node_builder(operation=RegisteredOperation, label=label)

        builder.set_resource_factory(_session_resource_factory)
        builder.set_input_description(StandardInputDescription())
        builder.set_handle(StandardOperationHandle)

        runner_director = _op.RunnerDirector(runner=_run, allow_duplicate=True)

        builder.set_runner_director(runner_director)
        builder.set_output_factory(_output_factory)
        builder.set_resource_manager(ResourceManager)

        # Note: we have not yet done any dispatching based on *resources*. We should
        # translate the resources provided into the form that the Context can subscribe to
        # using the dispatching resource_factory.
        assert isinstance(resources, _op.DataSourceCollection)
        for name, source in resources.items():
            builder.add_input(name, source)

        handle = builder.build()
        assert isinstance(handle, StandardOperationHandle)
        return handle

    def handle_type(self, context: gmxapi.abc.Context):
        return StandardOperationHandle

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
        # TODO: Normalize and consolidate the context-based dispatching.
        if source is None:
            # `source is None` implies source is coming from UI.
            if isinstance(target, _op.Context):
                # Return a factory that can bind to function call arguments to produce a DataSourceCollection.
                return ResourceFactory(target_context=target, source_context=source)
        if isinstance(source, _op.Context):
            # The source is a gmxapi.operation.Context when the operation is being evaluated through
            # a gmxapi.operation.ResourceManager. i.e. at run time.
            return ResourceFactory(target_context=target, source_context=source)
        if isinstance(source, ModuleObject):
            if isinstance(target, _op.Context):
                # We are creating a node in gmxapi.operation.Context from another gmxapi.simulation operation.
                # This means that we want to subscribe to the subcontext instead of the gmxapi.operation.Context.
                # In the first draft, though, we just access a special payload.
                # Return a factory that will consume *_simulation_input* and *parameters*
                # members of a received object.
                logger.info('Building mdrun operation from source {}'.format(source))

                def simulation_input_workaround(input, parameters: dict):
                    source = input
                    # TODO: Normalize mechanism for obtaining SimulationInput references.
                    if hasattr(source, 'output'):
                        source = input.output
                    assert hasattr(source, '_simulation_input')
                    assert hasattr(source, 'parameters')
                    logger.info('modify_input receiving input {}: {}'.format(source._simulation_input.name,
                                                                      source._simulation_input.description))
                    source_collection = _input.bind(_simulation_input=source._simulation_input,
                                                    parameters=parameters)
                    logger.info('modify_input input bound as source collection {}'.format(source_collection))
                    return source_collection

                return simulation_input_workaround

        raise gmxapi.exceptions.ValueError('No dispatching from {} context to {}'.format(source, target))


# TODO: This operation is intended to be able to compose complete simulation input. E.g:
#  def modify_input(input: SimulationInput, parameters=None, topology=None, conformation=None, simulation_state=None):
def modify_input(input,
                 parameters: dict,
                 label: str = None,
                 context=None):
    """Modify simulation input with data flow operations.

    Given simulation input *input*, override components of simulation input with
    additional arguments, such as *parameters*.
    """
    handle_context = context
    if handle_context is not None:
        raise gmxapi.exceptions.MissingImplementationError(
            'context must be None. This factory is only for the Python UI right now.')

    target_context = _op.current_context()
    assert isinstance(target_context, _op.Context)
    # Get a director that will create a node in the standard context.
    # TODO: For clarity, restructure code to use module-local helper functors that produce callables with clear
    #  interfaces. (May clash with the current ABC scheme.)
    node_director: StandardDirector = _op._get_operation_director(RegisteredOperation, context=target_context)
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

    # TODO: Fix this protocol to do dispatching in the correct place.
    # The source Context here is None (the handle Context). The resources themselves
    # may be from different Contexts, so we should dispatch at the add_resource
    # builder method, not here in the director client.
    resource_factory: ResourceFactory = node_director.resource_factory(source=source_context, target=target_context)
    resources: _op.DataSourceCollection = resource_factory(input=input, parameters=parameters)
    handle = node_director(resources=resources, label=label)
    return handle
