#
# This file is part of the GROMACS molecular simulation package.
#
# Copyright (c) 2019,2021,2022, by the GROMACS development team, led by
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

"""read_tpr operation module

Provides implementation classes and user interface for gmxapi.read_tpr.
"""

import inspect
import typing

import gmxapi
import gmxapi.abc
import gmxapi.exceptions
import gmxapi.operation as _op
import gmxapi.typing
from gmxapi.operation import InputCollectionDescription
from gmxapi.simulation.abc import ModuleObject

from . import fileio

# Initialize module-level logger
from gmxapi import logger as root_logger

logger = root_logger.getChild('read_tpr')
logger.info('Importing {}'.format(__name__))


#
# Interface classes and internal details
#

# TODO: The output of read_tpr and modify_input should be the same and part of the
#  simulation module specification. Such output is either a special type of output proxy
#  or Future.
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
    """Input and output run-time resources for a ReadTpr operation."""
    def __init__(self, tpr_filename, publisher: PublishingDataProxy):
        self.tpr_object = fileio.TprFile(filename=tpr_filename, mode='r')
        self.output = publisher


#
# Helpers
#


_input = _op.InputCollectionDescription(
    [('filename', inspect.Parameter('filename',
                                    inspect.Parameter.POSITIONAL_OR_KEYWORD,
                                    annotation=str))])


# TODO: Clarify. The actual input and output arguments passed are customized for this operation.
def _session_resource_factory(input: _op.InputPack, output: 'PublishingDataProxy',
                              **kwargs
                              ) -> SessionResources:
    """Translate resources from the gmxapi.operation Context to the ReadTpr implementation."""
    # TODO: Either get rid of **kwargs or clarify the roadmap and timeline for doing so.
    filename = input.kwargs['filename']
    return SessionResources(tpr_filename=filename, publisher=output)


def _standard_node_resource_factory(*args, **kwargs):
    """Translate Python UI input to the gmxapi.operation node builder inputs."""
    return _input.bind(*args, **kwargs)


def _run(resources: SessionResources):
    """Operation implementation in the gmxapi.operation module context."""
    # TODO: Implement for other source/target contexts. We don't always need to
    #  produce all outputs.
    with resources.tpr_object as fh:
        params = fh._tprFileHandle.params().extract()
        resources.output.parameters = params
        resources.output._simulation_input = fh.filename


#
# Implementation
#

# Note: we borrow the implementation from operation.ResourceManager for now,
# but in the future we want the implementations to either be decoupled or
# for implementations in a given context to be coupled to details that are clearly
# and explicitly related to that context. Right now, operation.ResourceManager
# is tied to the implementation of Contexts in gmxapi.operation, but that is not
# sufficiently clear and explicit.
class ResourceManager(gmxapi.operation.ResourceManager):
    """Manage resources for the ReadTpr operation in the gmxapi.operation contexts.

    Extends gmxapi.operation.ResourceManager to tolerate non-standard data payloads.
    Futures managed by this resource manager may contain additional attributes.
    """
    def future(self, name: str, description: _op.ResultDescription):
        tpr_future = super().future(name=name, description=description)
        return tpr_future

    def data(self) -> OutputDataProxy:
        return OutputDataProxy(self)

    def update_output(self):
        logger.debug('Updating output for {}.'.format(self.operation_id))
        super().update_output()
        for descriptor in _output_descriptors:
            name = descriptor._name
            if not self.is_done(name):
                raise gmxapi.exceptions.ApiError('Expected output {} not updated.'.format(name))

# TODO: Consider making Generic in source and target context type variables,
#  or leave unspecified and use generic function or pair of single_dispatch functions.
#  Need to know the right home for input_description, if needed.
class ResourceFactory(gmxapi.abc.ResourceFactory):
    """ReadTpr resource factory.

    Generic class for creating resources passed to read_tpr implementation details.
    Dispatching may occur based on the source and target Context of factory action.
    """
    def __init__(self, target_context, source_context):
        """Initialize an instance to support read_tpr action.

        Arguments:
            source_context: The source of the resources in the calling scope.
            target_context: The Context in which the factory product will be consumed.

        """
        self.source_context = source_context  # Determine input form of *create* method.
        self.target_context = target_context  # Context in which resources are consumed.

    # TODO: clean up the dispatching. What to do when input comes from multiple sources?
    # Use a typing overload or a single-dispatch functor for clearer input/result typing.
    def __call__(self, *args, **kwargs):
        """Create the resource product.

        Dispatch to appropriate factory functionality based on the *source_context*
        and *target_context*.
        """
        # context is used for dispatching and annotates the Context of the other arguments.
        # context == None implies inputs are from Python UI.
        if self.source_context is None:
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
        """Get an input description usable in the indicated Context.

        Arguments:
            context: Context for which to dispatch the generation of input description.

        Overrides gmxapi.abc.ResourceFactory for collaboration between this
        resource factory and a target Context.
        """
        # TODO: Clarify exposure and scope as gmxapi.abc.ResourceFactory is refined.
        # The expected use case is that a Context implementation may consult the
        # input_description when preparing input to a resource factory, or even
        # when determining in which Context a resource should be created or an
        # operation dispatched. It seems sensible that the ResourceFactory should
        # be able to describe its possible inputs. Instance functions are the
        # most common and simple Python detail, so it makes sense for this to be
        # an instance method instead of classmethod or staticmethod, until shown
        # otherwise.
        # Also note that ResourceFactory is likely to be a Generic class or a
        # container for composed functionality. Customizing method signatures
        # based on instance data is more brittle than behavior determined at the
        # class level. The behavior of this function is determined at the class
        # level to be a dispatcher, which may utilize instance data, but as an
        # implementation detail that is not the business of the caller.
        # This should probably be implemented in terms of the standard Python
        # functools.single_dispatch generic function, but it is cleaner if the
        # method itself is not generic beyond the level of typing overloads.
        if isinstance(context, _op.Context):
            return StandardInputDescription()
        raise gmxapi.exceptions.MissingImplementationError('No input description available for {} context'.format(context))


class StandardInputDescription(_op.InputDescription):
    """Provide the ReadTpr input description in gmxapi.operation Contexts."""

    # TODO: Improve fingerprinting.
    # If _make_uid can't make sufficiently unique IDs, use additional "salt":
    # _next_uid = 0

    @staticmethod
    def _make_uid(input) -> str:
        # TODO: Use input fingerprint for more useful identification.
        salt = hash(input)
        # If can't make sufficiently unique IDs, use additional "salt"
        # from a class data member. E.g.
        #     new_uid = 'read_tpr_{}_{}'.format(_next_uid, salt)
        #     _next_uid += 1
        new_uid = 'read_tpr_{}'.format(salt)
        return new_uid

    def signature(self) -> InputCollectionDescription:
        return _input

    def make_uid(self, input: _op.DataEdge) -> str:
        assert isinstance(input, _op.DataEdge)
        return self._make_uid(input)


class RegisteredOperation(_op.OperationImplementation, metaclass=_op.OperationMeta):
    """Provide the gmxapi compatible ReadTpr implementation."""

    # This is a class method to allow the class object to be used in gmxapi.operation._make_registry_key
    @classmethod
    def name(self) -> str:
        """Canonical name for the operation."""
        return 'read_tpr'

    @classmethod
    def namespace(self) -> str:
        """read_tpr is importable from the gmxapi module."""
        return 'gmxapi'

    @classmethod
    def director(cls, context: gmxapi.abc.Context) -> gmxapi.abc.OperationDirector:
        if isinstance(context, _op.Context):
            return StandardDirector(context)
        raise gmxapi.exceptions.MissingImplementationError(
            'No dispatcher for context {} of type {}'.format(context, type(context)))


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
    """Direct the instantiation of a read_tpr node in a gmxapi.operation Context.

    .. todo:: Compose this behavior in a more generic class.

    .. todo:: Describe where instances live.
    """
    def __init__(self, context: _op.Context):
        if not isinstance(context, _op.Context):
            raise gmxapi.exceptions.ValueError('StandardDirector requires a gmxapi.operation Context.')
        self.context = context

    def __call__(self, resources: _op.DataSourceCollection, label: str) -> StandardOperationHandle:
        builder = self.context.node_builder(operation=RegisteredOperation, label=label)

        builder.set_resource_factory(_session_resource_factory)
        builder.set_input_description(StandardInputDescription())
        builder.set_handle(StandardOperationHandle)

        runner_director = _op.RunnerDirector(
            runner=_run
        )

        builder.set_runner_director(runner_director)
        builder.set_output_factory(_output_factory)

        # Note: we have not yet done any dispatching based on *resources*. We should
        # translate the resources provided into the form that the Context can subscribe to
        # using the dispatching resource_factory. In the second draft, this operation
        # is dispatched to a SimulationModuleContext, which can be subscribed directly
        # to sources that are either literal filenames or gmxapi.simulation sources,
        # while standard Futures can be resolved in the standard context.
        #
        assert isinstance(resources, _op.DataSourceCollection)
        assert 'filename' in resources
        builder.add_input('filename', resources['filename'])

        handle = builder.build()
        assert isinstance(handle, StandardOperationHandle)
        return handle

    def handle_type(self, context: gmxapi.abc.Context):
        return StandardOperationHandle

    def resource_factory(self,
                         source: typing.Union[gmxapi.abc.Context, None],
                         target: gmxapi.abc.Context = None) -> ResourceFactory:
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
                return ResourceFactory(target_context=target, source_context=source)
        if isinstance(source, _op.Context):
            return ResourceFactory(target_context=target, source_context=source)
        raise gmxapi.exceptions.ValueError('No dispatching from {} context to {}'.format(source, target))


def read_tpr(filename, label: str = None, context=None):
    """

    Arguments:
        filename: input file name
        label: optional human-readable label with which to tag the new node
        context: Context in which to return a handle to the new node.
                 Use default (None) for Python scripting interface

    Returns:
        Reference (handle) to the new operation instance (node).

    """
    handle_context = context
    if handle_context is not None:
        raise gmxapi.exceptions.MissingImplementationError(
            'context must be None. This factory is only for the Python UI right now.')

    # 1. Handle node creation in the scripting interface.
    # When *context* input is None, dispatch to the current Context. Confirm that
    # it is a standard context from the gmxapi.operation module, translate the
    # input into that context, and create the node.

    # 2. Dispatch to SimulationModuleContext.
    # Operation is not fully implemented in gmxapi.operation context. When creating
    # a node in a gmxapi.operation context, dispatch and subscribe to a SimulationModuleContext,
    # in which

    # 3. Handle consumption in SimulationModuleContext.
    # When a consuming operation is native to the SimulationModuleContext,
    # detect that a richer interface is available and use it. Possible implementation:
    # Chain of Responsibility: subscribe() is serviced by the Context "closest"
    # to the source. subscriber is the most native compatible Context for the consuming operation
    # but decays to the context in which the handle is being created. subscribe()
    # can accept a callback function and an indication of the form of the message
    # to send (a consuming Context or ResourceFactory).

    # TODO: Other types of input, such as File placeholders.

    target_context = _op.current_context()
    assert isinstance(target_context, _op.Context)
    # Get a director that will create a node in the standard context.
    node_director = _op._get_operation_director(RegisteredOperation, context=target_context)
    assert isinstance(node_director, StandardDirector)
    # TODO: refine this protocol
    assert handle_context is None
    resource_factory = node_director.resource_factory(source=handle_context, target=target_context)
    resources = resource_factory(filename=filename)
    handle = node_director(resources=resources, label=label)
    # Note: One effect of the assertions above is to help the type checker infer
    # the return type of the handle. It is hard to convince the type checker that
    # the return value of the node builder is up-cast. We might be able to resolve
    # this by making both get_operation_director and ReadTprImplementation
    # generics of the handle type, or using the handle type as the key for
    # get_operation_director.
    return handle
