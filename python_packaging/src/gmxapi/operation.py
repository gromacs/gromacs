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

"""Define gmxapi-compliant Operations

Provide decorators and base classes to generate and validate gmxapi Operations.

Nodes in a work graph are created as instances of Operations. An Operation factory
accepts well-defined inputs as key word arguments. The object returned by such
a factory is a handle to the node in the work graph. It's ``output`` attribute
is a collection of the Operation's results.

function_wrapper(...) produces a wrapper that converts a function to an Operation
factory. The Operation is defined when the wrapper is called. The Operation is
instantiated when the factory is called. The function is executed when the Operation
instance is run.

The framework ensures that an Operation instance is executed no more than once.
"""

__all__ = ['computed_result',
           'function_wrapper',
           ]

import abc
import collections
import functools
import importlib.util
import inspect
import typing
import weakref
from contextlib import contextmanager

import gmxapi as gmx
from gmxapi import datamodel
from gmxapi import exceptions
from gmxapi import logger as root_logger

# Initialize module-level logger
logger = root_logger.getChild('operation')
logger.info('Importing {}'.format(__name__))


class ResultDescription:
    """Describe what will be returned when `result()` is called."""

    def __init__(self, dtype=None, width=1):
        assert isinstance(dtype, type)
        assert issubclass(dtype, valid_result_types)
        assert isinstance(width, int)
        self._dtype = dtype
        self._width = width

    @property
    def dtype(self) -> type:
        """node output type"""
        return self._dtype

    @property
    def width(self) -> int:
        """ensemble width"""
        return self._width


class OutputData(object):
    """Encapsulate the description and storage of a data output."""
    def __init__(self, name: str, description: ResultDescription):
        assert name != ''
        self._name = name
        assert isinstance(description, ResultDescription)
        self._description = description
        self._done = [False] * self._description.width
        self._data = [None] * self._description.width

    @property
    def name(self):
        """The name of the published output."""
        return self._name

    # TODO: Change to regular member function and add ensemble member arg.
    @property
    def done(self):
        """Ensemble completion status for this output."""
        return all(self._done)

    def data(self, member: int = None):
        """Access the raw data for localized output for the ensemble or the specified member."""
        if not self.done:
            raise exceptions.ApiError('Attempt to read before data has been published.')
        if self._data is None or None in self._data:
            raise exceptions.ApiError('Data marked "done" but contains null value.')
        if member is not None:
            return self._data[member]
        else:
            return self._data

    def set(self, value, member: int):
        """Set local data and mark as completed.

        Used internally to maintain the data store.
        """
        if self._description.dtype == datamodel.NDArray:
            self._data[member] = datamodel.ndarray(value)
        else:
            self._data[member] = self._description.dtype(value)
        self._done[member] = True

    def reset(self):
        """Reinitialize the data store.

        Note:
            This is a workaround until the execution model is more developed.
        
        Todo:
            Remove this method when all operation handles provide factories to
            reinstantiate on new Contexts and/or with new inputs.
        
        """
        self._done = [False] * self._description.width
        self._data = [None] * self._description.width


class EnsembleDataSource(object):
    """A single source of data with ensemble data flow annotations.

    Note that data sources may be Futures.
    """
    def __init__(self, source=None, width=1, dtype=None):
        self.source = source
        self.width = width
        self.dtype = dtype

    def node(self, member: int):
        return self.source[member]

    def reset(self):
        protocols = ('reset', '_reset')
        for protocol in protocols:
            if hasattr(self.source, protocol):
                getattr(self.source, protocol)()


ResultTypeVar = typing.TypeVar('ResultTypeVar', *(str, bool, int, float, dict, datamodel.NDArray))
valid_result_types = ResultTypeVar.__constraints__

SourceTypeVar = typing.TypeVar('SourceTypeVar',
                               *(str, bool, int, float, dict, datamodel.NDArray, EnsembleDataSource))
valid_source_types = SourceTypeVar.__constraints__


class DataSourceCollection(collections.OrderedDict):
    """Store and describe input data handles for an operation.

    When created from InputCollectionDescription.bind(), the DataSourceCollection
    has had default values inserted.
    """

    def __init__(self, **kwargs):
        """Initialize from key/value pairs of named data sources.

        Data sources may be any of the basic gmxapi data types, gmxapi Futures
        of those types, or gmxapi ensemble data bundles of the above.

        Note that the checking and conditioning could be moved to one or more
        creation functions. In conjunction with an InputCollectionDescription,
        such creation functions could make decisions about automatically shaping
        the data flow topology or making conversions of data type or shape.
        """
        named_data = []
        for name, value in kwargs.items():
            if not isinstance(name, str):
                raise exceptions.TypeError('Data must be named with str type.')
            # TODO: Encapsulate handling of proferred data sources to Context details.
            # Preprocessed input should be self-describing gmxapi data types. Structured
            # input must be recursively (depth-first) converted to gmxapi data types.
            # TODO: Handle gmxapi Futures stored as dictionary elements!
            if not isinstance(value, valid_source_types):
                if isinstance(value, collections.abc.Iterable):
                    # Iterables here are treated as arrays, but we do not have a robust typing system.
                    # Warning: In the initial implementation, the iterable may contain Future objects.
                    # TODO: (#2993) Revisit as we sort out data shape and Future protocol.
                    value = datamodel.ndarray(value)
                elif hasattr(value, 'result'):
                    # A Future object.
                    pass
                else:
                    raise exceptions.ApiError('Cannot process data source {}'.format(value))
            named_data.append((name, value))
        super().__init__(named_data)

    def __setitem__(self, key: str, value: SourceTypeVar) -> None:
        if not isinstance(key, str):
            raise exceptions.TypeError('Data must be named with str type.')
        if not isinstance(value, valid_source_types):
            if isinstance(value, collections.abc.Iterable):
                # Iterables here are treated as arrays, but we do not have a robust typing system.
                # Warning: In the initial implementation, the iterable may contain Future objects.
                # TODO: (#2993) Revisit as we sort out data shape and Future protocol.
                value = datamodel.ndarray(value)
            elif hasattr(value, 'result'):
                # A Future object.
                pass
            else:
                raise exceptions.ApiError('Cannot process data source {}'.format(value))
        super().__setitem__(key, value)

    def reset(self):
        """Reset all sources in the collection."""
        for source in self.values():
            if hasattr(source, 'reset'):
                source.reset()
            if hasattr(source, '_reset'):
                source._reset()


def computed_result(function):
    """Decorate a function to get a helper that produces an object with Result behavior.

    When called, the new function produces an ImmediateResult object.

    The new function has the same signature as the original function, but can accept
    gmxapi data proxies, assuming the provided proxy objects represent types
    compatible with the original signature.

    Calls to `result()` return the value that `function` would return when executed
    in the local context with the inputs fully resolved.

    The API does not specify when input data dependencies will be resolved
    or when the wrapped function will be executed. That is, ``@computed_result``
    functions may force immediate resolution of data dependencies and/or may
    be called more than once to satisfy dependent operation inputs.
    """

    # Attempt to inspect function signature. Introspection can fail, so we catch
    # the various exceptions. We re-raise as ApiError, indicating a bug, because
    # we should do more to handle this intelligently and provide better user
    # feedback.
    # TODO: Figure out what to do with exceptions where this introspection
    #  and rebinding won't work.
    # ref: https://docs.python.org/3/library/inspect.html#introspecting-callables-with-the-signature-object
    try:
        sig = inspect.signature(function)
    except TypeError as T:
        raise exceptions.ApiError('Can not inspect type of provided function argument.') from T
    except ValueError as V:
        raise exceptions.ApiError('Can not inspect provided function signature.') from V

    @functools.wraps(function)
    def new_function(*args, **kwargs):
        # The signature of the new function will accept abstractions
        # of whatever types it originally accepted. This wrapper must
        # * Create a mapping to the original call signature from `input`
        # * Add handling for typed abstractions in wrapper function.
        # * Process arguments to the wrapper function into `input`

        # 1. Inspect the return annotation to determine valid gmxapi type(s)
        # 2. Generate a Result object advertising the correct type, bound to the
        #    Input and implementing function.
        # 3. Transform the result() data to the correct type.

        # TODO: (FR3+) create a serializable data structure for inputs discovered
        #  from function introspection.

        for name, param in sig.parameters.items():
            assert not param.kind == param.POSITIONAL_ONLY
        bound_arguments = sig.bind(*args, **kwargs)
        wrapped_function = function_wrapper()(function)
        handle = wrapped_function(**bound_arguments.arguments)
        output = handle.output
        # TODO: Find a type hinting / generic way to assert output attributes.
        return output.data

    return new_function


class OutputCollectionDescription(collections.OrderedDict):
    def __init__(self, **kwargs):
        """Create the output description for an operation node from a dictionary of names and types."""
        outputs = []
        for name, flavor in kwargs.items():
            if not isinstance(name, str):
                raise exceptions.TypeError('Output descriptions are keyed by Python strings.')
            # Multidimensional outputs are explicitly NDArray
            if issubclass(flavor, (list, tuple)):
                flavor = datamodel.NDArray
            assert issubclass(flavor, valid_result_types)
            outputs.append((name, flavor))
        super().__init__(outputs)


class InputCollectionDescription(collections.OrderedDict):
    """Describe acceptable inputs for an Operation.

    Generally, an InputCollectionDescription is an aspect of the public API by
    which an Operation expresses its possible inputs. This class includes details
    of the Python package.

    Keyword Arguments:
        parameters : A sequence of named parameter descriptions.

    Parameter descriptions are objects containing an `annotation` attribute
    declaring the data type of the parameter and, optionally, a `default`
    attribute declaring a default value for the parameter.

    Instances can be used as an ordered map of parameter names to gmxapi data types.

    Analogous to inspect.Signature, but generalized for gmxapi Operations.
    Additional notable differences: typing is normalized at initialization, and
    the bind() method does not return an object that can be directly used as
    function input. The object returned from bind() is used to construct a data
    graph Edge for subsequent execution.
    """

    def __init__(self, parameters: typing.Iterable[typing.Tuple[str, inspect.Parameter]]):
        """Create the input description for an operation node from a dictionary of names and types."""
        inputs = []
        for name, param in parameters:
            if not isinstance(name, str):
                raise exceptions.TypeError('Input descriptions are keyed by Python strings.')
            # Multidimensional inputs are explicitly NDArray
            dtype = param.annotation
            if issubclass(dtype, collections.abc.Iterable) \
                    and not issubclass(dtype, (str, bytes, collections.abc.Mapping)):
                # TODO: we can relax this with some more input conditioning.
                if dtype != datamodel.NDArray:
                    raise exceptions.UsageError(
                        'Cannot accept input type {}. Sequence type inputs must use NDArray.'.format(param))
            assert issubclass(dtype, valid_result_types)
            if hasattr(param, 'kind'):
                disallowed = any([param.kind == param.POSITIONAL_ONLY,
                                  param.kind == param.VAR_POSITIONAL,
                                  param.kind == param.VAR_KEYWORD])
                if disallowed:
                    raise exceptions.ProtocolError(
                        'Cannot wrap function. Operations must have well-defined parameter names.')
                kind = param.kind
            else:
                kind = inspect.Parameter.POSITIONAL_OR_KEYWORD
            if hasattr(param, 'default'):
                default = param.default
            else:
                default = inspect.Parameter.empty
            inputs.append(inspect.Parameter(name, kind, default=default, annotation=dtype))
        super().__init__([(input.name, input.annotation) for input in inputs])
        self.signature = inspect.Signature(inputs)

    @staticmethod
    def from_function(function):
        """Inspect a function to be wrapped.

        Used internally by gmxapi.operation.function_wrapper()

            Raises:
                exceptions.ProtocolError if function signature cannot be determined to be valid.

            Returns:
                InputCollectionDescription for the function input signature.
        """
        # First, inspect the function.
        assert callable(function)
        signature = inspect.signature(function)
        # The function must have clear and static input schema
        # Make sure that all parameters have clear names, whether or not they are used in a call.
        for name, param in signature.parameters.items():
            disallowed = any([param.kind == param.POSITIONAL_ONLY,
                              param.kind == param.VAR_POSITIONAL,
                              param.kind == param.VAR_KEYWORD])
            if disallowed:
                raise exceptions.ProtocolError(
                    'Cannot wrap function. Operations must have well-defined parameter names.')
            if param.name == 'input':
                raise exceptions.ProtocolError('Function signature includes the (reserved) "input" keyword argument.')
        description = collections.OrderedDict()
        for param in signature.parameters.values():
            if param.name == 'output':
                # Wrapped functions may accept the output parameter to publish results, but
                # that is not part of the Operation input signature.
                continue
            if param.annotation == param.empty:
                if param.default == param.empty or param.default is None:
                    raise exceptions.ProtocolError('Could not infer parameter type for {}'.format(param.name))
                dtype = type(param.default)
                if isinstance(dtype, collections.abc.Iterable) \
                        and not isinstance(dtype, (str, bytes, collections.abc.Mapping)):
                    dtype = datamodel.NDArray
            else:
                dtype = param.annotation
            description[param.name] = param.replace(annotation=dtype)
        return InputCollectionDescription(description.items())

    def bind(self, *args, **kwargs) -> DataSourceCollection:
        """Create a compatible DataSourceCollection from provided arguments.

        Pre-process input and function signature to get named input arguments.

        This is a helper function to allow calling code to characterize the
        arguments in a Python function call with hints from the factory that is
        initializing an operation. Its most useful functionality is to allow a
        factory to accept positional arguments where named inputs are usually
        required. It also allows data sources to participate in multiple
        DataSourceCollections with minimal constraints.

        Note that the returned object has had data populated from any defaults
        described in the InputCollectionDescription.

        See wrapped_function_runner() and describe_function_input().
        """
        # For convenience, accept *args, but convert to **kwargs to pass to Operation.
        # Factory accepts an unadvertised `input` keyword argument that is used as a default kwargs dict.
        # If present, kwargs['input'] is treated as an input "pack" providing _default_ values.
        input_kwargs = collections.OrderedDict()
        if 'input' in kwargs:
            provided_input = kwargs.pop('input')
            if provided_input is not None:
                input_kwargs.update(provided_input)
        # `function` may accept an `output` keyword argument that should not be supplied to the factory.
        for key, value in kwargs.items():
            if key == 'output':
                raise exceptions.UsageError('Invalid keyword argument: output (reserved).')
            input_kwargs[key] = value
        try:
            bound_arguments = self.signature.bind_partial(*args, **input_kwargs)
        except TypeError as e:
            raise exceptions.UsageError('Could not bind operation parameters to function signature.') from e
        assert 'output' not in bound_arguments.arguments
        bound_arguments.apply_defaults()
        assert 'input' not in bound_arguments.arguments
        input_kwargs = collections.OrderedDict([pair for pair in bound_arguments.arguments.items()])
        if 'output' in input_kwargs:
            input_kwargs.pop('output')
        return DataSourceCollection(**input_kwargs)


class ProxyDataDescriptor(object):
    """Base class for data descriptors used in DataProxyBase subclasses.

    Subclasses should either not define __init__ or should call the base class
    __init__ explicitly: super().__init__(self, name, dtype)
    """

    def __init__(self, name: str, dtype: ResultTypeVar = None):
        self._name = name
        # TODO: We should not allow dtype==None, but we currently have a weak data
        #  model that does not allow good support of structured Futures.
        if dtype is not None:
            assert isinstance(dtype, type)
            assert issubclass(dtype, valid_result_types)
        self._dtype = dtype


class DataProxyMeta(type):
    # Key word arguments consumed by __prepare__
    _prepare_keywords = ('descriptors',)

    @classmethod
    def __prepare__(mcs, name, bases, descriptors: collections.abc.Mapping = None):
        """Allow dynamic sub-classing.

        DataProxy class definitions are collections of data descriptors. This
        class method allows subclasses to give the descriptor names and type(s)
        in the class declaration as arguments instead of as class attribute
        assignments.

            class MyProxy(DataProxyBase, descriptors={name: MyDescriptor() for name in datanames}): pass

        Note:
            If we are only using this metaclass for the __prepare__ hook by the
            time we require Python >= 3.6, we could reimplement __prepare__ as
            DataProxyBase.__init_subclass__ and remove this metaclass.
        """
        if descriptors is None:
            return {}
        else:
            return descriptors

    def __new__(cls, name, bases: typing.Iterable, namespace, **kwargs):
        for key in kwargs:
            if key not in DataProxyMeta._prepare_keywords:
                raise exceptions.ApiError('Unexpected class creation keyword: {}'.format(key))
        # See note about DataProxyBase._reserved.
        if '_reserved' not in namespace and not any(hasattr(base, '_reserved') for base in bases):
            raise exceptions.ApiError(
                'We currently expect DataProxy classes to provide a list of reserved attribute names.')
        for key in namespace:
            # Here we can check conformance with naming and typing rules.
            assert isinstance(key, str)
            if key.startswith('_'):
                # Skip non-public attributes.
                continue
            descriptor = namespace[key]
            # The purpose of the current data proxies is to serve as a limited namespace
            # containing only descriptors of a certain type. In the future, these proxies
            # may be flattened into a facet of a richer OperationHandle class
            # (this metaclass may become a decorator on an implementation class),
            # but for now we check that the class is being used according to the
            # documented framework. A nearer term update could be to restrict the
            # type of the data descriptor:
            # TODO: Use a member type of the derived cls (or a mix-in base) to specify a particular
            #  ProxyDataDescriptor subclass.
            # Also, see note about DataProxyBase._reserved
            if not isinstance(descriptor, ProxyDataDescriptor):
                if key not in namespace['_reserved'] and not any(key in getattr(base, '_reserved') for base in
                                                                 bases if hasattr(base, '_reserved')):
                    raise exceptions.ApiError('Unexpected data proxy attribute {}: {}'.format(key, repr(descriptor)))
            else:
                assert isinstance(descriptor, ProxyDataDescriptor)
                if not isinstance(descriptor._name, str) or descriptor._name == '':
                    descriptor._name = key
                else:
                    if descriptor._name != key:
                        raise exceptions.ApiError(
                            'Descriptor internal name {} does not match attribute name {}'.format(
                                descriptor._name, key))
        return type.__new__(cls, name, bases, namespace)

    # TODO: This keyword argument stripping is not necessary in more recent Python versions.
    # When Python minimum required version is increased, check if we can remove this.
    def __init__(cls, name, bases, namespace, **kwargs):
        for key in kwargs:
            if key not in DataProxyMeta._prepare_keywords:
                raise exceptions.ApiError('Unexpected class initialization keyword: {}'.format(key))
        super().__init__(name, bases, namespace)


class DataProxyBase(object, metaclass=DataProxyMeta):
    """Limited interface to managed resources.

    Inherit from DataProxy to specialize an interface to an ``instance``.
    In the derived class, either do not define ``__init__`` or be sure to
    initialize the super class (DataProxy) with an instance of the object
    to be proxied.

    A class deriving from DataProxyBase allows its instances to provide a namespace
    for proxies to named data by defining attributes that are data descriptors
    (subclasses of ProxyDataDescriptor).
    The ProxyDataDescriptors are accessed as attributes of the
    data proxy instance or by iterating on items(). Attributes that are not
    ProxyDataDescriptors are possible, but will not be returned by items() which
    is a necessary part of gmxapi execution protocol.

    Acts as an owning handle to the resources provide by ``instance``,
    preventing the reference count of ``instance`` from going to zero for the
    lifetime of the proxy object.

    When sub-classing DataProxyBase, data descriptors can be passed as a mapping
    to the ``descriptors`` key word argument in the class declaration. This
    allows data proxy subclasses to be easily defined dynamically.

        mydescriptors = {'foo': Publisher('foo', int), 'data': Publisher('data', float)}
        ...
        class MyDataProxy(DataProxyBase, descriptors=mydescriptors): pass
        assert hasattr(MyDataProxy, 'foo')

    """
    # This class attribute (which subclasses are free to replace to augment) is an
    # indication of a problem with the current data model. If we are allowing
    # reserved words that would otherwise be valid data names, there is not a
    # compelling reason for separate data proxy classes: we throw away the assertion
    # that we are preparing a clean namespace and we could have accomplished the
    # class responsibilities in the Operation handle with just descriptor classes.
    # If we want the clean namespace, we should figure out how to keep this interface
    # from growing and/or have some "hidden" internal interface.
    _reserved = ('ensemble_width', 'items')

    # This class can be expanded to be the attachment point for a metaclass for
    # data proxies such as PublishingDataProxy or OutputDataProxy, which may be
    # defined very dynamically and concisely as a set of Descriptors and a type()
    # call.
    # If development in this direction does not materialize, then this base
    # class is not very useful and should be removed.
    def __init__(self, instance: 'SourceResource', client_id: int = None):
        """Get partial ownership of a resource provider.

        Arguments:
            instance : resource-owning object
            client_id : identifier for client holding the resource handle (e.g. ensemble member id)

        If client_id is not provided, the proxy scope is for all clients.
        """
        # TODO: Decide whether _resource_instance is public or not.
        # Note: currently commonly needed for subclass implementations.
        self._resource_instance = instance
        # Developer note subclasses should handle self._client_identifier == None
        self._client_identifier = client_id

    @property
    def ensemble_width(self) -> int:
        return self._resource_instance.width()

    @classmethod
    def items(cls):
        """Generator for tuples of attribute name and descriptor instance.

        This almost certainly doesn't do quite what we want...
        """
        for name, value in cls.__dict__.items():
            if isinstance(value, ProxyDataDescriptor):
                yield name, value


class Publisher(ProxyDataDescriptor):
    """Data descriptor for write access to a specific named data resource.

    For a wrapped function receiving an ``output`` argument, provides the
    accessors for an attribute on the object passed as ``output``. Maps
    read and write access by the wrapped function to appropriate details of
    the resource manager.

    Used internally to implement settable attributes on PublishingDataProxy.
    Allows PublishingDataProxy to be dynamically defined in the scope of the
    operation.function_wrapper closure. Each named output is represented by
    an instance of Publisher in the PublishingDataProxy class definition for
    the operation.

    Ref: https://docs.python.org/3/reference/datamodel.html#implementing-descriptors

    Collaborations:
    Relies on implementation details of ResourceManager.
    """

    def __get__(self, instance: DataProxyBase, owner):
        if instance is None:
            # The current access has come through the class attribute of owner class
            return self
        resource_manager = instance._resource_instance
        client_id = instance._client_identifier
        # TODO: Fix API scope.
        # Either this class is a detail of the same implementation as ResourceManager,
        # or we need to enforce that instance._resource_instance provides _data (or equivalent)
        assert isinstance(resource_manager, ResourceManager)
        if client_id is None:
            return getattr(resource_manager._data, self._name)
        else:
            return getattr(resource_manager._data, self._name)[client_id]

    def __set__(self, instance: DataProxyBase, value):
        resource_manager = instance._resource_instance
        # TODO: Fix API scope.
        # Either this class is a detail of the same implementation as ResourceManager,
        # or we need to enforce that instance._resource_instance provides _data (or equivalent)
        assert isinstance(resource_manager, ResourceManager)
        client_id = instance._client_identifier
        resource_manager.set_result(name=self._name, value=value, member=client_id)

    def __repr__(self):
        return 'Publisher(name={}, dtype={})'.format(self._name, self._dtype.__qualname__)


def define_publishing_data_proxy(output_description) -> typing.Type[DataProxyBase]:
    """Returns a class definition for a PublishingDataProxy for the provided output description."""
    # This dynamic type creation hides collaborations with things like make_datastore.
    # We should encapsulate these relationships in Context details, explicit collaborations
    # between specific operations and Contexts, and in groups of Operation definition helpers.

    descriptors = collections.OrderedDict([(name, Publisher(name)) for name in output_description])

    class PublishingDataProxy(DataProxyBase, descriptors=descriptors):
        """Handler for write access to the `output` of an operation.

        Acts as a sort of PublisherCollection.
        """

    return PublishingDataProxy


class StaticSourceManager(object):
    """Provide the resource manager interface for local static data.

    Allow data transformations on the proxied resource.

    Keyword Args:
        proxied_data: A gmxapi supported data object.
        width: Size of (one-dimensional) shaped data produced by function.
        function: Transformation to perform on the managed data.

    The callable passed as ``function`` must accept a single argument. The
    argument will be an iterable when proxied_data represents an ensemble,
    or an object of the same type as proxied_data otherwise.
    """
    OutputDataProxyType = typing.TypeVar('OutputDataProxyType', bound=DataProxyBase)

    def __init__(self, *, name: str, proxied_data, width: int, function: typing.Callable):
        assert not isinstance(proxied_data, Future)
        if hasattr(proxied_data, 'width'):
            # Ensemble data source
            assert hasattr(proxied_data, 'source')
            self._result = function(proxied_data.source)
        else:
            self._result = function(proxied_data)
        if width > 1:
            if isinstance(self._result, (str, bytes)):
                # In this case, do not implicitly broadcast
                raise exceptions.ValueError('"function" produced data incompatible with "width".')
            else:
                if not isinstance(self._result, collections.abc.Iterable):
                    raise exceptions.DataShapeError(
                        'Expected iterable of size {} but "function" result is not iterable.')
            data = list(self._result)
            size = len(data)
            if len(data) != width:
                raise exceptions.DataShapeError(
                    'Expected iterable of size {} but "function" produced a {} of size {}'.format(width, type(data),
                                                                                                  size))
            dtype = type(data[0])
        else:
            if width != 1:
                raise exceptions.ValueError('width must be an integer 1 or greater.')
            dtype = type(self._result)
            if issubclass(dtype, (list, tuple)):
                dtype = datamodel.NDArray
                data = [datamodel.ndarray(self._result)]
            elif isinstance(self._result, collections.abc.Iterable):
                if not isinstance(self._result, (str, bytes)):
                    raise exceptions.ValueError(
                        'Expecting width 1 but "function" produced iterable type {}.'.format(type(self._result)))
                else:
                    dtype = str
                    data = [str(self._result)]
            else:
                data = [self._result]
        description = ResultDescription(dtype=dtype, width=width)
        self._data = OutputData(name=name, description=description)
        for member in range(width):
            self._data.set(data[member], member=member)

        output_collection_description = OutputCollectionDescription(**{name: dtype})
        self.output_data_proxy = define_output_data_proxy(output_description=output_collection_description)

    def is_done(self, name: str) -> bool:
        return True

    def get(self, name: str) -> 'OutputData':
        assert self._data.name == name
        return self._data

    def data(self) -> OutputDataProxyType:
        return self.output_data_proxy(self)

    def width(self) -> int:
        # TODO: It looks like the OutputData ResultDescription probably belongs
        #  in the public interface.
        return self._data._description.width

    def update_output(self):
        pass

    def reset(self):
        pass


class ProxyResourceManager(object):
    """Act as a resource manager for data managed by another resource manager.

    Allow data transformations on the proxied resource.

    Keyword Args:
        proxied_future: An object implementing the Future interface.
        width: Size of (one-dimensional) shaped data produced by function.
        function: Transformation to perform on the result of proxied_future.

    The callable passed as ``function`` must accept a single argument, which will
    be an iterable when proxied_future represents an ensemble, or an object of
    type proxied_future.description.dtype otherwise.
    """

    def __init__(self, proxied_future: 'Future', width: int, function: typing.Callable):
        self._done = False
        self._proxied_future = proxied_future
        self._width = width
        self.name = self._proxied_future.name
        self._result = None
        assert callable(function)
        self.function = function

    def width(self) -> int:
        return self._width

    def reset(self):
        self._done = False
        self._proxied_future._reset()
        self._result = None

    def is_done(self, name: str) -> bool:
        return self._done

    def get(self, name: str):
        if name != self.name:
            raise exceptions.ValueError('Request for unknown data.')
        if not self.is_done(name):
            raise exceptions.ProtocolError('Data not ready.')
        result = self.function(self._result)
        if self._width != 1:
            # TODO Fix this typing nightmare:
            #  ResultDescription should be fully knowable and defined when the resource manager is initialized.
            data = OutputData(name=self.name, description=ResultDescription(dtype=type(result[0]), width=self._width))
            for member, value in enumerate(result):
                data.set(value, member)
        else:
            data = OutputData(name=self.name, description=ResultDescription(dtype=type(result), width=self._width))
            data.set(result, 0)
        return data

    def update_output(self):
        self._result = self._proxied_future.result()
        self._done = True

    def data(self) -> DataProxyBase:
        raise exceptions.ApiError('ProxyResourceManager cannot yet manage a full OutputDataProxy.')


# TODO: Implement observer pattern for edge->node data flow.
# Step 0: implement subject interface subscribe()
# Step 1: implement subject interface get_state()
# Step 2: implement observer interface update()
# Step 3: implement subject interface notify()
# Step 4: implement observer hook to support substantial change in source that
#         invalidates downstream fingerprinting, such as a new subgraph iteration.
# class Observer(object):
#     """Abstract base class for data observers."""
#     def rebind(self, edge: DataEdge):
#         """Recreate the Operation at the consuming end of the DataEdge."""


class Future(object):
    """gmxapi data handle.

    Future is currently more than a Future right now. (should be corrected / clarified.)
    Future is also a facade to other features of the data provider.

    Future objects are the most atomic interface between Contexts. User scripts
    may hold Futures from which they extract data with result(). Operation output
    used as input for another Operation can be decomposed such that the Operation
    factory has only Future objects in its input.

    TODO: ``subscribe`` method allows consumers to bind as Observers.

    TODO: extract the abstract class for input inspection?
    Currently abstraction is handled through SourceResource subclassing.
    """

    def __init__(self, resource_manager, name: str, description: ResultDescription):
        self.name = name
        if not isinstance(description, ResultDescription):
            raise exceptions.ValueError('Need description of requested data.')
        self.resource_manager = resource_manager
        self.description = description

    def result(self):
        """Fetch data to the caller's Context.

        Returns an object of the concrete type specified according to
        the operation that produces this Result.

        Ensemble data are returned as a list. Scalar results or results from single
        member ensembles are returned as scalars.
        """
        self.resource_manager.update_output()
        assert self.resource_manager.is_done(self.name)
        # Return ownership of concrete data
        handle = self.resource_manager.get(self.name)

        # For intuitive use in non-ensemble cases, we represent data as bare scalars
        # when possible. It is easier for users to cast scalars to lists of length 1
        # than to introspect their own code to determine if a list of length 1 is
        # part of an ensemble or not. The data model will become clearer as we
        # develop more robust handling of multidimensional data and data flow topologies.
        # In the future,
        # we may distinguish between data of shape () and shape (1,), but we will need
        # to be careful with semantics. We are already starting to adopt a rule-of-thumb
        # that data objects assume the minimum dimensionality necessary unless told
        # otherwise, and we could make that a hard rule if it doesn't make other things
        # too difficult.
        if self.description.width == 1:
            return handle.data(member=0)
        else:
            return handle.data()

    def _reset(self):
        """Mark the Future "not done" to allow reexecution.

        Invalidates cached results, resets "done" markers in data sources, and
        triggers _reset recursively.

        Note: this is a hack that is inconsistent with the plan of unique mappings
        of inputs to outputs, but allows a quick prototype for looping operations.
        """
        self.resource_manager.reset()

    @property
    def dtype(self):
        return self.description.dtype

    def __getitem__(self, item):
        """Get a more limited view on the Future."""
        description = ResultDescription(dtype=self.dtype, width=self.description.width)
        # TODO: Use explicit typing when we have more thorough typing.
        description._dtype = None
        if self.description.width == 1:
            proxy = ProxyResourceManager(self,
                                         width=description.width,
                                         function=lambda value, key=item: value[key])
        else:
            proxy = ProxyResourceManager(self,
                                         width=description.width,
                                         function=lambda value, key=item:
                                         [subscriptable[key] for subscriptable in value])
        future = Future(proxy, self.name, description=description)
        return future


class OutputDescriptor(ProxyDataDescriptor):
    """Read-only data descriptor for proxied output access.

    Knows how to get a Future from the resource manager.
    """

    def __get__(self, proxy: DataProxyBase, owner):
        if proxy is None:
            # Access through class attribute of owner class
            return self
        result_description = ResultDescription(dtype=self._dtype, width=proxy.ensemble_width)
        return proxy._resource_instance.future(name=self._name, description=result_description)


def define_output_data_proxy(output_description: OutputCollectionDescription) -> typing.Type[DataProxyBase]:
    descriptors = {name: OutputDescriptor(name, description) for name, description in output_description.items()}

    class OutputDataProxy(DataProxyBase, descriptors=descriptors):
        """Handler for read access to the `output` member of an operation handle.

        Acts as a sort of ResultCollection.

        A ResourceManager creates an OutputDataProxy instance at initialization to
        provide the ``output`` property of an operation handle.
        """

    # Note: the OutputDataProxy has an inherent ensemble shape in the context
    # in which it is used, but that is an instance characteristic, not part of this type definition.
    # TODO: (FR5) The current tool does not support topology changing operations.
    return OutputDataProxy


# In the longer term, Contexts could provide metaclasses that allow transformation or dispatching
# of the basic aspects of the operation protocols between Contexts or from a result handle into a
# new context, based on some attribute or behavior in the result handle.


# Encapsulate the description of the input data flow.
PyFuncInput = collections.namedtuple('Input', ('args', 'kwargs', 'dependencies'))


class SinkTerminal(object):
    """Operation input end of a data edge.

    In addition to the information in an InputCollectionDescription, includes
    topological information for the Operation node (ensemble width).

    Collaborations: Required for creation of a DataEdge. Created with knowledge
    of a DataSourceCollection instance and a InputCollectionDescription.
    """

    # TODO: This clearly maps to a Builder pattern.
    # I think we want to get the sink terminal builder from a factory parameterized by InputCollectionDescription,
    # add data source collections, and then build the sink terminal for the data edge.
    def __init__(self, input_collection_description: InputCollectionDescription):
        """Define an appropriate data sink for a new operation node.

        Resolve data sources and input description to determine connectability,
        topology, and any necessary implicit data transformations.

        :param data_source_collection: Collection of offered input data.
        :param input_collection_description: Available inputs for Operation
        :return: Fully formed description of the Sink terminal for a data edge to be created.

        Collaborations: Execution Context implementation.
        """
        self.ensemble_width = 1
        self.inputs = input_collection_description

    def update_width(self, width: int):
        if not isinstance(width, int):
            try:
                width = int(width)
            except TypeError:
                raise exceptions.TypeError('Need an integer width > 0.')
        if width < 1:
            raise exceptions.ValueError('Nonsensical ensemble width: {}'.format(int(width)))
        if self.ensemble_width != 1:
            if width != self.ensemble_width:
                raise exceptions.ValueError(
                    'Cannot change ensemble width {} to width {}.'.format(self.ensemble_width, width))
        self.ensemble_width = width

    def update(self, data_source_collection: DataSourceCollection):
        """Update the SinkTerminal with the proposed data provider."""
        for name, sink_dtype in self.inputs.items():
            if name not in data_source_collection:
                # If/when we accept data from multiple sources, we'll need some additional sanity checking.
                if not hasattr(self.inputs.signature.parameters[name], 'default'):
                    raise exceptions.UsageError('No data or default for {}'.format(name))
            else:
                # With a single data source, we need data to be in the source or have a default
                assert name in data_source_collection
                assert issubclass(sink_dtype, valid_result_types)
                source = data_source_collection[name]
                if isinstance(source, sink_dtype):
                    continue
                else:
                    if isinstance(source, collections.abc.Iterable) and not isinstance(source, (
                            str, bytes, collections.abc.Mapping)):
                        assert isinstance(source, datamodel.NDArray)
                        if sink_dtype != datamodel.NDArray:
                            # Implicit scatter
                            self.update_width(len(source))


class DataEdge(object):
    """State and description of a data flow edge.

    A DataEdge connects a data source collection to a data sink. A sink is an
    input or collection of inputs of an operation (or fused operation). An operation's
    inputs may be fed from multiple data source collections, but an operation
    cannot be fully instantiated until all of its inputs are bound, so the DataEdge
    is instantiated at the same time the operation is instantiated because the
    required topology of a graph edge may be determined by the required topology
    of another graph edge.

    A data edge has a well-defined topology only when it is terminated by both
    a source and sink. Creation requires that a source collection is compared to
    a sink description.

    Calling code initiates edge creation by passing well-described data sources
    to an operation factory. The data sources may be annotated with explicit scatter
    or gather commands.

    The resource manager for the new operation determines the
    required shape of the sink to handle all of the offered input.

    Broadcasting
    and transformations of the data sources are then determined and the edge is
    established.

    At that point, the fingerprint of the input data at each operation
    becomes available to the resource manager for the operation. The fingerprint
    has sufficient information for the resource manager of the operation to
    request and receive data through the execution context.

    Instantiating operations and data edges implicitly involves collaboration with
    a Context instance. The state of a given Context or the availability of a
    default Context through a module function may affect the ability to instantiate
    an operation or edge. In other words, behavior may be different for connections
    being made in the scripting environment versus the running Session, and implementation
    details can determine whether or not new operations or data flow can occur in
    different code environments.
    """

    class ConstantResolver(object):
        def __init__(self, value):
            self.value = value

        def __call__(self, member=None):
            return self.value

    def __init__(self, source_collection: DataSourceCollection, sink_terminal: SinkTerminal):
        # Adapters are callables that transform a source and node ID to local data.
        # Every key in the sink has an adapter.
        self.adapters = {}
        self.source_collection = source_collection
        self.sink_terminal = sink_terminal
        for name in sink_terminal.inputs:
            if name not in source_collection:
                if hasattr(sink_terminal.inputs[name], 'default'):
                    self.adapters[name] = self.ConstantResolver(sink_terminal.inputs[name])
                else:
                    # TODO: Initialize with multiple DataSourceCollections?
                    raise exceptions.ValueError('No source or default for required input "{}".'.format(name))
            else:
                source = source_collection[name]
                sink = sink_terminal.inputs[name]
                if isinstance(source, (str, bool, int, float, dict)):
                    if issubclass(sink, (str, bool, int, float, dict)):
                        self.adapters[name] = self.ConstantResolver(source)
                    else:
                        assert issubclass(sink, datamodel.NDArray)
                        self.adapters[name] = self.ConstantResolver(datamodel.ndarray([source]))
                elif isinstance(source, datamodel.NDArray):
                    if issubclass(sink, datamodel.NDArray):
                        # TODO: shape checking
                        # Implicit broadcast may not be what is intended
                        self.adapters[name] = self.ConstantResolver(source)
                    else:
                        if source.shape[0] != sink_terminal.ensemble_width:
                            raise exceptions.ValueError(
                                'Implicit broadcast could not match array source to ensemble sink')
                        else:
                            self.adapters[name] = lambda member, source=source: source[member]
                elif hasattr(source, 'result'):
                    # Handle data futures...
                    # If the Future is part of an ensemble, result() will return a list.
                    # Otherwise, it will return a single object.
                    ensemble_width = source.description.width
                    if ensemble_width == 1:
                        self.adapters[name] = lambda member, source=source: source.result()
                    else:
                        self.adapters[name] = lambda member, source=source: source.result()[member]
                else:
                    assert isinstance(source, EnsembleDataSource)
                    self.adapters[name] = lambda member, source=source: source.node(member)

    def reset(self):
        self.source_collection.reset()

    def resolve(self, key: str, member: int):
        return self.adapters[key](member=member)

    def sink(self, node: int) -> dict:
        """Consume data for the specified sink terminal node.

        Run-time utility delivers data from the bound data source(s) for the
        specified terminal that was configured when the edge was created.

        Terminal node is identified by a member index number.

        Returns:
            A Python dictionary of the provided inputs as local data (not Future).
        """
        results = {}
        sink_ports = self.sink_terminal.inputs
        for key in sink_ports:
            results[key] = self.resolve(key, node)
        return results


class ResourceManager(object):
    """Provides data publication and subscription services.

        Owns the data published by the operation implementation or served to consumers.
        Mediates read and write access to the managed data streams.

        This ResourceManager implementation is defined in conjunction with a
        run-time definition of an Operation that wraps a Python callable (function).
        ResourceManager is instantiated with a reference to the callable.

        When the Operation is run, the resource manager prepares resources for the wrapped
        function. Inputs provided to the Operation factory are provided to the
        function as keyword arguments. The wrapped function publishes its output
        through the (additional) ``output`` key word argument. This argument is
        a short-lived resource, prepared by the ResourceManager, with writable
        attributes named in the call to function_wrapper().

        After the Operation has run and the outputs published, the data managed
        by the ResourceManager is marked "done."

        Protocols:

        The data() method produces a read-only collection of outputs named for
        the Operation when the Operation's ``output`` attribute is accessed.

        publishing_resources() can be called once during the ResourceManager lifetime
        to provide the ``output`` object for the wrapped function. (Used by update_output().)

        update_output() brings the managed output data up-to-date with the input
        when the Operation results are needed. If the Operation has not run, an
        execution session is prepared with input and output arguments for the
        wrapped Python callable. Output is publishable only during this session.

    TODO: This functionality should evolve to be a facet of Context implementations.
     There should be no more than one ResourceManager instance per work graph
     node in a Context. This will soon be at odds with letting the ResourceManager
     be owned by an operation instance handle.
    TODO: The publisher and data objects can be more strongly defined through
     interaction between the Context and clients.

    Design notes:

    The normative pattern for updating data is to execute a node in the work
    graph, passing Resources for an execution Session to an operation runner.
    The resources and runner are dependent on the implementation details of
    the operation and the execution context, so logical execution may look
    like the following.

        resource_builder = ResourcesBuilder()
        runner_builder = RunnerBuilder()
        input_resource_director = input_resource_factory.director(input)
        output_resource_director = publishing_resource_factory.director(output)
        input_resource_director(resource_builder, runner_builder)
        output_resource_director(resource_builder, runner_builder)
        resources = resource_builder.build()
        runner = runner_builder.build()
        runner(resources)

    Only the final line is intended to be literal. The preceding code, if it
    exists in entirety, may be spread across several code comments.

    TODO: Data should be pushed, not pulled.
    Early implementations executed operation code and extracted results directly.
    While we need to be able to "wait for" results and alert the data provider that
    we are ready for input, we want to defer execution management and data flow to
    the framework.
    """

    @contextmanager
    def __publishing_context(self, ensemble_member=0):
        """Get a context manager for resolving the data dependencies of this node.

        The returned object is a Python context manager (used to open a `with` block)
        to define the scope in which the operation's output can be published.
        'output' type resources can be published exactly once, and only while the
        publishing context is active. (See operation.function_wrapper())

        Used internally to implement ResourceManager.publishing_resources()

        Responsibilities of the context manager are to:
            * (TODO) Make sure dependencies are resolved.
            * Make sure outputs are marked 'done' when leaving the context.

        """

        # TODO:
        # if self._data.done():
        #     raise exceptions.ProtocolError('Resources have already been published.')

        # I don't think we want the OperationDetails to need to know about ensemble data,
        # (though the should probably be allowed to), so we may need a separate interface
        # for the resource manager with built-in scope-limiting to a single ensemble member.
        # Right now, one Operation handle owns one ResourceManager (which takes care of
        # the ensemble details), which owns one OperationDetails (which has no ensemble knowledge).
        # It is the responsibility of the calling code to make sure the PublishingDataProxy
        # gets used correctly.

        # ref: https://docs.python.org/3/library/contextlib.html#contextlib.contextmanager
        try:
            if not self._done[ensemble_member]:
                resource = self._operation.publishing_data_proxy(instance=weakref.proxy(self),
                                                                 client_id=ensemble_member)
                yield resource
        except Exception as e:
            message = 'Uncaught {} while providing output-publishing resources for {}.'.format(repr(e), self._operation)
            raise exceptions.ApiError(message) from e
        finally:
            self._done[ensemble_member] = True

    def __init__(self, *, source: DataEdge, operation):
        """Initialize a resource manager for the inputs and outputs of an operation.

        Arguments:
            operation : implementation details for a Python callable
            input_fingerprint : Uniquely identifiable input data description

        """
        # Note: This implementation assumes there is one ResourceManager instance per data source,
        # so we only stash the inputs and dependency information for a single set of resources.
        # TODO: validate input_fingerprint as its interface becomes clear.
        self._input_edge = source
        self.ensemble_width = self._input_edge.sink_terminal.ensemble_width
        self._operation = operation

        self._data = self._operation.make_datastore(self.ensemble_width)

        # We store a rereference to the publishing context manager implementation
        # in a data structure that can only produce one per Python interpreter
        # (using list.pop()).
        # TODO: reimplement as a data descriptor
        #  so that PublishingDataProxy does not need a bound circular reference.
        self.__publishing_resources = [self.__publishing_context]

        self._done = [False] * self.ensemble_width
        self.__operation_entrance_counter = 0

    def width(self) -> int:
        return self.ensemble_width

    def reset(self):
        self.__operation_entrance_counter = 0
        self._done = [False] * self.ensemble_width
        self.__publishing_resources = [self.__publishing_context]
        for data in self._data.values():
            data.reset()
        self._input_edge.reset()
        assert self.__operation_entrance_counter == 0

    def done(self, member=None):
        if member is None:
            return all(self._done)
        else:
            return self._done[member]

    def set_result(self, name, value, member: int):
        try:
            for item in value:
                # In this specification, it is antithetical to publish Futures.
                if hasattr(item, 'result'):
                    raise exceptions.ApiError('Operation produced Future instead of real output.')
        except TypeError:
            # Ignore when `item` is not iterable.
            pass
        self._data[name].set(value=value, member=member)

    def is_done(self, name):
        return self._data[name].done

    def get(self, name: str):
        """

        Raises exceptions.ProtocolError if requested data is not local yet.
        Raises exceptions.ValueError if data is requested for an unknown name.
        """
        if name not in self._data:
            raise exceptions.ValueError('Request for unknown data.')
        if not self.is_done(name):
            raise exceptions.ProtocolError('Data not ready.')
        assert isinstance(self._data[name], OutputData)
        return self._data[name]

    def update_output(self):
        """Bring the output of the bound operation up to date.

        Execute the bound operation once if and only if it has not
        yet been run in the lifetime of this resource manager.

        Used internally to implement Futures for the local operation
        associated with this resource manager.

        TODO: We need a different implementation for an operation whose output
         is served by multiple resource managers. E.g. an operation whose output
         is available across the ensemble, but which should only be executed on
         a single ensemble member.
        """
        # This code is not intended to be reentrant. We make a modest attempt to
        # catch unexpected reentrance, but this is not (yet) intended to be a thread-safe
        # resource manager implementation.
        # TODO: Handle checking just the ensemble members this resource manager is responsible for.
        # TODO: Replace with a managed observer pattern. Update once when input is available in the Context.
        if not self.done():
            self.__operation_entrance_counter += 1
            if self.__operation_entrance_counter > 1:
                raise exceptions.ProtocolError('Bug detected: resource manager tried to execute operation twice.')
            if not self.done():
                # Note! This is a detail of the ResourceManager in a SerialContext
                for i in range(self.ensemble_width):
                    with self.local_input(i) as input:
                        # Note: Resources are marked "done" by the publishing system
                        # before the following context manager finishes exiting.
                        with self.publishing_resources()(ensemble_member=i) as output:
                            # self._runner(*input.args, output=output, **input.kwargs)
                            ####
                            # Here we can make _runner a thing that accepts session resources, and
                            # is created by specializable builders. Separate out the expression of
                            # inputs.
                            #
                            # resource_builder = OperationDetails.ResourcesBuilder(context)
                            # runner_builder = OperationDetails.RunnerBuilder(context)
                            # input_resource_director = self._input_resource_factory.director(input)
                            # output_resource_director = self._publishing_resource_factory.director(output)
                            # input_resource_director(resource_builder, runner_builder)
                            # output_resource_director(resource_builder, runner_builder)
                            # resources = resource_builder.build()
                            # runner = runner_builder.build()
                            # runner(resources)
                            resources = self._operation.resource_director(input=input, output=output)
                            self._operation(resources)

    def future(self, name: str, description: ResultDescription):
        """Retrieve a Future for a named output.

        Provide a description of the expected result to check for compatibility or
        implicit topological conversion.

        TODO: (FR5+) Normalize this part of the interface between operation definitions and
         resource managers.
        """
        if not isinstance(name, str) or name not in self._data:
            raise exceptions.ValueError('"name" argument must name an output.')
        assert description is not None
        requested_dtype = description.dtype
        available_dtype = self._data[name]._description.dtype
        if requested_dtype != available_dtype:
            # TODO: framework to check for implicit conversions
            message = 'Requested Future of type {} is not compatible with available type {}.'
            message = message.format(requested_dtype, available_dtype)
            raise exceptions.ApiError(message)
        return Future(self, name, description)

    def data(self) -> DataProxyBase:
        """Get an adapter to the output resources to access results."""
        return self._operation.output_data_proxy(self)

    @contextmanager
    def local_input(self, member: int = None):
        """In an API session, get a handle to fully resolved locally available input data.

        Execution dependencies are resolved on creation of the context manager. Input data
        becomes available in the ``as`` object when entering the context manager, which
        becomes invalid after exiting the context manager. Resources allocated to hold the
        input data may be released when exiting the context manager.

        It is left as an implementation detail whether the context manager is reusable and
        under what circumstances one may be obtained.
        """
        # Localize data
        kwargs = self._input_edge.sink(node=member)
        assert 'input' not in kwargs

        # Check that we have real data
        for key, value in kwargs.items():
            assert not hasattr(value, 'result')
            assert not hasattr(value, 'run')
            value_list = []
            if isinstance(value, list):
                value_list = value
            if isinstance(value, datamodel.NDArray):
                value_list = value._values
            if isinstance(value, collections.abc.Mapping):
                value_list = value.values()
            assert not isinstance(value_list, (Future))
            assert not hasattr(value_list, 'result')
            assert not hasattr(value_list, 'run')
            for item in value_list:
                assert not hasattr(item, 'result')

        input_pack = collections.namedtuple('InputPack', ('kwargs'))(kwargs)

        # Prepare input data structure
        # Note: we use 'yield' instead of 'return' for the protocol expected by
        # the @contextmanager decorator
        yield input_pack

    def publishing_resources(self):
        """Get a context manager for resolving the data dependencies of this node.

        Use the returned object as a Python context manager.
        'output' type resources can be published exactly once, and only while the
        publishing context is active.

        Write access to publishing resources can be granted exactly once during the
        resource manager lifetime and conveys exclusive access.
        """
        return self.__publishing_resources.pop()


class PyFunctionRunnerResources(collections.UserDict):
    """Runtime resources for Python functions.

    Produced by a ResourceDirector for a particular Operation.
    """

    def output(self):
        if 'output' in self:
            return self['output']
        else:
            return None

    def input(self):
        return {key: value for key, value in self.items() if key != 'output'}


class PyFunctionRunner(abc.ABC):
    def __init__(self, *, function: typing.Callable, output_description: OutputCollectionDescription):
        assert callable(function)
        self.function = function
        self.output_description = output_description

    @abc.abstractmethod
    def __call__(self, resources: PyFunctionRunnerResources):
        self.function(output=resources.output(), **resources.input())


class CapturedOutputRunner(PyFunctionRunner):
    """Function runner that captures return value as output.data"""

    def __call__(self, resources: PyFunctionRunnerResources):
        resources['output'].data = self.function(**resources.input())


class OutputParameterRunner(PyFunctionRunner):
    """Function runner that uses output parameter to let function publish output."""

    def __call__(self, resources: PyFunctionRunnerResources):
        self.function(**resources)


def wrapped_function_runner(function, output_description: OutputCollectionDescription = None) -> PyFunctionRunner:
    """Get an adapter for a function to be wrapped.

    If the function does not accept a publishing data proxy as an `output`
    key word argument, the returned object has a `capture_output` attribute that
    must be re-assigned by the calling code before calling the runner. `capture_output`
    must be assigned to be a callable that will receive the output of the wrapped
    function.

    Returns:
        Callable with a signature `__call__(*args, **kwargs)` and no return value

    Collaborations:
        OperationDetails.resource_director assigns the `capture_output` member of the returned object.
    """
    assert callable(function)
    signature = inspect.signature(function)

    # Implementation note: this function dispatches an implementation with the
    # logic below. A better factoring would be a "chain of responsibility" in
    # which the concrete Runners would be tried in sequence and determine internally
    # whether to create a runner, raise an error, or defer.

    # Determine output details for proper dispatching.
    # First check for signature with output parameter.
    # TODO FR4: standardize typing
    if 'output' in signature.parameters:
        if not isinstance(output_description, OutputCollectionDescription):
            if not isinstance(output_description, collections.abc.Mapping):
                raise exceptions.UsageError(
                    'Function passes output through call argument, but output is not described.')
            return OutputParameterRunner(
                function=function,
                output_description=OutputCollectionDescription(**output_description))
        else:
            return OutputParameterRunner(function=function,
                                         output_description=output_description)
    # Next try output_description parameter or function return annotation.
    else:
        if isinstance(output_description, OutputCollectionDescription):
            return_type = output_description['data'].gmxapi_datatype
        elif output_description is not None:
            # output_description should be None for inferred output or
            # a singular mapping of the key 'data' to a gmxapi type.
            if not isinstance(output_description, collections.abc.Mapping) \
                    or set(output_description.keys()) != {'data'}:
                raise exceptions.ApiError(
                    'invalid output description for wrapped function: {}'.format(output_description))
            if signature.return_annotation != signature.empty:
                if signature.return_annotation != output_description['data']:
                    raise exceptions.ApiError(
                        'Wrapped function with return-value-capture provided with non-matching output description.')
            return_type = output_description['data']
        else:
            # Use return type inferred from function signature.
            return_type = signature.return_annotation
        if return_type == signature.empty or return_type is None:
            raise exceptions.ApiError('No return annotation or output_description for {}'.format(function))
        return CapturedOutputRunner(function=function,
                                    output_description=OutputCollectionDescription(data=return_type))


class OperationHandle(object):
    """Dynamically defined Operation handle.

    Define a gmxapi Operation for the functionality being wrapped by the enclosing code.

    An Operation type definition encapsulates description of allowed inputs
    of an Operation. An Operation instance represents a node in a work graph
    with uniquely fingerprinted inputs and well-defined output. The implementation
    of the operation is a collaboration with the resource managers resolving
    data flow for output Futures, which may depend on the execution context.
    """

    def __init__(self, resource_manager):
        """Initialization defines the unique input requirements of a work graph node.

        Initialization parameters map to the parameters of the wrapped function with
        addition(s) to support gmxapi data flow and deferred execution.

        If provided, an ``input`` keyword argument is interpreted as a parameter pack
        of base input. Inputs also present as standalone keyword arguments override
        values in ``input``.

        Inputs that are handles to gmxapi operations or outputs induce data flow
        dependencies that the framework promises to satisfy before the Operation
        executes and produces output.
        """
        # TODO: When the resource manager can be kept alive by an enclosing or
        #  module-level Context, convert to a weakref.
        self.__resource_manager = resource_manager
        # The unique identifier for the operation node allows the Context implementation
        # to manage the state of the handle. Reproducibility of node_uid is TBD, but
        # it must be unique in a Context where it references a different operation node.
        self.node_uid = None

    @property
    def output(self) -> DataProxyBase:
        # TODO: We can configure `output` as a data descriptor
        #  instead of a property so that we can get more information
        #  from the class attribute before creating an instance of OperationDetails.OutputDataProxy.
        # The C++ equivalence would probably be a templated free function for examining traits.
        return self.__resource_manager.data()

    def run(self):
        """Make a single attempt to resolve data flow conditions.

        This is a public method, but should not need to be called by users. Instead,
        just use the `output` data proxy for result handles, or force data flow to be
        resolved with the `result` methods on the result handles.

        `run()` may be useful to try to trigger computation (such as for remotely
        dispatched work) without retrieving results locally right away.

        `run()` is also useful internally as a facade to the Context implementation details
        that allow `result()` calls to ask for their data dependencies to be resolved.
        Typically, `run()` will cause results to be published to subscribing operations as
        they are calculated, so the `run()` hook allows execution dependency to be slightly
        decoupled from data dependency, as well as to allow some optimizations or to allow
        data flow to be resolved opportunistically. `result()` should not call `run()`
        directly, but should cause the resource manager / Context implementation to process
        the data flow graph.

        In one conception, `run()` can have a return value that supports control flow
        by itself being either runnable or not. The idea would be to support
        fault tolerance, implementations that require multiple iterations / triggers
        to complete, or looping operations.
        """
        self.__resource_manager.update_output()


# TODO: For outputs, distinguish between "results" and "events".
#  Both are published to the resource manager in the same way, but the relationship
#  with subscribers is potentially different.
def function_wrapper(output: dict = None):
    """Generate a decorator for wrapped functions with signature manipulation.

    New function accepts the same arguments, with additional arguments required by
    the API.

    The new function returns an object with an `output` attribute containing the named outputs.

    Example:
        @function_wrapper(output={'spam': str, 'foo': str})
        def myfunc(parameter: str = None, output=None):
            output.spam = parameter
            output.foo = parameter + ' ' + parameter

        operation1 = myfunc(parameter='spam spam')
        assert operation1.output.spam.result() == 'spam spam'
        assert operation1.output.foo.result() == 'spam spam spam spam'

    If 'output' is provided to the wrapper, a data structure will be passed to
    the wrapped functions with the named attributes so that the function can easily
    publish multiple named results. Otherwise, the `output` of the generated operation
    will just capture the return value of the wrapped function.
    """
    if output is not None and not isinstance(output, collections.abc.Mapping):
        raise exceptions.TypeError('If provided, `output` argument must be a mapping of data names to types.')

    # TODO: (FR5+) gmxapi operations need to allow a context-dependent way to generate an implementation with input.
    # This function wrapper reproduces the wrapped function's kwargs, but does not allow chaining a
    # dynamic `input` kwarg and does not dispatch according to a `context` kwarg. We should allow
    # a default implementation and registration of alternate implementations. We don't have to do that
    # with functools.singledispatch, but we could, if we add yet another layer to generate a wrapper
    # that takes the context as the first argument. (`singledispatch` inspects the first argument rather
    # that a named argument)

    # Implementation note: The closure of the current function is used to
    # dynamically define several classes that support the operation to be
    # created by the returned decorator.

    def decorator(function) -> typing.Callable:
        # Explicitly capture `function` and `output` references.
        provided_output_map = output

        # Note: Allow operations to be defined entirely in template headers to facilitate
        # compile-time optimization of fused operations. Consider what distinction, if any,
        # exists between a fused operation and a more basic operation. Probably it amounts
        # to aspects related to interaction with the Context that get combined in a fused
        # operation, such as the resource director, builder, etc.
        # Note that a ResourceManager holds a reference to an instance of OperationDetails,
        # though input signature needs to be well-defined at the type level.
        class OperationDetails(object):
            # Warning: function.__qualname__ is not rigorous since function may be in a local scope.
            # TODO: Improve base identifier.
            # Suggest registering directly in the Context instead of in this local class definition.
            __basename = function.__qualname__
            __last_uid = 0
            _input_signature_description = InputCollectionDescription.from_function(function)
            # TODO: Separate the class and instance logic for the runner.
            # Logically, the runner is a detail of a context-specific implementation class,
            # though the output is not generally fully knowable until an instance is initialized
            # for a certain input fingerprint.
            # Note: We are almost at a point where this class can be subsumed into two
            # possible return types for wrapped_function_runner, acting as an operation helper.
            _runner = wrapped_function_runner(function, provided_output_map)
            _output_description = _runner.output_description
            _output_data_proxy_type = define_output_data_proxy(_output_description)
            _publishing_data_proxy_type = define_publishing_data_proxy(_output_description)

            # TODO: This is a Context detail.
            def make_datastore(self, ensemble_width: int):
                """Create a container to hold the resources for this operation node.

                Used internally by the resource manager when setting up the node.
                Evolution of the C++ framework for creating the Operation SessionResources
                object will inform the future of this and the resource_director method.
                """
                datastore = {}
                for name, dtype in self.output_description().items():
                    assert isinstance(dtype, type)
                    result_description = ResultDescription(dtype, width=ensemble_width)
                    datastore[name] = OutputData(name=name, description=result_description)
                return datastore

            @classmethod
            def signature(cls) -> InputCollectionDescription:
                """Mapping of named inputs and input type.

                Used to determine valid inputs before an Operation node is created.
                """
                return cls._input_signature_description

            def output_description(self) -> OutputCollectionDescription:
                """Mapping of available outputs and types for an existing Operation node."""
                return self._output_description

            def publishing_data_proxy(self, *, instance, client_id: int):
                """Factory for Operation output publishing resources.

                Used internally when the operation is run with resources provided by instance."""
                assert isinstance(instance, ResourceManager)
                return self._publishing_data_proxy_type(instance=instance, client_id=client_id)

            def output_data_proxy(self, instance):
                """Get an object that can provide Futures for output data managed by instance."""
                assert isinstance(instance, ResourceManager)
                return self._output_data_proxy_type(instance=instance)

            def __call__(self, resources: PyFunctionRunnerResources):
                """Execute the operation with provided resources.

                Resources are prepared in an execution context with aid of resource_director()

                After the first call, output data has been published and is trivially
                available through the output_data_proxy()
                """
                self._runner(resources)

            @classmethod
            def resource_director(cls, *, input=None, output=None) -> PyFunctionRunnerResources:
                """a Director factory that helps build the Session Resources for the function.

                The Session launcher provides the director with all of the resources previously
                requested/negotiated/registered by the Operation. The director uses details of
                the operation to build the resources object required by the operation runner.

                For the Python Context, the protocol is for the Context to call the
                resource_director instance method, passing input and output containers.
                """
                resources = PyFunctionRunnerResources()
                resources.update(input.kwargs)
                resources.update({'output': output})

                # TODO: Remove this hack when we can better handle Futures of Containers and Future slicing.
                for name in resources:
                    if isinstance(resources[name], (list, tuple)):
                        resources[name] = datamodel.ndarray(resources[name])

                # Check data compatibility
                for name, value in resources.items():
                    if name != 'output':
                        expected = cls.signature()[name]
                        got = type(value)
                        if got != expected:
                            raise exceptions.TypeError('Expected {} but got {}.'.format(expected, got))
                return resources

            @classmethod
            def operation_director(cls, *args, context, label=None, **kwargs):
                """Dispatching Director for adding a work node.

                A Director for input of a particular sort knows how to reconcile
                input with the requirements of the Operation and Context node builder.
                The Director (using a less flexible / more standard interface)
                builds the operation node using a node builder provided by the Context.

                This is essentially the creation method, instead of __init__, but the
                object is created and owned by the framework, and the caller receives
                an OperationHandle instead of a reference to an instance of cls.

                # TODO: We need a way to compose this functionality for arbitrary Contexts.
                # That likely requires traits on the Contexts, and registration of Context
                # implementations. It seems likely that an Operation will register Director
                # implementations on import, and dispatching will be moved to the Context
                # implementations, which can either find an appropriate OperationDirector
                # or raise a compatibility error. To avoid requirements on import order of
                # Operations and Context implementations, we can change this to a non-abstract
                # dispatching method, requiring registration in the global gmxapi.context
                # module, or get rid of this method and use something like pkg_resources
                # "entry point" groups for independent registration of Directors and Contexts,
                # each annotated with relevant traits. E.g.:
                # https://setuptools.readthedocs.io/en/latest/setuptools.html#dynamic-discovery-of-services-and-plugins
                """
                data_source_collection = cls.signature().bind(*args, **kwargs)
                input_sink = SinkTerminal(cls.signature())
                input_sink.update(data_source_collection)
                edge = DataEdge(data_source_collection, input_sink)

                manager = ResourceManager(source=edge, operation=cls())
                handle = OperationHandle(resource_manager=manager)
                return handle


    # TODO: (FR4) Update annotations with gmxapi data types. E.g. return -> Future.
        @functools.wraps(function)
        def helper(*args, context=None, **kwargs):
            # Description of the Operation input (and output) occurs in the
            # decorator closure. By the time this factory is (dynamically) defined,
            # the OperationDetails and ResourceManager are well defined, but not
            # yet instantiated.
            # Inspection of the offered input occurs when this factory is called,
            # and OperationDetails, ResourceManager, and Operation are instantiated.

            # This calls a dispatching function that may not be able to reconcile the input
            # and Context capabilities. This is the place to handle various exceptions for
            # whatever reasons this reconciliation cannot occur.
            handle = OperationDetails.operation_director(*args, context=context, label=None, **kwargs)

            # TODO: NOW: The input fingerprint describes the provided input
            # as (a) ensemble input, (b) static, (c) future. By the time the
            # operation is instantiated, the topology of the node is known.
            # When compared to the InputCollectionDescription, the data compatibility
            # can be determined.

            return handle

        # to do: The factory itself needs to be able to register a factory with
        # the context that will be responsible for the Operation handle.
        # The factories need to be able to serve as dispatchers for themselves,
        # since an operation in one context may need to be reconstituted in a
        # different context.
        # The dispatching factory produces a director for a Context,
        # which will register a factory with the operation in that context.

        # The factory function has a DirectorFactory. Director instances talk to a NodeBuilder for a Context to
        # get handles to new operation nodes managed by the context. Part of that process includes registering
        # a DirectorFactory with the Context.
        return helper

    return decorator


@computed_result
def join_arrays(*, front: datamodel.NDArray = (), back: datamodel.NDArray = ()) -> datamodel.NDArray:
    """Operation that consumes two sequences and produces a concatenated single sequence.

    Note that the exact signature of the operation is not determined until this
    helper is called. Helper functions may dispatch to factories for different
    operations based on the inputs. In this case, the dtype and shape of the
    inputs determines dtype and shape of the output. An operation instance must
    have strongly typed output, but the input must be strongly typed on an
    object definition so that a Context can make runtime decisions about
    dispatching work and data before instantiating.
    # TODO: elaborate and clarify.
    # TODO: check type and shape.
    # TODO: figure out a better annotation.
    """
    # TODO: (FR4) Returned list should be an NDArray.
    if isinstance(front, (str, bytes)) or isinstance(back, (str, bytes)):
        raise exceptions.ValueError('Input must be a pair of lists.')
    assert isinstance(front, datamodel.NDArray)
    assert isinstance(back, datamodel.NDArray)
    new_list = list(front._values)
    new_list.extend(back._values)
    return new_list


Scalar = typing.TypeVar('Scalar')


def concatenate_lists(sublists: list = ()) -> Future:
    """Combine data sources into a single list.

    A trivial data flow restructuring operation.
    """
    if isinstance(sublists, (str, bytes)):
        raise exceptions.ValueError('Input must be a list of lists.')
    if len(sublists) == 0:
        return datamodel.ndarray([])
    else:
        return join_arrays(front=sublists[0], back=concatenate_lists(sublists[1:]))


def make_constant(value: Scalar) -> Future:
    """Provide a predetermined value at run time.

    This is a trivial operation that provides a (typed) value, primarily for
    internally use to manage gmxapi data flow.

    Accepts a value of any type. The object returned has a definite type and
    provides same interface as other gmxapi outputs. Additional constraints or
    guarantees on data type may appear in future versions.
    """
    dtype = type(value)
    source = StaticSourceManager(name='data', proxied_data=value, width=1, function=lambda x: x)
    description = ResultDescription(dtype=dtype, width=1)
    future = Future(source, 'data', description=description)
    return future
