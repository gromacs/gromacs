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

__all__ = [
    "computed_result",
    "function_wrapper",
]

import abc
import collections.abc
import contextlib
import functools
import inspect
import typing
import warnings
import weakref
from contextlib import contextmanager
from copy import copy

import gmxapi as gmx
import gmxapi.exceptions
from gmxapi import datamodel
from gmxapi import exceptions
from gmxapi import logger as root_logger
from gmxapi._transform import broadcast
from gmxapi._transform import identity
from gmxapi.abc import Future as FutureABC
from gmxapi.abc import MutableResource
from gmxapi.abc import NDArray
from gmxapi.abc import Node
from gmxapi.abc import OperationImplementation
from gmxapi.abc import OperationReference as AbstractOperationReference

# TODO: Relabel FutureABC as AbstractFuture, unless we can consolidate all gmxapi
#  Future types
#  into a single Protocol.
from gmxapi.datamodel import ArrayFuture
from gmxapi.exceptions import ApiError
from gmxapi.exceptions import DataShapeError
from gmxapi.exceptions import ProtocolError
from gmxapi.exceptions import TypeError as GmxapiTypeError
from gmxapi.exceptions import UsageError
from gmxapi.typing import _Context
from gmxapi.typing import Future as GenericFuture
from gmxapi.typing import ResultTypeVar
from gmxapi.typing import SourceTypeVar
from gmxapi.typing import valid_result_types
from gmxapi.typing import valid_source_types

# Initialize module-level logger
logger = root_logger.getChild("operation")
logger.info("Importing {}".format(__name__))

try:
    from mpi4py import MPI

    # TODO: Allow user to specify communicator.
    comm = MPI.COMM_WORLD
    rank_number = comm.Get_rank()
    comm_size = comm.Get_size()
except ImportError:
    comm = None
    rank_number = 0
    comm_size = 1
    MPI = None


class ResultDescription(typing.Generic[ResultTypeVar]):
    """Describe what will be returned when ``result()`` is called.

    Warning:
         *dtype* may be ``None`` if return type cannot be determined.
    """

    def __init__(self, dtype: typing.Type[ResultTypeVar] = None, width=1):
        assert isinstance(width, int)
        self._width = width

        if dtype is not None and not (
            isinstance(dtype, type) and issubclass(dtype, valid_result_types)
        ):
            raise ApiError(f"dtype {repr(dtype)} is not a valid result type.")
        self._dtype = typing.cast(typing.Type[ResultTypeVar], dtype)

        typename = getattr(
            self._dtype,
            "__qualname__",
            getattr(self._dtype, "__name__", repr(self._dtype)),
        )
        self.__representation = (
            f"{self.__class__.__name__}(dtype={typename}, width={width})"
        )

    @property
    def dtype(self) -> typing.Type[ResultTypeVar]:
        """node output type"""
        return self._dtype

    @property
    def width(self) -> int:
        """ensemble width"""
        return self._width

    def __repr__(self):
        return self.__representation

    def __eq__(self, other):
        return (
            self._dtype == typing.cast(ResultDescription, other).dtype
            and self._width == typing.cast(ResultDescription, other).width
        )


class OutputData(typing.Generic[ResultTypeVar]):
    """Encapsulate the description and storage of a data output.

    Allows a single signal to be issued when the object is "done".
    """

    _instance_count = 0
    _member_done: typing.List[bool]

    def __init__(
        self,
        name: str,
        description: ResultDescription[ResultTypeVar],
        done_callback: typing.Optional[typing.Callable[["OutputData"], None]] = None,
    ):
        assert name != ""
        self._name = name
        assert isinstance(description, ResultDescription)
        self._description = description
        # Warning: do not confuse OutputData._done with ResourceManager._done.
        self._member_done = [False] * self._description.width
        self._data = [None] * self._description.width
        self._reset_count = 0
        OutputData._instance_count += 1
        self._done_callback = done_callback

    @property
    def width(self):
        return self._description.width

    def __repr__(self):
        representation = (
            f'<{self.__class__.__name__}({self._description}) "{self._name}": ['
        )
        representation += ", ".join(
            repr(self._data[i]) if self._member_done[i] else "not done"
            for i in range(self._description.width)
        )
        representation += "]>"
        return representation

    @property
    def name(self):
        """The name of the published output."""
        return self._name

    def done(self, member=None):
        """Ensemble completion status for this output."""
        if member is not None:
            return self._member_done[member]
        else:
            return all(self._member_done)

    @typing.overload
    def data(self) -> typing.List[ResultTypeVar]:
        ...

    @typing.overload
    def data(self, member: int) -> ResultTypeVar:
        ...

    def data(self, member: int = None):
        """Access the raw data for localized output for the ensemble or the specified member."""
        if not self.done(member):
            raise exceptions.ApiError("Attempt to read before data has been published.")
        if member is not None:
            result = self._data[member]
        else:
            result = self._data
        # OutputData does not know its value type and cannot do a proper error check.
        # if result is None or None in result:
        #     raise exceptions.ApiError('Data marked "done" but contains null value.')
        return result

    def set(self, value: ResultTypeVar, member: int):
        """Set local data and mark as completed.

        Used internally to maintain the data store.
        """
        # Currently, we are not able to dynamically update the set of published outputs.
        # Allow `None` to indicate a null result.
        if self._member_done[member]:
            raise ApiError(
                f"{self}[{self.name}][{member}] is already set to {self._data[member]}!"
            )
        if value is None:
            self._data[member] = None
        elif self._description.dtype is not None and issubclass(
            self._description.dtype, gmxapi.abc.NDArray
        ):
            # Normalize sequence types to a gmxapi-specific type to aid in annotation and dispatching.
            # (Ultimately, we should use a gromacs-native data interface with universal shaped-type support,
            # but to simplify expression of new operations (e.g. function_wrapper) early gmxapi releases
            # use type annotations in familiar Python native types, where originally deemed feasible.)
            self._data[member] = datamodel.ndarray(value)
        else:
            if callable(self._description.dtype):
                try:
                    value = self._description.dtype(value)
                except TypeError as e:
                    logger.exception(
                        f"While trying to set {self} ({self._description}) to {value}: ",
                        exc_info=e,
                    )
                    raise e
            self._data[member] = value
        self._member_done[member] = True
        if all(self._member_done) and callable(self._done_callback):
            self._done_callback(self)

    def reset(self):
        """Reinitialize the data store.

        Note:
            This is a workaround until the execution model is more developed.

        Todo:
            Remove this method when all operation handles provide factories to
            reinstantiate on new Contexts and/or with new inputs.

        """
        self._member_done = [False] * self._description.width
        self._data = [None] * self._description.width
        self._reset_count += 1


class DataSourceCollection(collections.OrderedDict):
    """Store and describe input data handles for an operation.

    When created from InputCollectionDescription.bind(), the DataSourceCollection
    has had default values inserted.

    Note: DataSourceCollection is the input resource collection type for NodeBuilder.
    To use with the NodeBuilder interface, first create the collection, then
    iteratively add its resources.

    TODO: We should probably normalize the collection aspect of input. Is this a
          single input? Are all inputs added as collections? Should this be split
          to a standard iterable? Does a resource factory always produce a mapping?
    """

    def __init__(self, **kwargs):
        """Initialize from key/value pairs of named data sources.

        Data sources may be any of the basic gmxapi data types (or sequences thereof),
        or gmxapi Futures of those types.
        """
        super(DataSourceCollection, self).__init__()
        for name, value in kwargs.items():
            self[name] = value

    def __setitem__(self, key: str, item: SourceTypeVar) -> None:
        if not isinstance(key, str):
            raise exceptions.TypeError("Data must be named with str type.")

        # TODO(#3139): Encapsulate handling of provided data sources as a detail of the relationship
        #  between source and target Contexts.
        if isinstance(item, FutureABC):
            value = item
        else:
            # We check Mapping first, because a Mapping is a Sequence.
            if isinstance(item, collections.abc.Mapping):
                # Allow mappings of (possible) Futures to be automatically handled well.
                value = {}
                for name, obj in item.items():
                    value[str(name)] = obj
            elif isinstance(item, collections.abc.Sequence) and not isinstance(
                item, (str, bytes)
            ):
                try:
                    # Consider relaxing the requirements of ndarray, especially for the purposes
                    # of deducing ensemble inputs.
                    value = datamodel.ndarray(item)
                except (exceptions.ValueError, exceptions.TypeError) as e:
                    raise exceptions.TypeError(
                        "Iterable could not be converted to NDArray: {}".format(item)
                    ) from e
            elif isinstance(item, valid_source_types):
                value = item
            else:
                assert not isinstance(item, valid_source_types)
                # TODO: Support arbitrary types.
                raise exceptions.ApiError("Cannot process data source {}".format(item))
        super().__setitem__(key, value)

    @classmethod
    def _get_resettable(cls, obj):
        """Recursively (depth-first) yield possible data sources."""
        has_reset = lambda source: any(
            callable(getattr(source, method, None)) for method in ("reset", "_reset")
        )
        if has_reset(obj):
            yield obj
        elif isinstance(obj, collections.abc.Iterable):
            if isinstance(obj, collections.abc.Mapping):
                yield from cls._get_resettable(obj.values())
            elif not isinstance(obj, (str, bytes)):
                for element in obj:
                    if element is obj:
                        # Watch out for single-element iterables that return themselves.
                        break
                    yield from cls._get_resettable(element)

    def reset(self):
        """Reset all sources in the collection."""
        for source in self._get_resettable(self.values()):
            if hasattr(source, "reset"):
                source.reset()
            elif hasattr(source, "_reset"):
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
        raise exceptions.ApiError(
            "Can not inspect type of provided function argument."
        ) from T
    except ValueError as V:
        raise exceptions.ApiError("Can not inspect provided function signature.") from V

    wrapped_function = function_wrapper()(function)

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
                raise exceptions.TypeError(
                    "Output descriptions are keyed by Python strings."
                )
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

    # TODO: consider updating signature to *parameters or **parameters.
    def __init__(
        self, parameters: typing.Iterable[typing.Tuple[str, inspect.Parameter]]
    ):
        """Create the input description for an operation node from a dictionary of names and types."""
        inputs = []
        for name, param in parameters:
            if not isinstance(name, str):
                raise exceptions.TypeError(
                    "Input descriptions are keyed by Python strings."
                )
            # Multidimensional inputs are explicitly NDArray
            dtype = param.annotation
            if issubclass(dtype, collections.abc.Iterable) and not issubclass(
                dtype, (str, bytes, collections.abc.Mapping)
            ):
                # TODO: we can relax this with some more input conditioning.
                if dtype != datamodel.NDArray:
                    raise exceptions.UsageError(
                        f"Cannot accept input type {param}. Sequence type inputs must use NDArray."
                    )
            assert issubclass(dtype, valid_result_types)
            if hasattr(param, "kind"):
                disallowed = any(
                    [
                        param.kind == param.POSITIONAL_ONLY,
                        param.kind == param.VAR_POSITIONAL,
                        param.kind == param.VAR_KEYWORD,
                    ]
                )
                if disallowed:
                    raise exceptions.ProtocolError(
                        "Cannot wrap function. Operations must have well-defined parameter names."
                    )
                kind = param.kind
            else:
                kind = inspect.Parameter.POSITIONAL_OR_KEYWORD
            if hasattr(param, "default"):
                default = param.default
            else:
                default = inspect.Parameter.empty
            inputs.append(
                inspect.Parameter(name, kind, default=default, annotation=dtype)
            )
        super().__init__([(input.name, input.annotation) for input in inputs])
        self.signature = inspect.Signature(inputs)

    def __repr__(self):
        representation = f"<{self.__class__.__name__}: "
        representation += ", ".join(f"{key}={repr(val)}" for key, val in self.items())
        representation += ">"
        return representation

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
            disallowed = any(
                [
                    param.kind == param.POSITIONAL_ONLY,
                    param.kind == param.VAR_POSITIONAL,
                    param.kind == param.VAR_KEYWORD,
                ]
            )
            if disallowed:
                raise exceptions.ProtocolError(
                    "Cannot wrap function. Operations must have well-defined parameter names."
                )
            if param.name == "input":
                raise exceptions.ProtocolError(
                    'Function signature includes the (reserved) "input" keyword argument.'
                )
        description = collections.OrderedDict()
        for param in signature.parameters.values():
            if param.name == "output":
                # Wrapped functions may accept the output parameter to publish results, but
                # that is not part of the Operation input signature.
                continue
            if param.annotation == param.empty:
                if param.default == param.empty or param.default is None:
                    raise exceptions.ProtocolError(
                        f"Could not infer parameter type for {param.name}"
                    )
                dtype = type(param.default)
                if isinstance(dtype, collections.abc.Iterable) and not isinstance(
                    dtype, (str, bytes, collections.abc.Mapping)
                ):
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
        initializing an operation. Its most useful functionality is to  allows a
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
        # Additional kwargs (other than 'input' and 'output') override kwargs from 'input'.
        input_kwargs = collections.OrderedDict()
        if "input" in kwargs:
            provided_input = kwargs.pop("input")
            if provided_input is not None:
                input_kwargs.update(provided_input)
        # `function` may accept an `output` keyword argument that should not be supplied to the factory.
        for key, value in kwargs.items():
            if key == "output":
                raise exceptions.UsageError(
                    "Invalid keyword argument: output (reserved)."
                )
            input_kwargs[key] = value
        try:
            bound_arguments = self.signature.bind_partial(*args, **input_kwargs)
        except TypeError as e:
            raise exceptions.UsageError(
                "Could not bind operation parameters to function signature."
            ) from e
        assert "output" not in bound_arguments.arguments
        bound_arguments.apply_defaults()
        assert "input" not in bound_arguments.arguments
        input_kwargs = collections.OrderedDict(
            [pair for pair in bound_arguments.arguments.items()]
        )
        if "output" in input_kwargs:
            input_kwargs.pop("output")
        # Note that we have not done any type checking yet. Type and shape checking are part of
        # the "input binding" protocol, which consists of `NodeBuilder.add_input()` (which *may* be
        # implemented in terms of a `DataSourceCollection.__set_item__()`) and `NodeBuilder.build()`,
        # where checks occur in the DataEdge constructor using a DataSourceCollection and a
        # SinkTerminal (created from an InputCollectionDescription).
        return DataSourceCollection(**input_kwargs)


class ProxyDataDescriptor(object):
    """Base class for data descriptors used in DataProxyBase subclasses.

    Subclasses should either not define __init__ or should call the base class
    __init__ explicitly: super().__init__(self, name, dtype)
    """

    def __init__(self, name: str, dtype: ResultTypeVar = None):
        # Since Python 3.6, we have the opportunity to get the attribute name during
        # class creation using the __set_name__ hook.
        # TODO(4116): Simplify data descriptor set-up.
        self._name = name
        # TODO: We should not allow dtype==None, but we currently have a weak data
        #  model that does not allow good support of structured Futures.
        if dtype is not None:
            assert isinstance(dtype, type)
            assert issubclass(dtype, valid_result_types)
        self._dtype = dtype

    def __set_name__(self, owner, name):
        self._owner_name = owner.__qualname__
        if name != self._name:
            raise ApiError(
                f"Attribute name {self._owner_name}.{name} does not match "
                f"{self.__class__.__name__}._name: {self._name}"
            )

    def __repr__(self):
        return (
            f"<{self.__class__.__name__} "
            f'"{self._owner_name}.{self._name}", {self._dtype}>'
        )


class DataProxyMeta(abc.ABCMeta):
    # Key word arguments consumed by __prepare__
    _prepare_keywords = ("descriptors",)

    @classmethod
    def __prepare__(mcs, name, bases, descriptors: collections.abc.Mapping = None):
        """Allow dynamic sub-classing.

        DataProxy class definitions are collections of data descriptors. This
        class method allows subclasses to give the descriptor names and type(s)
        in the class declaration as arguments instead of as class attribute
        assignments.

            class MyProxy(DataProxyBase, descriptors={name: MyDescriptor() for name in datanames}): pass

        Note:
            Recent Python versions allow this to be replaced via ``__init_subclass__`` hook.
            See :issue:`4116`
        """
        if descriptors is None:
            return {}
        elif isinstance(descriptors, tuple):
            namespace = collections.OrderedDict([(d._name, d) for d in descriptors])
            return namespace
        else:
            assert isinstance(descriptors, collections.abc.Mapping)
            return descriptors

    def __new__(cls, name, bases: typing.Iterable, namespace, **kwargs):
        for key in kwargs:
            if key not in DataProxyMeta._prepare_keywords:
                raise exceptions.ApiError(
                    "Unexpected class creation keyword: {}".format(key)
                )
        # See note about DataProxyBase._reserved.
        if "_reserved" not in namespace and not any(
            hasattr(base, "_reserved") for base in bases
        ):
            raise exceptions.ApiError(
                "We currently expect DataProxy classes to provide a list of reserved attribute names."
            )
        for key in namespace:
            # Here we can check conformance with naming and typing rules.
            assert isinstance(key, str)
            if key.startswith("__"):
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
                if key not in namespace["_reserved"] and not any(
                    key in getattr(base, "_reserved")
                    for base in bases
                    if hasattr(base, "_reserved")
                ):
                    raise exceptions.ApiError(
                        "Unexpected data proxy attribute {}: {}".format(
                            key, repr(descriptor)
                        )
                    )
            else:
                assert isinstance(descriptor, ProxyDataDescriptor)
                if not isinstance(descriptor._name, str) or descriptor._name == "":
                    descriptor._name = key
                else:
                    if descriptor._name != key:
                        raise exceptions.ApiError(
                            "Descriptor internal name {} does not match attribute name {}".format(
                                descriptor._name, key
                            )
                        )
        return super().__new__(cls, name, bases, namespace)

    # TODO: This keyword argument stripping is not necessary in more recent Python versions.
    # When Python minimum required version is increased, check if we can remove this.
    def __init__(cls, name, bases, namespace, **kwargs):
        for key in kwargs:
            if key not in DataProxyMeta._prepare_keywords:
                raise exceptions.ApiError(
                    "Unexpected class initialization keyword: {}".format(key)
                )
        super().__init__(name, bases, namespace)

    # TODO: See if we can use __dir__ in the metaclass to help hint class attributes for better tab completion.
    #  Ref: https://ipython.readthedocs.io/en/stable/config/integrating.html#tab-completion
    # def __dir__(self) -> Iterable[str]:
    #     return super().__dir__()


class DataProxyBase(collections.abc.Mapping, metaclass=DataProxyMeta):
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
    _reserved = ("ensemble_width", "items", "_reserved")

    # This class can be expanded to be the attachment point for a metaclass for
    # data proxies such as PublishingDataProxy or OutputDataProxy, which may be
    # defined very dynamically and concisely as a set of Descriptors and a type()
    # call.
    # If development in this direction does not materialize, then this base
    # class is not very useful and should be removed.
    def __init__(self, instance: "SourceResource", client_id: int = None):
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
        # Collection is fixed by the time of instance creation, so cache it.
        self.__keys = tuple([key for key, _ in self.items()])
        self.__length = len(self.__keys)

    def __repr__(self):
        representation = (
            f"<{self.__class__.__name__} for {self._resource_instance} "
            f"client {self._client_identifier})>"
        )
        return representation

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

    def __getitem__(self, k):
        if hasattr(self, k):
            return getattr(self, k)

    def __len__(self):
        return self.__length

    def __iter__(self):
        for key in self.__keys:
            yield key


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
        return resource_manager._data.get(self._name).data(member=client_id)

    def __set__(self, instance: DataProxyBase, value):
        resource_manager = instance._resource_instance
        # TODO: Fix API scope.
        # Either this class is a detail of the same implementation as ResourceManager,
        # or we need to enforce that instance._resource_instance provides _data (or equivalent)
        assert isinstance(resource_manager, ResourceManager)
        client_id = instance._client_identifier
        resource_manager.set_result(name=self._name, value=value, member=client_id)

    def __repr__(self):
        return "{}(name={}, dtype={})".format(
            self.__class__.__name__, self._name, self._dtype.__qualname__
        )


def define_publishing_data_proxy(output_description) -> typing.Type[DataProxyBase]:
    """Returns a class definition for a PublishingDataProxy for the provided output description."""
    # This dynamic type creation hides collaborations with things like make_datastore.
    # We should encapsulate these relationships in Context details, explicit collaborations
    # between specific operations and Contexts, and in groups of Operation definition helpers.

    descriptors = collections.OrderedDict(
        [(name, Publisher(name)) for name in output_description]
    )

    # TODO: Instead of a local subclass and metaclass logic, we could follow the model
    #  of `dataclasses`.
    # At the very least, it could be helpful to distinguish DataProxies with Publishing
    # behavior versus Output behavior.
    class PublishingDataProxy(DataProxyBase, descriptors=descriptors):
        """Handler for write access to the `output` of an operation.

        Acts as a sort of PublisherCollection.
        """

    return PublishingDataProxy


# get symbols we can use to annotate input and output types more specifically.
_OutputDataProxyType = typing.TypeVar("_OutputDataProxyType", bound=DataProxyBase)
_PublishingDataProxyType = typing.TypeVar(
    "_PublishingDataProxyType", bound=DataProxyBase
)
# Currently, the ClientID type is an integer, but this may change.
ClientID = typing.NewType("ClientID", int)


class _Resources(typing.Generic[_PublishingDataProxyType]):
    pass


# TODO: Why generic in publishingdataproxytype?
class SourceResource(typing.Generic[_OutputDataProxyType, _PublishingDataProxyType]):
    """Resource Manager for a data provider.

    Supports Future instances in a particular context.
    """

    # Note: ResourceManager members not yet included:
    # future(), _data, set_result.

    def __init__(self):
        self._futures = weakref.WeakValueDictionary()

    # This might not belong here. Maybe separate out for a OperationHandleManager?
    @abc.abstractmethod
    def data(self) -> _OutputDataProxyType:
        """Get the output data proxy."""
        # Warning: this should probably be renamed, but "output_data_proxy" is already
        # a member in at least one derived class.
        ...

    @abc.abstractmethod
    def is_done(self, name: str) -> bool:
        return False

    @abc.abstractmethod
    def get(self, name: str) -> "OutputData":
        ...

    @abc.abstractmethod
    def update_output(self):
        """Bring the _data member up to date and local."""
        pass

    @abc.abstractmethod
    def reset(self):
        """Recursively reinitialize resources.

        Set the resource manager to its initialized state.
        All outputs are marked not "done".
        All inputs supporting the interface have ``_reset()`` called on them.
        """

    @abc.abstractmethod
    def width(self) -> int:
        """Ensemble width of the managed resources."""
        ...

    def future(self, name: str, description: ResultDescription) -> "Future":
        """Get a Future handle for managed data.

        Resource managers owned by subclasses of gmx.operation.Context provide
        this method to get references to output data.

        In addition to the interface described by gmx.abc.Future, returned objects
        provide the interface described by gmx.operation.Future.
        """
        key = (name, description.dtype, description.width)
        if key not in self._futures or self._futures[key] is None:
            _future = Future(self, name, description)
            self._futures[key] = _future
        return self._futures[key]


class StaticSourceManager(
    SourceResource[_OutputDataProxyType, _PublishingDataProxyType]
):
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

    def __init__(
        self, *, name: str, proxied_data, width: int, function: typing.Callable
    ):
        super().__init__()
        if isinstance(proxied_data, Future):
            raise ApiError(
                "Either StaticSourceManager is being misused, or a dispatching error has occurred for "
                f"{name}: {repr(proxied_data)}. "
                "StaticSourceManager is for managing local static data, "
                "but provided data is non-local or non-static. Please report bug."
            )
        self._result = function(proxied_data)
        # When width > 1, we expect `function(proxied_data)` to have produced a non-string
        # Iterable[dtype] of length `width`.
        if width > 1:
            if isinstance(self._result, (str, bytes)):
                # In this case, do not implicitly broadcast
                raise exceptions.ValueError(
                    '"function" produced data incompatible with "width".'
                )
            else:
                if not isinstance(self._result, collections.abc.Iterable):
                    raise exceptions.DataShapeError(
                        'Expected iterable of size {} but "function" result is not iterable.'
                    )
            data = list(self._result)
            size = len(data)
            if len(data) != width:
                raise exceptions.DataShapeError(
                    'Expected iterable of size {} but "function" produced a {} of size {}'.format(
                        width, type(data), size
                    )
                )
            dtype = type(data[0])
        else:
            if width != 1:
                raise exceptions.ValueError("width must be an integer 1 or greater.")
            dtype = type(self._result)
            if issubclass(dtype, (list, tuple)):
                dtype = datamodel.NDArray
                data = [datamodel.ndarray(self._result)]
            elif isinstance(self._result, collections.abc.Iterable):
                if not isinstance(self._result, (str, bytes, dict)):
                    raise exceptions.ValueError(
                        'Expecting width 1 but "function" produced iterable type {}.'.format(
                            type(self._result)
                        )
                    )
                else:
                    dtype = str
                    data = [str(self._result)]
            else:
                data = [self._result]
        if isinstance(None, dtype):
            dtype = None
        description = ResultDescription(dtype=dtype, width=width)
        self._data = OutputData(name=name, description=description)
        for member in range(width):
            self._data.set(data[member], member=member)

        output_collection_description = OutputCollectionDescription(**{name: dtype})
        self.output_data_proxy = define_output_data_proxy(
            output_description=output_collection_description
        )

        self.__representation = (
            f'<{self.__class__.__name__}: "{name}", '
            f"{repr(proxied_data)}, width={width}, function="
            f"{repr(function)}>"
        )

    def __repr__(self):
        return self.__representation

    def is_done(self, name: str) -> bool:
        return True

    def get(self, name: str) -> "OutputData":
        if name != self._data.name:
            raise KeyError(f'{repr(self)} does not hold "{name}"')
        return self._data

    def data(self) -> _OutputDataProxyType:
        return self.output_data_proxy(self)

    def width(self) -> int:
        # TODO: It looks like the OutputData ResultDescription probably belongs
        #  in the public interface.
        return self._data._description.width

    def update_output(self):
        pass

    def reset(self):
        pass

    def future(self, name: str, description: ResultDescription) -> "Future":
        return super().future(name, description=description)


class ProxyResourceManager(
    SourceResource[_OutputDataProxyType, _PublishingDataProxyType]
):
    """Act as a resource manager for data managed by another resource manager.

    Allow data transformations on the proxied resource.

    Keyword Args:
        proxied_future: An object implementing the Future interface.
        width: Size of (one-dimensional) shaped data produced by function.
        function: Transformation to perform on the result of proxied_future.
        dtype: (optional) Apply a check on the output data type.

    The callable passed as *function* must accept a single argument, which will
    be an iterable when proxied_future represents an ensemble, or an object of
    type proxied_future.description.dtype otherwise.

    If *width* is 1, *function* *should* produce a value with the intended output type.
    If *width* is greater than 1, *function* *should* produce a sequence with *width* elements.

    Warning:
        The GROMACS team has not yet approved a structured typing model for the API,
        so we cannot strongly assert source and target data type and shape. Data
        shape transformations may be ambiguous. If the heuristics for run-time
        transformations do not seem correct, please report bugs or comment at
        :issue:`2993`
    """

    def __init__(
        self,
        proxied_future: "Future",
        width: int,
        function: typing.Callable,
        dtype: typing.Optional[typing.Type] = None,
        # name: typing.Optional[str] = None
    ):
        super().__init__()
        self._proxied_result_available = False
        self._proxied_future = proxied_future
        self._width = width
        # if name is None:
        #     self.name = self._proxied_future.name
        # else:
        #     if not isinstance(name, str):
        #         raise GmxapiTypeError('*name* must be a string, or None.')
        #     self.name = name
        self.name = self._proxied_future.name
        self._result = None
        assert callable(function)
        self.function = function

        # Optional type hint for run time validation.
        self._dtype = dtype

        # Inferred data type of the transformed / proxied result.
        self._output_dtype = None

    def __repr__(self):
        status = "done" if self._done else "pending"
        dtype = self._output_dtype if self._output_dtype else self._dtype
        type_annotation = f" ({dtype.__name__})" if isinstance(dtype, type) else ""
        representation = (
            f"<{self.__class__.__name__} {self._width}{type_annotation} {status}>"
        )
        return representation

    def width(self) -> int:
        return self._width

    def reset(self):
        self._proxied_result_available = False
        self._proxied_future._reset()
        self._result = None
        self._output_dtype = None

    def is_done(self, name: str) -> bool:
        return self._proxied_result_available

    def _broadcast(
        self,
        result: ResultTypeVar,
        *,
        dtype: typing.Union[None, typing.Type[ResultTypeVar]],
    ):
        """Broadcast a scalar result to an ensemble."""
        if isinstance(dtype, type):
            if not isinstance(result, dtype):
                raise ApiError(f"{result} not compatible with {dtype}.")
        if dtype is None:
            dtype = type(result)
        data = OutputData(
            name=self.name,
            description=ResultDescription(dtype=dtype, width=self._width),
        )
        for member in range(self._width):
            data.set(result, member)
        return data

    def _map(self, result: collections.abc.Sequence, *, dtype):
        """Map elements of a result to the ensemble members."""
        if isinstance(result, str) or not isinstance(result, collections.abc.Sequence):
            raise ApiError(f"Expected array-like data. Got {result}.")
        data = OutputData(
            name=self.name,
            description=ResultDescription(dtype=dtype, width=self._width),
        )
        for member in range(self._width):
            data.set(result[member], member)
        return data

    def get(self, name: str):
        if name != self.name:
            raise exceptions.ValueError("Request for unknown data.")
        if not self.is_done(name):
            raise exceptions.ProtocolError("Data not ready.")
        result = self._result

        if self._width == 1:
            return self._broadcast(result, dtype=self._output_dtype)
        else:
            assert isinstance(result, collections.abc.Sequence)
            assert not isinstance(result, str)
            return self._map(result, dtype=self._output_dtype)

    def update_output(self):
        if not self._proxied_result_available:
            source = self._proxied_future.result()
            result = self.function(source)
            if self._dtype is None:
                if self._width == 1:
                    self._output_dtype = type(result)
                else:
                    assert self._width > 1
                    if isinstance(result, str) or not isinstance(
                        result, collections.abc.Sequence
                    ):
                        raise ApiError(
                            f"{self} must map a sequence of results to ensemble members. Got {result}."
                        )
                    elif len(result) != self._width:
                        raise DataShapeError(
                            f"Expected sequence of width {self._width}. Got {result}."
                        )
                    self._output_dtype = type(result[0])
            else:
                assert self._dtype is not None
                if isinstance(result, self._dtype):
                    self._output_dtype = type(result)
                else:
                    # We can try to implicitly broadcast as necessary, but let's try to leave that
                    # responsibility with the caller.
                    raise ApiError(f"Expected {self._dtype} but got {result}.")

            self._result = result
            self._proxied_result_available = True

    def data(self) -> _OutputDataProxyType:
        raise exceptions.ApiError(
            "ProxyResourceManager cannot yet manage a full OutputDataProxy."
        )

    def future(self, name: str, description: ResultDescription):
        return super().future(name, description=description)


class AbstractOperation(
    AbstractOperationReference, typing.Generic[_OutputDataProxyType]
):
    """Client interface to an operation instance (graph node).

    Note that this is a generic abstract class. Subclasses should provide a
    class subscript to help static type checkers.
    """

    @abc.abstractmethod
    def run(self):
        """Assert execution of an operation.

        After calling run(), the operation results are guaranteed to be available
        in the local context.
        """

    @property
    @abc.abstractmethod
    def output(self) -> _OutputDataProxyType:
        """Get a proxy collection to the output of the operation.

        Developer note: The 'output' property exists to isolate the namespace of
        output data from other operation handle attributes and we should consider
        whether it is actually necessary or helpful. To facilitate its possible
        future removal, do not enrich its interface beyond that of a collection
        of OutputDescriptor attributes.
        """
        ...


class OperationRegistryKey(typing.NamedTuple):
    """Helper class to define the key type for OperationRegistry."""

    namespace: str
    name: str

    def __hash__(self):
        return hash((self.namespace, self.name))

    def __eq__(self, other):
        """Test equivalence rather than identity.

        Note: Use `is` to test identity.
        """
        return other.namespace == self.namespace and other.name == self.name

    def __str__(self):
        return ".".join([self.namespace, self.name])

    def __repr__(self):
        return "{}(namespace={}, name={})".format(
            self.__qualname__, self.namespace, self.name
        )


def _make_registry_key(*args) -> OperationRegistryKey:
    """Normalize input to an OperationRegistryKey.

    Used to implement OperationRegistry.__getitem__(), which catches and converts
    the various exceptions this helper function can produce.
    """
    if len(args) > 1:
        return OperationRegistryKey(*args)
    else:
        if len(args) != 1:
            raise exceptions.UsageError(
                "Empty index value passed to OperationRegistry instance[]."
            )
        item = args[0]
    if isinstance(item, OperationRegistryKey):
        return item
    if isinstance(item, str):
        namespace, name = item.rsplit(sep=".", maxsplit=1)
        return OperationRegistryKey(namespace=namespace, name=name)
    # Item could be a class object or an instance object...
    if hasattr(item, "namespace") and hasattr(item, "name"):
        if callable(item.namespace):
            namespace = item.namespace()
        else:
            namespace = item.namespace
        if callable(item.name):
            name = item.name()
        else:
            name = item.name
        return OperationRegistryKey(namespace=namespace, name=name)
    raise exceptions.ValueError("Not a usable OperationRegistryKey: {}".format(item))


class OperationRegistry(collections.UserDict):
    """Helper class to map identifiers to Operation implementation instances.

    This is an implementation detail of gmxapi.operation and should not be used from
    outside of the package until a stable interface can be specified.
    """

    @typing.overload
    def __getitem__(self, item: OperationRegistryKey):
        ...

    @typing.overload
    def __getitem__(self, item: str):
        ...

    def __getitem__(self, *args):
        """Fetch the requested operation registrant.

        The getter is overloaded to support look-ups in multiple forms.

        The key can be given in the following forms.
        * As a period-delimited string of the form "namespace.operation".
        * As an OperationRegistryKey object.
        * As a sequence of positional arguments accepted by OperationRegistryKey.
        * As an object conforming to the OperationRegistryKey interface.
        """
        try:
            item = _make_registry_key(*args)
        except exceptions.Error as e:
            raise exceptions.TypeError(
                "Could not interpret key as a OperationRegistryKey."
            ) from e
        return self.data[item]


# Module level data store for locating operation implementations at run time.
# TODO: This may make sense as instance data of a root Context instance, but we
#  don't have a stable interface for communicating between Contexts yet.
# Alternatively, it could be held as class data in a RegisteredOperation class,
# but note that the "register" member function would behave less like the abc.ABCMeta
# support for "virtual subclassing" and more like the generic function machinery
# of e.g. functools.singledispatch.
_operation_registry = OperationRegistry()


def _register_operation(cls: typing.Type[OperationImplementation]):
    assert isinstance(cls, type)
    assert issubclass(cls, OperationImplementation)
    operation = _make_registry_key(cls)
    if operation in _operation_registry:
        full_name = str(operation)
        raise exceptions.ProtocolError(
            "Attempting to redefine operation {}.".format(full_name)
        )
    _operation_registry[operation] = cls


# TODO: replace with a generic function that we dispatch on so the type checker can infer a return type.
def _get_operation_director(operation, context: gmx.abc.Context):
    """

    :param operation:
    :param context:
    :return: gmxapi.abc.OperationDirector
    """
    registrant = _operation_registry[operation]
    director = registrant.director(context=context)
    return director


class InputDescription(abc.ABC):
    """Node input description for gmxapi.operation module.

    Provides the support needed to understand operation inputs in gmxapi.operation
    module Contexts.

    .. todo:: Implementation base class with heritable behavior and/or helpers to
              compose this functionality from more normative description of
              operation inputs. This will probably become a facet of the ResourceDirector
              when specialized for gmxapi.operation.Context.
    """

    @abc.abstractmethod
    def signature(self) -> InputCollectionDescription:
        """Mapping of named inputs and input type.

        Used to determine valid inputs before an Operation node is created.

        Collaborations:
            Related to the operation resource factory for this context.

        ..  todo::
            Better unification of this protocol, InputCollectionDescription, and
            resource factory.
            Note, also, that the *bind* method of the returned InputCollectionDescription
            serves as the resource factory for input to the node builder.
        """
        ...

    @abc.abstractmethod
    def make_uid(self, input: "DataEdge") -> str:
        """The unique identity of an operation node tags the output with respect to the input.

        Combines information on the Operation details and the input to uniquely
        identify the Operation node.

        Arguments:
            input : A (collection of) data source(s) that can provide Fingerprints.

        Used internally by the Context to manage ownership of data sources, to
        locate resources for nodes in work graphs, and to manage serialization,
        deserialization, and checkpointing of the work graph.

        The UID is a detail of the generic Operation that _should_ be independent
        of the Context details to allow the framework to manage when and where
        an operation is executed.

        TODO: We probably don't want to allow Operations to single-handedly determine their
         own uniqueness, but they probably should participate in the determination with the Context.

        TODO: Context implementations should be allowed to optimize handling of
         equivalent operations in different sessions or work graphs, but we do not
         yet guarantee that UIDs are globally unique!

        To be refined...
        """
        ...


class ConcreteInputDescription(InputDescription):
    """Simple composed InputDescription."""

    def __init__(
        self,
        input_signature: InputCollectionDescription,
        uid_helper: typing.Callable[["DataEdge"], str],
    ):
        self._input_signature_description = input_signature
        self._uid_helper = uid_helper

    def signature(self) -> InputCollectionDescription:
        return self._input_signature_description

    def make_uid(self, input: "DataEdge") -> str:
        return self._uid_helper(input)


class OperationMeta(abc.ABCMeta):
    """Metaclass to manage the definition of Operation implementation classes.

    Design Note:
        Note that this metaclass can be superseded by `__init_subclass__()`.
        See :issue:`4116`.

    """

    def __new__(meta, name, bases, class_dict):
        cls = super().__new__(meta, name, bases, class_dict)
        # Register subclasses, but not the base class.
        if (
            issubclass(cls, OperationImplementation)
            and cls is not OperationImplementation
        ):
            # TODO: Remove OperationDetailsBase layer and this extra check.
            # Note: we do not yet register the Operations built dynamically because we
            # don't have a clear definition of unique implementations yet. For instance,
            # a new OperationDetails class is defined for each call to gmx.join_arrays
            # TODO: Properly register and reuse Operations defined dynamically
            #  through function_wrapper (currently encompassed by OperationDetailsBase subclasses)
            if name != "OperationDetailsBase":
                if OperationDetailsBase not in bases:
                    _register_operation(cls)
        return cls


class OperationDetailsBase(
    OperationImplementation, InputDescription, metaclass=OperationMeta
):
    """Abstract base class for Operation details in this module's Python Context.

    Provides necessary interface for use with gmxapi.operation.ResourceManager.
    Separates the details of an Operation from those of the ResourceManager in
    a given Context.

    OperationDetails classes are almost stateless, serving mainly to compose implementation
    details. Instances (operation objects) provide the Context-dependent interfaces
    for a specific node in a work graph.

    OperationDetails subclasses are created dynamically by function_wrapper and
    make_operation.

    Developer note: when subclassing, note that the ResourceManager is responsible
    for managing Operation state. Do not add instance data members related to
    computation or output state.

    TODO: determine what is acceptable instance data and/or initialization information.
    Note that currently the subclass in function_wrapper has _no_ initialization input,
    but does not yet handle input-dependent output specification or node fingerprinting.
    It seems likely that instance initialization will require some characterization of
    supplied input, but nothing else. Even that much is not necessary if the instance
    is completely stateless, but that would require additional parameters to the member
    functions. However, an instance should be tied to a specific ResourceManager and
    Context, so weak references to these would be reasonable.
    """

    @abc.abstractmethod
    def output_description(self) -> OutputCollectionDescription:
        """Mapping of available outputs and types for an existing Operation node."""
        ...

    @abc.abstractmethod
    def publishing_data_proxy(
        self,
        *,
        instance: SourceResource[typing.Any, _PublishingDataProxyType],
        client_id,
    ) -> _PublishingDataProxyType:
        """Factory for Operation output publishing resources.

        Used internally when the operation is run with resources provided by instance."""
        ...

    @abc.abstractmethod
    def output_data_proxy(
        self, instance: SourceResource[_OutputDataProxyType, typing.Any]
    ) -> _OutputDataProxyType:
        """Get an object that can provide Futures for output data managed by instance."""
        ...

    @abc.abstractmethod
    def __call__(self, resources: _Resources):
        """Execute the operation with provided resources.

        Resources are prepared in an execution context with aid of resource_director()

        After the first call, output data has been published and is trivially
        available through the output_data_proxy()
        """
        ...

    @classmethod
    @abc.abstractmethod
    def resource_director(
        cls, *, input, output: _PublishingDataProxyType
    ) -> _Resources[_PublishingDataProxyType]:
        """a Director factory that helps build the Session Resources for the function.

        The Session launcher provides the director with all of the resources previously
        requested/negotiated/registered by the Operation. The director uses details of
        the operation to build the resources object required by the operation runner.

        For the Python Context, the protocol is for the Context to call the
        resource_director instance method, passing input and output containers.
        (See, for example, gmxapi.operation.PyFunctionRunnerResources)
        """
        ...

    # TODO: Don't run the director. Just return the correct callable.
    @classmethod
    def operation_director(
        cls, *args, context: "Context", label=None, **kwargs
    ) -> AbstractOperation:
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
        if not isinstance(context, Context):
            raise exceptions.UsageError("Context instance needed for dispatch.")
        # TODO: use Context characteristics rather than isinstance checks.
        if isinstance(context, ModuleContext):
            construct = OperationDirector(
                *args, operation_details=cls, context=context, label=label, **kwargs
            )
            return construct()
        elif isinstance(context, SubgraphContext):
            construct = OperationDirector(
                *args, operation_details=cls, context=context, label=label, **kwargs
            )
            return construct()
        else:
            raise exceptions.ApiError(
                "Cannot dispatch operation_director for context {}".format(context)
            )


# TODO(#3139): Implement observer pattern for edge->node data flow between Contexts or
#  ResourceManagers.
# (Lower level details would be expected to use a Futures interface, but we need ways to proxy
# between environments to set up fulfilment of the Futures.)
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


class Future(GenericFuture[ResultTypeVar]):
    """gmxapi data handle.

    Future is currently more than a Future right now. (should be corrected / clarified.)
    Future is also a facade to other features of the data provider.

    Future objects are the most atomic interface between Contexts. User scripts
    may hold Futures from which they extract data with result(). Operation output
    used as input for another Operation can be decomposed such that the Operation
    factory has only Future objects in its input.

    TODO(#3139): ``add_done_callback`` method allows consumers to bind as Observers.

    Currently abstraction is handled through SourceResource subclassing.

    Attributes:
        description (ResultDescription): Describes the result to be obtained from this Future.

    """

    description: ResultDescription
    name: str
    resource_manager: SourceResource

    def __init__(
        self,
        resource_manager: SourceResource,
        name: str,
        description: ResultDescription,
    ):
        self.name = name
        if not isinstance(description, ResultDescription):
            raise exceptions.ValueError("Need description of requested data.")
        self.description = description
        self.resource_manager = resource_manager
        # Note that we cannot confirm at this time that resource_manager can provide *name* because
        # SourceResource does not specify how implementations represent their output other than the
        # existence of a data() method. SourceResource.data() may lazily instantiate an OutputDataProxy
        # that only exposes output names as Future instances, which would cause an infinite recursion loop
        # if we tried to inspect here.

        # Deprecated. We should not "reset" futures, but reconstitute them, but we
        # need to move the data model to a subscription-based system so that we can
        # make Futures properly immutable and issue new ones across subgraph iterations.
        self._number_of_resets = 0

    def __repr__(self):
        return "<Future: name='{}', description={}>".format(self.name, self.description)

    def result(self) -> typing.Union[ResultTypeVar, typing.List[ResultTypeVar]]:
        """Fetch data to the caller's Context.

        Returns an object of the concrete type specified according to
        the operation that produces this Result.

        Ensemble data are returned as a list. Scalar results or results from single
        member ensembles are returned as scalars. If this behavior is confusing or
        problematic for you, please reopen https://gitlab.com/gromacs/gromacs/-/issues/3179
        or a new issue and join the discussion.
        """
        self.resource_manager.update_output()
        # Return ownership of concrete data
        handle: OutputData[ResultTypeVar] = self.resource_manager.get(self.name)
        assert isinstance(handle, OutputData)

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
        _width = self.description.width
        if _width == 1:
            result = handle.data(member=0)
        else:
            result = handle.data()
        # Normalize non-native sequences to a type that has general
        # serialization support and built-in dispatching.
        if isinstance(result, gmxapi.abc.NDArray):
            result = list(result)
        if (_width > 1) and (len(result) != _width):
            assert len(result) == 1
            # Implicitly broadcast
            result = [copy(result[0]) for _ in range(_width)]
        return result

    def _reset(self):
        """Mark the Future "not done" to allow reexecution.

        Invalidates cached results, resets "done" markers in data sources, and
        triggers _reset recursively.

        Note: this is a hack that is inconsistent with the plan of unique mappings
        of inputs to outputs, but allows a quick prototype for looping operations.
        """
        self._number_of_resets += 1
        self.resource_manager.reset()

    @property
    def dtype(self):
        return self.description.dtype

    def __getitem__(self, item):
        """Get a more limited view on a Future for subscriptable type.

        The result has the same ensemble dimensions as the original Future.

        For ensemble Futures, subscripting is applied to the object of each member, not to the
        ensemble array dimension. There is not currently a facility to hold Futures for specific
        non-local ensemble members. (Call *result()* to localize the ensemble data, then index
        locally.)

        Slicing is not supported.

        See Also:
            ArrayFuture
        """
        element_dtype = None

        if issubclass(self.description.dtype, collections.abc.Mapping):
            # TODO(#3130): Stronger typed gmxapi Mapping.
            element_dtype = None
        elif issubclass(self.description.dtype, NDArray):
            # TODO(#3130): Use proper data shaping instead of weakly typed containers if possible.
            # Try to extract dtype for fancier sequence-providers.
            if hasattr(self.description.dtype, "dtype"):
                if isinstance(self.description.dtype.dtype, type):
                    element_dtype = self.description.dtype.dtype
                elif callable(self.description.dtype.dtype):
                    element_dtype = self.description.dtype.dtype()

        description = ResultDescription(
            dtype=element_dtype, width=self.description.width
        )
        if description.width == 1:
            proxy = ProxyResourceManager(
                self,
                width=description.width,
                function=lambda value, key=item: value[key],
            )
        else:
            proxy = ProxyResourceManager(
                self,
                width=description.width,
                function=lambda value, key=item: [
                    subscriptable[key] for subscriptable in value
                ],
            )
        return proxy.future(self.name, description=description)


class OutputDataDescriptor(ProxyDataDescriptor):
    """Read-only data descriptor for proxied access to output data.

    Knows how to get a Future from the resource manager.
    """

    # TODO: Reconcile the internal implementation details with the visibility and
    #  usages of this class.

    def __get__(self, proxy: DataProxyBase, owner):
        if proxy is None:
            # Access through class attribute of owner class
            return self
        result_description = ResultDescription(
            dtype=self._dtype, width=proxy.ensemble_width
        )

        return proxy._resource_instance.future(
            name=self._name, description=result_description
        )


class MutableResourceDescriptor(ProxyDataDescriptor):
    """Accessor for rich binding interfaces.

    Allows operations to access resources beyond the scope of the current
    resource manager. Used by operations whose interactions are more complicated
    than standard typed data flow at the scope of the current Context.

    Instead of a Future interface, the returned object is a MutableResource with
    which a subscriber can collaborate with lower-level protocols.
    """

    def __get__(
        self, proxy: DataProxyBase, owner
    ) -> typing.Union[MutableResource, "MutableResourceDescriptor"]:
        if proxy is None:
            # Access through class attribute of owner class. We don't have a
            # specified use case for that, so allow inspection of the data
            # descriptor instance, itself.
            return self
        # TODO: implement.
        # The protocol for MD extension plugins requires that the simulation operation
        # subscribe to the plugin. Then the Context allows the plugin to access the
        # MdRunner interface as the simulation is launched.
        # The protocol for modify_input and for mdrun to consume the TPR payload
        # of read_tpr or modify_input should allow us to use the gmxapi 0.0.7
        # WorkSpec to configure and launch a simulation, which we can do by feeding
        # forward and building a fused operation at the mdrun node. The information
        # fed forward can just be references to the inputs and parameters of the
        # earlier operations, with annotations so that we know the intended behavior.


def define_output_data_proxy(
    output_description: OutputCollectionDescription,
) -> typing.Type[DataProxyBase]:
    descriptors = {
        name: OutputDataDescriptor(name, description)
        for name, description in output_description.items()
    }

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


# Encapsulate the description of the input data flow.
PyFuncInput = collections.namedtuple("Input", ("args", "kwargs", "dependencies"))


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

        :param input_collection_description: Available inputs for Operation
        :return: Fully formed description of the Sink terminal for a data edge to be created.

        Collaborations: Execution Context implementation.
        """
        self.ensemble_width = 1
        self.inputs = input_collection_description

    def __repr__(self):
        return f"<SinkTerminal: width={self.ensemble_width}, {self.inputs}>"

    def update_width(self, width: int):
        if not isinstance(width, int):
            try:
                width = int(width)
            except TypeError:
                raise exceptions.TypeError("Need an integer width > 0.")
        if width < 1:
            raise exceptions.ValueError(
                "Nonsensical ensemble width: {}".format(int(width))
            )
        elif width > 1:
            if self.ensemble_width != 1:
                if width != self.ensemble_width:
                    raise exceptions.ValueError(
                        "Cannot change ensemble width {} to width {}.".format(
                            self.ensemble_width, width
                        )
                    )
            self.ensemble_width = width
        else:
            assert width == 1
            # source data of width 1 is compatible with all consumer widths.

    def update(self, data_source_collection: DataSourceCollection):
        """Update the SinkTerminal with the proposed data provider.

        The primary purpose for this part of the protocol is to allow the sink to widen itself
        (assert that it is now part of an "ensemble" workflow) when implied by the shape of the input.
        """
        for name, sink_dtype in self.inputs.items():
            if name not in data_source_collection:
                # If/when we accept data from multiple sources, we'll need some additional sanity checking.
                if not hasattr(self.inputs.signature.parameters[name], "default"):
                    raise exceptions.UsageError(
                        "No data or default for {}".format(name)
                    )
            else:
                # With a single data source, we need data to be in the source or have a default
                assert name in data_source_collection
                assert issubclass(sink_dtype, valid_result_types)
                source = data_source_collection[name]
                logger.debug("Updating Sink for source {}: {}.".format(name, source))
                if isinstance(source, sink_dtype):
                    # Note: We do not currently have a way for sinks to advertise dimensionality.
                    # Check containers for nested ensemble input.
                    if issubclass(sink_dtype, collections.abc.Mapping):
                        assert isinstance(source, collections.abc.Mapping)
                        for value in source.values():
                            if hasattr(value, "description"):
                                source_description = typing.cast(
                                    ResultDescription, value.description
                                )
                                self.update_width(source_description.width)
                    elif issubclass(
                        sink_dtype, collections.abc.Sequence
                    ) and not issubclass(sink_dtype, str):
                        # This will need updating when Sinks can advertise greater than
                        # 1-dimensional input.
                        if (
                            all(
                                not isinstance(element, str)
                                and isinstance(
                                    element,
                                    (collections.abc.Sequence, Future, FutureABC),
                                )
                                for element in source
                            )
                            and len(source) > 0
                        ):
                            # If all elements of the list are lists or Futures, ensemble_width may be
                            # len(source) or may need to derive from one of the widths. Let's try to
                            # defer resolution of that ambiguity for the moment.
                            logger.warning(
                                "Ambiguous data shape: sending fully 2-dimensional "
                                f"input {source} to Array consumer."
                            )
                            # Try our best by assuming the outer dimension is an ensemble dimension.
                            width = len(source)
                            self.update_width(width)
                        else:
                            # If any element of the list is either a list or a Future with width
                            # greater than 1, call update for the wide element.
                            width = 0
                            for element in source:
                                if isinstance(
                                    element,
                                    (collections.abc.Sequence, FutureABC, Future),
                                ) and not isinstance(element, str):
                                    if hasattr(element, "description"):
                                        _width = typing.cast(
                                            ResultDescription, element.description
                                        ).width
                                    else:
                                        _width = len(element)
                                else:
                                    _width = 1
                                if width > 1:
                                    if _width > 1 and _width != width:
                                        raise DataShapeError(
                                            f"Inconsistent shape in data source {source}."
                                        )
                                else:
                                    width = max(width, _width)
                            if width:
                                self.update_width(width)
                    else:
                        logger.debug("Source matches sink. No update necessary.")
                else:
                    if isinstance(source, collections.abc.Iterable) and not isinstance(
                        source, (str, bytes, collections.abc.Mapping)
                    ):
                        assert isinstance(source, datamodel.NDArray)
                        if sink_dtype != datamodel.NDArray:
                            # Source is NDArray, but sink is not. Implicitly scatter.
                            self.update_width(len(source))
                        continue
                    if hasattr(source, "description"):
                        source_description = typing.cast(
                            ResultDescription, source.description
                        )
                        source_dtype = source_description.dtype
                        assert isinstance(sink_dtype, type)
                        # TODO: Handle typing of Future slices when we have a better data model.
                        if source_dtype is not None:
                            assert isinstance(source_dtype, type)
                            if not issubclass(source_dtype, sink_dtype):
                                raise exceptions.TypeError(
                                    "Expected {} but got {}.".format(
                                        sink_dtype, source_dtype
                                    )
                                )
                        source_width = source.description.width
                        self.update_width(source_width)


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

    This part of the data flow object model will be unnecessary with a more complete Future
    pattern and data model. Future slicing and broadcasting can be reconciled at binding time,
    and containers like DataSourceCollection and SinkTerminal can be reduced to Protocols or
    removed from the interface entirely.
    """

    # TODO: We can separate out these helpers to simplify DataEdge.
    class ConstantResolver(object):
        def __init__(self, value):
            self.value = value

        def __call__(self, member=None):
            return self.value

    class MappingResolver:
        """Adapter for resolving Mappings that combine Futures and static values.

        Wraps a dictionary like object to support data source resolution at run time.
        When called, produces a dictionary of localized values, calling ``value.result()``
        for values that are gmxapi Futures. For Futures with width greater than 1,
        results are sliced for the given ensemble member, if specified with the
        *member* key word argument.

        This class provides an explicit functor type to improve readability of the
        "adapters" dispatching in DataEdge.

        Note that gmxapi Mapping values do not support any sort of rigorous typing,
        pending #2993.
        """

        def __init__(self, mapping: collections.abc.Mapping):
            self.source = mapping

        def _resolve(self, member):
            for key, value in self.source.items():
                try:
                    if hasattr(value, "result") and callable(value.result):
                        # This logic allows ensemble Futures to be mixed with non-ensemble values.
                        # In the long run, the entire data structure (or View) should encapsulate the
                        # behavior for providing a consistent dimensionality, making appropriate
                        # updates internally during construction.
                        _result = value.result()
                        try:
                            width = value.description.width
                        except (TypeError, AttributeError):
                            width = 1
                        if width > 1 and member is not None:
                            # Make sure that list indices are valid.
                            try:
                                _result[member]
                            except Exception as e:
                                raise e
                            yield key, value.result()[member]
                        else:
                            yield key, _result
                    else:
                        yield key, value
                except Exception as e:
                    message = (
                        f'Trouble getting "{key}" for {self}. {value}.result() raised: '
                    )
                    warnings.warn(message + str(e))
                    logger.exception(message, exc_info=e)
                    # At some point, we should try to continue to allow debugging without immediate MPI_ABORT.
                    # We should also (https://gitlab.com/gromacs/gromacs/-/issues/3394) capture errors in a
                    # richer result representation so that the consumer (an operation-specific resource
                    # directory, runner director, or runner, or the client of a Future protocol) can take
                    # responsibility for interpreting the error and deciding what to do about it.
                    # Currently, however, we expect exceptions to propagate.
                    # yield key, None
                    raise e

        def __call__(self, member=None) -> dict:
            return dict(self._resolve(member))

    def __init__(
        self, source_collection: DataSourceCollection, sink_terminal: SinkTerminal
    ):
        """Reconcile the provided data source(s) with the advertised inputs.

        DataSourceCollection is constructed with consideration of the named inputs, but not the
        type or shape. (Currently) In NodeBuilder.build(), the sink is given a chance to check
        inputs and adjust its ensemble state by inspection of the DataSourceCollection via
        `input_sink.update(self.sources)`. The DataEdge is instantiated afterwards, so we should
        be able to use the sink ensemble width for some short-circuit logic and/or checking of
        invariants/assumptions.
        """

        # Adapters are callables that transform a source and node ID to local data.
        # Every key in the sink has an adapter.
        self.adapters = {}
        self.source_collection = source_collection
        self.sink_terminal = sink_terminal
        # For each input on the sink, dispatch data dependency resolution with help from source
        # and sink type (and shape).
        for name in sink_terminal.inputs:
            logger.debug(f"Resolving {name} for {repr(sink_terminal)}.")
            if name not in source_collection:
                if hasattr(sink_terminal.inputs[name], "default"):
                    self.adapters[name] = self.ConstantResolver(
                        sink_terminal.inputs[name]
                    )
                    logger.debug(
                        f"Using default value for {name} for {repr(sink_terminal)}."
                    )
                else:
                    raise exceptions.ValueError(
                        f"{repr(sink_terminal)}:{name} has no default, and no {name} in {source_collection}."
                    )
            else:
                source = source_collection[name]
                # Design note: The following dispatching is a bit awkward. Suggestion for future
                # implementation: Dispatch first on sink type, then source type.
                sink = sink_terminal.inputs[name]

                # Scalar -> (Scalar, Array[Scalar])
                if isinstance(source, (str, bool, int, float)):
                    logger.debug(
                        f"Input {name}:({sink.__name__}) provided by a local constant of "
                        f"type {type(source)}."
                    )
                    if issubclass(sink, (str, bool, int, float)):
                        self.adapters[name] = self.ConstantResolver(source)
                    else:
                        if issubclass(sink, (datamodel.NDArray, tuple, list)):
                            # Allow simple constants to implicitly convert to single-element arrays.
                            # The initial version of NDArray is a primitive local-only container.
                            self.adapters[name] = self.ConstantResolver(
                                datamodel.ndarray([source])
                            )
                        else:
                            raise UsageError(
                                f"Input {name} of type {sink} cannot accept data "
                                f"source {source}"
                            )
                # (Mapping, Array[Mapping], Future[Mapping]) -> dict
                elif issubclass(sink, dict):
                    if isinstance(source, collections.abc.Mapping):
                        if all(
                            isinstance(value, valid_source_types)
                            for value in source.values()
                        ):
                            logger.debug(
                                f"Input {name}:({sink.__name__}) provided by a local constant of "
                                f"type {type(source)}."
                            )
                            self.adapters[name] = self.ConstantResolver(source)
                        else:
                            self.adapters[name] = self.MappingResolver(source)
                    elif (
                        isinstance(source, collections.abc.Sequence)
                        and len(source) > 0
                        and isinstance(source[0], (collections.abc.Mapping, FutureABC))
                    ):
                        if len(source) == 1:
                            # Handle direct mapping or broadcast.
                            source = source[0]
                            if isinstance(source, collections.abc.Mapping):
                                self.adapters[name] = self.MappingResolver(source)
                            else:
                                # If the Future is part of an ensemble, result() will return a list.
                                # Otherwise, it will return a single object.
                                ensemble_width = source.description.width
                                # TODO: subscribe to futures so results can be pushed.
                                if ensemble_width == 1:
                                    self.adapters[
                                        name
                                    ] = lambda member, source=source: source.result()
                                else:
                                    self.adapters[
                                        name
                                    ] = lambda member, source=source: source.result()[
                                        member
                                    ]
                        else:
                            # The sink should already be updated before we get here.
                            assert sink_terminal.ensemble_width == len(source)
                            # Handle parallel mapping.
                            resolvers: typing.List[typing.Callable] = [object] * len(
                                source
                            )
                            for i, element in enumerate(source):
                                if isinstance(element, collections.abc.Mapping):
                                    resolvers[i] = self.MappingResolver(element)
                                else:
                                    source_width = element.description.width
                                    if source_width > 1:
                                        raise DataShapeError(
                                            f"Cannot nest {element} with width {source_width} "
                                            f"in parallel data edge with width {sink_terminal.ensemble_width}."
                                        )
                                    resolvers[
                                        i
                                    ] = lambda member, _source=element: _source.result()
                            self.adapters[name] = lambda member, _resolvers=tuple(
                                resolvers
                            ): _resolvers[member]()
                    elif hasattr(source, "result"):
                        # TODO: Simplify data model and type checking.
                        # Handle data futures...
                        logger.debug(
                            f"Input {name}:({sink.__name__}) provided by a Future "
                            f"{source.description}"
                        )
                        # If the Future is part of an ensemble, result() will return a list.
                        # Otherwise, it will return a single object.
                        ensemble_width = source.description.width
                        # TODO: subscribe to futures so results can be pushed.
                        if ensemble_width == 1:
                            self.adapters[
                                name
                            ] = lambda member, source=source: source.result()
                        else:
                            self.adapters[
                                name
                            ] = lambda member, source=source: source.result()[member]
                    else:
                        raise ApiError(
                            f"Input type {sink} cannot accept source {repr(source)}"
                        )
                # Array -> Array
                elif isinstance(source, datamodel.NDArray):
                    if issubclass(sink, datamodel.NDArray):
                        # TODO(#3136): shape checking
                        # Implicit broadcast may not be what is intended
                        if isinstance(source.dtype, type) and issubclass(
                            source.dtype, (list, tuple, datamodel.NDArray)
                        ):
                            # Try to treat array of arrays as an ensemble.
                            if len(source) == sink_terminal.ensemble_width:
                                self.adapters[
                                    name
                                ] = lambda member, _source=source.to_list(): _source[
                                    member
                                ]
                                continue
                            # Otherwise, hope that the sink is supposed to be a list of lists...
                        self.adapters[name] = self.ConstantResolver(source)
                    else:
                        if source.shape[0] != sink_terminal.ensemble_width:
                            raise exceptions.ValueError(
                                f"Implicit broadcast could not match array source {source} to "
                                f"ensemble sink {sink_terminal}"
                            )
                        else:
                            self.adapters[name] = lambda member, source=source: source[
                                member
                            ]
                # Future[Any] -> Any
                elif hasattr(source, "result"):
                    # TODO: Simplify data model and type checking.
                    # Handle data futures...
                    logger.debug(
                        f"Input {name}:({sink.__name__}) provided by a Future "
                        f"{source.description}"
                    )
                    # If the Future is part of an ensemble, result() will return a list.
                    # Otherwise, it will return a single object.
                    ensemble_width = source.description.width
                    # TODO: subscribe to futures so results can be pushed.
                    if ensemble_width == 1:
                        self.adapters[
                            name
                        ] = lambda member, source=source: source.result()
                    else:
                        self.adapters[
                            name
                        ] = lambda member, source=source: source.result()[member]
                else:
                    raise ApiError(
                        f"Input type {sink} cannot accept source {repr(source)}"
                    )
        for name in sink_terminal.inputs:
            assert name in self.adapters

    def __repr__(self):
        return "<DataEdge: source_collection={}, sink_terminal={}>".format(
            self.source_collection, self.sink_terminal
        )

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


class AbstractRunnerDirector(abc.ABC):
    allow_duplicate: bool = False

    @abc.abstractmethod
    def __call__(self, resources) -> typing.Callable[[], None]:
        ...


class PublishingManager(
    typing.Generic[_PublishingDataProxyType], contextlib.AbstractContextManager
):
    """Manages the output phase of task fulfilment for ResourceManager.

    * Generates the context managers for publishing operation results.
    * Provide a synchronization point in the otherwise flexible update_output() protocol.
    * Helps to ensure that each operation does not execute more than the expected number of times to
      provide output to consumers.

    Create one instance of PublishingManager for each ResourceManager.

    When entered as a Python context manager, the PublishingManager is "activated"
    and can provide context managers for PublishingDataProxies for each ensemble
    member exactly once (unless *reset*).

    When the main PublishingManager context is exited, a sequence of observers are contacted.
    The observers may be responsible for actually publishing the results of the operation,
    so the data state is not inspected until after the observers have been notified.

    After the visitors have been called, the main context manager checks whether
    all expected outputs have been published. If output is missing, and if no exceptions
    are already propagating, an exception is raised.
    """

    def __init__(
        self,
        resource_manager: "ResourceManager",
        publisher_factory: typing.Type[_PublishingDataProxyType],
    ):
        # Design Note: we could make this class into a Data Descriptor in order to
        # avoid storing a reference to the ResourceManager.
        self._manager = weakref.ref(resource_manager)

        assert callable(publisher_factory)
        self._factory: typing.Callable[
            [...], _PublishingDataProxyType
        ] = publisher_factory
        self._active = False
        self.publishing_data_proxies: typing.Dict[str, _PublishingDataProxyType] = {}
        self._observers: typing.List[typing.Callable[["PublishingManager"], None]] = []

    @contextmanager
    def publishing_context(self, ensemble_member=0) -> _PublishingDataProxyType:
        """Get a data proxy to which results may be published."""
        if not self._active:
            raise ApiError(
                "Before calling *publishing_context*, activate the PublishingManager by entering "
                "it as a Python context manager."
            )
        manager = self._manager()
        if manager is None:
            # Probably due to inappropriate use outside a ResourceManager instance.
            raise ApiError("ResourceManager does not exist.")

        if manager.done(ensemble_member):
            raise exceptions.ProtocolError(
                "Attempting to publish {}[{}] more than once.".format(
                    manager.operation_id, ensemble_member
                )
            )

        if ensemble_member in self.publishing_data_proxies:
            raise exceptions.ProtocolError(
                f"Publishing resources for {manager}[{ensemble_member}] already active."
            )

        try:
            resource = self._factory(instance=manager, client_id=ensemble_member)
        except Exception as e:
            logger.debug("Publishing context could not be created due to {}".format(e))
            raise e
        self.publishing_data_proxies[ensemble_member] = resource
        yield resource

    def reset(self):
        """Allow the manager to be used again.

        Generally, a ResourceManager executes an operation once to publish results,
        and PublishingManager helps enforce this by preventing context managers
        from acquiring a PublishingDataProxy for the same output more than once.
        In the gmxapi 0.1 while_loop implementation, ResourceManagers persist for
        all iterations of a subgraph and get reset for each iteration.

        We want to preserve the subscribed callbacks across iterations, so the
        PublishingManager gets reset, too.
        """
        if self._active:
            raise gmxapi.exceptions.ProtocolError(
                "Attempt to reset PublishingManager while still in use."
            )
        self.publishing_data_proxies.clear()
        self._observers = []

    def register_observer(self, cb: typing.Callable[["PublishingManager"], None]):
        self._observers.append(cb)

    def __enter__(self):
        assert len(self.publishing_data_proxies) == 0
        self._active = True
        return self

    def __exit__(self, __exc_type, __exc_value, __traceback):
        manager = self._manager()
        if __exc_value:
            logger.info(
                f"Exception occurred while updating {manager}. {self} exiting early."
            )
        else:
            for member in range(manager.ensemble_width):
                if member not in self.publishing_data_proxies:
                    logger.warning(
                        "PublishingManager has not been used for all ensemble members."
                    )
            assert all(
                isinstance(resource, DataProxyBase)
                for resource in self.publishing_data_proxies.values()
            )
            for observer in self._observers:
                try:
                    observer(self)
                except Exception as e:
                    logger.exception(
                        f"Exception while notifying PublishingManager observer {observer}",
                        exc_info=e,
                    )
            if not manager.done():
                logger.warning(
                    f"Publishing resources are being released, but {manager} does not appear done."
                )
        self._active = False
        # Allow __exc_value to be reraised:
        return False


class ResourceManager(SourceResource[_OutputDataProxyType, _PublishingDataProxyType]):
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

    TODO(#3139: Data should be pushed, not pulled.
    Early implementations executed operation code and extracted results directly.
    While we need to be able to "wait for" results and alert the data provider that
    we are ready for input, we want to defer execution management and data flow to
    the framework. Some forward motion is prototyped in a *add_done_callback* method.
    """

    @property
    def ensemble_width(self):
        return self._input_edge.sink_terminal.ensemble_width

    def __init__(
        self,
        *,
        source: DataEdge,
        operation_id,
        output_description: OutputCollectionDescription,
        output_data_proxy: typing.Type[_OutputDataProxyType],
        publishing_data_proxy: typing.Type[_PublishingDataProxyType],
        resource_factory,
        runner_director: AbstractRunnerDirector,
        output_context: "Context",
    ):
        """Initialize a resource manager for the inputs and outputs of an operation."""
        super().__init__()
        # Note: This implementation assumes there is one ResourceManager instance per data source,
        # so we only stash the inputs and dependency information for a single set of resources.
        # TODO: validate input_fingerprint as its interface becomes clear.
        self._input_edge = source

        self._base_operation_id = operation_id

        if isinstance(output_context, Context):
            self._output_context = output_context
        else:
            message = (
                "Provide an instance of gmxapi.operation.Context for output_context"
            )
            raise exceptions.UsageError(message)
        assert self._output_context is not None

        self._output_data_proxy = output_data_proxy
        assert self._output_data_proxy is not None
        assert callable(self._output_data_proxy)

        self._output_description = output_description
        assert self._output_description is not None

        self._runner_director = runner_director
        assert self._runner_director is not None
        # Discover from the RunnerDirector whether the task (whose resources we
        # are managing) should be executed on each rank where its results are
        # needed or whether we will have to synchronize its outputs across ranks
        # after executing in a single place.
        self._allow_duplicate = runner_director.allow_duplicate

        self._resource_factory = resource_factory
        assert self._resource_factory is not None

        self._data = DataStore(
            output_description=self._output_description,
            ensemble_width=self.ensemble_width,
            done_callback=self._receive_completion_signal,
        )

        # TODO: reimplement as a data descriptor
        #  so that PublishingManager (and PublishingDataProxy) do not need a bound circular reference.
        self.__publishing_resources = PublishingManager(
            resource_manager=self, publisher_factory=publishing_data_proxy
        )

        """Functions to call when the resource is brought up to date."""
        self._callbacks: typing.List[typing.Callable[[ResourceManager], None]] = []

        self.__operation_entrance_counter = 0
        self.__reset_counter = 0

    @property
    def operation_id(self):
        # Note that duplicate work on multiple ranks, if allowed, will have the same ID.
        return f"{self._base_operation_id}_i{self.__reset_counter}"

    def __repr__(self):
        representation = f"<{self.__class__.__qualname__} {self.operation_id}: width={self.width()}, "
        representation += ", ".join(
            f"{name}: {data}" for name, data in self._data.items()
        )
        representation += ">"
        return representation

    def _receive_completion_signal(self, obj: "DataStore"):
        # Finalize. Confirm that the runner calls succeeded in finalizing the work graph node.
        #
        # Note that ResourceManager subclasses may have
        # differently structured update_output that have an ensemble-wide
        # publishing_resources instead of using one for each ensemble member in sequence.
        # However, we need a general solution for ResourceManagers that are relying on
        # subscriptions to the ResourceManagers of other ranks for their results.
        # The providers need to run their callbacks, and the subscribers need to run a
        # corresponding set of callbacks in the same sequence at the same point in the data flow
        # resolution (then publish their results).
        #
        # Clearly, if we will continue to subclass ResourceManager and provide different algorithms for
        # bringing output up to date, we need to break up this function into the separate required
        # protocols versus extensible behaviors.

        if comm_size > 1:
            logger.debug(f"Synchronizing {self.operation_id} update_output().")

        self._run_callbacks()

        # This is probably overly conservative, but we can revisit for optimization
        # after the gmxapi 0.3.0 release.
        if comm_size > 1:
            comm.barrier()
            logger.debug("Passed barrier.")

    def width(self) -> int:
        return self.ensemble_width

    def reset(self):
        self.__operation_entrance_counter = 0
        # TODO: Update details of subscription relationships.
        # TODO: Update details of subscription relationships.
        self.__publishing_resources.reset()
        self._data.reset()
        self._input_edge.reset()
        self._callbacks = []
        assert self.__operation_entrance_counter == 0
        self.__reset_counter += 1

    def done(self, member=None):
        return all(data.done(member) for data in self._data.values())

    def is_done(self, name, member=None):
        return self._data[name].done(member=member)

    def set_result(self, name, value, member: int):
        if not isinstance(value, (str, bytes)):
            try:
                for item in value:
                    # In this specification, it is antithetical to publish Futures.
                    if hasattr(item, "result"):
                        raise exceptions.ApiError(
                            "Operation produced Future instead of real output."
                        )
            except TypeError:
                # Ignore when `item` is not iterable.
                pass
        # Note that data is not considered "done" until all members have been set.
        self._data[name].set(value=value, member=member)

    def _make_done_callback(
        self, func: typing.Callable[[OutputData], None], name: str, member: int
    ) -> typing.Callable[["ResourceManager"], None]:
        """Wrap a callback prototype for ResourceManager.add_done_callback.

        This shim helps us to independently evolve the ResourceManager interface,
        the Future interface, and the core callback functions that we need.

        We need to reconcile our notions of
        Future slicing and adopt a more conventional interface where `add_done_callback()`
        can be serviced by Future.
        """
        if member >= self.ensemble_width:
            raise gmxapi.exceptions.ValueError(
                f"No member {member} in {self} (width {self.ensemble_width})."
            )
        if name not in self._data:
            raise gmxapi.exceptions.ValueError(f'No output "{name}" in {self}.')
        return self._wrap_callback(func, name, self.operation_id)

    @staticmethod
    @functools.lru_cache()
    def _wrap_callback(
        func: typing.Callable[[OutputData], None], name: str, operation_id
    ) -> typing.Callable[["ResourceManager"], None]:
        """Wrap and cache a callback prototype for ResourceManager._make_done_callback.

        This is a separate static member to avoid potential circular reference problems with
        caching member functions, but we do want the ResourceManager identity to be
        part of the cache key.

        We would like ResourceManager.add_done_callback to be able to check whether
        specific callbacks have already been added, but we can't do that if we generate
        a new functools.partial object with each call to ResourceManager._make_done_callback.
        So this avoids making a new partially bound function if we have already done so.
        """

        def _filter(obj: ResourceManager):
            return func(obj._data[name])

        logger.debug(
            f"Created wrapper function {_filter} for {operation_id}[{name}] callback {func}."
        )

        return _filter

    def add_done_callback(self, func: typing.Callable[["ResourceManager"], None]):
        """Attaches the callable to the collection of managed resources.

        *func* will be called, with the ResourceManager as its only argument, when the
        operation finishes running.

        If the operation has already completed, *func* will be called immediately.

        Note that this signature is not yet normative.

        This is a step in the direction of callbacks for individual Futures.
        This implementation would have to accept *name* and *member* parameters
        in order to register callbacks to specific Futures. For the moment, we
        leave the responsibility with the caller to pick out the particular data
        it needs. The callback receives a reference to the ResourceManager itself,
        rather than a Future object. In the mean time, internally we use a helper
        functions to provide a bridge.

        """

        if func in self._callbacks:
            raise UsageError("Function {func} already registered in {self} callbacks.")
        self._callbacks.append(func)
        if self.done():
            try:
                func(self)
            except Exception as e:
                logger.exception(
                    "Exception during callbacks for {self}[{self.name}][{member}].",
                    exc_info=e,
                )

    def _run_callbacks(self):
        """Run all the scheduled callbacks when this operation is complete."""
        for cb in self._callbacks:
            try:
                cb(self)
            except Exception as e:
                logger.exception(
                    f"Exception in ResourceManager callback {cb}.", exc_info=e
                )

    def get(self, name: str) -> OutputData:
        """Get managed data by name.

        Raises exceptions.ProtocolError if requested data is not local yet.
        Raises exceptions.ValueError if data is requested for an unknown name.
        """
        if name not in self._data:
            raise exceptions.ValueError(
                f"Request for unknown data: {self.operation_id}:{name}"
            )
        if not self.is_done(name):
            raise exceptions.ProtocolError(f"{self.operation_id}:{name} not ready.")
        data = self._data[name]
        assert isinstance(data, OutputData)
        return data

    # TODO: Normalize. This is the no-argument, no-return callable member of an
    #  operation handle that dispatches to another Context (the operation implementation).
    # TODO: Allow update of a single ensemble member. As written, this always updates all ensemble members. That can
    #  be the default behavior, but we don't want to require non-local updates in all cases.
    def update_output(self):
        """Bring the output of the bound operation up to date.

        Execute the bound operation once if and only if it has not
        yet been run in the lifetime of this resource manager.

        Used internally to implement Futures for the local operation
        associated with this resource manager.

        Raises:
            exceptions.ApiError if operation runner fails to publish output.

        TODO: More comprehensive error handling for operations that fail to execute.

        When overriding this function, subclasses MUST use the PublishingManager
        context manager to ensure correct data flow.
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
            # Note: This check could also be encapsulated in a run_once decorator that
            # could even set a data descriptor to change behavior.
            self.__operation_entrance_counter += 1
            if self.__operation_entrance_counter > 1:
                raise exceptions.ProtocolError(
                    "Bug detected: resource manager tried to execute operation twice."
                )
            with self.publishing_resources() as publisher:
                # Note! This is a detail of the ResourceManager in a SerialContext
                for ensemble_member in range(self.ensemble_width):
                    # TODO: rewrite the following expression as a call to a resource factory.
                    # TODO: Consider whether the resource_factory behavior should be normalized
                    #  to always use `with` blocks to indicate the lifetime of a resource handle.
                    #  That implies that an operation handle can expire, but the operation handle
                    #  could be "yield"ed
                    #  from within the `with` block to keep the resource scope alive until the resulting
                    #  generator is exhausted. Not sure what that looks like or what the use case would be.
                    with self.local_input(ensemble_member) as input:
                        # Note: Resources are marked "done" by the publishing system
                        # before the following context manager finishes exiting.
                        with publisher.publishing_context(
                            ensemble_member=ensemble_member
                        ) as output:
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
                            #
                            # This resource factory signature might need to be inverted or broken up
                            # into a builder for consistency. I.e.
                            # option 1: Make the input and output resources with separate factories and add_resource on
                            # the runner builder.
                            # option 2: Pass resource_builder to input_director and then output_director.
                            error_message = (
                                "Got {} while executing {} for operation {}."
                            )

                            # output is a PublishingDataProxy
                            try:
                                # TODO: Provide runtime session resources with
                                #  consideration for what the operation actually needs
                                #  (see builder proposal above).
                                resources = self._resource_factory(
                                    input=input, output=output
                                )
                            except exceptions.TypeError as e:
                                message = error_message.format(
                                    e, self._resource_factory, self.operation_id
                                )
                                raise exceptions.ApiError(message) from e

                            runner = self._runner_director(resources)
                            try:
                                runner()
                            except Exception as e:
                                message = error_message.format(
                                    e, runner, self.operation_id
                                )
                                raise exceptions.ApiError(message) from e
                            # end `with publishing_context`
                        # end `with local_input`
                    # end `for ensemble_member`
                # end `with publisher`
            # PublishingManager will have made its observer callbacks now.

            # For gmxapi 0.3, we require that all implementations of ResourceManager.update_output()
            # also use PublishingManager instances to ensure synchronization of output publishing.
            # PublishingManager services callbacks when all runners have completed but before
            # asserting that the publishing_resources are "done". This is necessary as a bug fix,
            # to allow results to come through callbacks rather than from the runners, but it is not
            # necessarily a long-term feature.
            #
            # We could handle this better in a few different ways. We would either
            # need to communicate with a workflow manager to negotiate publishing to subscribers
            # (or just to receive a list of synchronizations or additional calls to make), or we would
            # need to allow asynchronous/concurrent finalization of tasks in `update_output`.
            # This is well beyond the scope of gmxapi 0.3.

        # We want to be sure that the ResourceManager will _receive_completion_signal
        # and call the registered callbacks immediately after this function ends.
        # Subclasses are advised to do the same until we have a more robust behavioral composition.
        if not self.done():
            message = (
                "update_output implementation failed to update all outputs for {}."
            )
            message = message.format(self.operation_id)
            raise exceptions.ApiError(message)

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
            message = (
                "Requested Future of type {} is not compatible with available type {}."
            )
            message = message.format(requested_dtype, available_dtype)
            raise exceptions.ApiError(message)
        return super().future(name, description)

    def data(self) -> _OutputDataProxyType:
        """Get an adapter to the output resources to access results."""
        # In the gmxapi 0.1 spec, the OutputDataProxy never binds to a specific
        # ensemble member (using the DataProxyBase *client_id* parameter).
        # This will probably not be changed before flattening the OutputDataProxy
        # into a set of Data Descriptors directly in the OperationHandle namespace
        # (remove the 'output' attribute: issue #3174), and should be coupled to
        # consideration of issue #3179.
        return self._output_data_proxy(self)

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
        # Localize data.
        # Implementation may call `result()`, triggering recursive `update_output()` for
        # dependencies.
        kwargs = self._input_edge.sink(node=member)
        assert "input" not in kwargs

        # Check that we have real data (not a Future or operation handle or something).
        for key, value in kwargs.items():
            assert not hasattr(value, "result")
            assert not hasattr(value, "run")
            value_list = []
            if isinstance(value, list):
                value_list = value
            if isinstance(value, datamodel.NDArray):
                value_list = value._values
            if isinstance(value, collections.abc.Mapping):
                value_list = value.values()
            assert not isinstance(value_list, Future)
            assert not hasattr(value_list, "result")
            assert not hasattr(value_list, "run")
            # Check one level into container objects to confirm that they aren't hiding
            # non-localized data.
            for item in value_list:
                if hasattr(item, "result"):
                    raise ApiError(
                        f"{repr(item)} in {repr(self._input_edge)} (key: {key}) "
                        f"unexpectedly looks like a Future."
                    )

        input_pack = InputPack(kwargs=kwargs)

        # Prepare input data structure
        # Note: we use 'yield' instead of 'return' for the protocol expected by
        # the @contextmanager decorator
        yield input_pack

    def publishing_resources(self) -> PublishingManager:
        """Get a context manager for resolving the data dependencies of this node.

        Activate the publishing resource factory by entering it as a Python context manager.

        Use the *publishing_context* method on the activated factory to get a publishing data proxy for a
        given ensemble_member.

        When these context managers are exited, the output publishing protocol is automatically satisfied.

        'output' type resources can be published exactly once, and only while the
        publishing context is active.

        Write access to publishing resources can be granted exactly once during the
        resource manager lifetime and conveys exclusive access. However, this state can be *reset*.
        """
        publishing_resource_factory: PublishingManager = self.__publishing_resources
        # The ensemble_width is certainly fixed by now, and this seems like a
        # reasonable place to establish subscription relationships between ranks.
        # self.publishing_resources() gets the (single-use) context manager
        # creation function (self.__publishing_context) that will take a member
        # id as an argument.
        if comm_size > 1 and not self._allow_duplicate:
            # owning_rank is the process where the work is actually executed. This is
            # hard-coded in the initial implementation. It is subject to optimization
            # and increased flexibility in future version.
            owning_rank = 0

            # Confirm that this code is executed everywhere we think it is.
            logger.debug(f"Synchronizing {self} on all ranks with a Comm.barrier.")
            # This is probably overly conservative, but we can revisit for optimization
            # after the gmxapi 0.3.0 release.
            comm.barrier()

            for member in range(self.ensemble_width):
                for name, output in self._data.items():
                    assert isinstance(output, OutputData)
                    assert output.name == name
                    # For each member of OutputData, either add a publishing callback for the Future or a
                    # subscriber callback to the publishing resource factory.
                    # Set up proxy subscribers for each member of each OutputData
                    if rank_number == owning_rank:
                        # Set up the callbacks for the owning rank to push updates.
                        for subscriber in [
                            rank for rank in range(comm_size) if rank != owning_rank
                        ]:
                            outputdata_cb = functools.partial(
                                _subscription_provider_callback,
                                ensemble_member=member,
                                subscriber=subscriber,
                            )
                            # Convert from a callback that takes OutputData to one that takes ResourceManager
                            cb = self._make_done_callback(
                                func=outputdata_cb, name=name, member=member
                            )
                            self.add_done_callback(func=cb)
                    else:
                        # We are on the subscribing rank.

                        publisher_cb: typing.Callable[
                            [_PublishingDataProxyType], None
                        ] = functools.partial(
                            _subscriber_callback,
                            name=name,
                            ensemble_member=member,
                            owner=owning_rank,
                        )
                        # Convert from a callback that takes OutputData to one that takes PublishingManager
                        cb = self._make_publisher_callback(
                            func=publisher_cb, name=name, member=member
                        )
                        publishing_resource_factory.register_observer(cb)

        return publishing_resource_factory

    def _make_publisher_callback(
        self,
        func: typing.Callable[[_PublishingDataProxyType], None],
        name: str,
        member: int,
    ):
        logger.debug(
            f"Wrapping publisher callback for {self} {name}[{member}]: {func}."
        )

        def callback(obj: PublishingManager):
            output = obj.publishing_data_proxies[member]
            return func(output)

        return callback


def _subscription_provider_callback(
    obj: OutputData, ensemble_member: int, subscriber: int
):
    """Support a subscription to the given element of the named output.

    We arrange for this to be called by the ResourceManager

    This function is a prototype: it holds the necessary behavior, but must be
    rewrapped into the correct signature.
    Partially bind the extra arguments to achieve the signature needed.
    """
    # This callback will be called just after ResourceManager.update_output()
    # when triggered by the successful publishing of all outputs on the
    # owning rank.
    result = obj.data(member=ensemble_member)
    try:
        result = result.to_list()
    except AttributeError:
        pass
    try:
        result = result.tolist()
    except AttributeError:
        pass
    if not isinstance(result, (str, bool, int, float, dict, list)):
        raise ApiError()
    logger.debug(f"Sending {obj.name}[{ensemble_member}]={result} to {subscriber}")
    comm.send(result, dest=subscriber)
    logger.debug(f"Sent {obj.name}[{ensemble_member}] to {subscriber}")


def _subscriber_callback(
    output: _PublishingDataProxyType, name: str, ensemble_member: int, owner: int
):
    """Receive a callback from the PublishingManager to be called before it is finalized.

    This function prototype assigns output in the way a runner would, but we use a hook
    in the PublishingManager to achieve correct synchronization with the ResourceManager
    "done" callbacks.
    """
    # The callback will be called when the publishing resources are being finalized
    # at the end of the ResourceManager.update_output() method.
    # Note that PublishingDataProxies hide the ensemble dimension.
    # We account for that here in wrapping code.
    logger.debug(f"Waiting for {name}[{ensemble_member}] from {owner}")
    result = comm.recv(source=owner)
    logger.debug(f"Received {name}[{ensemble_member}]={result}. Publishing locally.")
    setattr(output, name, result)


class PyFunctionRunnerResources(collections.UserDict):
    """Runtime resources for Python functions.

    Produced by a ResourceDirector for a particular Operation.
    """

    def output(self):
        if "output" in self:
            return self["output"]
        else:
            return None

    def input(self):
        return {key: value for key, value in self.items() if key != "output"}


class PyFunctionRunner(abc.ABC):
    def __init__(
        self,
        *,
        function: typing.Callable,
        output_description: OutputCollectionDescription,
    ):
        """Instantiate the wrapper for *function* and its parameters.

        *function* is a callable that does not accept PyFunctionRunnerResources directly.
        `wrapped_function_runner()` instantiates an appropriate subclass of
        PyFunctionRunner that will translate PyFunctionRunnerResources for *function* at
        run time.

        *function* is a callable that either returns a single value (to be captured as
        *data*) or which assigns output values to the attributes of an *output* keyword
        argument.
        """
        assert callable(function)
        self.function = function
        self.output_description = output_description

    @abc.abstractmethod
    def __call__(self, resources: PyFunctionRunnerResources):
        # Subclass should implement something like
        #     self.function(output=resources.output(), **resources.input())
        # and should not call super().__call__
        raise NotImplementedError


class CapturedOutputRunner(PyFunctionRunner):
    """Function runner that captures return value as output.data"""

    def __call__(self, resources: PyFunctionRunnerResources):
        resources["output"].data = self.function(**resources.input())


class OutputParameterRunner(PyFunctionRunner):
    """Function runner that uses output parameter to let function publish output."""

    def __call__(self, resources: PyFunctionRunnerResources):
        self.function(**resources)


class NoOpRunner(PyFunctionRunner):
    """Function runner that does not execute on the current rank.

    Alternative PyFunctionRunner for ranks where function should not be executed.

    Output will be synchronized separately through the ResourceManager facilities.
    """

    def __call__(self, resources: PyFunctionRunnerResources):
        logger.debug(f"null runner received {resources}.")


def wrapped_function_runner(
    function,
    output_description: OutputCollectionDescription = None,
    owning_rank: int = 0,
    allow_duplicate=False,
) -> PyFunctionRunner:
    """Get an adapter for a function to be wrapped.

    If the function does not accept a publishing data proxy as an `output`
    key word argument, the returned object has a `capture_output` attribute that
    must be re-assigned by the calling code before calling the runner. `capture_output`
    must be assigned to be a callable that will receive the output of the wrapped
    function.

    Returns:
        Callable with a signature `__call__(*args, **kwargs)` and no return value

    Collaborations:
        * OperationDetails.resource_director assigns the `capture_output` member of the
          returned object.
        * The publishing resources factory must be ready to make a compatible decision
          with respect to ensemble data flow for the current MPI rank.
    """
    # WARNING:  The hard-coded value for the owning rank is duplicated in the
    # ResourceManager.publishing_resources function. To deduplicate would require
    # either a global setting or more run-time dispatching logic. Currently,
    # this factory is called while dynamically defining the OperationDetails class for
    # function_wrapper. Another option would be to use additional information from
    # function_wrapper when calling this factory.

    assert callable(function)
    signature = inspect.signature(function)
    logger.debug(f"Creating runner for signature {signature}.")

    # Implementation note: this function dispatches an implementation with the
    # logic below. A better factoring would be a "chain of responsibility" in
    # which the concrete Runners would be tried in sequence and determine internally
    # whether to create a runner, raise an error, or defer.

    # Determine output details for proper dispatching.
    # First check for signature with output parameter.
    if "output" in signature.parameters:
        if not isinstance(output_description, OutputCollectionDescription):
            if not isinstance(output_description, collections.abc.Mapping):
                raise exceptions.UsageError(
                    "Function passes output through call argument, but output is not described."
                )
            runner = OutputParameterRunner(
                function=function,
                output_description=OutputCollectionDescription(**output_description),
            )
        else:
            runner = OutputParameterRunner(
                function=function, output_description=output_description
            )
    # Next try output_description parameter or function return annotation.
    else:
        if isinstance(output_description, OutputCollectionDescription):
            return_type = output_description["data"].gmxapi_datatype
        elif output_description is not None:
            # output_description should be None for inferred output or
            # a singular mapping of the key 'data' to a gmxapi type.
            if not isinstance(output_description, collections.abc.Mapping) or set(
                output_description.keys()
            ) != {"data"}:
                raise exceptions.ApiError(
                    "invalid output description for wrapped function: {}".format(
                        output_description
                    )
                )
            if signature.return_annotation != signature.empty:
                if signature.return_annotation != output_description["data"]:
                    raise exceptions.ApiError(
                        "Wrapped function with return-value-capture provided with non-matching output "
                        "description."
                    )
            return_type = output_description["data"]
        else:
            # Use return type inferred from function signature.
            return_type = signature.return_annotation
        if return_type == signature.empty or return_type is None:
            raise exceptions.ApiError(
                "No return annotation or output_description for {}".format(function)
            )
        runner = CapturedOutputRunner(
            function=function,
            output_description=OutputCollectionDescription(data=return_type),
        )
    if rank_number != owning_rank and not allow_duplicate:
        # Only do actual execution on the root rank. This rank will subscribe to
        # results from the root rank.
        # Replace the runner acquired above.
        runner = NoOpRunner(
            function=runner.function, output_description=runner.output_description
        )

    return runner


# TODO: Refactor in terms of reference to a node in a Context.
#  ResourceManager is an implementation detail of how the Context
#  manages a node.
class OperationHandle(AbstractOperation[_OutputDataProxyType]):
    """Generic Operation handle for dynamically defined operations.

    Define a gmxapi Operation for the functionality being wrapped by the enclosing code.

    An Operation type definition encapsulates description of allowed inputs
    of an Operation. An Operation instance represents a node in a work graph
    with uniquely fingerprinted inputs and well-defined output. The implementation
    of the operation is a collaboration with the resource managers resolving
    data flow for output Futures, which may depend on the execution context.
    """

    def __init__(
        self, resource_manager: SourceResource[_OutputDataProxyType, typing.Any]
    ):
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
        #  module-level Context, convert to a weakref or look up based on stored node uid.
        self.__resource_manager = resource_manager

    @property
    def node_uid(self):
        # The unique identifier for the operation node allows the Context implementation
        # to manage the state of the handle. Reproducibility of node_uid is TBD, but
        # it must be unique in a Context where it references a different operation node.
        uid = getattr(self.__resource_manager, "operation_id", None)
        return uid

    def __repr__(self):
        representation = f"<{self.__class__.__name__} ({self.__resource_manager})>"
        return representation

    @property
    def output(self) -> _OutputDataProxyType:
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
        # Note: `run()` is a synonym for `resolve` or `update` or whatever we choose
        #  to generically describe the request to bring a node up-to-date: i.e. the
        #  non-returning callable on the object produced by a director.
        self.__resource_manager.update_output()


class OperationPlaceholder(AbstractOperation):
    """Placeholder for Operation handle during subgraph definition."""

    def __init__(self, subgraph_resource_manager):
        ...

    def run(self):
        raise exceptions.UsageError(
            "This placeholder operation handle is not in an executable context."
        )

    @property
    def output(self):
        """Allow subgraph components to be connected without instantiating actual operations."""
        if not isinstance(current_context(), SubgraphContext):
            raise exceptions.UsageError("Invalid access to subgraph internals.")


_HandleType = typing.TypeVar("_HandleType", bound=gmx.abc.OperationReference)


class NodeBuilder(gmx.abc.NodeBuilder):
    """Add an operation node to be managed by a Context.

    The NodeBuilder interface implies minimal internal logic, and instead
    specifies the set of information that must or may be provided to construct
    a node.
    """

    def __init__(
        self, context: "Context", operation, label: typing.Optional[str] = None
    ):
        """Initialize the base NodeBuilder for gmxapi.operation module Nodes.

        TODO:
            Convert the *operation* argument to be the registrant in the Operation registry.
            Requires confirmation of conforming behavior for dynamically defined operations.

        """
        self.context = context
        self.label = label
        try:
            key = _make_registry_key(operation)
        except Exception as e:
            error = "Could not create an operation registry key from {}".format(
                operation
            )
            raise exceptions.ValueError(error) from e
        else:
            # TODO: sensibly handle dynamically defined operations.
            if key not in _operation_registry and not issubclass(
                operation, OperationDetailsBase
            ):
                error = "{} must be initialized with a registered operation. Got {}."
                raise exceptions.ValueError(
                    error.format(__class__.__qualname__, operation)
                )
        self.sources = DataSourceCollection()
        self._input_description = None
        self._resource_factory = None
        self._runner_director = None
        self._handle = None
        self._output_factory = None

        self._resource_manager = ResourceManager

    def set_input_description(self, input_description: InputDescription):
        self._input_description = input_description

    def set_resource_factory(self, factory):
        self._resource_factory = factory

    def set_runner_director(self, factory):
        self._runner_director = factory

    def set_handle(self, factory):
        self._handle = factory

    def set_output_factory(self, factory: "OutputFactory"):
        self._output_factory = factory

    def set_resource_manager(self, implementation: typing.Type[ResourceManager]):
        """Allow overriding the default ResourceManager implementation.

        This is a workaround until we figure out what parts of the ResourceManager
        could be composed, registered separately with the Context, or be customized
        through other dispatching. One likely part of the solution is for clients
        of the NodeBuilder to assert requirements of the Context.
        """
        assert issubclass(implementation, ResourceManager)
        self._resource_manager = implementation

    # TODO: Let the Context use the handle factory to produce a dynamically typed handle,
    #  and figure out the right way to annotate the return value.
    def build(self):
        """Create node and return a handle of the appropriate type."""

        # Check for the ability to instantiate operations.
        missing_details = list()
        for builder_resource in [
            "input_description",
            "resource_factory",
            "runner_director",
            "handle",
            "output_factory",
        ]:
            detail = "_" + builder_resource
            if getattr(self, detail, None) is None:
                missing_details.append(builder_resource)
        if len(missing_details) > 0:
            raise exceptions.UsageError(
                "Missing details needed for operation node: {}".format(
                    ", ".join(missing_details)
                )
            )

        assert hasattr(self._input_description, "signature")
        input_sink = SinkTerminal(self._input_description.signature())
        input_sink.update(self.sources)
        logger.debug("SinkTerminal configured: {}".format(SinkTerminal))
        edge = DataEdge(self.sources, input_sink)
        logger.debug(
            "Created data edge {} with Sink {}".format(edge, edge.sink_terminal)
        )
        # TODO: Fingerprinting: Each operation instance has unique output based on the unique input.
        #            input_data_fingerprint = edge.fingerprint()

        # Set up output proxy.
        assert hasattr(self._input_description, "make_uid")
        uid = self._input_description.make_uid(edge)
        # TODO: ResourceManager should fetch the relevant factories from the Context
        #  instead of getting an OperationDetails instance.
        output_data_proxy = self._output_factory.output_proxy()
        output_description = self._output_factory.output_description()
        publishing_data_proxy = self._output_factory.publishing_data_proxy()
        manager = self._resource_manager(
            output_context=self.context,
            source=edge,
            operation_id=uid,
            output_data_proxy=output_data_proxy,
            output_description=output_description,
            publishing_data_proxy=publishing_data_proxy,
            resource_factory=self._resource_factory,
            runner_director=self._runner_director,
        )
        self.context.work_graph[uid] = manager
        # TODO: Replace with a call to Node.handle()
        handle = self._handle(self.context.work_graph[uid])
        return handle

    def add_input(self, name, source):
        # TODO: We can move some input checking here as the data model matures.
        self.sources[name] = source


class InputPack(object):
    """Input container for data sources provided to resource factories.

    When gmxapi.operation Contexts provide run time inputs to operations,
    instances of this class are provided to the operation's registered
    Resource factory.

    Attributes:
        kwargs (dict): collection of named data sources.
    """

    def __init__(self, kwargs: typing.Mapping[str, SourceTypeVar]):
        self.kwargs = kwargs

    def __repr__(self):
        representation = f"{self.__class__.__name__}("
        representation += ", ".join(
            f"{key}={repr(value)}" for key, value in self.kwargs.items()
        )
        representation += ")"
        return representation


class Context(gmx.abc.Context):
    """API Context.

    All gmxapi data and operations are owned by a Context instance. The Context
    manages the details of how work is run and how data is managed.

    Context implementations are not required to inherit from gmxapi.context.Context,
    but this class definition serves to specify the current Context API.

    If subclassing is used to implement new Contexts, be sure to initialize the
    base class when providing a new __init__
    """

    def node(self, node_id) -> Node:
        if node_id in self.labels:
            return self.labels[node_id]
        elif node_id in self.work_graph:
            return self.work_graph[node_id]
        else:
            raise exceptions.ValueError(
                "Could not find a node identified by {}".format(node_id)
            )

    def __init__(self):
        self.operations = dict()
        self.labels = dict()
        self.work_graph = collections.OrderedDict()

    @abc.abstractmethod
    def node_builder(
        self, *, operation, label: typing.Optional[str] = None
    ) -> NodeBuilder:
        """Get a builder for a new work graph node.

        Nodes are elements of computational work, with resources and execution
        managed by the Context. The Context handles parallelism resources, data
        placement, work scheduling, and data flow / execution dependencies.

        This method is used by Operation director code and helper functions to
        add work to the graph.

        Arguments:
            operation: a registered gmxapi operation
            label: optional user-provided identifier to provide human-readable node locators.

        """
        ...

    # TODO: *node()* accessor.
    # @abc.abstractmethod
    # def node(self, node_identifier) -> AbstractOperation:
    #     ...


class ModuleNodeBuilder(NodeBuilder):
    """Builder for work nodes in gmxapi.operation.ModuleContext."""


class ModuleContext(Context):
    """Context implementation for the gmxapi.operation module."""

    __version__ = 0

    def node_builder(self, operation, label=None) -> NodeBuilder:
        """Get a builder for a new work node to add an operation in this context."""
        if label is not None:
            if label in self.labels:
                raise exceptions.ValueError("Label {} is already in use.".format(label))
            else:
                # The builder should update the labeled node when it is done.
                self.labels[label] = None

        return ModuleNodeBuilder(
            context=weakref.proxy(self), operation=operation, label=label
        )


# Context stack.
__current_context = [ModuleContext()]


def current_context() -> Context:
    """Get a reference to the currently active Context.

    The imported gmxapi.context module maintains some state for the convenience
    of the scripting environment. Internally, all gmxapi activity occurs under
    the management of an API Context, explicitly or implicitly. Some actions or
    code constructs will generate separate contexts or sub-contexts. This utility
    command retrieves a reference to the currently active Context.
    """
    return __current_context[-1]


def push_context(context) -> Context:
    """Enter a sub-context by pushing a context to the global context stack."""
    __current_context.append(context)
    return current_context()


def pop_context() -> Context:
    """Exit the current Context by popping it from the stack."""
    return __current_context.pop()


class OutputFactory(object):
    """Encapsulate the details of Operation output implementation in the gmxapi.operation Context.

    Currently, OutputFactory objects are containers that compose functionality
    with which to implement the required internal interface.
    """

    def __init__(
        self,
        *,
        output_proxy: typing.Callable[[SourceResource], _OutputDataProxyType],
        output_description: OutputCollectionDescription,
        publishing_data_proxy: typing.Callable[
            [SourceResource, ClientID], _PublishingDataProxyType
        ],
    ):
        """Package the output details for an operation.

        Arguments:
            output_proxy: factory to produce the *output* facet of an operation instance (node)
            output_description: fully formed output description
            publishing_data_proxy: factory to produce the run time output publishing resources

        """
        if not callable(output_proxy):
            raise exceptions.ValueError("output_proxy argument must be a callable.")
        if not callable(publishing_data_proxy):
            raise exceptions.ValueError(
                "publishing_data_proxy argument must be a callable."
            )
        if not isinstance(output_description, OutputCollectionDescription):
            raise exceptions.ValueError(
                "output_description must be an instance of "
                "gmxapi.operation.OutputCollectionDescription"
            )
        self._output_proxy = output_proxy
        self._output_description = output_description
        self._publishing_data_proxy = publishing_data_proxy

    def output_proxy(self) -> typing.Callable[[SourceResource], _OutputDataProxyType]:
        return self._output_proxy

    def output_description(self) -> OutputCollectionDescription:
        return self._output_description

    def publishing_data_proxy(
        self,
    ) -> typing.Callable[[SourceResource, ClientID], _PublishingDataProxyType]:
        return self._publishing_data_proxy


# TODO: Refactor in terms of gmx.abc.OperationDirector[_Op, gmx.operation.Context]
# Normalizing this OperationDirector may require other updates to the function_wrapper facilities.
class OperationDirector(object):
    """Direct the construction of an operation node in the gmxapi.operation module Context.

    Collaboration: used by OperationDetails.operation_director, which
    will likely dispatch to different implementations depending on
    requirements of work or context.
    """

    def __init__(
        self,
        *args,
        operation_details: typing.Type[OperationDetailsBase],
        context: Context,
        label=None,
        **kwargs,
    ):
        self.operation_details = operation_details
        self.context = weakref.proxy(context)
        self.args = args
        self.kwargs = kwargs
        self.label = label

    def __call__(self) -> AbstractOperation:
        builder = self.context.node_builder(
            operation=self.operation_details, label=self.label
        )

        builder.set_resource_factory(self.operation_details.resource_director)
        builder.set_input_description(self.operation_details)
        builder.set_handle(OperationHandle)

        # self.operation_details is a factory or class. We get an operation_details
        # instance here.
        operation_details = self.operation_details()

        # With gmxapi.abc.OperationDirector, this is the resource_factory member, which should be
        # composed instead of hard-coded to self.operation_details().signature().bind.
        node_input_factory = operation_details.signature().bind
        data_source_collection: DataSourceCollection = node_input_factory(
            *self.args, **self.kwargs
        )
        for name, source in data_source_collection.items():
            builder.add_input(name, source)

        allow_duplicate = getattr(operation_details, "_allow_duplicate", False)
        runner_director = RunnerDirector(
            runner=operation_details, allow_duplicate=allow_duplicate
        )

        builder.set_runner_director(runner_director)
        # Note that the call-backs held by OutputFactory cannot be annotated with
        # key word arguments under PEP 484, producing a weak warning in some cases.
        # We can consider more in the future how to balance call-back simplicity,
        # type annotation, and key word explicitness in helper functions like these.
        output_factory = OutputFactory(
            output_description=operation_details.output_description(),
            output_proxy=operation_details.output_data_proxy,
            publishing_data_proxy=operation_details.publishing_data_proxy,
        )
        builder.set_output_factory(output_factory)

        # The ResourceManager gets built in builder.build(). build() will only pass the
        # runner_director through to the ResourceManager runner_director argument. We
        # can replace the nested function definition with an instantiation of a
        # RunnerDirector object that we can set a *allow_duplicate* property on,
        # which we can get from operation_details. Then we can customize the creation of
        # call-backs while initializing the ResourceManager.
        handle = builder.build()
        return handle


class RunnerDirector(AbstractRunnerDirector):
    def __init__(self, runner: typing.Callable, allow_duplicate=False):
        self._runner = runner
        self.allow_duplicate = allow_duplicate

    def __call__(self, resources):
        runner = functools.partial(self._runner, resources)
        return runner


class DataStore(collections.OrderedDict, typing.Mapping[str, OutputData]):
    """A container to hold the resources for an operation node.

    Used internally by the resource manager when setting up the node.
    Evolution of the C++ framework for creating the Operation SessionResources
    object will inform the future of this and the resource_director method, but
    this data store is how the Context manages output data sources for resources
    that it manages.
    """

    _frozen = False
    _named_data_done: typing.MutableMapping[str, bool]

    def __init__(
        self,
        output_description: OutputCollectionDescription,
        ensemble_width: int,
        done_callback: typing.Optional[typing.Callable[["DataStore"], None]] = None,
    ):
        self._done_callback = done_callback

        elements = []
        for name, dtype in output_description.items():
            assert issubclass(dtype, valid_result_types)
            result_description = ResultDescription(dtype=dtype, width=ensemble_width)
            elements.append(
                (
                    name,
                    OutputData(
                        name=name,
                        description=result_description,
                        done_callback=self._mark_as_done,
                    ),
                )
            )

        collections.OrderedDict.__init__(self, elements)
        self._named_data_done = {name: False for name in self.keys()}
        self._frozen = True

    def _mark_as_done(self, key):
        """Mark the named element as done.

        Provided as a callback to the OutputData members. We collect the results
        and issue a signal of our own when the whole DataStore is done.
        """
        if isinstance(key, OutputData):
            key = key.name
        if self._named_data_done[key]:
            logger.warning(f'{self}[{key}] is being marked "done" more than once.')
        else:
            self._named_data_done[key] = True
            if all(self._named_data_done.values()) and callable(self._done_callback):
                self._done_callback(self)

    def __setitem__(self, k: str, v: OutputData) -> None:
        if self._frozen:
            raise UsageError("DataStore structure is immutable after initialization.")
        else:
            super().__setitem__(k, v)

    def __delitem__(self, v) -> None:
        if self._frozen:
            raise UsageError("DataStore structure is immutable after initialization.")
        else:
            super().__delitem__(v)

    def reset(self):
        self._named_data_done = {name: False for name in self.keys()}
        for output in self.values():
            if isinstance(output, OutputData):
                output.reset()


# TODO: For outputs, distinguish between "results" and "events".
#  Both are published to the resource manager in the same way, but the relationship
#  with subscribers is potentially different.
def function_wrapper(output: dict = None, allow_duplicate=False):
    # Suppress warnings in the example code.
    # noinspection PyUnresolvedReferences
    """Generate a decorator for wrapped functions with signature manipulation.

    New function accepts the same arguments, with additional arguments required by
    the API.

    The new function returns an object with an ``output`` attribute containing the named outputs.

    Example:

        >>> @function_wrapper(output={'spam': str, 'foo': str})
        ... def myfunc(parameter: str = None, output=None):
        ...    output.spam = parameter
        ...    output.foo = parameter + ' ' + parameter
        ...
        >>> operation1 = myfunc(parameter='spam spam')
        >>> assert operation1.output.spam.result() == 'spam spam'
        >>> assert operation1.output.foo.result() == 'spam spam spam spam'

    Arguments:
        output (dict): output names and types

    If ``output`` is provided to the wrapper, a data structure will be passed to
    the wrapped functions with the named attributes so that the function can easily
    publish multiple named results. Otherwise, the ``output`` of the generated operation
    will just capture the return value of the wrapped function.

    .. todo:: gmxapi typing stub file(s).
              The way this wrapper uses parameter annotations is not completely
              compatible with static type checking (PEP 484). If we decide to
              keep the convenience functionality by which operation details are
              inferred from parameter annotations, we should provide a separate
              stub file (.pyi) to support static type checking of the API.
    """

    if output is not None and not isinstance(output, collections.abc.Mapping):
        raise exceptions.TypeError(
            "If provided, `output` argument must be a mapping of data names to types."
        )

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
        class OperationDetails(OperationDetailsBase):
            # Warning: function.__qualname__ is not rigorous since function may be in a local scope.
            # TODO: Improve base identifier.
            # Suggest registering directly in the Context instead of in this local class definition.
            __basename = ".".join((str(function.__module__), function.__qualname__))
            __last_uid = 0
            _input_signature_description = InputCollectionDescription.from_function(
                function
            )
            # TODO: Separate the class and instance logic for the runner.
            # Logically, the runner is a detail of a context-specific implementation class,
            # though the output is not generally fully knowable until an instance is initialized
            # for a certain input fingerprint.
            # Note: We are almost at a point where this class can be subsumed into two
            # possible return types for wrapped_function_runner, acting as an operation helper.
            _runner = wrapped_function_runner(
                function=function,
                output_description=provided_output_map,
                allow_duplicate=allow_duplicate,
            )
            _output_description = _runner.output_description
            _output_data_proxy_type = define_output_data_proxy(_output_description)
            _publishing_data_proxy_type = define_publishing_data_proxy(
                _output_description
            )
            _SourceResource = SourceResource[
                _output_data_proxy_type, _publishing_data_proxy_type
            ]

            _allow_duplicate = allow_duplicate

            def __repr__(self):
                return f"<{self.__basename}(OperationDetailsBase)>"

            @classmethod
            def name(cls) -> str:
                return cls.__basename.split(".")[-1]

            @classmethod
            def namespace(cls) -> str:
                suffix = "." + cls.name()
                try:
                    index = cls.__basename.rindex(suffix)
                except ValueError:
                    index = None
                return cls.__basename[:index]

            @classmethod
            def director(cls, context: _Context):
                return cls.operation_director

            @classmethod
            def signature(cls) -> InputCollectionDescription:
                """Mapping of named inputs and input type.

                Used to determine valid inputs before an Operation node is created.

                Overrides OperationDetailsBase.signature() to provide an
                implementation for the bound operation.
                """
                return cls._input_signature_description

            def output_description(self) -> OutputCollectionDescription:
                """Mapping of available outputs and types for an existing Operation node.

                Overrides OperationDetailsBase.output_description() to provide an
                implementation for the bound operation.
                """
                return self._output_description

            def publishing_data_proxy(
                self, *, instance: SourceResource, client_id: int
            ) -> _publishing_data_proxy_type:
                """Factory for Operation output publishing resources.

                Used internally when the operation is run with resources provided by instance.

                Overrides OperationDetailsBase.publishing_data_proxy() to provide an
                implementation for the bound operation.
                """
                assert isinstance(instance, ResourceManager)
                return self._publishing_data_proxy_type(
                    instance=instance, client_id=client_id
                )

            def output_data_proxy(
                self, instance: _SourceResource
            ) -> _output_data_proxy_type:
                assert isinstance(instance, ResourceManager)
                return self._output_data_proxy_type(instance=instance)

            def __call__(self, resources: PyFunctionRunnerResources):
                """Execute the operation with provided resources.

                Resources are prepared in an execution context with aid of resource_director()

                After the first call, output data has been published and is trivially
                available through the output_data_proxy()

                Overrides OperationDetailsBase.__call__().
                """
                self._runner(resources)

            @classmethod
            def make_uid(cls, input):
                """The unique identity of an operation node tags the output with respect to the input.

                Combines information on the Operation details and the input to uniquely
                identify the Operation node.

                Arguments:
                    input : A (collection of) data source(s) that can provide Fingerprints.

                Used internally by the Context to manage ownership of data sources, to
                locate resources for nodes in work graphs, and to manage serialization,
                deserialization, and checkpointing of the work graph.

                The UID is a detail of the generic Operation that _should_ be independent
                of the Context details to allow the framework to manage when and where
                an operation is executed.

                Note:
                    This implementation creates a new identifier with every call, even if *input*
                    is the same, because we have not developed our Fingerprinting scheme in gmxapi 0.1+.

                Design notes on further refinement:
                    TODO: Operations should not single-handedly determine their own uniqueness
                    (but they should participate in the determination with the Context).

                    Context implementations should be allowed to optimize handling of
                    equivalent operations in different sessions or work graphs, but we do not
                    yet TODO: guarantee that UIDs are globally unique!

                    The UID should uniquely indicate an operation node based on that node's input.
                    We need input fingerprinting to identify equivalent nodes in a work graph
                    or distinguish nodes across work graphs.

                """
                uid = str(cls.__basename) + str(cls.__last_uid)
                cls.__last_uid += 1
                return uid

            @classmethod
            def resource_director(
                cls,
                *,
                input=None,
                output: _publishing_data_proxy_type = None,
            ) -> PyFunctionRunnerResources:
                """a Director factory that helps build the Session Resources for the function.

                The Session launcher provides the director with all resources previously
                requested/negotiated/registered by the Operation. The director uses details of
                the operation to build the resources object required by the operation runner.

                For the Python Context, the protocol is for the Context to call the
                resource_director instance method, passing input and output containers.

                Raises:
                    exceptions.TypeError if provided resource type does not match input signature.
                """
                resources = PyFunctionRunnerResources()
                resources.update(input.kwargs)
                if output:
                    resources["output"] = output

                # TODO: Remove this hack when we can better handle Futures of Containers and Future slicing.
                for name in resources:
                    if isinstance(resources[name], (list, tuple)):
                        resources[name] = datamodel.ndarray(resources[name])

                # Check data compatibility
                for name, value in resources.items():
                    if name != "output":
                        expected = cls.signature()[name]
                        got = type(value)
                        if got != expected:
                            raise exceptions.TypeError(
                                "Expected {} but got {} for {} resource {}.".format(
                                    expected, got, cls.__basename, name
                                )
                            )
                return resources

        # TODO: (FR4) Update annotations with gmxapi data types. E.g. return -> Future.
        @functools.wraps(function)
        def helper(*args, context=None, **kwargs):
            # Description of the Operation input (and output) occurs in the
            # decorator closure. By the time this factory is (dynamically) defined,
            # the OperationDetails and ResourceManager are well defined, but not
            # yet instantiated.
            # Inspection of the offered input occurs when this factory is called,
            # and OperationDetails, ResourceManager, and Operation are instantiated.

            # This operation factory is specialized for the default package Context.
            if context is None:
                context = current_context()
            else:
                raise exceptions.ApiError(
                    "Non-default context handling not implemented."
                )

            # This calls a dispatching function that may not be able to reconcile the input
            # and Context capabilities. This is the place to handle various exceptions for
            # whatever reasons this reconciliation cannot occur.
            handle = OperationDetails.operation_director(
                *args, context=context, label=None, **kwargs
            )

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


def ensure_future(
    source,
    description: typing.Optional[ResultDescription] = None,
    # name: typing.Optional[str] = None,
    _propagate_reset=True,
) -> Future:
    """Get a Future of the given description from a data source.

    If *source* is already a Future matching the provided description, return
    *source*. Otherwise, wrap *source* appropriately in a new Future instance.
    """
    # Our current convention is that unnamed data sources get called "data".
    # Set the default name of the returned Future.
    name = "data"

    # Note that the *reset* behavior for Futures is problematic and needs to be
    # replaced with a scheme by which an Operation can be duplicated or prototyped
    # or otherwise re-instantiated with new inputs in `while_loop`s, etc.
    #
    # Until then, we have to allow for Futures passed into a subgraph from outer
    # scopes to be excluded from the `reset` chain.
    if isinstance(source, Future):
        name = source.name
        source_description: ResultDescription = source.description

        if not _propagate_reset:
            logger.debug(
                f"Getting data event from {source}. Non-optimal data flow, pending work on #3139."
            )
            # WARNING: Initial implementation forces execution in order to convert Future to local constant
            # before re-wrapping. This is a potential performance issue, and interferes with evolution of
            # the execution dispatching, but is a necessary simplification for a minimal initial
            # implementation.
            #
            # If we follow up on #3139 to remove reset() and improve the prototyping of
            # new workflow nodes for a better Subgraph implementation, we can improve the situation. Please
            # raise discussion if this is important to you.
            logger.info(
                f"Localizing {source} to get a reference that will not be subject to Future.reset()."
            )
            source = source.result()
            if source_description.dtype is None:
                if source_description.width == 1:
                    _source_type = type(source)
                else:
                    assert source_description.width > 1
                    assert len(source) == source_description.width
                    _source_type = type(source[0])
                # Replace ResultDescription if we could infer something other than None (or NoneType)
                if _source_type is not None and not isinstance(None, _source_type):
                    source_description = ResultDescription(
                        dtype=_source_type, width=source_description.width
                    )
        else:
            # Future.reset() will be propagated to source, if applicable.
            # There may not be anything interesting to log yet, though.
            ...
    else:
        assert not isinstance(source, Future)
        # Resetting a constant has no effect, but for generality we do not forbid it.
        logger.debug(
            f"_propagate_reset={_propagate_reset} is ignored for non-Future data sources."
        )
        # Note that we do not currently provide a way to assert an ensemble
        # subgraph with pure constant variable initialization.
        _source_type = type(source) if source is not None else None
        source_description = ResultDescription(dtype=_source_type, width=1)

    if description is None:
        description = source_description

    target_dtype = description.dtype
    target_width = description.width
    if isinstance(source, Future):
        if target_width != source_description.width and source_description.width != 1:
            raise DataShapeError(f"{source} incompatible with {description}.")
        if target_dtype is not None:
            if source_description.dtype is not None and not issubclass(
                source_description.dtype, target_dtype
            ):
                raise GmxapiTypeError(f"{source} incompatible with {description}.")

        # We may be able to pass through or thinly wrap the Future.
        if description == source_description:
            return source

        if target_width == source_description.width:
            transform = identity
        elif target_width > 1:
            transform = functools.partial(broadcast, width=target_width)
        else:
            raise ProtocolError("Unknown bug. This line should not be reachable.")
        return ProxyResourceManager(
            proxied_future=source, width=target_width, function=transform
        ).future(name, description)

    else:
        assert not isinstance(source, Future)
        if target_width == 1:
            if target_dtype is not None and not isinstance(source, target_dtype):
                raise GmxapiTypeError(f"{source} incompatible with {description}.")
        if source_description.width == target_width:
            transform = identity
        elif source_description.width == 1:
            transform = functools.partial(broadcast, width=target_width)
        else:
            raise DataShapeError(f"Cannot map {source} to {description}.")
        future = StaticSourceManager(
            name=name, proxied_data=source, width=target_width, function=transform
        ).future(name, description)

    return future


class GraphVariableDescriptor(object):
    def __init__(self, name: str = None, dtype=None, default=None):
        self.name = name
        self.dtype = dtype
        self.default = default
        self.state = None

    @property
    def internal_name(self):
        try:
            return "_" + self.name
        except TypeError:
            return None

    def __get__(self, instance, owner):
        if instance is None:
            # Access is through the class attribute of the owning class.
            # Allows the descriptor itself to be inspected or reconfigured after
            # class definition.
            # TODO: There is probably some usage checking that can be performed here.
            return self
        try:
            value = getattr(instance, self.internal_name)
        except AttributeError:
            value = self.default
            # Lazily initialize the instance value from the class default.
            if value is not None:
                try:
                    setattr(instance, self.internal_name, value)
                except Exception as e:
                    message = (
                        "Could not assign default value to {} attribute of {}".format(
                            self.internal_name, instance
                        )
                    )
                    raise exceptions.ApiError(message) from e
        return value

    # Implementation note: If we have a version of the descriptor class with no `__set__` method,
    # it is a non-data descriptor that can be overridden by a data descriptor on an instance.
    # This could be one way to handle the polymorphism associated with state changes.
    def __set__(self, instance, value):
        if instance._editing:
            # Update the internal connections defining the subgraph.
            setattr(instance, self.internal_name, value)
        else:
            raise AttributeError("{} not assignable on {}".format(self.name, instance))


class GraphMeta(type):
    """Meta-class for gmxapi data flow graphs and subgraphs.

    Used to implement ``subgraph`` as GraphMeta.__new__(...).
    Also allows subgraphs to be defined as Python class definitions by inheriting
    from Subgraph or by using the ``metaclass=GraphMeta`` hint in the class
    statement arguments.

    The Python class creation protocol allows Subgraphs to be defined in as follows.

    See the Subgraph class documentation for customization of instances through
    the Python context manager protocol.
    """

    _prepare_keywords = ("variables",)

    @classmethod
    def __prepare__(mcs, name, bases, variables: typing.OrderedDict = None, **kwargs):
        """Prepare the class namespace.

        Keyword Args:
              variables: mapping of persistent graph variables to type / default value (optional)
        """
        # Python runs this before executing the class body of Subgraph or its
        # subclasses. This is our chance to handle key word arguments given in the
        # class declaration.

        if kwargs is not None:
            for keyword in kwargs:
                raise exceptions.UsageError(
                    "Unexpected key word argument: {}".format(keyword)
                )

        namespace = collections.OrderedDict()

        if variables is not None:
            if isinstance(variables, collections.abc.Mapping):
                for name, value in variables.items():
                    if isinstance(value, type):
                        dtype = value
                        if hasattr(value, "default"):
                            default = value.default
                        else:
                            default = None
                    else:
                        default = value
                        if hasattr(default, "dtype"):
                            dtype = default.dtype
                        else:
                            dtype = type(default)
                    namespace[name] = GraphVariableDescriptor(
                        name, default=default, dtype=dtype
                    )
                    # Note: we are not currently using the hook used by `inspect`
                    # to annotate with the class that defined the attribute.
                    # namespace[name].__objclass__ = mcs
                    assert not hasattr(namespace[name], "__objclass__")
            else:
                raise exceptions.ValueError(
                    '"variables" must be a mapping of graph variables to types or defaults.'
                )

        return namespace

    def __new__(cls, name, bases, namespace, **kwargs):
        for key in kwargs:
            if key not in GraphMeta._prepare_keywords:
                raise exceptions.ApiError(
                    "Unexpected class creation keyword: {}".format(key)
                )
        return type.__new__(cls, name, bases, namespace)

    # TODO: This is keyword argument stripping is not necessary in more recent Python versions.
    # When Python minimum required version is increased, check if we can remove this.
    def __init__(cls, name, bases, namespace, **kwargs):
        for key in kwargs:
            if key not in GraphMeta._prepare_keywords:
                raise exceptions.ApiError(
                    "Unexpected class initialization keyword: {}".format(key)
                )
        super().__init__(name, bases, namespace)


class SubgraphNodeBuilder(NodeBuilder):
    def __init__(
        self, context: "SubgraphContext", operation, label: typing.Optional[str] = None
    ):
        super().__init__(context, operation, label)

    def add_input(self, name: str, source):
        """Add an input resource for the Node under construction.

        Extends NodeBuilder.add_input()
        """
        # Inspect inputs.
        #  * Are they from outside the subgraph?
        #  * Subgraph variables?
        #  * Subgraph internal nodes?
        # Inputs from outside the subgraph are (provisionally) subgraph inputs.
        # Inputs that are subgraph variables or subgraph internal nodes mark operations that will need to be re-run.
        # For first implementation, let all operations be recreated, but we need to
        # manage the right input proxies.
        # For zeroeth implementation, try just tracking the entities that need a reset() method called.
        assert isinstance(self.context, SubgraphContext)
        if hasattr(source, "reset"):
            self.context.add_resetter(source.reset)
        elif hasattr(source, "_reset"):
            self.context.add_resetter(source._reset)
        super().add_input(name, source)

    def build(self) -> OperationPlaceholder:
        """Get a reference to the internal node in the subgraph definition.

        In the SubgraphContext, these handles cannot represent uniquely identifiable
        results. They are placeholders for relative positions in graphs that are
        not defined until the subgraph is being executed.

        Such references should be tracked and invalidated when exiting the
        subgraph context. Within the subgraph context, they are used to establish
        the recipe for updating the subgraph's outputs or persistent data during
        execution.
        """
        # Placeholder handles in the subgraph definition don't have real resource managers.
        # Check for the ability to instantiate operations.
        handle = super().build()
        # handle = OperationPlaceholder()
        return typing.cast(OperationPlaceholder, handle)


class SubgraphContext(Context):
    """Provide a Python context manager in which to set up a graph of operations.

    Allows operations to be configured without adding nodes to the global execution
    context.
    """

    def __init__(self):
        super().__init__()
        self.resetters = set()

    def node_builder(self, operation, label=None) -> NodeBuilder:
        if label is not None:
            if label in self.labels:
                raise exceptions.ValueError("Label {} is already in use.".format(label))
            else:
                # The builder should update the labeled node when it is done.
                self.labels[label] = None

        return SubgraphNodeBuilder(
            context=weakref.proxy(self), operation=operation, label=label
        )

    def add_resetter(self, function):
        assert callable(function)
        self.resetters.add(function)


class Subgraph(object, metaclass=GraphMeta):
    """

    When subclassing from Subgraph, aspects of the subgraph's data ports can be
    specified with keyword arguments in the class statement. Example::

        >>> class MySubgraph(Subgraph, variables={'int_with_default': 1, 'boolData': bool}): pass
        ...

    The key word *variables* is used in the class declaration to map the types
    of subgraph Variables with (optional) default values.

    Execution model:
        Subgraph execution must follow a well-defined protocol in order to sensibly
        resolve data Futures at predictable points. Note that subgraphs act as operations
        that can be automatically redefined at run time (in limited cases), such as
        to automatically form iteration chains to implement "while" loops. We refer
        to one copy / iteration / generation of a subgraph as an "element" below.

        When a subgraph element begins execution, each of its variables with an
        "updated" state from the previous iteration has the "updated" state moved
        to the new element's "initial" state, and the "updated" state is voided.

        Subgraph.next() appends an element of the subgraph to the chain. Subsequent
        calls to Subgraph.run() bring the new outputs up to date (and then call next()).
        Thus, for a subgraph with an output called ``is_converged``, calling
        ``while (not subgraph.is_converged): subgraph.run()`` has the desired effect,
        but is likely suboptimal for execution. Instead use the gmxapi while_loop.

        If the subgraph is not currently executing, it can be in one of two states:
        "editing" or "ready". In the "ready" state, class and instance Variables
        cannot be assigned, but can be read through the data descriptors. In this
        state, the descriptors only have a single state.

        If "editing," variables on the class object can be assigned to update the
        data flow defined for the subgraph. While "editing", reading or writing
        instance Variables is undefined behavior.

        When "editing" begins, each variable is readable as a proxy to the "initial"
        state in an element. Assignments while "editing" put variables in temporary
        states accessible only during editing.
        When "editing" finishes, the "updated" output data source of each
        variable for the element is set, if appropriate.

        # TODO: Use a get_context to allow operation factories or accessors to mark
        #  references for update/annotation when exiting the 'with' block.
    """


class SubgraphBuilder(object):
    """Helper for defining new Subgraphs.

    Manages a Python context in which to define a new Subgraph type. Can be used
    in a Python ``with`` construct exactly once to provide the Subgraph class body.
    When the ``with`` block is exited (or ``build()`` is called explicitly), a
    new type instance becomes available. Subsequent calls to SubgraphBuilder.__call__(self, ...)
    are dispatched to the Subgraph constructor.

    Outside of the ``with`` block, read access to data members is proxied to
    descriptors on the built Subgraph.

    Instances of SubgraphBuilder are returned by the ``subgraph()`` utility function.
    """

    def __init__(self, variables: collections.abc.Mapping):
        # Instances have different attribute access behavior depending whether the subgraph is actively
        # being edited or not, and we override the attribute accessors to provide a different view of the
        # namespace than what is really under the hood. This is confusing and needs to change. If we stay
        # with the current pattern, SubgraphBuilder would be sensibly split into a separate SubgraphManager
        # and SubgraphEditor class, at least. However, we should instead focus on reimplementing such that
        # Subgraph (sub)classes can be defined in terms of Data Descriptors (the class definition block
        # takes the place of the `with` block).
        self.__dict__.update(
            {
                "variables": dict(),
                "_editing": False,
                "_factory": None,
                "_fused_operation": None,
                "_initial_width": 1,
                "_internal_state": dict(),
                "_staging": collections.OrderedDict(),
                "_subgraph_context": None,
                "_subgraph_instance": None,
            }
        )
        # Establish the initial ensemble width. (If all variables have width 1,
        # we may widen the ensemble when creating the Subgraph instance.)
        for name, value in variables.items():
            if isinstance(value, Future):
                self._update_width(value.description.width)
            self.__dict__["variables"][name] = ensure_future(
                source=value, _propagate_reset=False
            )

        self._subgraph_instance = None

    def _update_width(self, width: int):
        if self._initial_width == 1:
            self._initial_width = width
        elif width != 1 and width != self._initial_width:
            raise DataShapeError(
                f"Cannot map source of width {width} to subgraph width {self._initial_width}."
            )
        return self._initial_width

    def __getattr__(self, item):
        # Note that obj.__getattr__ is called only if *item* is not an existing attribute.
        if self._editing:
            # This is where we are staging variables updates while initially defining the subgraph. The
            # instance has not yet been created.
            if item in self.variables:
                if item in self._staging:
                    logger.debug(
                        f"Read access to intermediate value of subgraph variable {item}"
                    )
                    return self._staging[item]
                else:
                    logger.debug("Read access to subgraph variable {}".format(item))
                    # Note that the Future returned here needs to represent a point in the internal data
                    # flow of the subgraph, but if the variable is assigned to later in the definition,
                    # the subgraph instance needs to make sure that this reference is updated during
                    # Subgraph.run() to point to the correct resource at the next iteration.
                    if item not in self._internal_state:
                        self._internal_state[item] = self.variables[item]
                    else:
                        assert self._internal_state[item] is self.variables[item]
                    return self._internal_state[item]
            else:
                raise AttributeError("Invalid attribute: {}".format(item))
        else:
            # The variable is being accessed outside the scope of the definition process (the `with`
            # block). Usually this represents either direct access by a user or internal usage by the
            # `while_loop` implementation.
            if not self._subgraph_instance:
                logger.debug(
                    f"Read access to initial value of {item} through before building subgraph."
                )
                return self.variables[item]
            # TODO: this is not quite the interface described for operation.Subgraph or gmxapi.subgraph...
            return SubgraphVariable(subgraph=self._subgraph_instance, item=item)
        # Conceivably, there is different preferred behavior for subgraph variable reading through the Builder
        #    * before starting to define the subgraph data flow
        #    * during subgraph definition (editing)
        #    * after finalizing the subgraph definition, but before executing
        #    * during execution of the while_loop wrapper
        #    * during execution of the lower-level Subgraph.run()
        #    * after Subgraph.run() has executed (at least once)
        #    * once the while_loop wrapper has finished
        # This is too much responsibility for one function. Most of these use cases are not currently well
        # specified. But they illustrate that we need a clearer distinction between the participants of
        # subgraph specification (its Builder, type, or metatype), instantiation of subgraph frames,
        # and use in iteration. Instead of stateful access logic, we should clarify distinct behavioral
        # protocols, as needed, and try to minimize misuse or accessibility of unspecified access in the
        # meantime.

    def __setattr__(self, key, value):
        """Part of the builder interface."""
        if key in self.__dict__:
            logger.debug(
                f"Replacing {self}.{key} with {value}. (Was {self.__dict__[key]}.)"
            )
            self.__dict__[key] = value
        else:
            if self._editing:
                self.add_update(key, value)
            else:
                raise exceptions.UsageError("Subgraph is not in an editable state.")

    def add_update(self, key, value):
        """Add a variable update to the internal subgraph.

        The new value must be consistent with the existing type, but is allowed to
        convert a single-element variable to an ensemble variable.
        """
        # To avoid potential confusion, do this check before triggering the conditionals
        # in __getattr__ or __setattr__.
        if not self._editing:
            raise exceptions.UsageError("Subgraph is not in an editable state.")
        # Else, stage the potential final value for the iteration.

        if key not in self.variables:
            raise AttributeError("No such attribute: {}".format(key))

        if not isinstance(self.variables[key], Future):
            raise ApiError(
                "SubgraphBuilder initialization should have normalized Futures for variables. "
                f"Found {key}: {self.variables[key]}"
            )

        target_description: ResultDescription = self.variables[key].description
        if not isinstance(value, Future):
            warnings.warn(
                f"Updating Subgraph.{key} with {value}, which is not a Future and will not be recalculated "
                f"at each iteration. If this is intentional, please contact the GROMACS team to help us "
                f"understand your use case, so we can replace this warning with a more appropriate check."
            )
        else:
            assert isinstance(value, Future)
            if (
                target_description.dtype is not None
                and value.description.dtype is not None
            ):
                if not issubclass(value.description.dtype, target_description.dtype):
                    raise ApiError(
                        f"Cannot update {key}: {target_description} from {value}."
                    )
            target_description._width = self._update_width(
                width=value.description.width
            )
        value = ensure_future(value, target_description)
        logger.debug(f"Staging subgraph update to {key} from {value}.")

        # Set a placeholder that we can update during iteration.
        # Long term, this is probably implemented with data descriptors
        # that will be moved to a new Subgraph type object.
        self._staging[key] = value
        self._staging.move_to_end(key)

    def __enter__(self):
        """Enter a Context managed by the subgraph to capture operation additions.

        Allows the internal data flow of the subgraph to be defined in the same
        manner as the default work graph while the Python context manager is active.

        The subgraph Context is activated when entering a ``with`` block and
        finalized at the end of the block.
        """
        # While editing the subgraph in the SubgraphContext, we need __get__ and __set__
        # data descriptor access on the Subgraph subclass type, but outside of the
        # context manager, the descriptor should be non-writable and more opaque,
        # while the instance should be readable, but not writable.
        self.__dict__["_editing"] = True
        # TODO: this probably needs to be configured with variables...
        self.__dict__["_subgraph_context"] = SubgraphContext()
        push_context(self._subgraph_context)
        return self

    def build(self):
        """Build the subgraph by defining some new operations.

        Examine the subgraph variables. Variables with handles to data sources
        in the SubgraphContext represent work that needs to be executed on
        a subgraph execution or iteration. Variables with handles to data sources
        outside of the subgraph represent work that needs to be executed once
        on only the first iteration to initialize the subgraph.

        Construct a factory for the fused operation that performs the work to
        update the variables on a single iteration.

        Construct a factory for the fused operation that performs the work to
        update the variables on subsequent iterations and which will be fed
        the outputs of a previous iteration. Both of the generated operations
        have the same output signature.

        Construct and wrap the generator function to recursively append work
        to the graph and update until condition is satisfied.

        TODO: Explore how to drop work from the graph once there are no more
         references to its output, including check-point machinery.
        """
        logger.debug("Finalizing subgraph definition.")

        # Updates to the subgraph output variables, as staged by assignments made while defining the subgraph.
        updates = self._staging

        # References (potentially) held by operations in the subgraph. These are references that were
        # initialized from input variable values, but which need to be updated for iterations after the
        # first, if the subgraph variable has an update staged in the subgraph specification.
        internal_state = self._internal_state

        # We can use the width inferred by the Builder to set the width for the
        # instance we will create.
        initial_width = self._initial_width

        # Whether we allow updates to the Subgraph after creation is a separate and uncoupled design
        # question, but the solution implemented in this closure assumes subgraph width and while_loop
        # width will remain *initial_width*

        # Establish the input data flow for the subgraph variables.
        inputs = collections.OrderedDict()
        for name, value in self.variables.items():
            # Variables *should* be initialized with clearly typed initial values or Futures.
            if isinstance(value, Future):
                dtype = value.description.dtype
            else:
                dtype = type(value)
            if dtype is type(None):
                dtype = None
            # Since Futures are not yet robustly typed, we may have to try to infer type from the staged
            # update.
            if dtype is None:
                update = updates.get(name, None)
                if isinstance(update, Future):
                    dtype = update.description.dtype
            logger.debug(
                "Inferred type {} for variable {}".format(
                    "None" if dtype is None else dtype.__class__.__name__, name
                )
            )
            description = ResultDescription(dtype=dtype, width=initial_width)

            # Get a placeholder that we can update during iteration.
            # Long term, this is probably implemented with data descriptors
            # that will be moved to a new Subgraph type object.
            inputs[name] = ensure_future(
                source=value, description=description, _propagate_reset=False
            )

        class Subgraph(object):
            # Note: Historically, it makes sense to distinguish fused operations as distinct types,
            # but there is not a technical reason why this class needs to be local to this closure.
            # Readability would likely be improved by refining the protocol for composing subgraph
            # definitions in terms of module-level classes.
            def __init__(
                self,
                *,
                input_variables: collections.OrderedDict,
                output_variables: collections.OrderedDict,
                state_variables: dict,
            ):
                # Current (instantaneous) values of subgraph variables.
                self.values = collections.OrderedDict(
                    [
                        (_name, _source.result())
                        for _name, _source in input_variables.items()
                    ]
                )
                self.value_width = {key: None for key in self.values}
                logger.debug(
                    "subgraph initialized with {}".format(
                        ", ".join(
                            [
                                "{}: {}".format(_name, _source)
                                for _name, _source in self.values.items()
                            ]
                        )
                    )
                )

                # Initialize the Futures supporting the output variables.
                self.futures = collections.OrderedDict(
                    [(_name, _source) for _name, _source in input_variables.items()]
                )

                # Stash the sequenced updates to the output variables, as gathered from assignments made
                # while editing the subgraph.
                self.output_variable_updates = collections.OrderedDict(
                    [(_name, _source) for _name, _source in output_variables.items()]
                )
                if len(self.output_variable_updates) > 0:
                    logger.debug("Subgraph updates staged:")
                    for _name, source in self.output_variable_updates.items():
                        logger.debug("    {} = {}".format(_name, source))

                # Operations in the subgraph are referencing subgraph variables that are updated
                # when the subgraph is run. Make sure these Futures get the same updates as the
                # Futures supporting the output variables, but only after the operation results are
                # consumed within the iteration.
                self.internal_references = {
                    _name: _future
                    for _name, _future in state_variables.items()
                    if _name in self.output_variable_updates
                }

                # Cache for representation. (initially "expired")
                self.__representation = None

            def _update_representation(self):
                # Developer note: expire the representation cache when updating *values* or
                # *futures* members.
                # Cached representation depends on self.values and self.update_sources. This feature is for
                # debugging, and these details are pretty fluid. We could, of course, use a fancier data
                # model with callbacks to expire the cache as needed, but that would be needless
                # obfuscation and would probably destroy the performance benefit of the cached
                # representation, anyway.
                representation = "<Subgraph: values={"
                representation += ", ".join(
                    f"{key}={repr(value)}" for key, value in self.values.items()
                )
                representation += "}, futures={"
                representation += ", ".join(
                    f"{key}={repr(value)}" for key, value in self.futures.items()
                )
                representation += "}>"
                self.__representation = representation
                return representation

            def __repr__(self):
                representation = self.__representation
                if not representation:
                    representation = self._update_representation()
                return representation

            def run(self):
                # Bring the stored Futures up to date for the current iteration.
                # For Futures that depend on subgraph-internal data flow, re-stage updates for the next
                # iteration.
                futures_updates = {}
                for name in self.output_variable_updates:
                    _future: Future = self.output_variable_updates[name]
                    source_description = _future.description
                    try:
                        logger.debug(
                            f"Getting update from {_future.resource_manager}[{_future.name}]."
                        )
                        result = _future.result()
                    except Exception as e:
                        logger.error(
                            f"Could not update {name}. Got exception: ", exc_info=e
                        )
                        raise e
                    logger.debug("Update: {} = {}".format(name, result))
                    self.value_width[name] = source_description.width
                    self.values[name] = result
                    # Expire the cached representation.
                    self.__representation = None

                    # Get a new ResourceManager for the Future serving the subgraph variable (for the next
                    # iteration).
                    # TODO: If the width changes, we need to serve the correct member.
                    # if source_description.width > 1 and self.futures[name].description.width == 1 or
                    #         self.internal_references[name].description.width == 1: transform = slice_by_rank
                    new_source = StaticSourceManager(
                        name=_future.name,
                        proxied_data=result,
                        width=source_description.width,
                        function=identity,
                    )
                    futures_updates[name] = {
                        "key": _future.name,
                        "description": _future.description,
                        "resource_manager": new_source,
                    }
                # Replace the data sources in the futures.
                for name in self.output_variable_updates:
                    attrs = futures_updates[name]
                    _output: Future = self.futures[name]
                    _output.name = attrs["key"]
                    _output.resource_manager = attrs["resource_manager"]
                    _output.description = attrs["description"]
                    # Expire the cached representation.
                    self.__representation = None
                    if name in self.internal_references:
                        _internal: Future = self.internal_references[name]
                        # The following assertions test the author's assumptions, not specified behavior.
                        assert _internal is not _output
                        assert _internal is not futures_updates[name]
                        # WARNING: Need to make sure things behave as expected if we change the width after
                        # initial binding.
                        if _internal.description != _output.description:
                            raise ProtocolError(
                                f"Attempting to update {_internal} with incompatible source: {_output}"
                            )
                        _internal.resource_manager = self.futures[name].resource_manager
                        _internal.name = self.futures[name].name
                for _future in self.output_variable_updates.values():
                    _future._reset()
                    self.__representation = None

        self._subgraph_instance = Subgraph(
            input_variables=inputs,
            output_variables=updates,
            state_variables=internal_state,
        )

        def factory():
            # Imitate the creation pattern for OperationHandles, but hold the (pre)built Subgraph in the
            # `build()` closure, rather than creating when the factory is called. This is a convenient and
            # minimal solution while the Subgraph class definition encapsulates the inputs. A more complete
            # Subgraph implementation should allow input descriptions to be specified in the definition,
            # but bound to actual inputs at instantiation to make Subgraphs more useful as general reusable
            # fused operations.
            return self._subgraph_instance

        return factory

    def __exit__(self, exc_type, exc_val, exc_tb):
        """End the Subgraph editing session and finalize the Subgraph build.

        After exiting, this instance forwards __call__() to a factory for an
        operation that carries out the work in the subgraph with inputs bound
        in the current context as defined by ``variables``.
        """
        self._factory = self.build()

        context = pop_context()
        assert context is self._subgraph_context
        self.__dict__["_editing"] = False
        # Returning False causes exceptions in the `with` block to be reraised.
        # Remember to switch this to return True if we want to transform or suppress
        # such an exception (we probably do).
        if exc_type is not None:
            logger.error(
                "Got exception {} while editing subgraph {}.".format(exc_val, self)
            )
            logger.debug("Subgraph exception traceback: \n{}".format(exc_tb))
        return False

    def __call__(self):
        # After build() has been called, SubgraphBuilder is callable and acts as a
        # factory for the subgraph instance.
        # TODO: the factory should return a proper OperationHandle.
        return self._factory()


class SubgraphVariable:
    """Provide access to the instantaneous value of a subgraph variable."""

    # Conundrum: For user access, we would expect a Future interface. To support the `while_loop`,
    # we expect to use the instance to get a peek at the value of the internal variable at the end of each
    # iteration.
    def __init__(self, subgraph, item):
        # TODO: Move nested Subgraph class to module level and add type checking.
        assert subgraph.__class__.__name__ == "Subgraph"
        self._subgraph = subgraph
        self._item = item

    def __call__(self, obj):
        """Condition callback.

        Support the while_loop *condition* callback. This protocol is a work in progress and should not be
        relied upon outside the `while_loop` implementation.
        """
        # This behavior is a work in progress. The original concept was akin to a Data Descriptor. The
        # while_loop would hold the Data Descriptor for a derived Subgraph class object and then provide
        # the instance at run time. This is probably still the scheme we want, but currently the while_loop
        # gets an attribute of an instance of a SubgraphBuilder (which is not a class object and not part
        # of the subgraph instance's MRO in any way).
        # Since we currently create a single subgraph instance for the subgraph builder, it seems redundant
        # to hold a reference to the owning subgraph at instantiation _and_ to accept the instance at the
        # time of reading. However, this is part of an evolution towards a model in which subgraphs look
        # more like regular Operations, with a base class instead of the current SubgraphBuilder,
        # and where the current class becomes a Data Descriptor for specific Subgraphs.
        if obj is not self._subgraph:
            raise ProtocolError(
                f"Received {obj} for a callback expecting {self._subgraph}."
            )
        return obj.values[self._item]


def while_loop(*, operation, condition, max_iteration=10):
    """Generate and run a chain of operations such that condition evaluates True.

    Returns and operation instance that acts like a single node in the current
    work graph, but which is a proxy to the operation at the end of a dynamically generated chain
    of operations. At run time, condition is evaluated for the last element in
    the current chain. If condition evaluates False, the chain is extended and
    the next element is executed. When condition evaluates True, the object
    returned by ``while_loop`` becomes a proxy for the last element in the chain.

    Equivalent to calling operation.while(condition), where available.

    Arguments:
        operation: a callable that produces an instance of an operation when called with no arguments.
        condition: a callable that accepts an object (returned by ``operation``) that returns a boolean.
        max_iteration: execute the loop no more than this many times (default 10)

    Warning:
        *max_iteration* is provided in part to minimize the cost of bugs in early
        versions of this software. The default value may be changed or
        removed on short notice.

    Warning:
        The protocol by which ``while_loop`` interacts with ``operation`` and ``condition``
        is very unstable right now. Please refer to this documentation when installing new
        versions of the package.

    Protocol:
        Warning:
            This protocol will be changed before the API is finalized.

        When called, ``while_loop`` calls ``operation`` without arguments
        and captures the return value for inspection of outputs.
        The object produced by ``operation()`` must have a ``reset``,
        a ``run`` method, and an ``output`` attribute. From gmxapi 0.1, an additional
        ``values`` attribute is examined to advertise the output members that will
        appear on the ``while_loop`` output.

        This is inspected
        to determine the output data proxy for the operation produced by the call
        to ``while_loop``. When that operation is called, it does the equivalent of

            while(condition(self._operation)):
                self._operation.reset()
                self._operation.run()

        Then, the output data proxy of ``self`` is updated with the results from
        self._operation.output.

    """
    # In the first implementation, Subgraph is NOT and OperationHandle.
    # if not isinstance(obj, AbstractOperationHandle):
    #     raise exceptions.UsageError(
    #     '"operation" key word argument must be a callable that produces an Operation handle.')
    # outputs = {}
    # for name, descriptor in obj.output.items():
    #     outputs[name] = descriptor._dtype

    # 1. Get the initial inputs.
    # 2. Initialize the subgraph with the initial inputs.
    # 3. Run the subgraph.
    # 4. Get the outputs.
    # 5. Initialize the subgraph with the outputs.
    # 6. Go to 3 if condition is not met.

    # TODO: Add type hint after moving Subgraph to module level.
    obj = operation()
    assert hasattr(obj, "values")
    # Note: Usually, we derive the ensemble width of an operation handle from the
    # inputs, but we have encapsulated the inputs in the subgraph and don't have
    # that information here. `while_loop` implicitly does an ensemble `gather`.
    # However, we can inspect the subgraph at this point to revise the output
    # descriptions.
    # TODO: Conditionally handle Subgraph versus standard OperationHandle.
    outputs = dict()
    for key, value in obj.futures.items():
        assert isinstance(value, Future)
        if value.description.width > 1:
            outputs[key] = NDArray
        else:
            outputs[key] = value.description.dtype

    def evaluate_condition():
        result = condition(obj)
        if isinstance(result, bool):
            return result
        elif isinstance(result, collections.abc.Sequence):
            if not isinstance(result, str):
                return all(result)
        else:
            raise ApiError("while_loop *condition* must produce boolean values.")

    @function_wrapper(output=outputs, allow_duplicate=True)
    def run_loop(output: OutputCollectionDescription):
        iteration = 0
        obj = operation()
        message = f"Created {repr(obj)} containing "
        values = ", ".join([f"{key}: {obj.values[key]}" for key in obj.values])
        if values:
            message += values
        else:
            message += "no values"
        logger.debug(message)
        logger.debug("Condition: {}".format(condition(obj)))

        while evaluate_condition():
            logger.debug("Running iteration {}".format(iteration))
            # WARNING: This is a dynamically generated task that can break logic based on
            # assumptions about the sequence of execution or resolution of data Futures,
            # if we aren't careful.
            obj.run()
            logger.debug(
                ", ".join(["{}: {}".format(key, obj.values[key]) for key in obj.values])
            )
            logger.debug("Condition: {}".format(condition(obj)))
            iteration += 1
            if iteration > max_iteration:
                logger.warning(
                    f"while_loop condition still True after {iteration} iterations. "
                    "If you are confident that the script is behaving correctly, "
                    f"but needs more than {max_iteration} iterations, "
                    "set a higher value for the *max_iteration* argument to `while_loop`"
                )
                break
        for name in outputs:
            setattr(output, name, obj.values[name])

    return run_loop


def subgraph(variables: collections.abc.Mapping = None):
    """Allow operations to be configured in a sub-context.

    The object returned functions as a Python context manager. When entering the
    context manager (the beginning of the ``with`` block), the object has an
    attribute for each of the named ``variables``. Reading from these variables
    gets a proxy for the initial value or its update from a previous loop iteration.
    At the end of the ``with`` block, any values or data flows assigned to these
    attributes become the output for an iteration.

    After leaving the ``with`` block, the variables are no longer assignable, but
    can be called as bound methods to get the current value of a variable.

    When the object is run, operations bound to the variables are ``reset`` and
    run to update the variables.
    """
    # Implementation note:
    # A Subgraph (type) has a subgraph context associated with it. The subgraph's
    # ability to capture operation additions is implemented in terms of the
    # subgraph context.
    logger.debug("Declare a new subgraph with variables {}".format(variables))

    return SubgraphBuilder(variables)


@function_wrapper()
def _join_arrays(
    *, front: datamodel.NDArray = (), back: datamodel.NDArray = ()
) -> datamodel.NDArray:
    # TODO: elaborate and clarify.
    # TODO: check type and shape.
    # TODO: figure out a better annotation.
    if isinstance(front, (str, bytes)) or isinstance(back, (str, bytes)):
        raise exceptions.ValueError("Input must be a pair of lists.")
    assert isinstance(front, datamodel.NDArray)
    assert isinstance(back, datamodel.NDArray)
    new_list = list(front._values)
    new_list.extend(back._values)
    return datamodel.NDArray(new_list)


def join_arrays(
    *, front: datamodel.NDArray = (), back: datamodel.NDArray = ()
) -> Future[datamodel.NDArray]:
    """Consumes two sequences and produces a concatenated single sequence.

    Note that the exact signature of the operation is not determined until this
    helper is called. Helper functions may dispatch to factories for different
    operations based on the inputs. In this case, the dtype and shape of the
    inputs determines dtype and shape of the output. An operation instance must
    have strongly typed output, but the input must be strongly typed on an
    object definition so that a Context can make runtime decisions about
    dispatching work and data before instantiating.
    """
    return _join_arrays(front=front, back=back).output.data


# TODO: Constrain
Scalar = typing.TypeVar("Scalar")


def concatenate_lists(
    sublists: typing.Sequence[typing.Sequence[Scalar]] = (),
) -> ArrayFuture[Scalar]:
    """Combine data sources into a single list.

    A trivial data flow restructuring helper.
    """
    if isinstance(sublists, (str, bytes)):
        raise exceptions.ValueError("Input must be a list of lists.")
    if len(sublists) == 0:
        return datamodel.ndarray([])
    else:
        # TODO: Fix the data model so that this can type-check properly.
        _front = sublists[0]
        _back = concatenate_lists(sublists[1:])
        return _join_arrays(
            front=_front, back=typing.cast(datamodel.NDArray, _back)
        ).output.data


def make_constant(value: Scalar) -> GenericFuture[Scalar]:
    """Provide a predetermined value at run time.

    This is a trivial operation that provides a (typed) value, primarily for
    internally use to manage gmxapi data flow.

    Accepts a value of any type. The object returned has a definite type and
    provides same interface as other gmxapi outputs. Additional constraints or
    guarantees on data type may appear in future versions.
    """
    dtype = type(value)
    source = StaticSourceManager(
        name="data", proxied_data=value, width=1, function=lambda x: x
    )
    description = ResultDescription(dtype=dtype, width=1)
    future = source.future("data", description=description)
    return future


def logical_not(value):
    """Boolean negation.

    If the argument is a gmxapi compatible Data or Future object, a new View or
    Future is created that proxies the boolean opposite of the input.

    If the argument is a callable, logical_not returns a wrapper function that
    produces the logical opposite of the result of the callable. If the callable
    produces a (non-string) sequence, the wrapper returns a list of the negated
    results of the callable.
    """
    if callable(value):

        def negate(arg):
            result = value(arg)
            if isinstance(result, collections.abc.Sequence) and not isinstance(
                result, str
            ):
                return list(not item for item in result)
            else:
                return not result

        return negate

    # TODO: Small data transformations like this don't need to be formal Operations.
    # This could be essentially a data annotation that affects the resolver in a
    # DataEdge. As an API detail, coding for different Contexts and optimizations
    # within those Context implementations could be simplified.
    operation = function_wrapper(output={"data": bool})(
        lambda data=bool(): not bool(data)
    )
    return operation(data=value).output.data
