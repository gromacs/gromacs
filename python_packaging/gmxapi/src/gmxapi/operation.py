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

    # make_operation is a utility used to create the `noop` module attribute
    # in the importable `test_support` module that copies its `data` input to its output.
    from test_support import noop
    op = noop(input={'data': True})
    assert op.output.data.extract() == True
    op = noop(input={'data': False})
    assert op.output.data.extract() == False

    # the `count` operation copies its `data` input to its output
    # and increments its `count` input to its output.
    from test_support import count
    op1 = count(input={'data': True})
    op2 = count(input=op1.output)
    op2.run()
    # TBD: how to introspect execution dependency or not-yet-executed status?
    # To allow introspection in testing, we might use module global data for operations executed in the local process.

    op1 = count(input={'data': False, 'count': 0})
    op2 = count(input=op1.output)
    assert op2.output.count.extract() == 2
"""

import abc
import functools
import inspect

import gmxapi.exceptions

class NamedResult(object):
    """Descriptor class for Result behavior on attributes.

    Stores actual data in the owning object, but ensures that the object
    accessed through the named port has the expected behavior.
    """
    def __init__(self, implementation, input):
        # The name of the public attribute in the owning class.
        self.name = None
        # The name of the private attribute in the owning class.
        self.internal_name = None
        self.implementation = Result(implementation, input)

    def __get__(self, instance, instance_type):
        # Note: `self` is the object providing an attribute accessor for an
        # `instance` of the class `instance_type`
        # (`instance` is None if attribute is accessed through the owning class).
        if instance is None: return self.implementation
        return self.implementation

# TODO: let Output be a Result, optionally containing a structure of Results
class AbstractResult(abc.ABC):
    """A Result serves as a proxy to or "future" for data produced by gmxapi operations.

    Result is the subset of the gmxapi Operation interface providing the output data proxy.

    An object implementing the Result interface has a `result()` method that forces resolution
    of data dependencies, triggering any necessary computation and communication, to return
    an object of the data type represented by the Result. If the Result is a container for other
    Results, nested Results are also resolved to local concrete data.
    """

    # TODO: A Result should defer implementation details to the Context or parent operation.
    # Note: A Result instance should not need to hold more than a weakref, if even that, to a
    # Context instance, and no reference to an Operation instance.

    #@property
    # @abc.abstractmethod
    @property
    def uid(self):
        """Get a unique identifier by which this result can be identified.

        Value should be distinct for logically different input parameters (that would be expected
        to produce non-equivalent results.

        Value must be equal for any situation in which the Context should or must reuse resources
        to acquire the result.

        Value should be equal for logically equivalent inputs that would be expected to produce
        equivalent results (to allow duplication of work or data to be avoided).

        Value may be equal but is not required to be equal in different environments or for different
        output parameters, in which different artifacts may be produced. The use case would be to
        allow greater data portability without explicit conversion between contexts.

        Can be used by a Context implementation to determine whether data is available locally,
        manage checkpointing and caching, map inputs to outputs, etc.

        The `uid` property is likely to be a constant string that is set on creation of the Result object.
        """
        return False

    # TODO:
    #"""An object implementing the Result interface (including a pure Result object) has
    #an `output` attribute that gets a pure Result copy of the Result handle.
    #"""

    # @property
    # @abc.abstractmethod
    # def output(self):
    #     return self.resulttype(self)

    # TODO: A Result needs to provide a way to get/validate a handle in the (new) Context of a subscriber.

    # TODO: We might reformulate Result as derived from concurrent.futures.Future
    # but it is currently conceptually dissimilar in that the availability of a Result
    # does not depend on having `submit`ed a task, but rather on describing tasks to be submitted.
    def result(self):
        """Get the represented data as a local object."""
        # Maybe a generic implementation would look like this.
        # TODO: Check for valid resource manager.
        # TODO: Check whether uid is registered with resource manager
        # TODO: do we need to provide a reference to this object and/or its class
        #       to allow Context to lazily create and perform the operation at call time?
        self._resource_manager.get(self.uid)

    @property
    @abc.abstractmethod
    def dtype(self):
        """Base data type of the result.

        Used to determine compatibility with the mapped inputs of consuming operations.
        """
        # At any point in time, the resource represented by this result may be in some abstract state
        # and may by an array of data sources that will be scattered to consumers.
        return None

    def _subscribe(self, subscriber):
        """Register interest in being provided with the result data when it is available.

        I don't think this is quite what we want. We want the consumer to register its interest
        through the resource management of the subscriber's Context, which can negotiate with
        the context of this object, probably by getting a Director from it.
        """
        pass

    @classmethod
    def _director(self, input, context):
        """Get an operation director appropriate for the input and context."""
        director = None
        return director

# Result scenarios:
# In (rough) order of increasing complexity:
# * stateless and reproducible locally: calculate when needed
# * stateful and reproducible locally: calculate as needed, but implementation needs to avoid
#   resource contention, race conditions, reentrancy issues.
# * deferred: need to allow resource manager to provide data as it becomes available.
# In the general case, then, the Result handle should
# 1. allow a consumer to register its interest in the result with its own resource manager
#    and allow itself to be provided with the result when it is available.
# 2. Allow the holder of the Result handle to request the data immediately, with the understanding
#    that the surrounding code is blocked on the request.
# Note that in case (1), the holder of the handle may not use the facility, especially if it will
# be using (2).

# TODO: This class can be removed for tidiness when more sophisticated classes are avialable.
# E.g. caching Results, ensemble-safe results.
class ImmediateResult(AbstractResult):
    """Simple Result obtainable with local computation.

    Operation and result are stateless and can be evaluated in any Context.
    """
    def __init__(self, implementation, input):
        """`implementation` is idempotent and may be called repeatedly without (additional) side effects."""
        # Retain input information for introspection.
        self.__input = input
        self.__cached_value = implementation(*input.args, **input.kwargs)
        # TODO: need a utility to resolve the base type of a value that may be a proxy object.
        self._dtype = None

    @property
    def dtype(self):
        return self._dtype

    def result(self):
        return self.__cached_value

class CachingResult(AbstractResult):
    """Simple Result obtainable with local computation.

    Operation and result are stateless and can be evaluated in any Context.
    """

    def __init__(self, implementation, input):
        """`implementation` is idempotent and may be called repeatedly without (additional) side effects."""
        self.__implementation = implementation
        self.__input = input
        self.__cached_value = None
        # TODO: need a utility to resolve the base type of a value that may be a proxy object.
        self._dtype = None

    @property
    def dtype(self):
        return self._dtype

    def result(self):
        # TODO: self.__cached_value should not be a resource owned by this object,
        # but a resource managed by the Context. This object should be just a proxy.
        if self.__cached_value is None:
            self.__cached_value = self.__implementation(self.__input)
        return self.__cached_value


class OutputMeta(type):
    """Output proxy metatype to annotate and enforce interface."""
    def __new__(meta, name, bases, class_dict):
        for key, value in class_dict.items():
            if isinstance(value, NamedResult):
                value.name = key
                value.internal_name = '_' + key
        cls = type.__new__(meta, name, bases, class_dict)
        return cls

class _OutputHandler(object, metaclass=OutputMeta):
    """Generic output attribute implementation."""

class Output(object):
    """Descriptor class for output proxy behavior."""
    # TODO: Output is a special case of Result that is a collection of Results.
    def __init__(self):
        self.name = None
        self.internal_name = None

    def __get__(self, instance, instance_type):
        if instance is None: return self
        return getattr(instance, self.internal_name, '')

    def __set__(self, instance, value):
        setattr(instance, self.internal_name, value)


def computed_result(function):
    """Decorate a function to get a helper that produces an object with Result behavior.
    """
    @functools.wraps(function)
    def new_function(*args, **kwargs):
        """When called, the new function produces an ImmediateResult object.

        The new function has the same signature as the original function, but can accept
        proxy objects (Result objects) for arguments if the provided proxy objects represent
        a type compatible with the original signature.

        The ImmediateResult object will be evaluated in the local Context when its `result()`
        method is called the first time.

        Calls to `result()` return the value that `function` would return when executed in
        the local context with the inputs fully resolved.
        """
        # The signature of the new function will accept abstractions of whatever types it originally accepted.
        # * Create a mapping to the original call signature from `input`
        # * Add handling for typed abstractions in wrapper function.
        # * Process arguments to the wrapper function into `input`
        sig = inspect.signature(function)
        # TODO: what to do with exceptions where this introspection and rebinding won't work?
        # ref: https://docs.python.org/3/library/inspect.html#introspecting-callables-with-the-signature-object

        # TODO: create a serializable data structure for inputs discovered from function introspection.

        # TODO: handle typed abstractions in input arguments
        input_pack = sig.bind(*args, **kwargs)

        result_object = ImmediateResult(function, input_pack)
        return result_object
    return new_function


# TODO: These ought to be proper gmxapi operations that can be deserialized or fetched to satisfy data flow.
@computed_result
def concatenate_lists(sublists=()):
    # TODO: Each sublist or sublist element could be a "future" handle; make sure input provider resolves that.
    # TODO: Returned list should be an NDArray.
    full_list = []
    for sublist in sublists:
        if sublist is not None:
            if isinstance(sublist, (str, bytes)):
                full_list.append(sublist)
            else:
                full_list.extend(sublist)
    return full_list

@computed_result
def make_constant(value):
    """Create a source of the provided value.

    Accepts a value of any type. The object returned has a definite type.
    """
    return type(value)(value)


def function_wrapper(output=()):
    """Generate a decorator for wrapped functions with signature manipulation.

    New function accepts the same arguments, with additional arguments required by
    the API.

    The new function returns an object with an `output` attribute containing the named outputs.

    Example:
        @function_wrapper(output=['spam', 'foo'])
        def myfunc(parameter=None, output=None):
            output.spam = parameter
            output.foo = parameter + ' ' + parameter

        operation1 = myfunc(parameter='spam spam')
        assert operation1.output.spam.result() == 'spam spam'
        assert operation1.output.foo.result() == 'spam spam spam spam'
    """
    # TODO: more flexibility?
    #     If 'output' is provided to the wrapper, a data structure will be passed to
    #     the wrapped functions with the named attributes so that the function can easily
    #     publish multiple named results. Otherwise, the `output` of the generated operation
    #     will just capture the return value of the wrapped function.
    # For now, this behavior is obtained with @computed_result

    # TODO: gmxapi operations need to allow a context-dependent way to generate an implementation with input.
    # This function wrapper reproduces the wrapped function's kwargs, but does not allow chaining a
    # dynamic `input` kwarg and does not dispatch according to a `context` kwarg. We should allow
    # a default implementation and registration of alternate implementations. We don't have to do that
    # with functools.singledispatch, but we could, if we add yet another layer to generate a wrapper
    # that takes the context as the first argument. (`singledispatch` inspects the first argument rather
    # that a named argument)

    output_names = list([name for name in output])
    def decorator(function):
        # TODO: add additional allowed arguments to signature.
        @functools.wraps(function)
        def new_helper(*args, **kwargs):

            # TODO: Should be an implementation detail of the context implementation.
            # The gmxapi Python package provides context implementations with ensemble management.
            # A simple operation should be able to easily get an OutputResource generator and/or
            # provide a module-specific implementation.
            # TODO: Dictionary-like assignment seems like a reasonable feature to add.
            class OutputResource(object):
                __slots__ = output_names

            #
            # class OutputResourceMeta(type):
            #     def __new__(mcs, name, bases, class_dict):
            #         for attr in output_names:
            #             # TODO: Replace None with an appropriate descriptor
            #             class_dict[attr] = None
            #         cls = type.__new__(mcs, name, bases, class_dict)
            #         return cls
            #
            # class OutputResource(object, metaclass=OutputResourceMeta):
            #     pass

            class DataProxyMeta(type):
                def __new__(mcs, name, bases, class_dict):
                    for attr in output_names:
                        # TODO: Replace None with an appropriate (read-only) Result descriptor
                        class_dict[attr] = None
                    cls = type.__new__(mcs, name, bases, class_dict)
                    return cls

            # TODO: Should be an implementation detail of the context implementation.
            # The gmxapi Python package provides context implementations with ensemble management.
            # A simple operation should be able to easily get an OutputResource generator and/or
            # provide a module-specific implementation.
            class OutputDataProxy(object, metaclass=DataProxyMeta):
                # TODO: we want some container behavior, in addition to the attributes...
                pass

            class Resources(object):
                # TODO: should be an implementation detail of the Context, valid only during session lifetime.
                __slots__ = ['output']

            class Operation(object):
                @property
                def _output_names(self):
                    return [name for name in output_names]

                def __init__(self, *args, **kwargs):
                    ##
                    # TODO: generate these as types at class level, not as instances at instance level
                    output_publishing_resource = OutputResource()
                    output_data_proxy = OutputDataProxy()
                    # TODO: output_data_proxy needs to have Result attributes now, not just after run.
                    for accessor in self._output_names:
                        setattr(output_publishing_resource, accessor, None)
                        setattr(output_data_proxy, accessor, None)
                    ##
                    self._resources = Resources()
                    self._resources.output = output_publishing_resource
                    self._output = output_data_proxy
                    self.input_args = tuple(args)
                    self.input_kwargs = {key: value for key, value in kwargs.items()}
                    # TODO: generalize
                    self.dependencies = []

                    # If present, kwargs['input'] is treated as an input "pack" providing _default_ values.
                    if 'input' in self.input_kwargs:
                        provided_input = self.input_kwargs.pop('input')
                        if provided_input is not None:
                            assert not 'input' in self.input_kwargs
                            # Try to determine what 'input' is.
                            # TODO: handling should be related to Context...
                            if hasattr(provided_input, 'run'):
                                self.dependencies.append(provided_input)
                            else:
                                # Assume a parameter pack is provided.
                                for key, value in provided_input.items():
                                    if key not in self.input_kwargs:
                                        self.input_kwargs[key] = value
                            assert not 'input' in self.input_kwargs
                    assert not 'input' in self.input_kwargs

                @property
                def output(self):
                    return self._output

                # TODO: This should be composed with help from the Context implementation.
                def run(self):
                    # TODO: take action only if outputs are not already done.
                    # TODO: make sure this gets run if outputs need to be satisfied for `result()`
                    for dependency in self.dependencies:
                        dependency.run()
                    args = []
                    for arg in self.input_args:
                        # TODO: be more rigorous...
                        if hasattr(arg, 'result'):
                            args.append(arg.result())
                        else:
                            args.append(arg)
                    kwargs = {}
                    for key, value in self.input_kwargs.items():
                        if hasattr(value, 'result'):
                            kwargs[key] = value.result()
                        else:
                            kwargs[key] = value

                    assert not 'input' in kwargs
                    function(*args, output=self._resources.output, **kwargs)
                    # TODO: Add publishing infrastructure to connect output resources to published output.
                    for name in output_names:
                        setattr(self._output, name, getattr(self._resources.output, name))

            operation = Operation(*args, **kwargs)
            return operation

        return new_helper

    return decorator

class Input(object):
    """Descriptor class for input proxy behavior."""
    def __init__(self):
        pass

    def __get__(self, instance, instance_type):
        return {}

class ValidateOperation(type):
    """Confirm that an Operation class implements the gmxapi Operation protocol."""
    def __new__(meta, name, bases, class_dict):
        # Don't validate the abstract Operation class
        if bases != (object,):
            if class_dict['run'] is None:
                raise gmxapi.exceptions.ApiError()
            else:
                class_dict['run'].name = 'run'
                class_dict['run'].internal_name = '_' + 'run'
            if class_dict['output'] is None:
                raise gmxapi.exceptions.ApiError()
        return type.__new__(meta, name, bases, class_dict)

class Operation(object, metaclass=ValidateOperation):
    run = None
    output = None

# Need 'make_operation' and 'make_operation_helper'. 'make_operation' returns something
# that can be called with (label=label, input=input, context=context)
# 'make_operation' returns something that can take additional arguments.
# Implement in terms of properties or descriptors and dispatch kwargs to setters.

def make_operation(implementation, input=None, output=None):
    """Utility to create a gmxapi-compliant Operation in an importable module.

    Args:
        implementation : importable class that will be instantiated to execute the operation.

    Keyword Arguments:
        input : list of parameters to accept as keyword arguments and relay to the constructor.
        output : list of outputs to map from instance properties or methods.

    Returns:
        A callable that produces objects supporting the gmxapi operation protocol.

    Example:
        # mymodule.py
    """
    operation_output = output
    if operation_output is None:
        operation_output = object()

    operation_input = input
    if operation_input is None:
        operation_input = []

    # Bind the signature of `implementation` with declared inputs...

    class NewOperation(Operation):
        # When the operation is `run`, the output is brought up-to-date with
        # available inputs and a reference is returned to a single output state
        run = Runner(implementation, input)
        # The output attribute gives (proxy) access to the output of the work graph node
        output = operation_output

    class NewOperationFactory(object):
        def __call__(self, *args, **kwargs):
            return NewOperationFactory.create(*args, **kwargs)

        @classmethod
        def create(cls, *args, **kwargs):
            # TODO: positional keyword arguments must all be parsed and munged
            return NewOperation()

    return NewOperationFactory()
