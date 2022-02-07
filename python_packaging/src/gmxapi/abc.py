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

"""Abstract base classes for gmxapi Python interfaces.

This module consolidates definitions of some basic interfaces in the gmxapi
Python package. These definitions are evolving and mostly for internal use, but
can be used to check compatibility with the gmxapi implementation details that
are not otherwise fully specified by the API.

Refer to [PEP 484](https://www.python.org/dev/peps/pep-0484) and to
[PEP 526](https://www.python.org/dev/peps/pep-0526) for Python background.

Developer Note:
    Type checking fails when accessing attributes of generic classes (see PEP 484).
    This affects some scoping choices.

Developer Note:
    It is worth noting some details regarding type erasure with Python generics.
    Generic type parameters (class subscripts) allow for annotation of dynamic
    types that are bound when the generic is instantiated. The subscript is not
    part of the run time class. This means that we can use subscripted generic
    classes for static type checking, but use for any other purpose is discouraged.
    Thus, type variable parameters to generics have less meaning that do template
    parameters in C++, and are orthogonal to subclassing. Note, though, that a
    subclass is not generic if the generic parameters are bound (explicitly
    with type subscripts or implicitly by omission (implying `typing.Any`).

    Also, type parameters
    (or types of parameters) to functions can be used to dispatch a generic
    function to a specific function.

    In other words: keep in mind the
    orthogonality of generic classes and base classes, and recognize that
    composed objects mimicking C++ template specializations are not distinguishable
    at the class level.

Note: This module overly specifies the API. As we figure out the relationships.
    As we clarify interactions, we should trim these specifications or migrate them
    to gmxapi.typing and focus on developing to functions and function interfaces
    rather than classes or types. In the mean time, these types help to illustrate
    the entities that can exist in the implementation. Just keep in mind that
    these ABCs describe the existence of "protocols" more than the details of
    the protocols, which may take longer to define.

..  todo:: Clarify protocol checks, ABCs, and mix-ins.

"""
# Note that the Python typing module defines generic classes in terms of abstract
# base classes defined in other modules (namely `collections`), but without
# actual inheritance. The ABCs are not intended to be instantiated, and the
# generics are very mangled objects that cannot be instantiated. However,
# user-defined subclasses of either may be instantiated.
#
# This restriction is probably due somewhat to implementation constraints in the
# typing module, but it represents a Separation of Concerns that we should
# consider borrowing in our model. Practically, this can mean using an abstract
# for run time checking and a generic for static type checking.
#
# There may not be a compelling reason to
# rigorously separate ABC and generic type, or to disallow instantiating a
# generic unless it is also an abstract. Note the distinction, though, between
# abstract generics and fully implemented generics. NDArray
# and Future are likely to be examples of fully implemented generics, while
# Context and various interface types are likely to have abstract generics.
#
# Use abstract base classes to define interfaces and relationships between
# interfaces. Use `typing` module machinery for static type checking. Use the
# presence or absence of an expected interface, or exception handling, for run
# time type checking. Use `isinstance` and `issubclass` checking against
# non-generic abstract base classes when the valid interface is complicated or
# unknown in the caller's context.
#
# Also, note that abstract base classes cannot assert that a simple data member
# is provided by a subclass, but class data members of Generic classes are
# interpreted by the type checker as instance members (unless annotated ClassVar).

import typing
from abc import ABC, abstractmethod
import collections.abc
from typing import Type, Callable


class NDArray(collections.abc.Sequence, ABC):
    """N-dimensional data interface.


    """
    # TODO: Fix the data model and ABC vs. type-checking conflicts so that we can
    #       recognize NDArrays and NDArray futures as compatible data in data flow operations.
    # We cannot do the following because of the forward reference to Future:
    # @classmethod
    # def __subclasshook__(cls, subclass):
    #     """Determine whether gmxapi should consider the provided type to be consistent with an NDArray."""
    #     if subclass is cls:
    #         return True
    #     if cls is NDArray:
    #         # For the purposes of gmxapi data flow, a Future[NDArray] is equivalent to an NDArray.
    #         if Future in subclass.__mro__:
    #             return issubclass(subclass.dtype, cls)
    #         else:
    #             def is_compatible(candidate):
    #                 if issubclass(candidate, collections.abc.Sequence) \
    #                         and not issubclass(candidate, (str, bytes)):
    #                     return True
    #                 return False
    #             any(is_compatible(base) for base in subclass.__mro__)
    #             return True
    #     return NotImplemented
NDArray.register(list)


# TODO: Define an enumeration.
SourceProtocol = typing.NewType('SourceProtocol', str)


class Resource(ABC):
    """gmxapi resource object interface.

    Resources are generally provided by Operation instances, but may be provided
    by the framework, or abstractly as the outputs of other scripts.

    A resource is owned in a specific Context. A resource handle may be created
    in a different Context than that in which the resource is owned. The Context
    instances negotiate the resource availability in the new Context.
    """
    # No public interface is yet defined for the Resource API.

    # Resource instances should have an attribute serving as a sequence of available
    # source protocols in decreasing order of preference.
    # TODO: Clarify. Define an enumeration with the allowed values.
    # TODO: Enforce. The existence of this attribute cannot be confirmed by the abc.ABC machinery.
    # TODO: Consider leaving this off of the ABC specification and make it an implementation
    #  detail of a single_dispatch function.
    # TODO: Similarly, find a way for the subscriber to have a chain of handlers.
    _gmxapi_source_protocol = typing.Sequence[SourceProtocol]


class Future(Resource):
    """Data source that may represent Operation output that does not yet exist.

    Futures represent "immutable resources," or fixed points in the data flow.
    """
    @property
    @abstractmethod
    def dtype(self) -> type:
        ...

    @abstractmethod
    def result(self) -> typing.Any:
        ...

    @classmethod
    def __subclasshook__(cls, C):
        if cls is Future:
            if any("result" in B.__dict__ and callable(B.result) for B in C.__mro__):
                return True
        return NotImplemented


class MutableResourceSubscriber(ABC):
    """Required interface for subscriber of a MutableResource.

    A MutableResource is bound to a collaborating object by passing a valid
    Subscriber to the resource's *subscribe()* method.
    """


class MutableResource(Resource):
    """An Operation interface that does not represent a fixed point in a data flow.

    Providers and consumers of mutable resources are more tightly coupled than
    data edge terminals and have additional machinery for binding at run time.
    Examples include the simulation plugin interface, binary payloads outside of
    the standard gmxapi types, and any operation interaction that the current
    context must defer to lower-level details.

    There is not yet a normative interface for MutableResources, but a consumer
    of MutableResources has chance to bind directly to the provider of a
    MutableResource without the mediation of a DataEdge. Accordingly, providers
    and consumers of MutableResources must be able to be instantiated in the
    same Context.
    """
    def subscribe(self, subscriber: MutableResourceSubscriber):
        """Create a dependency on this resource.

        Allows a gmxapi compatible operation to bind to this resource.
        The subscribing object will be provided with a lower-level registration
        interface at run time, as computing elements are being initialized.
        The nature of this registration interface is particular to the type of
        resource and its participants. See, for instance, the MD plugin binding
        interface.
        """


class OutputDataProxy(ABC):
    """A collection of Operation outputs.

    This abstract base class describes the interface to the output of an
    operation / work node.
    """
    # TODO: Specification.
    # Currently, the common aspect of OutputDataProxy is that a class has public
    # attributes that are exclusively OutputDescriptor instances, meaning that
    # getattr(instance, attr) returns a Future object. However, there are several
    # ways to implement getattr, and this sort of check in an ABC does not appear
    # to be common.
    #
    # We might choose to assert that all public attributes must be compatible data
    # descriptors, in conjunction with defining a more specific OutputDataProxy
    # metaclass, but this makes for a dubiously growing chain of data descriptors
    # we use for the output access.
    #
    # The data model might be cleaner if we move to something more
    # like a Collection or Mapping with more conventional getters, but we would
    # lose the ability to check type on individual elements. (However, the
    # typing of return values is not normally the defining aspect of an ABC.)
    #
    # Another alternative is to collapse the contents of the `output` attribute
    # into the Operation handle type, strongly define all handle types (so that
    # the type checker can identify the presence of attributes), and rely only
    # on type checking at the level of the data descriptors. (Dynamically defined
    # OutputDataProxy classes are the execption, rather than the rule.)
    #
    # We will need to consider the details of type checkers and syntax inspection
    # tools, like Jedi and MyPy, to make design choices that maximize API usability
    # and discoverability.


class OperationReference(ABC):
    """Client interface to an element of computational work already configured.

    An "instance" of an operation is assumed to be a node in a computational
    work graph, owned and managed by a Context. This class describes the
    interface of the reference held by a client once the node exists.

    The convergence of OperationReferences with Nodes as the results of the action
    of a Director implies that a Python user should also be able to "subscribe"
    to an operation handle (or its member resources). This could be a handy feature
    with which a user could register a call-back. Note that we will want to provide
    an optional way for the call-back (as with any subscriber) to assert a chain
    of prioritized Contexts to find the optimal path of subscription.
    """

    @abstractmethod
    def run(self):
        """Assert execution of an operation.

        After calling run(), the operation results are guaranteed to be available
        in the local context.
        """

    @property
    @abstractmethod
    def output(self) -> OutputDataProxy:
        """Get a proxy collection to the output of the operation.

        Developer note: The 'output' property exists to isolate the namespace of
        output data from other operation handle attributes and we should consider
        whether it is actually necessary or helpful. To facilitate its possible
        future removal, do not enrich its interface beyond that of a collection
        of OutputDescriptor attributes. The OutputDataProxy also serves as a Mapping,
        with keys matching the attributes. We may choose to keep only this aspect
        of the interface instead of trying to keep track of the set of attributes.
        """
        ...


class Fingerprint(ABC):
    """Unique global identifier for an Operation node.

    Represents the operation and operation inputs.

    No public interface.
    """


class OutputDescription(ABC):
    """

    There may not be a single OutputDescription base class, since the requirements
    are related to the Context implementation.
    """


class Edge(ABC):
    """Reference to the state and description of a data flow edge.

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

    A concrete Edge is a likely related to the consuming Context, and a single
    abstract base class may not be possible or appropriate. Possible Context-agnostic
    use cases for a global abstract Edge (along with Node) include topological aspects of data
    flow graphs or user-friendly inspection.
    """


class Node(ABC):
    """Object oriented interface to nodes configured in a Context.

    In gmxapi.operation Contexts, this functionality is implemented by subclasses
    of SourceResource.

    Likely additional interfaces for Node include subscribers(), label(), and
    (weak) reference helpers like context() and identifier().

    .. todo:: Converge. A Node is to a concrete Context what an operation handle is to the None (Python UI) Context.
    """
    @abstractmethod
    def handle(self, context: 'Context') -> OperationReference:
        """Get a reference to the Operation in the indicated Context.

        This is equivalent to the reference obtained from the helper function
        or factory that added the node if and only if the Contexts are the same.
        Otherwise, a new node is created in *context* that subscribes to the
        original node.
        """
        # Note that a member function like this is the same as dispatching a
        # director that translates from the Node's Context to *context*

    @abstractmethod
    def output_description(self, context: 'Context') -> OutputDescription:
        """Get a description of the output available from this node.

        Returns a subset of the information available through handle(), but
        without creating a subscription relationship. Allows data sources and
        consumers to determine compatibility and requirements for connecting
        two nodes.
        """

    @abstractmethod
    def input(self) -> Edge:
        """Describe the bound data sources.

        The returned object represents the data edge in the Context managing
        the node, though the data sources may be from other Contexts.
        """

    @abstractmethod
    def fingerprint(self) -> Fingerprint:
        """Uniquely identify this Node.

        Used internally to manage resources, check-point recovery, and messaging
        between Contexts. The fingerprint is dependent on the operation and the
        operation inputs, and is independent of the Context.

        Returns:
            Opaque identifier describing the unique output of this node.
        """

    @abstractmethod
    def operation(self) -> 'OperationImplementation':
        """Get a reference to the registered operation that produces nodes like this.
        """
        # Note that the uniqueness of a node is such that Node.operation() and
        # Node.input() used to get an OperationReference in the same Context
        # should result in a handle to the same Node.


class NodeBuilder(ABC):
    """Add an element of computational work to be managed by a gmxapi Context.

    A Node generally represents an instance of a registered Operation, but the
    only real requirement is that it contains sufficient information to run an
    operation and to direct the instantiation of an equivalent node in a
    different consumer Context.

    In the near future, Node and NodeBuilder will participate in check-pointing
    and in a serialization/deserialization scheme.

    .. todo:: As the NodeBuilder interface is minimized, we can look for a normative
              way to initialize a Generic NodeBuilder that supports the sorts of
              type inference and hinting we would like.
    """

    @abstractmethod
    def build(self) -> OperationReference:
        """Finalize the creation of the operation instance and get a reference."""
        ...

    @abstractmethod
    def set_input_description(self, input_description):
        """Add the details related to the operation input.

        Example: In gmxapi.operation, includes signature() and make_uid()

        .. todo:: This can probably be moved to an aspect of the resource factory.
        """
        ...

    @abstractmethod
    def set_output_factory(self, output_factory):
        """Set the factory that gives output description and resources for the Node.

        Output is not fully describable until the input is known and the Node is
        ready to be instantiated. This is the resource that can be used by the
        Context to finish completely describing the Node. The interface of the
        factory and of any object it produces is a lower level detail of the
        Context and Operation implementations in that Context.

        .. todo:: This can probably be moved to an aspect of the resource factory.
        """

    @abstractmethod
    def add_input(self, name: str, source):
        """Attach a client-provided data source to the named input.

        .. todo:: Generalize to add_resource().
        """
        ...

    @abstractmethod
    def set_handle(self, handle_builder):
        """Set the factory that gives a builder for a handle to the operation.

        .. todo:: Stabilize interface to handle_builder and move to an aspect of the
                  operation registrant.
        """

    @abstractmethod
    def set_runner_director(self, runner_builder):
        """Set the factory that gives a builder for the run-time callable.

        .. todo:: This should be a specialized Director obtained from the registrant
                  by the Context when translating a Node for execution.
        """

    @abstractmethod
    def set_resource_factory(self, factory: Callable):
        """Register a resource factory for the operation run-time resources.

        The factory will be called within the Context

        .. todo:: Along with the merged set_runner/Director, the resource factory
                  is the other core aspect of an operation implementation registrant
                  that the Context should fetch rather than receiving through the
                  NodeBuilder.
        """
        # The factory function takes input in the form the Context will provide it
        # and produces a resource object that will be passed to the callable that
        # implements the operation.
        assert callable(factory)
        ...


class Context(ABC):
    """API Context.

    A Context instance manages the details of the computing environment and
    provides for the allocation of resources.
    All gmxapi data and operations are owned by a Context instance.
    The Context manages the details of how work is run and how data is managed.

    Additionally, a concrete Context implementation determines some details of
    the interfaces used to manage operation execution and data flow. Thus, API
    calls may depend on multiple Contexts when, for instance, there is a source
    Context, a caller Context, and/or a destination Context. For Python data
    types and external interfaces to the gmxapi package (such as public function
    signatures) is equal to *None*.

    This abstract base class (ABC) defines the required interface of a Context
    implementation. Client code should use this ABC for type hints. Concrete
    implementations may, *but are not required*, to subclass from this ABC to
    help enforce compatibility.
    """
    @abstractmethod
    def node_builder(self, *, operation, label=None) -> NodeBuilder:
        """Get a builder for a new work graph node.

        Nodes are elements of computational work, with resources and execution
        managed by the Context. The Context handles parallelism resources, data
        placement, work scheduling, and data flow / execution dependencies.

        This method is used by Operation director code and helper functions to
        add work to the graph.
        """
        ...

    @abstractmethod
    def node(self, node_id) -> Node:
        """Get the indicated node from the Context.

        node_id may be an opaque identifier or a label used when the node was
        added.
        """


class ResourceFactory(ABC):
    """Packager for run time resources for a particular Operation in a particular Context.

    TODO: refine interface.
    """
    @abstractmethod
    def input_description(self, context: Context):
        """Get an input description in a form usable by the indicated Context."""


class OperationDirector(ABC):
    """Interface for Operation entry points.

    An operation director is instantiated for a specific operation and context
    (by a dispatching factory) to update the work managed by the context
    (add a computational element).
    """
    # TODO: How to handle subscriptions? Should the subscription itself be represented
    #  as a type of resource that is passed to __call__, or should the director be
    #  able to subscribe to resources as an alternative to passing with __call__?
    #  Probably the subscription is represented by a Future passed to __call__.

    # TODO: Annotate `resources`, whose validity is determined by both context and operation.
    @abstractmethod
    def __call__(self, resources, label: typing.Optional[str]):
        """Add an element of work (node) and return a handle to the client.

        Implements the client behavior in terms of the NodeBuilder interface
        for a NodeBuilder in the target Context. Return a handle to the resulting
        operation instance (node) that may be specialized to provide additional
        interface particular to the operation.
        """
        ...

    @abstractmethod
    def handle_type(self, context: Context) -> Type[OperationReference]:
        """Get the class used for operation references in this Context.

        Convenience function. May not be needed.
        """
        ...

    @abstractmethod
    def resource_factory(self,
                         source: typing.Union[Context, None],
                         target: typing.Optional[Context] = None) \
            -> typing.Union[ResourceFactory, typing.Callable]:
        """Get an appropriate resource factory.

        The ResourceFactory converts resources (in the form produced by the *source* Context)
        to the form consumed by the operation in the *target* Context.

        A *source* of None indicates that the source is an arbitrary Python function
        signature, or to try to detect appropriate dispatching. A *target* of
        None indicates that the Context of the Director instance is the target.

        As we merge the interface for a NodeBuilder and a RunnerBuilder, this will not need to be
        specified in multiple places, but it is not yet clear where the responsibility lies.
        Ultimately, the operation implementation should only need to provide a single-dispatching
        helper for finalizing input preparation in a target Context when building a node.
        Converting resources from one Context to another during edge creation is a
        collaboration between the resource type, the providing Context, and the consuming
        Context. It should be mediated by functions registered with the Context when the
        Operation is registered, and such functionality should be decoupled from the
        OperationDirector.

        Generally, the client should not need to call the resource_factory directly.
        The director should dispatch an appropriate factory. C++ versions need
        the resource_factory dispatcher to be available for reasons of compilation
        dependencies. Clients may want to use a specific resource_factory to explicitly
        control when and where data flow is resolved. Also, the resource_factory
        for the None Context can be used to set the function signature of the Python
        package helper function. As such, it may be appropriate to make it a "static"
        member function in Python.
        """
        ...


class OperationImplementation(ABC):
    """Essential interface of an Operation implementation.

    Describe the essential features of an Operation that can be registered with
    gmxapi to support building and executing work graphs in gmxapi compatible
    execution contexts.

    An Operation is usable in gmxapi when an OperationImplementation is registered
    with a valid identifier, consisting of a *namespace* and a *name*.

    Generally, the *namespace* is the module implementing the Operation
    and the *name* is a factory or helper importable from the same module that
    triggers the OperationDirector to configure a new instance of the Operation.

    The registered object must be able to describe its namespace and name, and
    to dispatch an OperationDirector appropriate for a given Context.
    """
    # TODO: Either OperationImplementations should be composed or subclasses should
    #  each be singletons. We can still consider that OperationReferences are instances
    #  of OperationImplementations or its subclasses, or that OperationReference classes
    #  have a class data member pointing to a single OperationImplementation instance
    #  or operation_registry value.

    # TODO: Consider a data descriptor and metaclass to validate the name and namespace.
    @classmethod
    @abstractmethod
    def name(self) -> str:
        """The name of the operation.

        Generally, this corresponds to a callable attribute of a Python module
        (named by namespace()) that acts as a factory for new operation instances.
        It is also used by Context implementations to locate code supporting
        the operation.
        """
        ...

    # TODO: Consider a data descriptor and metaclass to validate the name and namespace.
    @classmethod
    @abstractmethod
    def namespace(self) -> str:
        """The namespace of the operation.

        Generally, the namespace corresponds to a Python module importable in
        the execution environment.
        """

    # TODO: Allow this to be an instance method, and register instances.
    # Consider not storing an actual OperationImplementation in the registry.
    # Note, though, that if we want to automatically register on import via
    # base class (meta-class), the functionality must be in the class definition.
    @classmethod
    @abstractmethod
    def director(cls, context: Context) -> OperationDirector:
        """Factory to get an OperationDirector appropriate for the context."""
        ...
