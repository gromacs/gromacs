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
# https://github.com/kassonlab/gmxapi/blob/v0.0.7.4/src/gmx/workflow.py
# # https://github.com/kassonlab/gmxapi/blob/v0.0.7.4/LICENSE
"""
Provide workflow-level utilities and classes
============================================

Supports the implementation of operations in the gmxapi.simulation module.
"""

__all__ = ["from_tpr", "WorkSpec", "WorkElement"]

import collections
import warnings
import weakref
from typing import Text, Iterable, Set

import gmxapi as gmx
from gmxapi import exceptions

# Module-level logger
logger = gmx.logger.getChild("simulation.workflow")
logger.info("Importing gmx.workflow")

# Work specification version string.
workspec_version = "gmxapi_workspec_0_1"
logger.info("Using schema version {}.".format(workspec_version))


def to_utf8(input) -> bytes:
    """Return a utf8 encoded byte sequence of the Unicode ``input`` or its string representation.

    Returns:
         :py:bytes byte sequence.
    """
    if isinstance(input, str):
        value = input.encode("utf-8")
    elif isinstance(input, bytes):
        value = input
    else:
        try:
            string = str(input)
            value = string.encode("utf-8")
        except Exception as e:
            raise exceptions.ValueError(
                "Input cannot be interpreted as a UTF-8 compatible string."
            ) from e
    return value


def to_string(input) -> str:
    """Return a Unicode string representation of ``input``.

    If ``input`` or its string representation is not already a Unicode object, attempt to decode as utf-8.

    Returns a native string, decoding utf-8 encoded byte sequences if necessary.
    """
    if isinstance(input, str):
        value = input
    else:
        try:
            value = input.decode("utf-8")
        except Exception:
            try:
                value = str(input)
            except Exception as e:
                raise exceptions.ValueError(
                    "Cannot find a string representation of input."
                ) from e
    return value


class GmxMap(dict):
    """Utility/compatibility class to ensure consistent keys.

    Internally, converts all keys to native str for the current interpreter.
    """

    def keys(self):
        for key in dict.keys(self):
            if not isinstance(key, str):
                raise exceptions.ApiError(
                    "Invalid key type found: {} {}".format(key, type(key))
                )
            yield key

    def __getitem__(self, key):
        return super(GmxMap, self).__getitem__(str(key))

    def __setitem__(self, key, item):
        super(GmxMap, self).__setitem__(str(key), item)

    def __delitem__(self, key):
        super(GmxMap, self).__delitem__(str(key))


class WorkSpec(object):
    """
    Container of workflow elements with data dependency
    information and requirements for execution.

    An element cannot be added to a WorkSpec if it has dependencies that are not
    in the WorkSpec.

    Work is added to the specification by passing a WorkElement object to
    :py:func:`WorkSpec.add_element()`.
    Any dependencies in the WorkElement must already be specified in the target WorkSpec.

    When iterated over, a WorkSpec object returns WorkElement objects.
    WorkElement objects are yielded in a valid order to keep dependencies
    satisfied, but not necessarily the same order in which add_element()
    calls were originally made. In other words, the WorkSpec is a directed
    acyclic dependency graph, and its iterator returns nodes in an arbitrary
    but topologically correct order.

    The string representation of a WorkSpec object is a valid JSON serialized data object.

    The schema for version 0.1 of the specification allows data structures like
    the following.
    ::

        {
            'version': 'gmxapi_workspec_0_1',
            'elements':
            {
                'myinput':
                {
                    'namespace': 'gromacs',
                    'operation': 'load_tpr',
                    'params': {'input': ['tpr_filename1', 'tpr_filename2']}
                },
                'mydata':
                {
                    'namespace': 'gmxapi',
                    'operation': 'open_global_data_with_barrier',
                    'params': ['data_filename']
                },
                'mypotential':
                {
                    'namespace': 'myplugin',
                    'operation': 'create_mdmodule',
                    'params': {...},
                    'depends': ['mydata']
                },
                'mysim':
                {
                    'namespace': 'gmxapi',
                    'operation': 'md',
                    'depends': ['myinput', 'mypotential']
                }
            }
        }

    The first mapping (``version``) is required as shown. The ``elements`` map
    contains uniquely named elements specifying an operation, the operation's
    namespace, and parameters and dependencies of the operation for this element.
    ``depends`` is a sequence of string names of elements that are also in the
    work spec. ``params`` is a key-value map with string keys and values that
    are valid JSON data. ``namespace`` and ``operation`` are strings that the
    :py:class:`Context <gmx.context.Context>` can map to directors it uses to
    construct the session. Namespace ``gmxapi`` is reserved for operations
    specified by the API. Namespace ``gromacs`` is reserved for operations
    implemented as GROMACS adapters (versioned separately from gmxapi). The
    period character (".") has special meaning and should not be used in naming
    elements, namespaces, or operations.

    """

    def __init__(self):
        self.version = workspec_version
        self.elements = GmxMap()
        self.__context_weak_ref = None

    @property
    def _context(self):
        referent = None
        if self.__context_weak_ref is not None:
            referent = self.__context_weak_ref()
        return referent

    @_context.setter
    def _context(self, context):
        # We're moving towards having the context own the work, so the work should
        # not own the context.
        self.__context_weak_ref = weakref.ref(context)

    def __chase_deps(self, source_set: Set[str], name_list: Iterable[Text]):
        """Helper to recursively generate dependencies before dependents.

        Given a set of WorkElement objects and a list of element names, generate WorkElements for
        the members of name_list plus their dependencies in an order such that dependencies are
        guaranteed to occur before their dependent elements.

        For example, to sequence an entire work specification into a reasonable order for instantiation, use

            >>> workspec.__chase_deps(set(workspec.elements.keys()), list(workspec.elements.keys()))

        Note: as a member function of WorkSpec, we have access to the full WorkSpec elements data at all
        times, giving us extra flexibility in implementation and arguments.

        Args:
            source_set: a (super)set of element names from the current work spec (will be consumed)
            name_list: subset of *sources* to be sequenced

        Returns:
            Sequence of WorkElement objects drawn from the names in *source_set*

        Requires that WorkElements named in *name_list* and any elements on which
        they depend are all named in *source_list* and available in the current
        work spec.

        Warning: *source_set* is a reference to an object that is modified arbitrarily.
        The caller should not re-use the object after calling _chase_deps().
        (Make a copy first, if needed.)

        TODO: Separate out DAG topology operations from here and Context.__enter__()
        Our needs are simple enough that we probably don't need an external dependency
        like networkx...
        """
        # Recursively (depth-first) generate a topologically valid serialized DAG from source_set.
        assert isinstance(source_set, set)
        if isinstance(name_list, (str, bytes)):
            warnings.warn(
                "name_list appears to be a single name. Disambiguate a string by passing a list or tuple."
            )
        assert isinstance(name_list, collections.abc.Iterable)

        # Make a copy of name_list in case the input reference is being used elsewhere during
        # iteration, such as for source_set, which is modified during the loop.
        for name in tuple(name_list):
            assert isinstance(name, str)
            if name in source_set:
                source_set.remove(name)
                element = WorkElement.deserialize(
                    self.elements[name], name=name, workspec=self
                )
                dependencies = element.depends
                # items in element.depends are either element names or ensembles of element names.
                for item in dependencies:
                    if isinstance(item, (list, tuple, set)):
                        dependency_list = item
                    else:
                        if not isinstance(item, str):
                            raise exceptions.ValueError(
                                "Dependencies should be a string or sequence of strings. Got {}".format(
                                    type(item)
                                )
                            )
                        dependency_list = [item]
                    for dependency in dependency_list:
                        for recursive_dep in self.__chase_deps(
                            source_set, (dependency,)
                        ):
                            yield recursive_dep
                yield element
            else:
                # Note: The user is responsible for ensuring that source_set is complete.
                # Otherwise, we would need to maintain a list of elements previously yielded.
                pass

    def __iter__(self):
        source_set = set(self.elements.keys())
        for element in self.__chase_deps(source_set, source_set):
            yield element

    def __hash__(self):
        """Uniquely identify this work specification.

        Allows the spec to be used as a dictionary key in Python. Note that this hash is possibly dependent on the
        Python implementation. It is not part of the gmxapi specification and should not be used outside of a single
        invocation of a script.
        """
        # Hash the serialized elements, concatenated as a single string. Note that the order of elements and their
        # contents is not guaranteed, but should be consistent within a script invocation.
        return hash(to_string(self.serialize()))

    def add_element(self, element):
        """Add an element to a work specification if possible.

        Adding an element to a WorkSpec must preserve the validity of the workspec, which involves several checks.
        We do not yet check for element uniqueness beyond a string name.

        If an element is added that was previously in another WorkSpec, it must first be removed from the
        other WorkSpec.
        """
        if (
            hasattr(element, "namespace")
            and hasattr(element, "operation")
            and hasattr(element, "serialize")
        ):
            if (
                not hasattr(element, "name")
                or element.name is None
                or len(str(element.name)) < 1
            ):
                raise exceptions.UsageError(
                    "Only named elements may be added to a WorkSpec."
                )
            if element.name in self.elements:
                raise exceptions.UsageError(
                    "Elements in WorkSpec must be uniquely identifiable."
                )
            if hasattr(element, "depends"):
                for dependency in element.depends:
                    if not dependency in self.elements:
                        raise exceptions.UsageError(
                            "Element dependencies must already be specified before an Element may be added."
                        )
            # Okay, it looks like we have an element we can add
            if (
                hasattr(element, "workspec")
                and element.workspec is not None
                and element.workspec is not self
            ):
                raise exceptions.Error(
                    "Element must be removed from its current WorkSpec to be added to this WorkSpec, but element "
                    "removal is not yet implemented."
                )
            self.elements[element.name] = element.serialize()
            element.workspec = self
        else:
            raise exceptions.ValueError(
                "Provided object does not appear to be compatible with gmx.workflow.WorkElement."
            )
        logger.info("Added element {} to workspec.".format(element.name))

    def serialize(self):
        """Serialize the work specification in a form suitable to pass to any Context implementation.

        Serialization is performed with the JSON data serialization module.

        To simplify unique identification of work specifications, this function will also impose rules for reproducibility.

        1. All key-value maps are sorted alphanumerically by their string keys.
        2. Strings must consist of valid ASCII characters.
        3. Output is a byte sequence of the utf-8 encoded densely formatted JSON document.

        Returns:
            ``unicode`` object in Python 2, ``bytes`` object in Python 3

        Output of serialize() should be explicitly converted to a string before passing to a JSON deserializer.

            >>> my_object = my_workspec.serialize()
            >>> my_data_structure = json.loads(my_object.decode('utf-8'))
            >>> # or...
            >>> my_data_structure = json.loads(my_object, encoding='utf-8')

        """
        import json

        # Build the normalized dictionary
        dict_representation = {"version": self.version, "elements": {}}
        for name, element in [
            (e, json.loads(to_string(self.elements[e])))
            for e in sorted(self.elements.keys())
        ]:
            dict_representation["elements"][str(name)] = element
        serialization = json.dumps(
            dict_representation,
            ensure_ascii=True,
            sort_keys=True,
            separators=(",", ":"),
        )
        return serialization.encode("utf-8")

    @classmethod
    def deserialize(serialized):
        import json

        workspec = WorkSpec()
        dict_representation = json.loads(to_string(serialized))
        ver_in = dict_representation["version"]
        ver_out = workspec.version
        if ver_in != ver_out:
            message = "Expected work spec version {}. Got work spec version {}.".format(
                ver_out, ver_in
            )
            raise exceptions.ValueError(message)
        for element in dict_representation["elements"]:
            workspec.elements[element] = dict_representation["elements"][element]
        return workspec

    def uid(self):
        """Get a unique identifier for this work specification.

        Returns:
            hash value

        Generate a cryptographic hash of this work specification that is guaranteed to match that of another equivalent
        work specification. The returned string is a 64-character hexadecimal encoded SHA-256 hash digest of the
        serialized WorkSpec.

        The definition of equivalence is likely to evolve, but currently means a work spec of the
        same version with the same named elements containing the same operations, dependencies, and parameters, as
        represented in the serialized version of the work specification. Note that this does not include checks on the
        actual contents of input files or anything that does not appear in the work specification directly. Also, the
        hash is lossy, so it is remotely conceivable that two specs could have the same hash. The work specs
        should be compared before making any expensive decisions based on work spec equivalence, such as with hash(workspec).

        Element names probably shouldn't be included in the unique identifying information (so that we can optimize out
        duplicated artifacts), but they are. A future API specification may add unique identification to the elements...
        """
        # Get an alphanumeric string of the checksum of the serialized work spec. SHA-256 should require about 43 characters
        # of base64 to represent, which seems reasonable. We need to replace some of the base64 characters to make them
        # filesystem friendly, though. Hexadecimal may be more friendly, but would require 64 characters.
        import hashlib

        data = to_utf8(self.serialize())
        result = hashlib.sha256(data)
        return result.hexdigest()

    def __str__(self):
        """Generate string representation for str() or print().

        The string output should look like the abstract schema for gmxapi_workspec_1_0, but the exact
        format is unspecified and may change in future versions.

        For consistent JSON output, use WorkSpec.serialize().
        """
        import json

        string = to_string(self.serialize())
        data = json.loads(string)
        reserialized = json.dumps(data, indent=4, sort_keys=True)
        return str(reserialized)

    def __repr__(self):
        """Generate Pythonic representation for repr(workspec)."""
        return "gmx.workflow.WorkSpec()"


# A possible alternative name for WorkElement would be Operator, since there is a one-to-one
# mapping between WorkElements and applications of "operation"s. We need to keep in mind the
# sensible distinction between the WorkElement abstraction and the API objects and DAG nodes.
class WorkElement(object):
    """Encapsulate an element of a work specification."""

    def __init__(self, namespace="gmxapi", operation=None, params=None, depends=()):
        self._namespace = str(to_string(namespace))
        # We can add an operations submodule to validate these. E.g. self.operation = gmx.workflow.operations.normalize(operation)
        if operation is not None:
            self._operation = str(to_string(operation))
        else:
            raise exceptions.UsageError("Invalid argument type for operation.")

        # Note: Nothing currently prevents attribute updates by assignment after adding the element to a workspec,
        # but this protocol will be clarified with https://github.com/kassonlab/gmxapi/issues/92
        if params is None:
            self.params = GmxMap()
        elif isinstance(params, dict):
            self.params = GmxMap({to_string(name): params[name] for name in params})
        else:
            raise exceptions.UsageError(
                "If provided, params must be a dictionary of keyword arguments"
            )
        self.depends = []
        for d in depends:
            if isinstance(d, (list, tuple)):
                self.depends.append([str(name) for name in d])
            else:
                self.depends.append(str(d))

        # The Python class for work elements keeps a strong reference to a WorkSpec object containing its description
        self._name = None
        self._workspec = None

    @property
    def namespace(self):
        assert isinstance(self._namespace, str)
        return self._namespace

    @property
    def operation(self):
        assert isinstance(self._operation, str)
        return self._operation

    @property
    def name(self):
        assert isinstance(self._name, (str, type(None)))
        return self._name

    @name.setter
    def name(self, new_name):
        self._name = str(to_string(new_name))

    @property
    def workspec(self):
        return self._workspec

    @workspec.setter
    def workspec(self, input):
        self._workspec = input

    def add_dependency(self, element):
        """Add another element as a dependency.

        First move the provided element to the same WorkSpec, if not already here.
        Then, add to ``depends`` and update the WorkSpec.
        """

        def check_element(element):
            if element.workspec is None:
                self.workspec.add_element(element)
                assert element.workspec is self.workspec
                assert element.name in self.workspec.elements
            elif element.workspec is not self.workspec:
                raise exceptions.ApiError(
                    "Element will need to be moved to the same workspec."
                )
            return True

        if hasattr(element, "workspec") and hasattr(element, "name"):
            check_element(element)
            self.depends.append(element.name)
        else:
            assert isinstance(element, (list, tuple))
            self.depends.append(
                tuple([item.name for item in element if check_element(item)])
            )

        self.workspec.elements[self.name] = self.serialize()

    def serialize(self):
        """Create a byte sequence representation of the work element.

        The WorkElement class exists just to provide convenient handles in Python. The WorkSpec is not actually a
        container of WorkElement objects.

        Returns:
            Byte sequence of utf-8 encoded JSON document. May need to be decoded if needed as a (Unicode) string.

        """
        import json

        output_dict = {
            "namespace": self.namespace,
            "operation": self.operation,
            "params": self.params,
            "depends": self.depends,
        }
        serialization = json.dumps(output_dict, default=lambda x: x.to_list())
        return to_utf8(serialization)

    @classmethod
    def deserialize(cls, input, name=None, workspec=None):
        """Create a new WorkElement object from a serialized representation.

        Arguments:
            input: a serialized WorkElement
            name: new element name (optional) (deprecated)
            workspec: an existing workspec to attach this element to (optional)

        When subclasses become distinct, this factory function will need to do additional dispatching to create an object of the correct type.
        Alternatively, instead of subclassing, a slightly heavier single class may suffice, or more flexible duck typing might be better.
        """
        import json

        input_string = to_string(input)
        args = json.loads(input_string)
        element = cls(
            namespace=args["namespace"],
            operation=args["operation"],
            params=args["params"],
            depends=args["depends"],
        )
        if name is not None:
            element.name = name
            # This conditional is nested because we can only add named elements to a WorkSpec.
            if workspec is not None:
                element.workspec = workspec
                if element.name not in workspec.elements:
                    workspec.add_element(element)
        return element


class SharedDataElement(WorkElement):
    """Work element with MD-specific extensions.

    The schema may not need to be changed, but the API object may be expected to provide additional functionality.
    """

    def __init__(self, params, name=None):
        """Create a blank SharedDataElement representation.

        It may be appropriate to insist on creating objects of this type via helpers or factories, particularly if
        creation requires additional parameters.
        """
        self.args = params["args"]
        self.kwargs = params["kwargs"]
        super(SharedDataElement, self).__init__(
            namespace="gmxapi",
            operation="global_data",
            params={"args": self.args, "kwargs": self.kwargs},
        )
        self.name = name


def get_source_elements(workspec):
    """Get an iterator of the starting nodes in the work spec.

    Source elements have no dependencies and can be processed immediately. Elements with dependencies
    cannot be processed, instantiated, or added to a work spec until after their dependencies have been.

    Args:
        workspec : an existing work specification to analyze, such as by a Context implementation preparing to schedule work.

    Returns:
        iterator of gmx.workflow.WorkElement objects that may be processed without dependencies.

    This function is provided in the API to allow flexibility in how source elements are determined.
    """
    for name in workspec.elements:
        element_data = workspec.elements[name]
        element = WorkElement.deserialize(element_data)
        if len(element.depends) == 0:
            element.name = name
            element.workspec = workspec
            yield (element)


def from_tpr(input=None, **kwargs):
    """Create a WorkSpec from a (list of) tpr file(s).

    Generates a work specification based on the provided simulation input and returns a handle to the
    MD simulation element of the workflow. Key word arguments can override simulation behavior from
    ``input``.

    If the MD operation discovers artifacts from a previous simulation that was launched from the same input,
    the simulation resumes from the last checkpointed step. If ``append_output`` is set ``False``, existing
    artifacts are kept separate from new output with the standard file naming convention,
    and new output begins from the last checkpointed step, if any.

    Setting ``end_time`` redefines the end point of the simulation trajectory from what was provided in
    ``input``. It is equivalent to changing the number of steps requested in the MDP (or TPR) input, but
    the time is provided as picoseconds instead of a number of time steps.

    .. deprecated:: 0.0.7
        If ``steps=N`` is provided and N is an integer
        greater than or equal to 1, the MD operation advances the trajectory by ``N`` steps, regardless of the number
        of simulation steps specified in ``input`` or ``end_time``. For convenience, setting ``steps=None`` does not override
        ``input``.
        Note that when it is not ``None``, ``steps`` takes precedence over ``end_time`` and ``input``, but can still be
        superceded by a signal, such as if an MD plugin or other code has a simulation completion condition that occurs
        before ``N`` additional steps have run.

    Where key word arguments correspond to ``gmx mdrun`` command line options, the corresponding flags are noted below.

    Keyword Arguments:
        input (str): *Required* string or list of strings giving the filename(s) of simulation input
        append_output (bool): Append output for continuous trajectories if True, truncate existing output data if False. (default True)
        end_time (float): Specify the final time in the simulation trajectory, overriding input read from TPR.
        grid (tuple): Domain decomposition grid divisions (nx, ny, nz). (-dd)
        max_hours (float): Terminate after 0.99 times this many hours if simulation is still running. (-maxh)
        pme_ranks (int): number of separate ranks to be used for PME electrostatics. (-npme)
        threads_per_pme_rank (int): Number of OpenMP threads per PME rank. (-ntomp_pme)
        steps (int): Override input files and run for this many steps. (-nsteps; deprecated)
        threads (int): Total number of threads to start for thread-MPI GROMACS. (-nt)
        threads_per_rank (int): number of OpenMP threads to start per MPI rank. (-ntomp)
        tmpi (int): number of thread-MPI ranks to start. (-ntmpi)

    ..  versionchanged:: 0.1
        *pme_threads_per_rank* renamed to *threads_per_pme_rank*.

    Note that run time simulator arguments can also be provided through `gmxapi.mdrun` *runtime_args*.

    Returns:
        simulation member of a gmx.workflow.WorkSpec object

    Produces a WorkSpec with the following data::

        version: gmxapi_workspec_0_1
        elements:
            tpr_input:
                namespace: gromacs
                operation: load_tpr
                params: {'input': ['tpr_filename1', 'tpr_filename2', ...]}
            md_sim:
                namespace: gmxapi
                operation: md
                depends: ['tpr_input']
                params: {'kw1': arg1, 'kw2': arg2, ...}

    Bugs: version 0.0.6
        * There is not a way to programatically check the current step number on disk.
          See https://github.com/kassonlab/gmxapi/issues/56 and https://github.com/kassonlab/gmxapi/issues/85

    """
    import os
    from gmxapi.utility import config

    usage = "argument to from_tpr() should be a valid filename or list of filenames, followed by optional key word arguments."

    # Normalize to tuple input type.
    if isinstance(input, list) or isinstance(input, tuple):
        tpr_list = tuple([to_string(element) for element in input])
    else:
        try:
            tpr_list = (to_string(input),)
        except:
            raise exceptions.UsageError(usage)

    # Check for valid filenames
    for arg in tpr_list:
        if not (os.path.exists(arg) and os.path.isfile(arg)):
            arg_path = os.path.abspath(arg)
            raise exceptions.UsageError(usage + " Got {}".format(arg_path))

    # Prepare a message about the build configuration.
    # From gmxapi/CMakeLists.txt, this is "library", "tmpi", or None
    mpi_type = config().get("gmx_mpi_type")
    mpi_message = ""
    if not mpi_type:
        mpi_message = "Currently using GROMACS that was built with no MPI."
    elif mpi_type == "tmpi":
        mpi_message = "Currently using GROMACS that was built with internal thread-MPI."
    elif mpi_type == "library":
        mpi_message = (
            "Currently using GROMACS that was built against an external MPI library."
        )
    else:
        warnings.warn(
            f"The gmxapi Python package does not recognize the GROMACS library MPI type: {mpi_type}. "
            "Compatibility checks may be impaired."
        )

    # Note: These are runner parameters, not MD parameters, and should be in the call to gmx.run() instead of here.
    # Reference https://github.com/kassonlab/gmxapi/issues/95
    # Also note that it is not well defined whether or which arguments should be allowed to differ within
    # an ensemble. But to allow parallel trajectory continuation, we have to accept distinct `-cpi` arguments.
    # Note that this violates the traditional gmxapi assumption that all ranks see the same
    # instructions. The "md_sim" element of the gmxapi 0.0.7 workspec ends up being unique
    # to the rank that sees it.
    params = {}
    if "-noappend" in kwargs:
        kwargs["append_output"] = False
        del kwargs["-noappend"]
    for arg_key in kwargs:
        if arg_key == "grid" or arg_key == "dd":
            params["grid"] = tuple(kwargs[arg_key])
        elif arg_key == "pme_ranks" or arg_key == "npme":
            params["pme_ranks"] = int(kwargs[arg_key])
        elif arg_key == "threads" or arg_key == "nt":
            if mpi_type == "tmpi":
                params["threads"] = int(kwargs[arg_key])
            else:
                raise exceptions.ValueError(
                    f"'{arg_key}' argument requires thread-MPI GROMACS. " + mpi_message
                )
        elif arg_key == "tmpi" or arg_key == "ntmpi":
            if mpi_type == "tmpi":
                params["tmpi"] = int(kwargs[arg_key])
            else:
                raise exceptions.ValueError(
                    f"'{arg_key}' argument requires thread-MPI GROMACS. " + mpi_message
                )
        elif arg_key == "threads_per_rank" or arg_key == "ntomp":
            params["threads_per_rank"] = int(kwargs[arg_key])
        elif (
            arg_key == "pme_threads_per_rank"
            or arg_key == "threads_per_pme_rank"
            or arg_key == "ntomp_pme"
        ):
            # TODO: Remove this temporary accommodation.
            assert not gmx.version.api_is_at_least(0, 2)
            if arg_key == "pme_threads_per_rank":
                warnings.warn(
                    "Key word pme_threads_per_rank has been renamed to threads_per_pme_rank.",
                    DeprecationWarning,
                )
            params["threads_per_pme_rank"] = int(kwargs[arg_key])
        elif arg_key == "steps" or arg_key == "nsteps":
            if kwargs[arg_key] is None:
                # None means "don't override the input" which is indicated by a parameter value of -2 in GROMACS 2019
                steps = -2
            else:
                # Otherwise we require steps to be a positive integer
                try:
                    steps = int(kwargs[arg_key])
                    if steps < 1:
                        raise exceptions.ValueError("steps to run must be at least 1")
                except (TypeError, ValueError) as e:
                    # steps is not an integer.
                    raise exceptions.TypeError(
                        '"steps" could not be interpreted as an integer.'
                    )
                # The "nsteps" command line flag will be removed in GROMACS 2020
                # and so "steps" is deprecated in gmxapi 0.0.7
                warnings.warn(
                    "`steps` keyword argument is deprecated. Consider `end_time` instead.",
                    DeprecationWarning,
                )
            params["steps"] = steps
        elif arg_key == "max_hours" or arg_key == "maxh":
            params["max_hours"] = float(kwargs[arg_key])
        elif arg_key == "append_output":
            # Try not to encourage confusion with the `mdrun` `-noappend` flag, which would be a confusing double negative if represented as a bool.
            params["append_output"] = bool(kwargs[arg_key])
        elif arg_key == "end_time":
            params[arg_key] = float(kwargs[arg_key])
        else:
            # Other arguments are accepted without checking or conversion.
            params[arg_key] = kwargs[arg_key]

    # Create an empty WorkSpec
    workspec = WorkSpec()

    # Create and add the Element for the tpr file(s)
    inputelement = WorkElement(
        namespace="gromacs", operation="load_tpr", params={"input": tpr_list}
    )
    inputelement.name = "tpr_input"
    if inputelement.name not in workspec.elements:
        # Operations such as this need to be replaced with accessors or properties that can check the validity of the WorkSpec
        workspec.elements[inputelement.name] = inputelement.serialize()
        inputelement.workspec = workspec

    # Create and add the simulation element
    # We can add smarter handling of the `depends` argument, but it is only critical to check when adding the element
    # to a WorkSpec.
    mdelement = WorkElement(operation="md", depends=[inputelement.name], params=params)
    mdelement.name = "md_sim"
    # Check that the element has not already been added, but that its dependency has.
    workspec.add_element(mdelement)

    return mdelement
