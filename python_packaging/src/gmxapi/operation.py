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
"""

import functools
import inspect

__all__ = ['computed_result',
           'append_list',
           'concatenate_lists'
           ]

from gmxapi import exceptions


class ImmediateResult(object):
    """Data handle for a simple result.

    Instances of this class can be used to provide a gmxapi compatible data
    handle for trivial operations. Operation and result are stateless and can be
    evaluated in any Context.

    Used internally to implement the computed_result factory. The interface for
    this class will evolve as the gmxapi data model evolves. Generally, code
    providing gmxapi data sources should use one of the factories or decorators
    provided in the gmxapi.operation module rather than instantiating from this
    class directly.
    """

    def __init__(self, implementation=None, input=None):
        """Wrap a callable for a simple data source that does not need Future behavior.

        Provides a gmxapi compatible interface for data sources.

        Arguments:
            implementation : Python callable that consumes ``input`` and returns data
            input : object compatible with the call signature of ``implementation``

        ``input`` must have an ``args`` attribute and a ``kwargs`` attribute to be used as

            implementation(*input.args, **input.kwargs)

        Callers should not assume when or how often ``implementation`` could be called.
        Only suitable for function objects without side effects.
        """
        assert callable(implementation)
        assert hasattr(input, 'args')
        assert hasattr(input, 'kwargs')
        # Retain input information for introspection.
        self.__input = input

        self.__cached_value = implementation(*input.args, **input.kwargs)
        # TODO: (FR4) need a utility to resolve the base type of a value
        #  that may be a proxy object.
        self._dtype = type(self.__cached_value)

    @property
    def dtype(self):
        """The data type of the return value for the wrapped function."""
        return self._dtype

    def result(self):
        """Return value of the wrapped function."""
        return self.__cached_value


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

    @functools.wraps(function)
    def new_function(*args, **kwargs):
        # The signature of the new function will accept abstractions
        # of whatever types it originally accepted. This wrapper must
        # * Create a mapping to the original call signature from `input`
        # * Add handling for typed abstractions in wrapper function.
        # * Process arguments to the wrapper function into `input`

        sig = inspect.signature(function)
        # Note: Introspection could fail.
        # TODO: Figure out what to do with exceptions where this introspection
        #  and rebinding won't work.
        # ref: https://docs.python.org/3/library/inspect.html#introspecting-callables-with-the-signature-object

        # TODO: (FR3+) create a serializable data structure for inputs discovered
        #  from function introspection.

        # TODO: (FR4) handle typed abstractions in input arguments

        input_list = []
        for arg in args:
            if hasattr(arg, 'result'):
                input_list.append(arg.result())
            else:
                input_list.append(arg)
        input_dict = {}
        for name, value in kwargs.items():
            if hasattr(value, 'result'):
                input_dict[name] = value.result()
            else:
                input_dict[name] = value

        input_pack = sig.bind(*input_list, **input_dict)

        result_object = ImmediateResult(function, input_pack)
        return result_object

    return new_function


@computed_result
def append_list(a: list = (), b: list = ()):
    """Operation that consumes two lists and produces a concatenated single list."""
    # TODO: (FR3) Each sublist or sublist element could be a "future" handle;
    #  make sure input provider resolves that.
    # TODO: (FR4) Returned list should be an NDArray.
    if isinstance(a, (str, bytes)) or isinstance(b, (str, bytes)):
        raise exceptions.ValueError('Input must be a pair of lists.')
    return list(a) + list(b)


def concatenate_lists(sublists: list = ()):
    """Combine data sources into a single list.

    A trivial data flow restructuring operation
    """
    if isinstance(sublists, (str, bytes)):
        raise exceptions.ValueError('Input must be a list of lists.')
    if len(sublists) == 1:
        return make_constant(sublists[0])
    if len(sublists) == 2:
        return append_list(sublists[0], sublists[1])
    else:
        return append_list(sublists[0], concatenate_lists(sublists[1:]))


@computed_result
def make_constant(value):
    """Provide a predetermined value at run time.

    This is a trivial operation that provides a (typed) value, primarily for
    internally use to manage gmxapi data flow.

    Accepts a value of any type. The object returned has a definite type and
    provides same interface as other gmxapi outputs. Additional constraints or
    guarantees on data type may appear in future versions.
    """
    # TODO: (FR4+) Manage type compatibility with gmxapi data interfaces.
    return type(value)(value)
