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
"""
Provide command line operation.
"""
from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals

__all__ = ['commandline_operation']

from os import devnull
import subprocess

from gmxapi import exceptions
from gmxapi import logging
from gmxapi.operation import function_wrapper, concatenate_lists, make_constant
from gmxapi import util

# Module-level logger
logger = logging.getLogger(__name__)
logger.info('Importing gmxapi._commandline_operation')


#
# Wrap the defined function using a decorator that
#    * strips the `output` argument from the signature
#    * provides `output` to the inner function and
#    * returns the output object when called with the shorter signature.
#
@function_wrapper(output=['file', 'erroroutput', 'returncode'])
def cli(command=None, shell=None, output=None):
    """Execute a command line program in a subprocess.

    Configure an executable in a subprocess. Executes when run in an execution
    Context, as part of a work graph or via gmx.run(). Runs in the current
    working directory.

    Shell processing is not enabled, but can be considered for a future version.
    This means that shell expansions such as environment variables, globbing (`*`),
    and other special symbols (like `~` for home directory) are not available.
    This allows a simpler and more robust implementation, as well as a better
    ability to uniquely identify the effects of a command line operation. If you
    think this disallows important use cases, please let us know.

    Arguments:
         command : a tuple (or list) to be the subprocess arguments, including `executable`
         output : mapping of command line flags to output filename arguments
         shell : unused (provides forward-compatibility)

    Arguments are iteratively added to the command line with standard Python
    iteration, so you should use a tuple or list even if you have only one parameter.
    I.e. If you provide a string with `arguments="asdf"` then it will be passed as
    `... "a" "s" "d" "f"`. To pass a single string argument, `arguments=("asdf")`
    or `arguments=["asdf"]`.

    `input` and `output` should be a dictionary with string keys, where the keys
    name command line "flags" or options.

    Example:
        Execute a command named `exe` that takes a flagged option for file name
        (stored in a local Python variable `my_filename`) and an `origin` flag
        that uses the next three arguments to define a vector.

            >>> my_filename = "somefilename"
            >>> result = cli(('exe', '--origin', 1.0, 2.0, 3.0, '-f', my_filename), shell=False)
            >>> assert hasattr(result, 'file')
            >>> assert hasattr(result, 'erroroutput')
            >>> assert hasattr(result, 'returncode')

    Returns:
        A data structure with attributes for each of the results `file`, `erroroutput`, and `returncode`

    Result object attributes:
        * `file`: the mapping of CLI flags to filename strings resulting from the `output` kwarg
        * `erroroutput`: A string of error output (if any) if the process failed.
        * `returncode`: return code of the subprocess.

    """
    # Note: we could make provisions for stdio filehandles in a future version. E.g.
    # * STDOUT is available if a consuming operation is bound to `output.stdout`.
    # * STDERR is available if a consuming operation is bound to `output.stderr`.
    # * Otherwise, STDOUT and/or STDERR is(are) closed when command is called.
    #
    # Warning:
    #     Commands relying on STDIN cannot be used and is closed when command is called.

    # In the operation implementation, we expect the `shell` parameter to be intercepted by the
    # wrapper and set to False.
    if shell:
        raise exceptions.UsageError("Operation does not support shell processing.")

    if isinstance(command, (str, bytes)):
        command = [command]
    command = list([arg for arg in command])
    try:
        command[0] = util.which(command[0])
    except Exception:
        raise exceptions.ValueError('command argument could not be resolved to an executable file path.')

    # TODO: (FR9) can these be a responsibility of the code providing 'resources'?
    stdin = open(devnull, 'r')
    stdout = open(devnull, 'w')
    stderr = open(devnull, 'w')

    erroroutput = None
    logger.debug('executing subprocess')
    with stdin as null_in, stdout as null_out, stderr as null_err:
        try:
            returncode = subprocess.check_call(command,
                                               shell=shell,
                                               stdin=null_in,
                                               stdout=null_out,
                                               stderr=null_err)
        except subprocess.CalledProcessError as e:
            logger.info("commandline operation had non-zero return status when calling {}".format(e.cmd))
            erroroutput = e.output
            returncode = e.returncode
    # resources.output.erroroutput.publish(erroroutput)
    # resources.output.returncode.publish(returncode)
    # `publish` is descriptive, but redundant. Access to the output data handler is
    # assumed to coincide with publishing, and we assume data is published when the
    # handler is released. A class with a single `publish` method is overly complex
    # since we can just use the assignment operator.
    output.erroroutput = erroroutput
    output.returncode = returncode


# This doesn't need to be a formal operation. It can just generate the necessary data flow.
def filemap_to_flag_list(filemap=None):
    """Convert a map of command line flags and filenames to a list of command line arguments.

    Used to map inputs and outputs of command line tools to and from gmxapi data handles.
    User provides mappings of flags and filenames so that gmxapi can construct an
    executable command line.

    Primary use case is implicit. commandline_operation() generates this operation based on
    user input, and sends the output to cli()

    Arguments:
        filemap : key-value map of command line flags and filename arguments

    Returns:
        list of strings and/or gmxapi data references
    """
    flag_list = []
    if filemap is not None:
        for flag, value in filemap.items():
            flag_list.extend((flag, value))
    # TODO: (FR4) Should be a NDArray(dtype=str)
    return flag_list


# TODO: (FR4) Use generating function or decorator that can validate kwargs?
# TODO: (FR4) Outputs need to be fully formed and typed in the object returned from the helper (decorated function).
# Question: Should the implementing function be able to modify the work graph?
# Or should we have a restriction that implementation functions cannot call operations,
# and have more separation between the definition of the helper function and the definition
# of the implementation, allowing only helper functions to add to the work graph?
# At the moment, we don't have a way to prevent acquisition of new operation handles
# in runtime implementation functions, but we can discourage it for the moment and
# in the future we can check the current Context whenever getting a new operation
# handle to decide if it is allowed. Such a check could only be performed after the work is launched, of course.
def commandline_operation(executable=None, arguments=None, input_files=None, output_files=None, **kwargs):
    """Helper function to execute a subprocess in gmxapi data flow.

    Generate a chain of operations to process the named key word arguments and handle
    input/output data dependencies.

    Arguments:
        executable : name of an executable on the path
        arguments : list of positional arguments to insert at argv[1]
        input_files : mapping of command-line flags to input file names
        output_files : mapping of command-line flags to output file names

    Output:
        The output node of the resulting operation handle contains
        * `file`: the mapping of CLI flags to filename strings resulting from the `output` kwarg
        * `erroroutput`: A string of error output (if any) if the process failed.
        * `returncode`: return code of the subprocess.

    """
    command = concatenate_lists((executable,
                                 arguments,
                                 filemap_to_flag_list(input_files),
                                 filemap_to_flag_list(output_files)))
    shell = make_constant(False)
    cli_args = {'command': command,
                'shell': shell}
    cli_args.update(kwargs)
    return cli(**cli_args)
