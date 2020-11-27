#
# This file is part of the GROMACS molecular simulation package.
#
# Copyright (c) 2019,2020,2021, by the GROMACS development team, led by
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

__all__ = ['commandline_operation']

import functools
import os
import pathlib
import shutil
import subprocess

import gmxapi as gmx
from gmxapi import exceptions
from gmxapi import logger as root_logger
from gmxapi.datamodel import NDArray
from gmxapi.operation import OutputCollectionDescription

# Module-level logger
logger = root_logger.getChild('commandline')
logger.info('Importing {}'.format(__name__))


@functools.lru_cache()
def _config() -> dict:
    """Get the GROMACS configuration detected during installation.

    If this appears to be a useful function, it may become part of the regular
    interface, but it is currently unadvertised.
    """
    import json
    from importlib.resources import open_text
    with open_text('gmxapi', 'gmxconfig.json') as textfile:
        config = json.load(textfile)
    return config


@functools.lru_cache()
def cli_executable() -> pathlib.Path:
    """Report the installed GROMACS command line executable."""
    path = _config().get('gmx_executable', None)
    if path is not None:
        path = pathlib.Path(os.path.abspath(path))
        if path.is_file():
            return path
    raise exceptions.FeatureNotAvailableError('GROMACS installation unavailable.')


@functools.lru_cache()
def cli_bindir() -> pathlib.Path:
    """Report the installed GROMACS binary directory."""
    path = _config().get('gmx_bindir', None)
    if path is not None:
        path = pathlib.Path(os.path.abspath(path))
        if path.is_dir():
            return path
    raise exceptions.FeatureNotAvailableError('GROMACS installation unavailable.')


# Create an Operation that consumes a list and a boolean to produce a string and an integer.
#
# Wrap the defined function using a decorator that
#    * strips the `output` parameter from the signature
#    * provides `output` publishing proxy to the inner function and
#    * produce a result with attributes for
#       * file: mapping of output flags to output filenames
#       * stdout: process STDOUT
#       * stderr: porcess STDERR
#       * returncode: integer return code of wrapped command
#
# Note that the existence of the 'file' output map is expressed here, but
# the keys of the map are not implicit or set by the wrapped function.
# For the map to be non-empty, it must be defined before the resulting helper
# function is called.
#
# TODO: Operation returns the output object when called with the shorter signature.
#
@gmx.function_wrapper(output={'stdout': str,
                              'stderr': str,
                              'returncode': int})
def cli(command: NDArray, shell: bool, output: OutputCollectionDescription, stdin: str = ''):
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
         command: a tuple (or list) to be the subprocess arguments, including `executable`
         output: mapping of command line flags to output filename arguments
         shell: unused (provides forward-compatibility)
         stdin (str): String input to send to STDIN (terminal input) of the executable.

    Multi-line text sent to *stdin* should be joined into a single string
    (e.g. ``'\n'.join(list_of_strings) + '\n'``).
    If multiple strings are provided to *stdin*, gmxapi will assume an ensemble,
    and will run one operation for each provided string.

    Only string input (:py:func:str) to *stdin* is currently supported.
    If you have a use case that requires streaming input or binary input,
    please open an issue or contact the author(s).

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
            >>> assert hasattr(result, 'stdout')
            >>> assert hasattr(result, 'stderr')
            >>> assert hasattr(result, 'returncode')

    Returns:
        A data structure with attributes for each of the results `file`, `stdout`, `stderr`, and `returncode`

    Result object attributes:
        * `file`: the mapping of CLI flags to filename strings resulting from the `output` kwarg
        * `stdout`: A string mapping from process STDOUT.
        * `stderr`: A string mapping from process STDERR; it will be the
                    error output (if any) if the process failed.
        * `returncode`: return code of the subprocess.

    """
    # In the operation implementation, we expect the `shell` parameter to be intercepted by the
    # wrapper and set to False.
    if shell:
        raise exceptions.UsageError("Operation does not support shell processing.")

    if stdin == '':
        stdin = None

    if isinstance(command, (str, bytes)):
        command = [command]
    command = list([arg for arg in command])

    executable = shutil.which(command[0])
    if executable is None:
        executable = shutil.which(command[0], path=str(cli_bindir()))
    if executable is None:
        raise exceptions.ValueError('"{}" is not found or not executable.'.format(command[0]))
    command[0] = executable

    # TODO: (FR9) Can OS input/output filehandles be a responsibility of
    #  the code providing 'resources'?

    stdout = ''
    stderr = ''
    logger.debug('executing subprocess')
    try:
        completed_process = subprocess.run(command,
                                           shell=shell,
                                           input=stdin,
                                           check=True,
                                           stdout=subprocess.PIPE,
                                           stderr=subprocess.PIPE,
                                           encoding='utf-8',
                                           universal_newlines=True
                                           )
        returncode = completed_process.returncode
        # TODO: Resource management code should manage a safe data object for `output`.
        logger.debug('STDOUT:')
        if completed_process.stderr is not None:
            for line in completed_process.stdout.split('\n'):
                logger.debug(line)
        else:
            logger.debug('STDOUT is empty')
        logger.debug('STDERR:')
        if completed_process.stderr is not None:
            for line in completed_process.stderr.split('\n'):
                logger.debug(line)
        else:
            logger.debug('STDERR is empty')

        stdout = completed_process.stdout
        stderr = completed_process.stderr

    except subprocess.CalledProcessError as e:
        logger.info("commandline operation had non-zero return status"
                    "when calling {}".format(e.cmd))
        stdout = e.stdout
        stderr = e.stderr
        returncode = e.returncode

    # Publish outputs.
    output.stdout = stdout
    output.stderr = stderr
    output.returncode = returncode


# TODO: (FR4) Make this a formal operation to properly handle gmxapi data dependencies.
#  The consumer of this operation has an NDArray input. filemap may contain gmxapi data flow
#  aspects that we want the framework to handle for us.
def filemap_to_flag_list(filemap: dict = None):
    """Convert a map of command line flags and filenames to a list of command line arguments.

    Used to map inputs and outputs of command line tools to and from gmxapi data handles.
    User provides mappings of flags and filenames so that gmxapi can construct an
    executable command line.

    Primary use case is implicit. commandline_operation() instantiates this operation based on
    user input, and sends the output to cli()

    Arguments:
        filemap: key-value map of command line flags and filename arguments

    Returns:
        list of strings and/or gmxapi data references
    """
    result = []
    if filemap is not None:
        for key, value in filemap.items():
            # Note that the value may be a string, a list, an ndarray, or a future
            if not isinstance(value, (list, tuple, NDArray)):
                if hasattr(value, 'result') and value.dtype == NDArray:
                    pass
                elif hasattr(value, 'result') and value.dtype != NDArray:
                    # TODO: Fix this ugly hack when we have proper Future slicing and can make NDArray futures.
                    result_function = value.result
                    value.result = lambda function=result_function: [function()]
                else:
                    value = [value]
            result = gmx.join_arrays(front=result, back=gmx.join_arrays(front=[key], back=value))
    return result


# TODO: (FR4) Use generating function or decorator that can validate kwargs?
# TODO: (FR4) Outputs need to be fully formed and typed in the object returned
#  from the helper (decorated function).
def commandline_operation(executable=None,
                          arguments=(),
                          input_files: dict = None,
                          output_files: dict = None,
                          stdin: str = None,
                          **kwargs):
    """Helper function to define a new operation that executes a subprocess in gmxapi data flow.

    Define a new Operation for a particular executable and input/output parameter set.
    Generate a chain of operations to process the named key word arguments and handle
    input/output data dependencies.

    Arguments:
        executable: name of an executable on the path
        arguments: list of positional arguments to insert at ``argv[1]``
        input_files: mapping of command-line flags to input file names
        output_files: mapping of command-line flags to output file names
        stdin (str): String input to send to STDIN (terminal input) of the executable (optional).

    Multi-line text sent to *stdin* should be joined into a single string.
    E.g.::

        commandline_operation(..., stdin='\\n'.join(list_of_strings) + '\\n')

    If multiple strings are provided to *stdin*, gmxapi will assume an ensemble,
    and will run one operation for each provided string.

    Only string input (:py:func:`str`) to *stdin* is currently supported.
    If you have a use case that requires streaming input or binary input,
    please open an issue or contact the author(s).

    Output:
        The output node of the resulting operation handle contains

        * ``file``: the mapping of CLI flags to filename strings resulting from the ``output_files`` kwarg
        * ``stdout``: A string mapping from process STDOUT.
        * ``stderr``: A string mapping from process STDERR; it will be the
                      error output (if any) if the process failed.
        * ``returncode``: return code of the subprocess.

    """

    # Implementation details: When used in a script, this function returns an
    # instance of an operation. However, because of the dynamic specification of
    # inputs and outputs, each invocation may have the overhead of defining new
    # types to express the data flow topology, regardless of the executable.
    # If this overhead is problematic, consider exposing the intermediate step
    # at which the Operation is fully specified to facilitate reuse.

    ##
    # 1. Define a new operation with outputs from `cli()` plus `file` from `output_files`

    # output_files is essentially passed through, but we need assurance that results
    # will not be published until the rest of the operation has run (i.e. the cli() executable.)

    # Warning: decorating a local function like this is counter to the notion of Operations
    # as portable (importable, serializable/deserializable). The big picture here needs
    # some more consideration.
    # TODO: (NOW) Distinguish portable Operations from relocatable Futures.
    # There is nothing antithetical about objects implementing gmxapi data interfaces
    # that are only resolvable by a certain Context as long as that Context can convey
    # the results to another Context upon request. Re-instantiating Operations is
    # only one way of relocating Futures. In this case, though, the dynamic creation of
    # merged_ops doesn't seem right, and commandline_operation should probably be
    # a proper Operation.
    #
    # TODO: (FR4+) Characterize the `file` dictionary key type:
    #  explicitly sequences rather than maybe-string/maybe-sequence-of-strings
    @gmx.function_wrapper(output={'stdout': str,
                                  'stderr': str,
                                  'returncode': int,
                                  'file': dict})
    def merged_ops(stdout: str = None,
                   stderr: str = None,
                   returncode: int = None,
                   file: dict = None,
                   output: OutputCollectionDescription = None):
        assert stdout is not None
        assert stderr is not None
        assert returncode is not None
        assert file is not None
        assert output is not None
        output.returncode = returncode
        output.stdout = stdout
        output.stderr = stderr
        if returncode == 0:
            output.file = file
        else:
            output.file = {}

    ##
    # 2. Prepare data flow.

    if input_files is None:
        input_files = {}
    if output_files is None:
        output_files = {}
    if isinstance(arguments, (str, bytes)):
        arguments = [arguments]
    command = gmx.concatenate_lists([[executable],
                                     arguments,
                                     filemap_to_flag_list(input_files),
                                     filemap_to_flag_list(output_files)])
    shell = gmx.make_constant(False)
    cli_args = {'command': command,
                'shell': shell}
    cli_args.update(**kwargs)
    if stdin is not None:
        cli_args['stdin'] = str(stdin)

    ##
    # 3. Merge operations
    #
    # Note: Without a `label` argument, repeated calls to cli(**cli_args) should
    # produce references to the same unique resource. Creating this handle
    # separately should not be necessary, but we've got a way to go until we have the
    # fingerprinting and Context resource management we need for that.
    # TODO: ``label`` kwarg
    # TODO: input fingerprinting
    cli_result = cli(**cli_args)
    merged_result = merged_ops(stdout=cli_result.output.stdout,
                               stderr=cli_result.output.stderr,
                               returncode=cli_result.output.returncode,
                               file=output_files,
                               **kwargs)

    # Return an object with an OutputCollection granting access to outputs of
    # cli() and of output_files (as "file")
    return merged_result
