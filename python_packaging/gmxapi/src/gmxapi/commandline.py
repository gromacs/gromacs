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

"""
Provide command line operation.
"""

__all__ = ["commandline_operation"]

import collections.abc
import functools
import os
import pathlib
import shutil
import subprocess
import typing
from typing import Iterable
from typing import Union

import gmxapi as gmx
import gmxapi.operation
from gmxapi import exceptions
from gmxapi import logger as root_logger
from gmxapi.datamodel import NDArray
from gmxapi.exceptions import MissingImplementationError
from gmxapi.exceptions import UsageError
from gmxapi.operation import OutputCollectionDescription
from gmxapi.utility import config as _config

# Module-level logger
logger = root_logger.getChild("commandline")
logger.info("Importing {}".format(__name__))


@functools.lru_cache()
def cli_executable() -> pathlib.Path:
    """Report the installed GROMACS command line executable."""
    path = _config().get("gmx_executable", None)
    if path is not None:
        path = pathlib.Path(os.path.abspath(path))
        if path.is_file():
            return path
    raise exceptions.FeatureNotAvailableError("GROMACS installation unavailable.")


@functools.lru_cache()
def cli_bindir() -> pathlib.Path:
    """Report the installed GROMACS binary directory."""
    path = _config().get("gmx_bindir", None)
    if path is not None:
        path = pathlib.Path(os.path.abspath(path))
        if path.is_dir():
            return path
    raise exceptions.FeatureNotAvailableError("GROMACS installation unavailable.")


_seen = set()


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
@gmx.function_wrapper(
    output={
        "directory": str,
        "returncode": int,
        "stderr": str,
        "stdout": str,
    }
)
def cli(
    command: NDArray, shell: bool, env: dict, output, stdin: str = "", _exist_ok=True
):
    """Execute a command line program in a subprocess.

    Configure an executable in a subprocess. Executes when run in an execution
    Context, as part of a work graph or via gmx.run().

    Subprocesses should be run in separate working directories. The directory name is
    the responsibility of the workflow manager and execution manager (and derived from
    the operation id, tied to fingerprinting information).

    Until gmxapi can provide some support for workflow checkpointing, it is unclear
    how to react to the presence of working directories from previous invocations.
    Starting with gmxapi 0.3.0, users may (optionally) set *_exist_ok=False* to
    force an error when an existing working directory for the task is found.
    Otherwise, the command is executed (as with previous gmxapi releases), and
    the command is responsible for determining how to react to the presence of
    prior working files or output, if necessary.

    Shell processing is not enabled, but can be considered for a future version.
    This means that shell expansions such as environment variables, globbing (`*`),
    and other special symbols (like `~` for home directory) are not available.
    This allows a simpler and more robust implementation, as well as a better
    ability to uniquely identify the effects of a command line operation. If you
    think this disallows important use cases, please let us know.

    Arguments:
         env: The environment variables map for the subprocess.
         _exist_ok: If ``True`` (default), it is not an error for the
             working directory to exist already.
         command: a tuple (or list) to be the subprocess arguments, including *executable*
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

    Example:
        Execute a command named ``exe`` that takes a flagged option for file name
        (stored in a local Python variable `my_filename`) and an `origin` flag
        that uses the next three arguments to define a vector.

            >>> my_filename = "somefilename"
            >>> result = cli(('exe', '--origin', 1.0, 2.0, 3.0, '-f', my_filename), shell=False)
            >>> assert hasattr(result, 'file')
            >>> assert hasattr(result, 'stdout')
            >>> assert hasattr(result, 'stderr')
            >>> assert hasattr(result, 'returncode')

    Returns:
        A data structure with attributes for each of the results.

    Result object attributes:
        * file: the mapping of CLI flags to filename strings resulting from the `output` kwarg
        * returncode: return code of the subprocess.
        * stderr: A string mapping from process STDERR; it will be the
                    error output (if any) if the process failed.
        * stdout: A string mapping from process STDOUT.
    """
    # In the operation implementation, we expect the `shell` parameter to be intercepted by the
    # wrapper and set to False.
    if shell:
        raise exceptions.UsageError("Operation does not support shell processing.")

    if stdin == "":
        stdin = None

    if isinstance(command, (str, bytes)):
        command = [command]
    command = list([str(arg) for arg in command])

    executable = shutil.which(command[0])
    if executable is None:
        executable = shutil.which(command[0], path=str(cli_bindir()))
    if executable is None:
        raise exceptions.ValueError(
            '"{}" is not found or not executable.'.format(command[0])
        )
    command[0] = str(executable)

    # The SessionResources concept only exists in the form of the PublishingDataProxy provided as
    # *output* at run time.
    resource_manager: gmxapi.operation.SourceResource = typing.cast(
        gmxapi.operation.DataProxyBase, output
    )._resource_instance
    try:
        op_id = typing.cast(
            gmxapi.operation.ResourceManager, resource_manager
        ).operation_id
    except AttributeError as e:
        # Note that the SourceResource interface does not yet specify a means to get a canonical identifier
        # for the node in the work graph.
        raise MissingImplementationError(
            f"{resource_manager} does not provide required interface. gmxapi.commandline.cli() only "
            f"supports gmxapi.operation.ResourceManager."
        ) from e
    assert op_id

    ensemble_member = typing.cast(
        gmxapi.operation.DataProxyBase, output
    )._client_identifier
    uid = str(op_id)
    if ensemble_member:
        uid += f"_{ensemble_member}"
    assert uid not in _seen
    _seen.add(uid)

    # Note that `abspath` does not resolve symbolic links. (See `pathlib.Path.resolve`)
    # Some users may find it useful to preserve that abstraction.
    _cwd = os.path.abspath(uid)
    if os.path.exists(_cwd):
        message = f"Work directory {_cwd} already exists."
        if _exist_ok:
            logger.info(message)
        else:
            raise UsageError(message)
    os.makedirs(_cwd, exist_ok=_exist_ok)

    # TODO: (FR9) Can OS input/output filehandles be a responsibility of
    #  the code providing 'resources'?

    logger.debug("executing subprocess: %s", " ".join(str(word) for word in command))

    try:
        completed_process = subprocess.run(
            command,
            cwd=_cwd,
            env=env,
            shell=shell,
            input=stdin,
            check=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            encoding="utf-8",
            universal_newlines=True,
        )
        returncode = completed_process.returncode

        stdout = completed_process.stdout
        stderr = completed_process.stderr

    except subprocess.CalledProcessError as e:
        logger.info(f"Non-zero return status when calling {e.cmd} in {_cwd}")
        stdout = e.stdout
        stderr = e.stderr
        returncode = e.returncode

    if stderr is not None:
        if stdout is not None:
            logger.debug("STDOUT:")
            for line in stdout.split("\n"):
                logger.debug(line)
        else:
            logger.debug("STDOUT is empty")
        logger.debug("STDERR:")
        for line in stderr.split("\n"):
            logger.debug(line)
    else:
        logger.debug("STDERR is empty")

    # Publish outputs.
    output.directory = _cwd
    output.stdout = stdout
    output.stderr = stderr
    output.returncode = returncode


def _flatten_dict(mapping: dict):
    """Convert a dict of cli options to a sequence of strings.

    *mapping* values may be scalars or sequences. Non-string sequences will be converted to
    individual elements.

    Returns:
        Iterable of strings.
    """
    for key, value in mapping.items():
        yield str(key)
        if isinstance(value, collections.abc.Iterable) and not isinstance(
            value, (str, bytes)
        ):
            yield from [str(element) for element in value]
        else:
            yield value


@gmx.function_wrapper()
def filemap_to_flag_list(filemap: dict) -> list:
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
    assert isinstance(filemap, dict)
    result = list(_flatten_dict(filemap))
    return result


# TODO: (FR4) Use generating function or decorator that can validate kwargs?
# TODO: (FR4) Outputs need to be fully formed and typed in the object returned
#  from the helper (decorated function).
def commandline_operation(
    executable=None,
    arguments=(),
    input_files: Union[dict, Iterable[dict]] = None,
    output_files: Union[dict, Iterable[dict]] = None,
    stdin: Union[str, Iterable[str]] = None,
    env: Union[dict, Iterable[dict]] = None,
    **kwargs,
):
    """Helper function to define a new operation that executes a subprocess in gmxapi data flow.

    Define a new Operation for a particular executable and input/output parameter set.
    Generate a chain of operations to process the named key word arguments and handle
    input/output data dependencies.

    Arguments:
        env: Optional replacement for the environment variables seen by the subprocess.
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

    .. versionchanged:: 0.3.0
        output_files paths are converted to absolute paths at run time.

    If non-absolute paths are provided to *output_files*, paths are resolved relative to the
    working directory of the command instance (not relative to the working directory of the
    workflow script).

    By default, *executable* runs in a subprocess that inherits its environment
    and resources from the Python interpreter.
    (See https://docs.python.org/3/library/subprocess.html#subprocess.run)

    .. versionadded:: 0.3.1
        If specified, *env* replaces the default environment variables map seen by *executable* in the
        subprocess.

    In addition to controlling environment variables used for user-input, it may be
    necessary to adjust the environment to prevent the subprocess from inheriting variables that it
    should not. This is particularly relevant if the Python script is launched with ``mpiexec`` and
    then ``commandline_wrapper`` is used to launch an MPI-aware executable that may try to manage
    the MPI context. (Reference :issue:`4421`)

    When overriding the environment variables, don't forget to include basic variables
    like PATH that are necessary for the executable to run. `os.getenv` can help.
    E.g. ``commandline_operation(..., env={'PATH': os.getenv('PATH'), ...})``

    Output:
        The output node of the resulting operation handle contains

        * ``directory``: filesystem path that was used as the working directory for the subprocess
        * ``file``: the mapping of CLI flags to filename strings resulting from the ``output_files`` kwarg
        * ``returncode``: return code of the subprocess.
        * ``stderr``: A string mapping from process STDERR; it will be the
                      error output (if any) if the process failed.
        * ``stdout``: A string mapping from process STDOUT.

    .. versionchanged:: 0.3
        Subprocesses run in directories managed by gmxapi.
    .. versionadded:: 0.3
        The *directory* output.

    Working directory names are details of the gmxapi implementation; the naming scheme is not
    yet specified by the API, but is intended to be related to the operation ID.

    Note that re-executing a gmxapi script in the same directory will cause
    commands to be executed again in the same directories. If this presents a
    risk of data corruption (or just wasted compute cycles), you may include
    the key word argument ``_exist_ok=False`` to force an error. Please consider
    contacting the developers through any of the various GROMACS community
    channels to further discuss your use case.
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
    # TODO(3139): Distinguish portable Operations from relocatable Futures.
    # There is nothing antithetical about objects implementing gmxapi data interfaces
    # that are only resolvable by a certain Context as long as that Context can convey
    # the results to another Context upon request. Re-instantiating Operations is
    # only one way of relocating Futures. In this case, though, the dynamic creation of
    # merged_ops doesn't seem right, and commandline_operation should probably be
    # a proper Operation.
    #
    # TODO: (FR4+) Characterize the `file` dictionary key type:
    #  explicitly sequences rather than maybe-string/maybe-sequence-of-strings
    @gmx.function_wrapper(
        output={
            "directory": str,
            "file": dict,
            "returncode": int,
            "stderr": str,
            "stdout": str,
        }
    )
    def merged_ops(
        stdout: str = None,
        stderr: str = None,
        returncode: int = None,
        file: dict = None,
        directory: str = None,
        output: OutputCollectionDescription = None,
    ):
        assert stdout is not None
        assert stderr is not None
        assert returncode is not None
        assert file is not None
        assert output is not None
        output.returncode = returncode
        output.stdout = stdout
        output.stderr = stderr
        output.directory = directory
        file_map = file.copy()
        for key, value in file_map.items():
            if not os.path.isabs(value):
                value = os.path.abspath(os.path.join(directory, value))
            file_map[key] = value
        if returncode == 0:
            output.file = file_map
        else:
            output.file = {}

    ##
    # 2. Prepare data flow.

    if input_files is None:
        input_files = {}
    if output_files is None:
        output_files = {}
    try:
        executable = str(executable)
    except Exception as e:
        raise gmxapi.exceptions.TypeError(
            "gmxapi typing currently requires paths and names to be strings. *executable* argument is "
            f"{type(executable)}."
        )
    if isinstance(arguments, (str, bytes)):
        arguments = [arguments]
    input_options = filemap_to_flag_list(input_files).output.data
    output_options = filemap_to_flag_list(output_files).output.data
    command = gmx.concatenate_lists(
        [[executable], arguments, input_options, output_options]
    )
    shell = gmx.make_constant(False)

    if env is None:

        @gmx.function_wrapper(
            # allow_duplicate=True
        )
        def _env() -> dict:
            import os

            return dict(os.environ)

        env = _env().output.data
    cli_args = {"command": command, "shell": shell, "env": env}
    cli_args.update(**kwargs)
    if stdin is not None:
        cli_args["stdin"] = stdin

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
    merged_result = merged_ops(
        stdout=cli_result.output.stdout,
        stderr=cli_result.output.stderr,
        returncode=cli_result.output.returncode,
        directory=cli_result.output.directory,
        file=output_files,
        **kwargs,
    )

    # Return an object with an OutputCollection granting access to outputs of
    # cli() and of output_files (as "file")
    return merged_result
