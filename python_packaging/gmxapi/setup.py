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

# Note: most of this file is borrowed from
# https://github.com/pybind/cmake_example/commit/31bc276d91985c9bb94e2b4ec12f3fd528971f2c

"""Python setuptools script to build and install the gmxapi Python interface.

The `gmxapi` package requires an existing GROMACS installation, version 2020 or higher.
If the usual GROMACS environment variables are detected (such as from sourcing GMXRC),
the corresponding GROMACS installation will be used. Otherwise, provide CMake options
to the Python package installer.
Define ``gmxapi_ROOT`` to the GROMACS installation directory (``-Dgmxapi_ROOT=/path/to/gromacs``).
Help CMake find the GROMACS compiler toolchain with ``-C  /path/to/gromacs/share/cmake/gromacs/gromacs-hints.cmake``.

If this script is unable to locate GROMACS as above, it will check for a
GMXTOOLCHAINDIR environment variable, giving the location of a gromacs-toolchain.cmake file.
The toolchain file is not provided by GROMACS 2022 and above, so support for GMXTOOLCHAINDIR
is deprecated.

Example:
    GROMACS_DIR=/path/to/gromacs
    CMAKE_ARGS="-Dgmxapi_ROOT=$GROMACS_DIR -C $GROMACS_DIR/share/cmake/gromacs/gromacs-hints.cmake"
    pip install gmxapi

See https://manual.gromacs.org/current/gmxapi/userguide/install.html for more information.
"""

import os
import re
import subprocess
import sys
import typing
import warnings

from setuptools import setup, Extension
from setuptools.command.build_ext import build_ext
from setuptools.config import read_configuration

# Currently we have only partially migrated to PEP 517/518.
# It is unclear how we will access configuration details
# like options.package_data in the long run, but for now we
# can load the contents of setup.cfg into setup.py.
conf_dict = read_configuration(os.path.join(os.path.dirname(__file__), "setup.cfg"))

usage = __doc__[2:]

# gmxapi does not officially support Windows environments because GROMACS does not have automated testing
# infrastructure to verify correct functionality. However, we can try to be friendly or prepare for a possible future
# in which we can support more platforms.
# Convert distutils Windows platform specifiers to CMake -A arguments
PLAT_TO_CMAKE = {
    "win32": "Win32",
    "win-amd64": "x64",
    "win-arm32": "ARM",
    "win-arm64": "ARM64",
}


# A CMakeExtension needs a sourcedir instead of a file list.
# The name must be the _single_ output extension from the CMake build.
class CMakeExtension(Extension):
    def __init__(self, name, sourcedir=""):
        Extension.__init__(self, name, sources=[])
        self.sourcedir = os.path.abspath(sourcedir)


class CMakeBuild(typing.cast("setuptools.Command", build_ext)):
    """Derived distutils Command for build_extension.

    See https://github.com/pybind/cmake_example for the current version
    of the sample project from which this is borrowed.
    """

    def build_extension(self, ext):
        # The pybind11 package is only needed for `build_ext`, and may not be installed
        # in the caller's environment. By the time this function is called, though,
        # build dependencies will have been checked.
        import pybind11

        extdir = os.path.abspath(os.path.dirname(self.get_ext_fullpath(ext.name)))

        # required for auto-detection & inclusion of auxiliary "native" libs
        if not extdir.endswith(os.path.sep):
            extdir += os.path.sep

        debug = int(os.environ.get("DEBUG", 0)) if self.debug is None else self.debug
        cfg = "Debug" if debug else "Release"

        # CMake lets you override the generator - we need to check this.
        # Can be set with Conda-Build, for example.
        cmake_generator = os.environ.get("CMAKE_GENERATOR", "")

        cmake_args = [
            "-DCMAKE_LIBRARY_OUTPUT_DIRECTORY={}".format(extdir),
            "-DCMAKE_BUILD_TYPE={}".format(cfg),  # not used on MSVC, but no harm
        ]
        build_args = []
        # Adding CMake arguments set as environment variable
        # (needed e.g. to build for ARM OSx on conda-forge)
        if "CMAKE_ARGS" in os.environ:
            cmake_args += [item for item in os.environ["CMAKE_ARGS"].split(" ") if item]

        if self.compiler.compiler_type != "msvc":
            # Using Ninja-build since it a) is available as a wheel and b)
            # multithreads automatically. MSVC would require all variables be
            # exported for Ninja to pick it up, which is a little tricky to do.
            # Users can override the generator with CMAKE_GENERATOR in CMake
            # 3.15+.
            if not cmake_generator:
                try:
                    import ninja  # noqa: F401

                    cmake_args += ["-GNinja"]
                except ImportError:
                    pass

        else:

            # Single config generators are handled "normally"
            single_config = any(x in cmake_generator for x in {"NMake", "Ninja"})

            # CMake allows an arch-in-generator style for backward compatibility
            contains_arch = any(x in cmake_generator for x in {"ARM", "Win64"})

            # Specify the arch if using MSVC generator, but only if it doesn't
            # contain a backward-compatibility arch spec already in the
            # generator name.
            if not single_config and not contains_arch:
                cmake_args += ["-A", PLAT_TO_CMAKE[self.plat_name]]

            # Multi-config generators have a different way to specify configs
            if not single_config:
                cmake_args += [
                    "-DCMAKE_LIBRARY_OUTPUT_DIRECTORY_{}={}".format(cfg.upper(), extdir)
                ]
                build_args += ["--config", cfg]

        if sys.platform.startswith("darwin"):
            # Cross-compile support for macOS - respect ARCHFLAGS if set
            archs = re.findall(r"-arch (\S+)", os.environ.get("ARCHFLAGS", ""))
            if archs:
                cmake_args += ["-DCMAKE_OSX_ARCHITECTURES={}".format(";".join(archs))]

        # Set CMAKE_BUILD_PARALLEL_LEVEL to control the parallel build level
        # across all generators.
        if "CMAKE_BUILD_PARALLEL_LEVEL" not in os.environ:
            # self.parallel is a Python 3 only way to set parallel jobs by hand
            # using -j in the build_ext call, not supported by pip or PyPA-build.
            if hasattr(self, "parallel") and self.parallel:
                # CMake 3.12+ only.
                build_args += ["-j{}".format(self.parallel)]

        if not os.path.exists(self.build_temp):
            os.makedirs(self.build_temp)

        cmake_args = update_gromacs_client_cmake_args(cmake_args)
        try:
            hints_file_index = cmake_args.index("-C")
        except ValueError:
            pass
        else:
            hints_file = cmake_args[hints_file_index + 1]
            if not hints_file or not os.path.exists(hints_file):
                raise GmxapiInstallError(
                    f"CMake hints file {hints_file} does not exist."
                )
            # Insert before the end, in case there are trailing `-- ...` args.
            cmake_args.insert(hints_file_index, f"-D_setuppy_cmake_hints={hints_file}")

        has_pybind = False
        for arg in cmake_args:
            if arg.upper().startswith("-DPYBIND11_ROOT"):
                has_pybind = True
        if not has_pybind:
            pybind_root = pybind11.get_cmake_dir()
            if pybind_root:
                cmake_args.append(f"-Dpybind11_ROOT={pybind_root}")

        subprocess.check_call(
            ["cmake", ext.sourcedir] + cmake_args, cwd=self.build_temp
        )
        subprocess.check_call(
            ["cmake", "--build", "."] + build_args, cwd=self.build_temp
        )

    def run(self):
        """Build extensions in build directory, then copy if --inplace

        Extends setuptools.command.build_ext.run()
        """
        # Note that for installations in "develop" mode (or with `--in-place`),
        # setuptools sets `inplace=1` to cause build_ext.run() to copy built
        # extensions back to the source directory. The setuptools `develop`
        # command does not call `build_py` and does not provide a good
        # alternative insertion point for special handling of package data files
        # (such as the invocation of `build_ext` with `inplace=1`).
        # Since our data file is generated with the extension, this seems
        # like a reasonable place to extend the handling in develop mode.
        from distutils.file_util import copy_file

        super(CMakeBuild, self).run()
        if self.inplace:
            # We follow the model of build_ext.copy_extensions_to_source().
            build_py = self.get_finalized_command("build_py")
            for ext in self.extensions:
                fullname = self.get_ext_fullname(ext.name)
                modpath = fullname.split(".")
                package = ".".join(modpath[:-1])
                package_dir = build_py.get_package_dir(package)

                for filename in conf_dict["options"]["package_data"][package]:
                    src_filename = os.path.join(self.build_lib, *modpath[:-1], filename)
                    dest_filename = os.path.join(
                        package_dir, os.path.basename(filename)
                    )

                    copy_file(
                        src_filename,
                        dest_filename,
                        verbose=self.verbose,
                        dry_run=self.dry_run,
                    )


def _find_first_gromacs_suffix(directory):
    dir_contents = os.listdir(directory)
    for entry in dir_contents:
        if entry.startswith("gromacs"):
            return entry.strip("gromacs")


def update_gromacs_client_cmake_args(args: typing.Sequence[str]) -> typing.List[str]:
    """Try to convert information from command line environment to usable client CMake stuff.

    Normalize user input and automated GROMACS detection to produce a list of CMake arguments
    containing a ``-Dgmxapi_ROOT`` and, if possible, hints for generating the CMake build tool
    chain.

    Args:
        args: CMake args provided by CMakeBuild instance, including user input from CMAKE_ARGS.

    First, we must determine a single value of ``gmxapi_ROOT``.

    If ``CMAKE_ARGS`` contains a ``-C`` option, the value is checked for consistency with the
    ``gmxapi_ROOT`` path.

    If ``CMAKE_ARGS`` contains both ``-C`` and ``-DCMAKE_TOOLCHAIN_FILE`` options,
    both are passed along to CMake and we assume the user knows what they are doing.
    If ``CMAKE_ARGS`` contains neither a ``-C`` option, nor a ``-DCMAKE_TOOLCHAIN_FILE`` option,
    this script attempts to locate a ``gromacs-hints.cmake`` file in order to generate
    a ``-C`` option to add to the CMake arguments for CMakeBuild.

    If the script is unable to locate a ``gromacs-hints.cmake`` file, we fall back
    to the gmxapi 0.2 scheme for ``GMXTOOLCHAINDIR``.

    This function compartmentalizes details that are likely to evolve with issues
    https://gitlab.com/gromacs/gromacs/-/issues/3273
    and
    https://gitlab.com/gromacs/gromacs/-/issues/3279

    See linked issues for more discussion or to join in the conversation.
    """
    args = list(str(arg) for arg in args)

    # Existing input for gmxapi_root
    gmxapi_root_arg = False

    # Value of gmxapi_root for use in this function.
    gmxapi_root = None
    for arg in args:
        if arg.upper().startswith("-DGMXAPI_DIR") or arg.upper().startswith(
            "-DGMXAPI_ROOT"
        ):
            _dir = arg.split("=")[1]
            if _dir:
                gmxapi_root_arg = arg
                gmxapi_root = os.path.abspath(_dir)
            break

    # The CMake hints file, if available.
    # Note that the `-C` flag occurs as a distinct CLI argument from the value (the file path).
    cmake_hints = None
    for i, arg in enumerate(args):
        if arg.upper() == "-C":
            cmake_hints = args[i + 1]
            if cmake_hints:
                cmake_hints = os.path.abspath(cmake_hints)
            if not os.path.exists(cmake_hints):
                raise GmxapiInstallError(f"Hints file {cmake_hints} does not exist.")

    # CMake toolchain argument already present in args.
    toolchain_file_arg = None
    # The CMake toolchain file, if available.
    gmx_toolchain = None
    for arg in args:
        if arg.upper().startswith("-DCMAKE_TOOLCHAIN_FILE"):
            gmx_toolchain = arg.split("=")[1]
            if gmx_toolchain:
                toolchain_file_arg = arg
                gmx_toolchain = os.path.abspath(gmx_toolchain)
                if os.getenv("GMXTOOLCHAINDIR"):
                    warnings.warn(
                        "Overriding GMXTOOLCHAINDIR environment variable because"
                        "CMAKE_TOOLCHAIN_FILE CMake variable was specified."
                    )

    if (toolchain_file_arg or cmake_hints) and gmxapi_root_arg:
        # The user has provided enough information. No automatigical behavior required.
        return args
    # Else, try to find the needed configuration, and update args.

    # If toolchain file was not provided explicitly, look for the documented environment variable.
    if gmx_toolchain:
        gmx_toolchain_dir = os.path.dirname(gmx_toolchain)
    else:
        gmx_toolchain_dir = os.getenv("GMXTOOLCHAINDIR")

    # Final attempts to locate gmxapi_ROOT.
    if not gmxapi_root_arg:
        gmxapi_root = os.getenv("gmxapi_ROOT", default=os.getenv("gmxapi_DIR"))
    if not gmxapi_root:
        # Infer from GMXRC exports, if available.
        gmxapi_root = os.getenv("GROMACS_DIR")
    # Try to infer from other input.
    # Example: given /usr/local/gromacs/share/cmake/gromacs/gromacs-toolchain.cmake
    # or /usr/local/gromacs/share/cmake/gromacs/gromacs-hints.cmake,
    # we would want /usr/local/gromacs.
    # Note that we could point more directly to the gmxapi-config.cmake but,
    # so far, we have relied on CMake automatically looking into
    # <package>_DIR/share/cmake/<package>/ for such a file.
    if not gmxapi_root:
        if cmake_hints:
            gmxapi_root = os.path.join(os.path.dirname(cmake_hints), "..", "..", "..")
        elif gmx_toolchain_dir:
            gmxapi_root = os.path.join(gmx_toolchain_dir, "..", "..", "..")
    if (
        not gmxapi_root
        or not os.path.exists(gmxapi_root)
        or not os.path.isdir(gmxapi_root)
    ):
        print(usage)
        raise GmxapiInstallError("Please set a valid gmxapi_ROOT.")
    gmxapi_root = os.path.abspath(gmxapi_root)
    if not gmxapi_root_arg:
        args.append(f"-Dgmxapi_ROOT={gmxapi_root}")

    # If we need to keep guessing for cmake hints or toolchain files, make a final effort to find
    # the GROMACS-provided share/cmake/gromacs* directory.
    if not gmx_toolchain_dir:
        # Try to guess from standard GMXRC environment variables.
        if gmxapi_root:
            share_cmake = os.path.join(gmxapi_root, "share", "cmake")
            suffix = _find_first_gromacs_suffix(share_cmake)
            if suffix is not None:
                gmx_toolchain_dir = os.path.join(share_cmake, "gromacs" + suffix)

    if gmx_toolchain_dir and os.path.exists(gmx_toolchain_dir):
        suffix = os.path.basename(gmx_toolchain_dir).strip("gromacs")
        if not cmake_hints:
            cmake_hints = os.path.abspath(
                os.path.join(gmx_toolchain_dir, "gromacs-hints" + suffix + ".cmake")
            )
        if os.path.exists(cmake_hints):
            if "-C" not in args:
                args.extend(("-C", str(cmake_hints)))
        elif not toolchain_file_arg:
            # Only bother guessing a toolchain file if no hints file is available.
            if not gmx_toolchain:
                gmx_toolchain = os.path.abspath(
                    os.path.join(
                        gmx_toolchain_dir, "gromacs-toolchain" + suffix + ".cmake"
                    )
                )
            if os.path.exists(gmx_toolchain):
                args.append(f"-DCMAKE_TOOLCHAIN_FILE={gmx_toolchain}")

    return args


class GmxapiInstallError(Exception):
    """Error processing setup.py for gmxapi Python package."""


setup(
    ext_modules=[CMakeExtension("gmxapi._gmxapi")],
    cmdclass={"build_ext": CMakeBuild},
    zip_safe=False,
)
