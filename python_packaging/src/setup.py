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

# Note: most of this file is borrowed from
# https://github.com/pybind/cmake_example/commit/31bc276d91985c9bb94e2b4ec12f3fd528971f2c

"""Python setuptools script to build and install the gmxapi Python interface.

The `gmxapi` package requires an existing GROMACS installation, version 2020 or higher.
To specify the GROMACS installation to use, provide a GMXTOOLCHAINDIR
environment variable when running setup.py or `pip`.

Example:
    GMXTOOLCHAINDIR=/path/to/gromacs/share/cmake/gromacs pip install gmxapi

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


class CMakeBuild(build_ext):
    """Derived distutils Command for build_extension.

    See https://github.com/pybind/cmake_example for the current version
    of the sample project from which this is borrowed.
    """
    def build_extension(self, ext):
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

        update_gromacs_client_cmake_args(cmake_args)

        has_pybind = False
        for arg in cmake_args:
            if arg.upper().startswith('-DPYBIND11_ROOT'):
                has_pybind = True
        if not has_pybind:
            pybind_root = pybind11.get_cmake_dir()
            if pybind_root:
                cmake_args.append(f'-Dpybind11_ROOT={pybind_root}')

        subprocess.check_call(
            ["cmake", ext.sourcedir] + cmake_args, cwd=self.build_temp
        )
        subprocess.check_call(
            ["cmake", "--build", "."] + build_args, cwd=self.build_temp
        )


def _find_first_gromacs_suffix(directory):
    dir_contents = os.listdir(directory)
    for entry in dir_contents:
        if entry.startswith('gromacs'):
            return entry.strip('gromacs')


def update_gromacs_client_cmake_args(args: typing.List[str]):
    """Try to convert information from command line environment to usable client CMake stuff.

    This function compartmentalizes details that are likely to evolve with issues
    https://gitlab.com/gromacs/gromacs/-/issues/3273
    and
    https://gitlab.com/gromacs/gromacs/-/issues/3279

    See linked issues for more discussion or to join in the conversation.
    """
    has_gmxapi_dir = False
    gmxapi_DIR = None
    for arg in args:
        if arg.upper().startswith('-DGMXAPI_DIR'):
            gmxapi_DIR = arg.split('=')[1]
            if gmxapi_DIR:
                has_gmxapi_dir = True
            break
    if not has_gmxapi_dir:
        gmxapi_DIR = os.getenv('gmxapi_DIR')
    if not gmxapi_DIR:
        # Infer from GMXRC exports, if available.
        gmxapi_DIR = os.getenv('GROMACS_DIR')

    has_toolchain_file = False
    gmx_toolchain = None
    for arg in args:
        if arg.upper().startswith('-DCMAKE_TOOLCHAIN_FILE'):
            gmx_toolchain = arg.split('=')[1]
            if gmx_toolchain:
                has_toolchain_file = True

    if has_toolchain_file and has_gmxapi_dir:
        return

    gmx_toolchain_dir = os.getenv('GMXTOOLCHAINDIR')
    if gmx_toolchain:
        if gmx_toolchain_dir:
            warnings.warn('Overriding GMXTOOLCHAINDIR environment variable because CMAKE_TOOLCHAIN_FILE CMake '
                          'variable was specified.')
        gmx_toolchain_dir = os.path.dirname(gmx_toolchain)

    if gmx_toolchain_dir is None:
        # Try to guess from standard GMXRC environment variables.
        if gmxapi_DIR is not None:
            if os.path.exists(gmxapi_DIR) and os.path.isdir(gmxapi_DIR):
                share_cmake = os.path.join(gmxapi_DIR, 'share', 'cmake')
                suffix = _find_first_gromacs_suffix(share_cmake)
                if suffix is not None:
                    gmx_toolchain_dir = os.path.join(share_cmake, 'gromacs' + suffix)

    if gmx_toolchain_dir is None or not os.path.exists(gmx_toolchain_dir):
        print(usage)
        raise GmxapiInstallError(
            'Could not configure for GROMACS installation. '
            'Provide GMXTOOLCHAINDIR or CMAKE_TOOLCHAIN_FILE. '
            'See https://manual.gromacs.org/current/gmxapi/userguide/install.html'
        )

    suffix = os.path.basename(gmx_toolchain_dir).strip('gromacs')
    gmx_toolchain = os.path.abspath(os.path.join(gmx_toolchain_dir, 'gromacs-toolchain' + suffix + '.cmake'))

    if not gmxapi_DIR:
        # Example: given /usr/local/gromacs/share/cmake/gromacs/gromacs-toolchain.cmake,
        # we would want /usr/local/gromacs.
        # Note that we could point more directly to the gmxapi-config.cmake but,
        # so far, we have relied on CMake automatically looking into
        # <package>_DIR/share/cmake/<package>/ for such a file.
        # We would need a slightly different behavior for packages that link against
        # libgromacs directly, as sample_restraint currently does.
        gmxapi_DIR = os.path.join(os.path.dirname(gmx_toolchain), '..', '..', '..')

    gmxapi_DIR = os.path.abspath(gmxapi_DIR)

    if not os.path.exists(gmxapi_DIR) or not os.path.isdir(gmxapi_DIR):
        print(usage)
        raise GmxapiInstallError('Please set a valid gmxapi_DIR.')

    if gmxapi_DIR != os.path.commonpath([gmxapi_DIR, gmx_toolchain]):
        raise GmxapiInstallError('GROMACS toolchain file {} is not in gmxapi_DIR {}'.format(
            gmx_toolchain,
            gmxapi_DIR
        ))

    if not has_gmxapi_dir:
        args.append(f'-Dgmxapi_ROOT={gmxapi_DIR}')
    if not has_toolchain_file:
        args.append(f'-DCMAKE_TOOLCHAIN_FILE={gmx_toolchain}')


class GmxapiInstallError(Exception):
    """Error processing setup.py for gmxapi Python package."""

setup(
    ext_modules=[CMakeExtension("gmxapi._gmxapi")],
    cmdclass={"build_ext": CMakeBuild},
    zip_safe=False
)
