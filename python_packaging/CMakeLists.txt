#
# This file is part of the GROMACS molecular simulation package.
#
# Copyright (c) 2019,2020, by the GROMACS development team, led by
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

# Stub file for future higher-level CMake management. Ultimately, we want to
# drive building and testing from the repository-level CMake configuration, but
# in the initial packaging, tests are somewhat manual and Python package
# installation is driven through the setup.py file in src.
#
# Terminology note:
#  * A "Python package" is a directory that can be imported as a Python module
#    by the Python interpreter using the `import packagename` syntax.
#  * A "distribution archive" is a file in any of several distribution
#    formats (including archive files of specified directory structures)
#    that a Python package management tool or framework
#    (e.g. pypi.org, ``pip``, the ``setuptools`` module)
#    knows how to download, upload, or install
#    (building, if necessary, in the case of a source distribution).
#  * Here, the verb "to package" means to produce a package distribitution archive.
#    This terminology is used because
#     * The verb "to produce a package distribution archive" is awkward to conjugate.
#     * "To build a package" is otherwise ambiguous.
#  * Here, the meaning of "install" varies with context. In the context of CMake
#    or a build system like ``make``, "install" refers to the build system target
#    or command. In the context of a Python package or Python tool, "install"
#    refers to the managed installation of a Python package from a distribution
#    archive to a "site-packages" directory for the Python environment.
#    Note: a (valid and complete) Python does not need to be "installed" to be
#    usable, such as if the ``import`` and dependencies can be resolved with the
#    current ``PYTHONPATH`` environment variable or after ``sys.path.append()``.
#  * "To build a package" means to create a package from sources that may need
#    to be configured or compiled before the package can be imported. Reference
#    the ``distutils`` ``build`` and ``build_ext`` commands. For example, when
#    scikit-build is directed to build a package through the ``setuptools`` interface,
#    it invokes CMake to build **and install** a CMake library target for the C++
#    extension module in the Python package.
#
# For basic building and testing of the Python package, we don't necessarily
# need to create a distribution archive or install the Python package outside
# of the build tree. We can use a CMake ExternalProject to build and install the
# ``_gmxapi`` binary extension module from `python_packaging/src` into the main
# project build tree (note ``ExternalProject_add()`` option ``INSTALL_DIR``).
# If avoiding Python packaging tools in this way, the Python source files (``.py`` files)
# will have to be copied explicitly in additional CMake ``install(FILES)`` directives,
# and the ExternalProject installation directory will have to be added to the
# PYTHONPATH environment variable where the package needs to be available.
#
# Such an ExternalProject could be used as part of the process of packaging a
# distribution archive through arguments to the Python side of ``skbuild``,
# but coupling that machinery is not necessary. The packaging infrastructure may
# be simplest and most maintainable by using scikit-build with its default
# behavior (letting it maintain its own environment and artifacts).
#
# Note that a binary distribution needs to be built for a specific major and minor
# Python release, and is often sensitive to a particular Python distribution or
# packaging system. Instead, we can build `gmxapi` source distributions during
# the CMake build phase at the GROMACS project level. These source distributions
# can be installed into the GROMACS installation path or uploaded to, .e.g,
# PyPI.org so that the binary extension is built in the context of a specific
# (system or user) GROMACS installation and (user) Python environment. For a
# system-wide Python environment, the package needs to be built and installed
# (to the ``site-packages`` directory) for each supported Python interpreter.
# ``setup.py`` can just be invoked with different Python interpreters.
#
# To drive the packaging of a Python package distribution archive
# by a higher-level CMake configuration, a CMakeLists.txt file at this level would
# 1. configure: copy the contents of `src` to the build directory, potentially
#    configuring version info and CMakeToolchain for GROMACS
# 2. build: run (cd $CMAKE_CURRENT_BINARY_DIR/src && python setup.py sdist),
#    producing $CMAKE_CURRENT_BINARY_DIR/src/dist/gmxapi-0.1.0.dev0+fr7.tar.gz or similar
#    (add_custom_command(OUTPUT $sdist_name ...) where sdist_name is determined with
#    `python setup.py --full-name` and a choice of argument for `python setup.py sdist --formats`
# 3. install: copy the source distribution to $GROMACS_INSTALL_DIR/sdist/ or something
# 4. test:
#    0. copy src/test
#    1. require GROMACS library target to be present or importable
#    2. create and activate a Python venv in the build dir
#    3. `pip install <the_sdist_archive_file>` into the venv
#    4. run `pytest copy_of_src/test`

set(PYBIND11_PYTHON_VERSION "3")
add_subdirectory(src)
add_subdirectory(sample_restraint)
