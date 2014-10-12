#
# This file is part of the GROMACS molecular simulation package.
#
# Copyright (c) 2014,2015, by the GROMACS development team, led by
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

find_package(PythonLibrary REQUIRED)
find_package(NumPy REQUIRED)
find_package(SIP REQUIRED)

include(SIPMacros)
include(PythonMacros)

include_directories(
    ${PYTHON_INCLUDE_PATH}
    ${SIP_INCLUDE_DIR}
    ${GROMACS_INCLUDE_DIRS}
    ${CMAKE_SOURCE_DIR}/src/python/include
    ${NUMPY_INCLUDE_DIRS}
)

add_definitions(-DNPY_NO_DEPRECATED_API=NPY_1_7_API_VERSION)

set(SIP_INCLUDES ${CMAKE_BINARY_DIR} sip)
set(SIP_EXTRA_OPTIONS ${SIP_EXTRA_OPTIONS} -e -o)

file(GLOB common_files_sip sip/*.sip)

file(GLOB options_files_sip sip/options/*.sip)
set(SIP_EXTRA_FILES_DEPEND ${options_files_sip} ${common_files_sip})
add_sip_python_module(gromacs.Options sip/options/Options.sip libgromacs)

file(GLOB topology_files_sip sip/topology/*.sip)
set(SIP_EXTRA_FILES_DEPEND ${topology_files_sip} ${common_files_sip})
add_sip_python_module(gromacs.Topology
    sip/topology/Topology.sip libgromacs)

# This is needed for NumPy, which does not like code in many files
set(SIP_CONCAT_PARTS 1)

file(GLOB selection_files_sip sip/selection/*.sip)
set(SIP_EXTRA_FILES_DEPEND ${selection_files_sip} ${common_files_sip})
add_sip_python_module(gromacs.Selection
    sip/selection/Selection.sip libgromacs)

file(GLOB trajectoryanalysis_files_sip sip/trajectoryanalysis/*.sip)
set(SIP_EXTRA_FILES_DEPEND ${trajectoryanalysis_files_sip} ${common_files_sip})
add_sip_python_module(gromacs.TrajectoryAnalysis
    sip/trajectoryanalysis/TrajectoryAnalysis.sip libgromacs)

python_install(__init__.py ${PYTHON_SITE_PACKAGES_INSTALL_DIR}/gromacs)
