#
# This file is part of the GROMACS molecular simulation package.
#
# Copyright (c) 2013,2014,2015, by the GROMACS development team, led by
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

# TODO: Consider whether this could just be replaced by using
# GMX_DEVELOPER_BUILD to control this as well.
option(GMX_BUILD_MANUAL "Whether to try to configure to build the PDF manual" OFF)
mark_as_advanced(GMX_BUILD_MANUAL)

set(MANUAL_BUILD_IS_POSSIBLE ON)
set(MANUAL_BUILD_NOT_POSSIBLE_REASON)
if (NOT GMX_BUILD_MANUAL)
    # Make sure we only do detection of manual-building dependencies
    # when the user opted in for that.
    set(MANUAL_BUILD_IS_POSSIBLE OFF)
    set(MANUAL_BUILD_NOT_POSSIBLE_REASON "GMX_BUILD_MANUAL is not ON")
elseif (BUILD_IS_INSOURCE)
    # UseLATEX.cmake does not work with an in-source build
    set(MANUAL_BUILD_IS_POSSIBLE OFF)
    set(MANUAL_BUILD_NOT_POSSIBLE_REASON "the build is in-source")
else()
    include(UseLATEX.cmake)

    if(NOT PDFLATEX_COMPILER OR NOT IMAGEMAGICK_CONVERT)
        set(MANUAL_BUILD_IS_POSSIBLE OFF)
        set(MANUAL_BUILD_NOT_POSSIBLE_REASON "pdflatex or some other dependency (ImageMagick convert) is not available")

        # TODO Later, identify other dependencies like bibtex,
        # make_index, date, some graphics conversion program,
        # etc. Perhaps patch UseLATEX.cmake and contribute upstream.
    endif()
endif()

if(MANUAL_BUILD_IS_POSSIBLE)
    LATEX_GET_OUTPUT_PATH(output_dir)

    # Avoid generating a global "html" target that clashes with the
    # html target from program generaion of help text in share/html.
    set(LATEX2HTML_CONVERTER off)

    # This helps make the index look decent.
    set(MAKEINDEX_COMPILER_FLAGS -s ${CMAKE_CURRENT_SOURCE_DIR}/hfill.ist)
    # UseLATEX.cmake is hacked to call Anton Feenstra's subindex script
    # to improve the index generation. Don't replace UseLATEX.cmake
    # blindly from upstream!

    # Here follows various custom CMake "commands" that call our various
    # GROMACS helper scripts to make bits of the manual from the installed
    # binaries and/or source distribution, as appropriate.
    #
    # Only files that can be built should be listed in DEPENDS. Makefile
    # rules are generated for those files.

    # Finally, the command to build the manual.

    ADD_LATEX_DOCUMENT(gromacs.tex
        # Normal LaTeX \included files
        INPUTS algorithms.tex  defunits.tex    implement.tex  macros.tex     special.tex
        analyse.tex     files.tex       topology.tex
        averages.tex    forcefield.tex  gromacs.tex    intro.tex      programs.tex

        # CMake does variable replacement in these files
        CONFIGURE macros.tex

        BIBFILES monster.bib unpubl.bib
        IMAGE_DIRS plots
        DEFAULT_PDF
        USE_INDEX
        )
endif()

if(MANUAL_BUILD_IS_POSSIBLE)
    # An explicit "make manual" will still be required for the manual
    # to be built.
    add_custom_target(manual DEPENDS pdf)
else()
    add_custom_target(manual
        COMMAND ${CMAKE_COMMAND} -E echo
        "Cannot build PDF manual, because ${MANUAL_BUILD_NOT_POSSIBLE_REASON}"
        VERBATIM)
endif()

set(MANUAL_BUILD_IS_POSSIBLE ${MANUAL_BUILD_IS_POSSIBLE} PARENT_SCOPE)
