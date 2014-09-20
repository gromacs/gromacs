#
# This file is part of the GROMACS molecular simulation package.
#
# Copyright (c) 2013,2014, by the GROMACS development team, led by
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

if(${CMAKE_BINARY_DIR} STREQUAL ${CMAKE_SOURCE_DIR})
    # UseLATEX.cmake does not work with an in-source build
    message(STATUS "Cannot build the manual when building in the source directory")
    set(MANUAL_BUILD_IS_POSSIBLE OFF)
else()
    INCLUDE(UseLATEX.cmake)

    if(NOT PDFLATEX_COMPILER OR NOT IMAGEMAGICK_CONVERT)
        # No pdflatex was found, so don't build the manual.
        message(STATUS "A required dependency of the manual (pdflatex, ImageMagick's convert) was not found, so the manual build will not be available")
        set(MANUAL_BUILD_IS_POSSIBLE OFF)

        # TODO Later, identify other dependencies like bibtex,
        # make_index, date, some graphics conversion program,
        # etc. Perhaps patch UseLATEX.cmake and contribute upstream.
    else()
        set(MANUAL_BUILD_IS_POSSIBLE ON)
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
    #
    # The .mdp options section is now intended to be built directly
    # from a GROMACS source directory. This facilitates drafting
    # updates in that HTML file without copying files by hand.

    ADD_CUSTOM_COMMAND(OUTPUT ${output_dir}/mdp_opt.tex
        COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/mkmdp ARGS ${CMAKE_CURRENT_SOURCE_DIR}/../old-html
        DEPENDS mkmdp ${CMAKE_CURRENT_SOURCE_DIR}/../old-html/online/mdp_opt.html
        )

    # Finally, the command to build the manual.

    ADD_LATEX_DOCUMENT(gromacs.tex
        # Normal LaTeX \included files
        INPUTS algorithms.tex  defunits.tex    implement.tex  macros.tex     special.tex
        analyse.tex     files.tex       install.tex    topology.tex
        averages.tex    forcefield.tex  gromacs.tex    intro.tex      programs.tex

        # CMake does variable replacement in these files
        CONFIGURE macros.tex

        # These files we're responsible for creating in the
        # add_custom_targets() above. They should not be in the git
        # repository, or its directory, or the dependencies will not work
        # properly.
        DEPENDS mdp_opt.tex

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
    # TODO Arrange for the "make manual" target to explain that this can't
    # be done
endif()

set(MANUAL_BUILD_IS_POSSIBLE ${MANUAL_BUILD_IS_POSSIBLE} PARENT_SCOPE)
