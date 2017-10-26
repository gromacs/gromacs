#
# This file is part of the GROMACS molecular simulation package.
#
# Copyright (c) 2014,2015,2016,2017, by the GROMACS development team, led by
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
set(IMAGE_CONVERT_POSSIBLE ON)
set(MANUAL_BUILD_NOT_POSSIBLE_REASON)
set(NO_IMAGE_CONVERT_REASON)
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
        if (NOT IMAGEMAGICK_CONVERT)
            set(IMAGE_CONVERT_POSSIBLE OFF)
            set(NO_IMAGE_CONVERT_REASON "Can not convert files for online or pdf manual")
        endif()
        # TODO Later, identify other dependencies like bibtex,
        # make_index, date, some graphics conversion program,
        # etc. Perhaps patch UseLATEX.cmake and contribute upstream.
    endif()

endif()

if(MANUAL_BUILD_IS_POSSIBLE)
    LATEX_GET_OUTPUT_PATH(output_dir)

    # Avoid generating a global "html" target that clashes with the
    # html target from program generation of help text in share/html.
    set(LATEX2HTML_CONVERTER off)

    # This helps make the index look decent.
    set(MAKEINDEX_COMPILER_FLAGS -s ${CMAKE_CURRENT_SOURCE_DIR}/hfill.ist)
    # UseLATEX.cmake is hacked to call Anton Feenstra's subindex script
    # to improve the index generation. Don't replace UseLATEX.cmake
    # blindly from upstream!

    # An explicit "make manual" will still be required for the manual
    # to be built.
    ADD_LATEX_DOCUMENT(Gromacs.tex
        DEFAULT_PDF
        )


    add_custom_target(manual
        DEPENDS sphinx-create-texman
        DEPENDS pdf
        )
else()
        add_custom_target(manual
        COMMAND ${CMAKE_COMMAND} -E echo
        "Cannot build PDF manual, because ${MANUAL_BUILD_NOT_POSSIBLE_REASON}"
        VERBATIM)
endif()
if(IMAGE_CONVERT_POSSIBLE)
    # Need extra step to prepare the png files for input below
    add_custom_target(sphinx-image-conv-png
        DEPENDS sphinx-input sphinx-programs
        COMMAND
            mogrify -format png -density 300 -units PixelsPerInch *.eps *.pdf
        WORKING_DIRECTORY
            "${SPHINX_INPUT_DIR}/web-manual/plots/"
        COMMENT "Preparing png files for webpage"
        VERBATIM
        )
    # Need extra step to prepare the pdf files for input below
    add_custom_target(sphinx-image-conv-pdf
        DEPENDS sphinx-input sphinx-programs
        COMMAND
            mogrify -format pdf -density 300 -units PixelsPerInch *.eps
        WORKING_DIRECTORY
            "${SPHINX_INPUT_DIR}/web-manual/plots/"
        COMMENT "Preparing pdf files for print manual"
        VERBATIM
        )
else()
    # need dummy commands so that we don't break other things
    add_custom_target(sphinx-image-conv-png
        DEPENDS sphinx-input sphinx-programs
        COMMAND ${CMAKE_COMMAND} -E echo
        "Cannot convert images for web manual, skipping this part"
        VERBATIM
        )
    add_custom_target(sphinx-image-conv-pdf
        DEPENDS sphinx-input sphinx-programs
        COMMAND ${CMAKE_COMMAND} -E echo
        "Cannot convert images for pdf manual, skipping this part"
        VERBATIM
        )
endif()

set(MANUAL_BUILD_IS_POSSIBLE ${MANUAL_BUILD_IS_POSSIBLE} PARENT_SCOPE)
set(IMAGE_CONVERT_POSSIBLE ${IMAGE_CONVERT_POSSIBLE} PARENT_SCOPE)
