#
# This file is part of the GROMACS molecular simulation package.
#
# Copyright (c) 2018,2019, by the GROMACS development team, led by
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

# - Define macro to check if image conversion using ImageMagick convert
# actually works, because recent changes to it made it necessary to set
# a number of flags in /etc/ImageMagick-6/policy.xml. If a sample
# conversion fails due to this, the user is informed about this being
# a possible issue.
#
#  GMX_TEST_IMAGEMAGICK(VARIABLE)
#
#  VARIABLE will be set in the cache to either true or false
#  if convert is working or not.

function(GMX_TEST_IMAGEMAGICK VARIABLE)
    if(NOT ${ImageMagick_CONVERT_FOUND})
        MESSAGE(STATUS "No image conversion possible without ImageMagick")
        set(value_ OFF)
    elseif(NOT DEFINED ${VARIABLE})
        set(TEMPDIR "${CMAKE_CURRENT_BINARY_DIR}/imagemagicktmp")
        FILE(MAKE_DIRECTORY ${TEMPDIR})
        set(SAMPLE_INPUT "${CMAKE_CURRENT_SOURCE_DIR}/cmake/TestImageMagickConvert.pdf")
        set(SAMPLE_OUTPUT "${TEMPDIR}/TestImageMagickConvert.png")
        execute_process(
            COMMAND ${ImageMagick_convert_EXECUTABLE}  ${SAMPLE_INPUT} -antialias -quality 03 -quiet -pointsize 12 -density 1200 -units PixelsPerInch ${SAMPLE_OUTPUT}
            RESULT_VARIABLE TEST_OUTPUT
            OUTPUT_QUIET
            ERROR_QUIET
            )
        if (EXISTS ${SAMPLE_OUTPUT})
            set(value_ ON)
        else()
            if (GMX_BUILD_MANUAL)
                set(type_ "WARNING")
            else()
                set(type_ "STATUS")
            endif()
            MESSAGE(${type_} "Could not convert sample image, ImageMagick convert can not be used. A possible way to fix it can be found here: https://alexvanderbist.com/posts/2018/fixing-imagick-error-unauthorized")
            set(value_ OFF)
        endif()
        FILE(REMOVE_RECURSE ${TEMPDIR})
    endif()
    if(NOT DEFINED ${VARIABLE})
        set(${VARIABLE} ${value_} CACHE INTERNAL "Test if image conversion works")
        mark_as_advanced(${VARIABLE})
    endif()
endfunction()
