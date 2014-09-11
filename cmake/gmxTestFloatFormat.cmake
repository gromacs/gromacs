#
# This file is part of the GROMACS molecular simulation package.
#
# Copyright (c) 2009,2014, by the GROMACS development team, led by
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

# - Define macro to determine floating point format properties
#
#  GMX_TEST_FLOAT_FORMAT(FP_IEEE754 FP_BIG_ENDIAN_BYTE FP_BIG_ENDIAN_WORD)
#
#  The thee variables are set to true when:
#  FP_IEEE754          Floating-point numbers are stored in IEEE754 format
#  FP_BIG_ENDIAN_BYTE  The order of bytes in each FP word (4 bytes) is big endian
#  FP_BIG_ENDIAN_WORD  The order of FP words in double precision dwords (2 words) is big endian
#
#  On *most* platforms the two last tests will be the same as the integer endian,
#  big e.g. ARM processors have different byte/word order for floating-point storage,
#  so we figured it is a good idea to test both before relying on the format.

MACRO(GMX_TEST_FLOAT_FORMAT FP_IEEE754 FP_BIG_ENDIAN_BYTE FP_BIG_ENDIAN_WORD)
    IF(NOT DEFINED HAVE_${FP_IEEE754})
        MESSAGE(STATUS "Checking floating point format")

        TRY_COMPILE(HAVE_${FP_IEEE754} "${CMAKE_BINARY_DIR}"    
                    "${CMAKE_SOURCE_DIR}/cmake/TestFloatFormat.cpp"
                    COPY_FILE "${CMAKE_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/TestFloatFormat.bin")  

        if(HAVE_${FP_IEEE754})

            # dont match first/last letter because of string rounding errors :-)
            FILE(STRINGS "${CMAKE_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/TestFloatFormat.bin"
                 GMX_IEEE754_BB_BW LIMIT_COUNT 1 REGEX "ROMACS")
            FILE(STRINGS "${CMAKE_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/TestFloatFormat.bin"
                 GMX_IEEE754_BB_LW LIMIT_COUNT 1 REGEX "CSXGRO")
            FILE(STRINGS "${CMAKE_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/TestFloatFormat.bin"
                 GMX_IEEE754_LB_BW LIMIT_COUNT 1 REGEX "ORGXSC") 
            FILE(STRINGS "${CMAKE_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/TestFloatFormat.bin"
                 GMX_IEEE754_LB_LW REGEX "SCAMOR")
	    FILE(STRINGS "${CMAKE_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/TestFloatFormat.bin"
                 GMX_IEEE754 REGEX "GROMACS|CSXGRO|ORGXSC|SCAMOR")

            # OS X Universal binaries will contain both strings, set it to the host
            IF(GMX_IEEE754_BB_BW  AND  GMX_IEEE754_LB_LW)
                IF(CMAKE_SYSTEM_PROCESSOR MATCHES powerpc)
                    SET(GMX_IEEE754_BB_BW TRUE)
                    SET(GMX_IEEE754_LB_LW FALSE)
                ELSE()
                    SET(GMX_IEEE754_BB_BW FALSE)
                    SET(GMX_IEEE754_LB_LW TRUE)
                ENDIF()
                MESSAGE(STATUS "GMX_TEST_IEEE754_FORMAT found different results, consider setting CMAKE_OSX_ARCHITECTURES or CMAKE_TRY_COMPILE_OSX_ARCHITECTURES to one or no architecture !")
            ENDIF()

            IF(GMX_IEEE754)
                SET(${FP_IEEE754} 1 CACHE INTERNAL "Result of test for IEEE754 FP format" FORCE)
            ENDIF()

            IF(GMX_IEEE754_BB_BW  OR  GMX_IEEE754_BB_LW)
                SET(${FP_BIG_ENDIAN_BYTE} 1 CACHE INTERNAL "Result of test for big endian FP byte order" FORCE)
            ENDIF()

            IF(GMX_IEEE754_BB_BW  OR  GMX_IEEE754_LB_BW)
                SET(${FP_BIG_ENDIAN_WORD} 1 CACHE INTERNAL "Result of test for big endian FP word order" FORCE)
            ENDIF()

        endif()

        # just some informational output for the user
        if(GMX_IEEE754_BB_BW)
            MESSAGE(STATUS "Checking floating point format - IEEE754 (BE byte, BE word)")
        elseif(GMX_IEEE754_BB_LW)
            MESSAGE(STATUS "Checking floating point format - IEEE754 (BE byte, LE word)")
        elseif(GMX_IEEE754_LB_BW)
            MESSAGE(STATUS "Checking floating point format - IEEE754 (LE byte, BE word)")
        elseif(GMX_IEEE754_LB_LW)
            MESSAGE(STATUS "Checking floating point format - IEEE754 (LE byte, LE word)")
        else()
            MESSAGE(STATUS "Checking floating point format - unknown")
            MESSAGE(WARNING "Cannot detect your floating-point format. It is extremely unlikely to be anything else than IEEE754, but if we do not know the endian we need to rely on your OS providing the math functions erfd() and erfcd() rather than using our built-in ones.")
        endif()
    ENDIF()
ENDMACRO(GMX_TEST_FLOAT_FORMAT FP_IEEE754 FP_BIG_ENDIAN_BYTE FP_BIG_ENDIAN_WORD)



