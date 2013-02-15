#
# This file is part of the GROMACS molecular simulation package.
#
# Copyright (c) 2012,2013, by the GROMACS development team, led by
# David van der Spoel, Berk Hess, Erik Lindahl, and including many
# others, as listed in the AUTHORS file in the top-level source
# directory and at http://www.gromacs.org.
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
#
# - Define macro to check large file support
#
#  GMX_TEST_LARGE_FILES(VARIABLE)
#
#  VARIABLE will be set to true if off_t is 64 bits, and fseeko/ftello present.
#  This macro will also set defines necessary enable large file support, for instance
#  _LARGE_FILES
#  _LARGEFILE_SOURCE
#  _FILE_OFFSET_BITS 64  
#
#  However, it is YOUR job to make sure these defines are set in a cmakedefine so they
#  end up in a config.h file that is included in your source if necessary!

MACRO(GMX_TEST_LARGE_FILES VARIABLE)
    IF(NOT DEFINED ${VARIABLE})

        # On most platforms it is probably overkill to first test the flags for 64-bit off_t,
        # and then separately fseeko. However, in the future we might have 128-bit filesystems
        # (ZFS), so it might be dangerous to indiscriminately set e.g. _FILE_OFFSET_BITS=64.

        MESSAGE(STATUS "Checking for 64-bit off_t")

	# First check without any special flags
        TRY_COMPILE(FILE64_OK "${CMAKE_BINARY_DIR}"    
                    "${CMAKE_SOURCE_DIR}/cmake/TestFileOffsetBits.c")
        if(FILE64_OK)
	    MESSAGE(STATUS "Checking for 64-bit off_t - present")			
      	endif(FILE64_OK)

        if(NOT FILE64_OK)	
	    # Test with _FILE_OFFSET_BITS=64
            TRY_COMPILE(FILE64_OK "${CMAKE_BINARY_DIR}"
                        "${CMAKE_SOURCE_DIR}/cmake/TestFileOffsetBits.c"
                        COMPILE_DEFINITIONS "-D_FILE_OFFSET_BITS=64" )
            if(FILE64_OK)
	        MESSAGE(STATUS "Checking for 64-bit off_t - present with _FILE_OFFSET_BITS=64")
                set(_FILE_OFFSET_BITS 64 CACHE INTERNAL "64-bit off_t requires _FILE_OFFSET_BITS=64")
            endif(FILE64_OK)
        endif(NOT FILE64_OK)    

        if(NOT FILE64_OK)
            # Test with _LARGE_FILES
            TRY_COMPILE(FILE64_OK "${CMAKE_BINARY_DIR}"
                        "${CMAKE_SOURCE_DIR}/cmake/TestFileOffsetBits.c"
                        COMPILE_DEFINITIONS "-D_LARGE_FILES" )
            if(FILE64_OK)
                MESSAGE(STATUS "Checking for 64-bit off_t - present with _LARGE_FILES")
                set(_LARGE_FILES 1 CACHE INTERNAL "64-bit off_t requires _LARGE_FILES")
            endif(FILE64_OK)
        endif(NOT FILE64_OK)
	
        if(NOT FILE64_OK)
            # Test with _LARGEFILE_SOURCE
            TRY_COMPILE(FILE64_OK "${CMAKE_BINARY_DIR}"
                        "${CMAKE_SOURCE_DIR}/cmake/TestFileOffsetBits.c"
                        COMPILE_DEFINITIONS "-D_LARGEFILE_SOURCE" )
            if(FILE64_OK)
                MESSAGE(STATUS "Checking for 64-bit off_t - present with _LARGEFILE_SOURCE")
      		set(_LARGEFILE_SOURCE 1 CACHE INTERNAL "64-bit off_t requires _LARGEFILE_SOURCE")
            endif(FILE64_OK)
        endif(NOT FILE64_OK)

        if(NOT FILE64_OK)
            # now check for Windows stuff
            TRY_COMPILE(FILE64_OK "${CMAKE_BINARY_DIR}"
                        "${CMAKE_SOURCE_DIR}/cmake/TestWindowsFSeek.c")
            if(FILE64_OK)
                MESSAGE(STATUS "Checking for 64-bit off_t - present with _fseeki64")
                set(HAVE__FSEEKI64 1 CACHE INTERNAL "64-bit off_t requires _fseeki64")
            endif(FILE64_OK)
        endif(NOT FILE64_OK)

        if(NOT FILE64_OK)
            MESSAGE(STATUS "Checking for 64-bit off_t - not present")
        else(NOT FILE64_OK)

            # Set the flags we might have determined to be required above
            configure_file("${CMAKE_SOURCE_DIR}/cmake/TestLargeFiles.c.cmakein" 
                           "${CMAKE_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/TestLargeFiles.c")

            MESSAGE(STATUS "Checking for fseeko/ftello")            
	    # Test if ftello/fseeko are	available
	    TRY_COMPILE(FSEEKO_COMPILE_OK "${CMAKE_BINARY_DIR}"
                        "${CMAKE_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/TestLargeFiles.c")
	    if(FSEEKO_COMPILE_OK)
                MESSAGE(STATUS "Checking for fseeko/ftello - present")
            endif(FSEEKO_COMPILE_OK)

            if(NOT FSEEKO_COMPILE_OK)
                # glibc 2.2 neds _LARGEFILE_SOURCE for fseeko (but not 64-bit off_t...)
                TRY_COMPILE(FSEEKO_COMPILE_OK "${CMAKE_BINARY_DIR}"
                            "${CMAKE_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/TestLargeFiles.c"
                            COMPILE_DEFINITIONS "-D_LARGEFILE_SOURCE" )
                if(FSEEKO_COMPILE_OK)
                    MESSAGE(STATUS "Checking for fseeko/ftello - present with _LARGEFILE_SOURCE")
                    set(_LARGEFILE_SOURCE 1 CACHE INTERNAL "64-bit fseeko requires _LARGEFILE_SOURCE")
                endif(FSEEKO_COMPILE_OK)
            endif(NOT FSEEKO_COMPILE_OK)

        endif(NOT FILE64_OK)

	if(FSEEKO_COMPILE_OK)
            SET(${VARIABLE} 1 CACHE INTERNAL "Result of test for large file support" FORCE)
            set(HAVE_FSEEKO 1 CACHE INTERNAL "64bit fseeko is available" FORCE)
        else(FSEEKO_COMPILE_OK)
	    if (HAVE__FSEEKI64)
		SET(${VARIABLE} 1 CACHE INTERNAL "Result of test for large file support" FORCE)
		SET(HAVE__FSEEKI64 1 CACHE INTERNAL "Windows 64-bit fseek" FORCE)
	    else (HAVE__FSEEKI64)
                MESSAGE(STATUS "Checking for fseeko/ftello - not found")
                SET(${VARIABLE} 0 CACHE INTERNAL "Result of test for large file support" FORCE)
	    endif (HAVE__FSEEKI64)
        endif(FSEEKO_COMPILE_OK)

    ENDIF(NOT DEFINED ${VARIABLE})
ENDMACRO(GMX_TEST_LARGE_FILES VARIABLE)



