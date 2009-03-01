# - Find FFTW2
# Find the native FFTW2 includes and library, double precision
#
#  FFTW2_INCLUDES    - where to find [d]fftw.h
#  FFTW2_LIBRARIES   - List of libraries when using FFTW.
#  FFTW2_FOUND       - True if FFTW found.

if (FFTW2_INCLUDES)
  # Already in cache, be silent
  set (FFTW2_FIND_QUIETLY TRUE)
endif (FFTW2_INCLUDES)

set(FFTW2_FOUND 0)


foreach(fftw2_name dfftw fftw)
    string(TOUPPER ${fftw2_name} fftw2_uname)
    string(REPLACE "fftw" "rfftw" rfftw2_name ${fftw2_name})
    if(NOT FFTW2_FOUND)
        find_path (FFTW2_INCLUDES ${fftw2_name}.h)
 	find_library (CFFTW2_LIBRARIES  ${fftw2_name})	        
        find_library (RFFTW2_LIBRARIES ${rfftw2_name})
	TRY_COMPILE(FFTW2_FOUND "${CMAKE_BINARY_DIR}"
                    "${CMAKE_SOURCE_DIR}/cmake/TestFFTW2.c"
		    COMPILE_DEFINITIONS "-I${FFTW2_INCLUDES} -DDOUBLE -D${fftw2_uname}" )
    endif(NOT FFTW2_FOUND)
endforeach(fftw2_name dfftw fftw)

if(FFTW2_FOUND)
    set(FFTW2_LIBRARIES "${RFFTW2_LIBRARIES} ${CFFTW2_LIBRARIES}" CACHE STRING "Result of FFTW2 library check" FORCE)
else(FFTW2_FOUND)
    set(FFTW2_INCLUDES 0)
    set(FFTW2_LIBRARIES 0)
endif(FFTW2_FOUND)

# handle the QUIETLY and REQUIRED arguments and set FFTW_FOUND to TRUE if
# all listed variables are TRUE
include (FindPackageHandleStandardArgs)
find_package_handle_standard_args (FFTW2 DEFAULT_MSG FFTW2_LIBRARIES FFTW2_INCLUDES)

mark_as_advanced (RFFTW2_LIBRARIES CFFTW2_LIBRARIES FFTW2_LIBRARIES FFTW2_INCLUDES)


