# - Find FFTW2
# Find the native FFTW2 includes and library, double precision
#
#  FFTW2_INCLUDE_DIR    - where to find [d]fftw.h
#  FFTW2_LIBRARIES   - List of libraries when using FFTW.
#  FFTW2_FOUND       - True if FFTW found.
#
# The FFTW2F root installation directory can be provided in the FFTW2F_ROOT_DIR


if (FFTW2_INCLUDE_DIR)
  # Already in cache, be silent
  set (FFTW2_FIND_QUIETLY TRUE)
endif (FFTW2_INCLUDE_DIR)

set(FFTW2_FOUND 0)

file(TO_CMAKE_PATH "$ENV{FFTW2F_ROOT_DIR}" _env_FFTW2F_ROOT_DIR)

foreach(fftw2_name sfftw fftw)
    string(TOUPPER ${fftw2_name} fftw2_uname)
    string(REPLACE "fftw" "rfftw" rfftw2_name ${fftw2_name})
    if(NOT FFTW2_FOUND)
        find_path (FFTW2_INCLUDE_DIR
                    PATHS "${_env_FFTW2F_ROOT_DIR}/include"
                    ${fftw2_name}.h
                    CACHE STRING "Path single precision FFTW2 headers") 
 	find_library (CFFTW2_LIBRARIES  ${fftw2_name}
                        PATHS "${_env_FFTW2F_ROOT_DIR}/lib"
		        CACHE STRING "Single precision CFFTW2 libraries")
        find_library (RFFTW2_LIBRARIES ${rfftw2_name}
                        PATHS "${_env_FFTW2F_ROOT_DIR}/lib"
		        CACHE STRING "Single precision RFFTW2 libraries")
	TRY_COMPILE(FFTW2_FOUND "${CMAKE_BINARY_DIR}"
                    "${CMAKE_SOURCE_DIR}/cmake/TestFFTW2.c"
		    COMPILE_DEFINITIONS "-I${FFTW2_INCLUDE_DIR} -D${fftw2_uname}" )
    endif(NOT FFTW2_FOUND)
endforeach(fftw2_name sfftw fftw)

if(FFTW2_FOUND)
    set(FFTW2_LIBRARIES "${RFFTW2_LIBRARIES} ${CFFTW2_LIBRARIES}" CACHE STRING "Result of FFTW2 library check" FORCE)
else(FFTW2_FOUND)
    set(FFTW2_INCLUDE_DIR 0)
    set(FFTW2_LIBRARIES 0)
endif(FFTW2_FOUND)

# handle the QUIETLY and REQUIRED arguments and set FFTW_FOUND to TRUE if
# all listed variables are TRUE
include (FindPackageHandleStandardArgs)
set(__MSG "Could not find FFTW2F. Provide the fftw2 install directory in the FFTW2F_ROOT_DIR environment variable.")
find_package_handle_standard_args (FFTW2 ${__MSG} FFTW2_LIBRARIES FFTW2_INCLUDE_DIR)

mark_as_advanced (RFFTW2_LIBRARIES CFFTW2_LIBRARIES FFTW2_LIBRARIES FFTW2_INCLUDE_DIR)


