# - Find FFTW3F
# Find the native FFTW3 includes and library, single precision
#
#  FFTW3F_INCLUDE_DIR   where to find fftw3.h
#  FFTW3F_LIBRARIES     List of libraries when using FFTW.
#  FFTW3F_FOUND         True if FFTW found.
#  FFTW3F_THREADS       True if threaded version was found

if (FFTW3F_INCLUDE_DIR AND FFTW3F_LIBRARIES)
  # Already in cache, be silent
  set (FFTW3F_FIND_QUIETLY TRUE)
endif (FFTW3F_INCLUDE_DIR AND FFTW3F_LIBRARIES)

find_path (FFTW3F_INCLUDE_DIR fftw3.h 
	CACHE STRING "Path to headers for single precision FFTW3")

find_library (FFTW3F_LIBRARIES 
                NAMES fftw3f
                PATHS "${FFTW3F_INCLUDE_DIR}/../lib"
                CACHE STRING "Single precision FFTW3 libraries")

# if we are using threading try to find threaded fftw3
if(GMX_OPENMP)
    find_library (_FFTW3F_thread_libs
                    NAMES fftw3f_threads
                    PATHS "${FFTW3F_INCLUDE_DIR}/../lib")
    if(_FFTW3F_thread_libs)
        set(_threads True)
        list(APPEND FFTW3F_LIBRARIES ${_FFTW3F_thread_libs})
    else()
        message(WARNING "Threaded FFTW3F not found, FFT will run in serial (slower)!")
    endif()
endif()


# handle the QUIETLY and REQUIRED arguments and set FFTW_FOUND to TRUE if
# all listed variables are TRUE
include (FindPackageHandleStandardArgs)
find_package_handle_standard_args (FFTW3F DEFAULT_MSG FFTW3F_LIBRARIES FFTW3F_INCLUDE_DIR)

if(_threads)
    set(FFTW3F_THREADS CACHE STRING "Are we using threads FFTW3F?" "True")
else()
    set(FFTW3F_THREADS CACHE STRING "Are we using threads FFTW3F?" "False")
endif()

mark_as_advanced (FFTW3F_LIBRARIES FFTW3F_INCLUDE_DIR FFTW3F_THREADS)
