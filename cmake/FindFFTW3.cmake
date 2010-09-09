# - Find FFTW3
# Find the native FFTW3 includes and library, double precision
#
#  FFTW3_INCLUDE_DIR    where to find fftw3.h
#  FFTW3_LIBRARIES      List of libraries when using FFTW.
#  FFTW3_FOUND          True if FFTW found.
#  FFTW3_THREADS        True if threaded version was found

if (FFTW3_INCLUDE_DIR AND FFTW3_LIBRARIES)
  # Already in cache, be silent
  set (FFTW3_FIND_QUIETLY TRUE)
endif (FFTW3_INCLUDE_DIR AND FFTW3_LIBRARIES)

find_path (FFTW3_INCLUDE_DIR fftw3.h
		CACHE STRING "Path to headers for double precision FFTW3")

find_library (FFTW3_LIBRARIES 
                NAMES fftw3
                PATHS "${FFTW3_INCLUDE_DIR}/../lib"
                CACHE STRING "Double precision FFTW3 libraries")

# if we are using threading try to find threaded fftw3
if(GMX_OPENMP)
    find_library (_FFTW3_thread_libs
                    NAMES fftw3_threads
                    PATHS "${FFTW3F_INCLUDE_DIR}/../lib")
     if(_FFTW3_thread_libs)
        set(_theads True)
        list(APPEND FFTW3_LIBRARIES ${_FFTW3_thread_libs})
    else()
        message(WARNING "Threaded FFTW3F not found, FFT will run in serial (slower)!")
    endif()
endif()


# handle the QUIETLY and REQUIRED arguments and set FFTW_FOUND to TRUE if
# all listed variables are TRUE
include (FindPackageHandleStandardArgs)
find_package_handle_standard_args (FFTW3 DEFAULT_MSG FFTW3_LIBRARIES FFTW3_INCLUDE_DIR)

if(_threads)
    set(FFTW3_THREADS CACHE STRING "Are we using threads FFTW3?" "True")
else()
    set(FFTW3_THREADS CACHE STRING "Are we using threads FFTW3?" "False")
endif()

mark_as_advanced (FFTW3_LIBRARIES FFTW3_INCLUDE_DIR FFTW3_THREADS)


