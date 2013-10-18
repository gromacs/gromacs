# Custom build type "ThreadSanitizer", to be used for compiling GROMACS 
# with clang 3.4 or gcc 4.8 (currently pre-release) with ThreadSanitizer
# (aka "TSan") turned on, so that the tests can be run to detect data races.
#
# The main advantage of the clang version is that there can be a
# suppressions file that acts at compile time, though there is no use of 
# that yet. The main advantage of the gcc version is that it can be used
# with a OpenMP (if gomp is compiled with --disable-linux-futex).
#
# Unfortunately, out of the box Thread-MPI provokes several false
# positives. One example is that tMPI_Event_t contains an atomic int
# field "sync" that is initialized before thread spawn. During a
# collective, this is atomically updated by the source thread, and the
# update is observed by sink threads in tMPI_Event_wait, which do a
# yield wait when no change has occured. This means the read can
# happen before the write (by design, whether or not the read is
# atomic), but the surrounding logic prevents action until the write
# has happened. There is no way for the sink thread(s) to avoid
# reading until the write has happened - that is the point of the
# implementation.
#
# This ought to be able to be suppressed, but my attempts to apply
# suppressions on individual functions don't suppress reporting of the
# race event. Applying the suppression to the whole thread-MPI library
# might work, but seems to defeat the point. We want to be able to
# detect mis-use of the primitives provided by thread-MPI.
#
# This means there needs to be a way for this build type to trigger
# the use of the generic mutex-based fallback implementation within
# thread-MPI.
#
# Later, if a blacklist is needed, use something like
# "-fsanitize-blacklist=${CMAKE_SOURCE_DIR}/cmake/thread-sanitizer.supp"
# TODO find a better home for this and other suppression files
set(_flags "-O1 -g -fsanitize=thread")

foreach(_language C CXX)

    string(REPLACE "X" "+" _human_readable_language ${_language})

    if (CMAKE_${_language}_COMPILER_ID MATCHES "GNU")
        set(CMAKE_${_language}_FLAGS_THREADSANITIZER "${_flags} -pie -fPIE" CACHE STRING "${_human_readable_language} flags for thread sanitizer" FORCE)
    else()
        set(CMAKE_${_language}_FLAGS_THREADSANITIZER ${_flags} CACHE STRING "${_human_readable_language} flags for thread sanitizer" FORCE)
    endif()
    mark_as_advanced(CMAKE_${_language}_FLAGS_THREADSANITIZER)
    string(TOUPPER "${CMAKE_BUILD_TYPE}" _cmake_build_type)
    if (_cmake_build_type MATCHES THREADSANITIZER)
        set(TMPI_ATOMICS 0)
        if (NOT((CMAKE_${_language}_COMPILER_ID MATCHES "Clang" AND
                    CMAKE_${_language}_COMPILER_VERSION VERSION_GREATER 3.2.999)
             OR (CMAKE_${_language}_COMPILER_ID MATCHES "GNU" AND
                    CMAKE_${_language}_COMPILER_VERSION VERSION_GREATER 4.7.999)))
            message(FATAL_ERROR "The ThreadSanitizer build is only available with clang ${_human_readable_language} >=3.3 and gnu ${_human_readable_language} >=4.8.")
        endif()
    endif()

endforeach()
