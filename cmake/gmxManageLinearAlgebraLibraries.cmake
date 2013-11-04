# Helper macro for the function below. We treat BLAS and LAPACK the same
# way, so there's no reason to duplicate the logic for them
#
# MESSSAGE_TEXT variable is accumulated while checks for these
# libraries fail, but the message is only emitted if we are forced to
# fall back on the internal version.
#
# Arguments should be:
#     name                 "BLAS" or "LAPACK"
#     function_in_library  the name of a function to use in a linking test of that library
macro(manage_linear_algebra_library name function_in_library)
    set(_library_was_found 0)
    set(_find_quietly FALSE)
    set(_user_var_changed FALSE)
    if(NOT DEFINED GMX_EXTERNAL_${name} OR
       (GMX_EXTERNAL_${name} AND NOT "${GMX_${name}_USER}" STREQUAL "${GMX_${name}_USER_PREV}"))
        set(_user_var_changed TRUE)
    endif()
    if(DEFINED GMX_EXTERNAL_${name} AND NOT _user_var_changed)
        set(_find_quietly TRUE)
    endif()

    # We could consider printing status messages at the beginning and
    # end, which would require caching whether the previous provider
    # was user/MKL/external/internal. It's possible (we do it for
    # FFT), but the number of times the user changes these is pretty
    # low, so let's solve that one in master branch when we have
    # better CMake gear to support it.
    if(GMX_EXTERNAL_${name} OR NOT DEFINED GMX_EXTERNAL_${name})
        set(GMX_${name}_USER_PREV ${GMX_${name}_USER} CACHE INTERNAL
            "Previous value of GMX_${name}_USER (to detect changes)" FORCE)
        set(_message_text)
        # Check for user-specified libraries if external libraries have
        # been specified (which is the default).
        if(GMX_${name}_USER)
            set(_libraries_to_link ${GMX_${name}_USER})
            set(_library_was_found 1)

            if(NOT _find_quietly)
                set(CMAKE_REQUIRED_LIBRARIES ${GMX_${name}_USER})
                if(_user_var_changed)
                    unset(_${name}_user_works CACHE)
                endif()
                message(STATUS "Checking that user ${name} library ${GMX_${name}_USER} works")
                check_function_exists(${function_in_library} _${name}_user_works)
                if(NOT _${name}_user_works)
                    message(WARNING "GMX_${name}_USER library ${GMX_${name}_USER} was specified, but it may not provide ${name}. We are proceeding by assuming you know what you are doing and that linking F77-style to this library will work.")
                endif()

                if(HAVE_LIBMKL)
                    message(STATUS "MKL and GMX_${name}_USER were both specified. Using the latter for ${name}.")
                endif()
            endif()
        endif()

        if(NOT _library_was_found AND HAVE_LIBMKL)
            set(CMAKE_REQUIRED_LIBRARIES ${FFT_LINKER_FLAGS} ${FFT_LIBRARIES})
            # This may also not work correctly if the user changes
            # MKL_LIBRARIES after the first run. However,
            # MKL_LIBRARIES is only needed for icc version < 11, or
            # for trying to use MKL with a non-Intel compiler, and we
            # can live with that for now.
            check_function_exists(${function_in_library} _${name}_mkl_works)
            if(_${name}_mkl_works)
                # If we ever need to compile with MKL linear algebra and
                # not with FFT supplied by MKL, uncomment the next line
                # (and probably tweak other things).
#                list(APPEND LINEAR_ALGEBRA_LIBRARIES ${FFT_LINKER_FLAGS} ${FFT_LIBRARIES})
                set(_library_was_found 1)
            else()
                set(_message_text "Intel's MKL was specified, and it should provide ${name}, but it does not. ")
            endif()
        endif()

        # If detection of ${name} has never run, or none of the preceding
        # detection succeeded, try to detect ${name} in the CMake
        # detection paths, etc.
        if (NOT _library_was_found)
            set(${name}_FIND_QUIETLY _find_quietly)
            # Note that this finds all kinds of system libraries,
            # including Apple's Accelerate Framework (and perhaps MKL for
            # icc < 11).
            find_package(${name})
            if (${name}_FOUND)
                set(_libraries_to_link ${${name}_LIBRARIES})
                set(_library_was_found 1)
            endif()
        endif()

        if (NOT _library_was_found AND NOT _find_quietly)
            message("${_message_text}A ${name} library was not found by CMake in the paths available to it. Falling back on the GROMACS internal version of the ${name} library instead. This is fine for normal usage.")
        endif()
    endif()

    # Default behaviour is to try to use an external library, but fall
    # back on the internal one if none is found.
    set(GMX_EXTERNAL_${name} ${_library_was_found} CACHE BOOL "Use a ${name} library that is external to GROMACS if possible (ON), or the internal GROMACS one (OFF)")
    mark_as_advanced(GMX_EXTERNAL_${name})
    # Default behaviour is to use a library found on the system or in
    # GROMACS. The user must actively set GMX_${name}_USER if they
    # want to specify a library.
    set(GMX_${name}_USER "" CACHE BOOL "Use a ${name} library found on the system (OFF), or a ${name} library supplied by the user (any other value, which is a full path to that ${name} library)")
    mark_as_advanced(GMX_${name}_USER)

    if(GMX_EXTERNAL_${name})
        if (NOT _library_was_found)
            message(FATAL_ERROR "You have set GMX_EXTERNAL_${name}=ON to instruct GROMACS to use an external ${name} library, but no external library could be detected.")
        endif ()
        # Actually trigger linking.
        list(APPEND LINEAR_ALGEBRA_LIBRARIES ${_libraries_to_link})
    else()
        # Triggering the compilation of the internal version of the library is handled elsewhere.
    endif()
endmacro()

# The default behaviour is to try to detect an "external" BLAS and/or
# LAPACK, perhaps provided by a vendor, use those if found, and
# otherwise fall back on the GROMACS internal implementations of
# these. If the libraries are not in a standard location, the user can
# indicate a search path with CMAKE_PREFIX_PATH.
#
# However, if we are using icc+mkl (so a build command that includes
# -mkl), then it is probably painful to try to link some other BLAS or
# LAPACK. In that case, we use the BLAS & LAPACK provided by MKL. In
# principle, we could offer a more configurable behaviour if/when
# there is need to (say) use vendor BLAS with MKL for FFTs.
#
# If the vendor BLAS and/or LAPACK have abnormal library names, then
# the default searching procedure will fail (e.g. Redmine #771). The
# GMX_(BLAS|LAPACK)_USER variables can be used to indicate the correct
# libraries. If these do not work, a warning is emitted and we try to
# use them anyway, assuming the user knows what they are doing.

# Inputs:
#     GMX_EXTERNAL_BLAS     user input about whether to detect BLAS
#     GMX_EXTERNAL_LAPACK   user input about whether to detect LAPACK
#     HAVE_LIBMKL           true if the build will link to MKL
#     FFT_LINKER_FLAGS      used iff HAVE_MKL
#     FFT_LIBRARIES         used iff HAVE_MKL
#     GMX_BLAS_USER         user input for BLAS libraries to use
#     GMX_LAPACK_USER       user input for LAPACK libraries to use
#
# This function sets the following cache variables:
#     GMX_EXTERNAL_BLAS     according to whether external BLAS is being used
#     GMX_EXTERNAL_LAPACK   according to whether external LAPACK is being used
#     GMX_BLAS_USER         off = use a system library;
#                           any other value = full path to the library to use
#     GMX_LAPACK_USER       off = use a system library;
#                           any other value = full path to the library to use
#
# This function sets the following variables in its parent scope:
#     LINEAR_ALGEBRA_LIBRARIES  will be set as required to add libraries required for linear algebra
#
function(gmxManageLinearAlgebraLibraries)
    include(CheckFunctionExists)
    # Probably not necessary to unset, but let's be clear about usage.
    unset(LINEAR_ALGEBRA_LIBRARIES)

    manage_linear_algebra_library(BLAS dgemm_)
    set(BLAS_FIND_QUIETLY ON)
    manage_linear_algebra_library(LAPACK cheev_)

    # Propagate the new local value to the parent scope
    set(LINEAR_ALGEBRA_LIBRARIES "${LINEAR_ALGEBRA_LIBRARIES}" PARENT_SCOPE)
endfunction()

gmxManageLinearAlgebraLibraries()
