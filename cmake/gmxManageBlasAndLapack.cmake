# The default behaviour is to try to detect an "external" BLAS and/or
# LAPACK, use those if found, and otherwise fall back on the GROMACS
# internal implementations of these. If the libraries are not in a
# standard location, the user can indicate a search path with
# CMAKE_PREFIX_PATH.
#
# However, if we are using icc+mkl (so a build command that includes
# -mkl), then it is probably painful to try to link some other BLAS or
# LAPACK. In that case, we use the BLAS & LAPACK provided by
# MKL. Likewise, if we are on Apple and have the Accelerate Framework
# available. In principle, we could offer a more configurable
# behaviour if/when there is need to (say) use vendor BLAS with MKL
# for FFTs.
#
# If the vendor BLAS and/or LAPACK have abnormal library names, then
# the default searching procedure will fail (e.g. Redmine #771). The
# GMX_(BLAS|LAPACK)_USER variables can be used to indicate the correct
# libraries. If these do not work, a fatal error will result.

# Inputs:
#     GMX_EXTERNAL_BLAS     user input about whether to detect BLAS
#     GMX_EXTERNAL_LAPACK   user input about whether to detect LAPACK
#     HAVE_LIBMKL           true if the build will link to MKL
#     GMX_BLAS_USER         user input for BLAS libraries to use
#     GMX_LAPACK_USER       user input for LAPACK libraries to use
#
# This function sets the following cache variables (but does not use FORCE):
#     GMX_EXTERNAL_BLAS     according to whether detection took place
#     GMX_EXTERNAL_LAPACK   according to whether detection took place
#
# This function sets the following variables in its parent scope:
#     HAVE_EXTERNAL_BLAS    if an external BLAS will be used
#     HAVE_EXTERNAL_LAPACK  if an external LAPACK will be used
#     GMX_EXTRA_LIBRARIES   will be appended as required
#
function(gmxManageBlasAndLapack)
    include(CheckFunctionExists)

    if(APPLE)
        find_library(ACCELERATE_FRAMEWORK Accelerate)
        list(APPEND GMX_EXTRA_LIBRARIES ${ACCELERATE_FRAMEWORK})
    endif(APPLE)

    macro(manage_linear_algebra_library name function_in_library)

        if(HAVE_LIBMKL OR ACCELERATE_FRAMEWORK)

            # In either case, we have BLAS/LAPACK routines already organized
            # so no need to detect anything at all
            set(HAVE_EXTERNAL_${name} TRUE)

        else()

            # Set up default behaviour of doing detection and assuming
            # that an external library is not detected.
            set(DETECT_EXTERNAL_${name} TRUE)
            set(HAVE_EXTERNAL_${name} FALSE)

            # First, check for user-specified libraries
            if(GMX_${name}_USER)
                set(CMAKE_REQUIRED_LIBRARIES ${GMX_${name}_USER})
                check_function_exists(${function_in_library} ${name}_FOUND)
                if (NOT ${name}_FOUND)
                    message(FATAL_ERROR "GMX_${name}_USER library ${GMX_${name}_USER} was specified, but it does not provide ${name}. Cannot proceed.")
                endif()
                list(APPEND GMX_EXTRA_LIBRARIES ${GMX_${name}_USER})
                set(HAVE_EXTERNAL_${name} TRUE)
                set(DETECT_EXTERNAL_${name} FALSE)
            endif()
            mark_as_advanced(GMX_${name}_USER)
        
            # If the user has set GMX_EXTERNAL_(BLAS|LAPACK) on, then
            # we try to detect an external BLAS/LAPACK. If the user
            # has set these off, then we do not attempt to detect,
            # thus we fall back on the internal versions. In neither
            # case do we override the user's choice about detection
            # (i.e. no FORCE). If detection fails, a status message
            # is issued that the fall-back is in use.
            set(GMX_EXTERNAL_${name} ${DETECT_EXTERNAL_${name}} CACHE BOOL "Try to detect external ${name}, but fall back on built-in")
            mark_as_advanced(GMX_EXTERNAL_${name})

            if(GMX_EXTERNAL_${name})
                set(${name}_FIND_QUIETLY ON)
                find_package(${name})
                if (${name}_FOUND)
                    list(APPEND GMX_EXTRA_LIBRARIES ${${name}_LIBRARIES})
                    set(HAVE_EXTERNAL_${name} TRUE)
                else()
                    MESSAGE(STATUS "Detection of ${name} failed, using the GROMACS internal ${name} library instead.")
                endif()
            endif()

        endif()

        set(HAVE_EXTERNAL_${name} ${HAVE_EXTERNAL_${name}} PARENT_SCOPE)
    endmacro()

    manage_linear_algebra_library(BLAS dgemm)
    manage_linear_algebra_library(LAPACK cheev)

    # Propagate the new local value to the parent scope
    set(GMX_EXTRA_LIBRARIES "${GMX_EXTRA_LIBRARIES}" PARENT_SCOPE)
endfunction()

gmxManageBlasAndLapack()
