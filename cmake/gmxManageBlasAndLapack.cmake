# Helper macro for the function below. We treat BLAS and LAPACK the same
# way, so there's no reason to duplicate the logic for them
#
# Arguments should be:
#     name                 "BLAS" or "LAPACK"
#     function_in_library  the name of a function to use in a linking test of that library
macro(manage_linear_algebra_library name function_in_library)
    # Default behaviour is to try to use an external library, but fall
    # back on the internal one if none is found.
    set(GMX_EXTERNAL_${name}_DEFAULT ON)

    if(GMX_${name}_USER)
        # Check that the user-specified libraries work.

        if(HAVE_LIBMKL)
            message("MKL and GMX_${name}_USER were both specified. Using the latter.")
        endif()
        if(ACCELERATE_FRAMEWORK)
            message("The Apple Accelerate Framework is available and GMX_${name}_USER was specified. Using the latter.")
        endif()

        set(CMAKE_REQUIRED_LIBRARIES ${GMX_${name}_USER})
        check_function_exists(${function_in_library} ${name}_FOUND)
        if (NOT ${name}_FOUND)
            message(FATAL_ERROR "GMX_${name}_USER library ${GMX_${name}_USER} was specified, but it does not provide ${name}. Cannot proceed.")
        endif()

        list(APPEND LINEAR_ALGEBRA_LIBRARIES ${GMX_${name}_USER})
    else()
        # Detect libraries if necessary, then check they work.

        if(HAVE_LIBMKL OR ACCELERATE_FRAMEWORK)
            # In either case, we have BLAS/LAPACK routines already organized
            # so no need to detect anything at all
            if(HAVE_LIBMKL)
                set(LIBRARY_PROVIDING_${name} "Intel's MKL")
            endif()
            if(ACCELERATE_FRAMEWORK)
                set(LIBRARY_PROVIDING_${name} "Apple's Accelerate Framework")
            endif()

            check_function_exists(${function_in_library} ${name}_FOUND)
            if (NOT ${name}_FOUND)
                set(MESSAGE_TEXT "${LIBRARY_PROVIDING_${name}} was specified, and it should provide ${name}, but it does not. ")
            endif()
        endif()

        # If detection of ${name} has never run, or MKL/Accelerate
        # failed, try to detect ${name}.
        if (NOT ${name}_FOUND)
            set(${name}_FIND_QUIETLY ON)
            find_package(${name})
            if (${name}_FOUND)
                list(APPEND LINEAR_ALGEBRA_LIBRARIES ${${name}_LIBRARIES})
            else()
                set(MESSAGE_TEXT "${MESSAGE_TEXT}A ${name} library was not found by CMake in the paths available to it. ")
            endif()
        endif()

        if (NOT ${name}_FOUND AND GMX_EXTERNAL_${name})
            # Couldn't find a linear algebra library. If it was really
            # critical, then we'd issue a fatal error, but it
            # isn't. We set the default for the variable to trigger
            # the compilation of the GROMACS internal library instead.
            set(GMX_EXTERNAL_${name}_DEFAULT OFF)

            if(GMX_EXTERNAL_${name})
                # We can only get here if the user set this variable
                # on, and there's no way to satisfy their wishes. We
                # trigger falling back to the default internal
                # libraries, but we do not FORCE that setting in the
                # cache. That means the code will compile, and the
                # user will be warned each iteration of CMake that
                # their setting is not being respected. The down side
                # is that the user will see the value still ON in the
                # cache and might misunderstand, so the following
                # message tries to avoid that.
                set(MESSAGE_TEXT "${MESSAGE_TEXT}You have set GMX_EXTERNAL_${name}=ON to try to instruct GROMACS to use an external ${name} library, but something is inconsistent, so you are being ignored. ")
                set(GMX_EXTERNAL_${name} OFF)
            endif()

            message("${MESSAGE_TEXT}Falling back on the GROMACS internal version of the ${name} library instead.")
        endif()
    endif()

    set(GMX_EXTERNAL_${name} ${GMX_EXTERNAL_${name}_DEFAULT} CACHE BOOL "Use a ${name} library that is external to GROMACS (ON), or the internal GROMACS one (OFF)")
    mark_as_advanced(GMX_EXTERNAL_${name})
endmacro()

# The default behaviour is to try to detect an "external" BLAS and/or
# LAPACK, perhaps provided by a vendor, use those if found, and
# otherwise fall back on the GROMACS internal implementations of
# these. If the libraries are not in a standard location, the user can
# indicate a search path with CMAKE_PREFIX_PATH.
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
# This function sets the following variables:
#     GMX_EXTERNAL_BLAS     according to whether external BLAS is being used
#     GMX_EXTERNAL_LAPACK   according to whether external LAPACK is being used
# Their cached versions will normally be the same, except when the user
# specified ON and no library was actually available.
#
# This function sets the following variables in its parent scope:
#     LINEAR_ALGEBRA_LIBRARIES  will be appended as required to add libraries required for linear algebra
#
function(gmxManageBlasAndLapack)
    include(CheckFunctionExists)

    if(APPLE)
        find_library(ACCELERATE_FRAMEWORK Accelerate)
        list(APPEND LINEAR_ALGEBRA_LIBRARIES ${ACCELERATE_FRAMEWORK})
    endif(APPLE)

    manage_linear_algebra_library(BLAS dgemm)
    manage_linear_algebra_library(LAPACK cheev)

    # Propagate the new local value to the parent scope
    set(LINEAR_ALGEBRA_LIBRARIES "${LINEAR_ALGEBRA_LIBRARIES}" PARENT_SCOPE)
endfunction()

gmxManageBlasAndLapack()
