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
    # Default behaviour is to try to use an external library, but fall
    # back on the internal one if none is found.
    set(GMX_EXTERNAL_${name}_DEFAULT ON)

    # Look for user-specified libraries if external libraries have
    # been specified (or specified by default).
    if(GMX_${name}_USER AND (GMX_EXTERNAL_${name} OR NOT DEFINED GMX_EXTERNAL_${name}))
        set(CMAKE_REQUIRED_LIBRARIES ${GMX_${name}_USER})
        check_function_exists(${function_in_library} _${name}_user_works)
        if(NOT _${name}_user_works)
            message(WARNING "GMX_${name}_USER library ${GMX_${name}_USER} was specified, but it may not provide ${name}. We are proceeding by assuming you know what you are doing and that linking F77-style to this library will work.")
        endif()

        set(_${name}_libraries ${GMX_${name}_USER})
        set(GMX_${name}_FOUND 1)

        if(HAVE_LIBMKL)
            message(STATUS "MKL and GMX_${name}_USER were both specified. Using the latter for ${name}.")
        endif()
    endif()

    if(NOT GMX_${name}_FOUND AND HAVE_LIBMKL)
        set(CMAKE_REQUIRED_LIBRARIES ${FFT_LINKER_FLAGS} ${FFT_LIBRARIES})
        check_function_exists(${function_in_library} _${name}_mkl_works)
        if(_${name}_mkl_works)
            # If we ever need to compile with MKL linear algebra and
            # not with FFT supplied by MKL, uncomment the next line
            # (and probably tweak other things).
#            list(APPEND LINEAR_ALGEBRA_LIBRARIES ${FFT_LINKER_FLAGS} ${FFT_LIBRARIES})
    endif()
            set(GMX_${name}_FOUND 1)
        else()
            set(MESSAGE_TEXT "${MESSAGE_TEXT}Intel's MKL was specified, and it should provide ${name}, but it does not. ")
        endif()
    endif()

    # If detection of ${name} has never run, or none of the preceding
    # detection succeeded, try to detect ${name} in the CMake
    # detection paths, etc.
    if (NOT GMX_${name}_FOUND)
        set(${name}_FIND_QUIETLY ON)
        # Note that this finds all kinds of system libraries,
        # including Apple's Accelerate Framework (and perhaps MKL for
        # icc < 11).
        find_package(${name})
        if (${name}_FOUND)
            set(_${name}_libraries ${${name}_LIBRARIES})
            set(GMX_${name}_FOUND 1)
        endif()
    endif()

    if (NOT GMX_${name}_FOUND)
        set(MESSAGE_TEXT "${MESSAGE_TEXT}A ${name} library was not found by CMake in the paths available to it. ")

        # Couldn't find a linear algebra library. If it was really
        # critical, then we'd issue a fatal error, but it isn't.
        # Instead, we set the default for the variable so that it will
        # trigger the compilation of the GROMACS internal library.
        set(GMX_EXTERNAL_${name}_DEFAULT OFF)

        if(GMX_EXTERNAL_${name})
            # We can only get here if the user set this variable on,
            # and there's no way to satisfy their wishes. We trigger
            # falling back to the default internal libraries, but we
            # do not FORCE that setting in the cache. That means the
            # code will compile, and the user will be warned each
            # iteration of CMake that their setting is not being
            # respected. The down side is that the user will see the
            # value still ON in the cache and might misunderstand, so
            # the following message tries to avoid that.
            #
            # An alternative would be setting a variable of some other
            # name, which would later be used to determine the linker
            # line and/or whether the internal version of the
            # libraries was compiled.
            set(MESSAGE_TEXT "${MESSAGE_TEXT}You have set GMX_EXTERNAL_${name}=ON to try to instruct GROMACS to use an external ${name} library, but something is inconsistent, so you are being ignored. ")
            set(GMX_EXTERNAL_${name} OFF PARENT_SCOPE)
        endif()

        message("${MESSAGE_TEXT}Falling back on the GROMACS internal version of the ${name} library instead. This is fine for normal usage.")
    endif()

    set(GMX_EXTERNAL_${name} ${GMX_EXTERNAL_${name}_DEFAULT} CACHE BOOL "Use a ${name} library that is external to GROMACS if possible (ON), or the internal GROMACS one (OFF)")
    mark_as_advanced(GMX_EXTERNAL_${name})

    if(GMX_EXTERNAL_${name})
        # Actually trigger linking.
        list(APPEND LINEAR_ALGEBRA_LIBRARIES ${_${name}_libraries})
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
#     GMX_BLAS_USER         user input for BLAS libraries to use
#     GMX_LAPACK_USER       user input for LAPACK libraries to use
#
# This function sets the following cache variables:
#     GMX_EXTERNAL_BLAS     according to whether external BLAS is being used
#     GMX_EXTERNAL_LAPACK   according to whether external LAPACK is being used
# Their non-cached versions will normally be the same, except when the user
# specified ON and no library was actually available.
#
# This function sets the following variables in its parent scope:
#     LINEAR_ALGEBRA_LIBRARIES  will be appended as required to add libraries required for linear algebra
#
function(gmxManageBlasAndLapack)
    include(CheckFunctionExists)

    manage_linear_algebra_library(BLAS dgemm)
    manage_linear_algebra_library(LAPACK cheev)

    # Propagate the new local value to the parent scope
    set(LINEAR_ALGEBRA_LIBRARIES "${LINEAR_ALGEBRA_LIBRARIES}" PARENT_SCOPE)
endfunction()

gmxManageBlasAndLapack()
