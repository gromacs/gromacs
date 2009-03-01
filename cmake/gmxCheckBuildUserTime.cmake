# - Check the username performing the build, as well as date and time
#
# GMX_CHECK_BUILD_USER_TIME()
#
# No input arguments used for this macro.
#
# The output variables (strings) created are:
# (These could be 'unknown'/'anonymous' if detection fails)
#
# BUILD_TIME
# BUILD_USER 
# BUILD_MACHINE
#
macro(gmx_check_build_user_time)
    message(STATUS "Setting build user & time")
    if(CMAKE_HOST_UNIX)
        execute_process( COMMAND date     OUTPUT_VARIABLE BUILD_TIME OUTPUT_STRIP_TRAILING_WHITESPACE)
        execute_process( COMMAND whoami   OUTPUT_VARIABLE TMP_USER       OUTPUT_STRIP_TRAILING_WHITESPACE)
      	execute_process( COMMAND hostname OUTPUT_VARIABLE TMP_HOSTNAME   OUTPUT_STRIP_TRAILING_WHITESPACE)
        set(BUILD_USER    "@TMP_USER@\@@TMP_HOSTNAME@ [CMAKE]")
        execute_process( COMMAND uname -srm OUTPUT_VARIABLE BUILD_MACHINE OUTPUT_STRIP_TRAILING_WHITESPACE)
        message(STATUS "Setting build user & time - OK")
    else(CMAKE_HOST_UNIX)
        set(BUILD_TIME    "Unknown date")
        set(BUILD_USER    "Anonymous@unknown [CMAKE]")
        set(BUILD_MACHINE "@CMAKE_HOST_SYSTEM@ @CMAKE_HOST_SYSTEM_PROCESSOR@") 
        message(STATUS "Setting build user & time - not on Unix, using anonymous")
    endif(CMAKE_HOST_UNIX)
endmacro(gmx_check_build_user_time)
