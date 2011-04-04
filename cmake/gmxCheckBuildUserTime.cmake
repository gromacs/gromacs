# - Check the username performing the build, as well as date and time
#
# GMX_CHECK_BUILD_USER_TIME(BUILD_TIME BUILD_USER BUILD_MACHINE)
#
# The macro variables will be set to the user/machine used for configuration,
# or anonymous/unknown if it cannot be detected (windows)
#
# BUILD_TIME
# BUILD_USER 
# BUILD_MACHINE
#
macro(gmx_check_build_user_time BUILD_TIME BUILD_USER BUILD_MACHINE)
    IF(NOT DEFINED ${BUILD_USER})

    message(STATUS "Setting build user & time")
    if(CMAKE_HOST_UNIX)
        execute_process( COMMAND date     OUTPUT_VARIABLE TMP_TIME    OUTPUT_STRIP_TRAILING_WHITESPACE)
        execute_process( COMMAND whoami   OUTPUT_VARIABLE TMP_USER       OUTPUT_STRIP_TRAILING_WHITESPACE)
        execute_process( COMMAND hostname OUTPUT_VARIABLE TMP_HOSTNAME   OUTPUT_STRIP_TRAILING_WHITESPACE)
        set(${BUILD_USER}    "@TMP_USER@\@@TMP_HOSTNAME@ [CMAKE]" CACHE INTERNAL "Build user")
        set(${BUILD_TIME}    "@TMP_TIME@" CACHE INTERNAL "Build date & time")
        execute_process( COMMAND uname -srm OUTPUT_VARIABLE TMP_MACHINE OUTPUT_STRIP_TRAILING_WHITESPACE)
        set(${BUILD_MACHINE} "@TMP_MACHINE@" CACHE INTERNAL "Build host & architecture")
        message(STATUS "Setting build user & time - OK")
    else(CMAKE_HOST_UNIX)
        set(${BUILD_USER}    "Anonymous@unknown [CMAKE]" CACHE INTERNAL "Build user")
        set(${BUILD_TIME}    "Unknown date" CACHE INTERNAL "Build date & time")
        set(${BUILD_MACHINE} "@CMAKE_HOST_SYSTEM@ @CMAKE_HOST_SYSTEM_PROCESSOR@" CACHE INTERNAL "Build host & architecture") 
        message(STATUS "Setting build user & time - not on Unix, using anonymous")
    endif(CMAKE_HOST_UNIX)

    ENDIF(NOT DEFINED ${BUILD_USER})
endmacro(gmx_check_build_user_time BUILD_TIME BUILD_USER BUILD_MACHINE)

