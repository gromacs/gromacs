# This macro attempts to parse the version string of the C compiler in use.
# Currently supported are only compilers that accept "-dumpversion" argument:
# gcc, Intel Compiler (on Linux and Mac OS), Open64, EkoPath.
#
# C_COMPILER_VERSION    - version string of the current C compiler (CMAKE_C_COMPILER)
# CXX_COMPILER_VERSION  - version string of the current C++ compiler (CMAKE_CXX_COMPILER)
#
macro(get_compiler_version)
    if(NOT C_COMPILER_VERSION)
        execute_process(COMMAND ${CMAKE_C_COMPILER} -dumpversion
            RESULT_VARIABLE _cc_dumpversion_res
            OUTPUT_VARIABLE _cc_dumpversion_out
	    ERROR_VARIABLE  _cc_dumpversion_err
            OUTPUT_STRIP_TRAILING_WHITESPACE)

        if (${_cc_dumpversion_res} EQUAL 0)
            SET(C_COMPILER_VERSION ${_cc_dumpversion_out}
                CACHE STRING "C compiler version string" FORCE)
        else ()
            SET(C_COMPILER_VERSION ""
                CACHE STRING "C compiler version string not available" FORCE)
        endif ()
    endif()

    if(NOT CXX_COMPILER_VERSION AND CMAKE_CXX_COMPILER_LOADED)
        execute_process(COMMAND ${CMAKE_CXX_COMPILER} -dumpversion
            RESULT_VARIABLE _cxx_dumpversion_res
            OUTPUT_VARIABLE _cxx_dumpversion_out
	    ERROR_VARIABLE  _cxx_dumpversion_err
            OUTPUT_STRIP_TRAILING_WHITESPACE)

        if (${_cxx_dumpversion_res} EQUAL 0)
            SET(CXX_COMPILER_VERSION ${_cxx_dumpversion_out}
                CACHE STRING "C++ compiler version string" FORCE)
        else ()
            SET(CXX_COMPILER_VERSION ""
                CACHE STRING "C++ compiler version string not available" FORCE)
        endif ()
    endif ()

    if (NOT "${C_COMPILER_VERSION}" STREQUAL "${CXX_COMPILER_VERSION}" AND CMAKE_CXX_COMPILER_LOADED)
        message(WARNING "The version string of the C and C++ compilers does not match!")
    endif ()

    mark_as_advanced(C_COMPILER_VERSION CXX_COMPILER_VERSION)
endmacro()

macro(get_compiler_info LANGUAGE BUILD_COMPILER BUILD_FLAGS)
    execute_process(COMMAND ${CMAKE_${LANGUAGE}_COMPILER} --version
        RESULT_VARIABLE _exec_result
        OUTPUT_VARIABLE _compiler_version
        ERROR_VARIABLE  _compiler_version)
    # Try executing just the compiler command --version failed
    if(_exec_result)
        execute_process(COMMAND ${CMAKE_${LANGUAGE}_COMPILER}
            RESULT_VARIABLE _exec_result
            OUTPUT_VARIABLE _compiler_version
            ERROR_VARIABLE  _compiler_version)
    endif()
    if(NOT "${_compiler_version}" STREQUAL "")
        string(REGEX MATCH "[^\n]+" _compiler_version "${_compiler_version}")
    endif()

    set(${BUILD_COMPILER}
        "${CMAKE_${LANGUAGE}_COMPILER} ${CMAKE_${LANGUAGE}_COMPILER_ID} ${_compiler_version}")
    set(_build_flags "${CMAKE_${LANGUAGE}_FLAGS}")
    if(CMAKE_BUILD_TYPE STREQUAL "Debug")
        set(_build_flags "${_build_flags} ${CMAKE_${LANGUAGE}_FLAGS_DEBUG}")
    elseif(CMAKE_BUILD_TYPE STREQUAL "Release")
        set(_build_flags "${_build_flags} ${CMAKE_${LANGUAGE}_FLAGS_RELEASE}")
    elseif(CMAKE_BUILD_TYPE STREQUAL "MinSizeRel")
        set(_build_flags "${_build_flags} ${CMAKE_${LANGUAGE}_FLAGS_MINSIZEREL}")
    elseif(CMAKE_BUILD_TYPE STREQUAL "RelWithDebInfo")
        set(_build_flags "${_build_flags} ${CMAKE_${LANGUAGE}_FLAGS_RELWITHDEBINFO}")
    endif()
    set(${BUILD_FLAGS} ${_build_flags})
endmacro()
