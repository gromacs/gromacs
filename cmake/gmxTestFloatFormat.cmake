# - Define macro to determine floating point format properties
#
#  GMX_TEST_FLOAT_FORMAT(FP_IEEE754 FP_BIG_ENDIAN_BYTE FP_BIG_ENDIAN_WORD)
#
#  The thee variables are set to true when:
#  FP_IEEE754          Floating-point numbers are stored in IEEE754 format
#  FP_BIG_ENDIAN_BYTE  The order of bytes in each FP word (4 bytes) is big endian
#  FP_BIG_ENDIAN_WORD  The order of FP words in double precision dwords (2 words) is big endian
#
#  On *most* platforms the two last tests will be the same as the integer endian,
#  big e.g. ARM processors have different byte/word order for floating-point storage,
#  so we figured it is a good idea to test both before relying on the format.

MACRO(GMX_TEST_FLOAT_FORMAT FP_IEEE754 FP_BIG_ENDIAN_BYTE FP_BIG_ENDIAN_WORD)
    IF(NOT DEFINED HAVE_${FP_IEEE754})
        MESSAGE(STATUS "Checking floating point format")

        TRY_COMPILE(HAVE_${FP_IEEE754} "${CMAKE_BINARY_DIR}"    
                    "${CMAKE_SOURCE_DIR}/cmake/TestFloatFormat.c"
                    COPY_FILE "${CMAKE_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/TestFloatFormat.bin")  

        if(HAVE_${FP_IEEE754})

            # dont match first/last letter because of string rounding errors :-)
            FILE(STRINGS "${CMAKE_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/TestFloatFormat.bin"
                 GMX_IEEE754_BB_BW LIMIT_COUNT 1 REGEX "ROMACS")
            FILE(STRINGS "${CMAKE_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/TestFloatFormat.bin"
                 GMX_IEEE754_BB_LW LIMIT_COUNT 1 REGEX "CSXGRO")
            FILE(STRINGS "${CMAKE_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/TestFloatFormat.bin"
                 GMX_IEEE754_LB_BW LIMIT_COUNT 1 REGEX "ORGXSC") 
            FILE(STRINGS "${CMAKE_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/TestFloatFormat.bin"
                 GMX_IEEE754_LB_LW REGEX "SCAMOR")
	    FILE(STRINGS "${CMAKE_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/TestFloatFormat.bin"
                 GMX_IEEE754 REGEX "GROMACS|CSXGRO|ORGXSC|SCAMOR")

            # OS X Universal binaries will contain both strings, set it to the host
            IF(GMX_IEEE754_BB_BW  AND  GMX_IEEE754_LB_LW)
                IF(CMAKE_SYSTEM_PROCESSOR MATCHES powerpc)
                    SET(GMX_IEEE754_BB_BW TRUE)
                    SET(GMX_IEEE754_LB_LW FALSE)
                ELSE(CMAKE_SYSTEM_PROCESSOR MATCHES powerpc)
                    SET(GMX_IEEE754_BB_BW FALSE)
                    SET(GMX_IEEE754_LB_LW TRUE)
                ENDIF(CMAKE_SYSTEM_PROCESSOR MATCHES powerpc)
                MESSAGE(STATUS "GMX_TEST_IEEE754_FORMAT found different results, consider setting CMAKE_OSX_ARCHITECTURES or CMAKE_TRY_COMPILE_OSX_ARCHITECTURES to one or no architecture !")
            ENDIF(GMX_IEEE754_BB_BW  AND  GMX_IEEE754_LB_LW)

            IF(GMX_IEEE754)
                SET(${FP_IEEE754} 1 CACHE INTERNAL "Result of test for IEEE754 FP format" FORCE)
            ENDIF(GMX_IEEE754)

            IF(GMX_IEEE754_BB_BW  OR  GMX_IEEE754_BB_LW)
                SET(${FP_BIG_ENDIAN_BYTE} 1 CACHE INTERNAL "Result of test for big endian FP byte order" FORCE)
            ENDIF(GMX_IEEE754_BB_BW  OR  GMX_IEEE754_BB_LW)

            IF(GMX_IEEE754_BB_BW  OR  GMX_IEEE754_LB_BW)
                SET(${FP_BIG_ENDIAN_WORD} 1 CACHE INTERNAL "Result of test for big endian FP word order" FORCE)
            ENDIF(GMX_IEEE754_BB_BW  OR  GMX_IEEE754_LB_BW)

        endif(HAVE_${FP_IEEE754})

        # just some informational output for the user
        if(GMX_IEEE754_BB_BW)
            MESSAGE(STATUS "Checking floating point format - IEEE754 (BE byte, BE word)")
        elseif(GMX_IEEE754_BB_LW)
            MESSAGE(STATUS "Checking floating point format - IEEE754 (BE byte, LE word)")
        elseif(GMX_IEEE754_LB_BW)
            MESSAGE(STATUS "Checking floating point format - IEEE754 (LE byte, BE word)")
        elseif(GMX_IEEE754_LB_LW)
            MESSAGE(STATUS "Checking floating point format - IEEE754 (LE byte, LE word)")
        else(GMX_IEEE754_LB_LW)
            MESSAGE(STATUS "Checking floating point format - unknown")
        endif(GMX_IEEE754_BB_BW)
    ENDIF(NOT DEFINED HAVE_${FP_IEEE754})
ENDMACRO(GMX_TEST_FLOAT_FORMAT FP_IEEE754 FP_BIG_ENDIAN_BYTE FP_BIG_ENDIAN_WORD)



