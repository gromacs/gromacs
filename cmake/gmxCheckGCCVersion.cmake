# Check GCC version and if any of the 4.1.x family compiler suites is found
# quit the build system generating process. 
#
# The GCC 4.1.x compilers contain an optimization related bug which might 
# results in code that exhibits incorrect behaviour and often leads to 
# exploding systems or crashes. 
#
# For further details see e.g. 
# https://bugs.launchpad.net/ubuntu/+source/gcc-4.1/+bug/158799
#
# Szilard Pall (pszilard@cbr.su.se)
#

if(NOT GMX_DISABLE_GCC41_CHECK)

if(CMAKE_COMPILER_IS_GNUCC)
    # if we have -dumpversion flag use that, otherwise try the --version
    execute_process(COMMAND ${CMAKE_C_COMPILER} -dumpversion
        RESULT_VARIABLE _gcc_dumpversion_res
        OUTPUT_VARIABLE _gcc_dumpversion_out
        OUTPUT_STRIP_TRAILING_WHITESPACE)
    # if gcc returned with error the -dumpversion is not available 
    if(${_gcc_dumpversion_res} EQUAL 0)
        if(${_gcc_dumpversion_out} MATCHES ".*4\\.1\\.[0-9]+.*")
            message(FATAL_ERROR " The GCC compiler in use seems to belong to the 4.1.x 
                family (detected version: ${_gcc_dumpversion_out}). These compilers 
                contain an optimization related bug which might results in code that 
                exhibits incorrect behaviour and often leads to exploding systems or 
                crashes. To disable this check set GMX_DISABLE_GCC41_CHECK=YES.")
        endif()
    else()    
        message(WARNING " The GCC compiler in use does not support the -dumpversion flag. 
            Will attempt parsing the version from the \"gcc --version\" output.")        
        execute_process(COMMAND ${CMAKE_C_COMPILER} --version
            OUTPUT_VARIABLE _gcc_version_out
            OUTPUT_STRIP_TRAILING_WHITESPACE)            
        if("${_gcc_version_out}" MATCHES ".*4\\.1\\.[0-9]+.*")
            message(FATAL_ERROR " The GCC compiler in use seems to belong to the 4.1.x 
                family. These compiler  compilers contain an optimization related bug 
                which might results in code that exhibits incorrect behaviour and 
                often leads to exploding systems or crashes. To disable this check set 
                GMX_DISABLE_GCC41_CHECK=YES.")
        endif()
    endif()
endif()

endif()
