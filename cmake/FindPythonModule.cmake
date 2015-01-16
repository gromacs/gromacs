
# Adapted from code posted on cmake-users by Mark Moll
function(find_python_module module)
    string(TOUPPER ${module} module_upper)
    if(NOT PYTHONMODULE_${module_upper})
	if(ARGC GREATER 1 AND ARGV1 STREQUAL "REQUIRED")
	    set(${module}_FIND_REQUIRED TRUE)
	endif()
        if (NOT PYTHON_EXECUTABLE)
            message(STATUS "Cannot find python module ${module} because no python executable is known")
        else()
	    # A module's location is usually a directory, but for binary modules
	    # it's a .so file.
	    execute_process(COMMAND "${PYTHON_EXECUTABLE}" "-c" 
	        "import re, ${module}; print re.compile('/__init__.py.*').sub('',${module}.__file__)"
	        RESULT_VARIABLE _${module}_status 
	        OUTPUT_VARIABLE _${module}_location
	        ERROR_QUIET OUTPUT_STRIP_TRAILING_WHITESPACE)
        endif()
	if(NOT _${module}_status)
	    set(PYTHONMODULE_${module_upper} ${_${module}_location} CACHE STRING 
		"Location of Python module ${module}")
	endif()
    endif()
    find_package_handle_standard_args(PYTHONMODULE_${module} DEFAULT_MSG PYTHONMODULE_${module_upper})
endfunction()
