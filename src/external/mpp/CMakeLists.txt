project(mmpi CXX)

# States that CMake required version must be >= 2.6
cmake_minimum_required(VERSION 2.6)

# enable ctest support for testing
enable_testing()

# get code root directory (based on current file name path)
get_filename_component( mmpi_dir ${CMAKE_CURRENT_LIST_FILE} PATH )

set ( mmpi_include_dir ${mmpi_dir} )

set(CMAKE_MODULE_PATH ../../../cmake/)

# enable MPI
find_package(MPI REQUIRED)
set(CMAKE_REQUIRED_DEFINITIONS 	${MPI_COMPILE_FLAGS})
set(CMAKE_REQUIRED_INCLUDES 	${MPI_INCLUDE_PATH})
set(CMAKE_REQUIRED_LIBRARIES 	${MPI_LIBRARIES})

# enable C++0x support within gcc (if supported)
if (CMAKE_COMPILER_IS_GNUCXX)

	# add general flags
	add_definitions( -fshow-column )
	add_definitions( -Wall )

	# add flags for debug mode
	set (CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -g -O0")
	# set (CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -g -O3")

	# ENABLE PROFILING
	# add_definitions( -pg )
	# SET(CMAKE_EXE_LINKER_FLAGS -pg)

	# check for -std=c++0x
	include(CheckCXXCompilerFlag)
	check_cxx_compiler_flag( -std=c++11 CXX11_Support )
	if(CXX11_Support)
		set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
	else()
		message( "WARNING: --std=c++11 not supported by your compiler!" )
	endif()
endif()

set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${MPI_COMPILE_FLAGS}")

# avoid multiple import
if (NOT MEMORY_CHECK_SETUP)
	option(CONDUCT_MEMORY_CHECKS "Checks all test cases for memory leaks using valgrind if enabled." OFF)

	# add -all-valgrind target
	add_custom_target(valgrind)

	# define macro for adding tests
	macro ( add_unit_test case_name )
		
		# take value from environment variable
		set(USE_VALGRIND ${CONDUCT_MEMORY_CHECKS})

		# check whether there was a 2nd argument
		if(${ARGC} GREATER 1)
			# use last argument as a valgrind flag
			set(USE_VALGRIND ${ARG2})
		endif(${ARGC} GREATER 1)
		
		# add test case
		if(USE_VALGRIND)
			# no valgrind support in MSVC 
			if(NOT MSVC)
				# add valgrind as a test
				add_test(NAME valgrind_${case_name} 
					COMMAND valgrind
						--leak-check=full
						--show-reachable=no
						--track-fds=yes
						--error-exitcode=1
						#--log-file=${CMAKE_CURRENT_BINARY_DIR}/valgrind.log.${case_name}
						${CMAKE_CURRENT_BINARY_DIR}/ut_${case_name}
					WORKING_DIRECTORY
						${CMAKE_CURRENT_BINARY_DIR}
				)
			endif(NOT MSVC)
		else(USE_VALGRIND)
			# add normal test
			add_test(ut_${case_name} ut_${case_name})

			# + valgrind as a custom target (only of not explicitly prohibited)
			if ((NOT MSVC) AND ((NOT (${ARGC} GREATER 1)) OR (${ARG2})))
				add_custom_target(valgrind_${case_name} 
					COMMAND valgrind
						--leak-check=full
						--show-reachable=no
						--track-fds=yes
						--error-exitcode=1
						#--log-file=${CMAKE_CURRENT_BINARY_DIR}/valgrind.log.${case_name}
						${CMAKE_CURRENT_BINARY_DIR}/ut_${case_name}
					WORKING_DIRECTORY
						${CMAKE_CURRENT_BINARY_DIR}
				)
				add_dependencies(valgrind valgrind_${case_name})
			endif ((NOT MSVC) AND ((NOT (${ARGC} GREATER 1)) OR (${ARG2})))
		endif(USE_VALGRIND)
	endmacro(add_unit_test)
endif (NOT MEMORY_CHECK_SETUP)

# mark as defined
set(MEMORY_CHECK_SETUP OFF CACHE INTERNAL "Flag to avoid multiple setup" PARENT_SCOPE)

add_subdirectory(../gmock-1.6.0 gmock)

file(GLOB test_cases test/*.cc)
foreach ( case_file ${test_cases})
 	get_filename_component( case_name ${case_file} NAME_WE )

	add_executable(ut_${case_name} ${case_file} mpp.cpp test/util/printer.cc)
	include_directories(ut_${case_name} ${mmpi_include_dir})
	include_directories(ut_${case_name} ${MPI_INCLUDE_PATH})
	include_directories(ut_${case_name} ../gmock-1.6.0/gtest/include)
	include_directories(ut_${case_name} ../gmock-1.6.0/include)
	target_link_libraries(ut_${case_name} gmock)
	target_link_libraries(ut_${case_name} ${MPI_LIBRARIES})
	
	add_unit_test(${case_name})
	
endforeach(case_file)
