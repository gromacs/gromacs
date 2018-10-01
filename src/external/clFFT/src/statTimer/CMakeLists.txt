# ########################################################################
# Copyright 2013 Advanced Micro Devices, Inc.
# 
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
# 
# http://www.apache.org/licenses/LICENSE-2.0
# 
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
# ########################################################################


# List the names of the files to compile for the external client . . .
set( StatTimer.Source 	statisticalTimer.CPU.cpp
								statisticalTimer.GPU.cpp
								statisticalTimer.extern.cpp
								stdafx.cpp )

# Windows only uses dllmain
if( MSVC )
	set( StatTimer.Source ${StatTimer.Source} dllmain.cpp )
endif( )

set( StatTimer.Headers	statisticalTimer.h
								statisticalTimer.extern.h
								statisticalTimer.CPU.h
								statisticalTimer.GPU.h
								stdafx.h 
								targetver.h 
								../include/clFFT.h )
								
set( StatTimer.Files ${StatTimer.Source} ${StatTimer.Headers} )

# For a rainy day, add pre-compiled header support
#if( MSVC )
#	if (USE_MSVC_PCH)
	
#		set_source_files_properties(LungAnalysisPCH.cxx
#			PROPERTIES
#			COMPILE_FLAGS "/YcLungAnalysisPCH.h"
#			)
#		foreach( src_file ${UPMC_LA_SRCS} )
#			set_source_files_properties(
#				${src_file}
#				PROPERTIES
#				COMPILE_FLAGS "/YuLungAnalysisPCH.h"
#				)
#		endforeach( src_file ${UPMC_LA_SRCS} )
		
#		list(APPEND UPMC_LA_SRCS LungAnalysisPCH.cxx)
#		list(APPEND UPMC_LA_HDRS LungAnalysisPCH.h)

#	endif(USE_MSVC_PCH)
#endif (MSVC)

add_definitions( "/DSTATTIMER_EXPORTS" )

# Include standard OpenCL headers
include_directories( ${OpenCL_INCLUDE_DIRS} ${PROJECT_BINARY_DIR}/include ../include )

add_library( StatTimer ${StatTimer.Files} )
set_target_properties( StatTimer PROPERTIES VERSION ${CLFFT_VERSION} )
set_target_properties( StatTimer PROPERTIES SOVERSION ${CLFFT_SOVERSION} )
set_target_properties( StatTimer PROPERTIES RUNTIME_OUTPUT_DIRECTORY "${PROJECT_BINARY_DIR}/staging" )
target_link_libraries( StatTimer ${OpenCL_LIBRARIES} )

if( UNIX AND NOT APPLE )
	# This library dependency is brought in by the high precision timer available in linux
	target_link_libraries( StatTimer -lrt )
endif( )

# CPack configuration; include the executable into the package
install( TARGETS StatTimer
        EXPORT Library
        RUNTIME DESTINATION bin${SUFFIX_BIN}
        LIBRARY DESTINATION lib${SUFFIX_LIB}
        ARCHIVE DESTINATION lib${SUFFIX_LIB}/import
        )
