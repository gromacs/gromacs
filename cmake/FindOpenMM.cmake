# Find OpenMM library 
#
# Defines: 
#  OpenMM_FOUND     
#  OpenMM_ROOT_DIR
#  OpenMM_INCLUDE_DIR
#  OpenMM_LIBRARY_DIR
#  OpenMM_LIBRARIES
#  OpenMM_PLUGIN_DIR
#

# include(LibFindMacros) TODO get this: http://zi.fi/cmake/Modules/LibFindMacros.cmake

if(OpenMM_INCLUDE_DIR AND OpenMM_LIBRARY_DIR AND OpenMM_PLUGIN_DIR)
    set(OpenMM_FIND_QUIETLY)
endif()

if(IS_DIRECTORY "$ENV{OPENMM_ROOT_DIR}")
    set(OpenMM_ROOT_DIR $ENV{OPENMM_ROOT_DIR} CACHE PATH "OpenMM installation directory" FORCE)
else()
    if(IS_DIRECTORY "/usr/local/openmm")
        set(OpenMM_ROOT_DIR "/usr/local/openmm" CACHE PATH "OpenMM installation directory" FORCE)
    endif()
endif()

find_library(OpenMM_LIBRARIES
    NAMES OpenMM
    PATHS "${OpenMM_ROOT_DIR}/lib")

get_filename_component(OpenMM_LIBRARY_DIR 
    ${OpenMM_LIBRARIES} 
    PATH
    CACHE STRING "OpenMM library path")

find_path(OpenMM_INCLUDE_DIR 
    NAMES OpenMM.h 
    PATHS "${OpenMM_ROOT_DIR}/include" "${OpenMM_LIBRARY_DIR}/../include"
    CACHE STRING "OpenMM include directory")    

# if we did not manage to set the root dir at the beginning but found the 
# libs then set the ${OpenMM_LIBRARY_DIR}/.. as root
if(NOT IS_DIRECTORY ${OpenMM_ROOT_DIR})
    if (IS_DIRECTORY "${OpenMM_LIBRARY_DIR}/..") # just double-checking
        get_filename_component(OpenMM_ROOT_DIR 
            "${OpenMM_LIBRARY_DIR}/.." 
            ABSOLUTE
            CACHE PATH "OpenMM installation directory")
    endif()
endif()

if(NOT IS_DIRECTORY ${OpenMM_ROOT_DIR})
    message(FATAL_ERROR "Can't find OpenMM! Either set the OPENMM_ROOT_DIR environment " 
    "variable to the path where your OpenMM installation is located or install it to the default location (/usr/local/openmm)!")
endif()

if(NOT OpenMM_LIBRARY_DIR)
    message(FATAL_ERROR "Can't find OpenMM libraries. Check your OpenMM installation!")
endif()

# now we can be sure that we have the library dir
if(IS_DIRECTORY "${OpenMM_LIBRARY_DIR}/plugins")
    get_filename_component(OpenMM_PLUGIN_DIR
        "${OpenMM_LIBRARY_DIR}/plugins"
        ABSOLUTE
        CACHE PATH "OpenMM plugins directory")
else()
    message(WARNING "Could not detect the OpenMM plugin directory at the default location (${OpenMM_LIBRARY_DIR}/plugins)."
            "Check you OpenMM installation or manually set the OPENMM_PLUGIN_DIR environment variable!")
endif()

if(NOT OpenMM_INCLUDE_DIR)
    message(FATAL_ERROR "Can't find OpenMM includes. Check your OpenMM installation!")
endif()

include (FindPackageHandleStandardArgs)
find_package_handle_standard_args(OpenMM DEFAULT_MSG 
                                    OpenMM_ROOT_DIR
                                    OpenMM_LIBRARIES 
                                    OpenMM_LIBRARY_DIR 
                                    OpenMM_INCLUDE_DIR
                                    OpenMM_PLUGIN_DIR)
