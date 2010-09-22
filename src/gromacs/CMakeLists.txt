include_directories(${CMAKE_CURRENT_SOURCE_DIR}/..)

add_subdirectory(options)
add_subdirectory(selection)

file(GLOB GMXLIB_SOURCES */*.cpp */*/*.cpp)
file(GLOB_RECURSE NOT_GMXLIB_SOURCES */tests/*.cpp trajana//*.cpp)
list(REMOVE_ITEM GMXLIB_SOURCES ${NOT_GMXLIB_SOURCES})

add_library(gromacs ${GMXLIB_SOURCES})
target_link_libraries(gromacs gmx ${GMX_EXTRA_LIBRARIES} ${THREAD_LIB})
set_target_properties(gromacs PROPERTIES OUTPUT_NAME "gromacs${GMX_LIBS_SUFFIX}" SOVERSION ${SOVERSION} INSTALL_NAME_DIR "${LIB_INSTALL_DIR}")

install(TARGETS gromacs DESTINATION ${LIB_INSTALL_DIR})
