set(LIBGROMACS_SOURCES)

add_subdirectory(legacyheaders)
add_subdirectory(gmxlib)
add_subdirectory(mdlib)
add_subdirectory(gmxpreprocess)
add_subdirectory(analysisdata)
add_subdirectory(options)
add_subdirectory(selection)
add_subdirectory(trajectoryanalysis)
add_subdirectory(utility)

file(GLOB LIBGROMACS_HEADERS *.h)
install(FILES ${LIBGROMACS_HEADERS} DESTINATION ${INCL_INSTALL_DIR}/gromacs
        COMPONENT development)

# only fiddle with assembly kernels if we're not doing OpenMM build
if(NOT GMX_OPENMM) 
if(GMX_ASM_USEASM_NASM)
  enable_language(ASM_NASM)
  # if NASM is used, we need a special build command for windows...
  FOREACH(SRC ${GMX_SSEKERNEL_ASM_SRC})
    GET_FILENAME_COMPONENT(FILE_BASE ${SRC} NAME_WE)
    SET(OBJ ${CMAKE_CURRENT_BINARY_DIR}/${FILE_BASE}${CMAKE_C_OUTPUT_EXTENSION})

    ADD_CUSTOM_COMMAND(OUTPUT ${OBJ}
                       MAIN_DEPENDENCY ${SRC}
                       COMMAND ${CMAKE_ASM_NASM_COMPILER} -f ${CMAKE_ASM_NASM_OBJECT_FORMAT} -o ${OBJ} ${SRC})

    SET(ALL_ASM_OBJS ${ALL_ASM_OBJS} ${OBJ})
  ENDFOREACH(SRC ${GMX_SSEKERNEL_ASM_SRC})
  set(GMX_SSEKERNEL_ASM_SRC ${ALL_ASM_OBJS})
else(GMX_ASM_USEASM_NASM)
  enable_language(ASM-ATT)
  SET(CMAKE_ASM-ATT_COMPILER ${CMAKE_C_COMPILER})
  if(GMX_IA32_ASM)
    set_source_files_properties(${GMX_SSEKERNEL_ASM_SRC} PROPERTIES COMPILE_FLAGS "-c -m32")
  else()
    set_source_files_properties(${GMX_SSEKERNEL_ASM_SRC} PROPERTIES COMPILE_FLAGS "-c -m64")
  endif()
endif(GMX_ASM_USEASM_NASM)
endif(NOT GMX_OPENMM)

list(APPEND LIBGROMACS_SOURCES ${GMXLIB_SOURCES} ${GMX_SSEKERNEL_ASM_SRC} ${MDLIB_SOURCES})

# add target that generates version.c every time a make is run
# only do this if we generate the version
if (USE_VERSION_H)
    add_custom_target(gmx_version ALL
            COMMAND ${CMAKE_COMMAND} 
                -D GIT_EXECUTABLE="${GIT_EXECUTABLE}"
                -D GIT_VERSION="${GIT_VERSION}"
                -D PROJECT_VERSION="${PROJECT_VERSION}"
                -D PROJECT_SOURCE_DIR="${PROJECT_SOURCE_DIR}"
                -D VERSION_C_CMAKEIN="${CMAKE_CURRENT_SOURCE_DIR}/version.c.cmakein"
                -D VERSION_C_OUT="${CMAKE_CURRENT_BINARY_DIR}/version.c"
                -P ${CMAKE_SOURCE_DIR}/cmake/gmxGenerateVersionInfo.cmake 
            WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}/src/gmxlib 
            DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/version.c.cmakein
            COMMENT "Generating version information")
    list(APPEND LIBGROMACS_SOURCES ${CMAKE_CURRENT_BINARY_DIR}/version.c) # auto-generated
    set_source_files_properties(${CMAKE_CURRENT_BINARY_DIR}/version.c 
                                PROPERTIES GENERATED true)
endif (USE_VERSION_H)

add_library(libgromacs ${LIBGROMACS_SOURCES})
if (USE_VERSION_H)
    add_dependencies(libgromacs gmx_version)
endif (USE_VERSION_H)
target_link_libraries(libgromacs
                      ${GMX_EXTRA_LIBRARIES} ${FFT_LIBRARIES} ${XML_LIBRARIES}
                      ${THREAD_LIB})
set_target_properties(libgromacs PROPERTIES
                      OUTPUT_NAME "gromacs${GMX_LIBS_SUFFIX}"
                      SOVERSION ${SOVERSION}
                      INSTALL_NAME_DIR "${LIB_INSTALL_DIR}")

install(TARGETS libgromacs DESTINATION ${LIB_INSTALL_DIR} COMPONENT libraries)

configure_file(${CMAKE_CURRENT_SOURCE_DIR}/libgromacs.pc.cmakein
               ${CMAKE_CURRENT_BINARY_DIR}/libgromacs.pc @ONLY)
install(FILES ${CMAKE_CURRENT_BINARY_DIR}/libgromacs.pc
        DESTINATION ${LIB_INSTALL_DIR}/pkgconfig
        RENAME "libgromacs${GMX_LIBS_SUFFIX}.pc"
        COMPONENT development)
