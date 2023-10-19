file(GLOB cpu_source_files "${CMAKE_SOURCE_DIR}/tools/sz3/*.cpp")

foreach (SRCFILE IN LISTS cpu_source_files)
    get_filename_component(EXE ${SRCFILE} NAME_WLE)
    message(STATUS "Target = ${EXE}")
    add_executable(${EXE} ${EXE}.cpp)
    target_link_libraries(${EXE} PRIVATE SZ3)
    install(PROGRAMS ${CMAKE_CURRENT_BINARY_DIR}/${EXE}$<TARGET_FILE_SUFFIX:${EXE}> TYPE BIN)

endforeach ()

install(FILES testfloat_8_8_128.dat DESTINATION ${CMAKE_INSTALL_DATADIR}/SZ3)
