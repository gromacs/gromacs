#Create a single man page. Used by CreateManPage.cmake
file(TO_NATIVE_PATH ./${EXENAME} EXEPATH)
message(STATUS ${EXEPATH})
execute_process(COMMAND ${EXEPATH} -quiet -man nroff 
    RESULT_VARIABLE SUCCESS OUTPUT_QUIET ERROR_QUIET)
if(SUCCESS EQUAL 0)
    configure_file(${INFILE} ${OUTFILE})
    file(REMOVE ${INFILE})
endif()
