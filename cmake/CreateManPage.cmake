#Create a single man page. Used by CreateManPage.cmake
execute_process(COMMAND ./${EXENAME} -quiet -man nroff 
    RESULT_VARIABLE SUCCESS OUTPUT_QUIET ERROR_QUIET)
if(SUCCESS EQUAL 0)
    configure_file(${INFILE} ${OUTFILE})
    file(REMOVE ${INFILE})
endif()
