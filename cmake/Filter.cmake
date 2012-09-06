#Can be used to filter files at build time
#Usage: cmake  -DINFILE=... -DOUTFILE=... [variables to replace] -P Filter.cmake
configure_file(${INFILE} ${OUTFILE})
file(REMOVE ${INFILE})
