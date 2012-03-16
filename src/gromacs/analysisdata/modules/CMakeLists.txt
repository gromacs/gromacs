set(ANALYSISDATA_MODULES_PUBLIC_HEADERS
    average.h
    displacement.h
    histogram.h
    plot.h)
install(FILES ${ANALYSISDATA_MODULES_PUBLIC_HEADERS}
        DESTINATION ${INCL_INSTALL_DIR}/gromacs/analysisdata/modules
        COMPONENT development)
