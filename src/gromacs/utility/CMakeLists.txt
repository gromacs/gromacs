set(UTILITY_PUBLIC_HEADERS
    flags.h)
install(FILES ${UTILITY_PUBLIC_HEADERS}
        DESTINATION ${INCL_INSTALL_DIR}/gromacs/utility
        COMPONENT development)
