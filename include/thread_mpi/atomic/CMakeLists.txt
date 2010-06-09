# includes: Nothing to build, just installation
install(DIRECTORY . DESTINATION include/thread_mpi/atomic

        PATTERN "Makefile*" EXCLUDE
        PATTERN "CMake*" EXCLUDE
        PATTERN "cmake*" EXCLUDE
        PATTERN "*~" EXCLUDE
)


