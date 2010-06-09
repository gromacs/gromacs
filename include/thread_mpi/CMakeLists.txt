# includes: Nothing to build, just installation
install(DIRECTORY . DESTINATION include/thread_mpi

        PATTERN "Makefile*" EXCLUDE
        PATTERN "CMake*" EXCLUDE
        PATTERN "cmake*" EXCLUDE
        PATTERN "*~" EXCLUDE
)


