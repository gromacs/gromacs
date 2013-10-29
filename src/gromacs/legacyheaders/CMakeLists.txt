# includes: Nothing to build, just installation
file(GLOB ROOT_LEGACY_HEADERS       *.h)
file(GLOB THREAD_MPI_LEGACY_HEADERS thread_mpi/*.h)
file(GLOB THREAD_MPI_ATOMIC_LEGACY_HEADERS thread_mpi/atomic/*.h)
file(GLOB TYPES_LEGACY_HEADERS      types/*.h)
gmx_install_headers(legacyheaders                   ${ROOT_LEGACY_HEADERS})
gmx_install_headers(legacyheaders/thread_mpi        ${THREAD_MPI_LEGACY_HEADERS})
gmx_install_headers(legacyheaders/thread_mpi/atomic ${THREAD_MPI_ATOMIC_LEGACY_HEADERS})
gmx_install_headers(legacyheaders/types             ${TYPES_LEGACY_HEADERS})
