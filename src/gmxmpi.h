#ifdef GMX_LIB_MPI
/* MPI C++ binding is deprecated and can cause name conflicts (e.g. stdio/mpi seek) */
#define MPICH_SKIP_MPICXX 1
#define OMPI_SKIP_MPICXX 1
#include <mpi.h>
#endif
#ifdef GMX_THREAD_MPI
#include "tmpi.h"
#endif
