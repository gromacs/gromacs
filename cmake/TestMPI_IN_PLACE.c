
#include <mpi.h>

int main(void)
{
    void *buf;
    MPI_Allreduce(MPI_IN_PLACE, buf, 10, MPI_FLOAT, MPI_SUM, MPI_COMM_WORLD);
}


