
#define TMPI_CHECK_ATOMICS
#include "thread_mpi/atomic.h"

int main(void)
{
    int i;
    void *ptr;
    tMPI_Atomic_t some_atomic;
    tMPI_Atomic_ptr_t *some_atomic_ptr = NULL;
    tMPI_Spinlock_t some_spinlock;

    /* Make the compiler actually emit code for these functions, so
       that things like inability to emit inline assembly get
       tested. It is not expected that the code below can run. */
    tMPI_Atomic_memory_barrier();
    tMPI_Atomic_memory_barrier_acq();
    tMPI_Atomic_memory_barrier_rel();
    i = tMPI_Atomic_get(&some_atomic);
    tMPI_Atomic_set(&some_atomic, 0);
    ptr = tMPI_Atomic_ptr_get(some_atomic_ptr);
    tMPI_Atomic_ptr_set(some_atomic_ptr, ptr);
    tMPI_Atomic_add_return(&some_atomic, 0);
    tMPI_Atomic_fetch_add(&some_atomic, 0);
    tMPI_Atomic_cas(&some_atomic, 0, 1);
    tMPI_Atomic_ptr_cas(some_atomic_ptr, ptr, ptr);
    tMPI_Atomic_swap(&some_atomic, 0);
    tMPI_Atomic_ptr_swap(some_atomic_ptr, ptr);
    tMPI_Spinlock_init(&some_spinlock);
    tMPI_Spinlock_lock(&some_spinlock);
    tMPI_Spinlock_trylock(&some_spinlock);
    tMPI_Spinlock_unlock(&some_spinlock);
    tMPI_Spinlock_islocked(&some_spinlock);
    tMPI_Spinlock_wait(&some_spinlock);
return 0;
}
