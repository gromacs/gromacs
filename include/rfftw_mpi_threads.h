/* This is not part of fftw-2.1.2, but a hack to be able to
 * use several threads on each node with mpi
 * Erik 990726
 */
#ifndef RFFTW_MPI_THREADS_H
#define RFFTW_MPI_THREADS_H

#include <rfftw_mpi.h>
#include <fftw_threads.h>
#include <rfftw.h>

#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */

/***********************************************************************/

extern void rfftwnd_mpi_threads(int nthreads,
				    rfftwnd_mpi_plan p,
				    int n_fields,
				    fftw_real *local_data,
				    fftw_real *work,
				    fftwnd_mpi_output_order
				    output_order);

/***********************************************************************/

#ifdef __cplusplus
} /* extern "C" */
#endif /* __cplusplus */

#endif /* RFFTW_MPI_THREADS_H */
