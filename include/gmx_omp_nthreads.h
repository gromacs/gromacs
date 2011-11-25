#ifndef _GMX_OMP_NTHREADS_
#define _GMX_OMP_NTHREADS_

/*! Initializes the per-module thread count. It is compatible with tMPI, 
 *  thread-safety is ensured (for the features available with tMPI). 
 *  This function should only be caled once during the execution. */
void init_module_nthreads(t_commrec *cr);

/*! Returns the number of threads for domain decomposition. */
int gmx_omp_get_domdec_nthreads();

/*! Returns the number of threads for pair search. */
int gmx_omp_get_pairsearch_nthreads();

/*! Returns the number of threads for non-bonded force calculations. */
int gmx_omp_get_nonbonded_nthreads();

/*! Returns the number of threads for bonded force calculations. */
int gmx_omp_get_bonded_nthreads();

/*! Returns the number of threads for PME. */
int gmx_omp_get_pme_nthreads();

/*! Returns the number of threads for LINCS. */
int gmx_omp_get_lincs_nthreads();

/*! Returns the number of threads for update. */
int gmx_omp_get_update_nthreads();


#endif // _GMX_OMP_NTHREADS_
