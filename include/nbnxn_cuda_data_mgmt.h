#ifndef _NBNXN_CUDA_DATA_MGMT_H_
#define _NBNXN_CUDA_DATA_MGMT_H_

#include "typedefs.h" 
#include "types/nbnxn_cuda_types_ext.h"

#ifdef __cplusplus
extern "C" {
#endif

/*! Initializes the data structures related to CUDA nonbonded calculations. */
void nbnxn_cuda_init(FILE * /*fplog*/,
                     cu_nonbonded_t * /*p_cu_nb*/,
                     gmx_bool /*bDomDec*/);

/*! Initilizes simulation constant data. */
void nbnxn_cuda_init_const(FILE * /*fplog*/,
                           cu_nonbonded_t /*p_cu_nb*/,
                           const interaction_const_t * /*ic*/,
                           const nonbonded_verlet_t * /*nbv*/);

/*! Initilizes pair-list data for GPU, called at every neighbor search step. */
void nbnxn_cuda_init_pairlist(cu_nonbonded_t /*cu_nb*/,
                              const nbnxn_pairlist_t * /*h_nblist*/,
                              int /*iloc*/);

/*! Initilizes atom-data on the GPU, called at every neighbor search step. */
void nbnxn_cuda_init_atomdata(cu_nonbonded_t /*cu_nb*/,
                              const nbnxn_atomdata_t * /*atomdata*/);

/*! Re-generates the GPU Ewald force table and resets rlist - used with PME auto-tuning. */
void reset_gpu_rlist_ewaldtab(cu_nonbonded_t /*cu_nb*/,
                              const interaction_const_t * /*ic*/);

/*! Uploads shift vector to the GPU if the box is dynamic (otherwise just returns). */
void nbnxn_cuda_upload_shiftvec(cu_nonbonded_t /*cu_nb*/,
                                const nbnxn_atomdata_t * /*nbatom*/);

/*! Clears GPU outputs: nonbonded force, shift force and energy. */
void nbnxn_cuda_clear_outputs(cu_nonbonded_t /*cu_nb*/, int /*flags*/);

/*! Frees all GPU resources used for the nonbonded calculations. */
void nbnxn_cuda_free(FILE * /*fplog*/, cu_nonbonded_t /*cu_nb*/,
                     gmx_bool /*bDomDec*/);

/*! Returns the nonbonded GPU timings structure or NULL if GPU is not used 
   or timing is turned off. */
cu_timings_t * nbnxn_cuda_get_timings(cu_nonbonded_t /*cu_nb*/);

/*! Resets nonbonded GPU timings. */
void nbnxn_cuda_reset_timings(cu_nonbonded_t /*cu_nb*/);

/*! Calculates the minimum size of proximity lists to improve SM load balance 
    when executing the non-bonded kernels. */
int nbnxn_cuda_min_ci_balanced(cu_nonbonded_t /*cu_nb*/);

#ifdef __cplusplus
}
#endif

#endif // _NBNXN_CUDA_DATA_MGMT_H_
