#ifndef _NBNXN_CUDA_DATA_MGMT_H_
#define _NBNXN_CUDA_DATA_MGMT_H_

#include "typedefs.h" 
#include "types/nbnxn_cuda_types_ext.h"

#ifdef __cplusplus
extern "C" {
#endif

/*! Initializes the data structures related to CUDA nonbonded calculations. */
void cu_init_nonbonded(FILE * /*fplog*/,
                       cu_nonbonded_t * /*p_cu_nb*/,
                       gmx_bool /*bDomDec*/);

/*! Initilizes force-field (FIXME) related data. */
void cu_init_ff_data(FILE * /*fplog*/, 
                     cu_nonbonded_t /*p_cu_nb*/,
                     const interaction_const_t * /*ic*/,
                     const nonbonded_verlet_t * /*nbv*/);

/*! Initilizes pair-list data for GPU, called at every neighbor search step. */
void cu_init_nblist(cu_nonbonded_t /*cu_nb*/,
                    const nbnxn_pairlist_t * /*h_nblist*/,
                    int /*iloc*/);

/*! Initilizes atom-data on the GPU, called at every neighbor search step. */
void cu_init_atomdata(cu_nonbonded_t /*cu_nb*/,
                      const nbnxn_atomdata_t * /*atomdata*/);

/*! Re-genrates the GPU Ewald force table and resets rlist - used during 
 *  cut-off auto-tuning. */
void cu_reset_rlist_ewaldtab(cu_nonbonded_t /*cu_nb*/,
                             const interaction_const_t * /*ic*/);

/*! Uploads shift vector to the GPU if the box is dynamic (otherwise just returns). */
void cu_move_shift_vec(cu_nonbonded_t /*cu_nb*/, 
                       const nbnxn_atomdata_t * /*nbatom*/);

/*! Clears nonbonded force output array on the GPU. */
void cu_clear_nb_f_out(cu_nonbonded_t /*cu_nb*/);

/*! Clears nonbonded shift force output array and energy outputs on the GPU. */
void cu_clear_nb_e_fs_out(cu_nonbonded_t /*cu_nb*/);

/*! Frees all GPU resources used for the nonbonded calculations. */
void cu_free_nbdata(FILE * /*fplog*/, cu_nonbonded_t /*cu_nb*/,
                    gmx_bool /*bDomDec*/);

/*! Synchronizes the stream passed with the atomdata init stop event. */
void cu_synchstream_atomdata(cu_nonbonded_t /*cu_nb*/, int /*enbatATOMS*/);

/*! Returns the nonbonded GPU timings structure or NULL if GPU is not used 
 *  or timing is turned off. 
 */
cu_timings_t * cu_get_gpu_timings(cu_nonbonded_t /*cu_nb*/);

/*! Resets nonbonded GPU timings. */
void cu_reset_gpu_timings(cu_nonbonded_t /*cu_nb*/);

/*! Calculates the minimum size of prozimity lists to improve SM load balance 
    when executing the non-bonded kernels. */
int cu_calc_min_ci_balanced(cu_nonbonded_t /*cu_nb*/);

#ifdef __cplusplus
}
#endif

#endif // _NBNXN_CUDA_DATA_MGMT_H_
