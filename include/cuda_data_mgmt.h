#ifndef _CUDA_DATA_MGMT_H_
#define _CUDA_DATA_MGMT_H_

#include "typedefs.h" 
#include "cutypedefs_ext.h"

#ifdef __cplusplus
extern "C" {
#endif

void init_cudata_ff(FILE * /*fplog*/, 
                    cu_nonbonded_t * /*p_cu_nb*/,
                    const t_forcerec * /*fr*/);

void init_cudata_nblist(cu_nonbonded_t /*cu_nb*/,
                        const gmx_nblist_t * /*h_nblist*/,
                        gmx_bool /*nonLocal*/,
                        gmx_bool /*doStream*/);

void init_cudata_nblist(cu_nonbonded_t /*cu_nb*/,
                        const gmx_nblist_t * /*h_nblist*/,
                        gmx_bool /*nonLocal*/,
                        gmx_bool /*doStream*/);

void init_cudata_atoms(cu_nonbonded_t /*cu_nb*/,
                       const gmx_nb_atomdata_t * /*atomdata*/,
                       gmx_bool /*doStream*/);

void cu_move_shift_vec(cu_nonbonded_t /*cu_nb*/, 
                       const gmx_nb_atomdata_t * /*nbatom*/);

void cu_clear_nb_f_out(cu_nonbonded_t cu_nb);
void cu_clear_nb_e_fs_out(cu_nonbonded_t cu_nb);

void destroy_cudata(FILE * /*fplog*/, 
                    cu_nonbonded_t /*cu_nb*/);

void cu_blockwait_atomdata(cu_nonbonded_t /*cu_nb*/);
void cu_synchstream_atomdata(cu_nonbonded_t /*cu_nb*/, gmx_bool /*nonLocal*/);

cu_timings_t * get_gpu_timings(cu_nonbonded_t /*cu_nb*/);

void reset_gpu_timings(cu_nonbonded_t /*cu_nb*/);

int cu_upload_X(cu_nonbonded_t /*cu_nb*/, real * /*h_x*/);
int cu_download_F(real * /*h_f*/, cu_nonbonded_t /*cu_nb*/);

#ifdef __cplusplus
}
#endif

#endif // _CUDA_DATA_MGMT_H_
