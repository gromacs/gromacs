#ifndef GPU_NB_H
#define GPU_NB_H

#include "cutypedefs_ext.h"

#ifdef __cplusplus
extern "C" {
#endif

void cu_stream_nb(cu_nonbonded_t /*cu_nb*/, 
                  const gmx_nb_atomdata_t * /*nbdata*/,                  
                  int /*flags*/,
                  gmx_bool /*nonLocal*/,
                  gmx_bool /*sync*/);

void cu_copyback_nb_data(cu_nonbonded_t cu_nb,
                  const gmx_nb_atomdata_t *nbatom,
                  int flags,
                  gmx_bool nonLocal);

void cu_do_nb(cu_nonbonded_t /*cu_nb*/, rvec /*shiftvec*/[]);

gmx_bool cu_checkstat_nb(cu_nonbonded_t /*cu_nb*/, float * /*time*/);

void cu_blockwait_nb(cu_nonbonded_t /*cu_nb*/,
                     int /*flags*/,
                     gmx_bool /*nonLocal*/,
                     float * /*e_lj*/, float * /*e_el*/, rvec * /*fshift*/);

#ifdef __cplusplus
}
#endif

#endif
