#ifndef GPU_NB_H
#define GPU_NB_H

#include "cutypedefs_ext.h"

#ifdef __cplusplus
extern "C" {
#endif

void cu_stream_nb(cu_nonbonded_t /*cu_nb*/, 
                  const gmx_nb_atomdata_t * /*nbdata*/,                  
                  int /*flags*/,
                  int /*iloc*/);

void cu_copyback_nb_data(cu_nonbonded_t /*cu_nb*/,
                         const gmx_nb_atomdata_t * /*nbatom*/,
                         int /*flags*/,
                         int /*aloc*/);

void cu_do_nb(cu_nonbonded_t /*cu_nb*/, rvec /*shiftvec*/[]);

void cu_blockwait_nb(cu_nonbonded_t /*cu_nb*/,
                     int /*flags*/,
                     int /*aloc*/,
                     float * /*e_lj*/, float * /*e_el*/, rvec * /*fshift*/);

#ifdef __cplusplus
}
#endif

#endif
