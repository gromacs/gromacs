#ifndef GPU_NB_H
#define GPU_NB_H

#include "cutypedefs_ext.h"

#ifdef __cplusplus
extern "C" {
#endif

void cu_stream_nb(t_cudata /*d_data*/, 
                  const gmx_nb_atomdata_t * /*nbdata*/,                  
                  gmx_bool /*calc_ene*/,
                  gmx_bool /*sync*/);

void cu_do_nb(t_cudata /*d_data*/, rvec /*shiftvec*/[]);

gmx_bool cu_checkstat_nb(t_cudata /*d_data*/, float * /*time*/);

void cu_blockwait_nb(t_cudata /*d_data*/, gmx_bool /*calc_ene*/, 
                     float * /*e_lj*/, float * /*e_el*/);

#ifdef __cplusplus
}
#endif

#endif
