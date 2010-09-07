#ifndef GPU_NB_H
#define GPU_NB_H

#include "cutypedefs_ext.h"

#ifdef __cplusplus
extern "C" {
#endif

void cu_do_nb(t_cudata d_data, rvec x[], rvec f[]);

#ifdef __cplusplus
}
#endif

#endif
