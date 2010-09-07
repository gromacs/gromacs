#ifndef _GPU_DATA_H_
#define _GPU_DATA_H_

#include "typedefs.h" 
#include "cutypedefs_ext.h"

#ifdef __cplusplus
extern "C" {
#endif

void init_cudata(FILE * /*fplog*/, t_cudata * /*d_pdata*/,
        t_forcerec * /*fr*/, t_mdatoms * /*mdatoms*/,
        gmx_mtop_t * /*top_global*/ /*, gmx_localtop_t * */ /*top*/);

void destroy_cudata(FILE * /*fplog*/, t_cudata  /*d_data*/);

#ifdef __cplusplus
}
#endif

#endif // _GPU_DATA_H_
