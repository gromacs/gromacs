#ifndef _NB_VERLET_
#define _NB_VERLET_

#include "cutypedefs_ext.h"

#ifdef __cplusplus
extern "C" {
#endif

/* nonb-bonded data structure with Verlet-type cut-off */
typedef struct {
    /* The bounding box type neighbor searching data */
    gmx_nbsearch_t nbs;
    int            nnbl;
    gmx_nblist_t   **nbl;
    int            nnbl_nl;
    gmx_nblist_t   **nbl_nl;
    gmx_nb_atomdata_t *nbat;

    gmx_bool  useGPU;   /* use GPU acceleration */
    int kernel_type; /* non-bonded kernel - see enum above */

    /* GPU/CUDA nonbonded data structure */
    cu_nonbonded_t gpu_nb;  /*gpu_data*/
    int min_ci_balanced;    /* balance lists to this size - only used for 
                               the 8x8x8 CUDA kernels */
} nonbonded_verlet_t;

#ifdef __cplusplus
}
#endif

#endif /* _NB_VERLET_ */
