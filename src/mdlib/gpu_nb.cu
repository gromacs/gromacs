#include "stdlib.h"

#include "smalloc.h"

#include "cutypedefs.h"
#include "cudautils.h"

#include "gpu_nb.h"
#include "gpu_data.h"
#include "gpu_nb_kernels.h"


__global__ void __empty_kernel() {}

void cu_do_nb(t_cudata d_data, rvec x[], rvec f[])
{
    int     nb_blocks = (d_data->natoms % NB_DEFAULT_THREADS == 0 ? 
                d_data->natoms/NB_DEFAULT_THREADS : 
                d_data->natoms/NB_DEFAULT_THREADS + 1);
    dim3    dim_block(nb_blocks, 1, 1);
    dim3    dim_grid(NB_DEFAULT_THREADS, 1, 1);

    // printf("");

    /* transfer coordinates */
    upload_cudata(d_data->x, x, d_data->natoms*sizeof(*x));

    /* do the nonbonded calculations */
    k_calc_nb<<<dim_block, dim_grid>>>(d_data->f, d_data->x, d_data->natoms);
    CU_LAUNCH_ERR("k_calc_nb");

    /* get back the forces */
    download_cudata(f, d_data->f, d_data->natoms*sizeof(*f));
}
