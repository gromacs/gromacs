#include <stdlib.h>
#include <stdio.h>

#include "gmx_fatal.h"
#include "smalloc.h"

#include "cutypedefs.h"
#include "cudautils.h"
#include "gpu_data.h"

/*** CUDA MD Data operations ***/

void init_cudata(FILE *fplog, t_cudata *dp_data,
        t_forcerec *fr, t_mdatoms *mdatoms,
        gmx_mtop_t *top_global/*, gmx_localtop_t *top*/)
{
    t_cudata    d_data = NULL; 
    cudaError_t stat;
    int         natoms = top_global->natoms;
    int         ntypes = fr->ntype;;

    if (dp_data == NULL) return;
    
    snew(d_data, 1);

    d_data->ntypes  = ntypes;
    d_data->natoms  = natoms;
    
    d_data->eps_r = fr->epsilon_r;
    d_data->eps_rf = fr->epsilon_rf;   

    stat = cudaMalloc((void **)&d_data->f, natoms*sizeof(*(d_data->f)));
    CU_RET_ERR(stat, "cudaMalloc failed on d_data->f");
    stat = cudaMalloc((void **)&d_data->x, natoms*sizeof(*(d_data->x)));
    CU_RET_ERR(stat, "cudaMalloc failed on d_data->x");

    stat = cudaMalloc((void **)&d_data->atom_types, natoms*sizeof(*(d_data->atom_types)));
    CU_RET_ERR(stat, "cudaMalloc failed on d_data->atom_types"); 
    upload_cudata(d_data->atom_types, mdatoms->typeA, natoms*sizeof(*(d_data->atom_types)));

    stat = cudaMalloc((void **)&d_data->charges, ntypes*sizeof(*(d_data->charges)));
    CU_RET_ERR(stat, "cudaMalloc failed on d_data->charges");

    stat = cudaMalloc((void **)&d_data->masses, ntypes*sizeof(*(d_data->masses)));
    CU_RET_ERR(stat, "cudaMalloc failed on d_data->masses");   

    stat = cudaMalloc((void **)&d_data->nbfp, 2*ntypes*sizeof(*(d_data->nbfp)));
    CU_RET_ERR(stat, "cudaMalloc failed on d_data->masses"); 
    upload_cudata(d_data->nbfp, mdatoms->massT, 2*ntypes*sizeof(*(d_data->nbfp)));

    fprintf(fplog, "Initialized CUDA data structures.\n");

    *dp_data = d_data;
}

void destroy_cudata(FILE *fplog, t_cudata d_data)
{
    cudaError_t stat;

    if (d_data == NULL) return;

    stat = cudaFree(d_data->f);
    CU_RET_ERR(stat, "cudaFree failed on d_data->f");
    stat = cudaFree(d_data->x);   
    CU_RET_ERR(stat, "cudaFree failed on d_data->x");
    stat = cudaFree(d_data->atom_types);   
    CU_RET_ERR(stat, "cudaFree failed on d_data->atom_types");
    stat = cudaFree(d_data->charges);
    CU_RET_ERR(stat, "cudaFree failed on d_data->charges");
    stat = cudaFree(d_data->masses);
    CU_RET_ERR(stat, "cudaFree failed on d_data->masses");
    stat = cudaFree(d_data->nbfp);
    CU_RET_ERR(stat, "cudaFree failed on d_data->nbfp");

    fprintf(fplog, "Cleaned up CUDA data structures.\n");

}


/*** Nonbonded related CUDA data operations ***/

void upload_nb_cudata(t_cudata d_data, rvec h_x[])
{
    return;
}

void download_nb_cudata(t_cudata d_data, rvec h_f[])
{
    return;
}
