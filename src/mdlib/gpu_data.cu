#include <stdlib.h>
#include <stdio.h>

#include "gmx_fatal.h"
#include "smalloc.h"

#include "cutypedefs.h"
#include "cudautils.h"
#include "gpu_data.h"
#include "cupmalloc.h"

#define USE_CUDA_EVENT_BLOCKING_SYNC FALSE /* makes the CPU thread busy-wait! */
#define EWALD_COULOMB_FORCE_TABLE_SIZE 1536   /* size chosen such we do not run out of texture cache */

#define MY_PI               (3.1415926535897932384626433832795f)
#define TWO_OVER_SQRT_PI    (2.0f/sqrt(MY_PI))
    
#define TIME_GPU_TRANSFERS 1

#define NUM_NB_KERNELS 12

static const char * const nb_k1_names[NUM_NB_KERNELS] = 
{
    "_Z22k_calc_nb_RF_forces_1nPK12gmx_nbl_ci_tPK13gmx_nbl_sj4_tPK14gmx_nbl_excl_tPKiiPK6float4PKfPK6float3fffPSA_",
    "_Z25k_calc_nb_ewald_forces_1nPK12gmx_nbl_ci_tPK13gmx_nbl_sj4_tPK14gmx_nbl_excl_tPKiiPK6float4PKfPK6float3fffPSA_",
    "_Z26k_calc_nb_cutoff_forces_1nPK12gmx_nbl_ci_tPK13gmx_nbl_sj4_tPK14gmx_nbl_excl_tPKiiPK6float4PKfPK6float3fffPSA_",
    "_Z31k_calc_nb_RF_forces_energies_1nPK12gmx_nbl_ci_tPK13gmx_nbl_sj4_tPK14gmx_nbl_excl_tPKiiPK6float4PKfPK6float3fffffPfSI_PSA_",
    "_Z31k_calc_nb_RF_forces_prunenbl_1nPK12gmx_nbl_ci_tP13gmx_nbl_sj4_tPK14gmx_nbl_excl_tPKiiPK6float4PKfPK6float3ffffPS9_",
    "_Z34k_calc_nb_ewald_forces_energies_1nPK12gmx_nbl_ci_tPK13gmx_nbl_sj4_tPK14gmx_nbl_excl_tPKiiPK6float4PKfPK6float3fffffPfSI_PSA_",
    "_Z34k_calc_nb_ewald_forces_prunenbl_1nPK12gmx_nbl_ci_tP13gmx_nbl_sj4_tPK14gmx_nbl_excl_tPKiiPK6float4PKfPK6float3ffffPS9_",
    "_Z35k_calc_nb_cutoff_forces_energies_1nPK12gmx_nbl_ci_tPK13gmx_nbl_sj4_tPK14gmx_nbl_excl_tPKiiPK6float4PKfPK6float3fffffPfSI_PSA_",
    "_Z35k_calc_nb_cutoff_forces_prunenbl_1nPK12gmx_nbl_ci_tP13gmx_nbl_sj4_tPK14gmx_nbl_excl_tPKiiPK6float4PKfPK6float3ffffPS9_",
    "_Z40k_calc_nb_RF_forces_energies_prunenbl_1nPK12gmx_nbl_ci_tP13gmx_nbl_sj4_tPK14gmx_nbl_excl_tPKiiPK6float4PKfPK6float3ffffffPfSH_PS9_",
    "_Z43k_calc_nb_ewald_forces_energies_prunenbl_1nPK12gmx_nbl_ci_tP13gmx_nbl_sj4_tPK14gmx_nbl_excl_tPKiiPK6float4PKfPK6float3ffffffPfSH_PS9_",
    "_Z44k_calc_nb_cutoff_forces_energies_prunenbl_1nPK12gmx_nbl_ci_tP13gmx_nbl_sj4_tPK14gmx_nbl_excl_tPKiiPK6float4PKfPK6float3ffffffPfSH_PS9_"

};

static const char * const nb_k2_names[NUM_NB_KERNELS] = 
{
    "_Z22k_calc_nb_RF_forces_2nPK12gmx_nbl_ci_tPK13gmx_nbl_sj4_tPK14gmx_nbl_excl_tPKiiPK6float4PKfPK6float3fffPSA_",
    "_Z25k_calc_nb_ewald_forces_2nPK12gmx_nbl_ci_tPK13gmx_nbl_sj4_tPK14gmx_nbl_excl_tPKiiPK6float4PKfPK6float3fffPSA_",
    "_Z26k_calc_nb_cutoff_forces_2nPK12gmx_nbl_ci_tPK13gmx_nbl_sj4_tPK14gmx_nbl_excl_tPKiiPK6float4PKfPK6float3fffPSA_",
    "_Z31k_calc_nb_RF_forces_energies_2nPK12gmx_nbl_ci_tPK13gmx_nbl_sj4_tPK14gmx_nbl_excl_tPKiiPK6float4PKfPK6float3fffffPfSI_PSA_",
    "_Z31k_calc_nb_RF_forces_prunenbl_2nPK12gmx_nbl_ci_tP13gmx_nbl_sj4_tPK14gmx_nbl_excl_tPKiiPK6float4PKfPK6float3ffffPS9_",
    "_Z34k_calc_nb_ewald_forces_energies_2nPK12gmx_nbl_ci_tPK13gmx_nbl_sj4_tPK14gmx_nbl_excl_tPKiiPK6float4PKfPK6float3fffffPfSI_PSA_",
    "_Z34k_calc_nb_ewald_forces_prunenbl_2nPK12gmx_nbl_ci_tP13gmx_nbl_sj4_tPK14gmx_nbl_excl_tPKiiPK6float4PKfPK6float3ffffPS9_",
    "_Z35k_calc_nb_cutoff_forces_energies_2nPK12gmx_nbl_ci_tPK13gmx_nbl_sj4_tPK14gmx_nbl_excl_tPKiiPK6float4PKfPK6float3fffffPfSI_PSA_",
    "_Z35k_calc_nb_cutoff_forces_prunenbl_2nPK12gmx_nbl_ci_tP13gmx_nbl_sj4_tPK14gmx_nbl_excl_tPKiiPK6float4PKfPK6float3ffffPS9_",
    "_Z40k_calc_nb_RF_forces_energies_prunenbl_2nPK12gmx_nbl_ci_tP13gmx_nbl_sj4_tPK14gmx_nbl_excl_tPKiiPK6float4PKfPK6float3ffffffPfSH_PS9_",
    "_Z43k_calc_nb_ewald_forces_energies_prunenbl_2nPK12gmx_nbl_ci_tP13gmx_nbl_sj4_tPK14gmx_nbl_excl_tPKiiPK6float4PKfPK6float3ffffffPfSH_PS9_",
    "_Z44k_calc_nb_cutoff_forces_energies_prunenbl_2nPK12gmx_nbl_ci_tP13gmx_nbl_sj4_tPK14gmx_nbl_excl_tPKiiPK6float4PKfPK6float3ffffffPfSH_PS9_"
};

__device__ __global__ void k_empty(){}

/*** CUDA Data operations ***/

static void destroy_cudata_array(void * d_ptr, 
                                 int * n = NULL, int * nalloc = NULL);
static void realloc_cudata_array(void **d_dest, void *h_src, size_t type_size, 
                                 int *curr_size, int *curr_alloc_size, 
                                 int req_size, gmx_bool doStream);                                
static void tabulate_ewald_coulomb_force_r(t_cudata d_data);

/* 
  Tabulates the Ewald Coulomb force.
  Original idea: OpenMM 

 TODO 
    - replace smalloc with pmalloc
    - use double instead of float 
 */
static void tabulate_ewald_coulomb_force_r(t_cudata d_data)
{
    float       *ftmp;
    float       beta, r, x;
    int         i, tabsize;
    cudaError_t stat;
    
    cudaChannelFormatDesc   cd;
    const textureReference  *tex_coulomb_tab;

    beta        = d_data->ewald_beta;
    tabsize     = EWALD_COULOMB_FORCE_TABLE_SIZE;

    d_data->coulomb_tab_size   = tabsize;
    d_data->coulomb_tab_scale = (tabsize - 1) / sqrt(d_data->cutoff_sq);

    smalloc(ftmp, tabsize * sizeof(*ftmp)); 

    for (i = 1; i < tabsize; i++)
    {
        r       = i / d_data->coulomb_tab_scale;
        x       = r * beta;
        ftmp[i] = ((float) erfc(x) / r + beta * TWO_OVER_SQRT_PI * exp(-x * x)) / (r * r);
    }

    ftmp[0] = ftmp[1];

    stat = cudaMalloc((void **)&d_data->coulomb_tab, tabsize * sizeof(*d_data->coulomb_tab));
    CU_RET_ERR(stat, "cudaMalloc failed on d_data->coulomb_tab"); 
    upload_cudata(d_data->coulomb_tab, ftmp, tabsize * sizeof(*d_data->coulomb_tab));

    stat = cudaGetTextureReference(&tex_coulomb_tab, "tex_coulomb_tab");
    CU_RET_ERR(stat, "cudaGetTextureReference on tex_coulomb_tab failed");
    cd = cudaCreateChannelDesc<float>();
    stat = cudaBindTexture(NULL, tex_coulomb_tab, d_data->coulomb_tab, &cd, tabsize * sizeof(*d_data->coulomb_tab));
    CU_RET_ERR(stat, "cudaBindTexture on tex_coulomb_tab failed");

    sfree(ftmp);
}

/*
 Initilizes force-field related data (called only once, inthe beginning).
 */
void init_cudata_ff(FILE *fplog, 
                    t_cudata *dp_data,
                    const t_forcerec *fr)
{
    t_cudata            d_data = NULL;    
    cudaError_t         stat;
    gmx_nb_atomdata_t   *nbat;
    int                 ntypes, i, j;

    nbat = fr->nbat;
    ntypes = nbat->ntype;

    cudaChannelFormatDesc   cd;
    const textureReference  *tex_nbfp;

    int eventflags = ( USE_CUDA_EVENT_BLOCKING_SYNC ? cudaEventBlockingSync: cudaEventDefault );

    if (dp_data == NULL) return;
    
    snew(d_data, 1);

    d_data->time_transfers = TIME_GPU_TRANSFERS > 0; /* TODO fix this! */

    d_data->ewald_beta  = fr->ewaldcoeff;
    d_data->eps_r       = fr->epsilon_r;
    d_data->two_k_rf    = 2.0 * fr->k_rf;
    d_data->c_rf        = fr->c_rf;
    d_data->cutoff_sq   = fr->rcut_nsbox * fr->rcut_nsbox;
    d_data->rlist_sq    = fr->rlist_nsbox * fr->rlist_nsbox;
    
    if (fr->eeltype == eelCUT)
    {
        d_data->eeltype = cu_eelCUT;
    }
    else if (EEL_RF(fr->eeltype))
    {                
        d_data->eeltype = cu_eelRF;
    }
    else if ((EEL_PME(fr->eeltype) || fr->eeltype==eelEWALD))
    {
        d_data->eeltype = cu_eelEWALD;
    }
   else 
    {
        gmx_fatal(FARGS, "The requested electrostatics type is not implemented in the CUDA GPU accelerated kernels!");
    }

    if (d_data->eeltype == cu_eelEWALD)
    {
        tabulate_ewald_coulomb_force_r(d_data);
    }

    /* events for NB async ops */
    d_data->streamGPU = fr->streamGPU;    
    stat = cudaEventCreateWithFlags(&(d_data->start_nb), eventflags);
    CU_RET_ERR(stat, "cudaEventCreate on start_nb failed");
    stat = cudaEventCreateWithFlags(&(d_data->stop_nb), eventflags);
    CU_RET_ERR(stat, "cudaEventCreate on stop_nb failed");
    stat = cudaEventCreateWithFlags(&(d_data->start_atdat), eventflags);
    CU_RET_ERR(stat, "cudaEventCreate on start_atdat failed");
    stat = cudaEventCreateWithFlags(&(d_data->stop_atdat), eventflags);
    CU_RET_ERR(stat, "cudaEventCreate on stop_atdat failed");

    if (d_data->time_transfers)
    {
        stat = cudaEventCreateWithFlags(&(d_data->start_nb_h2d), eventflags);
        CU_RET_ERR(stat, "cudaEventCreate on start_nb_h2d failed");
        stat = cudaEventCreateWithFlags(&(d_data->stop_nb_h2d), eventflags);
        CU_RET_ERR(stat, "cudaEventCreate on stop_nb_h2d failed");

        stat = cudaEventCreateWithFlags(&(d_data->start_nb_d2h), eventflags);
        CU_RET_ERR(stat, "cudaEventCreate on start_nb_d2h failed");
        stat = cudaEventCreateWithFlags(&(d_data->stop_nb_d2h), eventflags);
        CU_RET_ERR(stat, "cudaEventCreate on stop_nb_d2h failed");
    }

    /* NB params */
    d_data->ntypes  = ntypes;
    stat = cudaMalloc((void **)&d_data->nbfp, 2*ntypes*ntypes*sizeof(*(d_data->nbfp)));
    CU_RET_ERR(stat, "cudaMalloc failed on d_data->nbfp"); 
    upload_cudata(d_data->nbfp, nbat->nbfp, 2*ntypes*ntypes*sizeof(*(d_data->nbfp)));

    stat = cudaGetTextureReference(&tex_nbfp, "tex_nbfp");
    CU_RET_ERR(stat, "cudaGetTextureReference on tex_nbfp failed");
    cd = cudaCreateChannelDesc<float>();
    stat = cudaBindTexture(NULL, tex_nbfp, d_data->nbfp, &cd, 2*ntypes*ntypes*sizeof(*(d_data->nbfp)));
    CU_RET_ERR(stat, "cudaBindTexture on tex_nbfp failed");

    stat = cudaMalloc((void**)&d_data->shift_vec, SHIFTS*sizeof(*d_data->shift_vec));
    CU_RET_ERR(stat, "cudaMalloc failed on d_data->shift_vec"); 
    d_data->shift_vec_copied = FALSE;

    stat = cudaMalloc((void**)&d_data->e_lj, sizeof(*d_data->e_lj));
    CU_RET_ERR(stat, "cudaMalloc failed on d_data->e_lj");
    stat = cudaMalloc((void**)&d_data->e_el, sizeof(*d_data->e_el));
    CU_RET_ERR(stat, "cudaMalloc failed on d_data->e_el");

    pmalloc((void**)&d_data->tmpdata.e_lj, sizeof(*d_data->tmpdata.e_lj));
    pmalloc((void**)&d_data->tmpdata.e_el, sizeof(*d_data->tmpdata.e_el));

    /* initilize timing structure */
    d_data->timings.nb_h2d_time = 0.0;
    d_data->timings.nb_d2h_time = 0.0;
    d_data->timings.nb_count = 0;
    d_data->timings.atomdt_h2d_total_time = 0.0;
    d_data->timings.atomdt_count = 0;
    for (i = 0; i < 2; i++)
    {
        for(j = 0; j < 2; j++)
        {
            d_data->timings.k_time[i][j].t = 0.0;
            d_data->timings.k_time[i][j].c = 0;
        }
    }

    /* initilize to NULL all data structures that might need reallocation 
       in init_cudata_atoms */
    d_data->xq      = NULL;
    d_data->f       = NULL;
    d_data->ci      = NULL;
    d_data->sj4     = NULL;
    d_data->excl    = NULL;

    /* size -1 just means that it has not been initialized yet */
    d_data->naps        = -1;
    d_data->natoms      = -1;
    d_data->nalloc      = -1;
    d_data->nci         = -1;
    d_data->ci_nalloc   = -1;
    d_data->nsj4        = -1;
    d_data->sj4_nalloc  = -1;
    d_data->nexcl       = -1;
    d_data->excl_nalloc = -1;

    d_data->prune_nbl = FALSE;

    *dp_data = d_data;

    if (fplog != NULL)
    {
        fprintf(fplog, "Initialized CUDA data structures.\n");
        
        printf("Initialized CUDA data structures.\n");
        fflush(stdout);
    }

    /* k_calc_nb_*_1 48/16 kB Shared/L1 */

    for (int i = 0; i < NUM_NB_KERNELS; i++)
    {
        stat = cudaFuncSetCacheConfig(nb_k1_names[i],  cudaFuncCachePreferShared);
        CU_RET_ERR(stat, "cudaFuncSetCacheConfig failed");
    }

    /* k_calc_nb_*_2 16/48 kB Shared/L1 */
    for (int i = 0; i < NUM_NB_KERNELS; i++)
    {
        stat = cudaFuncSetCacheConfig(nb_k2_names[i], cudaFuncCachePreferL1);
        CU_RET_ERR(stat, "cudaFuncSetCacheConfig failed");
    }

    k_empty<<<1, 512>>>();
}

/*
  Initilizes atom-data for the GPU, called at every neighbor search step. 
*/
void init_cudata_atoms(t_cudata d_data, 
                       const gmx_nb_atomdata_t *nbat, 
                       const gmx_nblist_t *nblist,
                       gmx_bool doStream)
{
    cudaError_t stat;
    char        sbuf[200];
    int         nalloc;
    int         natoms  = nbat->natoms;
   
    /* time async copy */
    stat = cudaEventRecord(d_data->start_atdat, 0);
    CU_RET_ERR(stat, "cudaEventRecord failed on d_data->start_atdat");

    if (d_data->naps < 0)
    {
        d_data->naps = nblist->naps;
    }
    else
    {
        if (d_data->naps != nblist->naps)
        {
            sprintf(sbuf, "In init_cudata_atoms: the #atoms per cell has changed (from %d to %d)",
                    d_data->naps, nblist->naps);            
            gmx_incons(sbuf);
        }
    }

    /* need to reallocate if we have to copy more atoms than the amount of space
       available and only allocate if we haven"t initilzed yet, i.e d_data->natoms == -1 */
    if (natoms > d_data->nalloc)
    {
        nalloc = natoms * 1.2 + 100;
    
        /* free up first if the arrays have already been initialized */
        if (d_data->nalloc != -1)
        {
            destroy_cudata_array(d_data->f, &d_data->natoms, &d_data->nalloc);
            destroy_cudata_array(d_data->xq);
            destroy_cudata_array(d_data->atom_types, &d_data->ntypes);             
        }
        
        stat = cudaMalloc((void **)&d_data->f, nalloc*sizeof(*(d_data->f)));
        CU_RET_ERR(stat, "cudaMalloc failed on d_data->f");                   
        stat = cudaMalloc((void **)&d_data->xq, nalloc*sizeof(*(d_data->xq)));
        CU_RET_ERR(stat, "cudaMalloc failed on d_data->xq");     

        stat = cudaMalloc((void **)&d_data->atom_types, nalloc*sizeof(*(d_data->atom_types)));
        CU_RET_ERR(stat, "cudaMalloc failed on d_data->atom_types"); 

        d_data->nalloc = nalloc;
    }
    /* XXX for the moment we just set all 8 values to the same value... 
       ATM not, we"ll do that later */    
    d_data->natoms = natoms;

    if(doStream)
    {
        upload_cudata_async(d_data->atom_types, nbat->type, natoms * sizeof(*d_data->atom_types), 0);      
    }
    else 
    {
        upload_cudata(d_data->atom_types, nbat->type, natoms * sizeof(*(d_data->atom_types)));
    
    }

    realloc_cudata_array((void **)&d_data->ci, nblist->ci, sizeof(*(d_data->ci)),
                         &d_data->nci, &d_data->ci_nalloc,
                         nblist->nci, doStream);

    realloc_cudata_array((void **)&d_data->sj4, nblist->sj4, sizeof(*(d_data->sj4)),
                         &d_data->nsj4, &d_data->sj4_nalloc,
                         nblist->nsj4, doStream);

    realloc_cudata_array((void **)&d_data->excl, nblist->excl, sizeof(*(d_data->excl)),
                         &d_data->nexcl, &d_data->excl_nalloc,
                         nblist->nexcl, doStream);

    stat = cudaEventRecord(d_data->stop_atdat, 0);
    CU_RET_ERR(stat, "cudaEventRecord failed on d_data->stop_atdat");

    d_data->timings.atomdt_count++;    

    /* pruning of the neighbor list needs to be done when a new list get uploaded */
    d_data->prune_nbl = TRUE;
}

void destroy_cudata(FILE *fplog, t_cudata d_data)
{
    cudaError_t stat;
    const textureReference  *tex;

    if (d_data == NULL) return;

    if (d_data->eeltype == cu_eelEWALD)
    {
        stat = cudaGetTextureReference(&tex, "tex_coulomb_tab");
        CU_RET_ERR(stat, "cudaGetTextureReference on tex_coulomb_tab failed");
        stat = cudaUnbindTexture(tex);
        CU_RET_ERR(stat, "cudaUnbindTexture failed on tex");
        destroy_cudata_array(d_data->coulomb_tab, &d_data->coulomb_tab_size);            
    }


    stat = cudaEventDestroy(d_data->start_nb);
    CU_RET_ERR(stat, "cudaEventDestroy failed on d_data->start_nb");
    stat = cudaEventDestroy(d_data->stop_nb);
    CU_RET_ERR(stat, "cudaEventDestroy failed on d_data->stop_nb");
    stat = cudaEventDestroy(d_data->start_atdat);
    CU_RET_ERR(stat, "cudaEventDestroy failed on d_data->start_atdat");
    stat = cudaEventDestroy(d_data->stop_atdat);
    CU_RET_ERR(stat, "cudaEventDestroy failed on d_data->stop_atdat");

    if (d_data->time_transfers)
    {
        stat = cudaEventDestroy(d_data->start_nb_h2d);
        CU_RET_ERR(stat, "cudaEventDestroy failed on d_data->start_nb_h2d");
        stat = cudaEventDestroy(d_data->stop_nb_h2d);
        CU_RET_ERR(stat, "cudaEventDestroy failed on d_data->stop_nb_h2d");

        stat = cudaEventDestroy(d_data->start_nb_d2h);
        CU_RET_ERR(stat, "cudaEventDestroy failed on d_data->start_nb_d2h");
        stat = cudaEventDestroy(d_data->stop_nb_d2h);
        CU_RET_ERR(stat, "cudaEventDestroy failed on d_data->stop_nb_d2h");
    }

    stat = cudaGetTextureReference(&tex, "tex_nbfp");
    CU_RET_ERR(stat, "cudaGetTextureReference on tex_nbfp failed");
    stat = cudaUnbindTexture(tex);
    CU_RET_ERR(stat, "cudaUnbindTexture failed on tex");
    destroy_cudata_array(d_data->nbfp);

    stat = cudaFree(d_data->shift_vec);
    CU_RET_ERR(stat, "cudaEventDestroy failed on d_data->shift_vec");

    stat = cudaFree(d_data->e_lj);
    CU_RET_ERR(stat, "cudaEventDestroy failed on d_data->e_lj");
    stat = cudaFree(d_data->e_el);
    CU_RET_ERR(stat, "cudaEventDestroy failed on d_data->e_el");

    destroy_cudata_array(d_data->f, &d_data->natoms, &d_data->nalloc);
    destroy_cudata_array(d_data->xq);
    destroy_cudata_array(d_data->atom_types, &d_data->ntypes);            

    destroy_cudata_array(d_data->ci, &d_data->nci, &d_data->ci_nalloc);
    destroy_cudata_array(d_data->sj4, &d_data->nsj4, &d_data->sj4_nalloc);
    destroy_cudata_array(d_data->excl, &d_data->nexcl, &d_data->excl_nalloc);

    stat = cudaThreadExit();
    CU_RET_ERR(stat, "cudaThreadExit failed");

    fprintf(fplog, "Cleaned up CUDA data structures.\n");
}

/* Frees the device memory pointed by d_ptr and resets the associated 
   size and allocation size variables to -1.
 */
static void destroy_cudata_array(void *d_ptr, int *n, int *nalloc)
{
    cudaError_t stat;
    
    if (d_ptr)
    {
        stat = cudaFree(d_ptr);
        CU_RET_ERR(stat, "cudaFree failed");
    }

    if (n)
    {        
        *n = -1;
    }

    if (nalloc)
    {
        *nalloc = -1;
    }
}

/* Reallocates the device memory pointed by d_ptr and copies the data from the 
   location pointed by h_src host-side pointer. Allocation is buffered and 
   therefor freeing is only needed if the previously allocated space is not 
   enough. 
 */
static void realloc_cudata_array(void **d_dest, void *h_src, size_t type_size, 
                                 int *curr_size, int *curr_alloc_size, 
                                 int req_size, gmx_bool doStream)
{
    cudaError_t stat;

    if (d_dest == NULL || req_size <= 0)
    {
        return;
    }

    /* reallocate only if the data does not fit = allocation size is smaller 
       than the current requested size */
    if (req_size > *curr_alloc_size)
    {
        /* only free if the array has already been initialized */
        if (*curr_alloc_size >= 0)
        {
            destroy_cudata_array(*d_dest, curr_size, curr_alloc_size);
        }

        *curr_alloc_size = 1.2 * req_size + 100;  /* TODO replace this with a fn pointer 
                                                     passed from outside */

        stat = cudaMalloc(d_dest, *curr_alloc_size * type_size);
        CU_RET_ERR(stat, "cudaMalloc failed in realloc_cudata_array");
    }

    /* size could have changed without actual reallocation */
    *curr_size = req_size;

    /* upload to device */
    if (h_src)
    {
        if(doStream)
        {
            upload_cudata_async(*d_dest, h_src, *curr_size * type_size, 0);
        }
        else 
        {
            upload_cudata(*d_dest, h_src,  *curr_size * type_size);
        }
    }
}

void cu_blockwait_atomdata(t_cudata d_data)
{   
    float t;
    cu_blockwait_event(d_data->stop_atdat, d_data->start_atdat, &t);
    d_data->timings.atomdt_h2d_total_time += t;
}

void cu_blockwait_atomdata_OLD(t_cudata d_data, float *time)
{    
    cudaError_t stat;     

    stat = cudaEventSynchronize(d_data->stop_atdat);
    CU_RET_ERR(stat, "the async trasfer of atomdata has failed");   

    stat = cudaEventElapsedTime(time, d_data->start_atdat, d_data->stop_atdat);
    CU_RET_ERR(stat, "cudaEventElapsedTime on start_atdat and stop_atdat failed");
}

/* GPU timerstruct  query functions */
gpu_times_t * get_gpu_times(t_cudata d_data)
{
    return &d_data->timings;
}

int cu_upload_X(t_cudata d_data, real *h_x) 
{
    return upload_cudata(d_data->xq, h_x, d_data->natoms*sizeof(*d_data->xq));
}

int cu_download_F(real *h_f, t_cudata d_data)
{
    return download_cudata(h_f, d_data->f, d_data->natoms*sizeof(*d_data->f));
}
