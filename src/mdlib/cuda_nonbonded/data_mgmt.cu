#include <stdlib.h>
#include <stdio.h>

#include "gmx_fatal.h"
#include "smalloc.h"

#include "cutypedefs.h"
#include "cudautils.h"
#include "cuda_data_mgmt.h"
#include "cupmalloc.h"

#define USE_CUDA_EVENT_BLOCKING_SYNC FALSE /* makes the CPU thread busy-wait! */
#define EWALD_COULOMB_FORCE_TABLE_SIZE 1536   /* size chosen such we do not run out of texture cache */

#define MY_PI               (3.1415926535897932384626433832795)
#define TWO_OVER_SQRT_PI    (2.0/sqrt(MY_PI))
    
#define TIME_GPU_TRANSFERS 1

#define NUM_NB_KERNELS 12

static const char * const nb_k1_names[NUM_NB_KERNELS] = 
{
    "_Z21k_calc_nb_RF_forces_111cu_atomdata12cu_nb_params9cu_nblist",
    "_Z24k_calc_nb_ewald_forces_111cu_atomdata12cu_nb_params9cu_nblist",
    "_Z25k_calc_nb_cutoff_forces_111cu_atomdata12cu_nb_params9cu_nblist",
    "_Z30k_calc_nb_RF_forces_energies_111cu_atomdata12cu_nb_params9cu_nblist",
    "_Z30k_calc_nb_RF_forces_prunenbl_111cu_atomdata12cu_nb_params9cu_nblist",
    "_Z33k_calc_nb_ewald_forces_energies_111cu_atomdata12cu_nb_params9cu_nblist",
    "_Z33k_calc_nb_ewald_forces_prunenbl_111cu_atomdata12cu_nb_params9cu_nblist",
    "_Z34k_calc_nb_cutoff_forces_energies_111cu_atomdata12cu_nb_params9cu_nblist",
    "_Z34k_calc_nb_cutoff_forces_prunenbl_111cu_atomdata12cu_nb_params9cu_nblist",
    "_Z39k_calc_nb_RF_forces_energies_prunenbl_111cu_atomdata12cu_nb_params9cu_nblist",
    "_Z42k_calc_nb_ewald_forces_energies_prunenbl_111cu_atomdata12cu_nb_params9cu_nblist",
    "_Z43k_calc_nb_cutoff_forces_energies_prunenbl_111cu_atomdata12cu_nb_params9cu_nblist"
};

static const char * const nb_k2_names[NUM_NB_KERNELS] = 
{
    "_Z21k_calc_nb_RF_forces_211cu_atomdata12cu_nb_params9cu_nblist",
    "_Z24k_calc_nb_ewald_forces_211cu_atomdata12cu_nb_params9cu_nblist",
    "_Z25k_calc_nb_cutoff_forces_211cu_atomdata12cu_nb_params9cu_nblist",
    "_Z30k_calc_nb_RF_forces_energies_211cu_atomdata12cu_nb_params9cu_nblist",
    "_Z30k_calc_nb_RF_forces_prunenbl_211cu_atomdata12cu_nb_params9cu_nblist",
    "_Z33k_calc_nb_ewald_forces_energies_211cu_atomdata12cu_nb_params9cu_nblist",
    "_Z33k_calc_nb_ewald_forces_prunenbl_211cu_atomdata12cu_nb_params9cu_nblist",
    "_Z34k_calc_nb_cutoff_forces_energies_211cu_atomdata12cu_nb_params9cu_nblist",
    "_Z34k_calc_nb_cutoff_forces_prunenbl_211cu_atomdata12cu_nb_params9cu_nblist",
    "_Z39k_calc_nb_RF_forces_energies_prunenbl_211cu_atomdata12cu_nb_params9cu_nblist",
    "_Z42k_calc_nb_ewald_forces_energies_prunenbl_211cu_atomdata12cu_nb_params9cu_nblist",
    "_Z43k_calc_nb_cutoff_forces_energies_prunenbl_211cu_atomdata12cu_nb_params9cu_nblist"
};

__device__ __global__ void k_empty_test(){}

/*** CUDA Data operations ***/
static void destroy_cudata_array(void * d_ptr, 
                                 int * n = NULL, int * nalloc = NULL);
static void realloc_cudata_array(void **d_dest, void *h_src, size_t type_size, 
                                 int *curr_size, int *curr_alloc_size, 
                                 int req_size, gmx_bool doStream);                                
static void tabulate_ewald_coulomb_force_r(cu_nb_params_t *nb_params);

/*! Tabulates the Ewald Coulomb force.
 *  Idea borrowed from OpenMM (https://simtk.org/home/openmm).
 */
static void tabulate_ewald_coulomb_force_r(cu_nb_params_t *nb_params)
{
    float       *ftmp;
    double      beta, r, x;
    int         i, tabsize;
    cudaError_t stat;    

    cudaChannelFormatDesc   cd;
    const textureReference  *tex_coulomb_tab;

    beta        = nb_params->ewald_beta;
    tabsize     = EWALD_COULOMB_FORCE_TABLE_SIZE;

    nb_params->coulomb_tab_size    = tabsize;
    nb_params->coulomb_tab_scale   = (tabsize - 1) / sqrt(nb_params->cutoff_sq);

    pmalloc((void**)&ftmp, tabsize*sizeof(*ftmp));

    for (i = 1; i < tabsize; i++)
    {
        r       = i / nb_params->coulomb_tab_scale;
        x       = r * beta;
        ftmp[i] = (float) ((erfc(x) / r + beta * TWO_OVER_SQRT_PI * exp(-x * x)) / (r * r));
    }

    ftmp[0] = ftmp[1];

    stat = cudaMalloc((void **)&nb_params->coulomb_tab, 
                      tabsize*sizeof(*nb_params->coulomb_tab));
    CU_RET_ERR(stat, "cudaMalloc failed on nb_params->coulomb_tab"); 
    upload_cudata(nb_params->coulomb_tab, ftmp, tabsize*sizeof(*nb_params->coulomb_tab));

    stat = cudaGetTextureReference(&tex_coulomb_tab, "tex_coulomb_tab");
    CU_RET_ERR(stat, "cudaGetTextureReference on tex_coulomb_tab failed");
    cd = cudaCreateChannelDesc<float>();
    stat = cudaBindTexture(NULL, tex_coulomb_tab, nb_params->coulomb_tab, 
                           &cd, tabsize*sizeof(*nb_params->coulomb_tab));
    CU_RET_ERR(stat, "cudaBindTexture on tex_coulomb_tab failed");

    pfree(ftmp);
}

/*! Initilizes the atomdata (XXX) data structure. */
void init_atomdata(cu_atomdata_t *ad, int ntypes)
{
    cudaError_t stat;

    ad->ntypes  = ntypes;
    stat = cudaMalloc((void**)&ad->shift_vec, SHIFTS*sizeof(*ad->shift_vec));
    CU_RET_ERR(stat, "cudaMalloc failed on ad->shift_vec"); 
    ad->shift_vec_copied = FALSE;

    stat = cudaMalloc((void**)&ad->e_lj, sizeof(*ad->e_lj));
    CU_RET_ERR(stat, "cudaMalloc failed on ad->e_lj");
    stat = cudaMalloc((void**)&ad->e_el, sizeof(*ad->e_el));
    CU_RET_ERR(stat, "cudaMalloc failed on ad->e_el");

    /* initilize to NULL poiters to data that is not allocated here and will
       need reallocation in init_cudata_atoms */
    ad->xq = NULL;
    ad->f  = NULL;

    /* size -1 indicates that the repective array hasn't been initialized yet */
    ad->natoms = -1;
    ad->nalloc = -1;
}

/*! Initilizes the nonbonded parameter data structure. */
void init_nb_params(cu_nb_params_t *nbp, const t_forcerec *fr)
{  
    cudaError_t             stat;
    cudaChannelFormatDesc   cd;
    const textureReference  *tex_nbfp;
    
    int ntypes = fr->nbat->ntype;
    
    nbp->ewald_beta  = fr->ewaldcoeff;
    nbp->eps_r       = fr->epsilon_r;
    nbp->two_k_rf    = 2.0 * fr->k_rf;
    nbp->c_rf        = fr->c_rf;
    nbp->cutoff_sq   = fr->rcut_nsbox * fr->rcut_nsbox;
    nbp->rlist_sq    = fr->rlist_nsbox * fr->rlist_nsbox;
    
    if (fr->eeltype == eelCUT)
    {
        nbp->eeltype = cu_eelCUT;
    }
    else if (EEL_RF(fr->eeltype))
    {                
        nbp->eeltype = cu_eelRF;
    }
    else if ((EEL_PME(fr->eeltype) || fr->eeltype==eelEWALD))
    {
        nbp->eeltype = cu_eelEWALD;
    }
   else 
    {
        gmx_fatal(FARGS, "The requested electrostatics type is not implemented in the CUDA GPU accelerated kernels!");
    }

    /* generate table for PME */
    if (nbp->eeltype == cu_eelEWALD)
    {
        tabulate_ewald_coulomb_force_r(nbp);
    }

    stat = cudaMalloc((void **)&nbp->nbfp, 2*ntypes*ntypes*sizeof(*(nbp->nbfp)));
    CU_RET_ERR(stat, "cudaMalloc failed on nbp->nbfp"); 
    upload_cudata(nbp->nbfp, fr->nbat->nbfp, 2*ntypes*ntypes*sizeof(*(nbp->nbfp)));

    stat = cudaGetTextureReference(&tex_nbfp, "tex_nbfp");
    CU_RET_ERR(stat, "cudaGetTextureReference on tex_nbfp failed");
    cd = cudaCreateChannelDesc<float>();
    stat = cudaBindTexture(NULL, tex_nbfp, nbp->nbfp, &cd, 2*ntypes*ntypes*sizeof(*(nbp->nbfp)));
    CU_RET_ERR(stat, "cudaBindTexture on tex_nbfp failed");
}

/*! Initilizes the neighborlist data structure. */
void init_nblist(cu_nblist_t *nbl)
{
    /* initilize to NULL poiters to data that is not allocated here and will
       need reallocation in init_cudata_atoms */
    nbl->ci      = NULL;
    nbl->sj4     = NULL;
    nbl->excl    = NULL;    
    
    /* size -1 indicates that the repective array hasn't been initialized yet */
    nbl->naps        = -1;
    nbl->nci         = -1;
    nbl->ci_nalloc   = -1;
    nbl->nsj4        = -1;
    nbl->sj4_nalloc  = -1;
    nbl->nexcl       = -1;
    nbl->excl_nalloc = -1;
    nbl->prune_nbl   = FALSE;
}

/*! Initilizes the timer data structure. */
static void init_timers(cu_timers_t *t)
{
    cudaError_t stat;
    /* XXX */ 
    int eventflags = ( USE_CUDA_EVENT_BLOCKING_SYNC ? cudaEventBlockingSync: cudaEventDefault );

    stat = cudaEventCreateWithFlags(&(t->start_nb), eventflags);
    CU_RET_ERR(stat, "cudaEventCreate on start_nb failed");
    stat = cudaEventCreateWithFlags(&(t->stop_nb), eventflags);
    CU_RET_ERR(stat, "cudaEventCreate on stop_nb failed");
    stat = cudaEventCreateWithFlags(&(t->start_atdat), eventflags);
    CU_RET_ERR(stat, "cudaEventCreate on start_atdat failed");
    stat = cudaEventCreateWithFlags(&(t->stop_atdat), eventflags);
    CU_RET_ERR(stat, "cudaEventCreate on stop_atdat failed");

    t->time_transfers = TIME_GPU_TRANSFERS > 0; /* XXX fix this! */

    if (t->time_transfers)
    {
        stat = cudaEventCreateWithFlags(&(t->start_nb_h2d), eventflags);
        CU_RET_ERR(stat, "cudaEventCreate on start_nb_h2d failed");
        stat = cudaEventCreateWithFlags(&(t->stop_nb_h2d), eventflags);
        CU_RET_ERR(stat, "cudaEventCreate on stop_nb_h2d failed");

        stat = cudaEventCreateWithFlags(&(t->start_nb_d2h), eventflags);
        CU_RET_ERR(stat, "cudaEventCreate on start_nb_d2h failed");
        stat = cudaEventCreateWithFlags(&(t->stop_nb_d2h), eventflags);
        CU_RET_ERR(stat, "cudaEventCreate on stop_nb_d2h failed");
    }
}

/*! Initilizes the timer data structure. */
static void init_timings(cu_timings_t *t)
{
    int i, j;

    t->nb_h2d_time = 0.0;
    t->nb_d2h_time = 0.0;
    t->nb_count    = 0;
    t->atomdt_h2d_total_time = 0.0;
    t->atomdt_count = 0;
    for (i = 0; i < 2; i++)
    {
        for(j = 0; j < 2; j++)
        {
            t->k_time[i][j].t = 0.0;
            t->k_time[i][j].c = 0;
        }
    }
}

/*! Initilizes force-field related data (called only once in the beginning).
 */
void init_cudata_ff(FILE *fplog, 
                    cu_nonbonded_t *p_cu_nb,
                    const t_forcerec *fr)
{
    cudaError_t     stat;
    cu_nonbonded_t  nb;

    if (p_cu_nb == NULL) return;
    
    snew(nb, 1); 
    snew(nb->atomdata, 1); 
    snew(nb->nb_params, 1); 
    snew(nb->nblist, 1); 
    snew(nb->timers, 1); 
    snew(nb->timings, 1); 

    init_atomdata(nb->atomdata, fr->nbat->ntype);
    init_nb_params(nb->nb_params, fr);
    init_nblist(nb->nblist);
    init_timers(nb->timers);
    init_timings(nb->timings);
    /* init tmpdata */
    pmalloc((void**)&nb->tmpdata.e_lj, sizeof(*nb->tmpdata.e_lj));
    pmalloc((void**)&nb->tmpdata.e_el, sizeof(*nb->tmpdata.e_el));

    nb->streamGPU   = fr->streamGPU;
    *p_cu_nb = nb;

    if (fplog != NULL)
    {
        fprintf(fplog, "Initialized CUDA data structures.\n");
    }

    /* k_calc_nb_*_1 48/16 kB Shared/L1 */
    for (int i = 0; i < NUM_NB_KERNELS; i++)
    {
        stat = cudaFuncSetCacheConfig(nb_k1_names[i],  cudaFuncCachePreferShared);
        // printf("--> %s\n", nb_k1_names[i]);
        CU_RET_ERR(stat, "cudaFuncSetCacheConfig failed");
    }

    /* k_calc_nb_*_2 16/48 kB Shared/L1 */
    for (int i = 0; i < NUM_NB_KERNELS; i++)
    {
        stat = cudaFuncSetCacheConfig(nb_k2_names[i], cudaFuncCachePreferL1);
        CU_RET_ERR(stat, "cudaFuncSetCacheConfig failed");
    }

    /* TODO: move this to gpu_utils module */
    k_empty_test<<<1, 512>>>();
    CU_LAUNCH_ERR_SYNC("test kernel");
}

/*! Initilizes atom-data and neighbor list for the GPU, called at every 
 *  neighbor search step. 
 */
void init_cudata_atoms(cu_nonbonded_t cu_nb, 
                       const gmx_nb_atomdata_t *nbat, 
                       const gmx_nblist_t *h_nblist,
                       gmx_bool doStream)
{
    cudaError_t stat;
    char        sbuf[200];
    int         nalloc;
    int         natoms  = nbat->natoms;
    cu_nblist_t *d_nblist   = cu_nb->nblist;
    cu_atomdata_t *d_atomd  = cu_nb->atomdata;
    cu_timers_t *timers     = cu_nb->timers;  
    
    /* time async copy */
    stat = cudaEventRecord(timers->start_atdat, 0);
    CU_RET_ERR(stat, "cudaEventRecord failed on timers->start_atdat");

    if (d_nblist->naps < 0)
    {
        d_nblist->naps = h_nblist->naps;
    }
    else
    {
        if (d_nblist->naps != h_nblist->naps)
        {
            sprintf(sbuf, "In init_cudata_atoms: the #atoms per cell has changed (from %d to %d)",
                    d_nblist->naps, h_nblist->naps);            
            gmx_incons(sbuf);
        }
    }

    /* need to reallocate if we have to copy more atoms than the amount of space
       available and only allocate if we haven"t initilzed yet, i.e d_nblist->natoms == -1 */
    if (natoms > d_atomd->nalloc)
    {
        nalloc = natoms * 1.2 + 100;
    
        /* free up first if the arrays have already been initialized */
        if (d_atomd->nalloc != -1)
        {
            destroy_cudata_array(d_atomd->f, &d_atomd->natoms, &d_atomd->nalloc);
            destroy_cudata_array(d_atomd->xq);
            destroy_cudata_array(d_atomd->atom_types, &d_atomd->ntypes);             
        }
        
        stat = cudaMalloc((void **)&d_atomd->f, nalloc*sizeof(*(d_atomd->f)));
        CU_RET_ERR(stat, "cudaMalloc failed on d_atomd->f");                   
        stat = cudaMalloc((void **)&d_atomd->xq, nalloc*sizeof(*(d_atomd->xq)));
        CU_RET_ERR(stat, "cudaMalloc failed on d_atomd->xq");     

        stat = cudaMalloc((void **)&d_atomd->atom_types, nalloc*sizeof(*(d_atomd->atom_types)));
        CU_RET_ERR(stat, "cudaMalloc failed on d_atomd->atom_types"); 

        d_atomd->nalloc = nalloc;
    }
    /* XXX for the moment we just set all 8 values to the same value... 
       ATM not, we"ll do that later */    
    d_atomd->natoms = natoms;

    if(doStream)
    {
        upload_cudata_async(d_atomd->atom_types, nbat->type, 
                            natoms*sizeof(*d_atomd->atom_types), 0);
    }
    else 
    {
        upload_cudata(d_atomd->atom_types, nbat->type, 
                      natoms*sizeof(*(d_atomd->atom_types)));
    
    }

    realloc_cudata_array((void **)&d_nblist->ci, h_nblist->ci, sizeof(*(d_nblist->ci)),
                         &d_nblist->nci, &d_nblist->ci_nalloc,
                         h_nblist->nci, doStream);

    realloc_cudata_array((void **)&d_nblist->sj4, h_nblist->sj4, sizeof(*(d_nblist->sj4)),
                         &d_nblist->nsj4, &d_nblist->sj4_nalloc,
                         h_nblist->nsj4, doStream);

    realloc_cudata_array((void **)&d_nblist->excl, h_nblist->excl, sizeof(*(d_nblist->excl)),
                         &d_nblist->nexcl, &d_nblist->excl_nalloc,
                         h_nblist->nexcl, doStream);

    stat = cudaEventRecord(timers->stop_atdat, 0);
    CU_RET_ERR(stat, "cudaEventRecord failed on timers->stop_atdat");

    cu_nb->timings->atomdt_count++;

    /* need to prune the neighbor list during the next step */
    d_nblist->prune_nbl = TRUE;
}

/*! Frees up all GPU resources used for the nonbonded calculations. */
void destroy_cudata(FILE *fplog, cu_nonbonded_t cu_nb)
{
    cudaError_t stat;
    const textureReference  *tex;
    cu_atomdata_t       *atomdata;
    cu_nb_params_t      *nb_params;
    cu_nblist_t         *nblist;
    cu_timers_t         *timers;

    atomdata    = cu_nb->atomdata;
    nb_params   = cu_nb->nb_params;
    nblist      = cu_nb->nblist;
    timers      = cu_nb->timers;

    if (cu_nb == NULL) return;

    if (nb_params->eeltype == cu_eelEWALD)
    {
        stat = cudaGetTextureReference(&tex, "tex_coulomb_tab");
        CU_RET_ERR(stat, "cudaGetTextureReference on tex_coulomb_tab failed");
        stat = cudaUnbindTexture(tex);
        CU_RET_ERR(stat, "cudaUnbindTexture failed on tex");
        destroy_cudata_array(nb_params->coulomb_tab, &nb_params->coulomb_tab_size);            
    }

    stat = cudaEventDestroy(timers->start_nb);
    CU_RET_ERR(stat, "cudaEventDestroy failed on timers->start_nb");
    stat = cudaEventDestroy(timers->stop_nb);
    CU_RET_ERR(stat, "cudaEventDestroy failed on timers->stop_nb");
    stat = cudaEventDestroy(timers->start_atdat);
    CU_RET_ERR(stat, "cudaEventDestroy failed on timers->start_atdat");
    stat = cudaEventDestroy(timers->stop_atdat);
    CU_RET_ERR(stat, "cudaEventDestroy failed on timers->stop_atdat");

    if (timers->time_transfers)
    {
        stat = cudaEventDestroy(timers->start_nb_h2d);
        CU_RET_ERR(stat, "cudaEventDestroy failed on timers->start_nb_h2d");
        stat = cudaEventDestroy(timers->stop_nb_h2d);
        CU_RET_ERR(stat, "cudaEventDestroy failed on timers->stop_nb_h2d");

        stat = cudaEventDestroy(timers->start_nb_d2h);
        CU_RET_ERR(stat, "cudaEventDestroy failed on timers->start_nb_d2h");
        stat = cudaEventDestroy(timers->stop_nb_d2h);
        CU_RET_ERR(stat, "cudaEventDestroy failed on timers->stop_nb_d2h");
    }

    stat = cudaGetTextureReference(&tex, "tex_nbfp");
    CU_RET_ERR(stat, "cudaGetTextureReference on tex_nbfp failed");
    stat = cudaUnbindTexture(tex);
    CU_RET_ERR(stat, "cudaUnbindTexture failed on tex");
    destroy_cudata_array(nb_params->nbfp);

    stat = cudaFree(atomdata->shift_vec);
    CU_RET_ERR(stat, "cudaEventDestroy failed on atomdata->shift_vec");

    stat = cudaFree(atomdata->e_lj);
    CU_RET_ERR(stat, "cudaEventDestroy failed on atomdata->e_lj");
    stat = cudaFree(atomdata->e_el);
    CU_RET_ERR(stat, "cudaEventDestroy failed on atomdata->e_el");

    destroy_cudata_array(atomdata->f, &atomdata->natoms, &atomdata->nalloc);
    destroy_cudata_array(atomdata->xq);
    destroy_cudata_array(atomdata->atom_types, &atomdata->ntypes);            

    destroy_cudata_array(nblist->ci, &nblist->nci, &nblist->ci_nalloc);
    destroy_cudata_array(nblist->sj4, &nblist->nsj4, &nblist->sj4_nalloc);
    destroy_cudata_array(nblist->excl, &nblist->nexcl, &nblist->excl_nalloc);

    stat = cudaThreadExit();
    CU_RET_ERR(stat, "cudaThreadExit failed");

    fprintf(fplog, "Cleaned up CUDA data structures.\n");
}

/*! Frees the device memory pointed by d_ptr and resets the associated 
 *  size and allocation size variables to -1.
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

/*! Reallocates the device memory pointed by d_ptr and copies the data from the 
 * location pointed by h_src host-side pointer. Allocation is buffered and 
 * therefor freeing is only needed if the previously allocated space is not 
 * enough. 
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

/*! Blocking waits until the atom data gets copied to the GPU and times the transfer.
 */
void cu_blockwait_atomdata(cu_nonbonded_t nb_cu)
{   
    float t;
    cu_blockwait_event(nb_cu->timers->stop_atdat, nb_cu->timers->start_atdat, &t);
    nb_cu->timings->atomdt_h2d_total_time += t;
}

/*! Returns the GPU timing structure or NULL if nb_cu is NULL. */
cu_timings_t * get_gpu_times(cu_nonbonded_t nb_cu)
{
    return nb_cu != NULL ? nb_cu->timings : NULL;
}

/*** Old stuff ***/
int cu_upload_X(cu_nonbonded_t nb_cu, real *h_x) 
{
    cu_atomdata_t *ad = nb_cu->atomdata;

    return upload_cudata(ad->xq, h_x, ad->natoms*sizeof(*ad->xq));
}

int cu_download_F(real *h_f, cu_nonbonded_t nb_cu)
{
    cu_atomdata_t *ad = nb_cu->atomdata;

    return download_cudata(h_f, ad->f, ad->natoms*sizeof(*ad->f));
}
