#include <stdlib.h>
#include <stdio.h>

#include "gmx_fatal.h"
#include "smalloc.h"
#include "force.h"

#include "cutypedefs.h"
#include "cudautils.h"
#include "cuda_data_mgmt.h"
#include "cupmalloc.h"

#define USE_CUDA_EVENT_BLOCKING_SYNC FALSE /* makes the CPU thread busy-wait! */
#define EWALD_COULOMB_FORCE_TABLE_SIZE (1536)   /* size chosen such we do not run out of texture cache */

#define MY_PI               (3.1415926535897932384626433832795)
#define TWO_OVER_SQRT_PI    (2.0/sqrt(MY_PI))
    
#define TIME_GPU_TRANSFERS 1

#define NUM_NB_KERNELS 12

static const char * const nb_k1_names[NUM_NB_KERNELS] = 
{
    "_Z21k_calc_nb_RF_forces_111cu_atomdata12cu_nb_params9cu_nblisti",
    "_Z24k_calc_nb_ewald_forces_111cu_atomdata12cu_nb_params9cu_nblisti",
    "_Z25k_calc_nb_cutoff_forces_111cu_atomdata12cu_nb_params9cu_nblisti",
    "_Z30k_calc_nb_RF_forces_energies_111cu_atomdata12cu_nb_params9cu_nblisti",
    "_Z30k_calc_nb_RF_forces_prunenbl_111cu_atomdata12cu_nb_params9cu_nblisti",
    "_Z33k_calc_nb_ewald_forces_energies_111cu_atomdata12cu_nb_params9cu_nblisti",
    "_Z33k_calc_nb_ewald_forces_prunenbl_111cu_atomdata12cu_nb_params9cu_nblisti",
    "_Z34k_calc_nb_cutoff_forces_energies_111cu_atomdata12cu_nb_params9cu_nblisti",
    "_Z34k_calc_nb_cutoff_forces_prunenbl_111cu_atomdata12cu_nb_params9cu_nblisti",
    "_Z39k_calc_nb_RF_forces_energies_prunenbl_111cu_atomdata12cu_nb_params9cu_nblisti",
    "_Z42k_calc_nb_ewald_forces_energies_prunenbl_111cu_atomdata12cu_nb_params9cu_nblisti",
    "_Z43k_calc_nb_cutoff_forces_energies_prunenbl_111cu_atomdata12cu_nb_params9cu_nblisti"
};

static const char * const nb_k2_names[NUM_NB_KERNELS] = 
{
    "_Z21k_calc_nb_RF_forces_211cu_atomdata12cu_nb_params9cu_nblisti",
    "_Z24k_calc_nb_ewald_forces_211cu_atomdata12cu_nb_params9cu_nblisti",
    "_Z25k_calc_nb_cutoff_forces_211cu_atomdata12cu_nb_params9cu_nblisti",
    "_Z30k_calc_nb_RF_forces_energies_211cu_atomdata12cu_nb_params9cu_nblisti",
    "_Z30k_calc_nb_RF_forces_prunenbl_211cu_atomdata12cu_nb_params9cu_nblisti",
    "_Z33k_calc_nb_ewald_forces_energies_211cu_atomdata12cu_nb_params9cu_nblisti",
    "_Z33k_calc_nb_ewald_forces_prunenbl_211cu_atomdata12cu_nb_params9cu_nblisti",
    "_Z34k_calc_nb_cutoff_forces_energies_211cu_atomdata12cu_nb_params9cu_nblisti",
    "_Z34k_calc_nb_cutoff_forces_prunenbl_211cu_atomdata12cu_nb_params9cu_nblisti",
    "_Z39k_calc_nb_RF_forces_energies_prunenbl_211cu_atomdata12cu_nb_params9cu_nblisti",
    "_Z42k_calc_nb_ewald_forces_energies_prunenbl_211cu_atomdata12cu_nb_params9cu_nblisti",
    "_Z43k_calc_nb_cutoff_forces_energies_prunenbl_211cu_atomdata12cu_nb_params9cu_nblisti"
};

__device__ __global__ void k_empty_test(){}

/*** CUDA Data operations ***/
static void destroy_cudata_array(void * d_ptr, 
                                 int * n = NULL, int * nalloc = NULL);
static void realloc_cudata_array(void **d_dest, void *h_src, size_t type_size, 
                                 int *curr_size, int *curr_alloc_size, 
                                 int req_size,
                                 cudaStream_t stream, gmx_bool doStream);
static void init_ewald_coulomb_force_table(cu_nb_params_t *nb_params);

/*! Tabulates the Ewald Coulomb force and initializes the related GPU resources. 
 */
static void init_ewald_coulomb_force_table(cu_nb_params_t *nb_params)
{
    float       *ftmp, *coul_tab;
    int         tabsize;
    double      tabscale;
    cudaError_t stat;

    tabsize     = EWALD_COULOMB_FORCE_TABLE_SIZE;
    tabscale    = (tabsize - 1) / sqrt(nb_params->cutoff_sq);

    pmalloc((void**)&ftmp, tabsize*sizeof(*ftmp));

    table_spline3_fill_ewald_force(ftmp, tabsize, 1/tabscale, nb_params->ewald_beta);

    stat = cudaMalloc((void **)&coul_tab, tabsize*sizeof(*coul_tab));
    CU_RET_ERR(stat, "cudaMalloc failed on coul_tab");
    upload_cudata(coul_tab, ftmp, tabsize*sizeof(*coul_tab));
    cu_bind_texture("tex_coulomb_tab", coul_tab, tabsize*sizeof(*coul_tab));

    nb_params->coulomb_tab          = coul_tab;
    nb_params->coulomb_tab_size     = tabsize;
    nb_params->coulomb_tab_scale    = tabscale;

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

    stat = cudaMalloc((void**)&ad->f_shift, SHIFTS*sizeof(*ad->f_shift));
    CU_RET_ERR(stat, "cudaMalloc failed on ad->f_shift");

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
    cudaError_t stat;
    int         ntypes, nnbfp; 

    ntypes = fr->nbat->ntype;
    
    nbp->ewald_beta  = fr->ewaldcoeff;
    nbp->eps_r       = fr->epsilon_r;
    nbp->two_k_rf    = 2.0 * fr->k_rf;
    nbp->c_rf        = fr->c_rf;
    nbp->cutoff_sq   = fr->rvdw * fr->rvdw;
    nbp->rlist_sq    = fr->rlist * fr->rlist;
    nbp->lj_shift    = (getenv("GMX_LJ_SHIFT") == NULL) ?
             0.0 : -1/(nbp->cutoff_sq * nbp->cutoff_sq * nbp->cutoff_sq);

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
        init_ewald_coulomb_force_table(nbp);
    }

    nnbfp = 2*ntypes*ntypes;
    stat = cudaMalloc((void **)&nbp->nbfp, nnbfp*sizeof(*nbp->nbfp));
    CU_RET_ERR(stat, "cudaMalloc failed on nbp->nbfp"); 
    upload_cudata(nbp->nbfp, fr->nbat->nbfp, nnbfp*sizeof(*nbp->nbfp));
    cu_bind_texture("tex_nbfp", nbp->nbfp, nnbfp*sizeof(*nbp->nbfp));
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
    stat = cudaEventCreateWithFlags(&(t->start_nb_nl), eventflags);
    CU_RET_ERR(stat, "cudaEventCreate on start_nb_nl failed");
    stat = cudaEventCreateWithFlags(&(t->stop_nb_nl), eventflags);
    CU_RET_ERR(stat, "cudaEventCreate on stop_nb_nl failed");

    stat = cudaEventCreateWithFlags(&(t->start_clear), eventflags);
    CU_RET_ERR(stat, "cudaEventCreate on start_clear failed");
    stat = cudaEventCreateWithFlags(&(t->stop_clear), eventflags);
    CU_RET_ERR(stat, "cudaEventCreate on stop_clear failed");

    stat = cudaEventCreateWithFlags(&(t->start_atdat), eventflags);
    CU_RET_ERR(stat, "cudaEventCreate on start_atdat failed");
    stat = cudaEventCreateWithFlags(&(t->stop_atdat), eventflags);
    CU_RET_ERR(stat, "cudaEventCreate on stop_atdat failed");
    stat = cudaEventCreateWithFlags(&(t->start_atdat_nl), eventflags);
    CU_RET_ERR(stat, "cudaEventCreate on start_atdat_nl failed");
    stat = cudaEventCreateWithFlags(&(t->stop_atdat_nl), eventflags);
    CU_RET_ERR(stat, "cudaEventCreate on stop_atdat_nl failed");

    stat = cudaStreamCreate(&t->nbstream);
    CU_RET_ERR(stat, "cudaStreamCreate on nbstream failed");
    stat = cudaStreamCreate(&t->nbstream_nl);
    CU_RET_ERR(stat, "cudaStreamCreate on nbstream_nl failed");

    t->time_transfers = TIME_GPU_TRANSFERS > 0; /* XXX fix this! */

    if (t->time_transfers)
    {
        stat = cudaEventCreateWithFlags(&(t->start_nb_h2d), eventflags);
        CU_RET_ERR(stat, "cudaEventCreate on start_nb_h2d failed");
        stat = cudaEventCreateWithFlags(&(t->stop_nb_h2d), eventflags);
        CU_RET_ERR(stat, "cudaEventCreate on stop_nb_h2d failed");
        stat = cudaEventCreateWithFlags(&(t->start_nb_h2d_nl), eventflags);
        CU_RET_ERR(stat, "cudaEventCreate on start_nb_h2d_nl failed");
        stat = cudaEventCreateWithFlags(&(t->stop_nb_h2d_nl), eventflags);
        CU_RET_ERR(stat, "cudaEventCreate on stop_nb_h2d_nl failed");

        stat = cudaEventCreateWithFlags(&(t->start_nb_d2h), eventflags);
        CU_RET_ERR(stat, "cudaEventCreate on start_nb_d2h failed");
        stat = cudaEventCreateWithFlags(&(t->stop_nb_d2h), eventflags);
        CU_RET_ERR(stat, "cudaEventCreate on stop_nb_d2h failed");
        stat = cudaEventCreateWithFlags(&(t->start_nb_d2h_nl), eventflags);
        CU_RET_ERR(stat, "cudaEventCreate on start_nb_d2h_nl failed");
        stat = cudaEventCreateWithFlags(&(t->stop_nb_d2h_nl), eventflags);
        CU_RET_ERR(stat, "cudaEventCreate on stop_nb_d2h_nl failed");
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
    snew(nb->nblist_nl, 1);
    snew(nb->timers, 1); 
    snew(nb->timings, 1); 

    init_atomdata(nb->atomdata, fr->nbat->ntype);
    init_nb_params(nb->nb_params, fr);
    init_nblist(nb->nblist);
    init_nblist(nb->nblist_nl);
    init_timers(nb->timers);
    init_timings(nb->timings);

    /* init tmpdata */
    pmalloc((void**)&nb->tmpdata.e_lj, sizeof(*nb->tmpdata.e_lj));
    pmalloc((void**)&nb->tmpdata.e_el, sizeof(*nb->tmpdata.e_el));
    pmalloc((void**)&nb->tmpdata.f_shift, SHIFTS * sizeof(*nb->tmpdata.f_shift));

    nb->streamGPU   = fr->streamGPU;
    *p_cu_nb = nb;

    if (debug)
    {
        fprintf(debug, "Initialized CUDA data structures.\n");
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

/*! Initilizes neighbor list on the GPU, called at every neighbor search step. 
 */
void init_cudata_nblist(cu_nonbonded_t cu_nb, 
                        const gmx_nblist_t *h_nblist,
                        gmx_bool nonLocal,
                        gmx_bool doStream)
{
    char        sbuf[STRLEN];

    cudaStream_t stream     = nonLocal ? cu_nb->timers->nbstream_nl : cu_nb->timers->nbstream;
    cu_nblist_t *d_nblist   = nonLocal ? cu_nb->nblist_nl : cu_nb->nblist;
    //cu_timers_t *timers   = cu_nb->timers;  // FIXME

    if (d_nblist->naps < 0)
    {
        d_nblist->naps = h_nblist->naps;
    }
    else
    {
        if (d_nblist->naps != h_nblist->naps)
        {
            sprintf(sbuf, "In init_cudata_nblist: the #atoms per cell has changed (from %d to %d)",
                    d_nblist->naps, h_nblist->naps);            
            gmx_incons(sbuf);
        }
    }

    realloc_cudata_array((void **)&d_nblist->ci, h_nblist->ci, sizeof(*(d_nblist->ci)),
                         &d_nblist->nci, &d_nblist->ci_nalloc,
                         h_nblist->nci,
                         stream, doStream);

    realloc_cudata_array((void **)&d_nblist->sj4, h_nblist->sj4, sizeof(*(d_nblist->sj4)),
                         &d_nblist->nsj4, &d_nblist->sj4_nalloc,
                         h_nblist->nsj4,
                         stream, doStream);

    realloc_cudata_array((void **)&d_nblist->excl, h_nblist->excl, sizeof(*(d_nblist->excl)),
                         &d_nblist->nexcl, &d_nblist->excl_nalloc,
                         h_nblist->nexcl, 
                         stream, doStream);

    /* need to prune the neighbor list during the next step */
    d_nblist->prune_nbl = TRUE;
}

void cu_move_shift_vec(cu_nonbonded_t cu_nb, 
                       const gmx_nb_atomdata_t *nbatom)
{
    cu_atomdata_t   *adat = cu_nb->atomdata;

    /* HtoD shift vec if we have a dynamic box */
    if (nbatom->dynamic_box || !adat->shift_vec_copied)
    {
        upload_cudata_async(adat->shift_vec, nbatom->shift_vec, SHIFTS * sizeof(*adat->shift_vec), 0);
        adat->shift_vec_copied = TRUE;
    }
}

/* FIXME put all the clear ops into a stream, otherwise it won't overlap with anything  */
void cu_clear_nb_outputs(cu_nonbonded_t cu_nb, 
                         const gmx_nb_atomdata_t *nbatom, // FIXME VEEERY dirty
                         int flags)
{
    cudaError_t stat;

    cu_atomdata_t   *adat = cu_nb->atomdata;
    cu_timers_t     *timers = cu_nb->timers;

    gmx_bool calc_ene   = flags & GMX_FORCE_VIRIAL;
    gmx_bool calc_fshift = flags & GMX_FORCE_VIRIAL;

    /* FIXME: this is not a clear OP! */
    cu_move_shift_vec(cu_nb, nbatom);

    stat = cudaEventRecord(timers->start_clear, 0);
    CU_RET_ERR(stat, "cudaEventRecord on start_clear falied");

    cudaMemsetAsync(adat->f, 0, adat->natoms * sizeof(*adat->f), 0);

    /* set the shift force output to 0 */
    if (calc_fshift)
    {
        cudaMemsetAsync(adat->f_shift, 0, SHIFTS * sizeof(*adat->f_shift), 0);
    }

    /* set energy outputs to 0 */
    if (calc_ene)
    {
        cudaMemsetAsync(adat->e_lj, 0, sizeof(*adat->e_lj), 0);
        cudaMemsetAsync(adat->e_el, 0, sizeof(*adat->e_el), 0);
    }

    cudaEventRecord(timers->stop_clear, 0);
    CU_RET_ERR(stat, "cudaEventRecord on stop_clear falied");

    /* block all future streams until this finishes */
    // XXX this is too restrictive stat = cudaStreamWaitEvent(NULL, timers->stop_clear 0);
}

/*! Initilizes atom-data on the GPU, called at every neighbor search step. 
 */
void init_cudata_atoms(cu_nonbonded_t cu_nb,
                       const gmx_nb_atomdata_t *nbat,
                       gmx_bool doStream)
{
    cudaError_t stat;
    int         nalloc;
    int         natoms  = nbat->natoms;

    cu_atomdata_t *d_atomd  = cu_nb->atomdata;
    cu_timers_t *timers     = cu_nb->timers;  // FIXME

    /* time async copy */
    stat = cudaEventRecord(timers->start_atdat, 0);
    CU_RET_ERR(stat, "cudaEventRecord failed on timers->start_atdat");

    /* need to reallocate if we have to copy more atoms than the amount of space
       available and only allocate if we haven't initilzed yet, i.e d_atomd->natoms == -1 */
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
    
    d_atomd->natoms = natoms;
    d_atomd->natoms_local = nbat->natoms_local;

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

    stat = cudaEventRecord(timers->stop_atdat, 0);
    CU_RET_ERR(stat, "cudaEventRecord failed on timers->stop_atdat");

    cu_nb->timings->atomdt_count++;
}

/*! Frees up all GPU resources used for the nonbonded calculations. */
void destroy_cudata(FILE *fplog, cu_nonbonded_t cu_nb)
{
    cudaError_t stat;
    cu_atomdata_t       *atomdata;
    cu_nb_params_t      *nb_params;
    cu_nblist_t         *nblist, *nblist_nl;
    cu_timers_t         *timers;

    atomdata    = cu_nb->atomdata;
    nb_params   = cu_nb->nb_params;
    nblist      = cu_nb->nblist;
    nblist_nl   = cu_nb->nblist_nl;
    timers      = cu_nb->timers;

    if (cu_nb == NULL) return;

    if (nb_params->eeltype == cu_eelEWALD)
    {
        cu_unbind_texture("tex_coulomb_tab");
        destroy_cudata_array(nb_params->coulomb_tab, &nb_params->coulomb_tab_size);            
    }

    stat = cudaEventDestroy(timers->start_nb);
    CU_RET_ERR(stat, "cudaEventDestroy failed on timers->start_nb");
    stat = cudaEventDestroy(timers->stop_nb);
    CU_RET_ERR(stat, "cudaEventDestroy failed on timers->stop_nb");
    stat = cudaEventDestroy(timers->start_nb_nl);
    CU_RET_ERR(stat, "cudaEventDestroy failed on timers->start_nb_nl");
    stat = cudaEventDestroy(timers->stop_nb_nl);
    CU_RET_ERR(stat, "cudaEventDestroy failed on timers->stop_nb_nl");

    stat = cudaEventDestroy(timers->start_clear);
    CU_RET_ERR(stat, "cudaEventDestroy failed on timers->start_clear");
    stat = cudaEventDestroy(timers->stop_clear);
    CU_RET_ERR(stat, "cudaEventDestroy failed on timers->stop_clear");

    stat = cudaEventDestroy(timers->start_atdat);
    CU_RET_ERR(stat, "cudaEventDestroy failed on timers->start_atdat");
    stat = cudaEventDestroy(timers->stop_atdat);
    CU_RET_ERR(stat, "cudaEventDestroy failed on timers->stop_atdat");
    stat = cudaEventDestroy(timers->start_atdat_nl);
    CU_RET_ERR(stat, "cudaEventDestroy failed on timers->start_atdat_nl");
    stat = cudaEventDestroy(timers->stop_atdat_nl);
    CU_RET_ERR(stat, "cudaEventDestroy failed on timers->stop_atdat_nl");

    stat = cudaStreamDestroy(timers->nbstream);
    CU_RET_ERR(stat, "cudaStreamDestroy failed on nbstream");
    stat = cudaStreamDestroy(timers->nbstream_nl);
    CU_RET_ERR(stat, "cudaStreamDestroy failed on nbstream_nl");

    if (timers->time_transfers)
    {
        stat = cudaEventDestroy(timers->start_nb_h2d);
        CU_RET_ERR(stat, "cudaEventDestroy failed on timers->start_nb_h2d");
        stat = cudaEventDestroy(timers->stop_nb_h2d);
        CU_RET_ERR(stat, "cudaEventDestroy failed on timers->stop_nb_h2d");
        stat = cudaEventDestroy(timers->start_nb_h2d_nl);
        CU_RET_ERR(stat, "cudaEventDestroy failed on timers->start_nb_h2d_nl");
        stat = cudaEventDestroy(timers->stop_nb_h2d_nl);
        CU_RET_ERR(stat, "cudaEventDestroy failed on timers->stop_nb_h2d_nl");

        stat = cudaEventDestroy(timers->start_nb_d2h);
        CU_RET_ERR(stat, "cudaEventDestroy failed on timers->start_nb_d2h");
        stat = cudaEventDestroy(timers->stop_nb_d2h);
        CU_RET_ERR(stat, "cudaEventDestroy failed on timers->stop_nb_d2h");
        stat = cudaEventDestroy(timers->start_nb_d2h_nl);
        CU_RET_ERR(stat, "cudaEventDestroy failed on timers->start_nb_d2h_nl");
        stat = cudaEventDestroy(timers->stop_nb_d2h_nl);
        CU_RET_ERR(stat, "cudaEventDestroy failed on timers->stop_nb_d2h_nl");
    }

    cu_unbind_texture("tex_nbfp");
    destroy_cudata_array(nb_params->nbfp);

    stat = cudaFree(atomdata->shift_vec);
    CU_RET_ERR(stat, "cudaEventDestroy failed on atomdata->shift_vec");
    stat = cudaFree(atomdata->f_shift);
    CU_RET_ERR(stat, "cudaEventDestroy failed on atomdata->f_shift");

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
    destroy_cudata_array(nblist_nl->ci, &nblist_nl->nci, &nblist_nl->ci_nalloc);
    destroy_cudata_array(nblist_nl->sj4, &nblist_nl->nsj4, &nblist_nl->sj4_nalloc);
    destroy_cudata_array(nblist_nl->excl, &nblist_nl->nexcl, &nblist->excl_nalloc);

    stat = cudaThreadExit();
    CU_RET_ERR(stat, "cudaThreadExit failed");

    if (debug)
    {
        fprintf(debug, "Cleaned up CUDA data structures.\n");
    }
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
                                 int req_size, 
                                 cudaStream_t stream, gmx_bool doStream)
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
            upload_cudata_async(*d_dest, h_src, *curr_size * type_size, stream);
        }
        else 
        {
            upload_cudata(*d_dest, h_src,  *curr_size * type_size);
        }
    }
}

void cu_move_xq(cu_nonbonded_t cu_nb, const gmx_nb_atomdata_t *nbat,
                gmx_bool nonLocal)
{
    cu_atomdata_t   *d_nbat = cu_nb->atomdata;
    cudaStream_t    stream = nonLocal ? cu_nb->timers->nbstream_nl :
                                        cu_nb->timers->nbstream;

    upload_cudata_async(d_nbat->xq, nbat->x,
                        d_nbat->natoms * sizeof(*d_nbat->xq), stream);
}

/*! Blocking waits until the atom data gets copied to the GPU and times the transfer.
 */
void cu_blockwait_atomdata(cu_nonbonded_t cu_nb)
{
    float t;
    cu_blockwait_event(cu_nb->timers->stop_atdat, cu_nb->timers->start_atdat, &t);
    cu_nb->timings->atomdt_h2d_total_time += t;
}

/*! Calculated the ellapsed time during atomdata transfer.
 */
void cu_time_atomdata(cu_nonbonded_t cu_nb)
{
    float t;
    cudaError_t stat;

    stat = cudaEventElapsedTime(&t, cu_nb->timers->start_atdat, cu_nb->timers->stop_atdat);
    CU_RET_ERR(stat, "cudaEventElapsedTime failed in cu_blockwait_nb");
    cu_nb->timings->atomdt_h2d_total_time += t;
}

/*! Synchronizes the respective stream with the atomdata init operation.
 */
void cu_synchstream_atomdata(cu_nonbonded_t cu_nb, gmx_bool nonLocal)
{
    cudaError_t stat;
    cudaStream_t stream = nonLocal ? cu_nb->timers->nbstream_nl : cu_nb->timers->nbstream;

    stat = cudaStreamWaitEvent(stream, cu_nb->timers->stop_atdat, 0);
    CU_RET_ERR(stat, "cudaStreamWaitEvent failed");
}

/*! Returns the GPU timing structure or NULL if cu_nb is NULL. */
cu_timings_t * get_gpu_timings(cu_nonbonded_t cu_nb)
{
    return cu_nb != NULL ? cu_nb->timings : NULL;
}

/*! Resets GPU timers. */
void reset_gpu_timings(cu_nonbonded_t cu_nb)
{
    init_timings(cu_nb->timings);
}

/*** Old stuff ***/
int cu_upload_X(cu_nonbonded_t cu_nb, real *h_x) 
{
    cu_atomdata_t *ad = cu_nb->atomdata;

    return upload_cudata(ad->xq, h_x, ad->natoms*sizeof(*ad->xq));
}

int cu_download_F(real *h_f, cu_nonbonded_t cu_nb)
{
    cu_atomdata_t *ad = cu_nb->atomdata;

    return download_cudata(h_f, ad->f, ad->natoms*sizeof(*ad->f));
}
