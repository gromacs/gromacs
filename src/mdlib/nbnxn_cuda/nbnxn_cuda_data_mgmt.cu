#include <stdlib.h>
#include <stdio.h>

#include "gmx_fatal.h"
#include "smalloc.h"
#include "force.h"
#include "types/nb_verlet.h"
#include "types/interaction_const.h"

#include "nbnxn_cuda_types.h"
#include "cudautils.cuh"
#include "nbnxn_cuda_data_mgmt.h"
#include "pmalloc_cuda.h"

#define USE_CUDA_EVENT_BLOCKING_SYNC FALSE  /* makes the CPU thread busy-wait! */
/* couldomb talble size chosen such that it fits along the NB params in the texture cache */
#define EWALD_COULOMB_FORCE_TABLE_SIZE (1536)

#define MY_PI               (3.1415926535897932384626433832795)
#define TWO_OVER_SQRT_PI    (2.0/sqrt(MY_PI))
    
#define TIME_GPU_TRANSFERS 1

#define NUM_NB_KERNELS 12

/*! v1 nonbonded kernel names with name mangling. */
static const char * const nb_k1_names[NUM_NB_KERNELS] = 
{
    "_Z12k_nbnxn_rf_111cu_atomdata12cu_nb_params9cu_nblisti",
    "_Z15k_nbnxn_ewald_111cu_atomdata12cu_nb_params9cu_nblisti",
    "_Z16k_nbnxn_cutoff_111cu_atomdata12cu_nb_params9cu_nblisti",
    "_Z17k_nbnxn_rf_ener_111cu_atomdata12cu_nb_params9cu_nblisti",
    "_Z18k_nbnxn_rf_prune_111cu_atomdata12cu_nb_params9cu_nblisti",
    "_Z20k_nbnxn_ewald_ener_111cu_atomdata12cu_nb_params9cu_nblisti",
    "_Z21k_nbnxn_ewald_prune_111cu_atomdata12cu_nb_params9cu_nblisti",
    "_Z21k_nbnxn_cutoff_ener_111cu_atomdata12cu_nb_params9cu_nblisti",
    "_Z22k_nbnxn_cutoff_prune_111cu_atomdata12cu_nb_params9cu_nblisti",
    "_Z23k_nbnxn_rf_ener_prune_111cu_atomdata12cu_nb_params9cu_nblisti",
    "_Z26k_nbnxn_ewald_ener_prune_111cu_atomdata12cu_nb_params9cu_nblisti",
    "_Z27k_nbnxn_cutoff_ener_prune_111cu_atomdata12cu_nb_params9cu_nblisti"
};

/*! v2 nonbonded kernel names with name mangling. */
static const char * const nb_k2_names[NUM_NB_KERNELS] = 
{
    "_Z12k_nbnxn_rf_211cu_atomdata12cu_nb_params9cu_nblisti",
    "_Z15k_nbnxn_ewald_211cu_atomdata12cu_nb_params9cu_nblisti",
    "_Z16k_nbnxn_cutoff_211cu_atomdata12cu_nb_params9cu_nblisti",
    "_Z17k_nbnxn_rf_ener_211cu_atomdata12cu_nb_params9cu_nblisti",
    "_Z18k_nbnxn_rf_prune_211cu_atomdata12cu_nb_params9cu_nblisti",
    "_Z20k_nbnxn_ewald_ener_211cu_atomdata12cu_nb_params9cu_nblisti",
    "_Z21k_nbnxn_ewald_prune_211cu_atomdata12cu_nb_params9cu_nblisti",
    "_Z21k_nbnxn_cutoff_ener_211cu_atomdata12cu_nb_params9cu_nblisti",
    "_Z22k_nbnxn_cutoff_prune_211cu_atomdata12cu_nb_params9cu_nblisti",
    "_Z23k_nbnxn_rf_ener_prune_211cu_atomdata12cu_nb_params9cu_nblisti",
    "_Z26k_nbnxn_ewald_ener_prune_211cu_atomdata12cu_nb_params9cu_nblisti",
    "_Z27k_nbnxn_cutoff_ener_prune_211cu_atomdata12cu_nb_params9cu_nblisti"
};


/*! Dummy kernel used for sanity check. */
__device__ __global__ void k_empty_test(){}

/*! Tabulates the Ewald Coulomb force and initializes the size/scale 
    and the table GPU array. If called with an already allocated table,
    it just re-uploads the table.
 */
static void init_ewald_coulomb_force_table(cu_nb_params_t *nbp)
{
    float       *ftmp, *coul_tab;
    int         tabsize;
    double      tabscale;
    cudaError_t stat;

    tabsize     = EWALD_COULOMB_FORCE_TABLE_SIZE;
    /* Subtract 2 iso 1 to avoid access out of range due to rounding */
    tabscale    = (tabsize - 2) / sqrt(nbp->rcoulomb_sq);

    pmalloc((void**)&ftmp, tabsize*sizeof(*ftmp));

    table_spline3_fill_ewald_lr(ftmp, NULL, tabsize, tableformatF,
                                1/tabscale, nbp->ewald_beta);

    /* If the table pointer == NULL the table is generated the first time =>
       the array pointer will be saved to nb_params and the texture is bound.
     */
    coul_tab = nbp->coulomb_tab;
    if (coul_tab == NULL)
    {
        stat = cudaMalloc((void **)&coul_tab, tabsize*sizeof(*coul_tab));
        CU_RET_ERR(stat, "cudaMalloc failed on coul_tab");

        nbp->coulomb_tab = coul_tab;
        cu_bind_texture("tex_coulomb_tab", coul_tab, tabsize*sizeof(*coul_tab));
    }

    cu_copy_H2D(coul_tab, ftmp, tabsize*sizeof(*coul_tab));

    nbp->coulomb_tab_size     = tabsize;
    nbp->coulomb_tab_scale    = tabscale;

    pfree(ftmp);
}


/*! Initilizes the atomdata structure. */
static void init_atomdata(cu_atomdata_t *ad, int ntypes)
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
       need reallocation in cu_init_atomdata */
    ad->xq = NULL;
    ad->f  = NULL;

    /* size -1 indicates that the repective array hasn't been initialized yet */
    ad->natoms = -1;
    ad->nalloc = -1;
}

/*! Initilizes the nonbonded parameter data structure. */
static void init_nb_params(cu_nb_params_t *nbp,
                           const interaction_const_t *ic,
                           const nonbonded_verlet_t *nbv)
{  
    cudaError_t stat;
    int         ntypes, nnbfp; 

    ntypes  = nbv->grp[0].nbat->ntype;
    
    nbp->ewald_beta = ic->ewaldcoeff;
    nbp->epsfac     = ic->epsfac;
    nbp->two_k_rf   = 2.0 * ic->k_rf;
    nbp->c_rf       = ic->c_rf;
    nbp->rvdw_sq    = ic->rvdw * ic->rvdw;
    nbp->rcoulomb_sq= ic->rcoulomb * ic->rcoulomb;
    nbp->rlist_sq   = ic->rlist * ic->rlist;
    nbp->lj_shift   = (getenv("GMX_LJ_SHIFT") == NULL) ?
             0.0 : -1/(nbp->rvdw_sq * nbp->rvdw_sq * nbp->rvdw_sq);

    if (ic->eeltype == eelCUT)
    {
        nbp->eeltype = cu_eelCUT;
    }
    else if (EEL_RF(ic->eeltype))
    {                
        nbp->eeltype = cu_eelRF;
    }
    else if ((EEL_PME(ic->eeltype) || ic->eeltype==eelEWALD))
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
        nbp->coulomb_tab = NULL;
        init_ewald_coulomb_force_table(nbp);
    }

    nnbfp = 2*ntypes*ntypes;
    stat = cudaMalloc((void **)&nbp->nbfp, nnbfp*sizeof(*nbp->nbfp));
    CU_RET_ERR(stat, "cudaMalloc failed on nbp->nbfp"); 
    cu_copy_H2D(nbp->nbfp, nbv->grp[0].nbat->nbfp, nnbfp*sizeof(*nbp->nbfp));
    cu_bind_texture("tex_nbfp", nbp->nbfp, nnbfp*sizeof(*nbp->nbfp));
}

void reset_gpu_rlist_ewaldtab(cu_nonbonded_t cu_nb,
                              const interaction_const_t *ic)
{
    cu_nb_params_t * nbp = cu_nb->nb_params;

    nbp->rlist_sq       = ic->rlist * ic->rlist;
    nbp->rcoulomb_sq    = ic->rcoulomb * ic->rcoulomb;
    nbp->ewald_beta     = ic->ewaldcoeff;

    init_ewald_coulomb_force_table(cu_nb->nb_params);
}

/*! Initilizes the neighborlist data structure. */
static void init_nblist(cu_nblist_t *nbl)
{
    /* initilize to NULL poiters to data that is not allocated here and will
       need reallocation in cu_init_atomdata */
    nbl->sci     = NULL;
    nbl->cj4     = NULL;
    nbl->excl    = NULL;    
    
    /* size -1 indicates that the repective array hasn't been initialized yet */
    nbl->na_c        = -1;
    nbl->nsci        = -1;
    nbl->sci_nalloc  = -1;
    nbl->ncj4        = -1;
    nbl->cj4_nalloc  = -1;
    nbl->nexcl       = -1;
    nbl->excl_nalloc = -1;
    nbl->prune_nbl   = FALSE;
}

/*! Initilizes the timer data structure. */
static void init_timers(cu_timers_t *t, gmx_bool bDomDec)
{
    cudaError_t stat;
    /* XXX */ 
    int eventflags = ( USE_CUDA_EVENT_BLOCKING_SYNC ? cudaEventBlockingSync: cudaEventDefault );

    stat = cudaEventCreateWithFlags(&(t->start_atdat), eventflags);
    CU_RET_ERR(stat, "cudaEventCreate on start_atdat failed");
    stat = cudaEventCreateWithFlags(&(t->stop_atdat), eventflags);
    CU_RET_ERR(stat, "cudaEventCreate on stop_atdat failed");

    /* The non-local counters/stream (second in the array) are needed only with DD. */
    for (int i = 0; i <= bDomDec ? 1 : 0; i++)
    {
        stat = cudaEventCreateWithFlags(&(t->start_nb_k[i]), eventflags);
        CU_RET_ERR(stat, "cudaEventCreate on start_nb_k failed");
        stat = cudaEventCreateWithFlags(&(t->stop_nb_k[i]), eventflags);
        CU_RET_ERR(stat, "cudaEventCreate on stop_nb_k failed");


        stat = cudaEventCreateWithFlags(&(t->start_nbl_h2d[i]), eventflags);
        CU_RET_ERR(stat, "cudaEventCreate on start_nbl_h2d failed");
        stat = cudaEventCreateWithFlags(&(t->stop_nbl_h2d[i]), eventflags);
        CU_RET_ERR(stat, "cudaEventCreate on stop_nbl_h2d failed");

        stat = cudaEventCreateWithFlags(&(t->start_nb_h2d[i]), eventflags);
        CU_RET_ERR(stat, "cudaEventCreate on start_nb_h2d failed");
        stat = cudaEventCreateWithFlags(&(t->stop_nb_h2d[i]), eventflags);
        CU_RET_ERR(stat, "cudaEventCreate on stop_nb_h2d failed");

        stat = cudaEventCreateWithFlags(&(t->start_nb_d2h[i]), eventflags);
        CU_RET_ERR(stat, "cudaEventCreate on start_nb_d2h failed");
        stat = cudaEventCreateWithFlags(&(t->stop_nb_d2h[i]), eventflags);
        CU_RET_ERR(stat, "cudaEventCreate on stop_nb_d2h failed");
    }
}

/*! Initilizes the timings data structure. */
static void init_timings(cu_timings_t *t)
{
    int i, j;

    t->nb_h2d_time = 0.0;
    t->nb_d2h_time = 0.0;
    t->nb_count    = 0;
    t->nbl_h2d_time = 0.0;
    t->nbl_h2d_count = 0;
    for (i = 0; i < 2; i++)
    {
        for(j = 0; j < 2; j++)
        {
            t->k_time[i][j].t = 0.0;
            t->k_time[i][j].c = 0;
        }
    }
}

void nbnxn_cuda_init(FILE *fplog,
                     cu_nonbonded_t *p_cu_nb,
                     gmx_bool bDomDec)
{
    cudaError_t stat;
    cu_nonbonded_t  nb;

    if (p_cu_nb == NULL) return;

    snew(nb, 1); 
    snew(nb->dev_info, 1);
    snew(nb->atomdata, 1); 
    snew(nb->nb_params, 1); 
    snew(nb->nblist[eintLocal], 1);
    if (bDomDec)
    {
        snew(nb->nblist[eintNonlocal], 1);
    }

    /* CUDA event timers don't work with multiple streams so 
       we have to disable timing with DD */
    nb->do_time = (!bDomDec && (getenv("GMX_DISABLE_CUDA_TIMING") == NULL));
    snew(nb->timers, 1); 
    snew(nb->timings, 1); 

    /* init tmpdata */
    pmalloc((void**)&nb->tmpdata.e_lj, sizeof(*nb->tmpdata.e_lj));
    pmalloc((void**)&nb->tmpdata.e_el, sizeof(*nb->tmpdata.e_el));
    pmalloc((void**)&nb->tmpdata.f_shift, SHIFTS * sizeof(*nb->tmpdata.f_shift));

    init_nblist(nb->nblist[eintLocal]);
    stat = cudaStreamCreate(&nb->stream[eintLocal]);
    CU_RET_ERR(stat, "cudaStreamCreate on stream[eintLocal] failed");
    if (bDomDec)
    {
        init_nblist(nb->nblist[eintNonlocal]);
        stat = cudaStreamCreate(&nb->stream[eintNonlocal]);
        CU_RET_ERR(stat, "cudaStreamCreate on stream[eintNonlocal] failed");
    }

    init_timers(nb->timers, bDomDec);
    init_timings(nb->timings);

    /* init device info */
    /* FIXME this should not be done here! */
    stat = cudaGetDevice(&nb->dev_info->dev_id);
    CU_RET_ERR(stat, "cudaGetDevice failed");
    stat = cudaGetDeviceProperties(&nb->dev_info->dev_prop, nb->dev_info->dev_id);
    CU_RET_ERR(stat, "cudaGetDeviceProperties failed");

    /* If ECC is enabled cudaStreamSynchronize introduces huge idling so we'll 
       switch to the (atmittedly fragile) memory polling waiting to preserve 
       performance. Event-timing also needs to be disabled. */
    nb->use_stream_sync = TRUE;
    if (nb->dev_info->dev_prop.ECCEnabled == 1 ||
        (getenv("GMX_NO_CUDA_STREAMSYNC") != NULL))
    {
        nb->use_stream_sync = FALSE;
        nb->do_time         = FALSE;
    }

    *p_cu_nb = nb;

    if (debug)
    {
        fprintf(debug, "Initialized CUDA data structures.\n");
    }

    /* k_nbnxn_*_1 48/16 kB Shared/L1 */
    for (int i = 0; i < NUM_NB_KERNELS; i++)
    {
        stat = cudaFuncSetCacheConfig(nb_k1_names[i],  cudaFuncCachePreferShared);
        CU_RET_ERR(stat, "cudaFuncSetCacheConfig failed");
    }

    /* k_nbnxn_*_2 16/48 kB Shared/L1 */
    for (int i = 0; i < NUM_NB_KERNELS; i++)
    {
        stat = cudaFuncSetCacheConfig(nb_k2_names[i], cudaFuncCachePreferL1);
        CU_RET_ERR(stat, "cudaFuncSetCacheConfig failed");
    }

    /* TODO: move this to gpu_utils module */
    k_empty_test<<<1, 512>>>();
    CU_LAUNCH_ERR_SYNC("dummy test kernel");
}

void nbnxn_cuda_init_const(FILE *fplogi,
                           cu_nonbonded_t cu_nb,
                           const interaction_const_t *ic,
                           const nonbonded_verlet_t *nbv)
{
    init_atomdata(cu_nb->atomdata, nbv->grp[0].nbat->ntype);
    init_nb_params(cu_nb->nb_params, ic, nbv);

    /* clear energy and shift force outputs */
    nbnxn_cuda_clear_e_fshift(cu_nb);
}

void nbnxn_cuda_init_pairlist(cu_nonbonded_t cu_nb,
                              const nbnxn_pairlist_t *h_nblist,
                              int iloc)
{
    char         sbuf[STRLEN];
    cudaError_t  stat;
    gmx_bool     do_time    = cu_nb->do_time;
    cudaStream_t stream     = cu_nb->stream[iloc];
    cu_nblist_t  *d_nblist  = cu_nb->nblist[iloc];

    if (d_nblist->na_c < 0)
    {
        d_nblist->na_c = h_nblist->na_c;
    }
    else
    {
        if (d_nblist->na_c != h_nblist->na_c)
        {
            sprintf(sbuf, "In cu_init_nblist: the #atoms per cell has changed (from %d to %d)",
                    d_nblist->na_c, h_nblist->na_c);
            gmx_incons(sbuf);
        }
    }

    if (do_time)
    {
        stat = cudaEventRecord(cu_nb->timers->start_nbl_h2d[iloc], stream);
        CU_RET_ERR(stat, "cudaEventRecord failed");
    }

    cu_realloc_buffered((void **)&d_nblist->sci, h_nblist->sci, sizeof(*d_nblist->sci),
                         &d_nblist->nsci, &d_nblist->sci_nalloc,
                         h_nblist->nsci,
                         stream, TRUE);

    cu_realloc_buffered((void **)&d_nblist->cj4, h_nblist->cj4, sizeof(*d_nblist->cj4),
                         &d_nblist->ncj4, &d_nblist->cj4_nalloc,
                         h_nblist->ncj4,
                         stream, TRUE);

    cu_realloc_buffered((void **)&d_nblist->excl, h_nblist->excl, sizeof(*d_nblist->excl),
                         &d_nblist->nexcl, &d_nblist->excl_nalloc,
                         h_nblist->nexcl, 
                         stream, TRUE);

    if (do_time)
    {
        stat = cudaEventRecord(cu_nb->timers->stop_nbl_h2d[iloc], stream);
        CU_RET_ERR(stat, "cudaEventRecord failed");
    }

    /* need to prune the neighbor list during the next step */
    d_nblist->prune_nbl = TRUE;
}

void nbnxn_cuda_upload_shiftvec(cu_nonbonded_t cu_nb,
                                const nbnxn_atomdata_t *nbatom)
{
    cu_atomdata_t   *adat = cu_nb->atomdata;

    /* only if we have a dynamic box */
    if (nbatom->dynamic_box || !adat->shift_vec_copied)
    {
        cu_copy_H2D_async(adat->shift_vec, nbatom->shift_vec, 
                          SHIFTS * sizeof(*adat->shift_vec), 0);
        adat->shift_vec_copied = TRUE;
    }
}


/*! Clears nonbonded force output array on the GPU - muldule internal 
    implementation that takes the number of atoms to clear the output for. */
static void nbnxn_cuda_clear_f(cu_nonbonded_t cu_nb, int natoms_clear)
{
    cudaError_t stat;
    cu_atomdata_t *adat = cu_nb->atomdata;

    stat = cudaMemsetAsync(adat->f, 0, natoms_clear * sizeof(*adat->f), 0);
    CU_RET_ERR(stat, "cudaMemsetAsync on f falied");
}

void nbnxn_cuda_clear_f(cu_nonbonded_t cu_nb)
{
    nbnxn_cuda_clear_f(cu_nb, cu_nb->atomdata->natoms);
}

void nbnxn_cuda_clear_e_fshift(cu_nonbonded_t cu_nb)
{
    cudaError_t stat;    
    cu_atomdata_t *adat = cu_nb->atomdata;

    stat = cudaMemsetAsync(adat->f_shift, 0, SHIFTS * sizeof(*adat->f_shift), 0);
    CU_RET_ERR(stat, "cudaMemsetAsync on f_shift falied");
    stat = cudaMemsetAsync(adat->e_lj, 0, sizeof(*adat->e_lj), 0);
    CU_RET_ERR(stat, "cudaMemsetAsync on e_lj falied");
    stat = cudaMemsetAsync(adat->e_el, 0, sizeof(*adat->e_el), 0);
    CU_RET_ERR(stat, "cudaMemsetAsync on e_el falied");
}

/* TODO: add gmx over_alloc call */
void nbnxn_cuda_init_atomdata(cu_nonbonded_t cu_nb,
                              const nbnxn_atomdata_t *nbat)
{
    cudaError_t stat;
    int         nalloc, natoms;
    gmx_bool    realloced = FALSE;
    gmx_bool    do_time     = cu_nb->do_time;
    cu_timers_t *timers     = cu_nb->timers;
    cu_atomdata_t *d_atomd  = cu_nb->atomdata;

    natoms = nbat->natoms;

    if (do_time)
    {
        /* time async copy */
        stat = cudaEventRecord(timers->start_atdat, 0);
        CU_RET_ERR(stat, "cudaEventRecord failed");
    }

    /* need to reallocate if we have to copy more atoms than the amount of space
       available and only allocate if we haven't initilzed yet, i.e d_atomd->natoms == -1 */
    if (natoms > d_atomd->nalloc)
    {
        nalloc = natoms * 1.2 + 100;
    
        /* free up first if the arrays have already been initialized */
        if (d_atomd->nalloc != -1)
        {
            cu_free_buffered(d_atomd->f, &d_atomd->natoms, &d_atomd->nalloc);
            cu_free_buffered(d_atomd->xq);
            cu_free_buffered(d_atomd->atom_types);
        }
        
        stat = cudaMalloc((void **)&d_atomd->f, nalloc*sizeof(*d_atomd->f));
        CU_RET_ERR(stat, "cudaMalloc failed on d_atomd->f");                   
        stat = cudaMalloc((void **)&d_atomd->xq, nalloc*sizeof(*d_atomd->xq));
        CU_RET_ERR(stat, "cudaMalloc failed on d_atomd->xq");     

        stat = cudaMalloc((void **)&d_atomd->atom_types, nalloc*sizeof(*d_atomd->atom_types));
        CU_RET_ERR(stat, "cudaMalloc failed on d_atomd->atom_types"); 

        d_atomd->nalloc = nalloc;
        realloced = TRUE;
    }

    d_atomd->natoms = natoms;
    d_atomd->natoms_local = nbat->natoms_local;

    /* need to clear GPU f output if realloc happened */
    if (realloced)
    {
        nbnxn_cuda_clear_f(cu_nb, nalloc);
    }

    cu_copy_H2D_async(d_atomd->atom_types, nbat->type,
                      natoms*sizeof(*d_atomd->atom_types), 0);

    if (do_time)
    {
        stat = cudaEventRecord(timers->stop_atdat, 0);
        CU_RET_ERR(stat, "cudaEventRecord failed");
    }
}

void nbnxn_cuda_free(FILE *fplog, cu_nonbonded_t cu_nb, gmx_bool bDomDec)
{
    cudaError_t     stat;
    cu_atomdata_t   *atomdata;
    cu_nb_params_t  *nb_params;
    cu_nblist_t     *nblist, *nblist_nl;
    cu_timers_t     *timers;

    atomdata    = cu_nb->atomdata;
    nb_params   = cu_nb->nb_params;
    nblist      = cu_nb->nblist[eintLocal];
    nblist_nl   = cu_nb->nblist[eintNonlocal];
    timers      = cu_nb->timers;

    if (cu_nb == NULL) return;

    if (nb_params->eeltype == cu_eelEWALD)
    {
        cu_unbind_texture("tex_coulomb_tab");
        cu_free_buffered(nb_params->coulomb_tab, &nb_params->coulomb_tab_size);
    }

    stat = cudaEventDestroy(timers->start_atdat);
    CU_RET_ERR(stat, "cudaEventDestroy failed on timers->start_atdat");
    stat = cudaEventDestroy(timers->stop_atdat);
    CU_RET_ERR(stat, "cudaEventDestroy failed on timers->stop_atdat");

    /* The non-local counters/stream (second in the array) are needed only with DD. */
    for (int i = 0; i <= bDomDec ? 1 : 0; i++)
    {

        stat = cudaEventDestroy(timers->start_nb_k[i]);
        CU_RET_ERR(stat, "cudaEventDestroy failed on timers->start_nb_k");
        stat = cudaEventDestroy(timers->stop_nb_k[i]);
        CU_RET_ERR(stat, "cudaEventDestroy failed on timers->stop_nb_k");

        stat = cudaEventDestroy(timers->start_nbl_h2d[i]);
        CU_RET_ERR(stat, "cudaEventDestroy failed on timers->start_nbl_h2d");
        stat = cudaEventDestroy(timers->stop_nbl_h2d[i]);
        CU_RET_ERR(stat, "cudaEventDestroy failed on timers->stop_nbl_h2d");

        stat = cudaStreamDestroy(cu_nb->stream[i]);
        CU_RET_ERR(stat, "cudaStreamDestroy failed on stream");

        stat = cudaEventDestroy(timers->start_nb_h2d[i]);
        CU_RET_ERR(stat, "cudaEventDestroy failed on timers->start_nb_h2d");
        stat = cudaEventDestroy(timers->stop_nb_h2d[i]);
        CU_RET_ERR(stat, "cudaEventDestroy failed on timers->stop_nb_h2d");

        stat = cudaEventDestroy(timers->start_nb_d2h[i]);
        CU_RET_ERR(stat, "cudaEventDestroy failed on timers->start_nb_d2h");
        stat = cudaEventDestroy(timers->stop_nb_d2h[i]);
        CU_RET_ERR(stat, "cudaEventDestroy failed on timers->stop_nb_d2h");
    }

    cu_unbind_texture("tex_nbfp");
    cu_free_buffered(nb_params->nbfp);

    stat = cudaFree(atomdata->shift_vec);
    CU_RET_ERR(stat, "cudaEventDestroy failed on atomdata->shift_vec");
    stat = cudaFree(atomdata->f_shift);
    CU_RET_ERR(stat, "cudaEventDestroy failed on atomdata->f_shift");

    stat = cudaFree(atomdata->e_lj);
    CU_RET_ERR(stat, "cudaEventDestroy failed on atomdata->e_lj");
    stat = cudaFree(atomdata->e_el);
    CU_RET_ERR(stat, "cudaEventDestroy failed on atomdata->e_el");

    cu_free_buffered(atomdata->f, &atomdata->natoms, &atomdata->nalloc);
    cu_free_buffered(atomdata->xq);
    cu_free_buffered(atomdata->atom_types, &atomdata->ntypes);

    cu_free_buffered(nblist->sci, &nblist->nsci, &nblist->sci_nalloc);
    cu_free_buffered(nblist->cj4, &nblist->ncj4, &nblist->cj4_nalloc);
    cu_free_buffered(nblist->excl, &nblist->nexcl, &nblist->excl_nalloc);
    if (bDomDec)
    {
        cu_free_buffered(nblist_nl->sci, &nblist_nl->nsci, &nblist_nl->sci_nalloc);
        cu_free_buffered(nblist_nl->cj4, &nblist_nl->ncj4, &nblist_nl->cj4_nalloc);
        cu_free_buffered(nblist_nl->excl, &nblist_nl->nexcl, &nblist->excl_nalloc);
    }

    stat = cudaThreadExit();
    CU_RET_ERR(stat, "cudaThreadExit failed");

    if (debug)
    {
        fprintf(debug, "Cleaned up CUDA data structures.\n");
    }
}

void cu_synchstream_atomdata(cu_nonbonded_t cu_nb, int iloc)
{
    cudaError_t stat;
    cudaStream_t stream = cu_nb->stream[iloc];

    stat = cudaStreamWaitEvent(stream, cu_nb->timers->stop_atdat, 0);
    CU_RET_ERR(stat, "cudaStreamWaitEvent failed");
}

cu_timings_t * nbnxn_cuda_get_timings(cu_nonbonded_t cu_nb)
{
    return (cu_nb != NULL && cu_nb->do_time) ? cu_nb->timings : NULL;
}

void nbnxn_cuda_reset_timings(cu_nonbonded_t cu_nb)
{
    if (cu_nb->do_time)
    {
        init_timings(cu_nb->timings);
    }
}

int nbnxn_cuda_min_ci_balanced(cu_nonbonded_t cu_nb)
{
    return cu_nb != NULL ? 
        GPU_MIN_CI_BALANCED_FACTOR*cu_nb->dev_info->dev_prop.multiProcessorCount : 0;
}


/****** FIXME Old stuff, mostly deprecated, remove before release  *****/

/* Upload asynchronously to the GPU the coordinate+charge array.
 * XXX not used  
 */
void cu_move_xq(cu_nonbonded_t cu_nb, const nbnxn_atomdata_t *nbat, int aloc)
{
    int iloc = -1; 

    /* determine interaction locality from atom locality 
       (needed for indexing timers/streams) */
    if (LOCAL_A(aloc))
    {
        iloc = eintLocal;
    }
    else if (NONLOCAL_A(aloc))
    {
        iloc = eintNonlocal;
    }
    else 
    {
        gmx_incons("Invalid atom locality passed (valid here is only local or nonlocal)");
    }

    cu_atomdata_t   *d_nbat = cu_nb->atomdata;
    cudaStream_t    stream  = cu_nb->stream[iloc];

    cu_copy_H2D_async(d_nbat->xq, nbat->x,
                        d_nbat->natoms * sizeof(*d_nbat->xq), stream);
}

/*! Waits until the atom data gets copied to the GPU and times the transfer.
 *  XXX not used  
 */
void cu_wait_atomdata(cu_nonbonded_t cu_nb)
{
    float t;
    cu_wait_event(cu_nb->timers->stop_atdat, cu_nb->timers->start_atdat, &t);
    cu_nb->timings->nbl_h2d_time += t;
}
