#include <stdlib.h>
#include <stdio.h>

#include "gmx_fatal.h"
#include "smalloc.h"
#include "force.h"
#include "nb_verlet.h"
#include "interaction_const.h"

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
                                 cudaStream_t stream, 
                                 gmx_bool doAsync);

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
    tabscale    = (tabsize - 1) / sqrt(nbp->rcoulomb_sq);

    pmalloc((void**)&ftmp, tabsize*sizeof(*ftmp));

    table_spline3_fill_ewald(ftmp, tabsize, tableformatF,
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

    upload_cudata(coul_tab, ftmp, tabsize*sizeof(*coul_tab));

    nbp->coulomb_tab_size     = tabsize;
    nbp->coulomb_tab_scale    = tabscale;

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
void init_nb_params(cu_nb_params_t *nbp,
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
    upload_cudata(nbp->nbfp, nbv->grp[0].nbat->nbfp, nnbfp*sizeof(*nbp->nbfp));
    cu_bind_texture("tex_nbfp", nbp->nbfp, nnbfp*sizeof(*nbp->nbfp));
}

void reset_cu_rlist_ewaldtab(cu_nonbonded_t cu_nb,
                             const interaction_const_t *ic)
{
    cu_nb_params_t * nbp = cu_nb->nb_params;

    nbp->rlist_sq       = ic->rlist * ic->rlist;
    nbp->rcoulomb_sq    = ic->rcoulomb * ic->rcoulomb;
    nbp->ewald_beta     = ic->ewaldcoeff;

    init_ewald_coulomb_force_table(cu_nb->nb_params);
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

/*! Initilizes the timer data structure. */
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

void init_cu_nonbonded(FILE *fplog,
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

    *p_cu_nb = nb;

    if (debug)
    {
        fprintf(debug, "Initialized CUDA data structures.\n");
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

    /* TODO: move this to gpu_utils module */
    k_empty_test<<<1, 512>>>();
    CU_LAUNCH_ERR_SYNC("test kernel");
}

/*! Initilizes force-field related data (called only once in the beginning).
 */
void init_cudata_ff(FILE *fplogi,
                    cu_nonbonded_t cu_nb,
                    const interaction_const_t *ic,
                    const nonbonded_verlet_t *nbv)
{
    init_atomdata(cu_nb->atomdata, nbv->grp[0].nbat->ntype);
    init_nb_params(cu_nb->nb_params, ic, nbv);

    /* clear energy and shift force outputs */
    cu_clear_nb_e_fs_out(cu_nb);
}

/*! Initilizes neighbor list on the GPU, called at every neighbor search step. 
 */
void init_cudata_nblist(cu_nonbonded_t cu_nb, 
                        const gmx_nblist_t *h_nblist,
                        int iloc)
{
    char sbuf[STRLEN];
    cudaError_t  stat;
    gmx_bool     do_time    = cu_nb->do_time;
    cudaStream_t stream     = cu_nb->stream[iloc];
    cu_nblist_t  *d_nblist  = cu_nb->nblist[iloc];

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

    if (do_time)
    {
        stat = cudaEventRecord(cu_nb->timers->start_nbl_h2d[iloc], stream);
        CU_RET_ERR(stat, "cudaEventRecord failed");
    }

    realloc_cudata_array((void **)&d_nblist->ci, h_nblist->ci, sizeof(*(d_nblist->ci)),
                         &d_nblist->nci, &d_nblist->ci_nalloc,
                         h_nblist->nci,
                         stream, TRUE);

    realloc_cudata_array((void **)&d_nblist->sj4, h_nblist->sj4, sizeof(*(d_nblist->sj4)),
                         &d_nblist->nsj4, &d_nblist->sj4_nalloc,
                         h_nblist->nsj4,
                         stream, TRUE);

    realloc_cudata_array((void **)&d_nblist->excl, h_nblist->excl, sizeof(*(d_nblist->excl)),
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

static void cu_clear_nb_f_out(cu_nonbonded_t cu_nb, int natoms_clear)
{
    cudaError_t stat;

    cu_atomdata_t   *adat = cu_nb->atomdata;

    stat = cudaMemsetAsync(adat->f, 0, natoms_clear * sizeof(*adat->f), 0);
    CU_RET_ERR(stat, "cudaMemsetAsync on f falied");
}

void cu_clear_nb_f_out(cu_nonbonded_t cu_nb)
{
    cu_clear_nb_f_out(cu_nb, cu_nb->atomdata->natoms);
}

void cu_clear_nb_e_fs_out(cu_nonbonded_t cu_nb)
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

/*! Initilizes atom-data on the GPU, called at every neighbor search step. 
    TODO: add gmx over_alloc call
 */
void init_cudata_atoms(cu_nonbonded_t cu_nb,
                       const gmx_nb_atomdata_t *nbat)
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
            destroy_cudata_array(d_atomd->f, &d_atomd->natoms, &d_atomd->nalloc);
            destroy_cudata_array(d_atomd->xq);
            destroy_cudata_array(d_atomd->atom_types);
        }
        
        stat = cudaMalloc((void **)&d_atomd->f, nalloc*sizeof(*(d_atomd->f)));
        CU_RET_ERR(stat, "cudaMalloc failed on d_atomd->f");                   
        stat = cudaMalloc((void **)&d_atomd->xq, nalloc*sizeof(*(d_atomd->xq)));
        CU_RET_ERR(stat, "cudaMalloc failed on d_atomd->xq");     

        stat = cudaMalloc((void **)&d_atomd->atom_types, nalloc*sizeof(*(d_atomd->atom_types)));
        CU_RET_ERR(stat, "cudaMalloc failed on d_atomd->atom_types"); 

        d_atomd->nalloc = nalloc;
        realloced = TRUE;
    }

    d_atomd->natoms = natoms;
    d_atomd->natoms_local = nbat->natoms_local;

    /* need to clear GPU f output if realloc happened */
    if (realloced)
    {
        cu_clear_nb_f_out(cu_nb, nalloc);
    }

    upload_cudata_async(d_atomd->atom_types, nbat->type,
                        natoms*sizeof(*d_atomd->atom_types), 0);

    if (do_time)
    {
        stat = cudaEventRecord(timers->stop_atdat, 0);
        CU_RET_ERR(stat, "cudaEventRecord failed");
    }
}

/*! Frees up all GPU resources used for the nonbonded calculations. */
void destroy_cudata(FILE *fplog, cu_nonbonded_t cu_nb, gmx_bool bDomDec)
{
    cudaError_t stat;
    cu_atomdata_t       *atomdata;
    cu_nb_params_t      *nb_params;
    cu_nblist_t         *nblist, *nblist_nl;
    cu_timers_t         *timers;

    atomdata    = cu_nb->atomdata;
    nb_params   = cu_nb->nb_params;
    nblist      = cu_nb->nblist[eintLocal];
    nblist_nl   = cu_nb->nblist[eintNonlocal];
    timers      = cu_nb->timers;

    if (cu_nb == NULL) return;

    if (nb_params->eeltype == cu_eelEWALD)
    {
        cu_unbind_texture("tex_coulomb_tab");
        destroy_cudata_array(nb_params->coulomb_tab, &nb_params->coulomb_tab_size);            
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
    if (bDomDec)
    {
        destroy_cudata_array(nblist_nl->ci, &nblist_nl->nci, &nblist_nl->ci_nalloc);
        destroy_cudata_array(nblist_nl->sj4, &nblist_nl->nsj4, &nblist_nl->sj4_nalloc);
        destroy_cudata_array(nblist_nl->excl, &nblist_nl->nexcl, &nblist->excl_nalloc);
    }

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
                                 cudaStream_t stream,
                                 gmx_bool doAsync)
{
    cudaError_t stat;

    if (d_dest == NULL || req_size < 0)
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

        *curr_alloc_size = 1.2 * req_size + 100;

        stat = cudaMalloc(d_dest, *curr_alloc_size * type_size);
        CU_RET_ERR(stat, "cudaMalloc failed in realloc_cudata_array");
    }

    /* size could have changed without actual reallocation */
    *curr_size = req_size;

    /* upload to device */
    if (h_src)
    {
        if (doAsync)
        {
            upload_cudata_async(*d_dest, h_src, *curr_size * type_size, stream);
        }
        else 
        {
            upload_cudata(*d_dest, h_src,  *curr_size * type_size);
        }
    }
}

/* XXX not used */
void cu_move_xq(cu_nonbonded_t cu_nb, const gmx_nb_atomdata_t *nbat,
                int aloc)
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

    upload_cudata_async(d_nbat->xq, nbat->x,
                        d_nbat->natoms * sizeof(*d_nbat->xq), stream);
}

/*! Blocking waits until the atom data gets copied to the GPU and times the transfer.
 * XXX not used  
 */
void cu_wait_atomdata(cu_nonbonded_t cu_nb)
{
    float t;
    cu_blockwait_event(cu_nb->timers->stop_atdat, cu_nb->timers->start_atdat, &t);
    cu_nb->timings->nbl_h2d_time += t;
}

/*! Synchronizes the respective stream with the atomdata init operation.
 */
void cu_synchstream_atomdata(cu_nonbonded_t cu_nb, int iloc)
{
    cudaError_t stat;
    cudaStream_t stream = cu_nb->stream[iloc];

    stat = cudaStreamWaitEvent(stream, cu_nb->timers->stop_atdat, 0);
    CU_RET_ERR(stat, "cudaStreamWaitEvent failed");
}

/*! Returns the GPU timings structure or NULL if cu_nb is NULL or timing 
    is turned off. 
 */
cu_timings_t * get_gpu_timings(cu_nonbonded_t cu_nb)
{
    return (cu_nb != NULL && cu_nb->do_time) ? cu_nb->timings : NULL;
}

/*! Resets GPU timers. */
void reset_gpu_timings(cu_nonbonded_t cu_nb)
{
    if (cu_nb->do_time)
    {
        init_timings(cu_nb->timings);
    }
}

int cu_calc_min_ci_balanced(cu_nonbonded_t cu_nb)
{
    return cu_nb != NULL ? 
        GPU_MIN_CI_BALANCED_FACTOR*cu_nb->dev_info->dev_prop.multiProcessorCount : 0;
}


/*** Old stuff ***/
/* XXX not used */
int cu_upload_X(cu_nonbonded_t cu_nb, real *h_x) 
{
    cu_atomdata_t *ad = cu_nb->atomdata;

    return upload_cudata(ad->xq, h_x, ad->natoms*sizeof(*ad->xq));
}

/* XXX not used */
int cu_download_F(real *h_f, cu_nonbonded_t cu_nb)
{
    cu_atomdata_t *ad = cu_nb->atomdata;

    return download_cudata(h_f, ad->f, ad->natoms*sizeof(*ad->f));
}
