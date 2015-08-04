/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2012,2013,2014,2015, by the GROMACS development team, led by
 * Mark Abraham, David van der Spoel, Berk Hess, and Erik Lindahl,
 * and including many others, as listed in the AUTHORS file in the
 * top-level source directory and at http://www.gromacs.org.
 *
 * GROMACS is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public License
 * as published by the Free Software Foundation; either version 2.1
 * of the License, or (at your option) any later version.
 *
 * GROMACS is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with GROMACS; if not, see
 * http://www.gnu.org/licenses, or write to the Free Software Foundation,
 * Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA.
 *
 * If you want to redistribute modifications to GROMACS, please
 * consider that scientific software is very special. Version
 * control is crucial - bugs must be traceable. We will be happy to
 * consider code for inclusion in the official distribution, but
 * derived work must not be called official GROMACS. Details are found
 * in the README & COPYING files - if they are missing, get the
 * official version at http://www.gromacs.org.
 *
 * To help us fund GROMACS development, we humbly ask that you cite
 * the research papers on the package. Check out http://www.gromacs.org.
 */
/*! \file
 *  \brief Define CUDA implementation of nbnxn_gpu_data_mgmt.h
 *
 *  \author Szilard Pall <pall.szilard@gmail.com>
 */
#include "gmxpre.h"

#include "config.h"

#include <assert.h>
#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>

#include <cuda_profiler_api.h>

#include "gromacs/gmxlib/cuda_tools/cudautils.cuh"
#include "gromacs/gmxlib/cuda_tools/pmalloc_cuda.h"
#include "gromacs/gmxlib/gpu_utils/gpu_utils.h"
#include "gromacs/legacyheaders/gmx_detect_hardware.h"
#include "gromacs/legacyheaders/typedefs.h"
#include "gromacs/legacyheaders/types/enums.h"
#include "gromacs/legacyheaders/types/force_flags.h"
#include "gromacs/legacyheaders/types/interaction_const.h"
#include "gromacs/mdlib/nb_verlet.h"
#include "gromacs/mdlib/nbnxn_consts.h"
#include "gromacs/mdlib/nbnxn_gpu_data_mgmt.h"
#include "gromacs/pbcutil/ishift.h"
#include "gromacs/timing/gpu_timing.h"
#include "gromacs/utility/basedefinitions.h"
#include "gromacs/utility/cstringutil.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/smalloc.h"

#include "nbnxn_cuda_types.h"

static bool bUseCudaEventBlockingSync = false; /* makes the CPU thread block */

/* This is a heuristically determined parameter for the Fermi, Kepler
 * and Maxwell architectures for the minimum size of ci lists by multiplying
 * this constant with the # of multiprocessors on the current device.
 * Since the maximum number of blocks per multiprocessor is 16, the ideal
 * count for small systems is 32 or 48 blocks per multiprocessor. Because
 * there is a bit of fluctuations in the generated block counts, we use
 * a target of 44 instead of the ideal value of 48.
 */
static unsigned int gpu_min_ci_balanced_factor = 44;

/* Functions from nbnxn_cuda.cu */
extern void nbnxn_cuda_set_cacheconfig(gmx_device_info_t *devinfo);
extern const struct texture<float, 1, cudaReadModeElementType> &nbnxn_cuda_get_nbfp_texref();
extern const struct texture<float, 1, cudaReadModeElementType> &nbnxn_cuda_get_nbfp_comb_texref();
extern const struct texture<float, 1, cudaReadModeElementType> &nbnxn_cuda_get_coulomb_tab_texref();

/* We should actually be using md_print_warn in md_logging.c,
 * but we can't include mpi.h in CUDA code.
 */
static void md_print_warn(FILE       *fplog,
                          const char *fmt, ...)
{
    va_list ap;

    if (fplog != NULL)
    {
        /* We should only print to stderr on the master node,
         * in most cases fplog is only set on the master node, so this works.
         */
        va_start(ap, fmt);
        fprintf(stderr, "\n");
        vfprintf(stderr, fmt, ap);
        fprintf(stderr, "\n");
        va_end(ap);

        va_start(ap, fmt);
        fprintf(fplog, "\n");
        vfprintf(fplog, fmt, ap);
        fprintf(fplog, "\n");
        va_end(ap);
    }
}


/* Fw. decl. */
static void nbnxn_cuda_clear_e_fshift(gmx_nbnxn_cuda_t *nb);

/* Fw. decl, */
static void nbnxn_cuda_free_nbparam_table(cu_nbparam_t            *nbparam,
                                          const gmx_device_info_t *dev_info);



#ifdef HAVE_CUDA_TEXOBJ_SUPPORT
static bool use_texobj(const gmx_device_info_t *dev_info)
{
    /* Only device CC >= 3.0 (Kepler and later) support texture objects */
    return (dev_info->prop.major >= 3);
}
#endif

/*! Tabulates the Ewald Coulomb force and initializes the size/scale
    and the table GPU array. If called with an already allocated table,
    it just re-uploads the table.
 */
static void init_ewald_coulomb_force_table(const interaction_const_t *ic,
                                           cu_nbparam_t              *nbp,
                                           const gmx_device_info_t   *dev_info)
{
    float       *coul_tab;
    cudaError_t  stat;

    if (nbp->coulomb_tab != NULL)
    {
        nbnxn_cuda_free_nbparam_table(nbp, dev_info);
    }

    stat = cudaMalloc((void **)&coul_tab, ic->tabq_size*sizeof(*coul_tab));
    CU_RET_ERR(stat, "cudaMalloc failed on coul_tab");

    nbp->coulomb_tab = coul_tab;

#ifdef HAVE_CUDA_TEXOBJ_SUPPORT
    /* Only device CC >= 3.0 (Kepler and later) support texture objects */
    if (use_texobj(dev_info))
    {
        cudaResourceDesc rd;
        memset(&rd, 0, sizeof(rd));
        rd.resType                  = cudaResourceTypeLinear;
        rd.res.linear.devPtr        = nbp->coulomb_tab;
        rd.res.linear.desc.f        = cudaChannelFormatKindFloat;
        rd.res.linear.desc.x        = 32;
        rd.res.linear.sizeInBytes   = ic->tabq_size*sizeof(*coul_tab);

        cudaTextureDesc td;
        memset(&td, 0, sizeof(td));
        td.readMode                 = cudaReadModeElementType;
        stat = cudaCreateTextureObject(&nbp->coulomb_tab_texobj, &rd, &td, NULL);
        CU_RET_ERR(stat, "cudaCreateTextureObject on coulomb_tab_texobj failed");
    }
    else
#endif  /* HAVE_CUDA_TEXOBJ_SUPPORT */
    {
        GMX_UNUSED_VALUE(dev_info);
        cudaChannelFormatDesc cd   = cudaCreateChannelDesc<float>();
        stat = cudaBindTexture(NULL, &nbnxn_cuda_get_coulomb_tab_texref(),
                               coul_tab, &cd,
                               ic->tabq_size*sizeof(*coul_tab));
        CU_RET_ERR(stat, "cudaBindTexture on coulomb_tab_texref failed");
    }

    cu_copy_H2D(coul_tab, ic->tabq_coul_F, ic->tabq_size*sizeof(*coul_tab));

    nbp->coulomb_tab_size     = ic->tabq_size;
    nbp->coulomb_tab_scale    = ic->tabq_scale;
}


/*! Initializes the atomdata structure first time, it only gets filled at
    pair-search. */
static void init_atomdata_first(cu_atomdata_t *ad, int ntypes)
{
    cudaError_t stat;

    ad->ntypes  = ntypes;
    stat        = cudaMalloc((void**)&ad->shift_vec, SHIFTS*sizeof(*ad->shift_vec));
    CU_RET_ERR(stat, "cudaMalloc failed on ad->shift_vec");
    ad->bShiftVecUploaded = false;

    stat = cudaMalloc((void**)&ad->fshift, SHIFTS*sizeof(*ad->fshift));
    CU_RET_ERR(stat, "cudaMalloc failed on ad->fshift");

    stat = cudaMalloc((void**)&ad->e_lj, sizeof(*ad->e_lj));
    CU_RET_ERR(stat, "cudaMalloc failed on ad->e_lj");
    stat = cudaMalloc((void**)&ad->e_el, sizeof(*ad->e_el));
    CU_RET_ERR(stat, "cudaMalloc failed on ad->e_el");

    /* initialize to NULL poiters to data that is not allocated here and will
       need reallocation in nbnxn_cuda_init_atomdata */
    ad->xq = NULL;
    ad->f  = NULL;

    /* size -1 indicates that the respective array hasn't been initialized yet */
    ad->natoms = -1;
    ad->nalloc = -1;
}

/*! Selects the Ewald kernel type, analytical on SM 3.0 and later, tabulated on
    earlier GPUs, single or twin cut-off. */
static int pick_ewald_kernel_type(bool                     bTwinCut,
                                  const gmx_device_info_t *dev_info)
{
    bool bUseAnalyticalEwald, bForceAnalyticalEwald, bForceTabulatedEwald;
    int  kernel_type;

    /* Benchmarking/development environment variables to force the use of
       analytical or tabulated Ewald kernel. */
    bForceAnalyticalEwald = (getenv("GMX_CUDA_NB_ANA_EWALD") != NULL);
    bForceTabulatedEwald  = (getenv("GMX_CUDA_NB_TAB_EWALD") != NULL);

    if (bForceAnalyticalEwald && bForceTabulatedEwald)
    {
        gmx_incons("Both analytical and tabulated Ewald CUDA non-bonded kernels "
                   "requested through environment variables.");
    }

    /* By default, on SM 3.0 and later use analytical Ewald, on earlier tabulated. */
    if ((dev_info->prop.major >= 3 || bForceAnalyticalEwald) && !bForceTabulatedEwald)
    {
        bUseAnalyticalEwald = true;

        if (debug)
        {
            fprintf(debug, "Using analytical Ewald CUDA kernels\n");
        }
    }
    else
    {
        bUseAnalyticalEwald = false;

        if (debug)
        {
            fprintf(debug, "Using tabulated Ewald CUDA kernels\n");
        }
    }

    /* Use twin cut-off kernels if requested by bTwinCut or the env. var.
       forces it (use it for debugging/benchmarking only). */
    if (!bTwinCut && (getenv("GMX_CUDA_NB_EWALD_TWINCUT") == NULL))
    {
        kernel_type = bUseAnalyticalEwald ? eelCuEWALD_ANA : eelCuEWALD_TAB;
    }
    else
    {
        kernel_type = bUseAnalyticalEwald ? eelCuEWALD_ANA_TWIN : eelCuEWALD_TAB_TWIN;
    }

    return kernel_type;
}

/*! Copies all parameters related to the cut-off from ic to nbp */
static void set_cutoff_parameters(cu_nbparam_t              *nbp,
                                  const interaction_const_t *ic)
{
    nbp->ewald_beta       = ic->ewaldcoeff_q;
    nbp->sh_ewald         = ic->sh_ewald;
    nbp->epsfac           = ic->epsfac;
    nbp->two_k_rf         = 2.0 * ic->k_rf;
    nbp->c_rf             = ic->c_rf;
    nbp->rvdw_sq          = ic->rvdw * ic->rvdw;
    nbp->rcoulomb_sq      = ic->rcoulomb * ic->rcoulomb;
    nbp->rlist_sq         = ic->rlist * ic->rlist;

    nbp->sh_lj_ewald      = ic->sh_lj_ewald;
    nbp->ewaldcoeff_lj    = ic->ewaldcoeff_lj;

    nbp->rvdw_switch      = ic->rvdw_switch;
    nbp->dispersion_shift = ic->dispersion_shift;
    nbp->repulsion_shift  = ic->repulsion_shift;
    nbp->vdw_switch       = ic->vdw_switch;
}

/*! Initializes the nonbonded parameter data structure. */
static void init_nbparam(cu_nbparam_t              *nbp,
                         const interaction_const_t *ic,
                         const nbnxn_atomdata_t    *nbat,
                         const gmx_device_info_t   *dev_info)
{
    cudaError_t stat;
    int         ntypes, nnbfp, nnbfp_comb;

    ntypes  = nbat->ntype;

    set_cutoff_parameters(nbp, ic);

    if (ic->vdwtype == evdwCUT)
    {
        switch (ic->vdw_modifier)
        {
            case eintmodNONE:
            case eintmodPOTSHIFT:
                nbp->vdwtype = evdwCuCUT;
                break;
            case eintmodFORCESWITCH:
                nbp->vdwtype = evdwCuFSWITCH;
                break;
            case eintmodPOTSWITCH:
                nbp->vdwtype = evdwCuPSWITCH;
                break;
            default:
                gmx_incons("The requested VdW interaction modifier is not implemented in the CUDA GPU accelerated kernels!");
                break;
        }
    }
    else if (ic->vdwtype == evdwPME)
    {
        if (ic->ljpme_comb_rule == ljcrGEOM)
        {
            assert(nbat->comb_rule == ljcrGEOM);
            nbp->vdwtype = evdwCuEWALDGEOM;
        }
        else
        {
            assert(nbat->comb_rule == ljcrLB);
            nbp->vdwtype = evdwCuEWALDLB;
        }
    }
    else
    {
        gmx_incons("The requested VdW type is not implemented in the CUDA GPU accelerated kernels!");
    }

    if (ic->eeltype == eelCUT)
    {
        nbp->eeltype = eelCuCUT;
    }
    else if (EEL_RF(ic->eeltype))
    {
        nbp->eeltype = eelCuRF;
    }
    else if ((EEL_PME(ic->eeltype) || ic->eeltype == eelEWALD))
    {
        /* Initially rcoulomb == rvdw, so it's surely not twin cut-off. */
        nbp->eeltype = pick_ewald_kernel_type(false, dev_info);
    }
    else
    {
        /* Shouldn't happen, as this is checked when choosing Verlet-scheme */
        gmx_incons("The requested electrostatics type is not implemented in the CUDA GPU accelerated kernels!");
    }

    /* generate table for PME */
    nbp->coulomb_tab = NULL;
    if (nbp->eeltype == eelCuEWALD_TAB || nbp->eeltype == eelCuEWALD_TAB_TWIN)
    {
        init_ewald_coulomb_force_table(ic, nbp, dev_info);
    }

    nnbfp      = 2*ntypes*ntypes;
    nnbfp_comb = 2*ntypes;

    stat  = cudaMalloc((void **)&nbp->nbfp, nnbfp*sizeof(*nbp->nbfp));
    CU_RET_ERR(stat, "cudaMalloc failed on nbp->nbfp");
    cu_copy_H2D(nbp->nbfp, nbat->nbfp, nnbfp*sizeof(*nbp->nbfp));


    if (ic->vdwtype == evdwPME)
    {
        stat  = cudaMalloc((void **)&nbp->nbfp_comb, nnbfp_comb*sizeof(*nbp->nbfp_comb));
        CU_RET_ERR(stat, "cudaMalloc failed on nbp->nbfp_comb");
        cu_copy_H2D(nbp->nbfp_comb, nbat->nbfp_comb, nnbfp_comb*sizeof(*nbp->nbfp_comb));
    }

#ifdef HAVE_CUDA_TEXOBJ_SUPPORT
    /* Only device CC >= 3.0 (Kepler and later) support texture objects */
    if (use_texobj(dev_info))
    {
        cudaResourceDesc rd;
        cudaTextureDesc  td;

        memset(&rd, 0, sizeof(rd));
        rd.resType                  = cudaResourceTypeLinear;
        rd.res.linear.devPtr        = nbp->nbfp;
        rd.res.linear.desc.f        = cudaChannelFormatKindFloat;
        rd.res.linear.desc.x        = 32;
        rd.res.linear.sizeInBytes   = nnbfp*sizeof(*nbp->nbfp);

        memset(&td, 0, sizeof(td));
        td.readMode                 = cudaReadModeElementType;
        stat = cudaCreateTextureObject(&nbp->nbfp_texobj, &rd, &td, NULL);
        CU_RET_ERR(stat, "cudaCreateTextureObject on nbfp_texobj failed");

        if (ic->vdwtype == evdwPME)
        {
            memset(&rd, 0, sizeof(rd));
            rd.resType                  = cudaResourceTypeLinear;
            rd.res.linear.devPtr        = nbp->nbfp_comb;
            rd.res.linear.desc.f        = cudaChannelFormatKindFloat;
            rd.res.linear.desc.x        = 32;
            rd.res.linear.sizeInBytes   = nnbfp_comb*sizeof(*nbp->nbfp_comb);

            memset(&td, 0, sizeof(td));
            td.readMode = cudaReadModeElementType;
            stat        = cudaCreateTextureObject(&nbp->nbfp_comb_texobj, &rd, &td, NULL);
            CU_RET_ERR(stat, "cudaCreateTextureObject on nbfp_comb_texobj failed");
        }
    }
    else
#endif /* HAVE_CUDA_TEXOBJ_SUPPORT */
    {
        cudaChannelFormatDesc cd = cudaCreateChannelDesc<float>();
        stat = cudaBindTexture(NULL, &nbnxn_cuda_get_nbfp_texref(),
                               nbp->nbfp, &cd, nnbfp*sizeof(*nbp->nbfp));
        CU_RET_ERR(stat, "cudaBindTexture on nbfp_texref failed");

        if (ic->vdwtype == evdwPME)
        {
            stat = cudaBindTexture(NULL, &nbnxn_cuda_get_nbfp_comb_texref(),
                                   nbp->nbfp_comb, &cd, nnbfp_comb*sizeof(*nbp->nbfp_comb));
            CU_RET_ERR(stat, "cudaBindTexture on nbfp_comb_texref failed");
        }
    }
}

/*! Re-generate the GPU Ewald force table, resets rlist, and update the
 *  electrostatic type switching to twin cut-off (or back) if needed. */
void nbnxn_gpu_pme_loadbal_update_param(const nonbonded_verlet_t    *nbv,
                                        const interaction_const_t   *ic)
{
    if (!nbv || nbv->grp[0].kernel_type != nbnxnk8x8x8_GPU)
    {
        return;
    }
    gmx_nbnxn_cuda_t *nb    = nbv->gpu_nbv;
    cu_nbparam_t     *nbp   = nb->nbparam;

    set_cutoff_parameters(nbp, ic);

    nbp->eeltype        = pick_ewald_kernel_type(ic->rcoulomb != ic->rvdw,
                                                 nb->dev_info);

    init_ewald_coulomb_force_table(ic, nb->nbparam, nb->dev_info);
}

/*! Initializes the pair list data structure. */
static void init_plist(cu_plist_t *pl)
{
    /* initialize to NULL pointers to data that is not allocated here and will
       need reallocation in nbnxn_gpu_init_pairlist */
    pl->sci     = NULL;
    pl->cj4     = NULL;
    pl->excl    = NULL;

    /* size -1 indicates that the respective array hasn't been initialized yet */
    pl->na_c        = -1;
    pl->nsci        = -1;
    pl->sci_nalloc  = -1;
    pl->ncj4        = -1;
    pl->cj4_nalloc  = -1;
    pl->nexcl       = -1;
    pl->excl_nalloc = -1;
    pl->bDoPrune    = false;
}

/*! Initializes the timer data structure. */
static void init_timers(cu_timers_t *t, bool bUseTwoStreams)
{
    cudaError_t stat;
    int         eventflags = ( bUseCudaEventBlockingSync ? cudaEventBlockingSync : cudaEventDefault );

    stat = cudaEventCreateWithFlags(&(t->start_atdat), eventflags);
    CU_RET_ERR(stat, "cudaEventCreate on start_atdat failed");
    stat = cudaEventCreateWithFlags(&(t->stop_atdat), eventflags);
    CU_RET_ERR(stat, "cudaEventCreate on stop_atdat failed");

    /* The non-local counters/stream (second in the array) are needed only with DD. */
    for (int i = 0; i <= (bUseTwoStreams ? 1 : 0); i++)
    {
        stat = cudaEventCreateWithFlags(&(t->start_nb_k[i]), eventflags);
        CU_RET_ERR(stat, "cudaEventCreate on start_nb_k failed");
        stat = cudaEventCreateWithFlags(&(t->stop_nb_k[i]), eventflags);
        CU_RET_ERR(stat, "cudaEventCreate on stop_nb_k failed");


        stat = cudaEventCreateWithFlags(&(t->start_pl_h2d[i]), eventflags);
        CU_RET_ERR(stat, "cudaEventCreate on start_pl_h2d failed");
        stat = cudaEventCreateWithFlags(&(t->stop_pl_h2d[i]), eventflags);
        CU_RET_ERR(stat, "cudaEventCreate on stop_pl_h2d failed");

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

/*! Initializes the timings data structure. */
static void init_timings(gmx_wallclock_gpu_t *t)
{
    int i, j;

    t->nb_h2d_t = 0.0;
    t->nb_d2h_t = 0.0;
    t->nb_c     = 0;
    t->pl_h2d_t = 0.0;
    t->pl_h2d_c = 0;
    for (i = 0; i < 2; i++)
    {
        for (j = 0; j < 2; j++)
        {
            t->ktime[i][j].t = 0.0;
            t->ktime[i][j].c = 0;
        }
    }
}

/*! Initializes simulation constant data. */
static void nbnxn_cuda_init_const(gmx_nbnxn_cuda_t               *nb,
                                  const interaction_const_t      *ic,
                                  const nonbonded_verlet_group_t *nbv_group)
{
    init_atomdata_first(nb->atdat, nbv_group[0].nbat->ntype);
    init_nbparam(nb->nbparam, ic, nbv_group[0].nbat, nb->dev_info);

    /* clear energy and shift force outputs */
    nbnxn_cuda_clear_e_fshift(nb);
}

void nbnxn_gpu_init(FILE                      *fplog,
                    gmx_nbnxn_cuda_t         **p_nb,
                    const gmx_gpu_info_t      *gpu_info,
                    const gmx_gpu_opt_t       *gpu_opt,
                    const interaction_const_t *ic,
                    nonbonded_verlet_group_t  *nbv_grp,
                    int                        my_gpu_index,
                    int                        /*rank*/,
                    gmx_bool                   bLocalAndNonlocal)
{
    cudaError_t       stat;
    gmx_nbnxn_cuda_t *nb;
    char              sbuf[STRLEN];
    bool              bStreamSync, bNoStreamSync, bTMPIAtomics, bX86, bOldDriver;
    int               cuda_drv_ver;

    assert(gpu_info);

    if (p_nb == NULL)
    {
        return;
    }

    snew(nb, 1);
    snew(nb->atdat, 1);
    snew(nb->nbparam, 1);
    snew(nb->plist[eintLocal], 1);
    if (bLocalAndNonlocal)
    {
        snew(nb->plist[eintNonlocal], 1);
    }

    nb->bUseTwoStreams = bLocalAndNonlocal;

    snew(nb->timers, 1);
    snew(nb->timings, 1);

    /* init nbst */
    pmalloc((void**)&nb->nbst.e_lj, sizeof(*nb->nbst.e_lj));
    pmalloc((void**)&nb->nbst.e_el, sizeof(*nb->nbst.e_el));
    pmalloc((void**)&nb->nbst.fshift, SHIFTS * sizeof(*nb->nbst.fshift));

    init_plist(nb->plist[eintLocal]);

    /* set device info, just point it to the right GPU among the detected ones */
    nb->dev_info = &gpu_info->gpu_dev[get_gpu_device_id(gpu_info, gpu_opt, my_gpu_index)];

    /* local/non-local GPU streams */
    stat = cudaStreamCreate(&nb->stream[eintLocal]);
    CU_RET_ERR(stat, "cudaStreamCreate on stream[eintLocal] failed");
    if (nb->bUseTwoStreams)
    {
        init_plist(nb->plist[eintNonlocal]);

        /* CUDA stream priority available in the CUDA RT 5.5 API.
         * Note that the device we're running on does not have to support
         * priorities, because we are querying the priority range which in this
         * case will be a single value.
         */
#if GMX_CUDA_VERSION >= 5050
        {
            int highest_priority;
            stat = cudaDeviceGetStreamPriorityRange(NULL, &highest_priority);
            CU_RET_ERR(stat, "cudaDeviceGetStreamPriorityRange failed");

            stat = cudaStreamCreateWithPriority(&nb->stream[eintNonlocal],
                                                cudaStreamDefault,
                                                highest_priority);
            CU_RET_ERR(stat, "cudaStreamCreateWithPriority on stream[eintNonlocal] failed");
        }
#else
        stat = cudaStreamCreate(&nb->stream[eintNonlocal]);
        CU_RET_ERR(stat, "cudaStreamCreate on stream[eintNonlocal] failed");
#endif
    }

    /* init events for sychronization (timing disabled for performance reasons!) */
    stat = cudaEventCreateWithFlags(&nb->nonlocal_done, cudaEventDisableTiming);
    CU_RET_ERR(stat, "cudaEventCreate on nonlocal_done failed");
    stat = cudaEventCreateWithFlags(&nb->misc_ops_and_local_H2D_done, cudaEventDisableTiming);
    CU_RET_ERR(stat, "cudaEventCreate on misc_ops_and_local_H2D_done failed");

    /* On GPUs with ECC enabled, cudaStreamSynchronize shows a large overhead
     * (which increases with shorter time/step) caused by a known CUDA driver bug.
     * To work around the issue we'll use an (admittedly fragile) memory polling
     * waiting to preserve performance. This requires support for atomic
     * operations and only works on x86/x86_64.
     * With polling wait event-timing also needs to be disabled.
     *
     * The overhead is greatly reduced in API v5.0 drivers and the improvement
     * is independent of runtime version. Hence, with API v5.0 drivers and later
     * we won't switch to polling.
     *
     * NOTE: Unfortunately, this is known to fail when GPUs are shared by (t)MPI,
     * ranks so we will also disable it in that case.
     */

    bStreamSync    = getenv("GMX_CUDA_STREAMSYNC") != NULL;
    bNoStreamSync  = getenv("GMX_NO_CUDA_STREAMSYNC") != NULL;

#ifdef TMPI_ATOMICS
    bTMPIAtomics = true;
#else
    bTMPIAtomics = false;
#endif

#ifdef GMX_TARGET_X86
    bX86 = true;
#else
    bX86 = false;
#endif

    if (bStreamSync && bNoStreamSync)
    {
        gmx_fatal(FARGS, "Conflicting environment variables: both GMX_CUDA_STREAMSYNC and GMX_NO_CUDA_STREAMSYNC defined");
    }

    stat = cudaDriverGetVersion(&cuda_drv_ver);
    CU_RET_ERR(stat, "cudaDriverGetVersion failed");

    bOldDriver = (cuda_drv_ver < 5000);

    if ((nb->dev_info->prop.ECCEnabled == 1) && bOldDriver)
    {
        /* Polling wait should be used instead of cudaStreamSynchronize only if:
         *   - ECC is ON & driver is old (checked above),
         *   - we're on x86/x86_64,
         *   - atomics are available, and
         *   - GPUs are not being shared.
         */
        bool bShouldUsePollSync = (bX86 && bTMPIAtomics &&
                                   (gmx_count_gpu_dev_shared(gpu_opt) < 1));

        if (bStreamSync)
        {
            nb->bUseStreamSync = true;

            /* only warn if polling should be used */
            if (bShouldUsePollSync)
            {
                md_print_warn(fplog,
                              "NOTE: Using a GPU with ECC enabled and CUDA driver API version <5.0, but\n"
                              "      cudaStreamSynchronize waiting is forced by the GMX_CUDA_STREAMSYNC env. var.\n");
            }
        }
        else
        {
            nb->bUseStreamSync = !bShouldUsePollSync;

            if (bShouldUsePollSync)
            {
                md_print_warn(fplog,
                              "NOTE: Using a GPU with ECC enabled and CUDA driver API version <5.0, known to\n"
                              "      cause performance loss. Switching to the alternative polling GPU wait.\n"
                              "      If you encounter issues, switch back to standard GPU waiting by setting\n"
                              "      the GMX_CUDA_STREAMSYNC environment variable.\n");
            }
            else
            {
                /* Tell the user that the ECC+old driver combination can be bad */
                sprintf(sbuf,
                        "NOTE: Using a GPU with ECC enabled and CUDA driver API version <5.0.\n"
                        "      A known bug in this driver version can cause performance loss.\n"
                        "      However, the polling wait workaround can not be used because\n%s\n"
                        "      Consider updating the driver or turning ECC off.",
                        (bX86 && bTMPIAtomics) ?
                        "      GPU(s) are being oversubscribed." :
                        "      atomic operations are not supported by the platform/CPU+compiler.");
                md_print_warn(fplog, sbuf);
            }
        }
    }
    else
    {
        if (bNoStreamSync)
        {
            nb->bUseStreamSync = false;

            md_print_warn(fplog,
                          "NOTE: Polling wait for GPU synchronization requested by GMX_NO_CUDA_STREAMSYNC\n");
        }
        else
        {
            /* no/off ECC, cudaStreamSynchronize not turned off by env. var. */
            nb->bUseStreamSync = true;
        }
    }

    /* CUDA timing disabled as event timers don't work:
       - with multiple streams = domain-decomposition;
       - with the polling waiting hack (without cudaStreamSynchronize);
       - when turned off by GMX_DISABLE_CUDA_TIMING.
     */
    nb->bDoTime = (!nb->bUseTwoStreams && nb->bUseStreamSync &&
                   (getenv("GMX_DISABLE_CUDA_TIMING") == NULL));

    if (nb->bDoTime)
    {
        init_timers(nb->timers, nb->bUseTwoStreams);
        init_timings(nb->timings);
    }

    /* set the kernel type for the current GPU */
    /* pick L1 cache configuration */
    nbnxn_cuda_set_cacheconfig(nb->dev_info);

    nbnxn_cuda_init_const(nb, ic, nbv_grp);

    *p_nb = nb;

    if (debug)
    {
        fprintf(debug, "Initialized CUDA data structures.\n");
    }
}

void nbnxn_gpu_init_pairlist(gmx_nbnxn_cuda_t       *nb,
                             const nbnxn_pairlist_t *h_plist,
                             int                     iloc)
{
    char          sbuf[STRLEN];
    cudaError_t   stat;
    bool          bDoTime    = nb->bDoTime;
    cudaStream_t  stream     = nb->stream[iloc];
    cu_plist_t   *d_plist    = nb->plist[iloc];

    if (d_plist->na_c < 0)
    {
        d_plist->na_c = h_plist->na_ci;
    }
    else
    {
        if (d_plist->na_c != h_plist->na_ci)
        {
            sprintf(sbuf, "In cu_init_plist: the #atoms per cell has changed (from %d to %d)",
                    d_plist->na_c, h_plist->na_ci);
            gmx_incons(sbuf);
        }
    }

    if (bDoTime)
    {
        stat = cudaEventRecord(nb->timers->start_pl_h2d[iloc], stream);
        CU_RET_ERR(stat, "cudaEventRecord failed");
    }

    cu_realloc_buffered((void **)&d_plist->sci, h_plist->sci, sizeof(*d_plist->sci),
                        &d_plist->nsci, &d_plist->sci_nalloc,
                        h_plist->nsci,
                        stream, true);

    cu_realloc_buffered((void **)&d_plist->cj4, h_plist->cj4, sizeof(*d_plist->cj4),
                        &d_plist->ncj4, &d_plist->cj4_nalloc,
                        h_plist->ncj4,
                        stream, true);

    cu_realloc_buffered((void **)&d_plist->excl, h_plist->excl, sizeof(*d_plist->excl),
                        &d_plist->nexcl, &d_plist->excl_nalloc,
                        h_plist->nexcl,
                        stream, true);

    if (bDoTime)
    {
        stat = cudaEventRecord(nb->timers->stop_pl_h2d[iloc], stream);
        CU_RET_ERR(stat, "cudaEventRecord failed");
    }

    /* need to prune the pair list during the next step */
    d_plist->bDoPrune = true;
}

void nbnxn_gpu_upload_shiftvec(gmx_nbnxn_cuda_t       *nb,
                               const nbnxn_atomdata_t *nbatom)
{
    cu_atomdata_t *adat  = nb->atdat;
    cudaStream_t   ls    = nb->stream[eintLocal];

    /* only if we have a dynamic box */
    if (nbatom->bDynamicBox || !adat->bShiftVecUploaded)
    {
        cu_copy_H2D_async(adat->shift_vec, nbatom->shift_vec,
                          SHIFTS * sizeof(*adat->shift_vec), ls);
        adat->bShiftVecUploaded = true;
    }
}

/*! Clears the first natoms_clear elements of the GPU nonbonded force output array. */
static void nbnxn_cuda_clear_f(gmx_nbnxn_cuda_t *nb, int natoms_clear)
{
    cudaError_t    stat;
    cu_atomdata_t *adat  = nb->atdat;
    cudaStream_t   ls    = nb->stream[eintLocal];

    stat = cudaMemsetAsync(adat->f, 0, natoms_clear * sizeof(*adat->f), ls);
    CU_RET_ERR(stat, "cudaMemsetAsync on f falied");
}

/*! Clears nonbonded shift force output array and energy outputs on the GPU. */
static void nbnxn_cuda_clear_e_fshift(gmx_nbnxn_cuda_t *nb)
{
    cudaError_t    stat;
    cu_atomdata_t *adat  = nb->atdat;
    cudaStream_t   ls    = nb->stream[eintLocal];

    stat = cudaMemsetAsync(adat->fshift, 0, SHIFTS * sizeof(*adat->fshift), ls);
    CU_RET_ERR(stat, "cudaMemsetAsync on fshift falied");
    stat = cudaMemsetAsync(adat->e_lj, 0, sizeof(*adat->e_lj), ls);
    CU_RET_ERR(stat, "cudaMemsetAsync on e_lj falied");
    stat = cudaMemsetAsync(adat->e_el, 0, sizeof(*adat->e_el), ls);
    CU_RET_ERR(stat, "cudaMemsetAsync on e_el falied");
}

void nbnxn_gpu_clear_outputs(gmx_nbnxn_cuda_t *nb, int flags)
{
    nbnxn_cuda_clear_f(nb, nb->atdat->natoms);
    /* clear shift force array and energies if the outputs were
       used in the current step */
    if (flags & GMX_FORCE_VIRIAL)
    {
        nbnxn_cuda_clear_e_fshift(nb);
    }
}

void nbnxn_gpu_init_atomdata(gmx_nbnxn_cuda_t              *nb,
                             const struct nbnxn_atomdata_t *nbat)
{
    cudaError_t    stat;
    int            nalloc, natoms;
    bool           realloced;
    bool           bDoTime   = nb->bDoTime;
    cu_timers_t   *timers    = nb->timers;
    cu_atomdata_t *d_atdat   = nb->atdat;
    cudaStream_t   ls        = nb->stream[eintLocal];

    natoms    = nbat->natoms;
    realloced = false;

    if (bDoTime)
    {
        /* time async copy */
        stat = cudaEventRecord(timers->start_atdat, ls);
        CU_RET_ERR(stat, "cudaEventRecord failed");
    }

    /* need to reallocate if we have to copy more atoms than the amount of space
       available and only allocate if we haven't initialized yet, i.e d_atdat->natoms == -1 */
    if (natoms > d_atdat->nalloc)
    {
        nalloc = over_alloc_small(natoms);

        /* free up first if the arrays have already been initialized */
        if (d_atdat->nalloc != -1)
        {
            cu_free_buffered(d_atdat->f, &d_atdat->natoms, &d_atdat->nalloc);
            cu_free_buffered(d_atdat->xq);
            cu_free_buffered(d_atdat->atom_types);
        }

        stat = cudaMalloc((void **)&d_atdat->f, nalloc*sizeof(*d_atdat->f));
        CU_RET_ERR(stat, "cudaMalloc failed on d_atdat->f");
        stat = cudaMalloc((void **)&d_atdat->xq, nalloc*sizeof(*d_atdat->xq));
        CU_RET_ERR(stat, "cudaMalloc failed on d_atdat->xq");

        stat = cudaMalloc((void **)&d_atdat->atom_types, nalloc*sizeof(*d_atdat->atom_types));
        CU_RET_ERR(stat, "cudaMalloc failed on d_atdat->atom_types");

        d_atdat->nalloc = nalloc;
        realloced       = true;
    }

    d_atdat->natoms       = natoms;
    d_atdat->natoms_local = nbat->natoms_local;

    /* need to clear GPU f output if realloc happened */
    if (realloced)
    {
        nbnxn_cuda_clear_f(nb, nalloc);
    }

    cu_copy_H2D_async(d_atdat->atom_types, nbat->type,
                      natoms*sizeof(*d_atdat->atom_types), ls);

    if (bDoTime)
    {
        stat = cudaEventRecord(timers->stop_atdat, ls);
        CU_RET_ERR(stat, "cudaEventRecord failed");
    }
}

static void nbnxn_cuda_free_nbparam_table(cu_nbparam_t            *nbparam,
                                          const gmx_device_info_t *dev_info)
{
    cudaError_t stat;

    if (nbparam->eeltype == eelCuEWALD_TAB || nbparam->eeltype == eelCuEWALD_TAB_TWIN)
    {
#ifdef HAVE_CUDA_TEXOBJ_SUPPORT
        /* Only device CC >= 3.0 (Kepler and later) support texture objects */
        if (use_texobj(dev_info))
        {
            stat = cudaDestroyTextureObject(nbparam->coulomb_tab_texobj);
            CU_RET_ERR(stat, "cudaDestroyTextureObject on coulomb_tab_texobj failed");
        }
        else
#endif
        {
            GMX_UNUSED_VALUE(dev_info);
            stat = cudaUnbindTexture(nbnxn_cuda_get_coulomb_tab_texref());
            CU_RET_ERR(stat, "cudaUnbindTexture on coulomb_tab_texref failed");
        }
        cu_free_buffered(nbparam->coulomb_tab, &nbparam->coulomb_tab_size);
    }
}

void nbnxn_gpu_free(gmx_nbnxn_cuda_t *nb)
{
    cudaError_t      stat;
    cu_atomdata_t   *atdat;
    cu_nbparam_t    *nbparam;
    cu_plist_t      *plist, *plist_nl;
    cu_timers_t     *timers;

    /* Stopping the nvidia profiler here allows us to eliminate the subsequent
       uninitialization API calls from the trace. */
    if (getenv("NVPROF_ID") != NULL)
    {
        stat = cudaProfilerStop();
        CU_RET_ERR(stat, "cudaProfilerStop failed");
    }

    if (nb == NULL)
    {
        return;
    }

    atdat       = nb->atdat;
    nbparam     = nb->nbparam;
    plist       = nb->plist[eintLocal];
    plist_nl    = nb->plist[eintNonlocal];
    timers      = nb->timers;

    nbnxn_cuda_free_nbparam_table(nbparam, nb->dev_info);

    stat = cudaEventDestroy(nb->nonlocal_done);
    CU_RET_ERR(stat, "cudaEventDestroy failed on timers->nonlocal_done");
    stat = cudaEventDestroy(nb->misc_ops_and_local_H2D_done);
    CU_RET_ERR(stat, "cudaEventDestroy failed on timers->misc_ops_and_local_H2D_done");

    if (nb->bDoTime)
    {
        stat = cudaEventDestroy(timers->start_atdat);
        CU_RET_ERR(stat, "cudaEventDestroy failed on timers->start_atdat");
        stat = cudaEventDestroy(timers->stop_atdat);
        CU_RET_ERR(stat, "cudaEventDestroy failed on timers->stop_atdat");

        /* The non-local counters/stream (second in the array) are needed only with DD. */
        for (int i = 0; i <= (nb->bUseTwoStreams ? 1 : 0); i++)
        {
            stat = cudaEventDestroy(timers->start_nb_k[i]);
            CU_RET_ERR(stat, "cudaEventDestroy failed on timers->start_nb_k");
            stat = cudaEventDestroy(timers->stop_nb_k[i]);
            CU_RET_ERR(stat, "cudaEventDestroy failed on timers->stop_nb_k");

            stat = cudaEventDestroy(timers->start_pl_h2d[i]);
            CU_RET_ERR(stat, "cudaEventDestroy failed on timers->start_pl_h2d");
            stat = cudaEventDestroy(timers->stop_pl_h2d[i]);
            CU_RET_ERR(stat, "cudaEventDestroy failed on timers->stop_pl_h2d");

            stat = cudaStreamDestroy(nb->stream[i]);
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
    }

#ifdef HAVE_CUDA_TEXOBJ_SUPPORT
    /* Only device CC >= 3.0 (Kepler and later) support texture objects */
    if (use_texobj(nb->dev_info))
    {
        stat = cudaDestroyTextureObject(nbparam->nbfp_texobj);
        CU_RET_ERR(stat, "cudaDestroyTextureObject on nbfp_texobj failed");
    }
    else
#endif
    {
        stat = cudaUnbindTexture(nbnxn_cuda_get_nbfp_texref());
        CU_RET_ERR(stat, "cudaUnbindTexture on nbfp_texref failed");
    }
    cu_free_buffered(nbparam->nbfp);

    if (nbparam->vdwtype == evdwCuEWALDGEOM || nbparam->vdwtype == evdwCuEWALDLB)
    {
#ifdef HAVE_CUDA_TEXOBJ_SUPPORT
        /* Only device CC >= 3.0 (Kepler and later) support texture objects */
        if (use_texobj(nb->dev_info))
        {
            stat = cudaDestroyTextureObject(nbparam->nbfp_comb_texobj);
            CU_RET_ERR(stat, "cudaDestroyTextureObject on nbfp_comb_texobj failed");
        }
        else
#endif
        {
            stat = cudaUnbindTexture(nbnxn_cuda_get_nbfp_comb_texref());
            CU_RET_ERR(stat, "cudaUnbindTexture on nbfp_comb_texref failed");
        }
        cu_free_buffered(nbparam->nbfp_comb);
    }

    stat = cudaFree(atdat->shift_vec);
    CU_RET_ERR(stat, "cudaFree failed on atdat->shift_vec");
    stat = cudaFree(atdat->fshift);
    CU_RET_ERR(stat, "cudaFree failed on atdat->fshift");

    stat = cudaFree(atdat->e_lj);
    CU_RET_ERR(stat, "cudaFree failed on atdat->e_lj");
    stat = cudaFree(atdat->e_el);
    CU_RET_ERR(stat, "cudaFree failed on atdat->e_el");

    cu_free_buffered(atdat->f, &atdat->natoms, &atdat->nalloc);
    cu_free_buffered(atdat->xq);
    cu_free_buffered(atdat->atom_types, &atdat->ntypes);

    cu_free_buffered(plist->sci, &plist->nsci, &plist->sci_nalloc);
    cu_free_buffered(plist->cj4, &plist->ncj4, &plist->cj4_nalloc);
    cu_free_buffered(plist->excl, &plist->nexcl, &plist->excl_nalloc);
    if (nb->bUseTwoStreams)
    {
        cu_free_buffered(plist_nl->sci, &plist_nl->nsci, &plist_nl->sci_nalloc);
        cu_free_buffered(plist_nl->cj4, &plist_nl->ncj4, &plist_nl->cj4_nalloc);
        cu_free_buffered(plist_nl->excl, &plist_nl->nexcl, &plist->excl_nalloc);
    }

    sfree(atdat);
    sfree(nbparam);
    sfree(plist);
    if (nb->bUseTwoStreams)
    {
        sfree(plist_nl);
    }
    sfree(timers);
    sfree(nb->timings);
    sfree(nb);

    if (debug)
    {
        fprintf(debug, "Cleaned up CUDA data structures.\n");
    }
}

void cu_synchstream_atdat(gmx_nbnxn_cuda_t *nb, int iloc)
{
    cudaError_t  stat;
    cudaStream_t stream = nb->stream[iloc];

    stat = cudaStreamWaitEvent(stream, nb->timers->stop_atdat, 0);
    CU_RET_ERR(stat, "cudaStreamWaitEvent failed");
}

gmx_wallclock_gpu_t * nbnxn_gpu_get_timings(gmx_nbnxn_cuda_t *nb)
{
    return (nb != NULL && nb->bDoTime) ? nb->timings : NULL;
}

void nbnxn_gpu_reset_timings(nonbonded_verlet_t* nbv)
{
    /* The NVPROF_ID environment variable is set by nvprof and indicates that
       mdrun is executed in the CUDA profiler.
       If nvprof was run is with "--profile-from-start off", the profiler will
       be started here. This way we can avoid tracing the CUDA events from the
       first part of the run. Starting the profiler again does nothing.
     */
    if (getenv("NVPROF_ID") != NULL)
    {
        cudaError_t stat;
        stat = cudaProfilerStart();
        CU_RET_ERR(stat, "cudaProfilerStart failed");
    }

    if (nbv->gpu_nbv && nbv->gpu_nbv->bDoTime)
    {
        init_timings(nbv->gpu_nbv->timings);
    }
}

int nbnxn_gpu_min_ci_balanced(gmx_nbnxn_cuda_t *nb)
{
    return nb != NULL ?
           gpu_min_ci_balanced_factor*nb->dev_info->prop.multiProcessorCount : 0;

}

gmx_bool nbnxn_gpu_is_kernel_ewald_analytical(const gmx_nbnxn_cuda_t *nb)
{
    return ((nb->nbparam->eeltype == eelCuEWALD_ANA) ||
            (nb->nbparam->eeltype == eelCuEWALD_ANA_TWIN));
}
