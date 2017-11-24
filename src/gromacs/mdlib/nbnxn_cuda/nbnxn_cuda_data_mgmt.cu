/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2012,2013,2014,2015,2016,2017, by the GROMACS development team, led by
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

#include <assert.h>
#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>

#include "gromacs/gpu_utils/cudautils.cuh"
#include "gromacs/gpu_utils/gpu_utils.h"
#include "gromacs/gpu_utils/pmalloc_cuda.h"
#include "gromacs/hardware/gpu_hw_info.h"
#include "gromacs/math/vectypes.h"
#include "gromacs/mdlib/force_flags.h"
#include "gromacs/mdlib/nb_verlet.h"
#include "gromacs/mdlib/nbnxn_consts.h"
#include "gromacs/mdlib/nbnxn_gpu_common_utils.h"
#include "gromacs/mdlib/nbnxn_gpu_data_mgmt.h"
#include "gromacs/mdtypes/interaction_const.h"
#include "gromacs/mdtypes/md_enums.h"
#include "gromacs/pbcutil/ishift.h"
#include "gromacs/timing/gpu_timing.h"
#include "gromacs/utility/basedefinitions.h"
#include "gromacs/utility/cstringutil.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/real.h"
#include "gromacs/utility/smalloc.h"

#include "nbnxn_cuda.h"
#include "nbnxn_cuda_types.h"

/* This is a heuristically determined parameter for the Fermi, Kepler
 * and Maxwell architectures for the minimum size of ci lists by multiplying
 * this constant with the # of multiprocessors on the current device.
 * Since the maximum number of blocks per multiprocessor is 16, the ideal
 * count for small systems is 32 or 48 blocks per multiprocessor. Because
 * there is a bit of fluctuations in the generated block counts, we use
 * a target of 44 instead of the ideal value of 48.
 */
static unsigned int gpu_min_ci_balanced_factor = 44;

/* Fw. decl. */
static void nbnxn_cuda_clear_e_fshift(gmx_nbnxn_cuda_t *nb);

/* Fw. decl, */
static void nbnxn_cuda_free_nbparam_table(cu_nbparam_t            *nbparam,
                                          const gmx_device_info_t *dev_info);

/*! \brief Return whether combination rules are used.
 *
 * \param[in]   pointer to nonbonded paramter struct
 * \return      true if combination rules are used in this run, false otherwise
 */
static inline bool useLjCombRule(const cu_nbparam_t  *nbparam)
{
    return (nbparam->vdwtype == evdwCuCUTCOMBGEOM ||
            nbparam->vdwtype == evdwCuCUTCOMBLB);
}

/*! \brief Initialized the Ewald Coulomb correction GPU table.

    Tabulates the Ewald Coulomb force and initializes the size/scale
    and the table GPU array. If called with an already allocated table,
    it just re-uploads the table.
 */
static void init_ewald_coulomb_force_table(const interaction_const_t *ic,
                                           cu_nbparam_t              *nbp,
                                           const gmx_device_info_t   *dev_info)
{
    if (nbp->coulomb_tab != NULL)
    {
        nbnxn_cuda_free_nbparam_table(nbp, dev_info);
    }

    nbp->coulomb_tab_scale = ic->tabq_scale;
    initParamLookupTable(nbp->coulomb_tab, nbp->coulomb_tab_texobj,
                         &nbnxn_cuda_get_coulomb_tab_texref(),
                         ic->tabq_coul_F, ic->tabq_size, dev_info);
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
                                  const interaction_const_t *ic,
                                  const NbnxnListParameters *listParams)
{
    nbp->ewald_beta        = ic->ewaldcoeff_q;
    nbp->sh_ewald          = ic->sh_ewald;
    nbp->epsfac            = ic->epsfac;
    nbp->two_k_rf          = 2.0 * ic->k_rf;
    nbp->c_rf              = ic->c_rf;
    nbp->rvdw_sq           = ic->rvdw * ic->rvdw;
    nbp->rcoulomb_sq       = ic->rcoulomb * ic->rcoulomb;
    nbp->rlistOuter_sq     = listParams->rlistOuter * listParams->rlistOuter;
    nbp->rlistInner_sq     = listParams->rlistInner * listParams->rlistInner;
    nbp->useDynamicPruning = listParams->useDynamicPruning;

    nbp->sh_lj_ewald       = ic->sh_lj_ewald;
    nbp->ewaldcoeff_lj     = ic->ewaldcoeff_lj;

    nbp->rvdw_switch       = ic->rvdw_switch;
    nbp->dispersion_shift  = ic->dispersion_shift;
    nbp->repulsion_shift   = ic->repulsion_shift;
    nbp->vdw_switch        = ic->vdw_switch;
}

/*! Initializes the nonbonded parameter data structure. */
static void init_nbparam(cu_nbparam_t              *nbp,
                         const interaction_const_t *ic,
                         const NbnxnListParameters *listParams,
                         const nbnxn_atomdata_t    *nbat,
                         const gmx_device_info_t   *dev_info)
{
    int         ntypes;

    ntypes  = nbat->ntype;

    set_cutoff_parameters(nbp, ic, listParams);

    /* The kernel code supports LJ combination rules (geometric and LB) for
     * all kernel types, but we only generate useful combination rule kernels.
     * We currently only use LJ combination rule (geometric and LB) kernels
     * for plain cut-off LJ. On Maxwell the force only kernels speed up 15%
     * with PME and 20% with RF, the other kernels speed up about half as much.
     * For LJ force-switch the geometric rule would give 7% speed-up, but this
     * combination is rarely used. LJ force-switch with LB rule is more common,
     * but gives only 1% speed-up.
     */
    if (ic->vdwtype == evdwCUT)
    {
        switch (ic->vdw_modifier)
        {
            case eintmodNONE:
            case eintmodPOTSHIFT:
                switch (nbat->comb_rule)
                {
                    case ljcrNONE:
                        nbp->vdwtype = evdwCuCUT;
                        break;
                    case ljcrGEOM:
                        nbp->vdwtype = evdwCuCUTCOMBGEOM;
                        break;
                    case ljcrLB:
                        nbp->vdwtype = evdwCuCUTCOMBLB;
                        break;
                    default:
                        gmx_incons("The requested LJ combination rule is not implemented in the CUDA GPU accelerated kernels!");
                        break;
                }
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

    /* set up LJ parameter lookup table */
    if (!useLjCombRule(nbp))
    {
        initParamLookupTable(nbp->nbfp, nbp->nbfp_texobj,
                             &nbnxn_cuda_get_nbfp_texref(),
                             nbat->nbfp, 2*ntypes*ntypes, dev_info);
    }

    /* set up LJ-PME parameter lookup table */
    if (ic->vdwtype == evdwPME)
    {
        initParamLookupTable(nbp->nbfp_comb, nbp->nbfp_comb_texobj,
                             &nbnxn_cuda_get_nbfp_comb_texref(),
                             nbat->nbfp_comb, 2*ntypes, dev_info);
    }
}

/*! Re-generate the GPU Ewald force table, resets rlist, and update the
 *  electrostatic type switching to twin cut-off (or back) if needed. */
void nbnxn_gpu_pme_loadbal_update_param(const nonbonded_verlet_t    *nbv,
                                        const interaction_const_t   *ic,
                                        const NbnxnListParameters   *listParams)
{
    if (!nbv || nbv->grp[0].kernel_type != nbnxnk8x8x8_GPU)
    {
        return;
    }
    gmx_nbnxn_cuda_t *nb    = nbv->gpu_nbv;
    cu_nbparam_t     *nbp   = nb->nbparam;

    set_cutoff_parameters(nbp, ic, listParams);

    nbp->eeltype        = pick_ewald_kernel_type(ic->rcoulomb != ic->rvdw,
                                                 nb->dev_info);

    init_ewald_coulomb_force_table(ic, nb->nbparam, nb->dev_info);
}

/*! Initializes the pair list data structure. */
static void init_plist(cu_plist_t *pl)
{
    /* initialize to NULL pointers to data that is not allocated here and will
       need reallocation in nbnxn_gpu_init_pairlist */
    pl->sci      = NULL;
    pl->cj4      = NULL;
    pl->imask    = NULL;
    pl->excl     = NULL;

    /* size -1 indicates that the respective array hasn't been initialized yet */
    pl->na_c           = -1;
    pl->nsci           = -1;
    pl->sci_nalloc     = -1;
    pl->ncj4           = -1;
    pl->cj4_nalloc     = -1;
    pl->nimask         = -1;
    pl->imask_nalloc   = -1;
    pl->nexcl          = -1;
    pl->excl_nalloc    = -1;
    pl->haveFreshList  = false;
}

/*! Initializes the timer data structure. */
static void init_timers(cu_timers_t *t, bool bUseTwoStreams)
{
    /* The non-local counters/stream (second in the array) are needed only with DD. */
    for (int i = 0; i <= (bUseTwoStreams ? 1 : 0); i++)
    {
        t->didPairlistH2D[i]  = false;
        t->didPrune[i]        = false;
        t->didRollingPrune[i] = false;
    }
}

/*! Initializes the timings data structure. */
static void init_timings(gmx_wallclock_gpu_nbnxn_t *t)
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
    t->pruneTime.c        = 0;
    t->pruneTime.t        = 0.0;
    t->dynamicPruneTime.c = 0;
    t->dynamicPruneTime.t = 0.0;
}

/*! Initializes simulation constant data. */
static void nbnxn_cuda_init_const(gmx_nbnxn_cuda_t               *nb,
                                  const interaction_const_t      *ic,
                                  const NbnxnListParameters      *listParams,
                                  const nbnxn_atomdata_t         *nbat)
{
    init_atomdata_first(nb->atdat, nbat->ntype);
    init_nbparam(nb->nbparam, ic, listParams, nbat, nb->dev_info);

    /* clear energy and shift force outputs */
    nbnxn_cuda_clear_e_fshift(nb);
}

void nbnxn_gpu_init(gmx_nbnxn_cuda_t         **p_nb,
                    const gmx_device_info_t   *deviceInfo,
                    const interaction_const_t *ic,
                    const NbnxnListParameters *listParams,
                    const nbnxn_atomdata_t    *nbat,
                    int                        /*rank*/,
                    gmx_bool                   bLocalAndNonlocal)
{
    cudaError_t       stat;
    gmx_nbnxn_cuda_t *nb;

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

    nb->timers = new cu_timers_t();
    snew(nb->timings, 1);

    /* init nbst */
    pmalloc((void**)&nb->nbst.e_lj, sizeof(*nb->nbst.e_lj));
    pmalloc((void**)&nb->nbst.e_el, sizeof(*nb->nbst.e_el));
    pmalloc((void**)&nb->nbst.fshift, SHIFTS * sizeof(*nb->nbst.fshift));

    init_plist(nb->plist[eintLocal]);

    /* set device info, just point it to the right GPU among the detected ones */
    nb->dev_info = deviceInfo;

    /* local/non-local GPU streams */
    stat = cudaStreamCreate(&nb->stream[eintLocal]);
    CU_RET_ERR(stat, "cudaStreamCreate on stream[eintLocal] failed");
    if (nb->bUseTwoStreams)
    {
        init_plist(nb->plist[eintNonlocal]);

        /* Note that the device we're running on does not have to support
         * priorities, because we are querying the priority range which in this
         * case will be a single value.
         */
        int highest_priority;
        stat = cudaDeviceGetStreamPriorityRange(NULL, &highest_priority);
        CU_RET_ERR(stat, "cudaDeviceGetStreamPriorityRange failed");

        stat = cudaStreamCreateWithPriority(&nb->stream[eintNonlocal],
                                            cudaStreamDefault,
                                            highest_priority);
        CU_RET_ERR(stat, "cudaStreamCreateWithPriority on stream[eintNonlocal] failed");
    }

    /* init events for sychronization (timing disabled for performance reasons!) */
    stat = cudaEventCreateWithFlags(&nb->nonlocal_done, cudaEventDisableTiming);
    CU_RET_ERR(stat, "cudaEventCreate on nonlocal_done failed");
    stat = cudaEventCreateWithFlags(&nb->misc_ops_and_local_H2D_done, cudaEventDisableTiming);
    CU_RET_ERR(stat, "cudaEventCreate on misc_ops_and_local_H2D_done failed");

    /* WARNING: CUDA timings are incorrect with multiple streams.
     *          This is the main reason why they are disabled by default.
     */
    // TODO: Consider turning on by default when we can detect nr of streams.
    nb->bDoTime = (getenv("GMX_ENABLE_GPU_TIMING") != NULL);

    if (nb->bDoTime)
    {
        init_timers(nb->timers, nb->bUseTwoStreams);
        init_timings(nb->timings);
    }

    /* set the kernel type for the current GPU */
    /* pick L1 cache configuration */
    nbnxn_cuda_set_cacheconfig(nb->dev_info);

    nbnxn_cuda_init_const(nb, ic, listParams, nbat);

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
    if (canSkipWork(nb, iloc))
    {
        return;
    }

    char          sbuf[STRLEN];
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
        nb->timers->pl_h2d[iloc].openTimingRegion(stream);
        nb->timers->didPairlistH2D[iloc] = true;
    }

    cu_realloc_buffered((void **)&d_plist->sci, h_plist->sci, sizeof(*d_plist->sci),
                        &d_plist->nsci, &d_plist->sci_nalloc,
                        h_plist->nsci,
                        stream, true);

    cu_realloc_buffered((void **)&d_plist->cj4, h_plist->cj4, sizeof(*d_plist->cj4),
                        &d_plist->ncj4, &d_plist->cj4_nalloc,
                        h_plist->ncj4,
                        stream, true);

    /* this call only allocates space on the device (no data is transferred) */
    cu_realloc_buffered((void **)&d_plist->imask, NULL, sizeof(*d_plist->imask),
                        &d_plist->nimask, &d_plist->imask_nalloc,
                        h_plist->ncj4*c_nbnxnGpuClusterpairSplit,
                        stream, true);

    cu_realloc_buffered((void **)&d_plist->excl, h_plist->excl, sizeof(*d_plist->excl),
                        &d_plist->nexcl, &d_plist->excl_nalloc,
                        h_plist->nexcl,
                        stream, true);

    if (bDoTime)
    {
        nb->timers->pl_h2d[iloc].closeTimingRegion(stream);
    }

    /* the next use of thist list we be the first one, so we need to prune */
    d_plist->haveFreshList = true;
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
                             const nbnxn_atomdata_t        *nbat)
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
        timers->atdat.openTimingRegion(ls);
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
            cu_free_buffered(d_atdat->lj_comb);
        }

        stat = cudaMalloc((void **)&d_atdat->f, nalloc*sizeof(*d_atdat->f));
        CU_RET_ERR(stat, "cudaMalloc failed on d_atdat->f");
        stat = cudaMalloc((void **)&d_atdat->xq, nalloc*sizeof(*d_atdat->xq));
        CU_RET_ERR(stat, "cudaMalloc failed on d_atdat->xq");
        if (useLjCombRule(nb->nbparam))
        {
            stat = cudaMalloc((void **)&d_atdat->lj_comb, nalloc*sizeof(*d_atdat->lj_comb));
            CU_RET_ERR(stat, "cudaMalloc failed on d_atdat->lj_comb");
        }
        else
        {
            stat = cudaMalloc((void **)&d_atdat->atom_types, nalloc*sizeof(*d_atdat->atom_types));
            CU_RET_ERR(stat, "cudaMalloc failed on d_atdat->atom_types");
        }

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

    if (useLjCombRule(nb->nbparam))
    {
        cu_copy_H2D_async(d_atdat->lj_comb, nbat->lj_comb,
                          natoms*sizeof(*d_atdat->lj_comb), ls);
    }
    else
    {
        cu_copy_H2D_async(d_atdat->atom_types, nbat->type,
                          natoms*sizeof(*d_atdat->atom_types), ls);
    }

    if (bDoTime)
    {
        timers->atdat.closeTimingRegion(ls);
    }
}

static void nbnxn_cuda_free_nbparam_table(cu_nbparam_t            *nbparam,
                                          const gmx_device_info_t *dev_info)
{
    if (nbparam->eeltype == eelCuEWALD_TAB || nbparam->eeltype == eelCuEWALD_TAB_TWIN)
    {
        destroyParamLookupTable(nbparam->coulomb_tab, nbparam->coulomb_tab_texobj,
                                &nbnxn_cuda_get_coulomb_tab_texref(), dev_info);
    }
}

void nbnxn_gpu_free(gmx_nbnxn_cuda_t *nb)
{
    cudaError_t      stat;
    cu_atomdata_t   *atdat;
    cu_nbparam_t    *nbparam;
    cu_plist_t      *plist, *plist_nl;

    if (nb == NULL)
    {
        return;
    }

    atdat       = nb->atdat;
    nbparam     = nb->nbparam;
    plist       = nb->plist[eintLocal];
    plist_nl    = nb->plist[eintNonlocal];

    nbnxn_cuda_free_nbparam_table(nbparam, nb->dev_info);

    stat = cudaEventDestroy(nb->nonlocal_done);
    CU_RET_ERR(stat, "cudaEventDestroy failed on timers->nonlocal_done");
    stat = cudaEventDestroy(nb->misc_ops_and_local_H2D_done);
    CU_RET_ERR(stat, "cudaEventDestroy failed on timers->misc_ops_and_local_H2D_done");

    delete nb->timers;
    if (nb->bDoTime)
    {
        /* The non-local counters/stream (second in the array) are needed only with DD. */
        for (int i = 0; i <= (nb->bUseTwoStreams ? 1 : 0); i++)
        {
            stat = cudaStreamDestroy(nb->stream[i]);
            CU_RET_ERR(stat, "cudaStreamDestroy failed on stream");
        }
    }

    if (!useLjCombRule(nb->nbparam))
    {
        destroyParamLookupTable(nbparam->nbfp, nbparam->nbfp_texobj,
                                &nbnxn_cuda_get_nbfp_texref(), nb->dev_info);

    }

    if (nbparam->vdwtype == evdwCuEWALDGEOM || nbparam->vdwtype == evdwCuEWALDLB)
    {
        destroyParamLookupTable(nbparam->nbfp_comb, nbparam->nbfp_comb_texobj,
                                &nbnxn_cuda_get_nbfp_comb_texref(), nb->dev_info);
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
    cu_free_buffered(atdat->lj_comb);

    cu_free_buffered(plist->sci, &plist->nsci, &plist->sci_nalloc);
    cu_free_buffered(plist->cj4, &plist->ncj4, &plist->cj4_nalloc);
    cu_free_buffered(plist->imask, &plist->nimask, &plist->imask_nalloc);
    cu_free_buffered(plist->excl, &plist->nexcl, &plist->excl_nalloc);
    if (nb->bUseTwoStreams)
    {
        cu_free_buffered(plist_nl->sci, &plist_nl->nsci, &plist_nl->sci_nalloc);
        cu_free_buffered(plist_nl->cj4, &plist_nl->ncj4, &plist_nl->cj4_nalloc);
        cu_free_buffered(plist_nl->imask, &plist_nl->nimask, &plist_nl->imask_nalloc);
        cu_free_buffered(plist_nl->excl, &plist_nl->nexcl, &plist->excl_nalloc);
    }

    sfree(atdat);
    sfree(nbparam);
    sfree(plist);
    if (nb->bUseTwoStreams)
    {
        sfree(plist_nl);
    }
    sfree(nb->timings);
    sfree(nb);

    if (debug)
    {
        fprintf(debug, "Cleaned up CUDA data structures.\n");
    }
}

//! This function is documented in the header file
gmx_wallclock_gpu_nbnxn_t *nbnxn_gpu_get_timings(gmx_nbnxn_cuda_t *nb)
{
    return (nb != nullptr && nb->bDoTime) ? nb->timings : nullptr;
}

void nbnxn_gpu_reset_timings(nonbonded_verlet_t* nbv)
{
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
