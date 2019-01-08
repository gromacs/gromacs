/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2019, by the GROMACS development team, led by
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
/*! \internal \file
 * \brief Common functions for the different NBNXN GPU implementations.
 *
 * \author Berk Hess <hess@kth.se>
 *
 * \ingroup module_nbnxm
 */

#include "gmxpre.h"

#include "gromacs/compat/make_unique.h"
#include "gromacs/domdec/domdec.h"
#include "gromacs/domdec/domdec_struct.h"
#include "gromacs/hardware/hw_info.h"
#include "gromacs/mdlib/gmx_omp_nthreads.h"
#include "gromacs/mdtypes/commrec.h"
#include "gromacs/mdtypes/forcerec.h"
#include "gromacs/mdtypes/inputrec.h"
#include "gromacs/nbnxm/atomdata.h"
#include "gromacs/nbnxm/gpu_data_mgmt.h"
#include "gromacs/nbnxm/nbnxm.h"
#include "gromacs/nbnxm/nbnxm_simd.h"
#include "gromacs/nbnxm/tuning.h"
#include "gromacs/nbnxm/utility.h"
#include "gromacs/simd/simd.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/logger.h"

#include "grid.h"
#include "internal.h"

/*! \brief Returns whether CPU SIMD support exists for the given inputrec
 *
 * If the return value is FALSE and fplog/cr != NULL, prints a fallback
 * message to fplog/stderr.
 */
static gmx_bool nbnxn_simd_supported(const gmx::MDLogger &mdlog,
                                     const t_inputrec    *ir)
{
    if (ir->vdwtype == evdwPME && ir->ljpme_combination_rule == eljpmeLB)
    {
        /* LJ PME with LB combination rule does 7 mesh operations.
         * This so slow that we don't compile SIMD non-bonded kernels
         * for that. */
        GMX_LOG(mdlog.warning).asParagraph().appendText("LJ-PME with Lorentz-Berthelot is not supported with SIMD kernels, falling back to plain C kernels");
        return FALSE;
    }

    return TRUE;
}

/*! \brief Returns the most suitable CPU kernel type and Ewald handling */
static void pick_nbnxn_kernel_cpu(const t_inputrec gmx_unused    *ir,
                                  int                            *kernel_type,
                                  int                            *ewald_excl,
                                  const gmx_hw_info_t gmx_unused &hardwareInfo)
{
    *kernel_type = nbnxnk4x4_PlainC;
    *ewald_excl  = ewaldexclTable;

#if GMX_SIMD
    {
#ifdef GMX_NBNXN_SIMD_4XN
        *kernel_type = nbnxnk4xN_SIMD_4xN;
#endif
#ifdef GMX_NBNXN_SIMD_2XNN
        *kernel_type = nbnxnk4xN_SIMD_2xNN;
#endif

#if defined GMX_NBNXN_SIMD_2XNN && defined GMX_NBNXN_SIMD_4XN
        /* We need to choose if we want 2x(N+N) or 4xN kernels.
         * This is based on the SIMD acceleration choice and CPU information
         * detected at runtime.
         *
         * 4xN calculates more (zero) interactions, but has less pair-search
         * work and much better kernel instruction scheduling.
         *
         * Up till now we have only seen that on Intel Sandy/Ivy Bridge,
         * which doesn't have FMA, both the analytical and tabulated Ewald
         * kernels have similar pair rates for 4x8 and 2x(4+4), so we choose
         * 2x(4+4) because it results in significantly fewer pairs.
         * For RF, the raw pair rate of the 4x8 kernel is higher than 2x(4+4),
         * 10% with HT, 50% without HT. As we currently don't detect the actual
         * use of HT, use 4x8 to avoid a potential performance hit.
         * On Intel Haswell 4x8 is always faster.
         */
        *kernel_type = nbnxnk4xN_SIMD_4xN;

#if !GMX_SIMD_HAVE_FMA
        if (EEL_PME_EWALD(ir->coulombtype) ||
            EVDW_PME(ir->vdwtype))
        {
            /* We have Ewald kernels without FMA (Intel Sandy/Ivy Bridge).
             * There are enough instructions to make 2x(4+4) efficient.
             */
            *kernel_type = nbnxnk4xN_SIMD_2xNN;
        }
#endif
        if (hardwareInfo.haveAmdZenCpu)
        {
            /* One 256-bit FMA per cycle makes 2xNN faster */
            *kernel_type = nbnxnk4xN_SIMD_2xNN;
        }
#endif  /* GMX_NBNXN_SIMD_2XNN && GMX_NBNXN_SIMD_4XN */


        if (getenv("GMX_NBNXN_SIMD_4XN") != nullptr)
        {
#ifdef GMX_NBNXN_SIMD_4XN
            *kernel_type = nbnxnk4xN_SIMD_4xN;
#else
            gmx_fatal(FARGS, "SIMD 4xN kernels requested, but GROMACS has been compiled without support for these kernels");
#endif
        }
        if (getenv("GMX_NBNXN_SIMD_2XNN") != nullptr)
        {
#ifdef GMX_NBNXN_SIMD_2XNN
            *kernel_type = nbnxnk4xN_SIMD_2xNN;
#else
            gmx_fatal(FARGS, "SIMD 2x(N+N) kernels requested, but GROMACS has been compiled without support for these kernels");
#endif
        }

        /* Analytical Ewald exclusion correction is only an option in
         * the SIMD kernel.
         * Since table lookup's don't parallelize with SIMD, analytical
         * will probably always be faster for a SIMD width of 8 or more.
         * With FMA analytical is sometimes faster for a width if 4 as well.
         * In single precision, this is faster on Bulldozer.
         */
#if GMX_SIMD_REAL_WIDTH >= 8 || \
        (GMX_SIMD_REAL_WIDTH >= 4 && GMX_SIMD_HAVE_FMA && !GMX_DOUBLE)
        /* On AMD Zen, tabulated Ewald kernels are faster on all 4 combinations
         * of single or double precision and 128 or 256-bit AVX2.
         */
        if (!hardwareInfo.haveAmdZenCpu)
        {
            *ewald_excl = ewaldexclAnalytical;
        }
#endif
        if (getenv("GMX_NBNXN_EWALD_TABLE") != nullptr)
        {
            *ewald_excl = ewaldexclTable;
        }
        if (getenv("GMX_NBNXN_EWALD_ANALYTICAL") != nullptr)
        {
            *ewald_excl = ewaldexclAnalytical;
        }

    }
#endif // GMX_SIMD
}

const char *lookup_nbnxn_kernel_name(int kernel_type)
{
    const char *returnvalue = nullptr;
    switch (kernel_type)
    {
        case nbnxnkNotSet:
            returnvalue = "not set";
            break;
        case nbnxnk4x4_PlainC:
            returnvalue = "plain C";
            break;
        case nbnxnk4xN_SIMD_4xN:
        case nbnxnk4xN_SIMD_2xNN:
#if GMX_SIMD
            returnvalue = "SIMD";
#else  // GMX_SIMD
            returnvalue = "not available";
#endif // GMX_SIMD
            break;
        case nbnxnk8x8x8_GPU: returnvalue    = "GPU"; break;
        case nbnxnk8x8x8_PlainC: returnvalue = "plain C"; break;

        case nbnxnkNR:
        default:
            gmx_fatal(FARGS, "Illegal kernel type selected");
    }
    return returnvalue;
};

/*! \brief Returns the most suitable kernel type and Ewald handling */
static void pick_nbnxn_kernel(const gmx::MDLogger &mdlog,
                              gmx_bool             use_simd_kernels,
                              const gmx_hw_info_t &hardwareInfo,
                              gmx_bool             bUseGPU,
                              EmulateGpuNonbonded  emulateGpu,
                              const t_inputrec    *ir,
                              int                 *kernel_type,
                              int                 *ewald_excl,
                              gmx_bool             bDoNonbonded)
{
    GMX_RELEASE_ASSERT(kernel_type, "Need a valid kernel_type pointer");

    *kernel_type = nbnxnkNotSet;
    *ewald_excl  = ewaldexclTable;

    if (emulateGpu == EmulateGpuNonbonded::Yes)
    {
        *kernel_type = nbnxnk8x8x8_PlainC;

        if (bDoNonbonded)
        {
            GMX_LOG(mdlog.warning).asParagraph().appendText("Emulating a GPU run on the CPU (slow)");
        }
    }
    else if (bUseGPU)
    {
        *kernel_type = nbnxnk8x8x8_GPU;
    }

    if (*kernel_type == nbnxnkNotSet)
    {
        if (use_simd_kernels &&
            nbnxn_simd_supported(mdlog, ir))
        {
            pick_nbnxn_kernel_cpu(ir, kernel_type, ewald_excl, hardwareInfo);
        }
        else
        {
            *kernel_type = nbnxnk4x4_PlainC;
        }
    }

    if (bDoNonbonded)
    {
        GMX_LOG(mdlog.info).asParagraph().appendTextFormatted(
                "Using %s %dx%d nonbonded short-range kernels",
                lookup_nbnxn_kernel_name(*kernel_type),
                nbnxn_kernel_to_cluster_i_size(*kernel_type),
                nbnxn_kernel_to_cluster_j_size(*kernel_type));

        if (nbnxnk4x4_PlainC == *kernel_type ||
            nbnxnk8x8x8_PlainC == *kernel_type)
        {
            GMX_LOG(mdlog.warning).asParagraph().appendTextFormatted(
                    "WARNING: Using the slow %s kernels. This should\n"
                    "not happen during routine usage on supported platforms.",
                    lookup_nbnxn_kernel_name(*kernel_type));
        }
    }
}

void init_nb_verlet(const gmx::MDLogger     &mdlog,
                    nonbonded_verlet_t     **nb_verlet,
                    gmx_bool                 bFEP_NonBonded,
                    const t_inputrec        *ir,
                    const t_forcerec        *fr,
                    const t_commrec         *cr,
                    const gmx_hw_info_t     &hardwareInfo,
                    const gmx_device_info_t *deviceInfo,
                    const gmx_mtop_t        *mtop,
                    matrix                   box)
{
    nonbonded_verlet_t *nbv;
    char               *env;

    nbv = new nonbonded_verlet_t();

    nbv->emulateGpu = ((getenv("GMX_EMULATE_GPU") != nullptr) ? EmulateGpuNonbonded::Yes : EmulateGpuNonbonded::No);
    nbv->bUseGPU    = deviceInfo != nullptr;

    GMX_RELEASE_ASSERT(!(nbv->emulateGpu == EmulateGpuNonbonded::Yes && nbv->bUseGPU), "When GPU emulation is active, there cannot be a GPU assignment");

    nbv->nbs             = nullptr;
    nbv->min_ci_balanced = 0;

    nbv->ngrp = (DOMAINDECOMP(cr) ? 2 : 1);
    for (int i = 0; i < nbv->ngrp; i++)
    {
        nbv->grp[i].nbl_lists.nnbl = 0;
        nbv->grp[i].kernel_type    = nbnxnkNotSet;

        if (i == 0) /* local */
        {
            pick_nbnxn_kernel(mdlog, fr->use_simd_kernels, hardwareInfo,
                              nbv->bUseGPU, nbv->emulateGpu, ir,
                              &nbv->grp[i].kernel_type,
                              &nbv->grp[i].ewald_excl,
                              fr->bNonbonded);
        }
        else /* non-local */
        {
            /* Use the same kernel for local and non-local interactions */
            nbv->grp[i].kernel_type = nbv->grp[0].kernel_type;
            nbv->grp[i].ewald_excl  = nbv->grp[0].ewald_excl;
        }
    }

    nbv->listParams = gmx::compat::make_unique<NbnxnListParameters>(ir->rlist);
    setupDynamicPairlistPruning(mdlog, ir, mtop, box, nbv->grp[0].kernel_type, fr->ic,
                                nbv->listParams.get());

    nbv->nbs = gmx::compat::make_unique<nbnxn_search>(DOMAINDECOMP(cr) ? &cr->dd->nc : nullptr,
                                                      DOMAINDECOMP(cr) ? domdec_zones(cr->dd) : nullptr,
                                                      bFEP_NonBonded,
                                                      gmx_omp_nthreads_get(emntPairsearch));

    for (int i = 0; i < nbv->ngrp; i++)
    {
        nbnxn_init_pairlist_set(&nbv->grp[i].nbl_lists,
                                nbnxn_kernel_pairlist_simple(nbv->grp[i].kernel_type),
                                /* 8x8x8 "non-simple" lists are ATM always combined */
                                !nbnxn_kernel_pairlist_simple(nbv->grp[i].kernel_type));
    }

    int      enbnxninitcombrule;
    if (fr->ic->vdwtype == evdwCUT &&
        (fr->ic->vdw_modifier == eintmodNONE ||
         fr->ic->vdw_modifier == eintmodPOTSHIFT) &&
        getenv("GMX_NO_LJ_COMB_RULE") == nullptr)
    {
        /* Plain LJ cut-off: we can optimize with combination rules */
        enbnxninitcombrule = enbnxninitcombruleDETECT;
    }
    else if (fr->ic->vdwtype == evdwPME)
    {
        /* LJ-PME: we need to use a combination rule for the grid */
        if (fr->ljpme_combination_rule == eljpmeGEOM)
        {
            enbnxninitcombrule = enbnxninitcombruleGEOM;
        }
        else
        {
            enbnxninitcombrule = enbnxninitcombruleLB;
        }
    }
    else
    {
        /* We use a full combination matrix: no rule required */
        enbnxninitcombrule = enbnxninitcombruleNONE;
    }

    nbv->nbat = new nbnxn_atomdata_t(nbv->bUseGPU ? gmx::PinningPolicy::PinnedIfSupported : gmx::PinningPolicy::CannotBePinned);
    int mimimumNumEnergyGroupNonbonded = ir->opts.ngener;
    if (ir->opts.ngener - ir->nwall == 1)
    {
        /* We have only one non-wall energy group, we do not need energy group
         * support in the non-bondeds kernels, since all non-bonded energy
         * contributions go to the first element of the energy group matrix.
         */
        mimimumNumEnergyGroupNonbonded = 1;
    }
    bool bSimpleList = nbnxn_kernel_pairlist_simple(nbv->grp[0].kernel_type);
    nbnxn_atomdata_init(mdlog,
                        nbv->nbat,
                        nbv->grp[0].kernel_type,
                        enbnxninitcombrule,
                        fr->ntype, fr->nbfp,
                        mimimumNumEnergyGroupNonbonded,
                        bSimpleList ? gmx_omp_nthreads_get(emntNonbonded) : 1);

    if (nbv->bUseGPU)
    {
        /* init the NxN GPU data; the last argument tells whether we'll have
         * both local and non-local NB calculation on GPU */
        nbnxn_gpu_init(&nbv->gpu_nbv,
                       deviceInfo,
                       fr->ic,
                       nbv->listParams.get(),
                       nbv->nbat,
                       cr->nodeid,
                       (nbv->ngrp > 1));

        if ((env = getenv("GMX_NB_MIN_CI")) != nullptr)
        {
            char *end;

            nbv->min_ci_balanced = strtol(env, &end, 10);
            if (!end || (*end != 0) || nbv->min_ci_balanced < 0)
            {
                gmx_fatal(FARGS, "Invalid value passed in GMX_NB_MIN_CI=%s, non-negative integer required", env);
            }

            if (debug)
            {
                fprintf(debug, "Neighbor-list balancing parameter: %d (passed as env. var.)\n",
                        nbv->min_ci_balanced);
            }
        }
        else
        {
            nbv->min_ci_balanced = nbnxn_gpu_min_ci_balanced(nbv->gpu_nbv);
            if (debug)
            {
                fprintf(debug, "Neighbor-list balancing parameter: %d (auto-adjusted to the number of GPU multi-processors)\n",
                        nbv->min_ci_balanced);
            }
        }

    }

    *nb_verlet = nbv;
}
