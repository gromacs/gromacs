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

#include "gromacs/domdec/domdec.h"
#include "gromacs/domdec/domdec_struct.h"
#include "gromacs/hardware/hw_info.h"
#include "gromacs/mdlib/gmx_omp_nthreads.h"
#include "gromacs/mdtypes/commrec.h"
#include "gromacs/mdtypes/forcerec.h"
#include "gromacs/mdtypes/inputrec.h"
#include "gromacs/nbnxm/gpu_data_mgmt.h"
#include "gromacs/nbnxm/nbnxm.h"
#include "gromacs/nbnxm/pairlist_tuning.h"
#include "gromacs/simd/simd.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/logger.h"

#include "atomdata.h"
#include "gpu_types.h"
#include "grid.h"
#include "nbnxm_geometry.h"
#include "nbnxm_simd.h"
#include "pairlist.h"
#include "pairlistset.h"
#include "pairlistsets.h"
#include "pairsearch.h"

namespace Nbnxm
{

/*! \brief Resources that can be used to execute non-bonded kernels on */
enum class NonbondedResource : int
{
    Cpu,
    Gpu,
    EmulateGpu
};

/*! \brief Returns whether CPU SIMD support exists for the given inputrec
 *
 * If the return value is FALSE and fplog/cr != NULL, prints a fallback
 * message to fplog/stderr.
 */
static gmx_bool nbnxn_simd_supported(const gmx::MDLogger& mdlog, const t_inputrec* ir)
{
    if (ir->vdwtype == evdwPME && ir->ljpme_combination_rule == eljpmeLB)
    {
        /* LJ PME with LB combination rule does 7 mesh operations.
         * This so slow that we don't compile SIMD non-bonded kernels
         * for that. */
        GMX_LOG(mdlog.warning)
                .asParagraph()
                .appendText(
                        "LJ-PME with Lorentz-Berthelot is not supported with SIMD kernels, falling "
                        "back to plain C kernels");
        return FALSE;
    }

    return TRUE;
}

/*! \brief Returns the most suitable CPU kernel type and Ewald handling */
static KernelSetup pick_nbnxn_kernel_cpu(const t_inputrec gmx_unused* ir,
                                         const gmx_hw_info_t gmx_unused& hardwareInfo)
{
    KernelSetup kernelSetup;

    if (!GMX_SIMD)
    {
        kernelSetup.kernelType         = KernelType::Cpu4x4_PlainC;
        kernelSetup.ewaldExclusionType = EwaldExclusionType::Table;
    }
    else
    {
#ifdef GMX_NBNXN_SIMD_4XN
        kernelSetup.kernelType = KernelType::Cpu4xN_Simd_4xN;
#endif
#ifdef GMX_NBNXN_SIMD_2XNN
        kernelSetup.kernelType = KernelType::Cpu4xN_Simd_2xNN;
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
        kernelSetup.kernelType = KernelType::Cpu4xN_Simd_4xN;

        if (!GMX_SIMD_HAVE_FMA && (EEL_PME_EWALD(ir->coulombtype) || EVDW_PME(ir->vdwtype)))
        {
            /* We have Ewald kernels without FMA (Intel Sandy/Ivy Bridge).
             * There are enough instructions to make 2x(4+4) efficient.
             */
            kernelSetup.kernelType = KernelType::Cpu4xN_Simd_2xNN;
        }

        if (hardwareInfo.haveAmdZen1Cpu)
        {
            /* One 256-bit FMA per cycle makes 2xNN faster */
            kernelSetup.kernelType = KernelType::Cpu4xN_Simd_2xNN;
        }
#endif /* GMX_NBNXN_SIMD_2XNN && GMX_NBNXN_SIMD_4XN */


        if (getenv("GMX_NBNXN_SIMD_4XN") != nullptr)
        {
#ifdef GMX_NBNXN_SIMD_4XN
            kernelSetup.kernelType = KernelType::Cpu4xN_Simd_4xN;
#else
            gmx_fatal(FARGS,
                      "SIMD 4xN kernels requested, but GROMACS has been compiled without support "
                      "for these kernels");
#endif
        }
        if (getenv("GMX_NBNXN_SIMD_2XNN") != nullptr)
        {
#ifdef GMX_NBNXN_SIMD_2XNN
            kernelSetup.kernelType = KernelType::Cpu4xN_Simd_2xNN;
#else
            gmx_fatal(FARGS,
                      "SIMD 2x(N+N) kernels requested, but GROMACS has been compiled without "
                      "support for these kernels");
#endif
        }

        /* Analytical Ewald exclusion correction is only an option in
         * the SIMD kernel.
         * Since table lookup's don't parallelize with SIMD, analytical
         * will probably always be faster for a SIMD width of 8 or more.
         * With FMA analytical is sometimes faster for a width if 4 as well.
         * In single precision, this is faster on Bulldozer.
         * On AMD Zen, tabulated Ewald kernels are faster on all 4 combinations
         * of single or double precision and 128 or 256-bit AVX2.
         */
        if (
#if GMX_SIMD
                (GMX_SIMD_REAL_WIDTH >= 8 || (GMX_SIMD_REAL_WIDTH >= 4 && GMX_SIMD_HAVE_FMA && !GMX_DOUBLE)) &&
#endif
                !hardwareInfo.haveAmdZen1Cpu)
        {
            kernelSetup.ewaldExclusionType = EwaldExclusionType::Analytical;
        }
        else
        {
            kernelSetup.ewaldExclusionType = EwaldExclusionType::Table;
        }
        if (getenv("GMX_NBNXN_EWALD_TABLE") != nullptr)
        {
            kernelSetup.ewaldExclusionType = EwaldExclusionType::Table;
        }
        if (getenv("GMX_NBNXN_EWALD_ANALYTICAL") != nullptr)
        {
            kernelSetup.ewaldExclusionType = EwaldExclusionType::Analytical;
        }
    }

    return kernelSetup;
}

const char* lookup_kernel_name(const KernelType kernelType)
{
    const char* returnvalue = nullptr;
    switch (kernelType)
    {
        case KernelType::NotSet: returnvalue = "not set"; break;
        case KernelType::Cpu4x4_PlainC: returnvalue = "plain C"; break;
        case KernelType::Cpu4xN_Simd_4xN:
        case KernelType::Cpu4xN_Simd_2xNN:
#if GMX_SIMD
            returnvalue = "SIMD";
#else  // GMX_SIMD
            returnvalue = "not available";
#endif // GMX_SIMD
            break;
        case KernelType::Gpu8x8x8: returnvalue = "GPU"; break;
        case KernelType::Cpu8x8x8_PlainC: returnvalue = "plain C"; break;

        default: gmx_fatal(FARGS, "Illegal kernel type selected");
    }
    return returnvalue;
};

/*! \brief Returns the most suitable kernel type and Ewald handling */
static KernelSetup pick_nbnxn_kernel(const gmx::MDLogger&     mdlog,
                                     gmx_bool                 use_simd_kernels,
                                     const gmx_hw_info_t&     hardwareInfo,
                                     const NonbondedResource& nonbondedResource,
                                     const t_inputrec*        ir,
                                     gmx_bool                 bDoNonbonded)
{
    KernelSetup kernelSetup;

    if (nonbondedResource == NonbondedResource::EmulateGpu)
    {
        kernelSetup.kernelType         = KernelType::Cpu8x8x8_PlainC;
        kernelSetup.ewaldExclusionType = EwaldExclusionType::DecidedByGpuModule;

        if (bDoNonbonded)
        {
            GMX_LOG(mdlog.warning)
                    .asParagraph()
                    .appendText("Emulating a GPU run on the CPU (slow)");
        }
    }
    else if (nonbondedResource == NonbondedResource::Gpu)
    {
        kernelSetup.kernelType         = KernelType::Gpu8x8x8;
        kernelSetup.ewaldExclusionType = EwaldExclusionType::DecidedByGpuModule;
    }
    else
    {
        if (use_simd_kernels && nbnxn_simd_supported(mdlog, ir))
        {
            kernelSetup = pick_nbnxn_kernel_cpu(ir, hardwareInfo);
        }
        else
        {
            kernelSetup.kernelType         = KernelType::Cpu4x4_PlainC;
            kernelSetup.ewaldExclusionType = EwaldExclusionType::Analytical;
        }
    }

    if (bDoNonbonded)
    {
        GMX_LOG(mdlog.info)
                .asParagraph()
                .appendTextFormatted("Using %s %dx%d nonbonded short-range kernels",
                                     lookup_kernel_name(kernelSetup.kernelType),
                                     IClusterSizePerKernelType[kernelSetup.kernelType],
                                     JClusterSizePerKernelType[kernelSetup.kernelType]);

        if (KernelType::Cpu4x4_PlainC == kernelSetup.kernelType
            || KernelType::Cpu8x8x8_PlainC == kernelSetup.kernelType)
        {
            GMX_LOG(mdlog.warning)
                    .asParagraph()
                    .appendTextFormatted(
                            "WARNING: Using the slow %s kernels. This should\n"
                            "not happen during routine usage on supported platforms.",
                            lookup_kernel_name(kernelSetup.kernelType));
        }
    }

    GMX_RELEASE_ASSERT(kernelSetup.kernelType != KernelType::NotSet
                               && kernelSetup.ewaldExclusionType != EwaldExclusionType::NotSet,
                       "All kernel setup parameters should be set here");

    return kernelSetup;
}

} // namespace Nbnxm

PairlistSets::PairlistSets(const PairlistParams& pairlistParams,
                           const bool            haveMultipleDomains,
                           const int             minimumIlistCountForGpuBalancing) :
    params_(pairlistParams),
    minimumIlistCountForGpuBalancing_(minimumIlistCountForGpuBalancing)
{
    localSet_ = std::make_unique<PairlistSet>(gmx::InteractionLocality::Local, params_);

    if (haveMultipleDomains)
    {
        nonlocalSet_ = std::make_unique<PairlistSet>(gmx::InteractionLocality::NonLocal, params_);
    }
}

namespace Nbnxm
{

/*! \brief Gets and returns the minimum i-list count for balacing based on the GPU used or env.var. when set */
static int getMinimumIlistCountForGpuBalancing(gmx_nbnxn_gpu_t* nbnxmGpu)
{
    int minimumIlistCount;

    if (const char* env = getenv("GMX_NB_MIN_CI"))
    {
        char* end;

        minimumIlistCount = strtol(env, &end, 10);
        if (!end || (*end != 0) || minimumIlistCount < 0)
        {
            gmx_fatal(FARGS,
                      "Invalid value passed in GMX_NB_MIN_CI=%s, non-negative integer required", env);
        }

        if (debug)
        {
            fprintf(debug, "Neighbor-list balancing parameter: %d (passed as env. var.)\n",
                    minimumIlistCount);
        }
    }
    else
    {
        minimumIlistCount = gpu_min_ci_balanced(nbnxmGpu);
        if (debug)
        {
            fprintf(debug,
                    "Neighbor-list balancing parameter: %d (auto-adjusted to the number of GPU "
                    "multi-processors)\n",
                    minimumIlistCount);
        }
    }

    return minimumIlistCount;
}

std::unique_ptr<nonbonded_verlet_t> init_nb_verlet(const gmx::MDLogger&     mdlog,
                                                   gmx_bool                 bFEP_NonBonded,
                                                   const t_inputrec*        ir,
                                                   const t_forcerec*        fr,
                                                   const t_commrec*         cr,
                                                   const gmx_hw_info_t&     hardwareInfo,
                                                   const gmx_device_info_t* deviceInfo,
                                                   const gmx_mtop_t*        mtop,
                                                   matrix                   box,
                                                   gmx_wallcycle*           wcycle)
{
    const bool emulateGpu = (getenv("GMX_EMULATE_GPU") != nullptr);
    const bool useGpu     = deviceInfo != nullptr;

    GMX_RELEASE_ASSERT(!(emulateGpu && useGpu),
                       "When GPU emulation is active, there cannot be a GPU assignment");

    NonbondedResource nonbondedResource;
    if (useGpu)
    {
        nonbondedResource = NonbondedResource::Gpu;
    }
    else if (emulateGpu)
    {
        nonbondedResource = NonbondedResource::EmulateGpu;
    }
    else
    {
        nonbondedResource = NonbondedResource::Cpu;
    }

    Nbnxm::KernelSetup kernelSetup = pick_nbnxn_kernel(mdlog, fr->use_simd_kernels, hardwareInfo,
                                                       nonbondedResource, ir, fr->bNonbonded);

    const bool haveMultipleDomains = (DOMAINDECOMP(cr) && cr->dd->nnodes > 1);

    PairlistParams pairlistParams(kernelSetup.kernelType, bFEP_NonBonded, ir->rlist,
                                  havePPDomainDecomposition(cr));

    setupDynamicPairlistPruning(mdlog, ir, mtop, box, fr->ic, &pairlistParams);

    int enbnxninitcombrule;
    if (fr->ic->vdwtype == evdwCUT
        && (fr->ic->vdw_modifier == eintmodNONE || fr->ic->vdw_modifier == eintmodPOTSHIFT)
        && getenv("GMX_NO_LJ_COMB_RULE") == nullptr)
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

    auto pinPolicy = (useGpu ? gmx::PinningPolicy::PinnedIfSupported : gmx::PinningPolicy::CannotBePinned);

    auto nbat = std::make_unique<nbnxn_atomdata_t>(pinPolicy);

    int mimimumNumEnergyGroupNonbonded = ir->opts.ngener;
    if (ir->opts.ngener - ir->nwall == 1)
    {
        /* We have only one non-wall energy group, we do not need energy group
         * support in the non-bondeds kernels, since all non-bonded energy
         * contributions go to the first element of the energy group matrix.
         */
        mimimumNumEnergyGroupNonbonded = 1;
    }
    nbnxn_atomdata_init(mdlog, nbat.get(), kernelSetup.kernelType, enbnxninitcombrule, fr->ntype,
                        fr->nbfp, mimimumNumEnergyGroupNonbonded,
                        (useGpu || emulateGpu) ? 1 : gmx_omp_nthreads_get(emntNonbonded));

    gmx_nbnxn_gpu_t* gpu_nbv                          = nullptr;
    int              minimumIlistCountForGpuBalancing = 0;
    if (useGpu)
    {
        /* init the NxN GPU data; the last argument tells whether we'll have
         * both local and non-local NB calculation on GPU */
        gpu_nbv = gpu_init(deviceInfo, fr->ic, pairlistParams, nbat.get(), cr->nodeid, haveMultipleDomains);

        minimumIlistCountForGpuBalancing = getMinimumIlistCountForGpuBalancing(gpu_nbv);
    }

    auto pairlistSets = std::make_unique<PairlistSets>(pairlistParams, haveMultipleDomains,
                                                       minimumIlistCountForGpuBalancing);

    auto pairSearch = std::make_unique<PairSearch>(
            ir->ePBC, EI_TPI(ir->eI), DOMAINDECOMP(cr) ? &cr->dd->nc : nullptr,
            DOMAINDECOMP(cr) ? domdec_zones(cr->dd) : nullptr, pairlistParams.pairlistType,
            bFEP_NonBonded, gmx_omp_nthreads_get(emntPairsearch), pinPolicy);

    return std::make_unique<nonbonded_verlet_t>(std::move(pairlistSets), std::move(pairSearch),
                                                std::move(nbat), kernelSetup, gpu_nbv, wcycle);
}

} // namespace Nbnxm

nonbonded_verlet_t::nonbonded_verlet_t(std::unique_ptr<PairlistSets>     pairlistSets,
                                       std::unique_ptr<PairSearch>       pairSearch,
                                       std::unique_ptr<nbnxn_atomdata_t> nbat_in,
                                       const Nbnxm::KernelSetup&         kernelSetup,
                                       gmx_nbnxn_gpu_t*                  gpu_nbv_ptr,
                                       gmx_wallcycle*                    wcycle) :
    pairlistSets_(std::move(pairlistSets)),
    pairSearch_(std::move(pairSearch)),
    nbat(std::move(nbat_in)),
    kernelSetup_(kernelSetup),
    wcycle_(wcycle),
    gpu_nbv(gpu_nbv_ptr)
{
    GMX_RELEASE_ASSERT(pairlistSets_, "Need valid pairlistSets");
    GMX_RELEASE_ASSERT(pairSearch_, "Need valid search object");
    GMX_RELEASE_ASSERT(nbat, "Need valid atomdata object");
}

nonbonded_verlet_t::~nonbonded_verlet_t()
{
    Nbnxm::gpu_free(gpu_nbv);
}
