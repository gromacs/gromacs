/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2019- The GROMACS Authors
 * and the project initiators Erik Lindahl, Berk Hess and David van der Spoel.
 * Consult the AUTHORS/COPYING files and https://www.gromacs.org for details.
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
 * https://www.gnu.org/licenses, or write to the Free Software Foundation,
 * Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA.
 *
 * If you want to redistribute modifications to GROMACS, please
 * consider that scientific software is very special. Version
 * control is crucial - bugs must be traceable. We will be happy to
 * consider code for inclusion in the official distribution, but
 * derived work must not be called official GROMACS. Details are found
 * in the README & COPYING files - if they are missing, get the
 * official version at https://www.gromacs.org.
 *
 * To help us fund GROMACS development, we humbly ask that you cite
 * the research papers on the package. Check out https://www.gromacs.org.
 */
/*! \internal \file
 * \brief Common functions for the different NBNXN GPU implementations.
 *
 * \author Berk Hess <hess@kth.se>
 *
 * \ingroup module_nbnxm
 */

#include "gmxpre.h"

#include <cstdio>
#include <cstdlib>

#include <algorithm>
#include <filesystem>
#include <memory>
#include <optional>
#include <utility>

#include "gromacs/domdec/domdec.h"
#include "gromacs/domdec/domdec_struct.h"
#include "gromacs/gpu_utils/hostallocator.h"
#include "gromacs/hardware/hw_info.h"
#include "gromacs/math/vectypes.h"
#include "gromacs/mdlib/calc_verletbuf.h"
#include "gromacs/mdlib/gmx_omp_nthreads.h"
#include "gromacs/mdtypes/commrec.h"
#include "gromacs/mdtypes/enerdata.h"
#include "gromacs/mdtypes/forcerec.h"
#include "gromacs/mdtypes/inputrec.h"
#include "gromacs/mdtypes/interaction_const.h"
#include "gromacs/mdtypes/md_enums.h"
#include "gromacs/nbnxm/atomdata.h"
#include "gromacs/nbnxm/gpu_data_mgmt.h"
#include "gromacs/nbnxm/nbnxm.h"
#include "gromacs/nbnxm/pairlist_tuning.h"
#include "gromacs/nbnxm/pairlistparams.h"
#include "gromacs/simd/simd.h"
#include "gromacs/topology/mtop_util.h"
#include "gromacs/utility/arrayref.h"
#include "gromacs/utility/basedefinitions.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/gmxassert.h"
#include "gromacs/utility/logger.h"
#include "gromacs/utility/real.h"

#include "exclusionchecker.h"
#include "freeenergydispatch.h"
#include "grid.h"
#include "nbnxm_geometry.h"
#include "nbnxm_simd.h"
#include "pairlist.h"
#include "pairlistset.h"
#include "pairlistsets.h"
#include "pairsearch.h"

struct gmx_mtop_t;
struct gmx_wallcycle;

namespace gmx
{
class DeviceStreamManager;
class ObservablesReducerBuilder;
struct NbnxmGpu;

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
static bool nbnxn_simd_supported(const MDLogger& mdlog, const t_inputrec& inputrec)
{
    if (inputrec.vdwtype == VanDerWaalsType::Pme && inputrec.ljpme_combination_rule == LongRangeVdW::LB)
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
static NbnxmKernelSetup pick_nbnxn_kernel_cpu(const t_inputrec gmx_unused& inputrec,
                                              const gmx_hw_info_t gmx_unused& hardwareInfo)
{
    NbnxmKernelSetup kernelSetup;

    if (!GMX_SIMD)
    {
        kernelSetup.kernelType         = NbnxmKernelType::Cpu4x4_PlainC;
        kernelSetup.ewaldExclusionType = EwaldExclusionType::Table;
    }
    else if (sc_haveNbnxmSimd4xmKernels && !sc_haveNbnxmSimd2xmmKernels)
    {
        kernelSetup.kernelType = NbnxmKernelType::Cpu4xN_Simd_4xN;
    }
    else if (!sc_haveNbnxmSimd4xmKernels && sc_haveNbnxmSimd2xmmKernels)
    {
        kernelSetup.kernelType = NbnxmKernelType::Cpu4xN_Simd_2xNN;
    }
    else
    {
        GMX_RELEASE_ASSERT(sc_haveNbnxmSimd4xmKernels && sc_haveNbnxmSimd2xmmKernels,
                           "Here both 4xM and 2xMM SIMD kernels should be supported");

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
        kernelSetup.kernelType = NbnxmKernelType::Cpu4xN_Simd_4xN;

        if (!GMX_SIMD_HAVE_FMA && (usingPmeOrEwald(inputrec.coulombtype) || usingLJPme(inputrec.vdwtype)))
        {
            /* We have Ewald kernels without FMA (Intel Sandy/Ivy Bridge).
             * There are enough instructions to make 2x(4+4) efficient.
             */
            kernelSetup.kernelType = NbnxmKernelType::Cpu4xN_Simd_2xNN;
        }

        if (hardwareInfo.haveAmdZen1Cpu)
        {
            /* One 256-bit FMA per cycle makes 2xNN faster */
            kernelSetup.kernelType = NbnxmKernelType::Cpu4xN_Simd_2xNN;
        }
    }

    if (getenv("GMX_NBNXN_SIMD_4XN") != nullptr)
    {
        if (sc_haveNbnxmSimd4xmKernels)
        {
            kernelSetup.kernelType = NbnxmKernelType::Cpu4xN_Simd_4xN;
        }
        else
        {
            gmx_fatal(FARGS,
                      "SIMD 4xN kernels requested, but GROMACS has been compiled without support "
                      "for these kernels");
        }
    }
    if (getenv("GMX_NBNXN_SIMD_2XNN") != nullptr)
    {
        if (sc_haveNbnxmSimd2xmmKernels)
        {
            kernelSetup.kernelType = NbnxmKernelType::Cpu4xN_Simd_2xNN;
        }
        else
        {
            gmx_fatal(FARGS,
                      "SIMD 2x(N+N) kernels requested, but GROMACS has been compiled without "
                      "support for these kernels");
        }
    }

    if (kernelSetup.kernelType == NbnxmKernelType::Cpu4xN_Simd_2xNN
        || kernelSetup.kernelType == NbnxmKernelType::Cpu4xN_Simd_4xN)
    {
        /* Analytical Ewald exclusion correction is only an option in
         * the SIMD kernel.
         * Since table lookup's don't parallelize with SIMD, analytical
         * will probably always be faster for a SIMD width of 8 or more.
         * With FMA analytical is sometimes faster for a width if 4 as well.
         * In single precision, this is faster on Bulldozer.
         * On AMD Zen, tabulated Ewald kernels are faster on all 4 combinations
         * of single or double precision and 128 or 256-bit AVX2.
         */
        MSVC_DIAGNOSTIC_IGNORE(6285) // Always zero because compile time constant
        if (
#if GMX_SIMD
                (GMX_SIMD_REAL_WIDTH >= 8 || (GMX_SIMD_REAL_WIDTH >= 4 && GMX_SIMD_HAVE_FMA && !GMX_DOUBLE)) &&
#endif
                !hardwareInfo.haveAmdZen1Cpu)
        {
            kernelSetup.ewaldExclusionType = EwaldExclusionType::Analytical;
        }
        MSVC_DIAGNOSTIC_RESET
        else { kernelSetup.ewaldExclusionType = EwaldExclusionType::Table; }
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

const char* nbnxmKernelTypeToName(const NbnxmKernelType kernelType)
{
    switch (kernelType)
    {
        case NbnxmKernelType::NotSet: return "not set";
        case NbnxmKernelType::Cpu4x4_PlainC: return "plain-C";
        case NbnxmKernelType::Cpu4xN_Simd_4xN: return "SIMD4xM";
        case NbnxmKernelType::Cpu4xN_Simd_2xNN: return "SIMD2xMM";
        case NbnxmKernelType::Gpu8x8x8: return "GPU";
        case NbnxmKernelType::Cpu8x8x8_PlainC: return "plain-C";

        default: gmx_fatal(FARGS, "Illegal kernel type selected");
    }
};

/*! \brief Returns the most suitable kernel type and Ewald handling */
static NbnxmKernelSetup pick_nbnxn_kernel(const gmx::MDLogger&     mdlog,
                                          gmx_bool                 use_simd_kernels,
                                          const gmx_hw_info_t&     hardwareInfo,
                                          const NonbondedResource& nonbondedResource,
                                          const t_inputrec&        inputrec)
{
    NbnxmKernelSetup kernelSetup;

    if (nonbondedResource == NonbondedResource::EmulateGpu)
    {
        kernelSetup.kernelType         = NbnxmKernelType::Cpu8x8x8_PlainC;
        kernelSetup.ewaldExclusionType = EwaldExclusionType::DecidedByGpuModule;

        GMX_LOG(mdlog.warning).asParagraph().appendText("Emulating a GPU run on the CPU (slow)");
    }
    else if (nonbondedResource == NonbondedResource::Gpu)
    {
        kernelSetup.kernelType         = NbnxmKernelType::Gpu8x8x8;
        kernelSetup.ewaldExclusionType = EwaldExclusionType::DecidedByGpuModule;
    }
    else
    {
        if (use_simd_kernels && nbnxn_simd_supported(mdlog, inputrec))
        {
            kernelSetup = pick_nbnxn_kernel_cpu(inputrec, hardwareInfo);
        }
        else
        {
            kernelSetup.kernelType         = NbnxmKernelType::Cpu4x4_PlainC;
            kernelSetup.ewaldExclusionType = EwaldExclusionType::Analytical;
        }
    }

    GMX_LOG(mdlog.info)
            .asParagraph()
            .appendTextFormatted("Using %s %dx%d nonbonded short-range kernels",
                                 nbnxmKernelTypeToName(kernelSetup.kernelType),
                                 sc_iClusterSize(kernelSetup.kernelType),
                                 sc_jClusterSize(kernelSetup.kernelType));

    if (NbnxmKernelType::Cpu4x4_PlainC == kernelSetup.kernelType
        || NbnxmKernelType::Cpu8x8x8_PlainC == kernelSetup.kernelType)
    {
        GMX_LOG(mdlog.warning)
                .asParagraph()
                .appendTextFormatted(
                        "WARNING: Using the slow %s kernels. This should\n"
                        "not happen during routine usage on common platforms.",
                        nbnxmKernelTypeToName(kernelSetup.kernelType));
    }

    GMX_RELEASE_ASSERT(kernelSetup.kernelType != NbnxmKernelType::NotSet
                               && kernelSetup.ewaldExclusionType != EwaldExclusionType::NotSet,
                       "All kernel setup parameters should be set here");

    return kernelSetup;
}

PairlistSets::PairlistSets(const PairlistParams& pairlistParams,
                           const bool            haveMultipleDomains,
                           const int             minimumIlistCountForGpuBalancing) :
    params_(pairlistParams), minimumIlistCountForGpuBalancing_(minimumIlistCountForGpuBalancing)
{
    localSet_ = std::make_unique<PairlistSet>(params_);

    if (haveMultipleDomains)
    {
        nonlocalSet_ = std::make_unique<PairlistSet>(params_);
    }
}

/*! \brief Gets and returns the minimum i-list count for balacing based on the GPU used or env.var. when set */
static int getMinimumIlistCountForGpuBalancing(NbnxmGpu* nbnxmGpu)
{
    if (const char* env = getenv("GMX_NB_MIN_CI"))
    {
        char* end = nullptr;

        int minimumIlistCount = strtol(env, &end, 10);
        if (!end || (*end != 0) || minimumIlistCount < 0)
        {
            gmx_fatal(
                    FARGS, "Invalid value passed in GMX_NB_MIN_CI=%s, non-negative integer required", env);
        }

        if (debug)
        {
            fprintf(debug, "Neighbor-list balancing parameter: %d (passed as env. var.)\n", minimumIlistCount);
        }
        return minimumIlistCount;
    }
    else
    {
        int minimumIlistCount = gpu_min_ci_balanced(nbnxmGpu);
        if (debug)
        {
            fprintf(debug,
                    "Neighbor-list balancing parameter: %d (auto-adjusted to the number of GPU "
                    "multi-processors)\n",
                    minimumIlistCount);
        }
        return minimumIlistCount;
    }
}

//! Returns the LJ combination rule choices for the LJ pair parameters
static std::optional<LJCombinationRule> chooseLJCombinationRule(const t_forcerec& forcerec)
{
    if (forcerec.ic->vdwtype == VanDerWaalsType::Cut
        && (forcerec.ic->vdw_modifier == InteractionModifiers::None
            || forcerec.ic->vdw_modifier == InteractionModifiers::PotShift)
        && getenv("GMX_NO_LJ_COMB_RULE") == nullptr)
    {
        /* Plain LJ cut-off: we can optimize with combination rules */
        return std::nullopt;
    }
    else if (forcerec.ic->vdwtype == VanDerWaalsType::Pme)
    { // NOLINT bugprone-branch-clone
        /* With LJ-PME the NBNxM module does not support combination rules for the pair parameters */
        return LJCombinationRule::None;
    }
    else
    {
        /* We use a full combination matrix: no rule required */
        return LJCombinationRule::None;
    }
}

//! Returns the LJ combination rule choices for the LJ PME-grid parameters
static LJCombinationRule chooseLJPmeCombinationRule(const t_forcerec& forcerec)
{
    if (forcerec.ic->vdwtype == VanDerWaalsType::Pme)
    {
        /* LJ-PME: we need to use a combination rule for the grid and none for the pairs */
        switch (forcerec.ljpme_combination_rule)
        {
            case LongRangeVdW::Geom: return LJCombinationRule::Geometric;
            case LongRangeVdW::LB: return LJCombinationRule::LorentzBerthelot;
            default: GMX_RELEASE_ASSERT(false, "Unhandled case");
        }
    }

    return LJCombinationRule::None;
}

std::unique_ptr<nonbonded_verlet_t> init_nb_verlet(const gmx::MDLogger& mdlog,
                                                   const t_inputrec&    inputrec,
                                                   const t_forcerec&    forcerec,
                                                   const t_commrec*     commrec,
                                                   const gmx_hw_info_t& hardwareInfo,
                                                   bool                 useGpuForNonbonded,
                                                   const gmx::DeviceStreamManager* deviceStreamManager,
                                                   const gmx_mtop_t&               mtop,
                                                   gmx::ObservablesReducerBuilder* observablesReducerBuilder,
                                                   gmx::ArrayRef<const gmx::RVec> coordinates,
                                                   matrix                         box,
                                                   gmx_wallcycle*                 wcycle)
{
    const bool emulateGpu = (getenv("GMX_EMULATE_GPU") != nullptr);

    GMX_RELEASE_ASSERT(!(emulateGpu && useGpuForNonbonded),
                       "When GPU emulation is active, there cannot be a GPU assignment");

    NonbondedResource nonbondedResource;
    if (useGpuForNonbonded)
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

    NbnxmKernelSetup kernelSetup = pick_nbnxn_kernel(
            mdlog, forcerec.use_simd_kernels, hardwareInfo, nonbondedResource, inputrec);

    const bool haveMultipleDomains = havePPDomainDecomposition(commrec);

    bool bFEP_NonBonded = (forcerec.efep != FreeEnergyPerturbationType::No)
                          && haveFepPerturbedNBInteractions(mtop);
    PairlistParams pairlistParams(
            kernelSetup.kernelType, bFEP_NonBonded, inputrec.rlist, haveMultipleDomains);

    const real effectiveAtomDensity = computeEffectiveAtomDensity(
            coordinates, box, std::max(inputrec.rcoulomb, inputrec.rvdw), commrec->mpi_comm_mygroup);

    setupDynamicPairlistPruning(mdlog, inputrec, mtop, effectiveAtomDensity, *forcerec.ic, &pairlistParams);

    if (EI_DYNAMICS(inputrec.eI))
    {
        printNbnxmPressureError(mdlog, inputrec, mtop, effectiveAtomDensity, pairlistParams);
    }

    auto pinPolicy = (useGpuForNonbonded ? gmx::PinningPolicy::PinnedIfSupported
                                         : gmx::PinningPolicy::CannotBePinned);

    int mimimumNumEnergyGroupNonbonded = inputrec.opts.ngener;
    if (inputrec.opts.ngener - inputrec.nwall == 1)
    {
        /* We have only one non-wall energy group, we do not need energy group
         * support in the non-bondeds kernels, since all non-bonded energy
         * contributions go to the first element of the energy group matrix.
         */
        mimimumNumEnergyGroupNonbonded = 1;
    }

    auto nbat = std::make_unique<nbnxn_atomdata_t>(
            pinPolicy,
            mdlog,
            kernelSetup.kernelType,
            chooseLJCombinationRule(forcerec),
            chooseLJPmeCombinationRule(forcerec),
            forcerec.nbfp,
            true,
            mimimumNumEnergyGroupNonbonded,
            (useGpuForNonbonded || emulateGpu) ? 1 : gmx_omp_nthreads_get(ModuleMultiThread::Nonbonded));

    if (forcerec.ic->vdwtype == VanDerWaalsType::Pme)
    {
        GMX_RELEASE_ASSERT(
                (forcerec.ljpme_combination_rule == LongRangeVdW::Geom
                 && nbat->params().ljCombinationRule == LJCombinationRule::Geometric)
                        || (forcerec.ljpme_combination_rule == LongRangeVdW::LB
                            && nbat->params().ljCombinationRule == LJCombinationRule::LorentzBerthelot),
                "nbat combination rule parameters should match those for LJ-PME");
    }

    NbnxmGpu* gpu_nbv                          = nullptr;
    int       minimumIlistCountForGpuBalancing = 0;
    if (useGpuForNonbonded)
    {
        /* init the NxN GPU data; the last argument tells whether we'll have
         * both local and non-local NB calculation on GPU */
        GMX_RELEASE_ASSERT(
                (deviceStreamManager != nullptr),
                "Device stream manager should be initialized in order to use GPU for non-bonded.");
        gpu_nbv = gpu_init(
                *deviceStreamManager, forcerec.ic.get(), pairlistParams, nbat.get(), haveMultipleDomains);

        minimumIlistCountForGpuBalancing = getMinimumIlistCountForGpuBalancing(gpu_nbv);
    }

    auto pairlistSets = std::make_unique<PairlistSets>(
            pairlistParams, haveMultipleDomains, minimumIlistCountForGpuBalancing);

    auto pairSearch = std::make_unique<PairSearch>(
            inputrec.pbcType,
            EI_TPI(inputrec.eI),
            haveDDAtomOrdering(*commrec) ? &commrec->dd->numCells : nullptr,
            haveDDAtomOrdering(*commrec) ? &getDomdecZones(*commrec->dd) : nullptr,
            pairlistParams.pairlistType,
            bFEP_NonBonded,
            gmx_omp_nthreads_get(ModuleMultiThread::Pairsearch),
            pinPolicy);

    std::unique_ptr<ExclusionChecker> exclusionChecker;
    if (inputrec.efep != FreeEnergyPerturbationType::No
        && (usingPmeOrEwald(inputrec.coulombtype) || usingLJPme(inputrec.vdwtype)))
    {
        exclusionChecker = std::make_unique<ExclusionChecker>(commrec, mtop, observablesReducerBuilder);
    }

    return std::make_unique<nonbonded_verlet_t>(std::move(pairlistSets),
                                                std::move(pairSearch),
                                                std::move(nbat),
                                                kernelSetup,
                                                std::move(exclusionChecker),
                                                gpu_nbv,
                                                wcycle);
}

nonbonded_verlet_t::nonbonded_verlet_t(std::unique_ptr<PairlistSets>     pairlistSets,
                                       std::unique_ptr<PairSearch>       pairSearch,
                                       std::unique_ptr<nbnxn_atomdata_t> nbat_in,
                                       const NbnxmKernelSetup&           kernelSetup,
                                       std::unique_ptr<ExclusionChecker> exclusionChecker,
                                       NbnxmGpu*                         gpu_nbv_ptr,
                                       gmx_wallcycle*                    wcycle) :
    pairlistSets_(std::move(pairlistSets)),
    pairSearch_(std::move(pairSearch)),
    nbat_(std::move(nbat_in)),
    kernelSetup_(kernelSetup),
    exclusionChecker_(std::move(exclusionChecker)),
    wcycle_(wcycle),
    gpuNbv_(gpu_nbv_ptr)
{
    GMX_RELEASE_ASSERT(pairlistSets_, "Need valid pairlistSets");
    GMX_RELEASE_ASSERT(pairSearch_, "Need valid search object");
    GMX_RELEASE_ASSERT(nbat_, "Need valid atomdata object");

    if (pairlistSets_->params().haveFep_)
    {
        freeEnergyDispatch_ = std::make_unique<FreeEnergyDispatch>(nbat_->params().numEnergyGroups);
    }
}

nonbonded_verlet_t::nonbonded_verlet_t(std::unique_ptr<PairlistSets>     pairlistSets,
                                       std::unique_ptr<PairSearch>       pairSearch,
                                       std::unique_ptr<nbnxn_atomdata_t> nbat_in,
                                       const NbnxmKernelSetup&           kernelSetup,
                                       NbnxmGpu*                         gpu_nbv_ptr) :
    pairlistSets_(std::move(pairlistSets)),
    pairSearch_(std::move(pairSearch)),
    nbat_(std::move(nbat_in)),
    kernelSetup_(kernelSetup),
    exclusionChecker_(),
    wcycle_(nullptr),
    gpuNbv_(gpu_nbv_ptr)
{
    GMX_RELEASE_ASSERT(pairlistSets_, "Need valid pairlistSets");
    GMX_RELEASE_ASSERT(pairSearch_, "Need valid search object");
    GMX_RELEASE_ASSERT(nbat_, "Need valid atomdata object");

    if (pairlistSets_->params().haveFep_)
    {
        freeEnergyDispatch_ = std::make_unique<FreeEnergyDispatch>(nbat_->params().numEnergyGroups);
    }
}

nonbonded_verlet_t::~nonbonded_verlet_t()
{
    gpu_free(gpuNbv_);
}

} // namespace gmx
