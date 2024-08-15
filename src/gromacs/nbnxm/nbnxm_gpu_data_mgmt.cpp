/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2012- The GROMACS Authors
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
 *  \brief Define common implementation of nbnxm_gpu_data_mgmt.h
 *
 *  \author Anca Hamuraru <anca@streamcomputing.eu>
 *  \author Dimitrios Karkoulis <dimitris.karkoulis@gmail.com>
 *  \author Teemu Virolainen <teemu@streamcomputing.eu>
 *  \author Szilárd Páll <pall.szilard@gmail.com>
 *  \author Artem Zhmurov <zhmurov@gmail.com>
 *
 *  \ingroup module_nbnxm
 */
#include "gmxpre.h"

#include "config.h"

#if GMX_GPU_CUDA
#    include "cuda/nbnxm_cuda_types.h"
#endif

#if GMX_GPU_OPENCL
#    include "opencl/nbnxm_ocl_types.h"
#endif

#if GMX_GPU_SYCL
#    include "sycl/nbnxm_sycl_types.h"
#endif

#include "gromacs/gpu_utils/device_stream_manager.h"
#include "gromacs/gpu_utils/gputraits.h"
#include "gromacs/gpu_utils/hostallocator.h"
#include "gromacs/hardware/device_information.h"
#include "gromacs/mdtypes/interaction_const.h"
#include "gromacs/mdtypes/simulation_workload.h"
#include "gromacs/nbnxm/gpu_common_utils.h"
#include "gromacs/nbnxm/gpu_data_mgmt.h"
#include "gromacs/nbnxm/gridset.h"
#include "gromacs/pbcutil/ishift.h"
#include "gromacs/timing/gpu_timing.h"
#include "gromacs/utility/cstringutil.h"
#include "gromacs/utility/exceptions.h"
#include "gromacs/utility/fatalerror.h"

#include "nbnxm_gpu.h"
#include "nbnxm_gpu_data_mgmt.h"
#include "pairlistsets.h"

namespace gmx
{

static inline void init_ewald_coulomb_force_table(const EwaldCorrectionTables& tables,
                                                  NBParamGpu*                  nbp,
                                                  const DeviceContext&         deviceContext)
{
    if (nbp->coulomb_tab)
    {
        destroyParamLookupTable(&nbp->coulomb_tab, &nbp->coulomb_tab_texobj);
    }

    nbp->coulomb_tab_scale = tables.scale;
    initParamLookupTable(
            &nbp->coulomb_tab, &nbp->coulomb_tab_texobj, tables.tableF.data(), tables.tableF.size(), deviceContext);
}

static bool useTabulatedEwaldByDefault(const DeviceInformation& deviceInfo)
{
    /* By default, use analytical Ewald except:
     *  - NVIDIA CC 7.0 and 8.0,
     *  - all AMD GPUs (although tested for gfx906 and 908 only).
     *
     * Note 1: this function does not handle OpenCL.
     * Note 2: for SYCL, the heuristics are taken from CUDA/HIP ports, and were only partially
     *         verified on oneAPI/hipSYCL themselves on AMD gfx 906/908/90a.
     */
#if GMX_GPU_CUDA
    return (deviceInfo.prop.major == 7 && deviceInfo.prop.minor == 0)
           || (deviceInfo.prop.major == 8 && deviceInfo.prop.minor == 0);
#elif GMX_GPU_SYCL
    switch (deviceInfo.deviceVendor)
    {
        case DeviceVendor::Amd: return true;
        case DeviceVendor::Nvidia:
        {
            const int major = deviceInfo.hardwareVersionMajor.value_or(-1);
            const int minor = deviceInfo.hardwareVersionMinor.value_or(-1);
            return ((major == 7 && minor == 0) || (major == 8 && minor == 0));
        }
        default: return false;
    }
#else
    GMX_UNUSED_VALUE(deviceInfo);
    return false;
#endif
}

static inline ElecType nbnxn_gpu_pick_ewald_kernel_type(const interaction_const_t& ic,
                                                        const DeviceInformation&   deviceInfo)
{
    bool bTwinCut = (ic.rcoulomb != ic.rvdw);

    /* Benchmarking/development environment variables to force the use of
       analytical or tabulated Ewald kernel. */
    const bool forceAnalyticalEwald = (getenv("GMX_GPU_NB_ANA_EWALD") != nullptr);
    const bool forceTabulatedEwald  = (getenv("GMX_GPU_NB_TAB_EWALD") != nullptr);
    const bool forceTwinCutoffEwald = (getenv("GMX_GPU_NB_EWALD_TWINCUT") != nullptr);

    if (forceAnalyticalEwald && forceTabulatedEwald)
    {
        gmx_incons(
                "Both analytical and tabulated Ewald GPU non-bonded kernels "
                "requested through environment variables.");
    }

    const bool c_useTabulatedEwaldDefault = useTabulatedEwaldByDefault(deviceInfo);
    bool       bUseAnalyticalEwald        = !c_useTabulatedEwaldDefault;
    if (forceAnalyticalEwald)
    {
        bUseAnalyticalEwald = true;
        if (debug)
        {
            fprintf(debug, "Using analytical Ewald GPU kernels\n");
        }
    }
    else if (forceTabulatedEwald)
    {
        bUseAnalyticalEwald = false;

        if (debug)
        {
            fprintf(debug, "Using tabulated Ewald GPU kernels\n");
        }
    }

    /* Use twin cut-off kernels if requested by bTwinCut or the env. var.
       forces it (use it for debugging/benchmarking only). */
    if (!bTwinCut && !forceTwinCutoffEwald)
    {
        return bUseAnalyticalEwald ? ElecType::EwaldAna : ElecType::EwaldTab;
    }
    else
    {
        return bUseAnalyticalEwald ? ElecType::EwaldAnaTwin : ElecType::EwaldTabTwin;
    }
}

static inline void set_cutoff_parameters(NBParamGpu*                nbp,
                                         const interaction_const_t& ic,
                                         const PairlistParams&      listParams)
{
    nbp->ewald_beta        = ic.ewaldcoeff_q;
    nbp->sh_ewald          = ic.sh_ewald;
    nbp->epsfac            = ic.epsfac;
    nbp->two_k_rf          = 2.0 * ic.reactionFieldCoefficient;
    nbp->c_rf              = ic.reactionFieldShift;
    nbp->rvdw_sq           = ic.rvdw * ic.rvdw;
    nbp->rcoulomb_sq       = ic.rcoulomb * ic.rcoulomb;
    nbp->rlistOuter_sq     = listParams.rlistOuter * listParams.rlistOuter;
    nbp->rlistInner_sq     = listParams.rlistInner * listParams.rlistInner;
    nbp->useDynamicPruning = listParams.useDynamicPruning;

    nbp->sh_lj_ewald   = ic.sh_lj_ewald;
    nbp->ewaldcoeff_lj = ic.ewaldcoeff_lj;

    nbp->rvdw_switch      = ic.rvdw_switch;
    nbp->dispersion_shift = ic.dispersion_shift;
    nbp->repulsion_shift  = ic.repulsion_shift;
    nbp->vdw_switch       = ic.vdw_switch;
}

GpuPairlistSorting::GpuPairlistSorting() {}

GpuPairlistSorting::~GpuPairlistSorting()
{
    if (nbnxmSortListsOnGpu())
    {
        freeDeviceBuffer(&scanTemporary);
        freeDeviceBuffer(&sciHistogram);
        freeDeviceBuffer(&sciOffset);
        freeDeviceBuffer(&sciCount);
        freeDeviceBuffer(&sciSorted);
    }
}


GpuPairlist::GpuPairlist() {}

GpuPairlist::~GpuPairlist()
{
    freeDeviceBuffer(&sci);
    freeDeviceBuffer(&cjPacked);
    freeDeviceBuffer(&imask);
    freeDeviceBuffer(&excl);
    freeDeviceBuffer(&d_rollingPruningPart);
}

static inline void init_timings(gmx_wallclock_gpu_nbnxn_t* t)
{
    t->nb_h2d_t = 0.0;
    t->nb_d2h_t = 0.0;
    t->nb_c     = 0;
    t->pl_h2d_t = 0.0;
    t->pl_h2d_c = 0;
    for (int i = 0; i < 2; i++)
    {
        for (int j = 0; j < 2; j++)
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

/*! \brief Initialize \p atomdata first time; it only gets filled at pair-search. */
static inline void initAtomdataFirst(NBAtomDataGpu*       atomdata,
                                     int                  numTypes,
                                     const DeviceContext& deviceContext,
                                     const DeviceStream&  localStream)
{
    atomdata->numTypes = numTypes;
    allocateDeviceBuffer(&atomdata->shiftVec, c_numShiftVectors, deviceContext);
    atomdata->shiftVecUploaded = false;

    allocateDeviceBuffer(&atomdata->fShift, c_numShiftVectors, deviceContext);
    allocateDeviceBuffer(&atomdata->eLJ, 1, deviceContext);
    allocateDeviceBuffer(&atomdata->eElec, 1, deviceContext);

    clearDeviceBufferAsync(&atomdata->fShift, 0, c_numShiftVectors, localStream);
    clearDeviceBufferAsync(&atomdata->eElec, 0, 1, localStream);
    clearDeviceBufferAsync(&atomdata->eLJ, 0, 1, localStream);

    /* initialize to nullptr pointers to data that is not allocated here and will
       need reallocation in later */
    atomdata->xq = nullptr;
    atomdata->f  = nullptr;

    /* size -1 indicates that the respective array hasn't been initialized yet */
    atomdata->numAtoms      = -1;
    atomdata->numAtomsAlloc = -1;
}

static inline VdwType nbnxmGpuPickVdwKernelType(const interaction_const_t& ic,
                                                LJCombinationRule          ljCombinationRule)
{
    if (ic.vdwtype == VanDerWaalsType::Cut)
    {
        switch (ic.vdw_modifier)
        {
            case InteractionModifiers::None:
            case InteractionModifiers::PotShift:
                switch (ljCombinationRule)
                {
                    case LJCombinationRule::None: return VdwType::Cut;
                    case LJCombinationRule::Geometric: return VdwType::CutCombGeom;
                    case LJCombinationRule::LorentzBerthelot: return VdwType::CutCombLB;
                    default:
                        GMX_THROW(InconsistentInputError(formatString(
                                "The requested LJ combination rule %s is not implemented in "
                                "the GPU accelerated kernels!",
                                enumValueToString(ljCombinationRule))));
                }
            case InteractionModifiers::ForceSwitch: return VdwType::FSwitch;
            case InteractionModifiers::PotSwitch: return VdwType::PSwitch;
            default:
                GMX_THROW(InconsistentInputError(
                        formatString("The requested VdW interaction modifier %s is not "
                                     "implemented in the GPU accelerated kernels!",
                                     enumValueToString(ic.vdw_modifier))));
        }
    }
    else if (ic.vdwtype == VanDerWaalsType::Pme)
    {
        if (ic.ljpme_comb_rule == LongRangeVdW::Geom)
        {
            GMX_RELEASE_ASSERT(
                    ljCombinationRule == LJCombinationRule::Geometric,
                    "Combination rules for long- and short-range interactions should match.");
            return VdwType::EwaldGeom;
        }
        else
        {
            GMX_RELEASE_ASSERT(
                    ljCombinationRule == LJCombinationRule::LorentzBerthelot,
                    "Combination rules for long- and short-range interactions should match.");
            return VdwType::EwaldLB;
        }
    }
    else
    {
        GMX_THROW(InconsistentInputError(formatString(
                "The requested VdW type %s is not implemented in the GPU accelerated kernels!",
                enumValueToString(ic.vdwtype))));
    }
}

static inline ElecType nbnxmGpuPickElectrostaticsKernelType(const interaction_const_t& ic,
                                                            const DeviceInformation&   deviceInfo)
{
    if (ic.eeltype == CoulombInteractionType::Cut)
    {
        return ElecType::Cut;
    }
    else if (usingRF(ic.eeltype))
    {
        return ElecType::RF;
    }
    else if ((usingPme(ic.eeltype) || ic.eeltype == CoulombInteractionType::Ewald))
    {
        return nbnxn_gpu_pick_ewald_kernel_type(ic, deviceInfo);
    }
    else
    {
        /* Shouldn't happen, as this is checked when choosing Verlet-scheme */
        GMX_THROW(InconsistentInputError(
                formatString("The requested electrostatics type %s is not implemented in "
                             "the GPU accelerated kernels!",
                             enumValueToString(ic.eeltype))));
    }
}

/*! \brief Initialize the nonbonded parameter data structure. */
static inline void initNbparam(NBParamGpu*                     nbp,
                               const interaction_const_t&      ic,
                               const PairlistParams&           listParams,
                               const nbnxn_atomdata_t::Params& nbatParams,
                               const DeviceContext&            deviceContext)
{
    const int numTypes = nbatParams.numTypes;

    set_cutoff_parameters(nbp, ic, listParams);

    nbp->vdwType  = nbnxmGpuPickVdwKernelType(ic, nbatParams.ljCombinationRule);
    nbp->elecType = nbnxmGpuPickElectrostaticsKernelType(ic, deviceContext.deviceInfo());

    if (ic.vdwtype == VanDerWaalsType::Pme)
    {
        GMX_ASSERT((ic.ljpme_comb_rule == LongRangeVdW::Geom
                    && nbatParams.ljCombinationRule == LJCombinationRule::Geometric)
                           || (ic.ljpme_comb_rule == LongRangeVdW::LB
                               && nbatParams.ljCombinationRule == LJCombinationRule::LorentzBerthelot),
                   "Combination rule mismatch!");
    }

    /* generate table for PME */
    if (nbp->elecType == ElecType::EwaldTab || nbp->elecType == ElecType::EwaldTabTwin)
    {
        GMX_RELEASE_ASSERT(ic.coulombEwaldTables, "Need valid Coulomb Ewald correction tables");
        init_ewald_coulomb_force_table(*ic.coulombEwaldTables, nbp, deviceContext);
    }

    /* set up LJ parameter lookup table */
    if (!useLjCombRule(nbp->vdwType))
    {
        static_assert(sizeof(decltype(nbp->nbfp)) == 2 * sizeof(decltype(*nbatParams.nbfp.data())),
                      "Mismatch in the size of host / device data types");
        initParamLookupTable(&nbp->nbfp,
                             &nbp->nbfp_texobj,
                             reinterpret_cast<const Float2*>(nbatParams.nbfp.data()),
                             numTypes * numTypes,
                             deviceContext);
    }

    /* set up LJ-PME parameter lookup table */
    if (ic.vdwtype == VanDerWaalsType::Pme)
    {
        static_assert(sizeof(decltype(nbp->nbfp_comb))
                              == 2 * sizeof(decltype(*nbatParams.nbfp_comb.data())),
                      "Mismatch in the size of host / device data types");
        initParamLookupTable(&nbp->nbfp_comb,
                             &nbp->nbfp_comb_texobj,
                             reinterpret_cast<const Float2*>(nbatParams.nbfp_comb.data()),
                             numTypes,
                             deviceContext);
    }
}

using GpuPairlistByLocality = EnumerationArray<InteractionLocality, std::unique_ptr<GpuPairlist>>;

static GpuPairlistByLocality initializeGpuLists(bool localAndNonLocal)
{
    GpuPairlistByLocality list;
    list[InteractionLocality::Local] = std::make_unique<GpuPairlist>();
    if (localAndNonLocal)
    {
        list[InteractionLocality::NonLocal] = std::make_unique<GpuPairlist>();
    }
    return list;
}

NbnxmGpu* gpu_init(const DeviceStreamManager& deviceStreamManager,
                   const interaction_const_t* ic,
                   const PairlistParams&      listParams,
                   const nbnxn_atomdata_t*    nbat,
                   const bool                 bLocalAndNonlocal)
{
    auto* nb           = new NbnxmGpu();
    nb->deviceContext_ = &deviceStreamManager.context();
    nb->atdat          = new NBAtomDataGpu;
    nb->nbparam        = new NBParamGpu;

    nb->plist = initializeGpuLists(bLocalAndNonlocal);

    nb->bUseTwoStreams = bLocalAndNonlocal;

    nb->timers = new GpuTimers();
    snew(nb->timings, 1);

    nb->bDoTime = decideGpuTimingsUsage();

    if (nb->bDoTime)
    {
        init_timings(nb->timings);
    }

    /* init nbst */
    changePinningPolicy(&nb->nbst.eLJ, PinningPolicy::PinnedIfSupported);
    changePinningPolicy(&nb->nbst.eElec, PinningPolicy::PinnedIfSupported);
    changePinningPolicy(&nb->nbst.fShift, PinningPolicy::PinnedIfSupported);

    nb->nbst.eLJ.resize(1);
    nb->nbst.eElec.resize(1);
    nb->nbst.fShift.resize(c_numShiftVectors);

    /* local/non-local GPU streams */
    GMX_RELEASE_ASSERT(deviceStreamManager.streamIsValid(DeviceStreamType::NonBondedLocal),
                       "Local non-bonded stream should be initialized to use GPU for non-bonded.");
    const DeviceStream& localStream = deviceStreamManager.stream(DeviceStreamType::NonBondedLocal);
    nb->deviceStreams[InteractionLocality::Local] = &localStream;
    // In general, it's not strictly necessary to use 2 streams for SYCL, since they are
    // out-of-order. But for the time being, it will be less disruptive to keep them.
    if (nb->bUseTwoStreams)
    {
        GMX_RELEASE_ASSERT(deviceStreamManager.streamIsValid(DeviceStreamType::NonBondedNonLocal),
                           "Non-local non-bonded stream should be initialized to use GPU for "
                           "non-bonded with domain decomposition.");
        nb->deviceStreams[InteractionLocality::NonLocal] =
                &deviceStreamManager.stream(DeviceStreamType::NonBondedNonLocal);
    }

    const nbnxn_atomdata_t::Params& nbatParams    = nbat->params();
    const DeviceContext&            deviceContext = *nb->deviceContext_;

    initNbparam(nb->nbparam, *ic, listParams, nbatParams, deviceContext);
    initAtomdataFirst(nb->atdat, nbatParams.numTypes, deviceContext, localStream);

    gpu_init_platform_specific(nb);

    if (debug)
    {
        fprintf(debug, "Initialized NBNXM GPU data structures.\n");
    }

    return nb;
}

void gpu_pme_loadbal_update_param(nonbonded_verlet_t* nbv, const interaction_const_t& ic)
{
    if (!nbv || !nbv->useGpu())
    {
        return;
    }
    NbnxmGpu*   nb  = nbv->gpuNbv();
    NBParamGpu* nbp = nb->nbparam;

    set_cutoff_parameters(nbp, ic, nbv->pairlistSets().params());

    nbp->elecType = nbnxn_gpu_pick_ewald_kernel_type(ic, nb->deviceContext_->deviceInfo());

    GMX_RELEASE_ASSERT(ic.coulombEwaldTables, "Need valid Coulomb Ewald correction tables");
    init_ewald_coulomb_force_table(*ic.coulombEwaldTables, nbp, *nb->deviceContext_);
}

void gpu_upload_shiftvec(NbnxmGpu* nb, const nbnxn_atomdata_t* nbatom)
{
    NBAtomDataGpu*      adat        = nb->atdat;
    const DeviceStream& localStream = *nb->deviceStreams[InteractionLocality::Local];

    /* only if we have a dynamic box */
    if (nbatom->bDynamicBox || !adat->shiftVecUploaded)
    {
        copyToDeviceBuffer(&adat->shiftVec,
                           asGenericFloat3Pointer(nbatom->shift_vec),
                           0,
                           c_numShiftVectors,
                           localStream,
                           GpuApiCallBehavior::Async,
                           nullptr);
        adat->shiftVecUploaded = true;
    }
}

//! This function is documented in the header file
void gpu_init_pairlist(NbnxmGpu* nb, const NbnxnPairlistGpu* h_plist, const InteractionLocality iloc)
{
    char sbuf[STRLEN];
    // Timing accumulation should happen only if there was work to do
    // because getLastRangeTime() gets skipped with empty lists later
    // which leads to the counter not being reset.
    bool                bDoTime      = (nb->bDoTime && !h_plist->sci.empty());
    const DeviceStream& deviceStream = *nb->deviceStreams[iloc];
    auto*               d_plist      = nb->plist[iloc].get();

    if (d_plist->numAtomsPerCluster < 0)
    {
        d_plist->numAtomsPerCluster = h_plist->na_ci;
    }
    else
    {
        if (d_plist->numAtomsPerCluster != h_plist->na_ci)
        {
            sprintf(sbuf,
                    "In init_plist: the #atoms per cell has changed (from %d to %d)",
                    d_plist->numAtomsPerCluster,
                    h_plist->na_ci);
            gmx_incons(sbuf);
        }
    }

    GpuTimers::Interaction& iTimers = nb->timers->interaction[iloc];

    if (bDoTime)
    {
        iTimers.pl_h2d.openTimingRegion(deviceStream);
        iTimers.didPairlistH2D = true;
    }

    // TODO most of this function is same in CUDA and OpenCL, move into the header
    const DeviceContext& deviceContext = *nb->deviceContext_;

    reallocateDeviceBuffer(
            &d_plist->sci, h_plist->sci.size(), &d_plist->numSci, &d_plist->sciAllocationSize, deviceContext);
    copyToDeviceBuffer(&d_plist->sci,
                       h_plist->sci.data(),
                       0,
                       h_plist->sci.size(),
                       deviceStream,
                       GpuApiCallBehavior::Async,
                       bDoTime ? iTimers.pl_h2d.fetchNextEvent() : nullptr);

    /* gpu sorting is only implemented for cuda */
    if (nbnxmSortListsOnGpu())
    {
        reallocateDeviceBuffer(&d_plist->sorting.sciSorted,
                               h_plist->sci.size(),
                               &d_plist->sorting.nsciSorted,
                               &d_plist->sorting.sciSortedNalloc,
                               deviceContext);
        copyToDeviceBuffer(&d_plist->sorting.sciSorted,
                           h_plist->sci.data(),
                           0,
                           h_plist->sci.size(),
                           deviceStream,
                           GpuApiCallBehavior::Async,
                           bDoTime ? iTimers.pl_h2d.fetchNextEvent() : nullptr);
        reallocateDeviceBuffer(&d_plist->sorting.sciCount,
                               h_plist->sci.size(),
                               &d_plist->sorting.nsciCounted,
                               &d_plist->sorting.sciCountedNalloc,
                               deviceContext);

        if (d_plist->sorting.nscanTemporary == -1)
        {
            reallocateDeviceBuffer(&d_plist->sorting.sciHistogram,
                                   c_sciHistogramSize + 1,
                                   &d_plist->sorting.nsciHistogram,
                                   &d_plist->sorting.sciHistogramNalloc,
                                   deviceContext);

            reallocateDeviceBuffer(&d_plist->sorting.sciOffset,
                                   c_sciHistogramSize,
                                   &d_plist->sorting.nsciOffset,
                                   &d_plist->sorting.sciOffsetNalloc,
                                   deviceContext);

            size_t scanTemporarySize = getExclusiveScanWorkingArraySize(d_plist, deviceStream);

            reallocateDeviceBuffer(&d_plist->sorting.scanTemporary,
                                   scanTemporarySize,
                                   &d_plist->sorting.nscanTemporary,
                                   &d_plist->sorting.scanTemporaryNalloc,
                                   deviceContext);
        }

        clearDeviceBufferAsync(&d_plist->sorting.sciHistogram, 0, c_sciHistogramSize, deviceStream);
    }

    reallocateDeviceBuffer(&d_plist->cjPacked,
                           h_plist->cjPacked.size(),
                           &d_plist->numPackedJClusters,
                           &d_plist->packedJClustersAllocationSize,
                           deviceContext);
    copyToDeviceBuffer(&d_plist->cjPacked,
                       h_plist->cjPacked.list_.data(),
                       0,
                       h_plist->cjPacked.size(),
                       deviceStream,
                       GpuApiCallBehavior::Async,
                       bDoTime ? iTimers.pl_h2d.fetchNextEvent() : nullptr);

    reallocateDeviceBuffer(&d_plist->imask,
                           h_plist->cjPacked.size() * c_nbnxnGpuClusterpairSplit,
                           &d_plist->numIMask,
                           &d_plist->iMaskAllocationSize,
                           deviceContext);

    reallocateDeviceBuffer(
            &d_plist->excl, h_plist->excl.size(), &d_plist->numExcl, &d_plist->exclAllocationSize, deviceContext);
    copyToDeviceBuffer(&d_plist->excl,
                       h_plist->excl.data(),
                       0,
                       h_plist->excl.size(),
                       deviceStream,
                       GpuApiCallBehavior::Async,
                       bDoTime ? iTimers.pl_h2d.fetchNextEvent() : nullptr);

    // Part index is managed on device, through a buffer containing a separate copy of the index
    // per block, to allow asynchronous incrementation on the GPU without CPU involvement.
    reallocateDeviceBuffer(&d_plist->d_rollingPruningPart,
                           d_plist->numSci,
                           &d_plist->d_rollingPruningPartSize,
                           &d_plist->d_rollingPruningPartAllocationSize,
                           *nb->deviceContext_);
    clearDeviceBufferAsync(&d_plist->d_rollingPruningPart, 0, d_plist->numSci, deviceStream);

    if (bDoTime)
    {
        iTimers.pl_h2d.closeTimingRegion(deviceStream);
    }

    /* need to prune the pair list during the next step */
    d_plist->haveFreshList = true;
}

void gpu_init_atomdata(NbnxmGpu* nb, const nbnxn_atomdata_t* nbat)
{
    bool                 bDoTime       = nb->bDoTime;
    GpuTimers*           timers        = bDoTime ? nb->timers : nullptr;
    NBAtomDataGpu*       atdat         = nb->atdat;
    const DeviceContext& deviceContext = *nb->deviceContext_;
    const DeviceStream&  localStream   = *nb->deviceStreams[InteractionLocality::Local];

    int  numAtoms  = nbat->numAtoms();
    bool realloced = false;

    if (bDoTime)
    {
        /* time async copy */
        timers->atdat.openTimingRegion(localStream);
    }

    /* need to reallocate if we have to copy more atoms than the amount of space
       available and only allocate if we haven't initialized yet, i.e atdat->natoms == -1 */
    if (numAtoms > atdat->numAtomsAlloc)
    {
        int numAlloc = over_alloc_small(numAtoms);

        /* free up first if the arrays have already been initialized */
        if (atdat->numAtomsAlloc != -1)
        {
            // Wait for the force-buffer clearing from the previous step
            localStream.synchronize();

            freeDeviceBuffer(&atdat->f);
            freeDeviceBuffer(&atdat->xq);
            if (useLjCombRule(nb->nbparam->vdwType))
            {
                freeDeviceBuffer(&atdat->ljComb);
            }
            else
            {
                freeDeviceBuffer(&atdat->atomTypes);
            }
        }


        allocateDeviceBuffer(&atdat->f, numAlloc, deviceContext);
        allocateDeviceBuffer(&atdat->xq, numAlloc, deviceContext);

        if (useLjCombRule(nb->nbparam->vdwType))
        {
            // Two Lennard-Jones parameters per atom
            allocateDeviceBuffer(&atdat->ljComb, numAlloc, deviceContext);
        }
        else
        {
            allocateDeviceBuffer(&atdat->atomTypes, numAlloc, deviceContext);
        }

        atdat->numAtomsAlloc = numAlloc;
        realloced            = true;
    }

    atdat->numAtoms      = numAtoms;
    atdat->numAtomsLocal = nbat->numLocalAtoms();

    /* need to clear GPU f output if realloc happened */
    if (realloced)
    {
        clearDeviceBufferAsync(&atdat->f, 0, atdat->numAtomsAlloc, localStream);
    }

    if (useLjCombRule(nb->nbparam->vdwType))
    {
        static_assert(
                sizeof(Float2) == 2 * sizeof(*nbat->params().lj_comb.data()),
                "Size of a pair of LJ parameters elements should be equal to the size of Float2.");
        copyToDeviceBuffer(&atdat->ljComb,
                           reinterpret_cast<const Float2*>(nbat->params().lj_comb.data()),
                           0,
                           numAtoms,
                           localStream,
                           GpuApiCallBehavior::Async,
                           bDoTime ? timers->atdat.fetchNextEvent() : nullptr);
    }
    else
    {
        static_assert(sizeof(int) == sizeof(*nbat->params().type.data()),
                      "Sizes of host- and device-side atom types should be the same.");
        copyToDeviceBuffer(&atdat->atomTypes,
                           nbat->params().type.data(),
                           0,
                           numAtoms,
                           localStream,
                           GpuApiCallBehavior::Async,
                           bDoTime ? timers->atdat.fetchNextEvent() : nullptr);
    }

    if (bDoTime)
    {
        timers->atdat.closeTimingRegion(localStream);
    }

    /* kick off the tasks enqueued above to ensure concurrency with the search */
    issueClFlushInStream(localStream);
}

void gpu_clear_outputs(NbnxmGpu* nb, bool computeVirial)
{
    NBAtomDataGpu*      adat        = nb->atdat;
    const DeviceStream& localStream = *nb->deviceStreams[InteractionLocality::Local];
    // Clear forces
    clearDeviceBufferAsync(&adat->f, 0, nb->atdat->numAtoms, localStream);
    // Clear shift force array and energies if the outputs were used in the current step
    if (computeVirial)
    {
        clearDeviceBufferAsync(&adat->fShift, 0, c_numShiftVectors, localStream);
        clearDeviceBufferAsync(&adat->eLJ, 0, 1, localStream);
        clearDeviceBufferAsync(&adat->eElec, 0, 1, localStream);
    }
    issueClFlushInStream(localStream);
}

//! This function is documented in the header file
gmx_wallclock_gpu_nbnxn_t* gpu_get_timings(NbnxmGpu* nb)
{
    return (nb != nullptr && nb->bDoTime) ? nb->timings : nullptr;
}

//! This function is documented in the header file
void gpu_reset_timings(nonbonded_verlet_t* nbv)
{
    if (nbv->gpuNbv() && nbv->gpuNbv()->bDoTime)
    {
        init_timings(nbv->gpuNbv()->timings);
    }
}

bool gpu_is_kernel_ewald_analytical(const NbnxmGpu* nb)
{
    return ((nb->nbparam->elecType == ElecType::EwaldAna)
            || (nb->nbparam->elecType == ElecType::EwaldAnaTwin));
}

void setupGpuShortRangeWorkLow(NbnxmGpu*                 nb,
                               const ListedForcesGpu*    listedForcesGpu,
                               const InteractionLocality iLocality)
{
    GMX_ASSERT(nb, "Need a valid nbnxn_gpu object");

    // There is short-range work if the pair list for the provided
    // interaction locality contains entries or if there is any
    // bonded work (as this is not split into local/nonlocal).
    nb->haveWork[iLocality] = ((nb->plist[iLocality]->numSci != 0)
                               || (listedForcesGpu != nullptr && listedForcesGpu->haveInteractions()));
}

bool haveGpuShortRangeWork(const NbnxmGpu* nb, const InteractionLocality interactionLocality)
{
    GMX_ASSERT(nb, "Need a valid nbnxn_gpu object");

    return nb->haveWork[interactionLocality];
}

/*! \brief
 * Launch asynchronously the download of nonbonded forces from the GPU
 * (and energies/shift forces if required).
 */
void gpu_launch_cpyback(NbnxmGpu*                nb,
                        struct nbnxn_atomdata_t* nbatom,
                        const StepWorkload&      stepWork,
                        const AtomLocality       atomLocality)
{
    GMX_ASSERT(nb, "Need a valid nbnxn_gpu object");

    /* determine interaction locality from atom locality */
    const InteractionLocality iloc = atomToInteractionLocality(atomLocality);
    GMX_ASSERT(iloc == InteractionLocality::Local
                       || (iloc == InteractionLocality::NonLocal && nb->bNonLocalStreamDoneMarked == false),
               "Non-local stream is indicating that the copy back event is enqueued at the "
               "beginning of the copy back function.");

    /* extract the data */
    NBAtomDataGpu*      adat         = nb->atdat;
    GpuTimers*          timers       = nb->timers;
    bool                bDoTime      = nb->bDoTime;
    const DeviceStream& deviceStream = *nb->deviceStreams[iloc];

    /* don't launch non-local copy-back if there was no non-local work to do */
    if ((iloc == InteractionLocality::NonLocal) && !haveGpuShortRangeWork(nb, iloc))
    {
        /* TODO An alternative way to signal that non-local work is
           complete is to use a clEnqueueMarker+clEnqueueBarrier
           pair. However, the use of bNonLocalStreamDoneMarked has the
           advantage of being local to the host, so probably minimizes
           overhead. Curiously, for NVIDIA OpenCL with an empty-domain
           test case, overall simulation performance was higher with
           the API calls, but this has not been tested on AMD OpenCL,
           so could be worth considering in future. */
        nb->bNonLocalStreamDoneMarked = false;
        return;
    }

    /* local/nonlocal offset and length used for xq and f */
    auto atomsRange = getGpuAtomRange(adat, atomLocality);

    /* beginning of timed D2H section */
    if (bDoTime)
    {
        timers->xf[atomLocality].nb_d2h.openTimingRegion(deviceStream);
    }

    /* With DD the local D2H transfer can only start after the non-local
       has been launched. */
    if (iloc == InteractionLocality::Local && nb->bNonLocalStreamDoneMarked)
    {
        nb->nonlocal_done.enqueueWaitEvent(deviceStream);
        nb->bNonLocalStreamDoneMarked = false;
    }

    /* DtoH f */
    if (!stepWork.useGpuFBufferOps)
    {
        static_assert(
                sizeof(*nbatom->outputBuffer(0).f.data()) == sizeof(float),
                "The host force buffer should be in single precision to match device data size.");
        copyFromDeviceBuffer(
                reinterpret_cast<Float3*>(nbatom->outputBuffer(0).f.data()) + atomsRange.begin(),
                &adat->f,
                atomsRange.begin(),
                atomsRange.size(),
                deviceStream,
                GpuApiCallBehavior::Async,
                bDoTime ? timers->xf[atomLocality].nb_d2h.fetchNextEvent() : nullptr);

        issueClFlushInStream(deviceStream);
    }

    /* After the non-local D2H is launched the nonlocal_done event can be
       recorded which signals that the local D2H can proceed. This event is not
       placed after the non-local kernel because we first need the non-local
       data back first. */
    if (iloc == InteractionLocality::NonLocal)
    {
        nb->nonlocal_done.markEvent(deviceStream);
        nb->bNonLocalStreamDoneMarked = true;
    }

    /* only transfer energies in the local stream */
    if (iloc == InteractionLocality::Local)
    {
        /* DtoH fshift when virial is needed */
        if (stepWork.computeVirial)
        {
            static_assert(
                    sizeof(*nb->nbst.fShift.data()) == sizeof(Float3),
                    "Sizes of host- and device-side shift vector elements should be the same.");
            copyFromDeviceBuffer(nb->nbst.fShift.data(),
                                 &adat->fShift,
                                 0,
                                 c_numShiftVectors,
                                 deviceStream,
                                 GpuApiCallBehavior::Async,
                                 bDoTime ? timers->xf[atomLocality].nb_d2h.fetchNextEvent() : nullptr);
        }

        /* DtoH energies */
        if (stepWork.computeEnergy)
        {
            static_assert(sizeof(*nb->nbst.eLJ.data()) == sizeof(float),
                          "Sizes of host- and device-side LJ energy terms should be the same.");
            copyFromDeviceBuffer(nb->nbst.eLJ.data(),
                                 &adat->eLJ,
                                 0,
                                 1,
                                 deviceStream,
                                 GpuApiCallBehavior::Async,
                                 bDoTime ? timers->xf[atomLocality].nb_d2h.fetchNextEvent() : nullptr);
            static_assert(sizeof(*nb->nbst.eElec.data()) == sizeof(float),
                          "Sizes of host- and device-side electrostatic energy terms should be the "
                          "same.");
            copyFromDeviceBuffer(nb->nbst.eElec.data(),
                                 &adat->eElec,
                                 0,
                                 1,
                                 deviceStream,
                                 GpuApiCallBehavior::Async,
                                 bDoTime ? timers->xf[atomLocality].nb_d2h.fetchNextEvent() : nullptr);
        }
    }

    if (bDoTime)
    {
        timers->xf[atomLocality].nb_d2h.closeTimingRegion(deviceStream);
    }
}

void nbnxnInsertNonlocalGpuDependency(NbnxmGpu* nb, const InteractionLocality interactionLocality)
{
    const DeviceStream& deviceStream = *nb->deviceStreams[interactionLocality];

    /* When we get here all misc operations issued in the local stream as well as
       the local xq H2D are done,
       so we record that in the local stream and wait for it in the nonlocal one.
       This wait needs to precede any PP tasks, bonded or nonbonded, that may
       compute on interactions between local and nonlocal atoms.
     */
    if (nb->bUseTwoStreams)
    {
        if (interactionLocality == InteractionLocality::Local)
        {
            nb->misc_ops_and_local_H2D_done.markEvent(deviceStream);
            issueClFlushInStream(deviceStream);
        }
        else
        {
            nb->misc_ops_and_local_H2D_done.enqueueWaitEvent(deviceStream);
        }
    }
}

/*! \brief Launch asynchronously the xq buffer host to device copy. */
void gpu_copy_xq_to_gpu(NbnxmGpu* nb, const nbnxn_atomdata_t* nbatom, const AtomLocality atomLocality)
{
    GMX_ASSERT(nb, "Need a valid nbnxn_gpu object");

    const InteractionLocality iloc = atomToInteractionLocality(atomLocality);

    NBAtomDataGpu*      adat         = nb->atdat;
    auto*               plist        = nb->plist[iloc].get();
    GpuTimers*          timers       = nb->timers;
    const DeviceStream& deviceStream = *nb->deviceStreams[iloc];

    const bool bDoTime = nb->bDoTime;

    /* Don't launch the non-local H2D copy if there is no dependent
       work to do: neither non-local nor other (e.g. bonded) work
       to do that has as input the nbnxn coordaintes.
       Doing the same for the local kernel is more complicated, since the
       local part of the force array also depends on the non-local kernel.
       So to avoid complicating the code and to reduce the risk of bugs,
       we always call the local local x+q copy (and the rest of the local
       work in nbnxn_gpu_launch_kernel().
     */
    if ((iloc == InteractionLocality::NonLocal) && !haveGpuShortRangeWork(nb, iloc))
    {
        plist->haveFreshList = false;

        // The event is marked for Local interactions unconditionally,
        // so it has to be released here because of the early return
        // for NonLocal interactions.
        nb->misc_ops_and_local_H2D_done.reset();

        return;
    }

    /* local/nonlocal offset and length used for xq and f */
    const auto atomsRange = getGpuAtomRange(adat, atomLocality);

    /* beginning of timed HtoD section */
    if (bDoTime)
    {
        timers->xf[atomLocality].nb_h2d.openTimingRegion(deviceStream);
    }

    /* HtoD x, q */
    GMX_ASSERT(nbatom->XFormat == nbatXYZQ,
               "The coordinates should be in xyzq format to copy to the Float4 device buffer.");
    copyToDeviceBuffer(&adat->xq,
                       reinterpret_cast<const Float4*>(nbatom->x().data()) + atomsRange.begin(),
                       atomsRange.begin(),
                       atomsRange.size(),
                       deviceStream,
                       GpuApiCallBehavior::Async,
                       nullptr);

    if (bDoTime)
    {
        timers->xf[atomLocality].nb_h2d.closeTimingRegion(deviceStream);
    }

    /* When we get here all misc operations issued in the local stream as well as
       the local xq H2D are done,
       so we record that in the local stream and wait for it in the nonlocal one.
       This wait needs to precede any PP tasks, bonded or nonbonded, that may
       compute on interactions between local and nonlocal atoms.
     */
    nbnxnInsertNonlocalGpuDependency(nb, iloc);
}


/* Initialization for X buffer operations on GPU. */
void nbnxn_gpu_init_x_to_nbat_x(const GridSet& gridSet, NbnxmGpu* gpu_nbv)
{
    const DeviceStream& localStream   = *gpu_nbv->deviceStreams[InteractionLocality::Local];
    const bool          bDoTime       = gpu_nbv->bDoTime;
    const int           maxNumColumns = gridSet.numColumnsMax();

    reallocateDeviceBuffer(&gpu_nbv->cxy_na,
                           maxNumColumns * gridSet.grids().size(),
                           &gpu_nbv->ncxy_na,
                           &gpu_nbv->ncxy_na_alloc,
                           *gpu_nbv->deviceContext_);
    reallocateDeviceBuffer(&gpu_nbv->cxy_ind,
                           maxNumColumns * gridSet.grids().size(),
                           &gpu_nbv->ncxy_ind,
                           &gpu_nbv->ncxy_ind_alloc,
                           *gpu_nbv->deviceContext_);

    for (unsigned int g = 0; g < gridSet.grids().size(); g++)
    {
        const Grid& grid = gridSet.grid(g);

        const int  numColumns      = grid.numColumns();
        const int* atomIndices     = gridSet.atomIndices().data();
        const int  atomIndicesSize = gridSet.atomIndices().size();
        const int* cxy_na          = grid.cxy_na().data();
        const int* cxy_ind         = grid.cxy_ind().data();

        auto* timerH2D = bDoTime ? &gpu_nbv->timers->xf[AtomLocality::Local].nb_h2d : nullptr;

        reallocateDeviceBuffer(&gpu_nbv->atomIndices,
                               atomIndicesSize,
                               &gpu_nbv->atomIndicesSize,
                               &gpu_nbv->atomIndicesSize_alloc,
                               *gpu_nbv->deviceContext_);

        if (atomIndicesSize > 0)
        {
            if (bDoTime)
            {
                timerH2D->openTimingRegion(localStream);
            }

            copyToDeviceBuffer(&gpu_nbv->atomIndices,
                               atomIndices,
                               0,
                               atomIndicesSize,
                               localStream,
                               GpuApiCallBehavior::Async,
                               bDoTime ? timerH2D->fetchNextEvent() : nullptr);

            if (bDoTime)
            {
                timerH2D->closeTimingRegion(localStream);
            }
        }

        if (numColumns > 0)
        {
            if (bDoTime)
            {
                timerH2D->openTimingRegion(localStream);
            }

            copyToDeviceBuffer(&gpu_nbv->cxy_na,
                               cxy_na,
                               maxNumColumns * g,
                               numColumns,
                               localStream,
                               GpuApiCallBehavior::Async,
                               bDoTime ? timerH2D->fetchNextEvent() : nullptr);

            if (bDoTime)
            {
                timerH2D->closeTimingRegion(localStream);
            }

            if (bDoTime)
            {
                timerH2D->openTimingRegion(localStream);
            }

            copyToDeviceBuffer(&gpu_nbv->cxy_ind,
                               cxy_ind,
                               maxNumColumns * g,
                               numColumns,
                               localStream,
                               GpuApiCallBehavior::Async,
                               bDoTime ? timerH2D->fetchNextEvent() : nullptr);

            if (bDoTime)
            {
                timerH2D->closeTimingRegion(localStream);
            }
        }
    }

    // The above data is transferred on the local stream but is a
    // dependency of the nonlocal stream (specifically the nonlocal X
    // buf ops kernel).  We therefore set a dependency to ensure
    // that the nonlocal stream waits on the local stream here.
    // This call records an event in the local stream:
    nbnxnInsertNonlocalGpuDependency(gpu_nbv, InteractionLocality::Local);
    // ...and this call instructs the nonlocal stream to wait on that event:
    nbnxnInsertNonlocalGpuDependency(gpu_nbv, InteractionLocality::NonLocal);
}

//! This function is documented in the header file
void gpu_free(NbnxmGpu* nb)
{
    if (nb == nullptr)
    {
        return;
    }

    /* Synchronize to make sure there are no leftover operations in the streams.
     * Ideally, there should not be any, but just in case we synchronize there
     * to avoid freeing buffers that can be still used. See #4519.
     *
     * In practice, that's not needed for CUDA since cudaFree synchronizes internally.
     * But explicitly waiting for tasks to complete before freeing the memory they use is logically
     * sound and should not have any performance impact.
     */
    nb->deviceStreams[InteractionLocality::Local]->synchronize();
    if (nb->deviceStreams[InteractionLocality::NonLocal])
    {
        nb->deviceStreams[InteractionLocality::NonLocal]->synchronize();
    }

    gpu_free_platform_specific(nb);

    delete nb->timers;
    sfree(nb->timings);

    NBAtomDataGpu* atdat   = nb->atdat;
    NBParamGpu*    nbparam = nb->nbparam;

    /* Free atdat */
    freeDeviceBuffer(&(nb->atdat->xq));
    freeDeviceBuffer(&(nb->atdat->f));
    freeDeviceBuffer(&(nb->atdat->eLJ));
    freeDeviceBuffer(&(nb->atdat->eElec));
    freeDeviceBuffer(&(nb->atdat->fShift));
    freeDeviceBuffer(&(nb->atdat->shiftVec));
    if (useLjCombRule(nb->nbparam->vdwType))
    {
        freeDeviceBuffer(&atdat->ljComb);
    }
    else
    {
        freeDeviceBuffer(&atdat->atomTypes);
    }

    /* Free nbparam */
    if (nbparam->elecType == ElecType::EwaldTab || nbparam->elecType == ElecType::EwaldTabTwin)
    {
        destroyParamLookupTable(&nbparam->coulomb_tab, &nbparam->coulomb_tab_texobj);
    }

    if (!useLjCombRule(nb->nbparam->vdwType))
    {
        destroyParamLookupTable(&nbparam->nbfp, &nbparam->nbfp_texobj);
    }

    if (nbparam->vdwType == VdwType::EwaldGeom || nbparam->vdwType == VdwType::EwaldLB)
    {
        destroyParamLookupTable(&nbparam->nbfp_comb, &nbparam->nbfp_comb_texobj);
    }

    freeDeviceBuffer(&nb->cxy_na);
    freeDeviceBuffer(&nb->cxy_ind);
    freeDeviceBuffer(&nb->atomIndices);

    delete atdat;
    delete nbparam;
    delete nb;

    if (debug)
    {
        fprintf(debug, "Cleaned up NBNXM GPU data structures.\n");
    }
}

NBAtomDataGpu* gpuGetNBAtomData(NbnxmGpu* nb)
{
    GMX_ASSERT(nb != nullptr, "nb pointer must be valid");
    return nb->atdat;
}

DeviceBuffer<RVec> gpu_get_f(NbnxmGpu* nb)
{
    GMX_ASSERT(nb != nullptr, "nb pointer must be valid");

    return nb->atdat->f;
}

} // namespace gmx
