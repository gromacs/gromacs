/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2020,2021, by the GROMACS development team, led by
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
 *  \brief
 *  Stubs of functions that must be defined by nbnxm sycl implementation.
 *
 *  \ingroup module_nbnxm
 */
#include "gmxpre.h"

#include "gromacs/gpu_utils/device_stream_manager.h"
#include "gromacs/gpu_utils/pmalloc.h"
#include "gromacs/hardware/device_information.h"
#include "gromacs/mdtypes/interaction_const.h"
#include "gromacs/nbnxm/atomdata.h"
#include "gromacs/nbnxm/gpu_data_mgmt.h"
#include "gromacs/nbnxm/nbnxm_gpu.h"
#include "gromacs/nbnxm/nbnxm_gpu_data_mgmt.h"
#include "gromacs/pbcutil/ishift.h"
#include "gromacs/utility/exceptions.h"

#include "nbnxm_sycl_types.h"

namespace Nbnxm
{

/*! \brief Initialize \p atomdata first time; it only gets filled at pair-search. */
static void initAtomdataFirst(NBAtomData*          atomdata,
                              int                  numTypes,
                              const DeviceContext& deviceContext,
                              const DeviceStream&  localStream)
{
    atomdata->numTypes = numTypes;
    allocateDeviceBuffer(&atomdata->shiftVec, SHIFTS, deviceContext);
    atomdata->shiftVecUploaded = false;

    allocateDeviceBuffer(&atomdata->fShift, SHIFTS, deviceContext);
    allocateDeviceBuffer(&atomdata->eLJ, 1, deviceContext);
    allocateDeviceBuffer(&atomdata->eElec, 1, deviceContext);

    clearDeviceBufferAsync(&atomdata->fShift, 0, SHIFTS, localStream);
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

/*! \brief Initialize the nonbonded parameter data structure. */
static void initNbparam(NBParamGpu*                     nbp,
                        const interaction_const_t&      ic,
                        const PairlistParams&           listParams,
                        const nbnxn_atomdata_t::Params& nbatParams,
                        const DeviceContext&            deviceContext)
{
    const int numTypes = nbatParams.numTypes;

    set_cutoff_parameters(nbp, &ic, listParams);

    nbp->vdwType  = nbnxmGpuPickVdwKernelType(&ic, nbatParams.ljCombinationRule);
    nbp->elecType = nbnxmGpuPickElectrostaticsKernelType(&ic, deviceContext.deviceInfo());

    /* generate table for PME */
    nbp->coulomb_tab = nullptr;
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

NbnxmGpu* gpu_init(const gmx::DeviceStreamManager& deviceStreamManager,
                   const interaction_const_t*      ic,
                   const PairlistParams&           listParams,
                   const nbnxn_atomdata_t*         nbat,
                   const bool                      bLocalAndNonlocal)
{
    auto* nb                              = new NbnxmGpu();
    nb->deviceContext_                    = &deviceStreamManager.context();
    nb->atdat                             = new NBAtomData;
    nb->nbparam                           = new NBParamGpu;
    nb->plist[InteractionLocality::Local] = new Nbnxm::gpu_plist;
    if (bLocalAndNonlocal)
    {
        nb->plist[InteractionLocality::NonLocal] = new Nbnxm::gpu_plist;
    }

    nb->bUseTwoStreams = bLocalAndNonlocal;

    nb->timers  = nullptr;
    nb->timings = nullptr;

    /* init nbst */
    pmalloc(reinterpret_cast<void**>(&nb->nbst.eLJ), sizeof(*nb->nbst.eLJ));
    pmalloc(reinterpret_cast<void**>(&nb->nbst.eElec), sizeof(*nb->nbst.eElec));
    pmalloc(reinterpret_cast<void**>(&nb->nbst.fShift), SHIFTS * sizeof(*nb->nbst.fShift));

    init_plist(nb->plist[InteractionLocality::Local]);

    /* local/non-local GPU streams */
    GMX_RELEASE_ASSERT(deviceStreamManager.streamIsValid(gmx::DeviceStreamType::NonBondedLocal),
                       "Local non-bonded stream should be initialized to use GPU for non-bonded.");
    const DeviceStream& localStream = deviceStreamManager.stream(gmx::DeviceStreamType::NonBondedLocal);
    nb->deviceStreams[InteractionLocality::Local] = &localStream;
    // In general, it's not strictly necessary to use 2 streams for SYCL, since they are
    // out-of-order. But for the time being, it will be less disruptive to keep them.
    if (nb->bUseTwoStreams)
    {
        init_plist(nb->plist[InteractionLocality::NonLocal]);

        GMX_RELEASE_ASSERT(deviceStreamManager.streamIsValid(gmx::DeviceStreamType::NonBondedNonLocal),
                           "Non-local non-bonded stream should be initialized to use GPU for "
                           "non-bonded with domain decomposition.");
        nb->deviceStreams[InteractionLocality::NonLocal] =
                &deviceStreamManager.stream(gmx::DeviceStreamType::NonBondedNonLocal);
    }

    nb->bDoTime = false;

    const nbnxn_atomdata_t::Params& nbatParams    = nbat->params();
    const DeviceContext&            deviceContext = *nb->deviceContext_;

    initNbparam(nb->nbparam, *ic, listParams, nbatParams, deviceContext);
    initAtomdataFirst(nb->atdat, nbatParams.numTypes, deviceContext, localStream);

    return nb;
}

void gpu_upload_shiftvec(NbnxmGpu* nb, const nbnxn_atomdata_t* nbatom)
{
    NBAtomData*         adat        = nb->atdat;
    const DeviceStream& localStream = *nb->deviceStreams[InteractionLocality::Local];

    /* only if we have a dynamic box */
    if (nbatom->bDynamicBox || !adat->shiftVecUploaded)
    {
        GMX_ASSERT(adat->shiftVec.elementSize() == sizeof(nbatom->shift_vec[0]),
                   "Sizes of host- and device-side shift vectors should be the same.");
        copyToDeviceBuffer(&adat->shiftVec,
                           reinterpret_cast<const Float3*>(nbatom->shift_vec.data()),
                           0,
                           SHIFTS,
                           localStream,
                           GpuApiCallBehavior::Async,
                           nullptr);
        adat->shiftVecUploaded = true;
    }
}

void gpu_free(NbnxmGpu* nb)
{
    if (nb == nullptr)
    {
        return;
    }

    NBAtomData* atdat   = nb->atdat;
    NBParamGpu* nbparam = nb->nbparam;

    if ((!nbparam->coulomb_tab)
        && (nbparam->elecType == ElecType::EwaldTab || nbparam->elecType == ElecType::EwaldTabTwin))
    {
        destroyParamLookupTable(&nbparam->coulomb_tab, nbparam->coulomb_tab_texobj);
    }

    if (!useLjCombRule(nb->nbparam->vdwType))
    {
        destroyParamLookupTable(&nbparam->nbfp, nbparam->nbfp_texobj);
    }

    if (nbparam->vdwType == VdwType::EwaldGeom || nbparam->vdwType == VdwType::EwaldLB)
    {
        destroyParamLookupTable(&nbparam->nbfp_comb, nbparam->nbfp_comb_texobj);
    }

    /* Free plist */
    auto* plist = nb->plist[InteractionLocality::Local];
    delete plist;
    if (nb->bUseTwoStreams)
    {
        auto* plist_nl = nb->plist[InteractionLocality::NonLocal];
        delete plist_nl;
    }

    /* Free nbst */
    pfree(nb->nbst.eLJ);
    nb->nbst.eLJ = nullptr;

    pfree(nb->nbst.eElec);
    nb->nbst.eElec = nullptr;

    pfree(nb->nbst.fShift);
    nb->nbst.fShift = nullptr;

    delete atdat;
    delete nbparam;
    delete nb;
}

int gpu_min_ci_balanced(NbnxmGpu* nb)
{
    // SYCL-TODO: Logic and magic values taken from OpenCL
    static constexpr unsigned int balancedFactor = 50;
    if (nb == nullptr)
    {
        return 0;
    }
    const cl::sycl::device device = nb->deviceContext_->deviceInfo().syclDevice;
    const int numComputeUnits     = device.get_info<cl::sycl::info::device::max_compute_units>();
    return balancedFactor * numComputeUnits;
}

} // namespace Nbnxm
