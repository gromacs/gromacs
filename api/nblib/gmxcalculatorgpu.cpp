/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2020- The GROMACS Authors
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
 * \brief Implements a force calculator based on GROMACS data structures.
 *
 * \author Victor Holanda <victor.holanda@cscs.ch>
 * \author Joe Jordan <ejjordan@kth.se>
 * \author Prashanth Kanduri <kanduri@cscs.ch>
 * \author Sebastian Keller <keller@cscs.ch>
 */
#include "nblib/gmxcalculatorgpu.h"

#include "gromacs/ewald/ewald_utils.h"
#include "gromacs/gpu_utils/device_stream_manager.h"
#include "gromacs/gpu_utils/devicebuffer.h"
#include "gromacs/mdtypes/enerdata.h"
#include "gromacs/mdtypes/interaction_const.h"
#include "gromacs/nbnxm/cuda/nbnxm_cuda_types.h"
#include "gromacs/nbnxm/gpu_data_mgmt.h"
#include "gromacs/nbnxm/nbnxm.h"
#include "gromacs/nbnxm/nbnxm_gpu.h"
#include "gromacs/nbnxm/pairlistset.h"
#include "gromacs/nbnxm/pairlistsets.h"
#include "gromacs/nbnxm/pairsearch.h"
#include "gromacs/utility/listoflists.h"
#include "gromacs/utility/range.h"

#include "nblib/exception.h"
#include "nblib/nbnxmsetuphelpers.h"
#include "nblib/pbc.hpp"

#include "gmxbackenddata.h"
#include "systemdescription.h"
#include "virials.h"

namespace nblib
{

class GmxNBForceCalculatorGpu::GpuImpl final
{
public:
    GpuImpl(gmx::ArrayRef<int>       particleTypeIdOfAllParticles,
            gmx::ArrayRef<real>      nonBondedParams,
            gmx::ArrayRef<real>      charges,
            gmx::ArrayRef<int64_t>   particleInteractionFlags,
            gmx::ArrayRef<int>       exclusionRanges,
            gmx::ArrayRef<int>       exclusionElements,
            const NBKernelOptions&   options,
            const DeviceInformation& deviceInfo);

    //! calculates a new pair list based on new coordinates (for every NS step)
    void updatePairlist(gmx::ArrayRef<gmx::RVec> coordinates, const Box& box);

    //! reorder coordinates into nbnxm ordering
    void reorder(gmx::ArrayRef<const gmx::RVec> coordinateInput, gmx::ArrayRef<real> coordinateOutput);

    void undoReorder(gmx::ArrayRef<const gmx::RVec> input, gmx::ArrayRef<gmx::RVec> output);

    //! Compute forces and return
    void compute(gmx::ArrayRef<const gmx::RVec> coordinateInput,
                 const Box&                     box,
                 gmx::ArrayRef<gmx::RVec>       forceOutput);

    void compute(DeviceBuffer<Float4> coordinateInput, const Box& box, DeviceBuffer<Float3> forceOutput);

    [[nodiscard]] std::size_t nbnxmBufferSize() const;

    //! returns the device context (duh!)
    const DeviceContext& deviceContext() const { return deviceStreamManager_->context(); }

    const DeviceStream& deviceStream() const
    {
        return *backend_.nbv_->gpuNbv()->deviceStreams[gmx::InteractionLocality::Local];
    }

private:
    //! \brief client-side provided system description data
    SystemDescription system_;

    //! \brief Gmx backend objects, employed for calculating the forces
    GmxBackendData backend_;

    //! Stream and context manager
    std::shared_ptr<gmx::DeviceStreamManager> deviceStreamManager_;
};

GmxNBForceCalculatorGpu::GpuImpl::GpuImpl(gmx::ArrayRef<int>       particleTypeIdOfAllParticles,
                                          gmx::ArrayRef<real>      nonBondedParams,
                                          gmx::ArrayRef<real>      charges,
                                          gmx::ArrayRef<int64_t>   particleInteractionFlags,
                                          gmx::ArrayRef<int>       exclusionRanges,
                                          gmx::ArrayRef<int>       exclusionElements,
                                          const NBKernelOptions&   options,
                                          const DeviceInformation& deviceInfo) :
    system_(SystemDescription(particleTypeIdOfAllParticles, nonBondedParams, charges, particleInteractionFlags)),
    // multiple energy groups not supported on GPUs
    backend_(GmxBackendData(options, 1, exclusionRanges, exclusionElements))
{
    backend_.simulationWork_ = createSimulationWorkloadGpu();

    // set DeviceInformation and create the DeviceStreamManager
    deviceStreamManager_ = createDeviceStreamManager(deviceInfo, backend_.simulationWork_);

    // create the nonbonded_verlet_t instance, corresponds roughly to init_nb_verlet in GROMACS
    backend_.nbv_ = createNbnxmGPU(system_.numParticleTypes_,
                                   options,
                                   system_.nonBondedParams_,
                                   backend_.interactionConst_,
                                   *deviceStreamManager_);
}

void GmxNBForceCalculatorGpu::GpuImpl::updatePairlist(gmx::ArrayRef<gmx::RVec> coordinates, const Box& box)
{
    if (coordinates.size() != system_.numParticles_)
    {
        throw InputException(
                "Coordinate array containing different number of entries than particles in the "
                "system");
    }

    const auto* legacyBox = box.legacyMatrix();
    system_.box_          = box;
    updateForcerec(&backend_.forcerec_, box.legacyMatrix());
    if (TRICLINIC(legacyBox))
    {
        throw InputException("Only rectangular unit-cells are supported here");
    }

    const rvec lowerCorner = { 0, 0, 0 };
    const rvec upperCorner = { legacyBox[dimX][dimX], legacyBox[dimY][dimY], legacyBox[dimZ][dimZ] };

    const real particleDensity = real(coordinates.size()) / det(legacyBox);

    // If particles are too far outside the box, the grid setup can fail
    // TODO: Call OMP function without passing unused arguments being passed as dummies
    put_atoms_in_box_omp(
            PbcType::Xyz, box.legacyMatrix(), false, box.legacyMatrix(), coordinates, coordinates, backend_.numThreads_);

    // nbnxn_put_on_grid(backend_.nbv_.get(),
    backend_.nbv_.get()->putAtomsOnGrid(legacyBox,
                                        0,
                                        lowerCorner,
                                        upperCorner,
                                        nullptr,
                                        { 0, int(coordinates.size()) },
                                        particleDensity,
                                        system_.particleInfo_,
                                        coordinates,
                                        0,
                                        nullptr);

    // construct the pairlist
    backend_.nbv_->constructPairlist(
            gmx::InteractionLocality::Local, backend_.exclusions_, 0, &backend_.nrnb_);

    // particle types, charges and info flags are stored in the nbnxm_atomdata_t member of
    // nonbonded_verlet_t. For now, we replicate the current GROMACS implementation to
    // add types, charges and infos to nbnxm_atomdata_t separately from the nonbondedParameters_.
    backend_.nbv_->setAtomProperties(
            system_.particleTypeIdOfAllParticles_, system_.charges_, system_.particleInfo_);

    // allocates f, xq, atom_types and lj_comb in the nbnxm_atomdata_gpu on the GPU
    // also uploads nbfp (LJ) parameters to the GPU
    Nbnxm::gpu_init_atomdata(backend_.nbv_->gpuNbv(), &(backend_.nbv_->nbat()));

    // this is needed to for atomdata_add_nbat_f_to_f, otherwise it thinks there's no work done
    // it just sets a flag to true in case the pair list is not empty
    backend_.nbv_->setupGpuShortRangeWork(nullptr, gmx::InteractionLocality::Local);
    backend_.updatePairlistCalled = true;
}

void GmxNBForceCalculatorGpu::GpuImpl::reorder(gmx::ArrayRef<const gmx::RVec> coordinateInput,
                                               gmx::ArrayRef<real>            coordinateOutput)
{
    if (!backend_.updatePairlistCalled)
    {
        throw InputException("reorder called without updating pairlist at least once");
    }

    // input size is actual number of particles
    assert(coordinateInput.size() == system_.numParticles_);
    // output size is equal to internal padded buffer size
    assert(coordinateOutput.size() == backend_.nbv_->nbat().x().size());

    backend_.nbv_->convertCoordinates(gmx::AtomLocality::Local, coordinateInput);
    gmx::ArrayRef<real> nbatXRef = backend_.nbv_->nbat().x();
    std::copy(nbatXRef.begin(), nbatXRef.begin() + coordinateOutput.ssize(), coordinateOutput.begin());
}

void GmxNBForceCalculatorGpu::GpuImpl::undoReorder(gmx::ArrayRef<const gmx::RVec> input,
                                                   gmx::ArrayRef<gmx::RVec>       output)
{
    if (!backend_.updatePairlistCalled)
    {
        throw InputException("undoReorder called without updating pairlist at least once");
    }

    assert(output.size() == system_.numParticles_);

    int fstride = backend_.nbv_->nbat().fstride;

    const real*         flattenedInput = reinterpret_cast<const real*>(input.data());
    gmx::ArrayRef<real> nbvBuffer(backend_.nbv_->nbat().out[0].f);

    // make sure the internal nbv buffer has enough space
    assert(input.size() * fstride == nbvBuffer.size());

    // Overwrite internal nbv state. We didn't like it anyway.
    std::copy(flattenedInput, flattenedInput + input.size() * fstride, nbvBuffer.data());

    // atomdata_add_nbat_f_to_f is additive, therefore we need to zero out the output first
    std::fill(output.begin(), output.end(), gmx::RVec{ 0, 0, 0 });
    backend_.nbv_->atomdata_add_nbat_f_to_f(gmx::AtomLocality::Local, output);
}

std::size_t GmxNBForceCalculatorGpu::GpuImpl::nbnxmBufferSize() const
{

    if (!backend_.updatePairlistCalled)
    {
        throw InputException("nbnxmBufferSize called without updating pairlist at least once");
    }

    return backend_.nbv_->nbat().x().size() / backend_.nbv_->nbat().xstride;
}

void GmxNBForceCalculatorGpu::GpuImpl::compute(gmx::ArrayRef<const gmx::RVec> coordinateInput,
                                               const Box&                     box,
                                               gmx::ArrayRef<gmx::RVec>       forceOutput)
{
    if (coordinateInput.size() != forceOutput.size())
    {
        throw InputException("coordinate array and force buffer size mismatch");
    }

    if (!backend_.updatePairlistCalled)
    {
        throw InputException("compute called without updating pairlist at least once");
    }

    // update the box if changed
    if (!(system_.box_ == box))
    {
        system_.box_ = box;
        updateForcerec(&backend_.forcerec_, box.legacyMatrix());
    }

    // needed for non-neighbor search steps
    backend_.nbv_->convertCoordinates(gmx::AtomLocality::Local, coordinateInput);

    // copy the coordinates to GPU buffer
    Nbnxm::gpu_copy_xq_to_gpu(backend_.nbv_->gpuNbv(), &(backend_.nbv_->nbat()), gmx::AtomLocality::Local);

    // clear force output in the NbnxmGpu owned buffer
    bool clearVirials = false;
    Nbnxm::gpu_clear_outputs(backend_.nbv_->gpuNbv(), clearVirials);

    // compute forces
    backend_.nbv_->dispatchNonbondedKernel(
            gmx::InteractionLocality::Local,
            backend_.interactionConst_,
            backend_.stepWork_,
            enbvClearFYes,
            backend_.forcerec_.shift_vec,
            backend_.enerd_.grpp.energyGroupPairTerms[backend_.forcerec_.haveBuckingham ? NonBondedEnergyTerms::BuckinghamSR
                                                                                        : NonBondedEnergyTerms::LJSR],
            backend_.enerd_.grpp.energyGroupPairTerms[NonBondedEnergyTerms::CoulombSR],
            &backend_.nrnb_);

    // copy force from the Gpu Nbnxm into the Cpu Nbnxm
    Nbnxm::gpu_launch_cpyback(
            backend_.nbv_->gpuNbv(), &(backend_.nbv_->nbat()), backend_.stepWork_, gmx::AtomLocality::Local);

    gmx::ArrayRef<gmx::RVec> nullShiftForce;
    Nbnxm::gpu_wait_finish_task(
            backend_.nbv_->gpuNbv(),
            backend_.stepWork_,
            gmx::AtomLocality::Local,
            backend_.enerd_.grpp.energyGroupPairTerms[NonBondedEnergyTerms::LJSR].data(),
            backend_.enerd_.grpp.energyGroupPairTerms[NonBondedEnergyTerms::CoulombSR].data(),
            nullShiftForce,
            nullptr);

    // copy non-bonded forces from the nbv-internal buffer to the client provided force output buffer
    backend_.nbv_->atomdata_add_nbat_f_to_f(gmx::AtomLocality::Local, forceOutput);
}

void GmxNBForceCalculatorGpu::GpuImpl::compute(DeviceBuffer<Float4> coordinateInput,
                                               const Box& /* box */,
                                               DeviceBuffer<Float3> forceOutput)
{
    // The current nbnxm implementation is not designed for dealing with externally provided
    // device buffers. Therefore we move the nbnxm owned buffers temporarily out of the way
    // and reroute the internal buffers to the externally provided ones for the duration
    // of the non-bonded kernel call.

    // save (pointers) of nbnxm owned buffers into temporaries
    auto tmpXq = backend_.nbv_->gpuNbv()->atdat->xq;
    auto tmpF  = backend_.nbv_->gpuNbv()->atdat->f;

    // set internal buffers to externally provided ones
    backend_.nbv_->gpuNbv()->atdat->xq = coordinateInput;
    backend_.nbv_->gpuNbv()->atdat->f  = forceOutput;

    backend_.nbv_->dispatchNonbondedKernel(
            gmx::InteractionLocality::Local,
            backend_.interactionConst_,
            backend_.stepWork_,
            enbvClearFYes,
            backend_.forcerec_.shift_vec,
            backend_.enerd_.grpp.energyGroupPairTerms[backend_.forcerec_.haveBuckingham ? NonBondedEnergyTerms::BuckinghamSR
                                                                                        : NonBondedEnergyTerms::LJSR],
            backend_.enerd_.grpp.energyGroupPairTerms[NonBondedEnergyTerms::CoulombSR],
            &backend_.nrnb_);

    gmx::ArrayRef<gmx::RVec> nullShiftForce;
    Nbnxm::gpu_wait_finish_task(
            backend_.nbv_->gpuNbv(),
            backend_.stepWork_,
            gmx::AtomLocality::Local,
            backend_.enerd_.grpp.energyGroupPairTerms[NonBondedEnergyTerms::LJSR].data(),
            backend_.enerd_.grpp.energyGroupPairTerms[NonBondedEnergyTerms::CoulombSR].data(),
            nullShiftForce,
            nullptr);

    // restore nbnxm owned buffers
    backend_.nbv_->gpuNbv()->atdat->xq = tmpXq;
    backend_.nbv_->gpuNbv()->atdat->f  = tmpF;
}


GmxNBForceCalculatorGpu::GmxNBForceCalculatorGpu(gmx::ArrayRef<int>  particleTypeIdOfAllParticles,
                                                 gmx::ArrayRef<real> nonBondedParams,
                                                 gmx::ArrayRef<real> charges,
                                                 gmx::ArrayRef<int64_t>   particleInteractionFlags,
                                                 gmx::ArrayRef<int>       exclusionRanges,
                                                 gmx::ArrayRef<int>       exclusionElements,
                                                 const NBKernelOptions&   options,
                                                 const DeviceInformation& deviceInfo)
{
    impl_ = std::make_unique<GpuImpl>(particleTypeIdOfAllParticles,
                                      nonBondedParams,
                                      charges,
                                      particleInteractionFlags,
                                      exclusionRanges,
                                      exclusionElements,
                                      options,
                                      deviceInfo);
}

GmxNBForceCalculatorGpu::~GmxNBForceCalculatorGpu() = default;

//! calculates a new pair list based on new coordinates (for every NS step)
void GmxNBForceCalculatorGpu::updatePairlist(gmx::ArrayRef<gmx::RVec> coordinates, const Box& box)
{
    impl_->updatePairlist(coordinates, box);
}

void GmxNBForceCalculatorGpu::reorder(gmx::ArrayRef<const gmx::RVec> coordinateInput,
                                      gmx::ArrayRef<real>            coordinateOutput)
{
    impl_->reorder(coordinateInput, coordinateOutput);
}

void GmxNBForceCalculatorGpu::undoReorder(gmx::ArrayRef<const gmx::RVec> input,
                                          gmx::ArrayRef<gmx::RVec>       output)
{
    impl_->undoReorder(input, output);
}

void GmxNBForceCalculatorGpu::compute(gmx::ArrayRef<const gmx::RVec> coordinateInput,
                                      const Box&                     box,
                                      gmx::ArrayRef<gmx::RVec>       forceOutput) const
{
    impl_->compute(coordinateInput, box, forceOutput);
}

void GmxNBForceCalculatorGpu::compute(DeviceBuffer<Float4> coordinateInput,
                                      const Box&           box,
                                      DeviceBuffer<Float3> forceOutput) const
{
    impl_->compute(coordinateInput, box, forceOutput);
}

std::size_t GmxNBForceCalculatorGpu::nbnxmBufferSize() const
{
    return impl_->nbnxmBufferSize();
}

const DeviceContext& GmxNBForceCalculatorGpu::deviceContext() const
{
    return impl_->deviceContext();
}

const DeviceStream& GmxNBForceCalculatorGpu::deviceStream() const
{
    return impl_->deviceStream();
}

void reorderScalarArray(GmxNBForceCalculatorGpu&  calculator,
                        gmx::ArrayRef<const real> input,
                        gmx::ArrayRef<real>       output)
{
    // insert two zeros after each value to get a compatible format
    std::vector<gmx::RVec> inputExpanded(input.size(), gmx::RVec{ 0, 0, 0 });

    for (std::size_t i = 0; i < input.size(); ++i)
    {
        inputExpanded[i][0] = input[i];
    }

    std::size_t       deviceBufferSize = calculator.nbnxmBufferSize();
    std::vector<real> outputExpanded(4 * deviceBufferSize);
    calculator.reorder(inputExpanded, outputExpanded);

    assert(output.size() == deviceBufferSize);

    // compact expanded output into the final output array
    for (std::size_t i = 0; i < deviceBufferSize; ++i)
    {
        output[i] = outputExpanded[4 * i];
    }
}


std::unique_ptr<GmxNBForceCalculatorGpu> setupGmxForceCalculatorGpu(const Topology&        topology,
                                                                    const NBKernelOptions& options,
                                                                    const DeviceInformation& deviceInfo)
{
    std::vector<real> nonBondedParameters = createNonBondedParameters(
            topology.getParticleTypes(), topology.getNonBondedInteractionMap());

    std::vector<int64_t> particleInteractionFlags = createParticleInfoAllVdw(topology.numParticles());

    return std::make_unique<GmxNBForceCalculatorGpu>(topology.getParticleTypeIdOfAllParticles(),
                                                     nonBondedParameters,
                                                     topology.getCharges(),
                                                     particleInteractionFlags,
                                                     topology.exclusionLists().ListRanges,
                                                     topology.exclusionLists().ListElements,
                                                     options,
                                                     deviceInfo);
}

} // namespace nblib
