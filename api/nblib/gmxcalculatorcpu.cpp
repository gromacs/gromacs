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
#include "nblib/gmxcalculatorcpu.h"

#include "gromacs/ewald/ewald_utils.h"
#include "gromacs/mdtypes/enerdata.h"
#include "gromacs/mdtypes/interaction_const.h"
#include "gromacs/nbnxm/atomdata.h"
#include "gromacs/nbnxm/nbnxm.h"
#include "gromacs/nbnxm/pairlistset.h"
#include "gromacs/nbnxm/pairlistsets.h"
#include "gromacs/nbnxm/pairsearch.h"
#include "gromacs/utility/listoflists.h"
#include "gromacs/utility/range.h"

#include "nblib/exception.h"
#include "nblib/nbnxmsetuphelpers.h"
#include "nblib/topology.h"
#include "nblib/tpr.h"

#include "gmxbackenddata.h"
#include "pbc.hpp"
#include "systemdescription.h"
#include "virials.h"

namespace nblib
{

class GmxNBForceCalculatorCpu::CpuImpl final
{
public:
    CpuImpl(gmx::ArrayRef<int>     particleTypeIdOfAllParticles,
            gmx::ArrayRef<real>    nonBondedParams,
            gmx::ArrayRef<real>    charges,
            gmx::ArrayRef<int64_t> particleInteractionFlags,
            gmx::ArrayRef<int>     exclusionRanges,
            gmx::ArrayRef<int>     exclusionElements,
            const NBKernelOptions& options);

    //! calculates a new pair list based on new coordinates (for every NS step)
    void updatePairlist(gmx::ArrayRef<gmx::RVec> coordinates, const Box& box);

    //! Compute forces and return
    void compute(gmx::ArrayRef<const gmx::RVec> coordinateInput,
                 const Box&                     box,
                 gmx::ArrayRef<gmx::RVec>       forceOutput);

    //! Compute forces and virial tensor
    void compute(gmx::ArrayRef<const gmx::RVec> coordinateInput,
                 const Box&                     box,
                 gmx::ArrayRef<gmx::RVec>       forceOutput,
                 gmx::ArrayRef<real>            virialOutput);

    //! Compute forces, virial tensor and potential energies
    void compute(gmx::ArrayRef<const gmx::RVec> coordinateInput,
                 const Box&                     box,
                 gmx::ArrayRef<gmx::RVec>       forceOutput,
                 gmx::ArrayRef<real>            virialOutput,
                 gmx::ArrayRef<real>            energyOutput);

private:
    //! \brief client-side provided system description data
    SystemDescription system_;

    //! \brief Gmx backend objects, employed for calculating the forces
    GmxBackendData backend_;
};

GmxNBForceCalculatorCpu::CpuImpl::CpuImpl(gmx::ArrayRef<int>     particleTypeIdOfAllParticles,
                                          gmx::ArrayRef<real>    nonBondedParams,
                                          gmx::ArrayRef<real>    charges,
                                          gmx::ArrayRef<int64_t> particleInteractionFlags,
                                          gmx::ArrayRef<int>     exclusionRanges,
                                          gmx::ArrayRef<int>     exclusionElements,
                                          const NBKernelOptions& options) :
    system_(SystemDescription(particleTypeIdOfAllParticles, nonBondedParams, charges, particleInteractionFlags)),
    backend_(GmxBackendData(options, findNumEnergyGroups(particleInteractionFlags), exclusionRanges, exclusionElements))
{
    // Set up non-bonded verlet in the backend
    backend_.nbv_ = createNbnxmCPU(system_.numParticleTypes_,
                                   options,
                                   findNumEnergyGroups(particleInteractionFlags),
                                   system_.nonBondedParams_);
}

void GmxNBForceCalculatorCpu::CpuImpl::updatePairlist(gmx::ArrayRef<gmx::RVec> coordinates, const Box& box)
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

    const real particleDensity = static_cast<real>(coordinates.size()) / det(legacyBox);

    // If particles are too far outside the box, the grid setup can fail
    put_atoms_in_box_omp(
            PbcType::Xyz, box.legacyMatrix(), false, nullptr, coordinates, {}, backend_.numThreads_);

    // Put particles on a grid based on bounds specified by the box
    backend_.nbv_->putAtomsOnGrid(legacyBox,
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

    backend_.nbv_->constructPairlist(
            gmx::InteractionLocality::Local, backend_.exclusions_, 0, &backend_.nrnb_);

    // Set Particle Types and Charges and VdW params
    backend_.nbv_->setAtomProperties(
            system_.particleTypeIdOfAllParticles_, system_.charges_, system_.particleInfo_);
    backend_.updatePairlistCalled = true;
}

void GmxNBForceCalculatorCpu::CpuImpl::compute(gmx::ArrayRef<const gmx::RVec> coordinateInput,
                                               const Box&                     box,
                                               gmx::ArrayRef<gmx::RVec>       forceOutput,
                                               gmx::ArrayRef<real>            virialOutput,
                                               gmx::ArrayRef<real>            energyOutput)
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

    bool computeVirial               = !virialOutput.empty();
    bool computeEnergies             = !energyOutput.empty();
    backend_.stepWork_.computeVirial = computeVirial;
    backend_.stepWork_.computeEnergy = computeEnergies;

    // update the coordinates in the backend
    backend_.nbv_->convertCoordinates(gmx::AtomLocality::Local, coordinateInput);

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

    backend_.nbv_->atomdata_add_nbat_f_to_f(gmx::AtomLocality::All, forceOutput);

    if (computeVirial)
    {
        // calculate shift forces and turn into an array ref
        std::vector<Vec3> shiftForcesVector(gmx::c_numShiftVectors, Vec3(0.0, 0.0, 0.0));
        nbnxn_atomdata_add_nbat_fshift_to_fshift(backend_.nbv_->nbat(), shiftForcesVector);
        auto shiftForcesRef = constArrayRefFromArray(shiftForcesVector.data(), shiftForcesVector.size());

        std::vector<Vec3> shiftVectorsArray(gmx::c_numShiftVectors);

        // copy shift vectors from ForceRec
        std::copy(backend_.forcerec_.shift_vec.begin(),
                  backend_.forcerec_.shift_vec.end(),
                  shiftVectorsArray.begin());

        computeVirialTensor(
                coordinateInput, forceOutput, shiftVectorsArray, shiftForcesRef, box, virialOutput);
    }

    // extract term energies (per interaction type)
    if (computeEnergies)
    {
        int nGroupPairs = backend_.enerd_.grpp.nener;
        if (int(energyOutput.size()) != int(NonBondedEnergyTerms::Count) * nGroupPairs)
        {
            throw InputException("Array size for energy output is wrong\n");
        }

        for (int eg = 0; eg < int(NonBondedEnergyTerms::Count); ++eg)
        {
            std::copy(begin(backend_.enerd_.grpp.energyGroupPairTerms[eg]),
                      end(backend_.enerd_.grpp.energyGroupPairTerms[eg]),
                      energyOutput.begin() + eg * nGroupPairs);
        }
    }
}

void GmxNBForceCalculatorCpu::CpuImpl::compute(gmx::ArrayRef<const gmx::RVec> coordinateInput,
                                               const Box&                     box,
                                               gmx::ArrayRef<gmx::RVec>       forceOutput)
{
    // compute forces and fill in force buffer
    compute(coordinateInput, box, forceOutput, gmx::ArrayRef<real>{}, gmx::ArrayRef<real>{});
}

void GmxNBForceCalculatorCpu::CpuImpl::compute(gmx::ArrayRef<const gmx::RVec> coordinateInput,
                                               const Box&                     box,
                                               gmx::ArrayRef<gmx::RVec>       forceOutput,
                                               gmx::ArrayRef<real>            virialOutput)
{
    // compute forces and fill in force buffer
    compute(coordinateInput, box, forceOutput, virialOutput, gmx::ArrayRef<real>{});
}

GmxNBForceCalculatorCpu::GmxNBForceCalculatorCpu(gmx::ArrayRef<int>  particleTypeIdOfAllParticles,
                                                 gmx::ArrayRef<real> nonBondedParams,
                                                 gmx::ArrayRef<real> charges,
                                                 gmx::ArrayRef<int64_t> particleInteractionFlags,
                                                 gmx::ArrayRef<int>     exclusionRanges,
                                                 gmx::ArrayRef<int>     exclusionElements,
                                                 const NBKernelOptions& options)
{
    if (options.useGpu)
    {
        throw InputException("Use GmxNBForceCalculatorGpu for GPU support");
    }

    impl_ = std::make_unique<CpuImpl>(particleTypeIdOfAllParticles,
                                      nonBondedParams,
                                      charges,
                                      particleInteractionFlags,
                                      exclusionRanges,
                                      exclusionElements,
                                      options);
}

GmxNBForceCalculatorCpu::~GmxNBForceCalculatorCpu() = default;

//! calculates a new pair list based on new coordinates (for every NS step)
void GmxNBForceCalculatorCpu::updatePairlist(gmx::ArrayRef<gmx::RVec> coordinates, const Box& box)
{
    impl_->updatePairlist(coordinates, box);
}

//! Compute forces and return
void GmxNBForceCalculatorCpu::compute(gmx::ArrayRef<const gmx::RVec> coordinateInput,
                                      const Box&                     box,
                                      gmx::ArrayRef<gmx::RVec>       forceOutput)
{
    impl_->compute(coordinateInput, box, forceOutput);
}

//! Compute forces and virial tensor
void GmxNBForceCalculatorCpu::compute(gmx::ArrayRef<const gmx::RVec> coordinateInput,
                                      const Box&                     box,
                                      gmx::ArrayRef<gmx::RVec>       forceOutput,
                                      gmx::ArrayRef<real>            virialOutput)
{
    impl_->compute(coordinateInput, box, forceOutput, virialOutput);
}

//! Compute forces, virial tensor and potential energies
void GmxNBForceCalculatorCpu::compute(gmx::ArrayRef<const gmx::RVec> coordinateInput,
                                      const Box&                     box,
                                      gmx::ArrayRef<gmx::RVec>       forceOutput,
                                      gmx::ArrayRef<real>            virialOutput,
                                      gmx::ArrayRef<real>            energyOutput)
{
    impl_->compute(coordinateInput, box, forceOutput, virialOutput, energyOutput);
}

std::unique_ptr<GmxNBForceCalculatorCpu> setupGmxForceCalculatorCpu(const Topology&        topology,
                                                                    const NBKernelOptions& options)
{
    std::vector<real> nonBondedParameters = createNonBondedParameters(
            topology.getParticleTypes(), topology.getNonBondedInteractionMap());

    std::vector<int64_t> particleInteractionFlags = createParticleInfoAllVdw(topology.numParticles());

    return std::make_unique<GmxNBForceCalculatorCpu>(topology.getParticleTypeIdOfAllParticles(),
                                                     nonBondedParameters,
                                                     topology.getCharges(),
                                                     particleInteractionFlags,
                                                     topology.exclusionLists().ListRanges,
                                                     topology.exclusionLists().ListElements,
                                                     options);
}

std::unique_ptr<GmxNBForceCalculatorCpu> setupGmxForceCalculatorCpu(TprReader& tprReader,
                                                                    const NBKernelOptions& options)
{
    return std::make_unique<GmxNBForceCalculatorCpu>(tprReader.particleTypeIdOfAllParticles_,
                                                     tprReader.nonbondedParameters_,
                                                     tprReader.charges_,
                                                     tprReader.particleInteractionFlags_,
                                                     tprReader.exclusionListRanges_,
                                                     tprReader.exclusionListElements_,
                                                     options);
}

} // namespace nblib
