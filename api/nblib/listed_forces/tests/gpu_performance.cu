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
 * \brief
 * This implements basic nblib utility tests
 *
 * \author Victor Holanda <victor.holanda@cscs.ch>
 * \author Joe Jordan <ejjordan@kth.se>
 * \author Prashanth Kanduri <kanduri@cscs.ch>
 * \author Sebastian Keller <keller@cscs.ch>
 */
#include <chrono>
#include <iostream>

#include "gromacs/gpu_utils/cudautils.cuh"
#include "gromacs/gpu_utils/gputraits.cuh"
#include "gromacs/hardware/device_information.h"
#include "gromacs/hardware/device_management.h"
#include "gromacs/listed_forces/listed_forces_gpu.h"
#include "gromacs/mdtypes/simulation_workload.h"
#include "gromacs/nbnxm/gpu_types_common.h"
#include "gromacs/topology/forcefieldparameters.h"

#include "nblib/listed_forces/gpu_interface.h"

#include "buildinfo.h"

// compiling tpr without a private impl in a cuda compilation unit will trigger some (unapplicable) warnings
#include "nblib/tpr.h"

using namespace nblib;

template<class T>
std::vector<T> compareForces(gmx::ArrayRef<const util::array<T, 3>> probe,
                             gmx::ArrayRef<const util::array<T, 3>> ref)
{
    if (probe.size() != ref.size())
    {
        throw InputException("compareForces input mismatch\n");
    }

    int            numParticles = probe.size();
    std::vector<T> relErr(numParticles);

    for (size_t i = 0; i < numParticles; ++i)
    {
        T num   = norm(probe[i] - ref[i]);
        T denom = norm(ref[i]);
        T error;
        if (denom == T(0))
        {
            error = (num == T(0)) ? T(0) : INFINITY;
        }
        else
        {
            error = num / denom;
        }

        relErr[i] = error;
    }

    std::sort(begin(relErr), end(relErr));

    return relErr;
}

thrust::host_vector<util::array<real, 4>> toXyzqNblib(gmx::ArrayRef<const gmx::RVec> xyz,
                                                      gmx::ArrayRef<const real>      q)
{
    thrust::host_vector<util::array<real, 4>> ret(xyz.size());

    for (int i = 0; i < xyz.size(); ++i)
    {
        ret[i] = util::array<real, 4>{ xyz[i][0], xyz[i][1], xyz[i][2], q[i] };
    }

    return ret;
}

thrust::host_vector<Float4> toXyzqGmx(gmx::ArrayRef<const gmx::RVec> xyz, gmx::ArrayRef<const real> q)
{
    thrust::host_vector<Float4> ret(xyz.size());

    for (int i = 0; i < xyz.size(); ++i)
    {
        ret[i] = Float4{ xyz[i][0], xyz[i][1], xyz[i][2], q[i] };
    }

    return ret;
}


int main(int argc, char* argv[])
{
    assert(canPerformDeviceDetection);
    std::vector<std::unique_ptr<DeviceInformation>> devices    = findDevices();
    const DeviceInformation&                        deviceInfo = *devices[0].get();

    const std::string filepath = argv[1];
    nblib::TprReader  tpr(filepath);
    std::cout << tpr.coordinates_.size() << std::endl;

    /* input definition *************************************/

    ListedInteractionData interactions = convertToNblibInteractions(*tpr.interactionDefinitions_);

    thrust::host_vector<util::array<real, 4>> xyzqNblib = toXyzqNblib(tpr.coordinates_, tpr.charges_);
    thrust::host_vector<Float4>               xyzqGmx   = toXyzqGmx(tpr.coordinates_, tpr.charges_);
    thrust::host_vector<util::array<real, 3>> forces(xyzqNblib.size(),
                                                     util::array<real, 3>{ 0.0, 0.0, 0.0 });
    thrust::host_vector<util::array<real, 3>> shiftForces(gmx::c_numShiftVectors,
                                                          util::array<real, 3>{ 0.0, 0.0, 0.0 });

    thrust::device_vector<util::array<real, 4>> d_xyzqNblib   = xyzqNblib;
    thrust::device_vector<Float4>               d_xyzqGmx     = xyzqGmx;
    thrust::device_vector<util::array<real, 3>> d_forces      = forces;
    thrust::device_vector<util::array<real, 3>> d_shiftForces = shiftForces;

    // for comparison
    thrust::host_vector<util::array<real, 3>> nblibForces, gmxForces;
    thrust::host_vector<util::array<real, 3>> nblibShiftf, gmxShiftf;

    nblib::Box box = tpr.getBox();

    /* convert to gmx *************************************/

    auto [idef, gmx_ffparams] = convertToGmxInteractions(interactions);
    printf("idef: bonds %d angles %d pdihs %d pairs %d\n",
           idef->il[F_BONDS].size() / 3,
           idef->il[F_ANGLES].size() / 4,
           idef->il[F_PDIHS].size() / 5,
           idef->il[F_LJ14].size() / 3);
    idef->ilsort = ilsortFE_SORTED;
    for (int fType = 0; fType < idef->numNonperturbedInteractions.size(); ++fType)
    {
        idef->numNonperturbedInteractions[fType] = idef->il[fType].size();
    }

    DeviceContext deviceContext(deviceInfo);
    DeviceStream  deviceStream(deviceContext, DeviceStreamPriority::High, true);

    auto wcycle = wallcycle_init(nullptr, 0, nullptr);

    float                   elecScale     = 1.0;
    util::array<real, 4>*   xqDeviceNblib = thrust::raw_pointer_cast(d_xyzqNblib.data());
    DeviceBuffer<gmx::RVec> forceDevice =
            reinterpret_cast<gmx::RVec*>(thrust::raw_pointer_cast(d_forces.data()));
    DeviceBuffer<gmx::RVec> fshiftDevice =
            reinterpret_cast<gmx::RVec*>(thrust::raw_pointer_cast(d_shiftForces.data()));
    std::vector<int> noReorder(tpr.coordinates_.size());
    std::iota(begin(noReorder), end(noReorder), 0);

    gmx_enerdata_t ener(1, 0);
    std::fill(ener.term.begin(), ener.term.end(), 0);
    ener.grpp.clear();

    /* compute with gmx facade *************************************/
    {
        gmx::ListedForcesNblibGpuImpl calculator(
                *gmx_ffparams, elecScale, deviceContext, deviceStream, wcycle.get());

        calculator.setPbc(PbcType::Xyz, box.legacyMatrix(), true);
        calculator.updateInteractionListsAndDeviceBuffers(
                noReorder, *idef, xqDeviceNblib, forceDevice, fshiftDevice);
        auto t0 = std::chrono::high_resolution_clock::now();
        calculator.launchKernel<true, true>();
        cudaDeviceSynchronize();
        auto   t1           = std::chrono::high_resolution_clock::now();
        double nblibElapsed = std::chrono::duration<double>(t1 - t0).count();

        calculator.launchEnergyTransfer();
        calculator.waitAccumulateEnergyTerms(&ener);
        // download and save force buffers
        forces      = d_forces;
        nblibForces = forces;
        shiftForces = d_shiftForces;
        nblibShiftf = shiftForces;

        printf("%f %f %f\n", forces[0][0], forces[0][1], forces[0][2]);
        printf("%f %f %f\n", forces[1][0], forces[1][1], forces[1][2]);
        printf("potentials, elapsed %f ms:\n", nblibElapsed * 1000);
        printf("  HarmonicBond %f\n", ener.term[FindIndex<HarmonicBondType, GmxToNblibMapping>{}]);
        printf("  HarmonicAngle %f\n", ener.term[FindIndex<HarmonicAngle, GmxToNblibMapping>{}]);
        printf("  ProperDihedral %f\n", ener.term[FindIndex<ProperDihedral, GmxToNblibMapping>{}]);
        printf("  ImproperDihedral %f\n", ener.term[FindIndex<ImproperDihedral, GmxToNblibMapping>{}]);
        printf("  VanDerWaals %f\n", ener.grpp.energyGroupPairTerms[NonBondedEnergyTerms::LJ14][0]);
        printf("  Coulomb %f\n", ener.grpp.energyGroupPairTerms[NonBondedEnergyTerms::Coulomb14][0]);
    }
    /* reset buffers *************************************/

    std::fill(ener.term.begin(), ener.term.end(), 0);
    ener.grpp.clear();
    thrust::fill(d_forces.begin(), d_forces.end(), util::array<real, 3>{ 0, 0, 0 });
    thrust::fill(d_shiftForces.begin(), d_shiftForces.end(), util::array<real, 3>{ 0, 0, 0 });
    Float4* xqDeviceGmx = thrust::raw_pointer_cast(d_xyzqGmx.data());

    /* compute with gmx implementation *************************************/
    {
        gmx::ListedForcesGpu gmxCalc(
                *gmx_ffparams, elecScale, deviceInfo, deviceContext, deviceStream, wcycle.get());
        gmxCalc.setPbc(PbcType::Xyz, box.legacyMatrix(), true);
        NBAtomDataGpu nbAtomDataGpu;
        nbAtomDataGpu.xq     = xqDeviceGmx;
        nbAtomDataGpu.f      = forceDevice;
        nbAtomDataGpu.fShift = fshiftDevice;
        gmxCalc.updateInteractionListsAndDeviceBuffers(noReorder, *idef, &nbAtomDataGpu);
        gmx::StepWorkload stepWork;
        stepWork.computeVirial = true;
        stepWork.computeEnergy = true;
        auto t2                = std::chrono::high_resolution_clock::now();
        gmxCalc.launchKernel(stepWork);
        cudaDeviceSynchronize();
        auto   t3         = std::chrono::high_resolution_clock::now();
        double gmxElapsed = std::chrono::duration<double>(t3 - t2).count();

        gmxCalc.launchEnergyTransfer();
        gmxCalc.waitAccumulateEnergyTerms(&ener);
        // download and save force buffers
        forces      = d_forces;
        gmxForces   = forces;
        shiftForces = d_shiftForces;
        gmxShiftf   = shiftForces;

        printf("\n");
        printf("%f %f %f\n", forces[0][0], forces[0][1], forces[0][2]);
        printf("%f %f %f\n", forces[1][0], forces[1][1], forces[1][2]);
        printf("gmx potentials, elapsed %f ms:\n", gmxElapsed * 1000);
        printf("  HarmonicBond %f\n", ener.term[FindIndex<HarmonicBondType, GmxToNblibMapping>{}]);
        printf("  HarmonicAngle %f\n", ener.term[FindIndex<HarmonicAngle, GmxToNblibMapping>{}]);
        printf("  ProperDihedral %f\n", ener.term[FindIndex<ProperDihedral, GmxToNblibMapping>{}]);
        printf("  ImproperDihedral %f\n", ener.term[FindIndex<ImproperDihedral, GmxToNblibMapping>{}]);
        printf("  VanDerWaals %f\n", ener.grpp.energyGroupPairTerms[NonBondedEnergyTerms::LJ14][0]);
        printf("  Coulomb %f\n", ener.grpp.energyGroupPairTerms[NonBondedEnergyTerms::Coulomb14][0]);
    }

    auto errors      = compareForces<real>(nblibForces, gmxForces);
    auto shiftErrors = compareForces<real>(nblibShiftf, gmxShiftf);

    printf("\n");
    printf("max relative force error: %f\n", errors.back());
    printf("max relative shift force error: %f\n", shiftErrors.back());
}
