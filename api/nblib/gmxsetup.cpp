/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2020, by the GROMACS development team, led by
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
 * \brief Translation layer to GROMACS data structures for force calculations.
 *
 * \author Victor Holanda <victor.holanda@cscs.ch>
 * \author Joe Jordan <ejjordan@kth.se>
 * \author Prashanth Kanduri <kanduri@cscs.ch>
 * \author Sebastian Keller <keller@cscs.ch>
 */
#include "nblib/gmxsetup.h"
#include "gromacs/ewald/ewald_utils.h"
#include "gromacs/gmxlib/nrnb.h"
#include "gromacs/mdlib/forcerec.h"
#include "gromacs/mdlib/gmx_omp_nthreads.h"
#include "gromacs/mdlib/rf_util.h"
#include "gromacs/mdtypes/forcerec.h"
#include "gromacs/mdtypes/interaction_const.h"
#include "gromacs/mdtypes/simulation_workload.h"
#include "gromacs/nbnxm/atomdata.h"
#include "gromacs/nbnxm/nbnxm.h"
#include "gromacs/nbnxm/nbnxm_simd.h"
#include "gromacs/nbnxm/pairlistset.h"
#include "gromacs/nbnxm/pairlistsets.h"
#include "gromacs/nbnxm/pairsearch.h"
#include "gromacs/pbcutil/pbc.h"
#include "gromacs/utility/logger.h"
#include "gromacs/utility/smalloc.h"
#include "nblib/exception.h"
#include "nblib/kerneloptions.h"
#include "nblib/particletype.h"
#include "nblib/simulationstate.h"

namespace nblib
{

//! Helper to translate between the different enumeration values.
static Nbnxm::KernelType translateBenchmarkEnum(const SimdKernels& kernel)
{
    int kernelInt = static_cast<int>(kernel);
    return static_cast<Nbnxm::KernelType>(kernelInt);
}

/*! \brief Checks the kernel setup
 *
 * Returns an error string when the kernel is not available.
 */
static void checkKernelSetup(const NBKernelOptions& options)
{
    if (options.nbnxmSimd >= SimdKernels::Count && options.nbnxmSimd == SimdKernels::SimdAuto)
    {
        throw InputException("Need a valid kernel SIMD type");
    }
    // Check SIMD support
    if ((options.nbnxmSimd != SimdKernels::SimdNo && !GMX_SIMD)
#ifndef GMX_NBNXN_SIMD_4XN
        || options.nbnxmSimd == SimdKernels::Simd4XM
#endif
#ifndef GMX_NBNXN_SIMD_2XNN
        || options.nbnxmSimd == SimdKernels::Simd2XMM
#endif
    )
    {
        throw InputException("The requested SIMD kernel was not set up at configuration time");
    }
}

NbvSetupUtil::NbvSetupUtil() : gmxForceCalculator_(std::make_unique<GmxForceCalculator>()) {}

void NbvSetupUtil::setExecutionContext(const NBKernelOptions& options)
{
    // Todo: find a more general way to initialize hardware
    gmx_omp_nthreads_set(emntPairsearch, options.numOpenMPThreads);
    gmx_omp_nthreads_set(emntNonbonded, options.numOpenMPThreads);
}

Nbnxm::KernelSetup NbvSetupUtil::getKernelSetup(const NBKernelOptions& options)
{
    checkKernelSetup(options);

    Nbnxm::KernelSetup kernelSetup;

    // The int enum options.nbnxnSimd is set up to match Nbnxm::KernelType + 1
    kernelSetup.kernelType = translateBenchmarkEnum(options.nbnxmSimd);
    // The plain-C kernel does not support analytical ewald correction
    if (kernelSetup.kernelType == Nbnxm::KernelType::Cpu4x4_PlainC)
    {
        kernelSetup.ewaldExclusionType = Nbnxm::EwaldExclusionType::Table;
    }
    else
    {
        kernelSetup.ewaldExclusionType = options.useTabulatedEwaldCorr
                                                 ? Nbnxm::EwaldExclusionType::Table
                                                 : Nbnxm::EwaldExclusionType::Analytical;
    }

    return kernelSetup;
}

void NbvSetupUtil::setParticleInfoAllVdv(const size_t numParticles)

{
    particleInfoAllVdw_.resize(numParticles);
    for (size_t particleI = 0; particleI < numParticles; particleI++)
    {
        SET_CGINFO_HAS_VDW(particleInfoAllVdw_[particleI]);
        SET_CGINFO_HAS_Q(particleInfoAllVdw_[particleI]);
    }
}

void NbvSetupUtil::setNonBondedParameters(const std::vector<ParticleType>& particleTypes,
                                          const NonBondedInteractionMap&   nonBondedInteractionMap)
{
    /* Todo: Refactor nbnxm to take nonbondedParameters_ directly
     *
     * initial self-handling of combination rules
     * size: 2*(numParticleTypes^2)
     */
    nonbondedParameters_.reserve(2 * particleTypes.size() * particleTypes.size());

    constexpr real c6factor  = 6.0;
    constexpr real c12factor = 12.0;

    for (const ParticleType& particleType1 : particleTypes)
    {
        for (const ParticleType& particleType2 : particleTypes)
        {
            nonbondedParameters_.push_back(
                    nonBondedInteractionMap.getC6(particleType1.name(), particleType2.name()) * c6factor);
            nonbondedParameters_.push_back(
                    nonBondedInteractionMap.getC12(particleType1.name(), particleType2.name()) * c12factor);
        }
    }
}

void NbvSetupUtil::setAtomProperties(const std::vector<int>&  particleTypeIdOfAllParticles,
                                     const std::vector<real>& charges)
{
    gmxForceCalculator_->nbv_->setAtomProperties(particleTypeIdOfAllParticles, charges, particleInfoAllVdw_);
}

//! Sets up and returns a Nbnxm object for the given options and system
void NbvSetupUtil::setupNbnxmInstance(const size_t numParticleTypes, const NBKernelOptions& options)
{
    const auto pinPolicy  = (options.useGpu ? gmx::PinningPolicy::PinnedIfSupported
                                           : gmx::PinningPolicy::CannotBePinned);
    const int  numThreads = options.numOpenMPThreads;
    // Note: the options and Nbnxm combination rule enums values should match
    const int combinationRule = static_cast<int>(options.ljCombinationRule);

    checkKernelSetup(options); // throws exception is setup is invalid

    Nbnxm::KernelSetup kernelSetup = getKernelSetup(options);

    PairlistParams pairlistParams(kernelSetup.kernelType, false, options.pairlistCutoff, false);
    Nbnxm::GridSet gridSet(PbcType::Xyz, false, nullptr, nullptr, pairlistParams.pairlistType,
                           false, numThreads, pinPolicy);
    auto           pairlistSets = std::make_unique<PairlistSets>(pairlistParams, false, 0);
    auto           pairSearch =
            std::make_unique<PairSearch>(PbcType::Xyz, false, nullptr, nullptr,
                                         pairlistParams.pairlistType, false, numThreads, pinPolicy);

    auto atomData = std::make_unique<nbnxn_atomdata_t>(pinPolicy);

    // Put everything together
    auto nbv = std::make_unique<nonbonded_verlet_t>(std::move(pairlistSets), std::move(pairSearch),
                                                    std::move(atomData), kernelSetup, nullptr, nullptr);

    // Needs to be called with the number of unique ParticleTypes
    nbnxn_atomdata_init(gmx::MDLogger(), nbv->nbat.get(), kernelSetup.kernelType, combinationRule,
                        numParticleTypes, nonbondedParameters_, 1, numThreads);

    gmxForceCalculator_->nbv_ = std::move(nbv);
}

//! Computes the Ewald splitting coefficient for Coulomb
static real ewaldCoeff(const real ewald_rtol, const real pairlistCutoff)
{
    return calc_ewaldcoeff_q(pairlistCutoff, ewald_rtol);
}

void NbvSetupUtil::setupStepWorkload(const NBKernelOptions& options)
{
    gmxForceCalculator_->stepWork_->computeForces          = true;
    gmxForceCalculator_->stepWork_->computeNonbondedForces = true;

    if (options.computeVirialAndEnergy)
    {
        gmxForceCalculator_->stepWork_->computeVirial = true;
        gmxForceCalculator_->stepWork_->computeEnergy = true;
    }
}

void NbvSetupUtil::setupInteractionConst(const NBKernelOptions& options)
{
    gmxForceCalculator_->interactionConst_->vdwtype      = evdwCUT;
    gmxForceCalculator_->interactionConst_->vdw_modifier = eintmodPOTSHIFT;
    gmxForceCalculator_->interactionConst_->rvdw         = options.pairlistCutoff;

    switch (options.coulombType)
    {
        case CoulombType::Pme: gmxForceCalculator_->interactionConst_->eeltype = eelPME; break;
        case CoulombType::Cutoff: gmxForceCalculator_->interactionConst_->eeltype = eelCUT; break;
        case CoulombType::ReactionField:
            gmxForceCalculator_->interactionConst_->eeltype = eelRF;
            break;
        case CoulombType::Count: throw InputException("Unsupported electrostatic interaction");
    }
    gmxForceCalculator_->interactionConst_->coulomb_modifier = eintmodPOTSHIFT;
    gmxForceCalculator_->interactionConst_->rcoulomb         = options.pairlistCutoff;
    // Note: values correspond to ic->coulomb_modifier = eintmodPOTSHIFT
    gmxForceCalculator_->interactionConst_->dispersion_shift.cpot =
            -1.0 / gmx::power6(gmxForceCalculator_->interactionConst_->rvdw);
    gmxForceCalculator_->interactionConst_->repulsion_shift.cpot =
            -1.0 / gmx::power12(gmxForceCalculator_->interactionConst_->rvdw);

    // These are the initialized values but we leave them here so that later
    // these can become options.
    gmxForceCalculator_->interactionConst_->epsilon_r  = 1.0;
    gmxForceCalculator_->interactionConst_->epsilon_rf = 1.0;

    /* Set the Coulomb energy conversion factor */
    if (gmxForceCalculator_->interactionConst_->epsilon_r != 0)
    {
        gmxForceCalculator_->interactionConst_->epsfac =
                ONE_4PI_EPS0 / gmxForceCalculator_->interactionConst_->epsilon_r;
    }
    else
    {
        /* eps = 0 is infinite dieletric: no Coulomb interactions */
        gmxForceCalculator_->interactionConst_->epsfac = 0;
    }

    calc_rffac(nullptr, gmxForceCalculator_->interactionConst_->epsilon_r,
               gmxForceCalculator_->interactionConst_->epsilon_rf,
               gmxForceCalculator_->interactionConst_->rcoulomb,
               &gmxForceCalculator_->interactionConst_->k_rf,
               &gmxForceCalculator_->interactionConst_->c_rf);

    if (EEL_PME_EWALD(gmxForceCalculator_->interactionConst_->eeltype))
    {
        // Ewald coefficients, we ignore the potential shift
        gmxForceCalculator_->interactionConst_->ewaldcoeff_q = ewaldCoeff(1e-5, options.pairlistCutoff);
        if (gmxForceCalculator_->interactionConst_->ewaldcoeff_q <= 0)
        {
            throw InputException("Ewald coefficient should be > 0");
        }
        gmxForceCalculator_->interactionConst_->coulombEwaldTables =
                std::make_unique<EwaldCorrectionTables>();
        init_interaction_const_tables(nullptr, gmxForceCalculator_->interactionConst_.get(), 0);
    }
}

void NbvSetupUtil::setupForceRec(const matrix& box)
{
    assert((gmxForceCalculator_->forcerec_ && "Forcerec not initialized"));
    gmxForceCalculator_->forcerec_->nbfp = nonbondedParameters_;
    snew(gmxForceCalculator_->forcerec_->shift_vec, numShiftVectors);
    calc_shifts(box, gmxForceCalculator_->forcerec_->shift_vec);
}

void NbvSetupUtil::setParticlesOnGrid(const std::vector<Vec3>& coordinates, const Box& box)
{
    gmxForceCalculator_->setParticlesOnGrid(particleInfoAllVdw_, coordinates, box);
}

void NbvSetupUtil::constructPairList(const gmx::ListOfLists<int>& exclusions)
{
    gmxForceCalculator_->nbv_->constructPairlist(gmx::InteractionLocality::Local, exclusions, 0,
                                                 gmxForceCalculator_->nrnb_.get());
}


std::unique_ptr<GmxForceCalculator> GmxSetupDirector::setupGmxForceCalculator(const SimulationState& system,
                                                                              const NBKernelOptions& options)
{
    NbvSetupUtil nbvSetupUtil;
    nbvSetupUtil.setExecutionContext(options);
    nbvSetupUtil.setNonBondedParameters(system.topology().getParticleTypes(),
                                        system.topology().getNonBondedInteractionMap());
    nbvSetupUtil.setParticleInfoAllVdv(system.topology().numParticles());

    nbvSetupUtil.setupInteractionConst(options);
    nbvSetupUtil.setupStepWorkload(options);
    nbvSetupUtil.setupNbnxmInstance(system.topology().getParticleTypes().size(), options);
    nbvSetupUtil.setParticlesOnGrid(system.coordinates(), system.box());
    nbvSetupUtil.constructPairList(system.topology().getGmxExclusions());
    nbvSetupUtil.setAtomProperties(system.topology().getParticleTypeIdOfAllParticles(),
                                   system.topology().getCharges());
    nbvSetupUtil.setupForceRec(system.box().legacyMatrix());

    return nbvSetupUtil.getGmxForceCalculator();
}

} // namespace nblib
