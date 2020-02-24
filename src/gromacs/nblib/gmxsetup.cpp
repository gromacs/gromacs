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
 * \brief
 * Implements the translation layer between the user scope and
 * GROMACS data structures for force calculations
 *
 * \author Victor Holanda <victor.holanda@cscs.ch>
 * \author Joe Jordan <ejjordan@kth.se>
 * \author Prashanth Kanduri <kanduri@cscs.ch>
 * \author Sebastian Keller <keller@cscs.ch>
 */

#include "simulationstate.h"
#include "gmxsetup.h"

#include "gromacs/compat/optional.h"
#include "gromacs/ewald/ewald_utils.h"
#include "gromacs/gmxlib/nrnb.h"
#include "gromacs/math/matrix.h"
#include "gromacs/math/units.h"
#include "gromacs/math/vec.h"
#include "gromacs/mdlib/gmx_omp_nthreads.h"
#include "gromacs/mdlib/rf_util.h"
#include "gromacs/mdtypes/forcerec.h"
#include "gromacs/mdtypes/mdatom.h"
#include "gromacs/nblib/atomtype.h"
#include "gromacs/nbnxm/atomdata.h"
#include "gromacs/nbnxm/nbnxm_simd.h"
#include "gromacs/nbnxm/pairlistset.h"
#include "gromacs/nbnxm/pairlistsets.h"
#include "gromacs/nbnxm/pairsearch.h"
#include "gromacs/pbcutil/ishift.h"
#include "gromacs/pbcutil/pbc.h"
#include "gromacs/simd/simd.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/logger.h"
#include "gromacs/utility/smalloc.h"

namespace nblib
{

namespace detail
{

static real combinationFunction(real v, real w, CombinationRule combinationRule)
{
    if (combinationRule == CombinationRule::Geometric)
    {
        return sqrt(v * w);
    }
    else
    {
        throw gmx::InvalidInputError("unknown LJ Combination rule specified\n");
    }
}

} // namespace detail

//! Helper to translate between the different enumeration values.
static Nbnxm::KernelType translateBenchmarkEnum(const BenchMarkKernels& kernel)
{
    int kernelInt = static_cast<int>(kernel);
    return static_cast<Nbnxm::KernelType>(kernelInt);
}


static real ewaldCoeff(const real ewald_rtol, const real pairlistCutoff)
{
    return calc_ewaldcoeff_q(pairlistCutoff, ewald_rtol);
}

/*! \brief Checks the kernel setup
 *
 * Returns an error string when the kernel is not available.
 */
static gmx::compat::optional<std::string> checkKernelSetup(const NBKernelOptions& options)
{
    GMX_RELEASE_ASSERT(options.nbnxmSimd < BenchMarkKernels::Count
                       && options.nbnxmSimd != BenchMarkKernels::SimdAuto,
                       "Need a valid kernel SIMD type");

    // Check SIMD support
    if ((options.nbnxmSimd != BenchMarkKernels::SimdNo && !GMX_SIMD)
        #ifndef GMX_NBNXN_SIMD_4XN
        || options.nbnxmSimd == BenchMarkKernels::Simd4XM
        #endif
        #ifndef GMX_NBNXN_SIMD_2XNN
        || options.nbnxmSimd == BenchMarkKernels::Simd2XMM
#endif
            )
    {
        return "the requested SIMD kernel was not set up at configuration time";
    }

    return {};
}


/*! \brief Returns the kernel setup
 */
static Nbnxm::KernelSetup getKernelSetup(const NBKernelOptions& options)
{
    auto messageWhenInvalid = checkKernelSetup(options);
    GMX_RELEASE_ASSERT(!messageWhenInvalid, "Need valid options");

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


//! Return an interaction constants struct with members used in the benchmark set appropriately
static interaction_const_t setupInteractionConst(const std::shared_ptr<NBKernelOptions> options)
{
    interaction_const_t ic;

    ic.vdwtype      = evdwCUT;
    ic.vdw_modifier = eintmodPOTSHIFT;
    ic.rvdw         = options->pairlistCutoff;

    switch (options->coulombType)
    {
        case BenchMarkCoulomb::Pme: ic.eeltype = eelPME; break;
        case BenchMarkCoulomb::Cutoff: ic.eeltype = eelCUT; break;
        case BenchMarkCoulomb::ReactionField: ic.eeltype = eelRF; break;
        case BenchMarkCoulomb::Count:
            GMX_THROW(gmx::InvalidInputError("Unsupported electrostatic interaction"));
    }
    ic.coulomb_modifier = eintmodPOTSHIFT;
    ic.rcoulomb         = options->pairlistCutoff;
    //! Note: values correspond to ic.coulomb_modifier = eintmodPOTSHIFT
    ic.dispersion_shift.cpot = -1.0 / gmx::power6(ic.rvdw);
    ic.repulsion_shift.cpot  = -1.0 / gmx::power12(ic.rvdw);

    // These are the initialized values but we leave them here so that later
    // these can become options.
    ic.epsilon_r  = 1.0;
    ic.epsilon_rf = 1.0;

    /* Set the Coulomb energy conversion factor */
    if (ic.epsilon_r != 0)
    {
        ic.epsfac = ONE_4PI_EPS0 / ic.epsilon_r;
    }
    else
    {
        /* eps = 0 is infinite dieletric: no Coulomb interactions */
        ic.epsfac = 0;
    }

    calc_rffac(nullptr, ic.epsilon_r, ic.epsilon_rf, ic.rcoulomb, &ic.k_rf, &ic.c_rf);

    if (EEL_PME_EWALD(ic.eeltype))
    {
        // Ewald coefficients, we ignore the potential shift
        ic.ewaldcoeff_q = ewaldCoeff(1e-5, options->pairlistCutoff);
        GMX_RELEASE_ASSERT(ic.ewaldcoeff_q > 0, "Ewald coefficient should be > 0");
        ic.coulombEwaldTables = std::make_unique<EwaldCorrectionTables>();
        init_interaction_const_tables(nullptr, &ic);
    }

    return ic;
}

NbvSetupUtil::NbvSetupUtil(SimulationState  system, const NBKernelOptions& options)
{
    system_ = std::make_shared<SimulationState>(system);
    options_ = std::make_shared<NBKernelOptions>(options);

    //! Todo: find a more general way to initialize hardware
    gmx_omp_nthreads_set(emntPairsearch, options.numThreads);
    gmx_omp_nthreads_set(emntNonbonded, options.numThreads);

    unpackTopologyToGmx();
}

void NbvSetupUtil::unpackTopologyToGmx()

{
    const Topology&              topology  = system_->topology();
    const std::vector<AtomType>& atomTypes = topology.getAtomTypes();

    size_t numAtoms = topology.numAtoms();

    //! Todo: Refactor nbnxm to take this (nonbondedParameters_) directly
    //!
    //! initial self-handling of combination rules
    //! size: 2*(numAtomTypes^2)
    nonbondedParameters_.reserve(2 * atomTypes.size() * atomTypes.size());
    for (const AtomType& atomType1 : atomTypes)
    {
        real c6_1  = atomType1.c6();
        real c12_1 = atomType1.c12();
        for (const AtomType& atomType2 : atomTypes)
        {
            real c6_2  = atomType2.c6();
            real c12_2 = atomType2.c12();

            real c6_combo  = detail::combinationFunction(c6_1, c6_2, CombinationRule::Geometric);
            real c12_combo = detail::combinationFunction(c12_1, c12_2, CombinationRule::Geometric);
            nonbondedParameters_.push_back(c6_combo);
            nonbondedParameters_.push_back(c12_combo);
        }
    }

    atomInfoAllVdw_.resize(numAtoms);
    for (size_t atomI = 0; atomI < numAtoms; atomI++)
    {
        SET_CGINFO_HAS_VDW(atomInfoAllVdw_[atomI]);
    }
}

//! Sets up and returns a Nbnxm object for the given options and system
std::unique_ptr<nonbonded_verlet_t> NbvSetupUtil::setupNbnxmInstance()
{
    const auto pinPolicy  = (options_->useGpu ? gmx::PinningPolicy::PinnedIfSupported
                                             : gmx::PinningPolicy::CannotBePinned);
    const int  numThreads = options_->numThreads;
    // Note: the options and Nbnxm combination rule enums values should match
    const int combinationRule = static_cast<int>(options_->ljCombinationRule);

    auto messageWhenInvalid = checkKernelSetup(*options_);
    if (messageWhenInvalid)
    {
        gmx_fatal(FARGS, "Requested kernel is unavailable because %s.", messageWhenInvalid->c_str());
    }

    Nbnxm::KernelSetup kernelSetup = getKernelSetup(*options_);

    PairlistParams pairlistParams(kernelSetup.kernelType, false, options_->pairlistCutoff, false);
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

    //! Needs to be called with the number of unique AtomTypes
    nbnxn_atomdata_init(gmx::MDLogger(), nbv->nbat.get(), kernelSetup.kernelType, combinationRule,
                        system_->topology().getAtomTypes().size(), nonbondedParameters_, 1, numThreads);

    matrix box_;
    gmx::fillLegacyMatrix(system_->box().matrix(), box_);

    GMX_RELEASE_ASSERT(!TRICLINIC(box_), "Only rectangular unit-cells are supported here");
    const rvec lowerCorner = { 0, 0, 0 };
    const rvec upperCorner = { box_[XX][XX], box_[YY][YY], box_[ZZ][ZZ] };

    const real atomDensity = system_->coordinates().size() / det(box_);

    nbnxn_put_on_grid(nbv.get(), box_, 0, lowerCorner, upperCorner, nullptr,
                      { 0, int(system_->coordinates().size()) }, atomDensity, atomInfoAllVdw_,
                      system_->coordinates(), 0, nullptr);

    t_nrnb nrnb;
    nbv->constructPairlist(gmx::InteractionLocality::Local, system_->topology().getGmxExclusions(), 0, &nrnb);

    t_mdatoms mdatoms;
    // We only use (read) the atom type and charge from mdatoms
    mdatoms.typeA   = const_cast<int*>(system_->topology().getAtomTypeIdOfAllAtoms().data());
    mdatoms.chargeA = const_cast<real*>(system_->topology().getCharges().data());
    nbv->setAtomProperties(mdatoms, atomInfoAllVdw_);

    return nbv;
}

std::unique_ptr<GmxForceCalculator> NbvSetupUtil::setupGmxForceCalculator()
{
    auto gmxForceCalculator_p = std::make_unique<GmxForceCalculator>(system_, options_);

    gmxForceCalculator_p->nbv_ = setupNbnxmInstance();

    // const PairlistSet& pairlistSet = nbv->pairlistSets().pairlistSet(gmx::InteractionLocality::Local);
    // const gmx::index numPairs = pairlistSet.natpair_ljq_ + pairlistSet.natpair_lj_ + pairlistSet.natpair_q_;
    // gmx_cycles_t cycles = gmx_cycles_read();

    matrix box_;
    gmx::fillLegacyMatrix(system_->box().matrix(), box_);

    gmxForceCalculator_p->forcerec_.nbfp  = nonbondedParameters_;
    snew(gmxForceCalculator_p->forcerec_.shift_vec, SHIFTS);
    calc_shifts(box_, gmxForceCalculator_p->forcerec_.shift_vec);

    put_atoms_in_box(PbcType::Xyz, box_, system_->coordinates());

    nbnxn_atomdata_t*                nbat = gmxForceCalculator_p->nbv_->nbat.get();
    gmxForceCalculator_p->verletForces_ = gmx::PaddedHostVector<gmx::RVec>(nbat->numAtoms(), gmx::RVec(0, 0, 0));

    return gmxForceCalculator_p;
}

GmxForceCalculator::GmxForceCalculator(const std::shared_ptr<SimulationState> system, const std::shared_ptr<NBKernelOptions> options) : enerd_(1, 0), verletForces_({})
{
    interactionConst_ = setupInteractionConst(options);

    gmx::fillLegacyMatrix(system->box().matrix(), box_);

    stepWork_.computeForces = true;
    if (options->computeVirialAndEnergy)
    {
        stepWork_.computeVirial = true;
        stepWork_.computeEnergy = true;
    }

    forcerec_.ntype = system->topology().getAtomTypes().size();
}

gmx::PaddedHostVector<gmx::RVec> GmxForceCalculator::compute()
{
    t_nrnb              nrnb = { 0 };

    nbv_->dispatchNonbondedKernel(gmx::InteractionLocality::Local, interactionConst_, stepWork_, enbvClearFNo,
                                 forcerec_, &enerd_, &nrnb);

    nbv_->atomdata_add_nbat_f_to_f(gmx::AtomLocality::All, verletForces_);

    return verletForces_;
}
} //namespace nblib
