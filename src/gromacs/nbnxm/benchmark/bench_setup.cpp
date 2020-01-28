/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2019,2020, by the GROMACS development team, led by
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
 * This file defines functions for setting up kernel benchmarks
 *
 * \author Berk Hess <hess@kth.se>
 * \ingroup module_nbnxm
 */

#include "gmxpre.h"

#include "bench_setup.h"

#include "gromacs/compat/optional.h"
#include "gromacs/gmxlib/nrnb.h"
#include "gromacs/mdlib/dispersioncorrection.h"
#include "gromacs/mdlib/force_flags.h"
#include "gromacs/mdlib/forcerec.h"
#include "gromacs/mdlib/gmx_omp_nthreads.h"
#include "gromacs/mdtypes/enerdata.h"
#include "gromacs/mdtypes/forcerec.h"
#include "gromacs/mdtypes/interaction_const.h"
#include "gromacs/mdtypes/mdatom.h"
#include "gromacs/mdtypes/simulation_workload.h"
#include "gromacs/nbnxm/atomdata.h"
#include "gromacs/nbnxm/gridset.h"
#include "gromacs/nbnxm/nbnxm.h"
#include "gromacs/nbnxm/nbnxm_simd.h"
#include "gromacs/nbnxm/pairlistset.h"
#include "gromacs/nbnxm/pairlistsets.h"
#include "gromacs/nbnxm/pairsearch.h"
#include "gromacs/pbcutil/ishift.h"
#include "gromacs/pbcutil/pbc.h"
#include "gromacs/simd/simd.h"
#include "gromacs/timing/cyclecounter.h"
#include "gromacs/utility/enumerationhelpers.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/logger.h"

#include "bench_system.h"

namespace Nbnxm
{

/*! \brief Checks the kernel setup
 *
 * Returns an error string when the kernel is not available.
 */
static gmx::compat::optional<std::string> checkKernelSetup(const KernelBenchOptions& options)
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

//! Helper to translate between the different enumeration values.
static KernelType translateBenchmarkEnum(const BenchMarkKernels& kernel)
{
    int kernelInt = static_cast<int>(kernel);
    return static_cast<KernelType>(kernelInt);
}

/*! \brief Returns the kernel setup
 */
static KernelSetup getKernelSetup(const KernelBenchOptions& options)
{
    auto messageWhenInvalid = checkKernelSetup(options);
    GMX_RELEASE_ASSERT(!messageWhenInvalid, "Need valid options");

    KernelSetup kernelSetup;

    // The int enum options.nbnxnSimd is set up to match KernelType + 1
    kernelSetup.kernelType = translateBenchmarkEnum(options.nbnxmSimd);
    // The plain-C kernel does not support analytical ewald correction
    if (kernelSetup.kernelType == KernelType::Cpu4x4_PlainC)
    {
        kernelSetup.ewaldExclusionType = EwaldExclusionType::Table;
    }
    else
    {
        kernelSetup.ewaldExclusionType = options.useTabulatedEwaldCorr ? EwaldExclusionType::Table
                                                                       : EwaldExclusionType::Analytical;
    }

    return kernelSetup;
}

//! Return an interaction constants struct with members used in the benchmark set appropriately
static interaction_const_t setupInteractionConst(const KernelBenchOptions& options)

{
    interaction_const_t ic;

    ic.vdwtype      = evdwCUT;
    ic.vdw_modifier = eintmodPOTSHIFT;
    ic.rvdw         = options.pairlistCutoff;

    ic.eeltype          = (options.coulombType == BenchMarkCoulomb::Pme ? eelPME : eelRF);
    ic.coulomb_modifier = eintmodPOTSHIFT;
    ic.rcoulomb         = options.pairlistCutoff;

    // Reaction-field with epsilon_rf=inf
    // TODO: Replace by calc_rffac() after refactoring that
    ic.k_rf = 0.5 * std::pow(ic.rcoulomb, -3);
    ic.c_rf = 1 / ic.rcoulomb + ic.k_rf * ic.rcoulomb * ic.rcoulomb;

    if (EEL_PME_EWALD(ic.eeltype))
    {
        // Ewald coefficients, we ignore the potential shift
        GMX_RELEASE_ASSERT(options.ewaldcoeff_q > 0, "Ewald coefficient should be > 0");
        ic.ewaldcoeff_q       = options.ewaldcoeff_q;
        ic.coulombEwaldTables = std::make_unique<EwaldCorrectionTables>();
        init_interaction_const_tables(nullptr, &ic);
    }

    return ic;
}

//! Sets up and returns a Nbnxm object for the given benchmark options and system
static std::unique_ptr<nonbonded_verlet_t> setupNbnxmForBenchInstance(const KernelBenchOptions& options,
                                                                      const gmx::BenchmarkSystem& system)
{
    const auto pinPolicy  = (options.useGpu ? gmx::PinningPolicy::PinnedIfSupported
                                           : gmx::PinningPolicy::CannotBePinned);
    const int  numThreads = options.numThreads;
    // Note: the options and Nbnxm combination rule enums values should match
    const int combinationRule = static_cast<int>(options.ljCombinationRule);

    auto messageWhenInvalid = checkKernelSetup(options);
    if (messageWhenInvalid)
    {
        gmx_fatal(FARGS, "Requested kernel is unavailable because %s.", messageWhenInvalid->c_str());
    }
    Nbnxm::KernelSetup kernelSetup = getKernelSetup(options);

    PairlistParams pairlistParams(kernelSetup.kernelType, false, options.pairlistCutoff, false);

    GridSet gridSet(PbcType::Xyz, false, nullptr, nullptr, pairlistParams.pairlistType, false,
                    numThreads, pinPolicy);

    auto pairlistSets = std::make_unique<PairlistSets>(pairlistParams, false, 0);

    auto pairSearch =
            std::make_unique<PairSearch>(PbcType::Xyz, false, nullptr, nullptr,
                                         pairlistParams.pairlistType, false, numThreads, pinPolicy);

    auto atomData = std::make_unique<nbnxn_atomdata_t>(pinPolicy);

    // Put everything together
    auto nbv = std::make_unique<nonbonded_verlet_t>(std::move(pairlistSets), std::move(pairSearch),
                                                    std::move(atomData), kernelSetup, nullptr, nullptr);

    nbnxn_atomdata_init(gmx::MDLogger(), nbv->nbat.get(), kernelSetup.kernelType, combinationRule,
                        system.numAtomTypes, system.nonbondedParameters, 1, numThreads);

    t_nrnb nrnb;

    GMX_RELEASE_ASSERT(!TRICLINIC(system.box), "Only rectangular unit-cells are supported here");
    const rvec lowerCorner = { 0, 0, 0 };
    const rvec upperCorner = { system.box[XX][XX], system.box[YY][YY], system.box[ZZ][ZZ] };

    gmx::ArrayRef<const int> atomInfo;
    if (options.useHalfLJOptimization)
    {
        atomInfo = system.atomInfoOxygenVdw;
    }
    else
    {
        atomInfo = system.atomInfoAllVdw;
    }

    const real atomDensity = system.coordinates.size() / det(system.box);

    nbnxn_put_on_grid(nbv.get(), system.box, 0, lowerCorner, upperCorner, nullptr,
                      { 0, int(system.coordinates.size()) }, atomDensity, atomInfo,
                      system.coordinates, 0, nullptr);

    nbv->constructPairlist(gmx::InteractionLocality::Local, system.excls, 0, &nrnb);

    t_mdatoms mdatoms;
    // We only use (read) the atom type and charge from mdatoms
    mdatoms.typeA   = const_cast<int*>(system.atomTypes.data());
    mdatoms.chargeA = const_cast<real*>(system.charges.data());
    nbv->setAtomProperties(mdatoms, atomInfo);

    return nbv;
}

//! Add the options instance to the list for all requested kernel SIMD types
static void expandSimdOptionAndPushBack(const KernelBenchOptions&        options,
                                        std::vector<KernelBenchOptions>* optionsList)
{
    if (options.nbnxmSimd == BenchMarkKernels::SimdAuto)
    {
        bool addedInstance = false;
#ifdef GMX_NBNXN_SIMD_4XN
        optionsList->push_back(options);
        optionsList->back().nbnxmSimd = BenchMarkKernels::Simd4XM;
        addedInstance                 = true;
#endif
#ifdef GMX_NBNXN_SIMD_2XNN
        optionsList->push_back(options);
        optionsList->back().nbnxmSimd = BenchMarkKernels::Simd2XMM;
        addedInstance                 = true;
#endif
        if (!addedInstance)
        {
            optionsList->push_back(options);
            optionsList->back().nbnxmSimd = BenchMarkKernels::SimdNo;
        }
    }
    else
    {
        optionsList->push_back(options);
    }
}

//! Sets up and runs the requested benchmark instance and prints the results
//
// When \p doWarmup is true runs the warmup iterations instead
// of the normal ones and does not print any results
static void setupAndRunInstance(const gmx::BenchmarkSystem& system,
                                const KernelBenchOptions&   options,
                                const bool                  doWarmup)
{
    // Generate an, accurate, estimate of the number of non-zero pair interactions
    const real atomDensity = system.coordinates.size() / det(system.box);
    const real numPairsWithinCutoff =
            atomDensity * 4.0 / 3.0 * M_PI * std::pow(options.pairlistCutoff, 3);
    const real numUsefulPairs = system.coordinates.size() * 0.5 * (numPairsWithinCutoff + 1);

    std::unique_ptr<nonbonded_verlet_t> nbv = setupNbnxmForBenchInstance(options, system);

    // We set the interaction cut-off to the pairlist cut-off
    interaction_const_t ic = setupInteractionConst(options);

    t_nrnb nrnb = { 0 };

    gmx_enerdata_t enerd(1, 0);

    gmx::StepWorkload stepWork;
    stepWork.computeForces = true;
    if (options.computeVirialAndEnergy)
    {
        stepWork.computeVirial = true;
        stepWork.computeEnergy = true;
    }

    const gmx::EnumerationArray<BenchMarkKernels, std::string> kernelNames = { "auto", "no", "4xM",
                                                                               "2xMM" };

    const gmx::EnumerationArray<BenchMarkCombRule, std::string> combruleNames = { "geom.", "LB",
                                                                                  "none" };

    if (!doWarmup)
    {
        fprintf(stdout, "%-7s %-4s %-5s %-4s ",
                options.coulombType == BenchMarkCoulomb::Pme ? "Ewald" : "RF",
                options.useHalfLJOptimization ? "half" : "all",
                combruleNames[options.ljCombinationRule].c_str(), kernelNames[options.nbnxmSimd].c_str());
    }

    // Run pre-iteration to avoid cache misses
    for (int iter = 0; iter < options.numPreIterations; iter++)
    {
        nbv->dispatchNonbondedKernel(gmx::InteractionLocality::Local, ic, stepWork, enbvClearFYes,
                                     system.forceRec, &enerd, &nrnb);
    }

    const int numIterations = (doWarmup ? options.numWarmupIterations : options.numIterations);
    const PairlistSet& pairlistSet = nbv->pairlistSets().pairlistSet(gmx::InteractionLocality::Local);
    const gmx::index numPairs = pairlistSet.natpair_ljq_ + pairlistSet.natpair_lj_ + pairlistSet.natpair_q_;
    gmx_cycles_t cycles = gmx_cycles_read();
    for (int iter = 0; iter < numIterations; iter++)
    {
        // Run the kernel without force clearing
        nbv->dispatchNonbondedKernel(gmx::InteractionLocality::Local, ic, stepWork, enbvClearFNo,
                                     system.forceRec, &enerd, &nrnb);
    }
    cycles = gmx_cycles_read() - cycles;
    if (!doWarmup)
    {
        const double dCycles = static_cast<double>(cycles);
        if (options.cyclesPerPair)
        {
            fprintf(stdout, "%10.3f %10.4f %8.4f %8.4f\n", cycles * 1e-6,
                    dCycles / options.numIterations * 1e-6, dCycles / (options.numIterations * numPairs),
                    dCycles / (options.numIterations * numUsefulPairs));
        }
        else
        {
            fprintf(stdout, "%10.3f %10.4f %8.4f %8.4f\n", dCycles * 1e-6,
                    dCycles / options.numIterations * 1e-6, options.numIterations * numPairs / dCycles,
                    options.numIterations * numUsefulPairs / dCycles);
        }
    }
}

void bench(const int sizeFactor, const KernelBenchOptions& options)
{
    // We don't want to call gmx_omp_nthreads_init(), so we init what we need
    gmx_omp_nthreads_set(emntPairsearch, options.numThreads);
    gmx_omp_nthreads_set(emntNonbonded, options.numThreads);

    const gmx::BenchmarkSystem system(sizeFactor);

    real minBoxSize = norm(system.box[XX]);
    for (int dim = YY; dim < DIM; dim++)
    {
        minBoxSize = std::min(minBoxSize, norm(system.box[dim]));
    }
    if (options.pairlistCutoff > 0.5 * minBoxSize)
    {
        gmx_fatal(FARGS, "The cut-off should be shorter than half the box size");
    }

    std::vector<KernelBenchOptions> optionsList;
    if (options.doAll)
    {
        KernelBenchOptions                        opt = options;
        gmx::EnumerationWrapper<BenchMarkCoulomb> coulombIter;
        for (auto coulombType : coulombIter)
        {
            opt.coulombType = coulombType;
            for (int halfLJ = 0; halfLJ <= 1; halfLJ++)
            {
                opt.useHalfLJOptimization = (halfLJ == 1);

                gmx::EnumerationWrapper<BenchMarkCombRule> combRuleIter;
                for (auto combRule : combRuleIter)
                {
                    if (combRule == BenchMarkCombRule::RuleNone)
                    {
                        continue;
                    }
                    opt.ljCombinationRule = combRule;

                    expandSimdOptionAndPushBack(opt, &optionsList);
                }
            }
        }
    }
    else
    {
        expandSimdOptionAndPushBack(options, &optionsList);
    }
    GMX_RELEASE_ASSERT(!optionsList.empty(), "Expect at least on benchmark setup");

#if GMX_SIMD
    if (options.nbnxmSimd != BenchMarkKernels::SimdNo)
    {
        fprintf(stdout, "SIMD width:           %d\n", GMX_SIMD_REAL_WIDTH);
    }
#endif
    fprintf(stdout, "System size:          %zu atoms\n", system.coordinates.size());
    fprintf(stdout, "Cut-off radius:       %g nm\n", options.pairlistCutoff);
    fprintf(stdout, "Number of threads:    %d\n", options.numThreads);
    fprintf(stdout, "Number of iterations: %d\n", options.numIterations);
    fprintf(stdout, "Compute energies:     %s\n", options.computeVirialAndEnergy ? "yes" : "no");
    if (options.coulombType != BenchMarkCoulomb::ReactionField)
    {
        fprintf(stdout, "Ewald excl. corr.:    %s\n",
                options.nbnxmSimd == BenchMarkKernels::SimdNo || options.useTabulatedEwaldCorr
                        ? "table"
                        : "analytical");
    }
    printf("\n");

    if (options.numWarmupIterations > 0)
    {
        setupAndRunInstance(system, optionsList[0], true);
    }

    fprintf(stdout, "Coulomb LJ   comb. SIMD    Mcycles  Mcycles/it.   %s\n",
            options.cyclesPerPair ? "cycles/pair" : "pairs/cycle");
    fprintf(stdout, "                                                total    useful\n");

    for (const auto& optionsInstance : optionsList)
    {
        setupAndRunInstance(system, optionsInstance, false);
    }
}

} // namespace Nbnxm
