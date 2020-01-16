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
 * Implements nblib kernel setup options
 *
 * \author Berk Hess <hess@kth.se>
 * \author Victor Holanda <victor.holanda@cscs.ch>
 * \author Joe Jordan <ejjordan@kth.se>
 * \author Prashanth Kanduri <kanduri@cscs.ch>
 * \author Sebastian Keller <keller@cscs.ch>
 */
#include "gmxpre.h"

#include "nbkerneloptions.h"

#include "gromacs/compat/optional.h"
#include "gromacs/ewald/ewald_utils.h"
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
#include "gromacs/pbcutil/mshift.h"
#include "gromacs/pbcutil/pbc.h"
#include "gromacs/timing/cyclecounter.h"
#include "gromacs/topology/topology.h"
#include "gromacs/utility/enumerationhelpers.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/logger.h"

namespace nblib
{

/*! \brief Checks the kernel setup
 *
 * Returns an error string when the kernel is not available.
 */
static gmx::compat::optional<std::string> checkKernelSetup(const NBKernelOptions &options)
{
    GMX_RELEASE_ASSERT(options.nbnxmSimd < BenchMarkKernels::Count &&
                       options.nbnxmSimd != BenchMarkKernels::SimdAuto, "Need a valid kernel SIMD type");

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
static Nbnxm::KernelType translateBenchmarkEnum(const BenchMarkKernels &kernel)
{
    int kernelInt = static_cast<int>(kernel);
    return static_cast<Nbnxm::KernelType>(kernelInt);
}

/*! \brief Returns the kernel setup
 */
static Nbnxm::KernelSetup getKernelSetup(const NBKernelOptions &options)
{
    auto messageWhenInvalid = checkKernelSetup(options);
    GMX_RELEASE_ASSERT(!messageWhenInvalid, "Need valid options");

    Nbnxm::KernelSetup kernelSetup;

    //The int enum options.nbnxnSimd is set up to match Nbnxm::KernelType + 1
    kernelSetup.kernelType         = translateBenchmarkEnum(options.nbnxmSimd);
    // The plain-C kernel does not support analytical ewald correction
    if (kernelSetup.kernelType == Nbnxm::KernelType::Cpu4x4_PlainC)
    {
        kernelSetup.ewaldExclusionType = Nbnxm::EwaldExclusionType::Table;
    }
    else
    {
        kernelSetup.ewaldExclusionType = options.useTabulatedEwaldCorr ? Nbnxm::EwaldExclusionType::Table : Nbnxm::EwaldExclusionType::Analytical;
    }

    return kernelSetup;
}

static real ewaldCoeff(const real ewald_rtol, const real pairlistCutoff)
{
    return calc_ewaldcoeff_q(pairlistCutoff, ewald_rtol);
}

//! Return an interaction constants struct with members used in the benchmark set appropriately
interaction_const_t setupInteractionConst(const NBKernelOptions &options)
{
    interaction_const_t ic;

    ic.vdwtype          = evdwCUT;
    ic.vdw_modifier     = eintmodPOTSHIFT;
    ic.rvdw             = options.pairlistCutoff;

    ic.eeltype          = (options.coulombType == BenchMarkCoulomb::Pme ? eelPME : eelRF);
    ic.coulomb_modifier = eintmodPOTSHIFT;
    ic.rcoulomb         = options.pairlistCutoff;

    // Reaction-field with epsilon_rf=inf
    // TODO: Replace by calc_rffac() after refactoring that
    ic.k_rf             = 0.5*std::pow(ic.rcoulomb, -3);
    ic.c_rf             = 1/ic.rcoulomb + ic.k_rf*ic.rcoulomb*ic.rcoulomb;

    if (EEL_PME_EWALD(ic.eeltype))
    {
        // Ewald coefficients, we ignore the potential shift
        ic.ewaldcoeff_q = ewaldCoeff(1e-5, options.pairlistCutoff);
        GMX_RELEASE_ASSERT(ic.ewaldcoeff_q > 0, "Ewald coefficient should be > 0");
        ic.coulombEwaldTables = std::make_unique<EwaldCorrectionTables>();
        init_interaction_const_tables(nullptr, &ic);
    }

    return ic;
}

//! Sets up and returns a Nbnxm object for the given options and system
std::unique_ptr<nonbonded_verlet_t>
setupNbnxmInstance(const NBKernelOptions   &options,
                   NBKernelSystem          &system)
{
    const auto         pinPolicy       = (options.useGpu ? gmx::PinningPolicy::PinnedIfSupported : gmx::PinningPolicy::CannotBePinned);
    const int          numThreads      = options.numThreads;
    // Note: the options and Nbnxm combination rule enums values should match
    const int          combinationRule = static_cast<int>(options.ljCombinationRule);

    auto               messageWhenInvalid = checkKernelSetup(options);
    if (messageWhenInvalid)
    {
        gmx_fatal(FARGS, "Requested kernel is unavailable because %s.",
                  messageWhenInvalid->c_str());
    }
    Nbnxm::KernelSetup        kernelSetup = getKernelSetup(options);

    PairlistParams            pairlistParams(kernelSetup.kernelType, false, options.pairlistCutoff, false);

    Nbnxm::GridSet            gridSet(epbcXYZ, false, nullptr, nullptr, pairlistParams.pairlistType, false, numThreads, pinPolicy);

    auto                      pairlistSets = std::make_unique<PairlistSets>(pairlistParams, false, 0);

    auto                      pairSearch   = std::make_unique<PairSearch>(epbcXYZ, false, nullptr, nullptr,
                                                                          pairlistParams.pairlistType,
                                                                          false, numThreads, pinPolicy);

    auto atomData     = std::make_unique<nbnxn_atomdata_t>(pinPolicy);

    // Put everything together
    auto nbv = std::make_unique<nonbonded_verlet_t>(std::move(pairlistSets),
                                                    std::move(pairSearch),
                                                    std::move(atomData),
                                                    kernelSetup,
                                                    nullptr,
                                                    nullptr);

    nbnxn_atomdata_init(gmx::MDLogger(),
                        nbv->nbat.get(), kernelSetup.kernelType,
                        combinationRule, system.numAtoms, system.nonbondedParameters,
                        1, numThreads);

    t_nrnb nrnb;

    GMX_RELEASE_ASSERT(!TRICLINIC(system.box), "Only rectangular unit-cells are supported here");
    const rvec               lowerCorner = { 0, 0, 0 };
    const rvec               upperCorner = {
        system.box[XX][XX],
        system.box[YY][YY],
        system.box[ZZ][ZZ]
    };

    gmx::ArrayRef<const int> atomInfo;
    atomInfo = system.atomInfoAllVdw;

    const real atomDensity = system.coordinates.size()/det(system.box);

    nbnxn_put_on_grid(nbv.get(),
                      system.box, 0, lowerCorner, upperCorner,
                      nullptr, {0, int(system.coordinates.size())}, atomDensity,
                      atomInfo, system.coordinates,
                      0, nullptr);

    nbv->constructPairlist(gmx::InteractionLocality::Local,
                           &system.excls, 0, &nrnb);

    t_mdatoms mdatoms;
    // We only use (read) the atom type and charge from mdatoms
    mdatoms.typeA   = const_cast<int *>(system.atomTypes.data());
    mdatoms.chargeA = const_cast<real *>(system.charges.data());
    nbv->setAtomProperties(mdatoms, atomInfo);

    return nbv;
}

//! Add the options instance to the list for all requested kernel SIMD types
//! TODO This should be refactored so that if SimdAuto is set only one kernel
//!      layout is chosen.
//! TODO This should be refactored to only return the desired kernel layout
static void
expandSimdOptionAndPushBack(const NBKernelOptions        &options,
                            std::vector<NBKernelOptions> *optionsList)
{
    if (options.nbnxmSimd == BenchMarkKernels::SimdAuto)
    {
        bool addedInstance = false;
#ifdef GMX_NBNXN_SIMD_4XN
        optionsList->push_back(options);
        optionsList->back().nbnxmSimd = BenchMarkKernels::Simd4XM;
        addedInstance = true;
#endif
#ifdef GMX_NBNXN_SIMD_2XNN
        optionsList->push_back(options);
        optionsList->back().nbnxmSimd = BenchMarkKernels::Simd2XMM;
        addedInstance = true;
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



void nbKernel(NBKernelSystem        &system,
              const NBKernelOptions &options,
              const bool            &printTimings)
{
    // We don't want to call gmx_omp_nthreads_init(), so we init what we need
    gmx_omp_nthreads_set(emntPairsearch, options.numThreads);
    gmx_omp_nthreads_set(emntNonbonded, options.numThreads);

    real                       minBoxSize = norm(system.box[XX]);
    for (int dim = YY; dim < DIM; dim++)
    {
        minBoxSize = std::min(minBoxSize, norm(system.box[dim]));
    }
    if (options.pairlistCutoff > 0.5*minBoxSize)
    {
        gmx_fatal(FARGS, "The cut-off should be shorter than half the box size");
    }

    std::vector<NBKernelOptions> optionsList;
    expandSimdOptionAndPushBack(options, &optionsList);
    GMX_RELEASE_ASSERT(!optionsList.empty(), "Expect at least one benchmark setup");

    // setupAndRunInstance(system, optionsList[0], printTimings);
}

} // namespace nblib
