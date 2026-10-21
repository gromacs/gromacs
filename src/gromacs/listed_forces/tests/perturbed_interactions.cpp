/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2026- The GROMACS Authors
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
 * Tests for perturbed interactions with ListedForces interaction selections
 *
 * With multiple time stepping, the listed interactions are split
 * over several ListedForces objects, each computing a subset of
 * the interaction groups. These tests check that, for free energy
 * calculations with perturbed interactions, the sum over complementary
 * selections gives the same forces, energies, dH/dlambda and
 * foreign-lambda energies as a single object computing all listed
 * interactions.
 *
 * \ingroup module_listed_forces
 */
#include "gmxpre.h"

#include <cstdint>

#include <array>
#include <ostream>
#include <string>
#include <vector>

#include <gtest/gtest.h>

#include "gromacs/gmxlib/nrnb.h"
#include "gromacs/listed_forces/listed_forces.h"
#include "gromacs/math/paddedvector.h"
#include "gromacs/mdtypes/enerdata.h"
#include "gromacs/mdtypes/fcdata.h"
#include "gromacs/mdtypes/forceoutput.h"
#include "gromacs/mdtypes/forcerec.h"
#include "gromacs/mdtypes/inputrec.h"
#include "gromacs/mdtypes/md_enums.h"
#include "gromacs/mdtypes/simulation_workload.h"
#include "gromacs/pbcutil/ishift.h"
#include "gromacs/pbcutil/pbc.h"
#include "gromacs/topology/forcefieldparameters.h"
#include "gromacs/topology/idef.h"
#include "gromacs/topology/ifunc.h"
#include "gromacs/topology/topsort.h"
#include "gromacs/utility/arrayref.h"
#include "gromacs/utility/enumerationhelpers.h"
#include "gromacs/utility/real.h"
#include "gromacs/utility/vectypes.h"

#include "testutils/testasserts.h"

namespace gmx
{
namespace test
{
namespace
{

//! Shorthand for the interaction groups of ListedForces
using InteractionGroup = ListedForces::InteractionGroup;

//! The number of atoms in the test system, a linear chain
constexpr int c_numAtoms = 6;

//! The current lambda value for all components
constexpr real c_lambda = 0.3;

//! The foreign lambda values for all components
const std::vector<double> c_foreignLambdas = { 0.0, 0.5, 1.0 };

/*! \brief Force-field parameters and FE-sorted interaction definitions for a linear chain
 *
 * The chain has bonds, angles and proper dihedrals. Each type has both
 * perturbed and non-perturbed interactions, interleaved so the free-energy
 * sorting moves interactions.
 */
class TestSystem
{
public:
    TestSystem() : idef_(ffparams_)
    {
        t_iparams params = {};

        // Bonds: type 0 non-perturbed, type 1 perturbed
        params.harmonic = { 0.15, 4.0e5, 0.15, 4.0e5 };
        addType(InteractionFunction::Bonds, params);
        params.harmonic = { 0.14, 3.0e5, 0.16, 2.0e5 };
        addType(InteractionFunction::Bonds, params);

        // Angles: type 2 non-perturbed, type 3 perturbed
        params.harmonic = { 109.5, 400.0, 109.5, 400.0 };
        addType(InteractionFunction::Angles, params);
        params.harmonic = { 105.0, 350.0, 115.0, 450.0 };
        addType(InteractionFunction::Angles, params);

        // Proper dihedrals: type 4 non-perturbed, type 5 perturbed
        params.pdihs = { 0.0, 5.0, 3, 0.0, 5.0 };
        addType(InteractionFunction::ProperDihedrals, params);
        params.pdihs = { 30.0, 8.0, 2, 60.0, 2.0 };
        addType(InteractionFunction::ProperDihedrals, params);

        // Interleave non-perturbed and perturbed interactions
        for (int i = 0; i + 1 < c_numAtoms; i++)
        {
            idef_.il[InteractionFunction::Bonds].push_back(i % 2, std::array<int, 2>{ i, i + 1 });
        }
        for (int i = 0; i + 2 < c_numAtoms; i++)
        {
            idef_.il[InteractionFunction::Angles].push_back(2 + i % 2,
                                                            std::array<int, 3>{ i, i + 1, i + 2 });
        }
        for (int i = 0; i + 3 < c_numAtoms; i++)
        {
            idef_.il[InteractionFunction::ProperDihedrals].push_back(
                    4 + (i + 1) % 2, std::array<int, 4>{ i, i + 1, i + 2, i + 3 });
        }

        const std::vector<int32_t> atomInfo(c_numAtoms, 0);
        gmx_sort_ilist_fe(&idef_, atomInfo);
    }

    //! Returns the FE-sorted interaction definitions
    const InteractionDefinitions& idef() const { return idef_; }
    //! Returns the force-field parameters
    const gmx_ffparams_t& ffparams() const { return ffparams_; }

private:
    //! Adds an interaction type with parameters \p params
    void addType(InteractionFunction ftype, const t_iparams& params)
    {
        ffparams_.functype.push_back(ftype);
        ffparams_.iparams.push_back(params);
    }

    //! The force-field parameters, idef_ refers to these
    gmx_ffparams_t ffparams_;
    //! The interaction definitions
    InteractionDefinitions idef_;
};

//! Output of the listed-forces calculation
struct ListedOutput
{
    //! Constructor
    ListedOutput(const EnumerationArray<FreeEnergyPerturbationCouplingType, std::vector<double>>* allLambdas) :
        force(c_numAtoms, { 0.0_real, 0.0_real, 0.0_real }),
        shiftForce(c_numShiftVectors, { 0.0_real, 0.0_real, 0.0_real }),
        enerd(1, allLambdas)
    {
    }

    //! The forces
    PaddedVector<RVec> force;
    //! The shift forces
    std::vector<RVec> shiftForce;
    //! The energies, dH/dlambda and foreign-lambda terms
    gmx_enerdata_t enerd;
};

/*! \brief Computes all listed forces in \p output using one ListedForces object
 * per selection in \p selections
 *
 * The ListedForces objects are set up as mdrun does with multiple time stepping.
 */
void computeListedForces(const TestSystem&                                  system,
                         ArrayRef<const ListedForces::InteractionSelection> selections,
                         const t_lambda&                                    fepvals,
                         ListedOutput*                                      output)
{
    t_forcerec fr;
    fr.pbcType = PbcType::No;
    fr.bMolPBC = false;

    t_nrnb nrnb;

    std::vector<ListedForces> listedForces;
    for (const auto& selection : selections)
    {
        listedForces.emplace_back(system.ffparams(), 1, 1, 1, selection, nullptr, nullptr, nullptr);
    }

    // Coordinates of a twisted chain
    const PaddedVector<RVec> x   = { { 1.382, 1.573, 1.482 }, { 1.281, 1.559, 1.596 },
                                     { 1.292, 1.422, 1.663 }, { 1.189, 1.407, 1.775 },
                                     { 1.221, 1.302, 1.881 }, { 1.105, 1.243, 1.947 } };
    matrix                   box = { { 3.0, 0.0, 0.0 }, { 0.0, 3.0, 0.0 }, { 0.0, 0.0, 3.0 } };

    t_disresdata disres = {};
    t_fcdata     fcdata;
    fcdata.disres = &disres;

    const std::vector<real> lambda(static_cast<int>(FreeEnergyPerturbationCouplingType::Count), c_lambda);
    const std::vector<real>           charges(c_numAtoms, 0.0_real);
    const std::vector<unsigned short> cENER(c_numAtoms, 0);

    StepWorkload stepWork;
    stepWork.computeForces       = true;
    stepWork.computeListedForces = true;
    stepWork.computeVirial       = true;
    stepWork.computeEnergy       = true;
    stepWork.computeDhdl         = true;

    ForceWithShiftForces forceWithShiftForces(output->force.arrayRefWithPadding(), true, output->shiftForce);
    ForceOutputs forceOutputs(forceWithShiftForces, false, ForceWithVirial(output->force, true));

    for (auto& lf : listedForces)
    {
        lf.setup(system.idef(), c_numAtoms, false, {});
        lf.calculate(nullptr,
                     box,
                     x.constArrayRefWithPadding(),
                     {},
                     &fcdata,
                     nullptr,
                     &forceOutputs,
                     &fr,
                     nullptr,
                     &output->enerd,
                     &nrnb,
                     lambda,
                     charges,
                     charges,
                     {},
                     cENER,
                     0,
                     nullptr,
                     stepWork);
    }

    output->enerd.foreignLambdaTerms.finalizePotentialContributions(output->enerd.dvdl_lin, lambda, fepvals);
}

//! Parameters: the groups assigned to the second (slow) level, as with mts-level2-forces
struct SecondLevelGroups
{
    //! The interaction groups in the second level
    std::vector<InteractionGroup> groups;
};

//! Returns a name for the groups in \p secondLevelGroups, e.g. "PairsDihedralsAngles"
std::string secondLevelGroupsName(const SecondLevelGroups& secondLevelGroups)
{
    std::string name;
    for (const auto group : secondLevelGroups.groups)
    {
        switch (group)
        {
            case InteractionGroup::Pairs: name += "Pairs"; break;
            case InteractionGroup::Dihedrals: name += "Dihedrals"; break;
            case InteractionGroup::Angles: name += "Angles"; break;
            default: name += "Rest"; break;
        }
    }
    return name;
}

//! Prints \p secondLevelGroups readably in GoogleTest messages
void PrintTo(const SecondLevelGroups& secondLevelGroups, std::ostream* os)
{
    *os << secondLevelGroupsName(secondLevelGroups);
}

//! Returns the test name for the second-level groups in \p info
std::string secondLevelGroupsTestName(const ::testing::TestParamInfo<SecondLevelGroups>& info)
{
    return secondLevelGroupsName(info.param);
}

//! Test fixture, parametrized by the groups in the second level
class PerturbedInteractionSelectionTest : public ::testing::TestWithParam<SecondLevelGroups>
{
};

TEST_P(PerturbedInteractionSelectionTest, SplitMatchesAll)
{
    TestSystem system;
    ASSERT_EQ(system.idef().ilsort, ilsortFE_SORTED);

    t_lambda fepvals;
    fepvals.n_lambda = gmx::ssize(c_foreignLambdas);
    for (auto& foreignLambdas : fepvals.all_lambda)
    {
        foreignLambdas = c_foreignLambdas;
    }

    // As mdrun with multiple time stepping: the first level gets
    // the groups not in the second level and the rest group
    const SecondLevelGroups&           secondLevelGroups = GetParam();
    ListedForces::InteractionSelection secondLevel;
    for (const auto group : secondLevelGroups.groups)
    {
        secondLevel.set(static_cast<int>(group));
    }
    ListedForces::InteractionSelection firstLevel = ~secondLevel;
    ASSERT_TRUE(firstLevel.test(static_cast<int>(InteractionGroup::Rest)));

    const std::vector<ListedForces::InteractionSelection> selectionsAll = {
        ListedForces::interactionSelectionAll()
    };
    const std::vector<ListedForces::InteractionSelection> selectionsSplit = { firstLevel, secondLevel };

    ListedOutput outputAll(&fepvals.all_lambda);
    computeListedForces(system, selectionsAll, fepvals, &outputAll);
    ListedOutput outputSplit(&fepvals.all_lambda);
    computeListedForces(system, selectionsSplit, fepvals, &outputSplit);

    const auto tolerance = relativeToleranceAsFloatingPoint(100.0, GMX_REAL_EPS * 100);

    for (int a = 0; a < c_numAtoms; a++)
    {
        for (int d = 0; d < DIM; d++)
        {
            EXPECT_REAL_EQ_TOL(outputAll.force[a][d], outputSplit.force[a][d], tolerance)
                    << "force on atom " << a << " dim " << d;
        }
    }
    for (int s = 0; s < c_numShiftVectors; s++)
    {
        for (int d = 0; d < DIM; d++)
        {
            EXPECT_REAL_EQ_TOL(outputAll.shiftForce[s][d], outputSplit.shiftForce[s][d], tolerance)
                    << "shift force " << s << " dim " << d;
        }
    }
    for (const auto ftype :
         { InteractionFunction::Bonds, InteractionFunction::Angles, InteractionFunction::ProperDihedrals })
    {
        EXPECT_REAL_EQ_TOL(outputAll.enerd.term[ftype], outputSplit.enerd.term[ftype], tolerance)
                << "energy of " << interaction_function[ftype].longname;
    }
    // The free-energy terms are stored in double, but are sums of terms computed in real
    EXPECT_REAL_EQ_TOL(outputAll.enerd.dvdl_nonlin[FreeEnergyPerturbationCouplingType::Bonded],
                       outputSplit.enerd.dvdl_nonlin[FreeEnergyPerturbationCouplingType::Bonded],
                       tolerance);

    // The foreign-lambda energies only contain the perturbed interactions
    const auto energiesAll   = outputAll.enerd.foreignLambdaTerms.energies();
    const auto energiesSplit = outputSplit.enerd.foreignLambdaTerms.energies();
    ASSERT_EQ(energiesAll.size(), energiesSplit.size());
    for (Index i = 0; i < energiesAll.ssize(); i++)
    {
        EXPECT_REAL_EQ_TOL(energiesAll[i], energiesSplit[i], tolerance)
                << "foreign-lambda energy index " << i;
    }
}

INSTANTIATE_TEST_SUITE_P(WithMts,
                         PerturbedInteractionSelectionTest,
                         ::testing::Values(SecondLevelGroups{ { InteractionGroup::Dihedrals } },
                                           SecondLevelGroups{ { InteractionGroup::Angles } },
                                           SecondLevelGroups{ { InteractionGroup::Pairs,
                                                                InteractionGroup::Dihedrals,
                                                                InteractionGroup::Angles } }),
                         secondLevelGroupsTestName);

} // namespace
} // namespace test
} // namespace gmx
