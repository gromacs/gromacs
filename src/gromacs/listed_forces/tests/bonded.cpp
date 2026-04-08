/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2018- The GROMACS Authors
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
 * Implements test of bonded force routines
 *
 *
 * \todo We should re-work this. For example, a harmonic bond
 * has so few computations that force components should be
 * accurate to a small *and computed* relative error.
 *
 * \author David van der Spoel <david.vanderspoel@icm.uu.se>
 * \ingroup module_listed_forces
 */
#include "gmxpre.h"

#include "gromacs/listed_forces/bonded.h"

#include <cmath>
#include <cstdint>

#include <algorithm>
#include <iterator>
#include <limits>
#include <memory>
#include <ostream>
#include <string>
#include <tuple>
#include <vector>

#include <gtest/gtest.h>

#include "gromacs/gpu_utils/capabilities.h"
#include "gromacs/listed_forces/listed_forces.h"
#include "gromacs/math/paddedvector.h"
#include "gromacs/mdtypes/md_enums.h"
#include "gromacs/pbcutil/ishift.h"
#include "gromacs/pbcutil/pbc.h"
#include "gromacs/topology/idef.h"
#include "gromacs/topology/ifunc.h"
#include "gromacs/utility/arrayref.h"
#include "gromacs/utility/real.h"
#include "gromacs/utility/strconvert.h"
#include "gromacs/utility/vec.h"
#include "gromacs/utility/vectypes.h"

#include "testutils/refdata.h"
#include "testutils/test_device.h"
#include "testutils/test_hardware_environment.h"
#include "testutils/testasserts.h"

#include "bondedtestdata.h"
#include "bondedtestrunners.h"

namespace gmx
{
namespace test
{
namespace
{

class ListedForcesTest :
    public ::testing::TestWithParam<std::tuple<iListInput, PaddedVector<RVec>, PbcType>>
{
public:
    //! Before any test is run, work out whether any compatible GPUs exist.
    static std::vector<std::unique_ptr<IBondedTestRunner>> getRunners()
    {
        std::vector<std::unique_ptr<IBondedTestRunner>> runners;
        runners.push_back(std::make_unique<BondedHostTestRunner>());
        if (GpuConfigurationCapabilities::Bonded)
        {
            for (const auto& testDevice : getTestHardwareEnvironment()->getTestDeviceList())
            {
                runners.push_back(std::make_unique<BondedDeviceTestRunner>(*testDevice));
            }
        }
        return runners;
    }

protected:
    matrix                 box_;
    t_pbc                  pbc_;
    PaddedVector<RVec>     x_;
    PbcType                pbcType_;
    iListInput             input_;
    TestReferenceData      refData_;
    TestReferenceChecker   checker_;
    FloatingPointTolerance shiftForcesTolerance_ = defaultRealTolerance();
    ListedForcesTest() : checker_(refData_.rootChecker())
    {
        input_   = std::get<0>(GetParam());
        x_       = std::get<1>(GetParam());
        pbcType_ = std::get<2>(GetParam());
        clear_mat(box_);
        box_[0][0] = box_[1][1] = box_[2][2] = 1.5;
        set_pbc(&pbc_, pbcType_, box_);
        // We need quite specific tolerances here since angle functions
        // etc. are not very precise and reproducible.
        test::FloatingPointTolerance tolerance(test::FloatingPointTolerance(
                input_.ftoler, input_.dtoler, 1.0e-6, 1.0e-12, 10000, 100, false));
        checker_.setDefaultTolerance(tolerance);
        // The SIMD acos() is only accurate to 2-3 ULP, so the angles
        // computed by it and the non-SIMD code paths (that use
        // std::acos) differ by enough to require quite large
        // tolerances for the shift forces in mixed precision.
        float singleShiftForcesAbsoluteTolerance =
                ((input_.ftype == InteractionFunction::Polarization)
                                 || (input_.ftype == InteractionFunction::AnharmonicPolarization)
                                 || (IS_ANGLE(input_.ftype.value()))
                         ? 5e-3
                         : 5e-5);
        // Note that std::numeric_limits isn't required by the standard to
        // have an implementation for uint64_t(!) but this is likely to
        // work because that type is likely to be a typedef for one of
        // the other numerical types that happens to be 64-bits wide.
        shiftForcesTolerance_ = FloatingPointTolerance(singleShiftForcesAbsoluteTolerance,
                                                       1e-8,
                                                       1e-6,
                                                       1e-12,
                                                       std::numeric_limits<uint64_t>::max(),
                                                       std::numeric_limits<uint64_t>::max(),
                                                       false);
    }
    void testOneIfunc(TestReferenceChecker* checker, const std::vector<t_iatom>& iatoms, const real lambda)
    {
        SCOPED_TRACE(std::string("Testing PBC type: ") + c_pbcTypeNames[pbcType_]);

        std::vector<BondedKernelFlavor> flavors = { BondedKernelFlavor::ForcesAndVirialAndEnergy };
        if (!input_.fep || lambda == 0)
        {
            flavors.push_back(BondedKernelFlavor::ForcesSimdWhenAvailable);
        }

        for (const auto& runner : getRunners())
        {
            for (const auto flavor : flavors)
            {
                if (!runner->supportsFlavor(flavor, input_, lambda))
                {
                    continue;
                }
                SCOPED_TRACE("Testing on " + runner->hardwareDescription()
                             + " with bonded kernel flavor: " + c_bondedKernelFlavorStrings[flavor]);
                OutputQuantities output;
                runner->run(input_, x_, pbc_, lambda, iatoms, flavor, &output);
                EXPECT_TRUE((input_.fep || (output.dvdlambda == 0.0)))
                        << "dvdlambda was " << output.dvdlambda;
                checkOutput(checker, output, flavor);
                auto shiftForcesChecker = checker->checkCompound("Shift-Forces", "Shift-forces");
                if (computeVirial(flavor))
                {
                    shiftForcesChecker.setDefaultTolerance(shiftForcesTolerance_);
                    shiftForcesChecker.checkVector(output.fshift[c_centralShiftIndex], "Central");
                }
                else
                {
                    shiftForcesChecker.disableUnusedEntriesCheck();
                }
            }
        }
    }
    void testIfunc()
    {
        EXPECT_TRUE(input_.ftype.has_value() && input_.ftype < InteractionFunction::Count);
        TestReferenceChecker thisChecker =
                checker_.checkCompound("FunctionType", interaction_function[input_.ftype.value()].name)
                        .checkCompound("FEP", (input_.fep ? "Yes" : "No"));
        std::vector<t_iatom> iatoms;
        fillIatoms(input_.ftype, &iatoms);
        if (input_.fep)
        {
            const int numLambdas = 3;
            for (int i = 0; i < numLambdas; ++i)
            {
                const real lambda       = i / (numLambdas - 1.0);
                auto       valueChecker = thisChecker.checkCompound("Lambda", toString(lambda));
                testOneIfunc(&valueChecker, iatoms, lambda);
            }
        }
        else
        {
            testOneIfunc(&thisChecker, iatoms, 0.0);
        }
    }
};

TEST_P(ListedForcesTest, Ifunc)
{
    testIfunc();
}

//! Function types for testing bonds. Add new terms at the end.
std::vector<iListInput> c_InputBonds = {
    { iListInput(2e-6F, 1e-8).setHarmonic(InteractionFunction::Bonds, 0.15, 500.0) },
    { iListInput(2e-6F, 1e-8).setHarmonic(InteractionFunction::Bonds, 0.15, 500.0, 0.17, 400.0) },
    { iListInput(1e-4F, 1e-8).setHarmonic(InteractionFunction::GROMOS96Bonds, 0.15, 50.0) },
    { iListInput().setHarmonic(InteractionFunction::GROMOS96Bonds, 0.15, 50.0, 0.17, 40.0) },
    { iListInput().setCubic(0.16, 50.0, 2.0) },
    { iListInput(2e-6F, 1e-8).setMorse(0.15, 50.0, 2.0, 0.17, 40.0, 1.6) },
    { iListInput(2e-6F, 1e-8).setMorse(0.15, 30.0, 2.7) },
    { iListInput().setFene(0.4, 5.0) }
};

//! Constants for Quartic Angles
const real cQuarticAngles[5] = { 1.1, 2.3, 4.6, 7.8, 9.2 };

//! Function types for testing angles. Add new terms at the end.
std::vector<iListInput> c_InputAngles = {
    { iListInput(2e-3, 1e-8).setHarmonic(InteractionFunction::Angles, 100.0, 50.0) },
    { iListInput(2e-3, 1e-8).setHarmonic(InteractionFunction::Angles, 100.15, 50.0, 95.0, 30.0) },
    { iListInput(8e-3, 1e-8).setHarmonic(InteractionFunction::GROMOS96Angles, 100.0, 50.0) },
    { iListInput(8e-3, 1e-8).setHarmonic(InteractionFunction::GROMOS96Angles, 100.0, 50.0, 95.0, 30.0) },
    { iListInput().setLinearAngle(50.0, 0.4) },
    { iListInput().setLinearAngle(50.0, 0.4, 40.0, 0.6) },
    { iListInput(2e-6, 1e-8).setCrossBondBonds(0.8, 0.7, 45.0) },
    { iListInput(3e-6, 1e-8).setCrossBondAngles(0.8, 0.7, 0.3, 45.0) },
    { iListInput(2e-2, 1e-8).setUreyBradley(950.0, 46.0, 0.3, 5.0) },
    { iListInput(2e-2, 1e-8).setUreyBradley(100.0, 45.0, 0.3, 5.0, 90.0, 47.0, 0.32, 7.0) },
    { iListInput(2e-3, 1e-8).setQuarticAngles(87.0, cQuarticAngles) }
};

//! Constants for Ryckaert-Bellemans A
const real rbcA[NR_RBDIHS] = { -5.35, 13.6, 8.4, -16.7, 0.3, 12.4 };

//! Constants for Ryckaert-Bellemans B
const real rbcB[NR_RBDIHS] = { -6.35, 12.6, 8.1, -10.7, 0.9, 15.4 };

//! Constants for Ryckaert-Bellemans without FEP
const real rbc[NR_RBDIHS] = { -7.35, 13.6, 8.4, -16.7, 1.3, 12.4 };

//! Function types for testing dihedrals. Add new terms at the end.
std::vector<iListInput> c_InputDihs = {
    { iListInput(5e-4, 1e-8).setPDihedrals(InteractionFunction::ProperDihedrals, -100.0, 10.0, 2, -80.0, 20.0) },
    { iListInput(1e-4, 1e-8).setPDihedrals(InteractionFunction::ProperDihedrals, -105.0, 15.0, 2) },
    { iListInput(2e-4, 1e-8).setHarmonic(InteractionFunction::ImproperDihedrals, 100.0, 50.0) },
    { iListInput(2e-4, 1e-8).setHarmonic(InteractionFunction::ImproperDihedrals, 100.15, 50.0, 95.0, 30.0) },
    { iListInput(4e-4, 1e-8).setRbDihedrals(rbcA, rbcB) },
    { iListInput(4e-4, 1e-8).setRbDihedrals(rbc) }
};

//! Function types for testing polarization. Add new terms at the end.
std::vector<iListInput> c_InputPols = {
    { iListInput(2e-5, 1e-8).setPolarization(0.12) },
    { iListInput(2e-3, 1e-8).setAnharmPolarization(0.0013, 0.02, 1235.6) },
    { iListInput(1.4e-3, 1e-8).setTholePolarization(0.26, 0.07, 0.09) },
    { iListInput(2e-3, 1e-8).setWaterPolarization(0.001, 0.0012, 0.0016, 0.095, 0.15, 0.02) },
};

//! Function types for testing polarization. Add new terms at the end.
std::vector<iListInput> c_InputRestraints = {
    { iListInput(1e-4, 1e-8).setPDihedrals(InteractionFunction::AngleRestraints, -100.0, 10.0, 2, -80.0, 20.0) },
    { iListInput(1e-4, 1e-8).setPDihedrals(InteractionFunction::AngleRestraints, -105.0, 15.0, 2) },
    { iListInput(1e-4, 1e-8).setPDihedrals(InteractionFunction::AngleZAxisRestraints, -100.0, 10.0, 2, -80.0, 20.0) },
    { iListInput(1e-4, 1e-8).setPDihedrals(InteractionFunction::AngleZAxisRestraints, -105.0, 15.0, 2) },
    { iListInput(2e-3, 1e-8).setHarmonic(InteractionFunction::RestrictedBendingPotential, 100.0, 50.0) },
    { iListInput(2e-3, 1e-8).setHarmonic(InteractionFunction::RestrictedBendingPotential, 100.0, 50.0, 110.0, 45.0) }
};

//! Function types for testing bond with zero length, has zero reference length to make physical sense.
std::vector<iListInput> c_InputBondsZeroLength = {
    { iListInput().setHarmonic(InteractionFunction::Bonds, 0.0, 500.0) },
};

//! Function types for testing angles with zero angle, has zero reference angle to make physical sense.
std::vector<iListInput> c_InputAnglesZeroAngle = {
    { iListInput(2e-3, 1e-8).setHarmonic(InteractionFunction::Angles, 0.0, 50.0) },
};

} // namespace
} // namespace test

//! Print an RVec to \c os
static void PrintTo(const RVec& value, std::ostream* os)
{
    *os << value[XX] << " " << value[YY] << " " << value[ZZ] << std::endl;
}

//! Print a padded vector of RVec to \c os
static void PrintTo(const PaddedVector<RVec>& vector, std::ostream* os)
{
    if (vector.empty())
    {
        *os << "Empty vector" << std::endl;
    }
    else
    {
        *os << "Vector of RVec containing:" << std::endl;
        std::for_each(vector.begin(), vector.end(), [os](const RVec& v) { PrintTo(v, os); });
    }
}

namespace test
{
namespace
{

/*! \brief Coordinates for testing
 *
 * Taken from a butane molecule, so we have some
 * normal-sized bonds and angles to test.
 *
 * \todo Test also some weirder values */
std::vector<PaddedVector<RVec>> c_coordinatesForTests = {
    { { 1.382, 1.573, 1.482 }, { 1.281, 1.559, 1.596 }, { 1.292, 1.422, 1.663 }, { 1.189, 1.407, 1.775 } }
};

//! Coordinates for testing bonds with zero length
std::vector<PaddedVector<RVec>> c_coordinatesForTestsZeroBondLength = {
    { { 0.0, 0.0, 0.0 }, { 0.0, 0.0, 0.0 }, { 0.005, 0.0, 0.1 }, { 0.005, 0.0, 0.1 } }
};

//! Coordinates for testing bonds with zero length
std::vector<PaddedVector<RVec>> c_coordinatesForTestsZeroAngle = {
    { { 0.005, 0.0, 0.1 }, { 0.0, 0.0, 0.0 }, { 0.005, 0.0, 0.1 }, { 0.5, 0.18, 0.22 } }
};

//! PBC values for testing
std::vector<PbcType> c_pbcForTests = { PbcType::No, PbcType::XY, PbcType::Xyz };

INSTANTIATE_TEST_SUITE_P(Bond,
                         ListedForcesTest,
                         ::testing::Combine(::testing::ValuesIn(c_InputBonds),
                                            ::testing::ValuesIn(c_coordinatesForTests),
                                            ::testing::ValuesIn(c_pbcForTests)));

INSTANTIATE_TEST_SUITE_P(Angle,
                         ListedForcesTest,
                         ::testing::Combine(::testing::ValuesIn(c_InputAngles),
                                            ::testing::ValuesIn(c_coordinatesForTests),
                                            ::testing::ValuesIn(c_pbcForTests)));

INSTANTIATE_TEST_SUITE_P(Dihedral,
                         ListedForcesTest,
                         ::testing::Combine(::testing::ValuesIn(c_InputDihs),
                                            ::testing::ValuesIn(c_coordinatesForTests),
                                            ::testing::ValuesIn(c_pbcForTests)));

INSTANTIATE_TEST_SUITE_P(Polarize,
                         ListedForcesTest,
                         ::testing::Combine(::testing::ValuesIn(c_InputPols),
                                            ::testing::ValuesIn(c_coordinatesForTests),
                                            ::testing::ValuesIn(c_pbcForTests)));

INSTANTIATE_TEST_SUITE_P(Restraints,
                         ListedForcesTest,
                         ::testing::Combine(::testing::ValuesIn(c_InputRestraints),
                                            ::testing::ValuesIn(c_coordinatesForTests),
                                            ::testing::ValuesIn(c_pbcForTests)));

INSTANTIATE_TEST_SUITE_P(BondZeroLength,
                         ListedForcesTest,
                         ::testing::Combine(::testing::ValuesIn(c_InputBondsZeroLength),
                                            ::testing::ValuesIn(c_coordinatesForTestsZeroBondLength),
                                            ::testing::ValuesIn(c_pbcForTests)));

INSTANTIATE_TEST_SUITE_P(AngleZero,
                         ListedForcesTest,
                         ::testing::Combine(::testing::ValuesIn(c_InputAnglesZeroAngle),
                                            ::testing::ValuesIn(c_coordinatesForTestsZeroAngle),
                                            ::testing::ValuesIn(c_pbcForTests)));

} // namespace

} // namespace test

} // namespace gmx
