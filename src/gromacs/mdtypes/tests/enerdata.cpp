/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2022- The GROMACS Authors
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
 * Tests for the \p ForeignLambdaTerms class
 *
 * \author berk Hess <hess@kth.se>
 * \ingroup module_mdtypes
 */
#include "gmxpre.h"

#include "gromacs/mdtypes/enerdata.h"

#include <array>
#include <string>
#include <tuple>
#include <vector>

#include <gmock/gmock.h>
#include <gtest/gtest.h>

#include "gromacs/mdtypes/inputrec.h"
#include "gromacs/mdtypes/md_enums.h"
#include "gromacs/utility/arrayref.h"
#include "gromacs/utility/basedefinitions.h"
#include "gromacs/utility/enumerationhelpers.h"
#include "gromacs/utility/real.h"

#include "testutils/testasserts.h"

namespace gmx
{
namespace test
{
namespace
{

constexpr int c_numLambdas = 5;

const std::array<double, c_numLambdas> dhdlVdw  = { 1.0, 2.0, 4.0, 8.0, 16.0 };
const std::array<double, c_numLambdas> dhdlCoul = { 32.0, 64.0, 128.0, 256.0, 512.0 };

const gmx::EnumerationArray<FreeEnergyPerturbationCouplingType, double> dhdlLinearZero = { 0.0, 0.0,
                                                                                           0.0, 0.0,
                                                                                           0.0, 0.0 };

t_lambda makeFepvals(const std::vector<double>& lambdaVdw, const std::vector<double>& lambdaCoul)
{
    t_lambda fepvals;

    fepvals.n_lambda = gmx::ssize(lambdaVdw);

    for (auto couplingType : gmx::EnumerationArray<FreeEnergyPerturbationCouplingType, real>::keys())
    {
        fepvals.all_lambda[couplingType].resize(fepvals.n_lambda);
    }
    fepvals.all_lambda[FreeEnergyPerturbationCouplingType::Vdw]  = lambdaVdw;
    fepvals.all_lambda[FreeEnergyPerturbationCouplingType::Coul] = lambdaCoul;

    for (auto& sep : fepvals.separate_dvdl)
    {
        sep = false;
    }

    return fepvals;
}

} // namespace

// Check that the rate check catches a setup with different rates
TEST(ForeingLambdaTermsDhdl, RateCheckWorks)
{
    const std::vector<double> lamVdw{ 0.0, 0.5 };
    const std::vector<double> lamCoul{ 0.0, 0.25 };

    const t_lambda fepvals = makeFepvals(lamVdw, lamCoul);

    ASSERT_EQ(fepLambdasChangeAtSameRate(fepvals.all_lambda), false);
}

TEST(ForeingLambdaTermsDhdl, AllLinear)
{
    const std::vector<double> lamSet{ 0.0, 0.25, 0.5, 0.75, 1.0 };

    const t_lambda fepvals = makeFepvals(lamSet, lamSet);

    ASSERT_EQ(fepLambdasChangeAtSameRate(fepvals.all_lambda), true);

    ForeignLambdaTerms foreignLambdaTerms(&fepvals.all_lambda);

    for (int i = 0; i < c_numLambdas; i++)
    {
        foreignLambdaTerms.accumulate(1 + i, FreeEnergyPerturbationCouplingType::Vdw, 0.0, dhdlVdw[i]);
        foreignLambdaTerms.accumulate(1 + i, FreeEnergyPerturbationCouplingType::Coul, 0.0, dhdlCoul[i]);
    }

    foreignLambdaTerms.finalizePotentialContributions(dhdlLinearZero, {}, fepvals);

    std::vector<double> dummy;
    std::vector<double> dhdl;
    std::tie(dummy, dhdl) = foreignLambdaTerms.getTerms(nullptr);

    for (int i = 0; i < c_numLambdas; i++)
    {
        EXPECT_FLOAT_EQ(dhdl[i], dhdlVdw[i] + dhdlCoul[i]);
    }
}

TEST(ForeingLambdaTermsDhdl, AllLinearNegative)
{
    const std::vector<double> lamSet{ 1.0, 0.75, 0.5, 0.25, 0.0 };

    const t_lambda fepvals = makeFepvals(lamSet, lamSet);

    EXPECT_EQ(fepLambdasChangeAtSameRate(fepvals.all_lambda), true);

    ForeignLambdaTerms foreignLambdaTerms(&fepvals.all_lambda);

    for (int i = 0; i < c_numLambdas; i++)
    {
        foreignLambdaTerms.accumulate(1 + i, FreeEnergyPerturbationCouplingType::Vdw, 0.0, dhdlVdw[i]);
        foreignLambdaTerms.accumulate(1 + i, FreeEnergyPerturbationCouplingType::Coul, 0.0, dhdlCoul[i]);
    }

    foreignLambdaTerms.finalizePotentialContributions(dhdlLinearZero, {}, fepvals);

    std::vector<double> dummy;
    std::vector<double> dhdl;
    std::tie(dummy, dhdl) = foreignLambdaTerms.getTerms(nullptr);

    for (int i = 0; i < c_numLambdas; i++)
    {
        EXPECT_FLOAT_EQ(dhdl[i], -dhdlVdw[i] - dhdlCoul[i]);
    }
}

TEST(ForeingLambdaTermsDhdl, SeparateVdwCoul)
{
    const std::vector<double> lamVdw{ 0.0, 0.5, 1.0, 1.0, 1.0 };
    const std::vector<double> lamCoul{ 0.0, 0.0, 0.0, 0.5, 1.0 };

    const t_lambda fepvals = makeFepvals(lamVdw, lamCoul);

    ASSERT_EQ(fepLambdasChangeAtSameRate(fepvals.all_lambda), true);

    ForeignLambdaTerms foreignLambdaTerms(&fepvals.all_lambda);

    for (int i = 0; i < c_numLambdas; i++)
    {
        foreignLambdaTerms.accumulate(1 + i, FreeEnergyPerturbationCouplingType::Vdw, 0.0, dhdlVdw[i]);
        foreignLambdaTerms.accumulate(1 + i, FreeEnergyPerturbationCouplingType::Coul, 0.0, dhdlCoul[i]);
    }

    foreignLambdaTerms.finalizePotentialContributions(dhdlLinearZero, {}, fepvals);

    std::vector<double> dummy;
    std::vector<double> dhdl;
    std::tie(dummy, dhdl) = foreignLambdaTerms.getTerms(nullptr);

    for (int i = 0; i < 2; i++)
    {
        EXPECT_FLOAT_EQ(dhdl[i], dhdlVdw[i]);
    }
    EXPECT_FLOAT_EQ(dhdl[2], (dhdlVdw[2] + dhdlCoul[2]) * 0.5);
    for (int i = 3; i < 5; i++)
    {
        EXPECT_FLOAT_EQ(dhdl[i], dhdlCoul[i]);
    }
}

} // namespace test
} // namespace gmx
