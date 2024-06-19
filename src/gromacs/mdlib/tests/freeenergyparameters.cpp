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
 *
 * \brief Tests routines in freeenergyparameters.h .
 *
 * \author Christian Blau <blau@kth.se>
 */

#include "gmxpre.h"

#include "gromacs/mdlib/freeenergyparameters.h"

#include <cstdint>

#include <algorithm>
#include <string>
#include <vector>

#include <gmock/gmock.h>
#include <gtest/gtest.h>

#include "gromacs/mdtypes/inputrec.h"
#include "gromacs/utility/arrayref.h"
#include "gromacs/utility/enumerationhelpers.h"
#include "gromacs/utility/real.h"

#include "testutils/testasserts.h"
#include "testutils/testmatchers.h"

enum class FreeEnergyPerturbationCouplingType : int;


namespace gmx
{
namespace test
{
namespace
{


/*! \brief Parameters that will vary from test to test.
 */
struct FreeEnergyParameterTestParameters
{
    //! current state of lambda in the simulation, -1 if not set
    int currentLambdaState = -1;
    //! Fractional value of lambda to start from, -1 if not set
    double initLambda = -1;
    //! The initial number of the state, -1 if not set
    int initFepState = -1;
    //! Change of lambda per time step (fraction of (0.1)
    double deltaLambda = 0;
    //! number of lambda entries
    int nLambda = 0;
    //! the current simulation step
    int64_t step = 0;
    //! the expected lambda at the current simulation step
    gmx::EnumerationArray<FreeEnergyPerturbationCouplingType, real> expectedLambdas = { -1, -1, -1,
                                                                                        -1, -1, -1,
                                                                                        -1 };
};


/*! \brief Sets of parameters on which to run the tests.
 */
const FreeEnergyParameterTestParameters freeEnergyParameterSets[] = {
    // no parameters set at all
    { -1, -1, -1, 0, 1, 0, { -1, -1, -1, -1, -1, -1, -1 } },
    // setting current lambda state to 0, no other variables set, using eftpNR * [0.8] as lambda state matrix
    { 0, -1, -1, 0, 1, 1, { 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8 } },
    // setting current lambda state to 1, no other variables set, using eftpNR * [0.2,0.8] as lambda state matrix
    { 1, -1, -1, 0, 2, 1, { 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8 } },
    // test that non-zero deltaLambda trumps current lambda state
    // setting current lambda state to 1, using deltaLambda 0.1 and setp 0 and eftpNR * [0.2,0.8] as lambda state matrix
    { 1, -1, -1, 0.1, 2, 0, { 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2 } },
    // test that interpolating lambda values starting
    // from lambda = 0, deltaLambda = 0.1, step = 10 results in values at end-state
    { 1, 0, -1, 0.1, 2, 10, { 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8 } },
    // interpolating half the way
    { 1, 0, -1, 0.1, 2, 5, { 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5 } },
    // starting from end-state 1 and move lambda half-way with negative deltaLambda = -0.1
    { -1, -1, 1, -0.1, 2, 5, { 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5 } },
    // starting from end-state 1 and move one step backwards with negative deltaLambda = -0.1
    { -1, -1, 1, -0.1, 2, 1, { 0.74, 0.74, 0.74, 0.74, 0.74, 0.74, 0.74 } },
    // three lambda states, the last two equal ([0.2,0.8,0.8]), move forward with deltaLambda = 0.1
    { -1, -1, 0, 0.1, 3, 0, { 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2 } },
    { -1, -1, 0, 0.1, 3, 3, { 0.56, 0.56, 0.56, 0.56, 0.56, 0.56, 0.56 } },
    { -1, -1, 0, 0.1, 3, 7, { 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8 } },
    { -1, -1, 0, 0.1, 3, 8, { 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8 } },
    { -1, -1, 0, 0.1, 3, 10, { 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8 } },
    // three lambda states, the last two equal ([0.2,0.8,0.8]), move backwards
    { -1, -1, 2, -0.1, 3, 1, { 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8 } },
    { -1, -1, 2, -0.1, 3, 2, { 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8 } },
    { -1, -1, 2, -0.1, 3, 3, { 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8 } },
    { -1, -1, 2, -0.1, 3, 5, { 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8 } },
    { -1, -1, 2, -0.1, 3, 7, { 0.56, 0.56, 0.56, 0.56, 0.56, 0.56, 0.56 } },
    { -1, -1, 2, -0.1, 3, 8, { 0.44, 0.44, 0.44, 0.44, 0.44, 0.44, 0.44 } },
    { -1, -1, 2, -0.1, 3, 10, { 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2 } },
    // three lambda states, the last two equal ([0.2,0.8,0.8]), start in middle state, move backwards
    { -1, -1, 1, -0.1, 3, 0, { 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8 } },
    { -1, -1, 1, -0.1, 3, 2, { 0.56, 0.56, 0.56, 0.56, 0.56, 0.56, 0.56 } },
    { -1, -1, 1, -0.1, 3, 3, { 0.44, 0.44, 0.44, 0.44, 0.44, 0.44, 0.44 } },
};

/*! \brief Test for setting free energy parameters.
 */
class FreeEnergyParameterTest : public ::testing::TestWithParam<FreeEnergyParameterTestParameters>
{
public:
    //! Fills in the required FEP values for testing from the test parameters
    t_lambda getFepVals()
    {
        t_lambda fepvals;
        fepvals.init_fep_state             = GetParam().initFepState;
        fepvals.init_lambda_without_states = GetParam().initLambda;
        fepvals.delta_lambda               = GetParam().deltaLambda;
        std::fill(fepvals.all_lambda.begin(),
                  fepvals.all_lambda.end(),
                  defaultLambdaArrayForTest_[GetParam().nLambda]);
        fepvals.n_lambda = GetParam().nLambda;
        return fepvals;
    }

private:
    //! a set of default lambda arrays for different lengths
    std::vector<std::vector<double>> defaultLambdaArrayForTest_ = { {}, { 0.8 }, { 0.2, 0.8 }, { 0.2, 0.8, 0.8 } };
};

TEST_P(FreeEnergyParameterTest, CorrectLambdas)
{
    EXPECT_THAT(GetParam().expectedLambdas,
                Pointwise(RealEq(defaultRealTolerance()),
                          currentLambdas(GetParam().step, getFepVals(), GetParam().currentLambdaState)));
}

INSTANTIATE_TEST_SUITE_P(WithParameters, FreeEnergyParameterTest, ::testing::ValuesIn(freeEnergyParameterSets));

} // namespace

} // namespace test

} // namespace gmx
