/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2015- The GROMACS Authors
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
 * Tests for functionality of the electric field module.
 *
 * \author David van der Spoel <david.vanderspoel@icm.uu.se>
 * \ingroup module_applied_forces
 */
#include "gmxpre.h"

#include "gromacs/applied_forces/electricfield.h"

#include <memory>
#include <string>
#include <vector>

#include <gtest/gtest.h>

#include "gromacs/gmxlib/network.h"
#include "gromacs/math/paddedvector.h"
#include "gromacs/math/vectypes.h"
#include "gromacs/mdlib/forcerec.h"
#include "gromacs/mdtypes/commrec.h"
#include "gromacs/mdtypes/enerdata.h"
#include "gromacs/mdtypes/forceoutput.h"
#include "gromacs/mdtypes/iforceprovider.h"
#include "gromacs/mdtypes/imdmodule.h"
#include "gromacs/mdtypes/imdpoptionprovider.h"
#include "gromacs/mdtypes/inputrec.h"
#include "gromacs/mdtypes/mdatom.h"
#include "gromacs/options/options.h"
#include "gromacs/options/treesupport.h"
#include "gromacs/utility/arrayref.h"
#include "gromacs/utility/basedefinitions.h"
#include "gromacs/utility/gmxassert.h"
#include "gromacs/utility/keyvaluetreebuilder.h"
#include "gromacs/utility/keyvaluetreetransform.h"
#include "gromacs/utility/real.h"
#include "gromacs/utility/smalloc.h"
#include "gromacs/utility/stringcompare.h"
#include "gromacs/utility/stringutil.h"

#include "testutils/testasserts.h"

namespace gmx
{
namespace test
{
namespace
{

/********************************************************************
 * ElectricFieldTest
 */

class ElectricFieldTest : public ::testing::Test
{
public:
    static void test(int dim, real E0, real omega, real t0, real sigma, real expectedValue)
    {
        // Make the electric field module
        auto module = createElectricFieldModule();

        // Fill the module as if from .mdp inputs
        {
            const char* dimXYZ[3] = { "x", "y", "z" };
            GMX_RELEASE_ASSERT(dim >= 0 && dim < DIM, "Dimension should be 0, 1 or 2");

            KeyValueTreeBuilder mdpValues;
            mdpValues.rootObject().addValue(formatString("electric-field-%s", dimXYZ[dim]),
                                            formatString("%g %g %g %g", E0, omega, t0, sigma));

            KeyValueTreeTransformer transform;
            transform.rules()->addRule().keyMatchType("/", StringCompareType::CaseAndDashInsensitive);
            module->mdpOptionProvider()->initMdpTransform(transform.rules());
            auto    result = transform.transform(mdpValues.build(), nullptr);
            Options moduleOptions;
            module->mdpOptionProvider()->initMdpOptions(&moduleOptions);
            assignOptionsFromKeyValueTree(&moduleOptions, result.object(), nullptr);
        }

        // Prepare a ForceProviderInput
        std::vector<real> chargeA{ 1 };
        t_commrec         cr;
        matrix            boxDummy = { { 0, 0, 0 }, { 0, 0, 0 }, { 0, 0, 0 } };
        ForceProviderInput forceProviderInput({}, gmx::ssize(chargeA), chargeA, {}, 0.0, 0, boxDummy, cr);

        // Prepare a ForceProviderOutput
        PaddedVector<RVec>  f = { { 0, 0, 0 } };
        ForceWithVirial     forceWithVirial(f, true);
        gmx_enerdata_t      enerdDummy(1, nullptr);
        ForceProviderOutput forceProviderOutput(&forceWithVirial, &enerdDummy);

        // Use the ForceProviders to calculate forces
        ForceProviders forceProviders;
        module->initForceProviders(&forceProviders);
        forceProviders.calculateForces(forceProviderInput, &forceProviderOutput);

        FloatingPointTolerance tolerance(relativeToleranceAsFloatingPoint(1.0, 0.005));
        EXPECT_REAL_EQ_TOL(f[0][dim], expectedValue, tolerance);
    }
};

TEST_F(ElectricFieldTest, Static)
{
    test(0, 1, 0, 0, 0, 96.4853363);
}

TEST_F(ElectricFieldTest, Oscillating)
{
    test(0, 1, 5, 0.2, 0, 96.4853363);
}

TEST_F(ElectricFieldTest, Pulsed)
{
    test(0, 1, 5, 0.5, 1, -68.215782);
}

} // namespace
} // namespace test
} // namespace gmx
