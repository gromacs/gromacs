/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2015,2016,2017, by the GROMACS development team, led by
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
 * Tests for functionality of the "angle" trajectory analysis module.
 *
 * \author David van der Spoel <david.vanderspoel@icm.uu.se>
 * \ingroup module_applied_forces
 */
#include "gmxpre.h"

#include "gromacs/applied-forces/electricfield.h"

#include <gtest/gtest.h>

#include "gromacs/gmxlib/network.h"
#include "gromacs/math/vec.h"
#include "gromacs/mdlib/forcerec.h"
#include "gromacs/mdtypes/forceoutput.h"
#include "gromacs/mdtypes/iforceprovider.h"
#include "gromacs/mdtypes/imdmodule.h"
#include "gromacs/mdtypes/imdpoptionprovider.h"
#include "gromacs/mdtypes/inputrec.h"
#include "gromacs/mdtypes/mdatom.h"
#include "gromacs/options/options.h"
#include "gromacs/options/treesupport.h"
#include "gromacs/utility/arrayref.h"
#include "gromacs/utility/keyvaluetreebuilder.h"
#include "gromacs/utility/keyvaluetreetransform.h"
#include "gromacs/utility/real.h"
#include "gromacs/utility/smalloc.h"
#include "gromacs/utility/stringcompare.h"
#include "gromacs/utility/stringutil.h"

#include "testutils/testasserts.h"

namespace
{

/********************************************************************
 * ElectricFieldTest
 */

class ElectricFieldTest : public ::testing::Test
{
    public:
        void test(int  dim,
                  real E0,
                  real omega,
                  real t0,
                  real sigma,
                  real expectedValue)
        {
            gmx::test::FloatingPointTolerance tolerance(
                    gmx::test::relativeToleranceAsFloatingPoint(1.0, 0.005));
            auto                              module(gmx::createElectricFieldModule());

            // Prepare MDP inputs
            const char *dimXYZ[3] = { "x", "y", "z" };
            GMX_RELEASE_ASSERT(dim >= 0 && dim < DIM, "Dimension should be 0, 1 or 2");

            gmx::KeyValueTreeBuilder     mdpValues;
            mdpValues.rootObject().addValue(gmx::formatString("E%s", dimXYZ[dim]),
                                            gmx::formatString("1 %g 0", E0));
            mdpValues.rootObject().addValue(gmx::formatString("E%s-t", dimXYZ[dim]),
                                            gmx::formatString("3 %g 0 %g 0 %g 0", omega, t0, sigma));

            gmx::KeyValueTreeTransformer transform;
            transform.rules()->addRule()
                .keyMatchType("/", gmx::StringCompareType::CaseAndDashInsensitive);
            module->mdpOptionProvider()->initMdpTransform(transform.rules());
            auto         result = transform.transform(mdpValues.build(), nullptr);
            gmx::Options moduleOptions;
            module->mdpOptionProvider()->initMdpOptions(&moduleOptions);
            gmx::assignOptionsFromKeyValueTree(&moduleOptions, result.object(), nullptr);

            ForceProviders forceProviders;
            module->initForceProviders(&forceProviders);

            t_mdatoms            md;
            PaddedRVecVector     f = { { 0, 0, 0 } };
            gmx::ForceWithVirial forceWithVirial(f, true);
            md.homenr = 1;
            snew(md.chargeA, md.homenr);
            md.chargeA[0] = 1;

            t_commrec  *cr = init_commrec();
            forceProviders.calculateForces(cr, &md, nullptr, 0, nullptr, &forceWithVirial);
            done_commrec(cr);

            EXPECT_REAL_EQ_TOL(f[0][dim], expectedValue, tolerance);
            sfree(md.chargeA);
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
