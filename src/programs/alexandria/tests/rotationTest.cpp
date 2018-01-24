/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2017, by the GROMACS development team, led by
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
 * Implements test of regression code
 *
 * \author David van der Spoel <david.vanderspoel@icm.uu.se>
 * \ingroup module_alexandria
 */
 
#include <cmath>
#include <cstdlib>
#include <gtest/gtest.h>

#include "testutils/cmdlinetest.h"
#include "testutils/refdata.h"
#include "testutils/testasserts.h"
#include "testutils/testfilemanager.h"

#include "gromacs/math/functions.h"

#include "programs/alexandria/regression.h"

class RotationTest : public gmx::test::CommandLineTestBase
{
    protected:
        gmx::test::TestReferenceChecker checker_;

        //init set tolecrance
        RotationTest () : checker_(this->rootChecker())
        {
            auto tolerance = gmx::test::relativeToleranceAsFloatingPoint(1.0, 5e-2);
            checker_.setDefaultTolerance(tolerance);
        }

    double RMSD(tensor a, tensor b)
        {
            double rmsd = 0;
            for (int i = 0; i < DIM; i++)
            {
                for (int j = 0; j < DIM; j++)
                {
                    rmsd += gmx::square(a[i][j] - b[i][j]);
                }
            } 
            rmsd /= (DIM*DIM);
            return rmsd;
        }
    
    void testRotation(tensor p, tensor q)
        {
            double rmsd;
            tensor rotated_p;
            kabsch_rotation(p, q, rotated_p);
            rmsd = RMSD(rotated_p, q);
            checker_.checkDouble(rmsd, "rmsd");
        }
};

TEST_F (RotationTest, Rotate_p_to_q)
{
    tensor p = {{0, 0, 0}, {7, 3, 5}, {0, 0, 0}};
    tensor q = {{0, 0, 0}, {5, 3, 7}, {0, 0, 0}};
    testRotation(p, q);
}








