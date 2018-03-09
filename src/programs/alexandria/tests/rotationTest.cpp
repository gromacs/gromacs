/*
 * This source file is part of the Alexandria program.
 *
 * Copyright (C) 2014-2018 
 *
 * Developers:
 *             Mohammad Mehdi Ghahremanpour, 
 *             Paul J. van Maaren, 
 *             David van der Spoel (Project leader)
 *
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, 
 * Boston, MA  02110-1301, USA.
 */
 
/*! \internal \brief
 * Implements part of the alexandria program.
 * \author Mohammad Mehdi Ghahremanpour <mohammad.ghahremanpour@icm.uu.se>
 * \author David van der Spoel <david.vanderspoel@icm.uu.se>
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








