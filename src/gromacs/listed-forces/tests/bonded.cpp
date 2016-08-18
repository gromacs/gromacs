/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2016, by the GROMACS development team, led by
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
 * Implements test of bonded force routines
 *
 * \author David van der Spoel <david.vanderspoel@icm.uu.se>
 * \ingroup module_listed-forces
 */
#include "gmxpre.h"

#include "gromacs/listed-forces/bonded.h"

#include <cmath>

#include <memory>

#include <gtest/gtest.h>

#include "gromacs/pbcutil/ishift.h"
#include "gromacs/pbcutil/pbc.h"
#include "gromacs/utility/fatalerror.h"

#include "testutils/testasserts.h"
#include "testutils/testfilemanager.h"

namespace gmx
{
namespace
{

class BondedTest : public ::testing::Test
{
    protected:
        rvec   xi, xj, xk, xl;
        matrix box;
        gmx::test::FloatingPointTolerance tolerance;
        // Use erefdataCreateMissing for creating new files
        BondedTest( ) : tolerance(gmx::test::defaultRealTolerance())
        {
            xi[0] = xi[1] = xi[2] = 0;
            xj[0] = xj[1] = 0; xj[2] = 1;
            xk[0] = 0; xk[1] = xk[2] = 1;
            xl[0] = xl[1] = xl[2] = 1;
            clear_mat(box);
            box[0][0] = box[1][1] = box[2][2] = 1.5;
        }

        // Static initiation, only run once every test.
        static void SetUpTestCase()
        {
        }

        static void TearDownTestCase()
        {
        }

        void testBondAngle(int epbc)
        {
            rvec  r_ij, r_kj;
            real  costh, angle;
            int   t1, t2;
            t_pbc pbc;

            set_pbc(&pbc, epbc, box);
            angle = bond_angle(xi, xj, xk, &pbc,
                               r_ij, r_kj, &costh,
                               &t1, &t2);
            EXPECT_REAL_EQ_TOL(M_PI/2.0, angle, tolerance);
            EXPECT_REAL_EQ_TOL(0, costh, tolerance);
            switch (epbc)
            {
                case epbcNONE:
                    EXPECT_EQ(CENTRAL, t1);
                    EXPECT_EQ(CENTRAL, t2);
                    break;
                case epbcXY:
                    EXPECT_EQ(CENTRAL, t1);
                    EXPECT_EQ(17, t2);
                    break;
                case epbcXYZ:
                    EXPECT_EQ(37, t1);
                    EXPECT_EQ(17, t2);
                    break;
                default:
                    gmx_fatal(FARGS, "Invalid epbc %d", epbc);
            }
        }

        void testDihedralAngle(int epbc)
        {
            rvec  r_ij, r_kj, r_kl, m, n;
            real  costh, angle;
            int   t1, t2, t3;
            t_pbc pbc;

            set_pbc(&pbc, epbc, box);
            angle = dih_angle(xi, xj, xk, xl, &pbc,
                              r_ij, r_kj, r_kl, m, n, &costh,
                              &t1, &t2, &t3);
            switch (epbc)
            {
                case epbcNONE:
                    EXPECT_REAL_EQ_TOL(-M_PI/2.0, angle, tolerance);
                    EXPECT_REAL_EQ_TOL(-1, costh, tolerance);
                    EXPECT_EQ(CENTRAL, t1);
                    EXPECT_EQ(CENTRAL, t2);
                    EXPECT_EQ(CENTRAL, t3);
                    break;
                case epbcXY:
                    EXPECT_REAL_EQ_TOL(-M_PI/2.0, angle, tolerance);
                    EXPECT_REAL_EQ_TOL(-1, costh, tolerance);
                    EXPECT_EQ(CENTRAL, t1);
                    EXPECT_EQ(17, t2);
                    EXPECT_EQ(CENTRAL+1, t3);
                    break;
                case epbcXYZ:
                    EXPECT_REAL_EQ_TOL(M_PI/2.0, angle, tolerance);
                    EXPECT_REAL_EQ_TOL(1, costh, tolerance);
                    EXPECT_EQ(37, t1);
                    EXPECT_EQ(17, t2);
                    EXPECT_EQ(CENTRAL+1, t3);
                    break;
                default:
                    gmx_fatal(FARGS, "Invalid epbc %d", epbc);
            }
        }

};

TEST_F (BondedTest, BondAnglePbcNone)
{
    testBondAngle(epbcNONE);
}

TEST_F (BondedTest, BondAnglePbcXy)
{
    testBondAngle(epbcXY);
}

TEST_F (BondedTest, BondAnglePbcXyz)
{
    testBondAngle(epbcXYZ);
}

TEST_F (BondedTest, DihedralAnglePbcNone)
{
    testDihedralAngle(epbcNONE);
}

TEST_F (BondedTest, DihedralAnglePbcXy)
{
    testDihedralAngle(epbcXY);
}

TEST_F (BondedTest, DihedarlAnglePbcXyz)
{
    testDihedralAngle(epbcXYZ);
}

}

}
