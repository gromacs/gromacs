/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2016,2017, by the GROMACS development team, led by
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

#include "gromacs/math/units.h"
#include "gromacs/pbcutil/ishift.h"
#include "gromacs/pbcutil/pbc.h"

#include "testutils/refdata.h"
#include "testutils/testasserts.h"
#include "testutils/testfilemanager.h"

namespace gmx
{
namespace
{

//! Number of atoms used in these tests.
#define NATOMS 4

class BondedTest : public ::testing::Test
{
    protected:
        rvec   x[NATOMS];
        matrix box;
        test::TestReferenceData           refData_;
        test::TestReferenceChecker        checker_;
        BondedTest( ) :
            checker_(refData_.rootChecker())
        {
            test::FloatingPointTolerance tolerance(test::relativeToleranceAsFloatingPoint(1.0, 1e-6));
            checker_.setDefaultTolerance(tolerance);
            clear_rvecs(NATOMS, x);
            x[1][2] = 1;
            x[2][1] = x[2][2] = 1;
            x[3][0] = x[3][1] = x[3][2] = 1;

            clear_mat(box);
            box[0][0] = box[1][1] = box[2][2] = 1.5;
        }

        void testBondAngle(int epbc)
        {
            rvec  r_ij, r_kj;
            real  cosine_angle, angle;
            int   t1, t2;
            t_pbc pbc;

            set_pbc(&pbc, epbc, box);
            angle = bond_angle(x[0], x[1], x[2], &pbc,
                               r_ij, r_kj, &cosine_angle,
                               &t1, &t2);
            checker_.checkReal(angle, "angle");
            checker_.checkReal(cosine_angle, "cosine_angle");
            checker_.checkInteger(t1, "t1");
            checker_.checkInteger(t2, "t2");
        }

        void testDihedralAngle(int epbc)
        {
            rvec  r_ij, r_kj, r_kl, m, n;
            real  angle;
            int   t1, t2, t3;
            t_pbc pbc;

            set_pbc(&pbc, epbc, box);
            angle = dih_angle(x[0], x[1], x[2], x[3], &pbc,
                              r_ij, r_kj, r_kl, m, n,
                              &t1, &t2, &t3);

            checker_.checkReal(angle, "angle");
            checker_.checkInteger(t1, "t1");
            checker_.checkInteger(t2, "t2");
            checker_.checkInteger(t3, "t3");
        }

        void testIfunc(int                         ftype,
                       const std::vector<t_iatom> &iatoms,
                       const t_iparams             iparams[],
                       int                         epbc)
        {
            real  lambda    = 0;
            real  dvdlambda = 0;
            rvec4 f[NATOMS];
            for (int i = 0; i < NATOMS; i++)
            {
                for (int j = 0; j < 4; j++)
                {
                    f[i][j] = 0;
                }
            }
            rvec  fshift[N_IVEC];
            clear_rvecs(N_IVEC, fshift);
            t_pbc pbc;
            set_pbc(&pbc, epbc, box);
            int   ddgatindex = 0;
            real  energy     = interaction_function[ftype].ifunc(iatoms.size(),
                                                                 iatoms.data(),
                                                                 iparams,
                                                                 x, f, fshift,
                                                                 &pbc,
                                                                 /* const struct t_graph *g */ nullptr,
                                                                 lambda, &dvdlambda,
                                                                 /* const struct t_mdatoms *md */ nullptr,
                                                                 /* struct t_fcdata *fcd */ nullptr,
                                                                 &ddgatindex
#ifdef BUILD_WITH_FDA
                                                                 , nullptr
#endif
                                                                );
            checker_.checkReal(energy, interaction_function[ftype].longname);
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

TEST_F (BondedTest, DihedralAnglePbcXyz)
{
    testDihedralAngle(epbcXYZ);
}

TEST_F (BondedTest, IfuncBondsPbcNo)
{
    std::vector<t_iatom> iatoms = { 0, 0, 1, 0, 1, 2, 0, 2, 3 };
    t_iparams            iparams;
    iparams.harmonic.rA  = iparams.harmonic.rB  = 0.8;
    iparams.harmonic.krA = iparams.harmonic.krB = 50;
    testIfunc(F_BONDS, iatoms, &iparams, epbcNONE);
}

TEST_F (BondedTest, IfuncBondsPbcXy)
{
    std::vector<t_iatom> iatoms = { 0, 0, 1, 0, 1, 2, 0, 2, 3 };
    t_iparams            iparams;
    iparams.harmonic.rA  = iparams.harmonic.rB  = 0.8;
    iparams.harmonic.krA = iparams.harmonic.krB = 50;
    testIfunc(F_BONDS, iatoms, &iparams, epbcXY);
}

TEST_F (BondedTest, IfuncBondsPbcXyz)
{
    std::vector<t_iatom> iatoms = { 0, 0, 1, 0, 1, 2, 0, 2, 3 };
    t_iparams            iparams;
    iparams.harmonic.rA  = iparams.harmonic.rB  = 0.8;
    iparams.harmonic.krA = iparams.harmonic.krB = 50;
    testIfunc(F_BONDS, iatoms, &iparams, epbcXYZ);
}

TEST_F (BondedTest, IfuncAnglesPbcNo)
{
    std::vector<t_iatom> iatoms = { 0, 0, 1, 2, 0, 1, 2, 3 };
    t_iparams            iparams;
    real                 k = 50;
    iparams.harmonic.rA  = iparams.harmonic.rB  = 100;
    iparams.harmonic.krA = iparams.harmonic.krB = k;
    testIfunc(F_ANGLES, iatoms, &iparams, epbcNONE);
}

TEST_F (BondedTest, IfuncAnglesPbcXy)
{
    std::vector<t_iatom> iatoms  = { 0, 0, 1, 2, 0, 1, 2, 3 };
    t_iparams            iparams;
    real                 k = 50;
    iparams.harmonic.rA  = iparams.harmonic.rB  = 100;
    iparams.harmonic.krA = iparams.harmonic.krB = k;
    testIfunc(F_ANGLES, iatoms, &iparams, epbcXY);
}

TEST_F (BondedTest, IfuncAnglesPbcXYZ)
{
    std::vector<t_iatom> iatoms = { 0, 0, 1, 2, 0, 1, 2, 3 };
    t_iparams            iparams;
    real                 k = 50;
    iparams.harmonic.rA  = iparams.harmonic.rB  = 100;
    iparams.harmonic.krA = iparams.harmonic.krB = k;
    testIfunc(F_ANGLES, iatoms, &iparams, epbcXYZ);
}

TEST_F (BondedTest, IfuncProperDihedralsPbcNo)
{
    std::vector<t_iatom> iatoms = { 0, 0, 1, 2, 3 };
    t_iparams            iparams;
    iparams.pdihs.phiA = iparams.pdihs.phiB = -100;
    iparams.pdihs.cpA  = iparams.pdihs.cpB  = 10;
    iparams.pdihs.mult = 1;
    testIfunc(F_PDIHS, iatoms, &iparams, epbcNONE);
}

TEST_F (BondedTest, IfuncProperDihedralsPbcXy)
{
    std::vector<t_iatom> iatoms = { 0, 0, 1, 2, 3 };
    t_iparams            iparams;
    iparams.pdihs.phiA = iparams.pdihs.phiB = -100;
    iparams.pdihs.cpA  = iparams.pdihs.cpB  = 10;
    iparams.pdihs.mult = 1;
    testIfunc(F_PDIHS, iatoms, &iparams, epbcXY);
}

TEST_F (BondedTest, IfuncProperDihedralsPbcXyz)
{
    std::vector<t_iatom> iatoms  = { 0, 0, 1, 2, 3 };
    t_iparams            iparams;
    iparams.pdihs.phiA = iparams.pdihs.phiB = -100;
    iparams.pdihs.cpA  = iparams.pdihs.cpB  = 10;
    iparams.pdihs.mult = 1;
    testIfunc(F_PDIHS, iatoms, &iparams, epbcXYZ);
}

}

}
