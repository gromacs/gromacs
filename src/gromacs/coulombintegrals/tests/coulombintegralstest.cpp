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
 * Implements test of Coulomb routines
 *
 * \author David van der Spoel <david.vanderspoel@icm.uu.se>
 * \ingroup module_correlationfunctions
 */


#include <math.h>
#include <gtest/gtest.h>

#include "gromacs/coulombintegrals/coulombintegrals.h"

#include "testutils/cmdlinetest.h"
#include "testutils/refdata.h"
#include "testutils/testasserts.h"
#include "testutils/testfilemanager.h"

namespace
{

enum ChargeDistribution  {
    ecdSlater,
    ecdGaussian, 
    ecdNR
};

typedef struct {
    ChargeDistribution cd;
    const char        *name;
} t_chargedistribution;

t_chargedistribution chargedistributions[ecdNR] = {
    { ecdSlater,   "Slater"},
    { ecdGaussian, "Gaussian"}
};

static const char *getChargeDistributionName(ChargeDistribution cd)
{
    for (auto i = 0; i < ecdNR; i++)
    {
        if (cd == chargedistributions[i].cd)
        {
            return chargedistributions[i].name;
        }
    }
    return nullptr;
}

class CoulombTest : public gmx::test::CommandLineTestBase
{
    protected:
        gmx::test::TestReferenceChecker checker_;

        //init set tolecrance
        CoulombTest () : checker_(this->rootChecker())
        {
            #if !HAVE_LIBCLN
            gmx::test::FloatingPointTolerance tolerance = gmx::test::relativeToleranceAsFloatingPoint(1.0, 1e-13);
            checker_.setDefaultTolerance(tolerance);
            #endif
        }

        // Static initiation, only run once every test.
        static void SetUpTestCase()
        {
        }

        void testCoulomb(ChargeDistribution cd, 
                         int                irow, 
                         int                jrow, 
                         double             izeta, 
                         double             jzeta, 
                         double             lambda)
        {
            std::vector<double> coulomb;
            std::vector<double> force;
            std::vector<double> ncoulomb;
            std::vector<double> nforce;
            for(int i = 0; i < 9; i++)
            {
                double r = lambda*(i+1);
                
                switch (cd)
                {
                case ecdGaussian:
                    coulomb.push_back(Coulomb_GG(r,  izeta, jzeta));
                    force.push_back(-DCoulomb_GG(r,  izeta, jzeta));
                    ncoulomb.push_back(Nuclear_GG(r, izeta));
                    nforce.push_back(-DNuclear_GG(r, izeta));
                    break;
                case ecdSlater:
                    coulomb.push_back(Coulomb_SS(r,  irow, jrow, izeta, jzeta));
                    force.push_back(-DCoulomb_SS(r,  irow, jrow, izeta, jzeta));
                    ncoulomb.push_back(Nuclear_SS(r, irow, izeta));
                    nforce.push_back(-DNuclear_SS(r, irow, izeta));
                    break;
                default:
                    break;
                }
                
                if (lambda == 0)
                    break;
            }
            const char *name = getChargeDistributionName(cd);
            char buf[256];
            snprintf(buf, sizeof(buf), "Potential-%s", name);
            checker_.checkSequence(coulomb.begin(), coulomb.end(), buf);
            snprintf(buf, sizeof(buf), "Force-%s", name);
            checker_.checkSequence(force.begin(), force.end(), buf);
            snprintf(buf, sizeof(buf), "NuclearPotential-%s", name);
            checker_.checkSequence(ncoulomb.begin(), ncoulomb.end(), buf);
            snprintf(buf, sizeof(buf), "NuclearForce-%s", name);
            checker_.checkSequence(nforce.begin(), nforce.end(), buf);
        }
        

        static void TearDownTestCase()
        {
        }

};

TEST_F (CoulombTest, Gaussian)
{
    testCoulomb(ecdGaussian, 1, 1, 5.0, 8.0, 0.1);
}

TEST_F (CoulombTest, GaussianSameXi)
{
    testCoulomb(ecdGaussian, 1, 1, 7.0, 7.0, 0.1);
}

TEST_F (CoulombTest, SlaterRow1Row2)
{
    testCoulomb(ecdSlater, 1, 2, 12.0, 16.0, 0.1);
}

TEST_F (CoulombTest, SlaterRow1Row1)
{
    testCoulomb(ecdSlater, 1, 1, 3.0, 16.0, 0.1);
}

TEST_F (CoulombTest, SlaterRow1Row3)
{
    testCoulomb(ecdSlater, 1, 3, 12.0, 17.0, 0.1);
}

TEST_F (CoulombTest, SlaterRow2Row2)
{
    testCoulomb(ecdSlater, 2, 2, 5.0, 16.0, 0.1);
}

TEST_F (CoulombTest, SlaterRow2Row3)
{
    testCoulomb(ecdSlater, 2, 3, 12.0, 9.0, 0.1);
}

TEST_F (CoulombTest, SlaterRow3Row3)
{
    testCoulomb(ecdSlater, 3, 3, 11.0, 8.0, 0.1);
}

TEST_F (CoulombTest, GaussianLargeDistance)
{
    testCoulomb(ecdGaussian, 1, 1, 5.0, 8.0, 2);
}

TEST_F (CoulombTest, SlaterLargeDistance)
{
    testCoulomb(ecdSlater, 3, 3, 11.0, 8.0, 2);
}

TEST_F (CoulombTest, GaussianZeroDistance)
{
    testCoulomb(ecdGaussian, 1, 1, 5.0, 8.0, 0);
}

TEST_F (CoulombTest, SlaterZeroDistance)
{
    testCoulomb(ecdSlater, 3, 3, 11.0, 8.0, 0);
}

}

