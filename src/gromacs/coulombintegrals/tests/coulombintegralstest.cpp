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
#include <config.h>

#include <math.h>

#include <gtest/gtest.h>

#include "gromacs/coulombintegrals/coulombintegrals.h"

#include "testutils/cmdlinetest.h"
#include "testutils/refdata.h"
#include "testutils/testasserts.h"
#include "testutils/testfilemanager.h"

namespace gmx
{

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

/*! \brief Utility to do the real testing
 *
 * \param[in] cd    Charge distribution
 * \param[in] irow  Row number for atom i
 * \param[in] jrow  Row number for atom j
 * \param[in] xi    Distribution width atom i (may be 0)
 * \param[in] xj    Distribution width atom j (may be 0)
 * \param[in] checker The checker data structure
 */
void testCoulomb(ChargeDistribution          cd,
                 int                         irow,
                 int                         jrow,
                 double                      izeta,
                 double                      jzeta,
                 gmx::test::TestReferenceChecker *checker)
{
    std::vector<double> coulomb;
    std::vector<double> force;
    std::vector<double> ncoulomb;
    std::vector<double> nforce;
    for(int i = 0; i <= 5; i++)
    {
        double r = 0.2*i;
        
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
    }
    const char *name = getChargeDistributionName(cd);
    char buf[256];
    snprintf(buf, sizeof(buf), "Potential-%s", name);
    checker->checkSequence(coulomb.begin(), coulomb.end(), buf);
    snprintf(buf, sizeof(buf), "Force-%s", name);
    checker->checkSequence(force.begin(), force.end(), buf);
    snprintf(buf, sizeof(buf), "NuclearPotential-%s", name);
    checker->checkSequence(ncoulomb.begin(), ncoulomb.end(), buf);
    snprintf(buf, sizeof(buf), "NuclearForce-%s", name);
    checker->checkSequence(nforce.begin(), nforce.end(), buf);
}

class SlaterTest : public ::testing::TestWithParam<std::tuple<std::tuple<int, int>, std::tuple<double, double> > >
{
    protected:
        int                        irow_;
        int                        jrow_;
        double                     xi_;
        double                     xj_;
        test::TestReferenceData    refData_;
        test::TestReferenceChecker checker_;
        
        SlaterTest () : checker_(refData_.rootChecker())
        {
            double toler     = 1e-3;
            auto   tolerance = gmx::test::relativeToleranceAsFloatingPoint(1.0, toler);
            checker_.setDefaultTolerance(tolerance);
            irow_ = std::get<0>(std::get<0>(GetParam()));
            jrow_ = std::get<1>(std::get<0>(GetParam()));
            xi_   = std::get<0>(std::get<1>(GetParam()));
            xj_   = std::get<1>(std::get<1>(GetParam()));
        }
        void runTest()
        {
            testCoulomb(ecdSlater, irow_, jrow_, xi_, xj_, &checker_);
        }
};

class GaussianTest : public ::testing::TestWithParam<std::tuple<double, double> >
{
    protected:
        double                     xi_;
        double                     xj_;
        test::TestReferenceData    refData_;
        test::TestReferenceChecker checker_;
        
        GaussianTest () : checker_(refData_.rootChecker())
        {
            double toler     = 1e-5;
            auto   tolerance = gmx::test::relativeToleranceAsFloatingPoint(1.0, toler);
            checker_.setDefaultTolerance(tolerance);
            xi_   = std::get<0>(GetParam());
            xj_   = std::get<1>(GetParam());
        }
        void runTest()
        {
            testCoulomb(ecdGaussian, 0, 0, xi_, xj_, &checker_);
        }

};

TEST_P (GaussianTest, All)
{
    runTest();
}

TEST_P (SlaterTest, All)
{
    runTest();
}

//! Rows for Slater tests
std::vector<std::tuple<int, int>> c_rows =
{
    { 1, 1 }, { 1, 2 }, { 2, 2 },
    { 1, 3 }, { 2, 3 }, { 3, 3 }
};

//! xi and xj for tests
std::vector<std::tuple<double, double> > c_xi =
{
    {  5.6,   5.7  },
    {  5.7,   5.71 },
    {  5.91,  5.9  },
    { 15.8,  16.0  },
    {  6.1,   6.6  },
    { 22.3,  22.4  },
    { 34.6,  34.5  }
};

INSTANTIATE_TEST_CASE_P(Xi, SlaterTest, ::testing::Combine(::testing::ValuesIn(c_rows), ::testing::ValuesIn(c_xi)));

INSTANTIATE_TEST_CASE_P(Xi, GaussianTest, ::testing::ValuesIn(c_xi));

//! integer xi and xj for tests
std::vector<std::tuple<double, double> > c_xiInteger =
{
    {  3.0,  4.0 },
    { 17.0, 18.0 },
    { 25.0, 26.0 },
    { 29.0, 28.0 },
    { 30.0, 29.0 },
    { 31.0, 33.0 },
    { 37.0, 38.0 }
};

INSTANTIATE_TEST_CASE_P(IntegerXi, SlaterTest, ::testing::Combine(::testing::ValuesIn(c_rows), ::testing::ValuesIn(c_xiInteger)));

} // namespace

} // namespace gmx

