/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2015, by the GROMACS development team, led by
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
 * Tests various corners of tables implementations.
 *
 * The main point of these tests is to check that different functions
 * using table routines compile and work as intended.
 *
 * \author David van der Spoel <david.vanderspoel@icm.uu.se>
 * \ingroup module_tables
 */
#include "gmxpre.h"

#include "gromacs/tables/longrangetables.h"

#include <vector>

#include <gtest/gtest.h>

#include "gromacs/legacyheaders/names.h"
#include "gromacs/legacyheaders/types/enums.h"
#include "gromacs/math/calculate-ewald-splitting-coefficient.h"
#include "gromacs/tables/splinetable.h"
#include "gromacs/utility/snprintf.h"

#include "testutils/refdata.h"
#include "testutils/testasserts.h"
#include "testutils/testfilemanager.h"

namespace gmx
{

namespace
{

class TableTest : public ::testing::Test
{
    protected:
        test::TestReferenceData        refData_;
        test::TestReferenceChecker     checker_;
        TableTest( )
            : checker_(refData_.rootChecker())
        {
            checker_.setDefaultTolerance(test::relativeToleranceAsFloatingPoint(1, 1e-3));
        }

        void testTable(interaction_potential_function pot, const double *param)
        {
            std::vector<double> ewald;
            double              r, e, f;
            for (unsigned int i = 0; (i < 20); i++)
            {
                r = i*0.005;
                (*pot)(param, r, &e, &f);
                ewald.push_back(e);
            }
            for (unsigned int i = 1; (i <= 20); i++)
            {
                r = i*0.1;
                (*pot)(param, r, &e, &f);
                ewald.push_back(e);
            }
            checker_.checkSequenceArray(ewald.size(), &ewald[0], "ewald");
        }
};

TEST_F(TableTest, QEwald)
{
    const double beta = 3.0;
    testTable(v_q_ewald_lr, &beta);
}

TEST_F(TableTest, LJEwald)
{
    const double beta = 3.0;
    testTable(v_lj_ewald_lr, &beta);
}

class TabscaleTest : public ::testing::Test
{
    protected:
        test::TestReferenceData        refData_;
        test::TestReferenceChecker     checker_;
        TabscaleTest( )
            : checker_(refData_.rootChecker())
        {
            checker_.setDefaultTolerance(test::relativeToleranceAsFloatingPoint(1, 1e-3));
        }

        void testTabscale(double ewaldcoeff_q,
                          double ewaldcoeff_lj)
        {
            for (int eeltype = 0; (eeltype < eelNR); eeltype++)
            {
                for (int vdwtype = 0; (vdwtype < evdwNR); vdwtype++)
                {
                    if (EEL_PME_EWALD(eeltype) ||
                        EVDW_PME(vdwtype))
                    {
                        for (int rc = 8; (rc <= 16); rc += 2)
                        {
                            double rcoulomb = 0.1*rc;
                            for (int rv = 8; (rv <= 16); rv += 2)
                            {
                                double                     rvdw = 0.1*rv;
                                char                       buf[256];
                                LongRangeInteractionTables it(eeltype,
                                                              ewaldcoeff_q,
                                                              rcoulomb,
                                                              &ewaldcoeff_q,
                                                              v_q_ewald_lr,
                                                              NULL,
                                                              NULL,
                                                              vdwtype,
                                                              ewaldcoeff_lj,
                                                              rvdw,
                                                              &ewaldcoeff_lj,
                                                              v_lj_ewald_lr,
                                                              NULL,
                                                              NULL,
                                                              0,
                                                              0.0,
                                                              1.0);
                                snprintf(buf, sizeof(buf),
                                         "tabSpacing_Coulomb-%s-%g_LJ-%s-%g",
                                         eel_names[eeltype],
                                         rcoulomb,
                                         evdw_names[vdwtype],
                                         rvdw);
                                checker_.checkDouble(it.tabSpacing(), buf);
                            }
                        }
                    }
                }
            }
        }
};

TEST_F(TabscaleTest, SameEwaldCoeff)
{
    double ewaldcoeff_q  = 2;
    double ewaldcoeff_lj = 2;
    testTabscale(ewaldcoeff_q, ewaldcoeff_lj);
}

TEST_F(TabscaleTest, DifferentEwaldCoeff)
{
    double ewaldcoeff_q  = 2;
    double ewaldcoeff_lj = 3;
    testTabscale(ewaldcoeff_q, ewaldcoeff_lj);
}

//! Short cut for testing vector of double.
static void my_check(test::TestReferenceChecker  checker,
                     const std::vector<double>  &v,
                     const char                 *label)
{
    std::vector<double>::const_iterator vend = v.begin()+20;
    if (vend > v.end())
    {
        vend = v.end();
    }
    checker.checkSequence(v.begin(), vend, label);
}

class ArrayTest : public ::testing::Test
{
    protected:
        test::TestReferenceData        refData_;
        test::TestReferenceChecker     checker_;
        ArrayTest( )
            : checker_(refData_.rootChecker())
        {
            checker_.setDefaultTolerance(test::relativeToleranceAsFloatingPoint(1, 1e-3));
        }

        void testArray(int    eeltype,
                       double q_tol,
                       double rcoulomb,
                       int    vdwtype,
                       double lj_tol,
                       double rvdw)
        {
            double                     ewaldcoeff_q  = calc_ewaldcoeff_q(rcoulomb, q_tol);
            double                     ewaldcoeff_lj = calc_ewaldcoeff_q(rvdw, lj_tol);
            LongRangeInteractionTables it(eeltype,
                                          ewaldcoeff_q,
                                          rcoulomb,
                                          &ewaldcoeff_q,
                                          v_q_ewald_lr,
                                          NULL,
                                          NULL,
                                          vdwtype,
                                          ewaldcoeff_lj,
                                          rvdw,
                                          &ewaldcoeff_lj,
                                          v_lj_ewald_lr,
                                          NULL,
                                          NULL,
                                          0,
                                          0,
                                          1);
            const SplineTable *coulomb = it.coulomb();
            my_check(checker_, coulomb->F(), "coulombF");
            my_check(checker_, coulomb->V(), "coulombV");
            my_check(checker_, coulomb->FDV0(), "coulombFDV0");
            const SplineTable *vdw = it.vanderWaals();
            my_check(checker_, vdw->F(), "vdwF");
            my_check(checker_, vdw->V(), "vdwV");
            my_check(checker_, vdw->FDV0(), "vdwFDV0");
        }
};

TEST_F(ArrayTest, All)
{
    double q_tol    = 1e-3;
    double lj_tol   = 1e-2;
    double rcoulomb = 0.3;
    double rvdw     = 0.4;
    testArray(eelPME, q_tol, rcoulomb, evdwPME, lj_tol, rvdw);
}

} // namespace

} // namespace gmx
