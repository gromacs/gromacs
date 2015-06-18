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

#include "gromacs/tables/tables.h"

#include <vector>

#include <gtest/gtest.h>

#include "gromacs/legacyheaders/names.h"
#include "gromacs/legacyheaders/types/enums.h"
#include "gromacs/math/calculate-ewald-splitting-coefficient.h"
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
            double              dr   = 0.05;

            std::vector<double> ewald;
            for (unsigned int i = 1; (i <= 20); i++)
            {
                ewald.push_back((*pot)(param, i*dr));
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

        void testTabscale(real ewaldcoeff_q,
                          real ewaldcoeff_lj)
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
                            real rcoulomb = 0.1*rc;
                            for (int rv = 8; (rv <= 16); rv += 2)
                            {
                                real              rvdw = 0.1*rv;
                                char              buf[256];
                                InteractionTables it(eeltype,
                                                     ewaldcoeff_q,
                                                     rcoulomb,
                                                     vdwtype,
                                                     ewaldcoeff_lj,
                                                     rvdw);
                                snprintf(buf, sizeof(buf),
                                         "tabScale_Coulomb-%s-%g_LJ-%s-%g",
                                         eel_names[eeltype],
                                         rcoulomb,
                                         evdw_names[vdwtype],
                                         rvdw);
                                checker_.checkDouble(it.tabScale(), buf);
                            }
                        }
                    }
                }
            }
        }
};

TEST_F(TabscaleTest, DifferentEwaldCoeff)
{
    real ewaldcoeff_q  = 1e-4;
    real ewaldcoeff_lj = 1e-5;
    testTabscale(ewaldcoeff_q, ewaldcoeff_lj);
}

TEST_F(TabscaleTest, SameEwaldCoeff)
{
    real ewaldcoeff_q  = 1e-4;
    real ewaldcoeff_lj = 1e-4;
    testTabscale(ewaldcoeff_q, ewaldcoeff_lj);
}

static void my_check(test::TestReferenceChecker  checker,
                     const std::vector<real>    &v,
                     const char                 *label)
{
    checker.checkSequence(v.begin(), v.end(), label);
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

        void testArray(int  eeltype,
                       real q_tol,
                       real rcoulomb,
                       int  vdwtype,
                       real lj_tol,
                       real rvdw)
        {
            real              ewaldcoeff_q  = calc_ewaldcoeff_q(rcoulomb, q_tol);
            real              ewaldcoeff_lj = calc_ewaldcoeff_q(rvdw, lj_tol);
            InteractionTables it(eeltype,
                                 ewaldcoeff_q,
                                 rcoulomb,
                                 vdwtype,
                                 ewaldcoeff_lj,
                                 rvdw);
            double beta_q      = ewaldcoeff_q;
            it.fillEwaldSpline3(&beta_q, etCoulomb, v_q_ewald_lr);
            my_check(checker_, it.coulombF(), "coulombF");
            my_check(checker_, it.coulombV(), "coulombV");
            my_check(checker_, it.coulombFDV0(), "coulombFDV0");

            double beta_lj = ewaldcoeff_lj;
            it.fillEwaldSpline3(&beta_lj, etLJ, v_lj_ewald_lr);
            my_check(checker_, it.vdwF(), "vdwF");
            my_check(checker_, it.vdwV(), "vdwV");
            my_check(checker_, it.vdwFDV0(), "vdwFDV0");
        }
};

TEST_F(ArrayTest, All)
{
    real q_tol    = 1e-3;
    real lj_tol   = 1e-2;
    real rcoulomb = 0.3;
    real rvdw     = 0.4;
    testArray(eelPME, q_tol, rcoulomb, evdwPME, lj_tol, rvdw);
}

} // namespace

} // namespace gmx
