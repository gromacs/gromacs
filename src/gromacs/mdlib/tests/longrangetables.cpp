/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2015,2016, by the GROMACS development team, led by
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

#include "gromacs/mdlib/longrangetables.h"

#include <vector>

#include <gtest/gtest.h>

#include "gromacs/math/calculate-ewald-splitting-coefficient.h"
#include "gromacs/mdtypes/md_enums.h"
#include "gromacs/tables/splinetable.h"
#include "gromacs/utility/snprintf.h"

#include "testutils/refdata.h"
#include "testutils/testasserts.h"
#include "testutils/testfilemanager.h"

namespace gmx
{

namespace
{

class LongRangeCorrectionFunctionTest : public ::testing::Test
{
    protected:
        test::TestReferenceData        refData_;
        test::TestReferenceChecker     checker_;
        LongRangeCorrectionFunctionTest( )
            : checker_(refData_.rootChecker())
        {
            checker_.setDefaultTolerance(test::relativeToleranceAsFloatingPoint(1, 1e-5));
        }

        void testTable(interactionPotentialFunction pot, const std::vector<double> &param)
        {
            std::vector<double> f_ewald, v_ewald;
            double              r, e, f, dx = 0.001;

            for (unsigned int i = 0; (i <= 250); i++)
            {
                r = i*dx;
                (*pot)(param, r, &e, &f);
                v_ewald.push_back(e);
                f_ewald.push_back(f);
            }
            checker_.checkSequenceArray(v_ewald.size(), &v_ewald[0], "LR-Corr-V");
            checker_.checkSequenceArray(f_ewald.size(), &f_ewald[0], "LR-Corr-F");
        }
};

TEST_F(LongRangeCorrectionFunctionTest, QEwald)
{
    std::vector<double> beta;
    beta.push_back(3.0);
    testTable(v_q_ewald_lr, beta);
}

TEST_F(LongRangeCorrectionFunctionTest, LJEwald)
{
    std::vector<double> beta;
    beta.push_back(3.0);
    testTable(v_lj_ewald_lr, beta);
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
                                double rvdw = 0.1*rv;
                                char   buf[256];
                                double scale = longRangeCorrectionTableScale(eeltype,
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
                                checker_.checkDouble(scale, buf);
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

} // namespace

} // namespace gmx
