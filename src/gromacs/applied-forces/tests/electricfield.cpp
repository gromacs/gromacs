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
 * Tests for functionality of the "angle" trajectory analysis module.
 *
 * \author David van der Spoel <david.vanderspoel@icm.uu.se>
 * \ingroup module_applied_forces
 */
#include "gmxpre.h"

#include "gromacs/applied-forces/electricfield.h"

#include <gtest/gtest.h>

#include "gromacs/fileio/readinp.h"
#include "gromacs/fileio/warninp.h"
#include "gromacs/gmxlib/network.h"
#include "gromacs/math/vec.h"
#include "gromacs/mdlib/forcerec.h"
#include "gromacs/mdrunutility/mdmodules.h"
#include "gromacs/mdtypes/forcerec.h"
#include "gromacs/mdtypes/inputrec.h"
#include "gromacs/mdtypes/mdatom.h"
#include "gromacs/utility/cstringutil.h"
#include "gromacs/utility/real.h"
#include "gromacs/utility/smalloc.h"

#include "testutils/testasserts.h"

namespace
{

/********************************************************************
 * ElectricFieldTest
 */

class ElectricFieldTest : public ::testing::Test
{
    public:
        ElectricFieldTest() {}

        void test(bool newMdp,
                  int  dim,
                  real E0,
                  real omega,
                  real t0,
                  real sigma,
                  real expectedValue)
        {
            gmx::test::FloatingPointTolerance tolerance(
                    gmx::test::relativeToleranceAsFloatingPoint(1.0, 0.005));
            gmx::MDModules                    module;
            t_inputrec *inputrec = module.inputrec();
            // Prepare MDP inputs
            int         ninp = 0;
            t_inpfile  *inp  = NULL;
            char        name[128];
            char        value[128];
            const char *dimXYZ[3] = { "X", "Y", "Z" };
            GMX_RELEASE_ASSERT((dim >= 0 && dim < 3), "Dimension should be 0, 1 or 2");
            t_inpfile   tmp;
            tmp.bObsolete = FALSE;
            tmp.bSet      = TRUE;
            tmp.inp_count = 0;
            tmp.count     = 0;
            tmp.name      = NULL;
            tmp.value     = NULL;
            if (newMdp)
            {
                ninp = 1;
                snew(inp, ninp);
                inp[0] = tmp;
                snprintf(name, sizeof(name), "ElectricField-%s", dimXYZ[dim]);
                snprintf(value, sizeof(value), "%g %g %g %g", E0, omega, t0, sigma);
                inp[0].name      = gmx_strdup(name);
                inp[0].value     = gmx_strdup(value);

            }
            else
            {
                ninp = 2;
                snew(inp, ninp);
                inp[0] = tmp;
                snprintf(name, sizeof(name), "E%s", dimXYZ[dim]);
                snprintf(value, sizeof(value), "1 %g 0", E0);
                inp[0].name      = gmx_strdup(name);
                inp[0].value     = gmx_strdup(value);

                inp[1]       = tmp;
                inp[1].count = 1;
                snprintf(name, sizeof(name), "E%s-t", dimXYZ[dim]);
                snprintf(value, sizeof(value), "3 %g 0 %g 0 %g 0", omega, t0, sigma);
                inp[1].name      = gmx_strdup(name);
                inp[1].value     = gmx_strdup(value);
            }

            warninp_t  wi = init_warning(TRUE, 0);

            inputrec->efield->readMdp(&ninp, &inp, wi);
            t_mdatoms md;
            rvec      f[1];
            clear_rvec(f[0]);
            md.homenr = 1;
            snew(md.chargeA, md.homenr);
            md.chargeA[0] = 1;

            t_commrec  *cr       = init_commrec();
            t_forcerec *forcerec = mk_forcerec();
            inputrec->efield->initForcerec(forcerec);
            forcerec->efield->calculateForces(cr, &md, f, 0);

            EXPECT_REAL_EQ_TOL(f[0][dim], expectedValue, tolerance);
            for (int i = 0; i < ninp; i++)
            {
                sfree(inp[i].name);
                sfree(inp[i].value);
            }
            done_warning(wi, 0, "no file", 0);
        }
};

TEST_F(ElectricFieldTest, NewStatic)
{
    test(true, 0, 1, 0, 0, 0, 96.4853363);
}

TEST_F(ElectricFieldTest, NewOscillating)
{
    test(true, 0, 1, 5, 0.2, 0, 96.4853363);
}

TEST_F(ElectricFieldTest, NewPulsed)
{
    test(true, 0, 1, 5, 0.5, 1, -68.215782);
}

TEST_F(ElectricFieldTest, OldStatic)
{
    test(false, 0, 1, 0, 0, 0, 96.4853363);
}

TEST_F(ElectricFieldTest, OldOscillating)
{
    test(false, 0, 1, 5, 0.2, 0, 96.4853363);
}

TEST_F(ElectricFieldTest, OldPulsed)
{
    test(false, 0, 1, 5, 0.5, 1, -68.215782);
}

} // namespace
