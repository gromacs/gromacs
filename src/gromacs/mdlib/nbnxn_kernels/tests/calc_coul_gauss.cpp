/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2018, by the GROMACS development team, led by
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

#include "gmxpre.h"

#include <math.h>

#include <algorithm>
#include <iomanip>
#include <iostream>

#include <gtest/gtest.h>

#include "gromacs/math/functions.h"
#include "gromacs/math/vectypes.h"
#include "gromacs/mdlib/forcerec.h"
#include "gromacs/mdlib/nbnxn_pairlist.h"
#include "gromacs/mdlib/nbnxn_kernels/nbnxn_kernel_ref.h"
#include "gromacs/mdtypes/interaction_const.h"
#include "gromacs/pbcutil/ishift.h"
#include "gromacs/utility/smalloc.h"

#include "testutils/cmdlinetest.h"
#include "testutils/refdata.h"
#include "testutils/testasserts.h"
#include "testutils/testfilemanager.h"

namespace
{

class CoulombTest : public gmx::test::CommandLineTestBase
{

    protected:
        gmx::test::TestReferenceChecker checker_;
        //init set tolerance
        CoulombTest () : checker_(this->rootChecker())
        {
            gmx::test::FloatingPointTolerance tolerance = gmx::test::relativeToleranceAsFloatingPoint(1.0, 5e-5);
            checker_.setDefaultTolerance(tolerance);
        }
        // Static initiation, only run once every test.
        static void SetUpTestCase()
        {
        }
        void testCoulomb(real zeta1,
                         real zeta2,
                         real lambda)
        {
            std::cout << "1 \n";
            nbnxn_pairlist_t *nbl = new nbnxn_pairlist_t();

            nbl->nci = 1;
            std::cout << "2 \n";
            snew(nbl->ci, nbl->nci);

            /* nbnxn_ci_t */
            nbl->ci[0].ci           = 0;
            nbl->ci[0].shift        = 534;
            nbl->ci[0].cj_ind_start = 0;
            nbl->ci[0].cj_ind_end   = 1;
            std::cout << "3 \n";
            snew(nbl->cj, nbl->nci);
            /* nbnxn_cj_t */
            nbl->cj[0].cj   = 0;
            nbl->cj[0].excl = 2254;
            std::cout << "4 \n";
            nbnxn_atomdata_t *nbat = new nbnxn_atomdata_t();
            nbat->ntype = 3;
            std::cout << "5 \n";
            snew(nbat->q, 4);
            nbat->q[0] = -1.000;
            nbat->q[1] = 1.000;
            nbat->q[2] = 0.000;
            nbat->q[3] = 0.000;
            std::cout << "6 \n";
            snew(nbat->type, 4);
            nbat->type[0] = 0;
            nbat->type[1] = 1;
            nbat->type[2] = 2;
            nbat->type[3] = 2;
            std::cout << "7 \n";
            snew(nbat->x, 12);
            nbat->x[0]  = 0.000;
            nbat->x[1]  = 5.000;
            nbat->x[2]  = 10.000;
            nbat->x[3]  = 0.000;
            nbat->x[4]  = 5.000;
            nbat->x[5]  = 10.000;
            nbat->x[6]  = -1000000.000;
            nbat->x[7]  = -1000000.000;
            nbat->x[8]  = -1000000.000;
            nbat->x[9]  = -1000000.000;
            nbat->x[10] = -1000000.000;
            nbat->x[11] = -1000000.000;
            std::cout << "8 \n";
            snew(nbat->nbfp, nbat->ntype*nbat->ntype*2);
            std::cout << "9 \n";
            interaction_const_t *ic = new interaction_const_t();

            ic->rcoulomb              = 1.100000;
            ic->epsfac                = 138.935455;
            ic->repulsion_shift.cpot  = 0.000000;
            ic->dispersion_shift.cpot = 0.000000;
            ic->ewaldcoeff_q          = 3.424148;
            std::cout << "10 \n";
            real  zeta[] = {zeta1, zeta2};
            std::cout << "11 \n";
            real *zeta_matrix = make_zeta_matrix(2, zeta);
            std::cout << "12 \n";
            /* Now copy the zeta_matrix */
            snew(nbat->zeta_matrix, 13);
            for (int i = 0; i < 4; i++)
            {
                nbat->zeta_matrix[i] = zeta_matrix[i];
            }
            std::cout << "13 \n";
            std::vector<double> force;
            std::vector<double> potential;

            for (int i = 0; i < 9; i++)
            {
                /* Reset the x-coordinate of particle 2 */
                nbat->x[3] = lambda*(i+1);
                real f[12];
                memset(f, 0.0, sizeof(f));
                real fshift[12];
                memset(fshift, 0.0, sizeof(fshift));
                real Vvdw = 0.0;
                rvec shift_vec[SHIFTS];
                memset(shift_vec, 0.0, sizeof(shift_vec));
                real Vc = 0.0;
                /* Call the function to be tested */
                nbnxn_kernel_ElecGauss_VdwLJ_VF_ref(nbl, nbat, ic, shift_vec, f, fshift, &Vvdw, &Vc);
                force.push_back(f[0]);
                potential.push_back(Vc);
            }
            std::cout << "14 \n";
            char buf[256];
            snprintf(buf, sizeof(buf), "Potential");
            checker_.checkSequence(potential.begin(), potential.end(), buf);
            snprintf(buf, sizeof(buf), "Force");
            checker_.checkSequence(force.begin(), force.end(), buf);
            std::cout << "15 \n";
            /* Free memory */
            sfree(nbl->ci);
            sfree(nbl->cj);
            sfree(nbat->q);
            sfree(nbat->x);
            sfree(nbat->type);
            sfree(nbat->nbfp);
            sfree(nbat->zeta_matrix);
            std::cout << "16 \n";
            delete nbl;
            delete nbat;
            delete ic;
            std::cout << "17 \n";
            std::vector<double>().swap(force);
            std::vector<double>().swap(potential);
            std::cout << "18 \n";
            sfree(zeta_matrix);
        }
        static void TearDownTestCase()
        {
        }

};
TEST_F (CoulombTest, Gaussian)
{
    std::cout << "for Gaussian, calling testCoulomb \n";
    testCoulomb(70.0, 90.0, 0.1);
    std::cout << "called tesCoulomb for Gaussian \n";
}
TEST_F (CoulombTest, GaussianSmallDistance)
{
    std::cout << "for GaussianSmallDistance, calling testCoulomb \n";
    testCoulomb(70.0, 90.0, 0.001);
    std::cout << "called testCoulomb for GaussianSmallDistance\n";
}
TEST_F (CoulombTest, GaussianSameZi)
{
    std::cout << "for GaussianSameZi, calling testCoulomb\n";
    testCoulomb(70.0, 70.0, 0.001);
    std::cout << "called testCoulomb for GaussianSameZi\n";
}

}  // namespace
