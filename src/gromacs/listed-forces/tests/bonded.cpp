/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2018,2019, by the GROMACS development team, led by
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
#include <unordered_map>

#include <gtest/gtest.h>

#include "gromacs/listed-forces/listed-forces.h"
#include "gromacs/math/units.h"
#include "gromacs/math/vec.h"
#include "gromacs/math/vectypes.h"
#include "gromacs/pbcutil/ishift.h"
#include "gromacs/pbcutil/pbc.h"
#include "gromacs/topology/idef.h"
#include "gromacs/utility/strconvert.h"

#include "testutils/refdata.h"
#include "testutils/testasserts.h"

namespace gmx
{
namespace
{

//! Number of atoms used in these tests.
constexpr int c_numAtoms = 4;

/*! \brief Output from bonded kernels
 *
 * \todo Later this might turn into the actual output struct. */
struct OutputQuantities
{
    //! Energy of this interaction
    real  energy         = 0;
    //! Derivative with respect to lambda
    real  dvdlambda      = 0;
    //! Shift vectors
    rvec  fshift[N_IVEC] = {{0}};
    //! Forces
    rvec4 f[c_numAtoms]  = {{0}};
};

/*! \brief Utility to check the output from bonded tests
 *
 * \param[in] checker Reference checker
 * \param[in] output  The output from the test to check
 */
void checkOutput(test::TestReferenceChecker *checker,
                 const OutputQuantities     &output)
{
    checker->checkReal(output.energy, "Epot ");
    // Should still be zero if not doing FEP, so may as well test it.
    checker->checkReal(output.dvdlambda, "dVdlambda ");
    checker->checkVector(output.fshift[CENTRAL], "Central shift forces");
    checker->checkSequence(std::begin(output.f), std::end(output.f), "Forces");
}

/*! \brief Input structure for listed forces tests
 */
struct iListInput
{
    //! Function type
    int       ftype;
    //! Tolerance for float evaluation
    float     ftoler;
    //! Tolerance for double evaluation
    double    dtoler;
    //! Do free energy perturbation?
    bool      fep;
    //! Interaction parameters
    t_iparams iparams;
};

/*! \brief Utility to fill iatoms struct
 *
 * \param[in]  ftype  Function type
 * \param[out] iatoms Pointer to iatoms struct
 */
void fillIatoms(int ftype, std::vector<t_iatom> *iatoms)
{
    std::unordered_map<int, std::vector<int> > ia =
    { { 2, { 0, 0, 1, 0, 1, 2, 0, 2, 3 } },
      { 3, { 0, 0, 1, 2, 0, 1, 2, 3 } },
      { 4, { 0, 0, 1, 2, 3 } }};
    EXPECT_TRUE(ftype >= 0 && ftype < F_NRE);
    int nral = interaction_function[ftype].nratoms;
    for (auto &i : ia[nral])
    {
        iatoms->push_back(i);
    }
}

class ListedForcesTest : public ::testing::TestWithParam<std::tuple<iListInput, std::vector<gmx::RVec>, int> >
{
    protected:
        matrix                     box_;
        t_pbc                      pbc_;
        std::vector<gmx::RVec>     x_;
        int                        epbc_;
        iListInput                 input_;
        test::TestReferenceData    refData_;
        test::TestReferenceChecker checker_;
        ListedForcesTest( ) :
            checker_(refData_.rootChecker())
        {
            input_ = std::get<0>(GetParam());
            x_     = std::get<1>(GetParam());
            epbc_  = std::get<2>(GetParam());
            clear_mat(box_);
            box_[0][0] = box_[1][1] = box_[2][2] = 1.5;
            set_pbc(&pbc_, epbc_, box_);

            // We need quite specific tolerances here since angle functions
            // etc. are not very precise and reproducible.
            test::FloatingPointTolerance tolerance(test::FloatingPointTolerance(input_.ftoler, 1.0e-6,
                                                                                input_.dtoler, 1.0e-12,
                                                                                10000, 100, false));
            checker_.setDefaultTolerance(tolerance);
        }
        void testOneIfunc(test::TestReferenceChecker *checker,
                          const std::vector<t_iatom> &iatoms,
                          const real                  lambda)
        {
            SCOPED_TRACE(std::string("Testing PBC ") + epbc_names[epbc_]);
            std::vector<int> ddgatindex = { 0, 1, 2, 3 };
            OutputQuantities output;
            output.energy = bondedFunction(input_.ftype)(iatoms.size(),
                                                         iatoms.data(),
                                                         &input_.iparams,
                                                         as_rvec_array(x_.data()),
                                                         output.f, output.fshift,
                                                         &pbc_,
                                                         /* const struct t_graph *g */ nullptr,
                                                         lambda, &output.dvdlambda,
                                                         /* const struct t_mdatoms *md */ nullptr,
                                                         /* struct t_fcdata *fcd */ nullptr,
                                                         ddgatindex.data());
            // Internal consistency test of both test input
            // and bonded functions.
            EXPECT_TRUE((input_.fep || (output.dvdlambda == 0.0)));
            checkOutput(checker, output);
        }
        void testIfunc()
        {
            test::TestReferenceChecker thisChecker =
                checker_.checkCompound("FunctionType",
                                       interaction_function[input_.ftype].name).checkCompound("FEP", (input_.fep ? "Yes" : "No"));
            std::vector<t_iatom> iatoms;
            fillIatoms(input_.ftype, &iatoms);
            if (input_.fep)
            {
                const int numLambdas = 3;
                for (int i = 0; i < numLambdas; ++i)
                {
                    const real lambda       = i / (numLambdas - 1.0);
                    auto       valueChecker = thisChecker.checkCompound("Lambda", toString(lambda));
                    testOneIfunc(&valueChecker, iatoms, lambda);
                }
            }
            else
            {
                testOneIfunc(&thisChecker, iatoms, 0.0);
            }
        }
};

TEST_P (ListedForcesTest, Ifunc)
{
    testIfunc();
}

//! Default absolute floating point tolerance
const float  ftol = 1e-6;
//! Default absolute floating point tolerance
const double dtol = 1e-8;
//! Function types for testing bonds. Add new terms at the end.
std::vector<iListInput> c_InputBonds =
{
    { F_BONDS, ftol, dtol, false,
      { .harmonic.rA = 0.15, .harmonic.krA = 500.0, .harmonic.rB = 0.15, .harmonic.krB = 500.0 } },
    { F_BONDS, ftol, dtol, true,
      { .harmonic.rA = 0.15, .harmonic.krA = 500.0, .harmonic.rB = 0.17, .harmonic.krB = 400.0 } },
    { F_G96BONDS, ftol, dtol, false,
      { .harmonic.rA = 0.15, .harmonic.krA = 50.0, .harmonic.rB = 0.15, .harmonic.krB = 50.0 } },
    { F_G96BONDS, ftol, dtol, true,
      { .harmonic.rA = 0.15, .harmonic.krA = 50.0, .harmonic.rB = 0.17, .harmonic.krB = 40.0 } }, 
    { F_CUBICBONDS, ftol, dtol, false,
      { .cubic.b0 = 0.16, .cubic.kb = 50.0, .cubic.kcub = 2.0 } },
    { F_MORSE, ftol, dtol, true,
      { .morse.b0A = 0.15, .morse.cbA = 50.0, .morse.betaA = 2.0,
        .morse.b0B = 0.17, .morse.cbB = 40.0, .morse.betaB = 1.6 } },
    { F_MORSE, ftol, dtol, false,
      { .morse.b0A = 0.15, .morse.cbA = 50.0, .morse.betaA = 2.0,
        .morse.b0B = 0.15, .morse.cbB = 50.0, .morse.betaB = 2.0 } },
    { F_FENEBONDS, ftol, dtol, false, 
      { .fene.bm = 0.4, .fene.kb = 5.0 } }
};

//! Function types for testing angles. Add new terms at the end.
std::vector<iListInput> c_InputAngles =
{
    { F_ANGLES, ftol, dtol, false,
      { .harmonic.rA = 100.0, .harmonic.krA = 50.0, .harmonic.rB = 100.0, .harmonic.krB = 50.0 } },
    { F_ANGLES, ftol, dtol, true,
      { .harmonic.rA = 100.0, .harmonic.krA = 50.0, .harmonic.rB = 95.0, .harmonic.krB = 30.0 } },
    { F_G96ANGLES, ftol, dtol, false,
      { .harmonic.rA = 100.0, .harmonic.krA = 50.0, .harmonic.rB = 100.0, .harmonic.krB = 50.0 } },
    { F_G96ANGLES, ftol, dtol, true,
      { .harmonic.rA = 100.0, .harmonic.krA = 50.0, .harmonic.rB = 95.0, .harmonic.krB = 30.0 } },
    { F_LINEAR_ANGLES, ftol, dtol, true,
      { .linangle.klinA = 50.0, .linangle.aA = 0.4, .linangle.klinB = 40.0, .linangle.aB = 0.6 } },
    { F_LINEAR_ANGLES, ftol, dtol, false,
      { .linangle.klinA = 50.0, .linangle.aA = 0.5, .linangle.klinB = 50.0, .linangle.aB = 0.5 } },
    { F_CROSS_BOND_BONDS, ftol, dtol, false,
      { .cross_bb.r1e = 0.8, .cross_bb.r2e = 0.7, .cross_bb.krr = 45.0 } },
    { F_CROSS_BOND_ANGLES, ftol, dtol, false,
      { .cross_bb.r1e = 0.8, .cross_bb.r2e = 0.7, .cross_bb.krr = 45.0 } },
    { F_UREY_BRADLEY, ftol, dtol, true,
      { .u_b.thetaA = 100.0, .u_b.kthetaA = 45.0, .u_b.r13A = 0.3,  .u_b.kUBA = 5.0,
        .u_b.thetaB =  90.0, .u_b.kthetaB = 47.0, .u_b.r13B = 0.32, .u_b.kUBB = 7.0 } },
    { F_UREY_BRADLEY, ftol, dtol, false,
      { .u_b.thetaA =  90.0, .u_b.kthetaA = 47.0, .u_b.r13A = 0.31, .u_b.kUBA = 6.0,
        .u_b.thetaB =  90.0, .u_b.kthetaB = 47.0, .u_b.r13B = 0.31, .u_b.kUBB = 6.0 } },
    { F_QUARTIC_ANGLES, ftol, dtol, false,
      { .qangle.theta = 87.0, .qangle.c = { 1.1, 2.2, 3.3, 4.5, 6.7 } } }
};
  
//! Function types for testing dihedrals. Add new terms at the end.
std::vector<iListInput> c_InputDihs =
{
    { F_PDIHS, ftol, dtol, true,
      { .pdihs.phiA = -100.0, .pdihs.cpA = 10, .pdihs.mult = 2,
        .pdihs.phiB =  -80.0, .pdihs.cpB = 20 } },
    { F_PDIHS, ftol, dtol, false,
      { .pdihs.phiA = -90.0, .pdihs.cpA = 30, .pdihs.mult = 2,
        .pdihs.phiB = -90.0, .pdihs.cpB = 30 } },
    { F_IDIHS, ftol, dtol, true,
      { .harmonic.rA =  0.0, .harmonic.krA = 5.0, .harmonic.rB = 35.5, .harmonic.krB = 10.0 } },
    { F_IDIHS, ftol, dtol, false,
      { .harmonic.rA = 10.0, .harmonic.krA = 7.0, .harmonic.rB = 10.0, .harmonic.krB =  7.0 } },
    { F_RBDIHS, ftol, dtol, true,
      { .rbdihs.rbcA = { -5.35, 13.6, 8.4, -16.7, 0.3, 12.4 },
        .rbdihs.rbcB = { -6.35, 12.6, 8.1, -10.7, 0.9, 15.4 } } },
    { F_RBDIHS, ftol, dtol, false,
      { .rbdihs.rbcA = { -7.35, 13.6, 8.4, -16.7, 1.3, 12.4 },
        .rbdihs.rbcB = { -7.35, 13.6, 8.4, -16.7, 1.3, 12.4 } } },
    { F_CBTDIHS, ftol, dtol, true,
      { .cbtdihs.cbtcA = { -5.35, 13.6, 8.4, -16.7, 0.3, 12.4 },
        .cbtdihs.cbtcB = { -6.35, 12.6, 8.1, -10.7, 0.9, 15.4 } } },
    { F_CBTDIHS, ftol, dtol, false,
      { .cbtdihs.cbtcA = { -7.35, 13.6, 8.4, -16.7, 1.3, 12.4 },
        .cbtdihs.cbtcB = { -7.35, 13.6, 8.4, -16.7, 1.3, 12.4 } } },
};
 
//! Coordinates for testing
std::vector<std::vector<gmx::RVec> > c_coordinatesForTests =
{
    {{  0.0, 0.0, 0.0 }, { 0.0, 0.0, 0.2 }, { 0.005, 0.0, 0.1 }, { -0.001, 0.1, 0.0 }},
    {{  0.5, 0.0, 0.0 }, { 0.5, 0.0, 0.15 }, { 0.5, 0.07, 0.22 }, { 0.5, 0.18, 0.22 }},
    {{ -0.1143, -0.0282, 0.0 }, { 0.0, 0.0434, 0.0 }, { 0.1185, -0.0138, 0.0 }, { -0.0195, 0.1498, 0.0 }}
};

//! PBC values for testing
std::vector<int> c_pbcForTests = { epbcNONE, epbcXY, epbcXYZ };

INSTANTIATE_TEST_CASE_P(Bond, ListedForcesTest, ::testing::Combine(::testing::ValuesIn(c_InputBonds), ::testing::ValuesIn(c_coordinatesForTests), ::testing::ValuesIn(c_pbcForTests)));

INSTANTIATE_TEST_CASE_P(Angle, ListedForcesTest, ::testing::Combine(::testing::ValuesIn(c_InputAngles), ::testing::ValuesIn(c_coordinatesForTests), ::testing::ValuesIn(c_pbcForTests)));

INSTANTIATE_TEST_CASE_P(Dihedral, ListedForcesTest, ::testing::Combine(::testing::ValuesIn(c_InputDihs), ::testing::ValuesIn(c_coordinatesForTests), ::testing::ValuesIn(c_pbcForTests)));

}  // namespace

}  // namespace gmx
