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

/*! \brief Utility to to fill iparams struct
 *
 * \param[in]  ftype Function type
 * \param[in]  fep   Whether or not to use free energy
 * \param[out] ip    Pointer to iparams struct
 */
void fillIparams(int ftype, bool fep, t_iparams *ip)
{
    switch (ftype)
    {
        case F_BONDS:
        case F_G96BONDS:
            if (fep)
            {
                ip->harmonic = { 0.8, 50.0, 0.6, 40.0 };
            }
            else
            {
                ip->harmonic = { 0.8, 50.0, 0.8, 50.0 };
            }
            break;
        case F_MORSE:
            ip->morse = { 0.8, 50.0, 2.0, 0.8, 50.0, 2.0 };
            break;
        case F_CUBICBONDS:
            ip->cubic = { 0.8, 50.0, 2.0 };
            break;
        case F_FENEBONDS:
            ip->fene = { 1.01, 5.0 };
            break;
        case F_ANGLES:
            if (fep)
            {
                ip->harmonic = { 100, 50.0, 95.0, 30.0 };
            }
            else
            {
                ip->harmonic = { 100, 50.0, 100.0, 50.0 };
            }
            break;
        case F_PDIHS:
            if (fep)
            {
                ip->pdihs = { -100.0, 10, 1, -80, 20 };
            }
            else
            {
                ip->pdihs = { -100.0, 10, 1, -100, 10 };
            }
            break;
        case  F_IDIHS:
            if (fep)
            {
                ip->harmonic = { 0.0, 5.0, 0.0, 5.0 };
            }
            else
            {
                ip->harmonic = { 0.0, 5.0, 35.5, 10.0 };
            }
            break;
        default:
            break;
    }
}

/*! \brief Utility to to fill iatoms struct
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
    int nral = interaction_function[ftype].nratoms;
    for (auto &i : ia[nral])
    {
        iatoms->push_back(i);
    }
}

class ListedForcesTest : public ::testing::TestWithParam<std::tuple<std::vector<gmx::RVec>, int, int, bool> >
{
    protected:
        matrix                     box_;
        t_pbc                      pbc_;
        std::vector<gmx::RVec>     x_;
        int                        epbc_;
        int                        ftype_;
        bool                       fep_;
        test::TestReferenceData    refData_;
        test::TestReferenceChecker checker_;
        ListedForcesTest( ) :
            checker_(refData_.rootChecker())
        {
            // We need quite specific tolerances here since angle functions
            // etc. are not very precise and reproducible.
            test::FloatingPointTolerance tolerance(test::FloatingPointTolerance(1.0e-4, 1.0e-6,
                                                                                1.0e-6, 1.0e-12,
                                                                                1000000000, 10000, false));
            checker_.setDefaultTolerance(tolerance);
            x_     = std::get<0>(GetParam());
            epbc_  = std::get<1>(GetParam());
            ftype_ = std::get<2>(GetParam());
            fep_   = std::get<3>(GetParam());
            clear_mat(box_);
            box_[0][0] = box_[1][1] = box_[2][2] = 1.5;
            set_pbc(&pbc_, epbc_, box_);

        }

        void testBondAngle()
        {
            rvec  r_ij, r_kj;
            real  cosine_angle, angle;
            int   t1, t2;

            angle = bond_angle(x_[0], x_[1], x_[2], &pbc_,
                               r_ij, r_kj, &cosine_angle,
                               &t1, &t2);
            checker_.checkReal(angle, "angle");
            checker_.checkReal(cosine_angle, "cosine_angle");
            checker_.checkInteger(t1, "t1");
            checker_.checkInteger(t2, "t2");
        }

        void testDihedralAngle()
        {
            rvec  r_ij, r_kj, r_kl, m, n;
            real  angle;
            int   t1, t2, t3;

            angle = dih_angle(x_[0], x_[1], x_[2], x_[3], &pbc_,
                              r_ij, r_kj, r_kl, m, n,
                              &t1, &t2, &t3);

            checker_.checkReal(angle, "angle");
            checker_.checkInteger(t1, "t1");
            checker_.checkInteger(t2, "t2");
            checker_.checkInteger(t3, "t3");
        }
        void testOneIfunc(test::TestReferenceChecker *checker,
                          const std::vector<t_iatom> &iatoms,
                          const t_iparams            &iparams,
                          const real                  lambda)
        {
            SCOPED_TRACE(std::string("Testing PBC ") + epbc_names[epbc_]);
            int                           ddgatindex = 0;
            OutputQuantities              output;
            output.energy = bondedFunction(ftype_)(iatoms.size(),
                                                   iatoms.data(),
                                                   &iparams,
                                                   as_rvec_array(x_.data()),
                                                   output.f, output.fshift,
                                                   &pbc_,
                                                   /* const struct t_graph *g */ nullptr,
                                                   lambda, &output.dvdlambda,
                                                   /* const struct t_mdatoms *md */ nullptr,
                                                   /* struct t_fcdata *fcd */ nullptr,
                                                   &ddgatindex);
            checkOutput(checker, output);
        }
        void testIfunc()
        {
            test::TestReferenceChecker thisChecker =
                checker_.checkCompound("FunctionType",
                                       interaction_function[ftype_].name).checkCompound("FEP", (fep_ ? "Yes" : "No"));
            t_iparams            iparams = {{0}};
            fillIparams(ftype_, fep_, &iparams);
            std::vector<t_iatom> iatoms;
            fillIatoms(ftype_, &iatoms);
            if (fep_)
            {
                const int numLambdas = 3;
                for (int i = 0; i < numLambdas; ++i)
                {
                    const real lambda       = i / (numLambdas - 1.0);
                    auto       valueChecker = thisChecker.checkCompound("Lambda", toString(lambda));
                    testOneIfunc(&valueChecker, iatoms, iparams, lambda);
                }
            }
            else
            {
                testOneIfunc(&thisChecker, iatoms, iparams, 0.0);
            }
        }
};

TEST_P (ListedForcesTest, Ifunc)
{
    testIfunc();
}

//! Function types for testing
std::vector<int> c_fTypes =
{
    F_BONDS, F_G96BONDS, F_CUBICBONDS, F_MORSE, F_FENEBONDS,
    F_ANGLES,
    F_PDIHS, F_IDIHS
};

//! Coordinates for testing
std::vector<std::vector<gmx::RVec> > c_coordinatesForTests =
{
    {{  0.0, 0.0, 0.0 }, { 0.0, 0.0, 1.0 }, { 0.0, 1.0, 1.0 }, { 1.0, 1.0, 1.0 }},
    {{  0.5, 0.0, 0.0 }, { 0.5, 0.0, 0.15 }, { 0.5, 0.07, 0.22 }, { 0.5, 0.18, 0.22 }},
    {{ -0.1143, -0.0282, 0.0 }, { 0.0, 0.0434, 0.0 }, { 0.1185, -0.0138, 0.0 }, { -0.0195, 0.1498, 0.0 }}
};

//! PBC values for testing
std::vector<int> c_pbcForTests = { epbcNONE, epbcXY, epbcXYZ };

//! FEP on or off
std::vector<bool> c_Fep = { false, true };

INSTANTIATE_TEST_CASE_P(Ifunc, ListedForcesTest, ::testing::Combine(::testing::ValuesIn(c_coordinatesForTests), ::testing::ValuesIn(c_pbcForTests), ::testing::ValuesIn(c_fTypes), ::testing::ValuesIn(c_Fep)));
}  // namespace

}  // namespace gmx
