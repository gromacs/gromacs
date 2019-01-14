/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2016,2017,2018,2019, by the GROMACS development team, led by
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
    // TODO This is pretty inefficient, perhaps we should just check
    // whether the central box fields have their values.
    checker->checkSequence(std::begin(output.fshift), std::end(output.fshift), "ShiftForces");
    checker->checkSequence(std::begin(output.f), std::end(output.f), "Forces");
}

class BondedTest : public ::testing::TestWithParam<std::tuple<std::vector<gmx::RVec>, int> >
{
    protected:
        matrix                     box_;
        t_pbc                      pbc_;
        std::vector<gmx::RVec>     x_;
        int                        epbc_;
        test::TestReferenceData    refData_;
        test::TestReferenceChecker checker_;
        BondedTest( ) :
            checker_(refData_.rootChecker())
        {
            // We need quite specific tolerances here since angle functions
            // etc. are not very precise and reproducible.
            test::FloatingPointTolerance tolerance(test::FloatingPointTolerance(1.0e-4, 2.0e-7,
                                                                                1.0e-6, 1.0e-12,
                                                                                1000000000, 10000, false));
            checker_.setDefaultTolerance(tolerance);
            x_ = std::get<0>(GetParam());
            clear_mat(box_);
            box_[0][0] = box_[1][1] = box_[2][2] = 1.5;
            epbc_      = std::get<1>(GetParam());
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
        void testIfunc(test::TestReferenceChecker *checker,
                       const int                   ftype,
                       const std::vector<t_iatom> &iatoms,
                       const t_iparams            &iparams,
                       const real                  lambda)
        {
            SCOPED_TRACE(std::string("Testing PBC ") + epbc_names[epbc_]);

            int                           ddgatindex = 0;
            OutputQuantities              output;
            output.energy = bondedFunction(ftype)(iatoms.size(),
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
};

TEST_P (BondedTest, BondAngle)
{
    testBondAngle();
}

TEST_P (BondedTest, DihedralAngle)
{
    testDihedralAngle();
}

TEST_P (BondedTest, IfuncBonds)
{
    std::vector<t_iatom> iatoms = { 0, 0, 1, 0, 1, 2, 0, 2, 3 };
    t_iparams            iparams;
    iparams.harmonic.rA  = iparams.harmonic.rB  = 0.8;
    iparams.harmonic.krA = iparams.harmonic.krB = 50.0;
    const real lambda = 0.0;
    testIfunc(&checker_, F_BONDS, iatoms, iparams, lambda);
}

TEST_P (BondedTest, IfuncAngles)
{
    std::vector<t_iatom> iatoms = { 0, 0, 1, 2, 0, 1, 2, 3 };
    t_iparams            iparams;
    real                 k = 50.0;
    iparams.harmonic.rA  = iparams.harmonic.rB  = 100.0;
    iparams.harmonic.krA = iparams.harmonic.krB = k;
    const real lambda = 0.0;
    testIfunc(&checker_, F_ANGLES, iatoms, iparams, lambda);
}

TEST_P (BondedTest, IfuncProperDihedrals)
{
    std::vector<t_iatom> iatoms = { 0, 0, 1, 2, 3 };
    t_iparams            iparams;
    iparams.pdihs.phiA = iparams.pdihs.phiB = -100.0;
    iparams.pdihs.cpA  = iparams.pdihs.cpB  = 10.0;
    iparams.pdihs.mult = 1;
    const real lambda = 0.0;
    testIfunc(&checker_, F_PDIHS, iatoms, iparams, lambda);
}

TEST_P (BondedTest, IfuncImproperDihedrals)
{
    std::vector<t_iatom> iatoms = { 0, 0, 1, 2, 3 };
    t_iparams            iparams;
    iparams.harmonic.rA  = iparams.harmonic.rB  = 0.0;
    iparams.harmonic.krA = iparams.harmonic.krB = 5.0;
    const real lambda = 0.0;
    testIfunc(&checker_, F_IDIHS, iatoms, iparams, lambda);
}

TEST_P (BondedTest, IfuncImproperDihedralsFEP)
{
    std::vector<t_iatom> iatoms = { 0, 0, 1, 2, 3 };
    t_iparams            iparams;
    iparams.harmonic.rA  = iparams.harmonic.rB  = 0.0;
    iparams.harmonic.krA = iparams.harmonic.krB = 5.0;
    iparams.harmonic.rB  = 35.5;
    iparams.harmonic.krB = 10.0;

    const int numLambdas = 3;
    for (int i = 0; i < numLambdas; ++i)
    {
        const real lambda       = i / (numLambdas - 1.0);
        auto       valueChecker = checker_.checkCompound("Lambda", toString(lambda));
        testIfunc(&valueChecker, F_IDIHS, iatoms, iparams, lambda);
    }
}

//! Coordinates for testing
std::vector<std::vector<gmx::RVec> > c_coordinatesForTests =
{
    {{  0.0, 0.0, 0.0 }, { 0.0, 0.0, 1.0 }, { 0.0, 1.0, 1.0 }, { 1.0, 1.0, 1.0 }},
    {{  0.5, 0.0, 0.0 }, { 0.5, 0.0, 0.15 }, { 0.5, 0.07, 0.22 }, { 0.5, 0.18, 0.22 }},
    {{ -0.1143, -0.0282, 0.0 }, { 0.0, 0.0434, 0.0 }, { 0.1185, -0.0138, 0.0 }, { -0.0195, 0.1498, 0.0 }}
};

//! PBC values for testing
std::vector<int> c_pbcForTests = { epbcNONE, epbcXY, epbcXYZ };

INSTANTIATE_TEST_CASE_P(ForPbcValues, BondedTest, ::testing::Combine(::testing::ValuesIn(c_coordinatesForTests), ::testing::ValuesIn(c_pbcForTests)));
}  // namespace

}  // namespace gmx
