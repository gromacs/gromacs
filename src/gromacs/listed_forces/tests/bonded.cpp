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
 * \ingroup module_listed_forces
 */
#include "gmxpre.h"

#include "gromacs/listed_forces/bonded.h"

#include <cmath>

#include <memory>
#include <unordered_map>

#include <gtest/gtest.h>

#include "gromacs/listed_forces/listed_forces.h"
#include "gromacs/math/units.h"
#include "gromacs/math/vec.h"
#include "gromacs/math/vectypes.h"
#include "gromacs/mdtypes/mdatom.h"
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
    public:
        //! Function type
        int       ftype = -1;
        //! Tolerance for float evaluation
        float     ftoler = 1e-6;
        //! Tolerance for double evaluation
        double    dtoler = 1e-8;
        //! Do free energy perturbation?
        bool      fep = false;
        //! Interaction parameters
        t_iparams iparams = {{ 0 }};

        //! Constructor
        iListInput() {}

        /*! \brief Constructor with tolerance
         *
         * \param[in] ftol Single precision tolerance
         * \param[in] dtol Double precision tolerance
         */
        iListInput(float ftol, double dtol)
        {
            ftoler = ftol;
            dtoler = dtol;
        }
        /*! \brief Set parameters for harmonic potential
         *
         * Free energy perturbation is turned on when A
         * and B parameters are different.
         * \param[in] ft  Function type
         * \param[in] rA  Equilibrium value A
         * \param[in] krA Force constant A
         * \param[in] rB  Equilibrium value B
         * \param[in] krB Force constant B
         * \return The structure itself.
         */
        iListInput setHarmonic(int ft, real rA, real krA, real rB, real krB)
        {
            iparams.harmonic.rA  = rA;
            iparams.harmonic.rB  = rB;
            iparams.harmonic.krA = krA;
            iparams.harmonic.krB = krB;
            ftype                = ft;
            fep                  = (rA != rB || krA != krB);
            return *this;
        }
        /*! \brief Set parameters for harmonic potential
         *
         * \param[in] ft  Function type
         * \param[in] rA  Equilibrium value
         * \param[in] krA Force constant
         * \return The structure itself.
         */
        iListInput setHarmonic(int ft, real rA, real krA)
        {
            return setHarmonic(ft, rA, krA, rA, krA);
        }
        /*! \brief Set parameters for cubic potential
         *
         * \param[in] b0   Equilibrium bond length
         * \param[in] kb   Harmonic force constant
         * \param[in] kcub Cubic force constant
         * \return The structure itself.
         */
        iListInput setCubic(real b0, real kb, real kcub)
        {
            ftype              = F_CUBICBONDS;
            iparams.cubic.b0   = b0;
            iparams.cubic.kb   = kb;
            iparams.cubic.kcub = kcub;
            return *this;
        }
        /*! \brief Set parameters for morse potential
         *
         * Free energy perturbation is turned on when A
         * and B parameters are different.
         * \param[in] b0A   Equilibrium value A
         * \param[in] cbA   Force constant A
         * \param[in] betaA Steepness parameter A
         * \param[in] b0B   Equilibrium value B
         * \param[in] cbB   Force constant B
         * \param[in] betaB Steepness parameter B
         * \return The structure itself.
         */
        iListInput setMorse(real b0A, real cbA, real betaA, real b0B, real cbB, real betaB)
        {
            ftype               = F_MORSE;
            iparams.morse.b0A   = b0A;
            iparams.morse.cbA   = cbA;
            iparams.morse.betaA = betaA;
            iparams.morse.b0B   = b0B;
            iparams.morse.cbB   = cbB;
            iparams.morse.betaB = betaB;
            fep                 = (b0A != b0B || cbA != cbB || betaA != betaB);
            return *this;
        }
        /*! \brief Set parameters for morse potential
         *
         * \param[in] b0A   Equilibrium value
         * \param[in] cbA   Force constant
         * \param[in] betaA Steepness parameter
         * \return The structure itself.
         */
        iListInput setMorse(real b0A, real cbA, real betaA)
        {
            return setMorse(b0A, cbA, betaA, b0A, cbA, betaA);
        }
        /*! \brief Set parameters for fene potential
         *
         * \param[in] bm Equilibrium bond length
         * \param[in] kb Force constant
         * \return The structure itself.
         */
        iListInput setFene(real bm, real kb)
        {
            ftype           = F_FENEBONDS;
            iparams.fene.bm = bm;
            iparams.fene.kb = kb;
            return *this;
        }
        /*! \brief Set parameters for linear angle potential
         *
         * Free energy perturbation is turned on when A
         * and B parameters are different.
         * \param[in] klinA Force constant A
         * \param[in] aA    The position of the central atom A
         * \param[in] klinB Force constant B
         * \param[in] aB    The position of the central atom B
         * \return The structure itself.
         */
        iListInput setLinearAngle(real klinA, real aA, real klinB, real aB)
        {
            ftype                  = F_LINEAR_ANGLES;
            iparams.linangle.klinA = klinA;
            iparams.linangle.aA    = aA;
            iparams.linangle.klinB = klinB;
            iparams.linangle.aB    = aB;
            fep                    = (klinA != klinB || aA != aB);
            return *this;
        }
        /*! \brief Set parameters for linear angle potential
         *
         * \param[in] klinA Force constant
         * \param[in] aA    The position of the central atom
         * \return The structure itself.
         */
        iListInput setLinearAngle(real klinA, real aA)
        {
            return setLinearAngle(klinA, aA, klinA, aA);
        }
        /*! \brief Set parameters for Urey Bradley potential
         *
         * Free energy perturbation is turned on when A
         * and B parameters are different.
         * \param[in] thetaA  Equilibrium angle A
         * \param[in] kthetaA Force constant A
         * \param[in] r13A    The distance between i and k atoms A
         * \param[in] kUBA    The force constant for 1-3 distance A
         * \param[in] thetaB  Equilibrium angle B
         * \param[in] kthetaB Force constant B
         * \param[in] r13B    The distance between i and k atoms B
         * \param[in] kUBB    The force constant for 1-3 distance B
         * \return The structure itself.
         */
        iListInput setUreyBradley(real thetaA, real kthetaA, real r13A, real kUBA,
                                  real thetaB, real kthetaB, real r13B, real kUBB)
        {
            ftype               = F_UREY_BRADLEY;
            iparams.u_b.thetaA  = thetaA;
            iparams.u_b.kthetaA = kthetaA;
            iparams.u_b.r13A    = r13A;
            iparams.u_b.kUBA    = kUBA;
            iparams.u_b.thetaB  =  thetaB;
            iparams.u_b.kthetaB = kthetaB;
            iparams.u_b.r13B    = r13B;
            iparams.u_b.kUBB    = kUBB;
            fep                 = (thetaA != thetaB || kthetaA != kthetaB || r13A != r13B || kUBA != kUBB);
            return *this;
        }
        /*! \brief Set parameters for Urey Bradley potential
         *
         * \param[in] thetaA  Equilibrium angle
         * \param[in] kthetaA Force constant
         * \param[in] r13A    The distance between i and k atoms
         * \param[in] kUBA    The force constant for 1-3 distance
         * \return The structure itself.
         */
        iListInput setUreyBradley(real thetaA, real kthetaA, real r13A, real kUBA)
        {
            return setUreyBradley(thetaA, kthetaA, r13A,  kUBA,
                                  thetaA, kthetaA, r13A,  kUBA);
        }
        /*! \brief Set parameters for Cross Bond Bonds potential
         *
         * \param[in] r1e  First bond length i-j
         * \param[in] r2e  Second bond length i-k
         * \param[in] krr  The force constant
         * \return The structure itself.
         */
        iListInput setCrossBondBonds(real r1e, real r2e, real krr)
        {
            ftype                = F_CROSS_BOND_BONDS;
            iparams.cross_bb.r1e = r1e;
            iparams.cross_bb.r2e = r2e;
            iparams.cross_bb.krr = krr;
            return *this;
        }
        /*! \brief Set parameters for Cross Bond Angles potential
         *
         * \param[in] r1e  First bond length i-j
         * \param[in] r2e  Second bond length j-k
         * \param[in] r3e  Third bond length i-k
         * \param[in] krt  The force constant
         * \return The structure itself.
         */
        iListInput setCrossBondAngles(real r1e, real r2e, real r3e, real krt)
        {
            ftype                = F_CROSS_BOND_ANGLES;
            iparams.cross_ba.r1e = r1e;
            iparams.cross_ba.r2e = r2e;
            iparams.cross_ba.r3e = r3e;
            iparams.cross_ba.krt = krt;
            return *this;
        }
        /*! \brief Set parameters for Quartic Angles potential
         *
         * \param[in] theta Angle
         * \param[in] c     Array of parameters
         * \return The structure itself.
         */
        iListInput setQuarticAngles(real theta, const real c[5])
        {
            ftype                = F_QUARTIC_ANGLES;
            iparams.qangle.theta = theta;
            iparams.qangle.c[0]  = c[0];
            iparams.qangle.c[1]  = c[1];
            iparams.qangle.c[2]  = c[2];
            iparams.qangle.c[3]  = c[3];
            iparams.qangle.c[4]  = c[4];
            return *this;
        }
        /*! \brief Set parameters for proper dihedrals potential
         *
         * Free energy perturbation is turned on when A
         * and B parameters are different.
         * \param[in] ft   Function type
         * \param[in] phiA Dihedral angle A
         * \param[in] cpA  Force constant A
         * \param[in] mult Multiplicity of the angle
         * \param[in] phiB Dihedral angle B
         * \param[in] cpB  Force constant B
         * \return The structure itself.
         */
        iListInput setPDihedrals(int ft, real phiA, real cpA, int mult, real phiB, real cpB)
        {
            ftype              = ft;
            iparams.pdihs.phiA = phiA;
            iparams.pdihs.cpA  = cpA;
            iparams.pdihs.phiB = phiB;
            iparams.pdihs.cpB  = cpB;
            iparams.pdihs.mult = mult;
            fep                = (phiA != phiB || cpA != cpB);
            return *this;
        }
        /*! \brief Set parameters for proper dihedrals potential
         *
         * \param[in] ft   Function type
         * \param[in] phiA Dihedral angle
         * \param[in] cpA  Force constant
         * \param[in] mult Multiplicity of the angle
         * \return The structure itself.
         */
        iListInput setPDihedrals(int ft, real phiA, real cpA, int mult)
        {
            return setPDihedrals(ft, phiA, cpA, mult, phiA, cpA);
        }
        /*! \brief Set parameters for Ryckaert-Bellemans dihedrals potential
         *
         * Free energy perturbation is turned on when A
         * and B parameters are different.
         * \param[in] rbcA Force constants A
         * \param[in] rbcB Force constants B
         * \return The structure itself.
         */
        iListInput setRbDihedrals(const real rbcA[NR_RBDIHS], const real rbcB[NR_RBDIHS])
        {
            ftype = F_RBDIHS;
            fep   = false;
            for (int i = 0; i < NR_RBDIHS; i++)
            {
                iparams.rbdihs.rbcA[i] = rbcA[i];
                iparams.rbdihs.rbcB[i] = rbcB[i];
                fep                    = fep || (rbcA[i] != rbcB[i]);
            }
            return *this;
        }
        /*! \brief Set parameters for Ryckaert-Bellemans dihedrals potential
         *
         * \param[in] rbc Force constants
         * \return The structure itself.
         */
        iListInput setRbDihedrals(const real rbc[NR_RBDIHS])
        {
            return setRbDihedrals(rbc, rbc);
        }
        /*! \brief Set parameters for Polarization
         *
         * \param[in] alpha Polarizability
         * \return The structure itself.
         */
        iListInput setPolarization(real alpha)
        {
            ftype                  = F_POLARIZATION;
            fep                    = false;
            iparams.polarize.alpha = alpha;
            return *this;
        }
        /*! \brief Set parameters for Anharmonic Polarization
         *
         * \param[in] alpha Polarizability (nm^3)
         * \param[in] drcut The cut-off distance (nm) after which the
         *                  fourth power kicks in
         * \param[in] khyp  The force constant for the fourth power
         * \return The structure itself.
         */
        iListInput setAnharmPolarization(real alpha,
                                         real drcut,
                                         real khyp)
        {
            ftype                         = F_ANHARM_POL;
            fep                           = false;
            iparams.anharm_polarize.alpha = alpha;
            iparams.anharm_polarize.drcut = drcut;
            iparams.anharm_polarize.khyp  = khyp;
            return *this;
        }
        /*! \brief Set parameters for Thole Polarization
         *
         * \param[in] a      Thole factor
         * \param[in] alpha1 Polarizability 1 (nm^3)
         * \param[in] alpha2 Polarizability 2 (nm^3)
         * \param[in] rfac   Distance factor
         * \return The structure itself.
         */
        iListInput setTholePolarization(real a,
                                        real alpha1,
                                        real alpha2,
                                        real rfac)
        {
            ftype                = F_THOLE_POL;
            fep                  = false;
            iparams.thole.a      = a;
            iparams.thole.alpha1 = alpha1;
            iparams.thole.alpha2 = alpha2;
            iparams.thole.rfac   = rfac;
            return *this;
        }
        /*! \brief Set parameters for Water Polarization
         *
         * \param[in] alpha_x Polarizability X (nm^3)
         * \param[in] alpha_y Polarizability Y (nm^3)
         * \param[in] alpha_z Polarizability Z (nm^3)
         * \param[in] rOH     Oxygen-Hydrogen distance
         * \param[in] rHH     Hydrogen-Hydrogen distance
         * \param[in] rOD     Oxygen-Dummy distance
         * \return The structure itself.
         */
        iListInput setWaterPolarization(real alpha_x,
                                        real alpha_y,
                                        real alpha_z,
                                        real rOH,
                                        real rHH,
                                        real rOD)
        {
            ftype             = F_WATER_POL;
            fep               = false;
            iparams.wpol.al_x = alpha_x;
            iparams.wpol.al_y = alpha_y;
            iparams.wpol.al_z = alpha_z;
            iparams.wpol.rOH  = rOH;
            iparams.wpol.rHH  = rHH;
            iparams.wpol.rOD  = rOD;
            return *this;
        }
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
      { 4, { 0, 0, 1, 2, 3 } },
      { 5, { 0, 0, 1, 2, 3, 0 } } };
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
            std::vector<int>  ddgatindex = { 0, 1, 2, 3 };
            std::vector<real> chargeA    = { 1.5, -2.0, 1.5, -1.0 };
            t_mdatoms         mdatoms    = {0};
            mdatoms.chargeA              = chargeA.data();
            OutputQuantities  output;
            output.energy = bondedFunction(input_.ftype) (iatoms.size(),
                                                          iatoms.data(),
                                                          &input_.iparams,
                                                          as_rvec_array(x_.data()),
                                                          output.f, output.fshift,
                                                          &pbc_,
                                                          /* const struct t_graph *g */ nullptr,
                                                          lambda, &output.dvdlambda,
                                                          &mdatoms,
                                                          /* struct t_fcdata * */ nullptr,
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

//! Function types for testing bonds. Add new terms at the end.
std::vector<iListInput> c_InputBonds =
{
    { iListInput().setHarmonic(F_BONDS, 0.15, 500.0) },
    { iListInput(2e-6f, 1e-8).setHarmonic(F_BONDS, 0.15, 500.0, 0.17, 400.0) },
    { iListInput(1e-4f, 1e-8).setHarmonic(F_G96BONDS, 0.15, 50.0) },
    { iListInput().setHarmonic(F_G96BONDS, 0.15, 50.0, 0.17, 40.0) },
    { iListInput().setCubic(0.16, 50.0, 2.0) },
    { iListInput(2e-6f, 1e-8).setMorse(0.15, 50.0, 2.0, 0.17, 40.0, 1.6) },
    { iListInput(2e-6f, 1e-8).setMorse(0.15, 30.0, 2.7) },
    { iListInput().setFene(0.4, 5.0) }
};

//! Constants for Quartic Angles
const real cQuarticAngles[5] = { 1.1, 2.3, 4.6, 7.8, 9.2 };

//! Function types for testing angles. Add new terms at the end.
std::vector<iListInput> c_InputAngles =
{
    { iListInput(2e-3, 1e-8).setHarmonic(F_ANGLES, 100.0, 50.0) },
    { iListInput(2e-3, 1e-8).setHarmonic(F_ANGLES, 100.15, 50.0, 95.0, 30.0) },
    { iListInput(8e-3, 1e-8).setHarmonic(F_G96ANGLES, 100.0, 50.0) },
    { iListInput(8e-3, 1e-8).setHarmonic(F_G96ANGLES, 100.0, 50.0, 95.0, 30.0) },
    { iListInput().setLinearAngle(50.0, 0.4) },
    { iListInput().setLinearAngle(50.0, 0.4, 40.0, 0.6) },
    { iListInput(2e-6, 1e-8).setCrossBondBonds(0.8, 0.7, 45.0) },
    { iListInput(2e-6, 1e-8).setCrossBondAngles(0.8, 0.7, 0.3, 45.0) },
    { iListInput(2e-2, 1e-8).setUreyBradley(950.0, 46.0, 0.3, 5.0) },
    { iListInput(2e-2, 1e-8).setUreyBradley(100.0, 45.0, 0.3, 5.0, 90.0, 47.0, 0.32, 7.0) },
    { iListInput(2e-3, 1e-8).setQuarticAngles(87.0, cQuarticAngles ) }
};

//! Constants for Ryckaert-Bellemans A
const real rbcA[NR_RBDIHS] = { -5.35, 13.6, 8.4, -16.7, 0.3, 12.4 };

//! Constants for Ryckaert-Bellemans B
const real rbcB[NR_RBDIHS] = { -6.35, 12.6, 8.1, -10.7, 0.9, 15.4 };

//! Constants for Ryckaert-Bellemans without FEP
const real rbc[NR_RBDIHS] = { -7.35, 13.6, 8.4, -16.7, 1.3, 12.4 };

//! Function types for testing dihedrals. Add new terms at the end.
std::vector<iListInput> c_InputDihs =
{
    { iListInput(1e-4, 1e-8).setPDihedrals(F_PDIHS, -100.0, 10.0, 2, -80.0, 20.0) },
    { iListInput(1e-4, 1e-8).setPDihedrals(F_PDIHS, -105.0, 15.0, 2) },
    { iListInput(2e-4, 1e-8).setHarmonic(F_IDIHS, 100.0, 50.0) },
    { iListInput(2e-4, 1e-8).setHarmonic(F_IDIHS, 100.15, 50.0, 95.0, 30.0) },
    { iListInput(4e-4, 1e-8).setRbDihedrals(rbcA, rbcB) },
    { iListInput(4e-4, 1e-8).setRbDihedrals(rbc) }
};

//! Function types for testing polarization. Add new terms at the end.
std::vector<iListInput> c_InputPols =
{
    { iListInput(2e-5, 1e-8).setPolarization(0.12) },
    { iListInput(1.7e-3, 1e-8).setAnharmPolarization(0.0013, 0.02, 1235.6) },
    { iListInput(1.4e-3, 1e-8).setTholePolarization(0.26, 0.07, 0.09, 1.6) },
    { iListInput().setWaterPolarization(0.001, 0.0012, 0.0016, 0.095, 0.15, 0.02) },
};

//! Function types for testing polarization. Add new terms at the end.
std::vector<iListInput> c_InputRestraints =
{
    { iListInput(1e-4, 1e-8).setPDihedrals(F_ANGRES, -100.0, 10.0, 2, -80.0, 20.0) },
    { iListInput(1e-4, 1e-8).setPDihedrals(F_ANGRES, -105.0, 15.0, 2) },
    { iListInput(1e-4, 1e-8).setPDihedrals(F_ANGRESZ, -100.0, 10.0, 2, -80.0, 20.0) },
    { iListInput(1e-4, 1e-8).setPDihedrals(F_ANGRESZ, -105.0, 15.0, 2) },
    { iListInput(2e-3, 1e-8).setHarmonic(F_RESTRANGLES, 100.0, 50.0) },
    { iListInput(2e-3, 1e-8).setHarmonic(F_RESTRANGLES, 100.0, 50.0, 110.0, 45.0) }
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

INSTANTIATE_TEST_CASE_P(Polarize, ListedForcesTest, ::testing::Combine(::testing::ValuesIn(c_InputPols), ::testing::ValuesIn(c_coordinatesForTests), ::testing::ValuesIn(c_pbcForTests)));

INSTANTIATE_TEST_CASE_P(Restraints, ListedForcesTest, ::testing::Combine(::testing::ValuesIn(c_InputRestraints), ::testing::ValuesIn(c_coordinatesForTests), ::testing::ValuesIn(c_pbcForTests)));

}  // namespace

}  // namespace gmx
