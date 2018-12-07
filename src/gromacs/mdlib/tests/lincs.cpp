/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2004, The GROMACS development team.
 * Copyright (c) 2013,2014,2015,2016,2017,2018, by the GROMACS development team, led by
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
 * \brief LINCS tests.
 *
 * \author Artem Zhmurov <zhmurov@gmail.com>
 * \ingroup module_mdlib
 */

#include "gmxpre.h"

#include <assert.h>

#include <cmath>

#include <algorithm>
#include <vector>

#include <gtest/gtest.h>

#include "gromacs/gmxlib/nonbonded/nonbonded.h"
#include "gromacs/math/vec.h"
#include "gromacs/math/vectypes.h"
#include "gromacs/mdlib/constr.h"
#include "gromacs/mdlib/gmx_omp_nthreads.h"
#include "gromacs/mdlib/lincs.h"
#include "gromacs/mdtypes/commrec.h"
#include "gromacs/mdtypes/inputrec.h"
#include "gromacs/mdtypes/mdatom.h"
#include "gromacs/pbcutil/pbc.h"
#include "gromacs/topology/block.h"
#include "gromacs/topology/idef.h"
#include "gromacs/topology/ifunc.h"
#include "gromacs/topology/topology.h"
#include "gromacs/utility/smalloc.h"
#include "gromacs/utility/stringutil.h"
#include "gromacs/utility/unique_cptr.h"

#include "gromacs/mdlib/tests/lincs_data.h"
#include "testutils/refdata.h"
#include "testutils/testasserts.h"

namespace gmx
{
namespace test
{

/*! \brief Test fixture for LINCS
 *
 * The fixture uses three simple and one real-life system for testing.
 * Simple systems include two atoms, connected with one constrain, three
 * atoms, connected consequently and three atoms, forming a triangle of
 * constraints. The real life system is a CVW tri-peptide with constraints
 * both on bonds, containing hydrogens, and on all bonds.
 *
 */
class LincsTest : public ::testing::Test
{
    public:

        /*! \brief
         * Method used to initialize atoms, constraints and LINCS parameters.
         *
         * \param[in] natom          Number of atoms in the system.
         * \param[in] masses         Atom masses. Size of this vector should be equal to natom.
         * \param[in] constraints    List of constraints, organized in tripples of integers.
         *                           First integer is the index of type for a constrain, second
         *                           and third are the indices of constrained atoms. The types
         *                           of constraints shoulb be sequantial (which is the way they
         *                           normally are in GROMACS).
         * \param[in] constraintsR0  Target values for bond lengths for bonds of each type. The
         *                           size of this vector should be equal to the total number of
         *                           unique types in constraints vector.
         * \param[in] firstType      The index of first interaction type used for constraints
         *                           (i.e. the lowest number that appear as type in constraints
         *                           vector). This is not fail-safe, since technically the
         *                           numbering of the types may nod be sequantial, but (i) this
         *                           is the way things are in GROMACS, (ii) everything is fine
         *                           in the examples provided for the test.
         * \param[in] epbc           Type of PBC (epbcXYZ / epbcNONE are used in this test,
         *                           see src/gromacs/pbcutil/pbc.h for details).
         * \param[in] box            The PBC box.
         * \param[in] nLincsIter     Number of iterations used to compute the inverse matrix.
         * \param[in] nProjOrder     The order for algorithm that adjusts the direction of the
         *                           bond after constraints are applied.
         * \param[in] LincsWarnAngle The value for the change in bond angle after which the
         *                           program will issue a warning.
         * \param[in] init_t         Initial time.
         * \param[in] delta_t        Timestep.
         *
         */
        void initSystem(int natom, std::vector<real> masses,
                        std::vector<int> constraints, std::vector<real> constraintsR0, int firstType,
                        int epbc, const matrix box,
                        int nLincsIter, int nProjOrder, real LincsWarnAngle, real init_t, real delta_t)
        {

            this->n = natom;   // Number of atoms

            invmass.resize(n); // Vector of inverce masses

            for (int i = 0; i < n; i++)
            {
                invmass.at(i) = 1.0/masses.at(i);
            }

            this->constraints   = constraints;   // Constraints indexes (in type-i-j format)
            this->constraintsR0 = constraintsR0; // Euilibrium distances for each type of constraint
            this->firstType     = firstType;     // The first type index used for constraints

            this->epbc = epbc;
            for (int i = 0; i < DIM; i++)
            {
                for (int j = 0; j < DIM; j++)
                {
                    this->box[i][j] = box[i][j]; // Periodic box
                }
            }

            set_pbc(&pbc, epbc, box);

            this->invdt  = 1.0/delta_t;
            this->wangle = LincsWarnAngle; // Change in angle between bond starting from which warning will be issued

            // Communication record
            cr.nnodes = 1;
            cr.dd     = nullptr;

            // Input record - data that usualy comes from configurtion file (.mdp)
            ir.nLincsIter     = nLincsIter;
            ir.nProjOrder     = nProjOrder;
            ir.efep           = 0;
            ir.LincsWarnAngle = LincsWarnAngle;
            ir.init_t         = init_t;
            ir.delta_t        = delta_t;
            ir.eI             = 0;

            // MD atoms data
            md.nMassPerturbed = 0;
            md.lambda         = 0.0;
            md.invmass        = invmass.data();
            md.nr             = n;
            md.homenr         = n;

            // No multisim
            ms = nullptr;

            // Constraints and their parameters in old data format
            idef.il[F_CONSTR].nr = constraints.size();

            snew(idef.il[F_CONSTR].iatoms, constraints.size());
            maxtype = 0;
            for (unsigned i = 0; i < constraints.size(); i++)
            {
                if (i % 3 == 0)
                {
                    if (maxtype < constraints.at(i))
                    {
                        maxtype = constraints.at(i);
                    }
                }
                idef.il[F_CONSTR].iatoms[i] = constraints.at(i);
            }
            snew(idef.iparams, maxtype + 1);
            for (unsigned i = 0; i < constraints.size()/3; i++)
            {
                idef.iparams[constraints.at(3*i)].constr.dA = constraintsR0.at(constraints.at(3*i) - firstType);
                idef.iparams[constraints.at(3*i)].constr.dB = constraintsR0.at(constraints.at(3*i) - firstType);
            }


            // Constraints and their parameters in new data format
            InteractionList ilist;
            ilist.iatoms.resize(constraints.size());
            for (unsigned i = 0; i < constraints.size(); i++)
            {
                ilist.iatoms.at(i) = constraints.at(i);
            }

            gmx_moltype_t moltype;
            moltype.atoms.nr           = n;
            moltype.ilist.at(F_CONSTR) = ilist;
            mtop.moltype.push_back(moltype);

            gmx_molblock_t molblock;
            molblock.type = 0;
            molblock.nmol = 1;
            mtop.molblock.push_back(molblock);

            mtop.natoms = n;
            mtop.ffparams.iparams.assign(idef.iparams, idef.iparams + (maxtype+1)*sizeof(t_iparams));
        }

        /*! \brief
         * Test for the final length of constraints and final coordinates.
         *
         * This method initializes and applies LINCS constraints. Two main test are performed on the results:
         * 1. If the final length of all the constraints is equal to the target lenght with provided tolerance.
         * 2. If the final coordinates correspond to the reference set of coordinates. This test is optional.
         *
         * \param[in] coordinates                Initial (unconstrained) coordinates.
         * \param[in] tolerance                  Allowed tolerance in final lengths.
         * \param[in] chechFinalCoordinates      Should final coordinates be tested againts reference coords.
         * \param[in] finalCoordinates           The reference set of coordinates.
         * \param[in] finalCoordinatesTolerance  Tolerance for the coordinaes test.
         */
        void testConstraints(const std::vector<RVec> coordinates, FloatingPointTolerance tolerance,
                             bool chechFinalCoordinates,
                             std::vector<RVec> finalCoordinates, FloatingPointTolerance finalCoordinatesTolerance)
        {
            // Set the number of omp threads
            gmx_omp_nthreads_set(emntLINCS, 1);

            // Set the coordinates in convinient format
            std::vector<RVec> x_vector(coordinates);
            std::vector<RVec> xprime_vector(coordinates);
            std::vector<RVec> min_proj_vector(coordinates);

            rvec             *x        = as_rvec_array(x_vector.data());
            rvec             *xprime   = as_rvec_array(xprime_vector.data());
            rvec             *min_proj = as_rvec_array(min_proj_vector.data());

            // Velocities (not used)
            rvec *v;
            snew(v, n);

            // Make blocka structure for faster LINCS setup
            std::vector<t_blocka> at2con_mt;
            at2con_mt.reserve(mtop.moltype.size());
            for (const gmx_moltype_t &moltype : mtop.moltype)
            {
                // This function is in constr.cpp
                at2con_mt.push_back(make_at2con(moltype,
                                                mtop.ffparams.iparams,
                                                flexibleConstraintTreatment(EI_DYNAMICS(ir.eI))));
            }

            // Log file may be here. Redirect it somewhere usefull?
            FILE* log = nullptr;

            // Disable free energy and virial computation
            real * dvdlambda     = nullptr;
            bool   computeVirial = false;
            tensor vir_r_m_dr;
            int    maxwarn         = 100;
            int    warncount_lincs = 0;

            // Initialie LINCS
            lincsd = init_lincs(log, mtop,
                                nflexcon, at2con_mt,
                                false,
                                ir.nLincsIter, ir.nProjOrder);

            set_lincs(idef, md, EI_DYNAMICS(ir.eI), &cr, lincsd);

            // Evaluate constraints
            bool bOK = constrain_lincs(false, ir, 0, lincsd, md,
                                       &cr,
                                       *ms,
                                       x, xprime, min_proj,
                                       box, &pbc, md.lambda, dvdlambda,
                                       invdt, v,
                                       computeVirial, vir_r_m_dr,
                                       gmx::ConstraintVariable::Positions, &nrnb,
                                       maxwarn, &warncount_lincs);
            EXPECT_TRUE(bOK) << "Something went wrong in LINCS: RMS of at least one atom is larger than 50%% of target length.";
            EXPECT_EQ(warncount_lincs, 0) << "There were warnings in LINCS.";

            // Test if all the constraints are satisfied
            for (unsigned c = 0; c < constraints.size()/3; c++)
            {
                real r0 = constraintsR0.at(constraints.at(3*c) - firstType);
                int  i  = constraints.at(3*c + 1);
                int  j  = constraints.at(3*c + 2);
                rvec v0, v1;
                real d0, d1;
                if (this->epbc == epbcXYZ)
                {
                    pbc_dx_aiuc(&pbc, x[i], x[j], v0);
                    pbc_dx_aiuc(&pbc, xprime[i], xprime[j], v1);
                }
                else
                {
                    rvec_sub(x[i], x[j], v0);
                    rvec_sub(xprime[i], xprime[j], v1);
                }
                d0 = norm(v0);
                d1 = norm(v1);
                EXPECT_REAL_EQ_TOL(d1, r0, tolerance) << "rij = " << d1 << ", which is not equal to r0 = " << r0
                << " for constraint #" << c << ", between atoms " << i << " and " << j
                << " (before lincs rij was " << d0 << ").";
            }

            // Check if coordinates after constraining correspond to the reference set
            if (chechFinalCoordinates)
            {
                for (int i = 0; i < n; i++)
                {
                    for (int d = 0; d < DIM; d++)
                    {
                        EXPECT_REAL_EQ_TOL(finalCoordinates[i][d], xprime[i][d], finalCoordinatesTolerance) <<
                        "Coordinates after LINCS constrains were applied differ from these in the reference set for atom #" << i << ".";
                    }
                }
            }

            sfree(v);
        }


        /*! \brief
         * Test for final length of the constraints.
         *
         * Same as previous, but only constrained length test is performed.
         *
         * \param[in] coordinates                Initial (unconstrained) coordinates.
         * \param[in] tolerance                  Allowed tolerance in final lengths.
         *
         */
        void testConstraints(const std::vector<RVec> coordinates, FloatingPointTolerance tolerance)
        {
            testConstraints(coordinates, tolerance, false, coordinates, tolerance);
        }

    private:
        int               n;
        gmx_mtop_t        mtop;
        Lincs            *lincsd;
        std::vector<real> invmass;
        int               epbc;
        matrix            box;
        t_pbc             pbc;
        t_commrec         cr;
        t_inputrec        ir;
        t_idef            idef;
        t_mdatoms         md;
        gmx_multisim_t   *ms;
        t_nrnb            nrnb;
        real              invdt    = 1.0/0.001;
        real              wangle   = 30.0;
        int               nflexcon = 0;
        int               maxtype  = 0;
        std::vector<int>  constraints;
        std::vector<real> constraintsR0;
        int               firstType;
};

/*
 *
 * Tests that the interatomic distance along constraints correspond to the reference distances.
 * Performed on basic systems of one, two and three constaints.
 *
 */

/*
 * For simple bond.
 */

TEST_F(LincsTest, ConstraintsOneBond)
{

    FloatingPointTolerance tolerance = relativeToleranceAsFloatingPoint(0.1, 0.00001);

    initSystem(2, c_oneBondMasses,
               c_oneBondConstraints, c_oneBondConstraintsR0, c_oneBondConstraintsFirstType,
               epbcNONE, c_infinitesimalBox,
               1, 4, real(30.0), real(0.0), real(0.001));
    std::vector<RVec> coordinates(c_oneBondCoordinates);
    real              distMod = 10.0; // Distortion
    coordinates[0][0] += 0.10*distMod;
    coordinates[0][1] -= 0.20*distMod;
    coordinates[0][2] += 0.05*distMod;
    coordinates[1][0] -= 0.20*distMod;
    coordinates[1][1] += 0.05*distMod;
    coordinates[1][2] -= 0.10*distMod;
    testConstraints(coordinates, tolerance);
}

/*
 * For two consequentive bonds.
 */

TEST_F(LincsTest, ConstraintsTwoBonds)
{

    FloatingPointTolerance tolerance = relativeToleranceAsFloatingPoint(0.1, 0.00001);

    initSystem(3, c_twoBondsMasses,
               c_twoBondsConstraints, c_twoBondsConstraintsR0, c_twoBondsConstraintsFirstType,
               epbcNONE, c_infinitesimalBox,
               1, 4, real(30.0), real(0.0), real(0.001));
    std::vector<RVec> coordinates(c_twoBondsCoordinates);
    real              distMod = 0.1; // Distortion
    coordinates[0][0] += 0.10*distMod;
    coordinates[0][1] -= 0.20*distMod;
    coordinates[0][2] += 0.05*distMod;
    coordinates[1][0] -= 0.20*distMod;
    coordinates[1][1] += 0.25*distMod;
    coordinates[1][2] -= 0.20*distMod;
    coordinates[2][0] += 0.25*distMod;
    coordinates[2][1] -= 0.10*distMod;
    coordinates[2][2] -= 0.05*distMod;
    testConstraints(coordinates, tolerance);
}

/*
 * For triangle of bonds.
 */
TEST_F(LincsTest, ConstraintsTriangle)
{

    FloatingPointTolerance tolerance = relativeToleranceAsFloatingPoint(0.1, 0.001);

    initSystem(3, c_triangleMasses,
               c_triangleConstraints, c_triangleConstraintsR0, c_triangleConstraintsFirstType,
               epbcNONE, c_infinitesimalBox,
               2, 8, real(30.0), real(0.0), real(0.001));
    std::vector<RVec> coordinates(c_triangleCoordinates);
    real              distMod = 0.05; // Distortion
    coordinates[0][0] += 0.10*distMod;
    coordinates[0][1] -= 0.20*distMod;
    coordinates[0][2] += 0.05*distMod;
    coordinates[1][0] -= 0.25*distMod;
    coordinates[1][1] += 0.20*distMod;
    coordinates[1][2] -= 0.25*distMod;
    coordinates[2][0] += 0.20*distMod;
    coordinates[2][1] -= 0.05*distMod;
    coordinates[2][2] -= 0.10*distMod;
    testConstraints(coordinates, tolerance);
}


/*
 *
 * Tests of both final interatomic distances and final coordinates.
 * Peformed on the peptide Cys-Val-Trp, with constraints along HBonds
 * or AllBonds.
 *
 */

/*
 * All bonds that involve hydrogen atom (HBonds) are constrained.
 */
TEST_F(LincsTest, ConstraintsCVWHBondsNone)
{

    FloatingPointTolerance tolerance            = relativeToleranceAsFloatingPoint(0.1, 0.00001);
    FloatingPointTolerance coordinatesTolerance = absoluteTolerance(0.000001);

    initSystem(54, c_cvwMasses,
               c_cvwHBondsConstraints, c_cvwHBondsConstraintsR0, c_cvwHBondsConstraintsFirstType,
               epbcXYZ, c_cvwBox,
               1, 4, real(30.0), real(0.0), real(0.001));
    testConstraints(c_cvwInitialCoordinates, tolerance,
                    true, c_cvwFinalCoordinatesHBonds, coordinatesTolerance);
}

/*
 * All bonds are constrained.
 */
TEST_F(LincsTest, ConstraintsCVWAllBondsNone)
{

    FloatingPointTolerance tolerance            = relativeToleranceAsFloatingPoint(0.1, 0.01);
    FloatingPointTolerance coordinatesTolerance = absoluteTolerance(0.000001);

    initSystem(54, c_cvwMasses,
               c_cvwAllBondsConstraints, c_cvwAllBondsConstraintsR0, c_cvwAllBondsConstraintsFirstType,
               epbcXYZ, c_cvwBox,
               1, 4, real(30.0), real(0.0), real(0.001));
    testConstraints(c_cvwInitialCoordinates, tolerance,
                    true, c_cvwFinalCoordinatesAllBonds, coordinatesTolerance);
}

} // namespace
} // namespace gmx
