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
/*! \internal \file
 * \brief SHAKE and LINCS tests.
 *
 * \author Artem Zhmurov <zhmurov@gmail.com>
 * \ingroup module_mdlib
 */

#include "gmxpre.h"

#include "gromacs/mdlib/constr.h"

#include <assert.h>

#include <cmath>

#include <algorithm>
#include <vector>

#include <gtest/gtest.h>

#include "gromacs/fileio/gmxfio.h"
#include "gromacs/gmxlib/nonbonded/nonbonded.h"
#include "gromacs/math/vec.h"
#include "gromacs/math/vectypes.h"
#include "gromacs/mdlib/gmx_omp_nthreads.h"
#include "gromacs/mdlib/lincs.h"
#include "gromacs/mdlib/shake.h"
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

#include "testutils/refdata.h"
#include "testutils/testasserts.h"

#include "constrdata.h"

namespace gmx
{
namespace test
{

/*! \brief Test fixture for SHAKE and LINCS
 *
 * The fixture uses three simple and one real-life system for testing.
 * Simple systems include:
 * 1.1. Two atoms, connected with one constrain.
 * 1.2. Three atoms, connected consequently with two constraints
 * 1.3. Four atoms, connected by two independent constraints.
 * 1.4. Three atoms, connected by three constraints in a triangle.
 * 1.5. Four atoms, connected by three consequentive constraints.
 * The real-life system is a Cys-Val-Trp peptide with constraints
 * on:
 * 2.1. Bonds, containing hydrogens (SHAKE and LINCS)
 * 2.2. All bonds (SHAKE and LINCS).
 * 2.3. Angles with hydrogens and all bonds (SHAKE).
 * 2.4. All angles and all bonds (disabled).
 *
 * For all systems, the final length of the constraints is tested ageinst the
 * reference values. For systems 2.1 and 2.2, the exact final coordinates of
 * all atoms are also tested against the reference values.
 */
class ConstraintsTest : public ::testing::Test
{
    public:

        /*! \brief
         * Cleaning up the memory.
         */
        void TearDown() override
        {
            sfree(idef.il[F_CONSTR].iatoms);
            sfree(idef.iparams);
            sfree(x);
            sfree(xprime);
            sfree(xprime2);
            sfree(v);
            gmx_fio_close(log);
        }

        /*! \brief
         * Method used to initialize atoms, constraints and their parameters.
         *
         * This method constructs stubs for all the data structures, required to initialize
         * and apply LINCS and SHAKE constraints. The constraints data is also saved as a field
         * to check if the target distances are achieved after the constraints are applied.
         *
         * \param[in] natom          Number of atoms in the system.
         * \param[in] masses         Atom masses. Size of this vector should be equal to natom.
         * \param[in] constraints    List of constraints, organized in tripples of integers.
         *                           First integer is the index of type for a constrain, second
         *                           and third are the indices of constrained atoms. The types
         *                           of constraints should be sequential but not nesseseraly
         *                           start from zero (which is the way they normally are in
         *                           GROMACS).
         * \param[in] constraintsR0  Target values for bond lengths for bonds of each type. The
         *                           size of this vector should be equal to the total number of
         *                           unique types in constraints vector.
         * \param[in] firstType      The index of first interaction type used for constraints
         *                           (i.e. the lowest number that appear as type in constraints
         *                           vector). This is not fail-safe, since technically the
         *                           numbering of the types may not be sequantial, but (i) this
         *                           is the way things are in GROMACS, (ii) everything is fine
         *                           in the examples provided for the test.
         * \param[in] epbc           Type of PBC (epbcXYZ / epbcNONE are used in this test,
         *                           see src/gromacs/pbcutil/pbc.h for details).
         * \param[in] box            The PBC box.
         * \param[in] init_t         Initial time.
         * \param[in] delta_t        Timestep.
         * \param[in] coordinates    Initial (unconstrained) coordinates.
         *
         */
        void initSystem(int natom, std::vector<real> masses,
                        std::vector<int> constraints, std::vector<real> constraintsR0, int firstType,
                        int epbc, const matrix box,
                        real init_t, real delta_t,
                        std::vector<RVec> coordinates)

        {

            this->n = natom;   // Number of atoms

            invmass.resize(n); // Vector of inverce masses

            for (int i = 0; i < n; i++)
            {
                invmass.at(i) = 1.0/masses.at(i);
            }

            // Saving constraints to check if they are satisfied after algorithm was applied
            this->constraints   = constraints;   // Constraints indexes (in type-i-j format)
            this->constraintsR0 = constraintsR0; // Euilibrium distances for each type of constraint
            this->firstType     = firstType;     // The first type index used for constraints

            // PBC initialization
            this->epbc = epbc;
            for (int i = 0; i < DIM; i++)
            {
                for (int j = 0; j < DIM; j++)
                {
                    this->box[i][j] = box[i][j]; // Periodic box
                }
            }
            set_pbc(&pbc, epbc, box);

            this->invdt  = 1.0/delta_t; // Inverse timestep

            // Communication record
            cr.nnodes = 1;
            cr.dd     = nullptr;

            // Input record - data that usualy comes from configurtion file (.mdp)
            ir.efep           = 0;
            ir.init_t         = init_t;
            ir.delta_t        = delta_t;
            ir.eI             = 0;

            // MD atoms data
            md.nMassPerturbed = 0;
            md.lambda         = 0.0;
            md.invmass        = invmass.data();
            md.nr             = n;
            md.homenr         = n;

            // Constraints and their parameters in old data format (local topology)
            for (int i = 0; i < F_NRE; i++)
            {
                idef.il[i].nr = 0;
            }
            idef.il[F_CONSTR].nr   = constraints.size();

            snew(idef.il[F_CONSTR].iatoms, constraints.size());
            int maxType = 0;
            for (unsigned i = 0; i < constraints.size(); i++)
            {
                if (i % 3 == 0)
                {
                    if (maxType < constraints.at(i))
                    {
                        maxType = constraints.at(i);
                    }
                }
                idef.il[F_CONSTR].iatoms[i] = constraints.at(i);
            }
            snew(idef.iparams, maxType + 1);
            for (unsigned i = 0; i < constraints.size()/3; i++)
            {
                idef.iparams[constraints.at(3*i)].constr.dA = constraintsR0.at(constraints.at(3*i) - firstType);
                idef.iparams[constraints.at(3*i)].constr.dB = constraintsR0.at(constraints.at(3*i) - firstType);
            }


            // Constraints and their parameters in new data format (global topology)
            InteractionList ilist;
            ilist.iatoms.resize(constraints.size());
            for (unsigned i = 0; i < constraints.size(); i++)
            {
                ilist.iatoms.at(i) = constraints.at(i);
            }
            InteractionList ilistEmpty;
            ilistEmpty.iatoms.resize(0);

            gmx_moltype_t moltype;
            moltype.atoms.nr             = n;
            moltype.ilist.at(F_CONSTR)   = ilist;
            moltype.ilist.at(F_CONSTRNC) = ilistEmpty;
            mtop.moltype.push_back(moltype);

            gmx_molblock_t molblock;
            molblock.type = 0;
            molblock.nmol = 1;
            mtop.molblock.push_back(molblock);

            mtop.natoms = n;
            mtop.ffparams.iparams.resize(maxType + 1);
            for (int i = 0; i <= maxType; i++)
            {
                mtop.ffparams.iparams.at(i) = idef.iparams[i];
            }
            mtop.bIntermolecularInteractions = false;

            // Log file may be here. Redirect it somewhere usefull?
            log = gmx_fio_open("constraintstest.log", "w");

            // Set the coordinates in convinient format
            snew(x, n);
            snew(xprime, n);
            snew(xprime2, n);

            for (int i = 0; i < n; i++)
            {
                for (int j = 0; j < DIM; j++)
                {
                    x[i][j]      = coordinates.at(i)[j];
                    xprime[i][j] = coordinates.at(i)[j];
                }
            }

            // Velocities (not used)
            snew(v, n);


        }

        /*! \brief
         * This method initializes and applies LINCS constraints.
         *
         * \param[in] nLincsIter     Number of iterations used to compute the inverse matrix.
         * \param[in] nProjOrder     The order for algorithm that adjusts the direction of the
         *                           bond after constraints are applied.
         * \param[in] LincsWarnAngle The value for the change in bond angle after which the
         *                           program will issue a warning.
         */
        void applyConstraintsLincs(int nLincsIter, int nProjOrder, real LincsWarnAngle)
        {

            Lincs                *lincsd;
            int                   maxwarn         = 100;
            int                   warncount_lincs = 0;
            ir.nLincsIter     = nLincsIter;
            ir.nProjOrder     = nProjOrder;
            ir.LincsWarnAngle = LincsWarnAngle;
            gmx_omp_nthreads_set(emntLINCS, 1);

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
            // Initialie LINCS
            lincsd = init_lincs(gmx_fio_getfp(log), mtop,
                                nflexcon, at2con_mt,
                                false,
                                ir.nLincsIter, ir.nProjOrder);
            set_lincs(idef, md, EI_DYNAMICS(ir.eI), &cr, lincsd);

            // Evaluate constraints
            bool bOK = constrain_lincs(false, ir, 0, lincsd, md,
                                       &cr,
                                       nullptr,
                                       x, xprime, xprime2,
                                       box, &pbc, md.lambda, dvdlambda,
                                       invdt, v,
                                       computeVirial, vir_r_m_dr,
                                       gmx::ConstraintVariable::Positions, &nrnb,
                                       maxwarn, &warncount_lincs);
            EXPECT_TRUE(bOK) << "Something went wrong: LINCS returned false value.";
            EXPECT_EQ(warncount_lincs, 0) << "There were warnings in LINCS.";
            for (unsigned int i = 0; i < mtop.moltype.size(); i++)
            {
                sfree(at2con_mt.at(i).index);
                sfree(at2con_mt.at(i).a);
            }
            done_lincs(lincsd);
        }

        /*! \brief
         * Initialize and apply SHAKE constraints.
         *
         * \param[in] shake_tol      Target tolerance for SHAKE.
         * \param[in] bShakeSOR      Use successive over-relaxation method for SHAKE iterative process.
         *                           The general formula is:
         *                              x_n+1 = (1-omega)*x_n + omega*f(x_n),
         *                           where omega = 1 if SOR is off and may be < 1 if SOR is on.
         *
         */
        void applyConstraintsShake(real shake_tol, gmx_bool bShakeSOR)
        {

            ir.shake_tol      = shake_tol;
            ir.bShakeSOR      = bShakeSOR;

            shakedata* shaked = shake_init();
            make_shake_sblock_serial(shaked, &idef, md);
            bool       bOK = constrain_shake(
                        gmx_fio_getfp(log),
                        shaked,
                        invmass.data(),
                        idef,
                        ir,
                        x,
                        xprime,
                        xprime2,
                        &nrnb,
                        md.lambda,
                        dvdlambda,
                        invdt,
                        v,
                        computeVirial,
                        vir_r_m_dr,
                        false,
                        gmx::ConstraintVariable::Positions);
            EXPECT_TRUE(bOK) << "Something went wrong: SHAKE returned false value.";
            done_shake(shaked);  // Not yet implemented. We will have memory leaks without this.
        }

        /*! \brief
         * The test on the final length of constrained bonds.
         *
         * Goes through all the constraints and checks if the final length of all the constraints is equal
         * to the target lenght with provided tolerance.
         *
         * \param[in] tolerance                  Allowed tolerance in final lengths.
         */

        void checkConstrained(FloatingPointTolerance tolerance)
        {

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
                EXPECT_REAL_EQ_TOL(r0, d1, tolerance) << "rij = " << d1 << ", which is not equal to r0 = " << r0
                << " for constraint #" << c << ", between atoms " << i << " and " << j
                << " (before lincs rij was " << d0 << ").";
            }
        }

        /*! \brief
         * The test on the final coordinates.
         *
         * Goes through all atoms and checks if the final positions correspond to the
         * provided reference set of coordinates.
         *
         * \param[in] finalCoordinates           The reference set of coordinates.
         * \param[in] finalCoordinatesTolerance  Tolerance for the coordinaes test.
         */
        void checkFinalCoordinates(std::vector<RVec> finalCoordinates, FloatingPointTolerance finalCoordinatesTolerance)
        {
            for (int i = 0; i < n; i++)
            {
                for (int d = 0; d < DIM; d++)
                {
                    EXPECT_REAL_EQ_TOL(finalCoordinates.at(i)[d], xprime[i][d], finalCoordinatesTolerance) <<
                    "Coordinates after constrains were applied differ from these in the reference set for atom #" << i << ".";
                }
            }
        }


    private:
        int                   n;                            // Number of atoms
        gmx_mtop_t            mtop;                         // Topology
        std::vector<real>     invmass;                      // Inverse masses
        int                   epbc;                         // PBC used (epbcNONE/epbcXYZ)
        matrix                box;                          // PBC box
        t_pbc                 pbc;                          // PBC object
        t_commrec             cr;                           // Communication record
        t_inputrec            ir;                           // Input record (info that usually in .mdp file)
        t_idef                idef;                         // Local topology
        t_mdatoms             md;                           // MD atoms
        t_nrnb                nrnb;

        real                  invdt         = 1.0/0.001;    // Inverse timestep
        int                   nflexcon      = 0;            // Number of flexible constraints
        real                 *dvdlambda     = nullptr;      // Free energy computation stub (not tested)
        bool                  computeVirial = false;        // If the virials should be conputed (not tested)
        tensor                vir_r_m_dr;                   // Virials stuff
        rvec                 *x;                            // Coordinates before constraints are applied
        rvec                 *xprime;                       // Coordinates after constraints are applied
        rvec                 *xprime2;                      // Intermediate set of coordinates used by LINCS and
                                                            // SHAKE for different purposes
        rvec                 *v;                            // Velocities (can also be constrained, this is not tested)
        t_fileio             *log;                          // log file (now set to "constraintstest.log")

        // Fields to store constraints data for testing
        int                   firstType;                    // First interaction type used for constraints
        std::vector<int>      constraints;                  // Constraints data (type1-i1-j1-type2-i2-j2-...)
        std::vector<real>     constraintsR0;                // Target lengths for al the constrain
};

/*
 *
 * Tests that the interatomic distance along constraints correspond to the reference distances.
 * Performed on basic systems of up to four constaints.
 *
 */

/*
 * Simple bond test for SHAKE
 */
TEST_F(ConstraintsTest, ShakeOneBond)
{

    FloatingPointTolerance tolerance = relativeToleranceAsFloatingPoint(0.1, 0.0001);
    initSystem(2, c_oneBondMasses,
               c_oneBondConstraints, c_oneBondConstraintsR0, c_oneBondConstraintsFirstType,
               epbcNONE, c_infinitesimalBox,
               real(0.0), real(0.001), c_oneBondCoordinates);
    applyConstraintsShake(0.0001, false);
    checkConstrained(tolerance);
}

/*
 * Simple bond test for LINCS
 */
TEST_F(ConstraintsTest, LincsOneBond)
{

    FloatingPointTolerance tolerance = relativeToleranceAsFloatingPoint(0.1, 0.0000002);
    initSystem(2, c_oneBondMasses,
               c_oneBondConstraints, c_oneBondConstraintsR0, c_oneBondConstraintsFirstType,
               epbcNONE, c_infinitesimalBox,
               real(0.0), real(0.001), c_oneBondCoordinates);
    applyConstraintsLincs(1, 4, real(30.0));
    checkConstrained(tolerance);
}

/*
 * Two disjoint bonds test for SHAKE
 */
TEST_F(ConstraintsTest, ShakeTwoDisjointBonds)
{

    FloatingPointTolerance tolerance = relativeToleranceAsFloatingPoint(0.1, 0.000001);
    initSystem(4, c_twoDJBondsMasses,
               c_twoDJBondsConstraints, c_twoDJBondsConstraintsR0, c_twoDJBondsConstraintsFirstType,
               epbcNONE, c_infinitesimalBox,
               real(0.0), real(0.001), c_twoDJBondsCoordinates);

    applyConstraintsShake(0.000001, true);
    checkConstrained(tolerance);
}

/*
 * Two disjoint bonds test for LINCS
 */
TEST_F(ConstraintsTest, LincsTwoDisjointBonds)
{

    FloatingPointTolerance tolerance = relativeToleranceAsFloatingPoint(0.1, 0.0000002);
    initSystem(4, c_twoDJBondsMasses,
               c_twoDJBondsConstraints, c_twoDJBondsConstraintsR0, c_twoDJBondsConstraintsFirstType,
               epbcNONE, c_infinitesimalBox,
               real(0.0), real(0.001), c_twoDJBondsCoordinates);

    applyConstraintsLincs(1, 4, real(30.0));
    checkConstrained(tolerance);
}

/*
 * Two consequentive constraints test for SHAKE.
 */
TEST_F(ConstraintsTest, ShakeTwoBonds)
{

    FloatingPointTolerance tolerance = relativeToleranceAsFloatingPoint(0.1, 0.0001);
    initSystem(3, c_twoBondsMasses,
               c_twoBondsConstraints, c_twoBondsConstraintsR0, c_twoBondsConstraintsFirstType,
               epbcNONE, c_infinitesimalBox,
               real(0.0), real(0.001), c_twoBondsCoordinates);

    applyConstraintsShake(0.0001, true);
    checkConstrained(tolerance);
}

/*
 * Two consequentive constraints test for LINCS.
 */
TEST_F(ConstraintsTest, LincsTwoBonds)
{

    FloatingPointTolerance tolerance = relativeToleranceAsFloatingPoint(0.1, 0.000001);
    initSystem(3, c_twoBondsMasses,
               c_twoBondsConstraints, c_twoBondsConstraintsR0, c_twoBondsConstraintsFirstType,
               epbcNONE, c_infinitesimalBox,
               real(0.0), real(0.001), c_twoBondsCoordinates);

    applyConstraintsLincs(1, 4, real(30.0));
    checkConstrained(tolerance);
}

/*
 * Triangle of constraints test for SHAKE
 */
TEST_F(ConstraintsTest, ShakeTriangleOfBonds)
{

    FloatingPointTolerance tolerance = relativeToleranceAsFloatingPoint(0.1, 0.0001);
    initSystem(3, c_triangleMasses,
               c_triangleConstraints, c_triangleConstraintsR0, c_triangleConstraintsFirstType,
               epbcNONE, c_infinitesimalBox,
               real(0.0), real(0.001), c_triangleCoordinates);
    applyConstraintsShake(0.0001, true);
    checkConstrained(tolerance);
}

/*
 * Triangle of constraints test for LINCS
 */
TEST_F(ConstraintsTest, LincsTriangleOfBonds)
{

    FloatingPointTolerance tolerance = relativeToleranceAsFloatingPoint(0.1, 0.0001);
    initSystem(3, c_triangleMasses,
               c_triangleConstraints, c_triangleConstraintsR0, c_triangleConstraintsFirstType,
               epbcNONE, c_infinitesimalBox,
               real(0.0), real(0.001), c_triangleCoordinates);
    applyConstraintsLincs(4, 4, real(30.0));
    checkConstrained(tolerance);
}

/*
 * Three consecutive constraints test for SHAKE.
 */
TEST_F(ConstraintsTest, ShakeThreeConsequativeConstraints)
{

    FloatingPointTolerance tolerance = relativeToleranceAsFloatingPoint(0.1, 0.001);
    initSystem(4, c_threeBondsMasses,
               c_threeBondsConstraints, c_threeBondsConstraintsR0, c_threeBondsConstraintsFirstType,
               epbcNONE, c_infinitesimalBox,
               real(0.0), real(0.001), c_threeBondsCoordinates);
    applyConstraintsShake(0.0001, true);
    checkConstrained(tolerance);
}

/*
 * Three consecutive constraints test for LINCS.
 */
TEST_F(ConstraintsTest, LincsThreeConsequativeConstraints)
{

    FloatingPointTolerance tolerance = relativeToleranceAsFloatingPoint(0.1, 0.001);
    initSystem(4, c_threeBondsMasses,
               c_threeBondsConstraints, c_threeBondsConstraintsR0, c_threeBondsConstraintsFirstType,
               epbcNONE, c_infinitesimalBox,
               real(0.0), real(0.001), c_threeBondsCoordinates);
    applyConstraintsLincs(4, 4, real(30.0));
    checkConstrained(tolerance);
}


/*
 *
 * Tests of both final interatomic distances and final coordinates.
 * Peformed on the peptide Cys-Val-Trp, with constraints along HBonds,
 * AllBonds (for both SHAKE and LINCS), HAngles and AllAlngles (for
 * SHAKE only).
 */

/*
 * All bonds that involve hydrogen atom (HBonds) are constrained by SHAKE. Both final length
 * of each constraint and final coordinates for each atom are tested.
 */
TEST_F(ConstraintsTest, ShakeCysValTrpHBonds)
{

    FloatingPointTolerance tolerance            = relativeToleranceAsFloatingPoint(0.1, 0.0001);
    FloatingPointTolerance coordinatesTolerance = absoluteTolerance(0.0001);

    initSystem(54, c_cvwMasses,
               c_cvwHBondsConstraints, c_cvwHBondsConstraintsR0, c_cvwHBondsConstraintsFirstType,
               epbcXYZ, c_cvwBox,
               real(0.0), real(0.001), c_cvwInitialCoordinates);
    applyConstraintsShake(0.0001, true);
    checkConstrained(tolerance);
    checkFinalCoordinates(c_cvwFinalCoordinatesHBonds, coordinatesTolerance);
}

/*
 * All bonds that involve hydrogen atom (HBonds) are constrained by LINCS. Both final length
 * of each constraint and final coordinates for each atom are tested.
 */
TEST_F(ConstraintsTest, LincsCysValTrpHBonds)
{

    FloatingPointTolerance tolerance            = relativeToleranceAsFloatingPoint(0.1, 0.00001);
    FloatingPointTolerance coordinatesTolerance = absoluteTolerance(0.0001);

    initSystem(54, c_cvwMasses,
               c_cvwHBondsConstraints, c_cvwHBondsConstraintsR0, c_cvwHBondsConstraintsFirstType,
               epbcXYZ, c_cvwBox,
               real(0.0), real(0.001), c_cvwInitialCoordinates);
    applyConstraintsLincs(1, 4, real(30.0));
    checkConstrained(tolerance);
    checkFinalCoordinates(c_cvwFinalCoordinatesHBonds, coordinatesTolerance);
}

/*
 * All bonds are constrained by SHAKE. Both final length of each constraint and final coordinates
 * for each atom are tested.
 */
TEST_F(ConstraintsTest, ShakeCysValTrpAllBonds)
{

    FloatingPointTolerance tolerance            = relativeToleranceAsFloatingPoint(0.1, 0.0001);
    FloatingPointTolerance coordinatesTolerance = absoluteTolerance(0.001);

    initSystem(54, c_cvwMasses,
               c_cvwAllBondsConstraints, c_cvwAllBondsConstraintsR0, c_cvwAllBondsConstraintsFirstType,
               epbcXYZ, c_cvwBox,
               real(0.0), real(0.001), c_cvwInitialCoordinates);
    applyConstraintsShake(0.0001, true);
    checkConstrained(tolerance);
    checkFinalCoordinates(c_cvwFinalCoordinatesAllBonds, coordinatesTolerance);
}

/*
 * All bonds are constrained by LINCS. Both final length of each constraint and final coordinates
 * for each atom are tested.
 */
TEST_F(ConstraintsTest, LincsCysValTrpAllBonds)
{

    FloatingPointTolerance tolerance            = relativeToleranceAsFloatingPoint(0.1, 0.0001);
    FloatingPointTolerance coordinatesTolerance = absoluteTolerance(0.001);

    initSystem(54, c_cvwMasses,
               c_cvwAllBondsConstraints, c_cvwAllBondsConstraintsR0, c_cvwAllBondsConstraintsFirstType,
               epbcXYZ, c_cvwBox,
               real(0.0), real(0.001), c_cvwInitialCoordinates);
    applyConstraintsLincs(4, 4, real(30.0));
    checkConstrained(tolerance);
    checkFinalCoordinates(c_cvwFinalCoordinatesAllBonds, coordinatesTolerance);
}

/*
 * H angles and all bonds are constrained by SHAKE. Only final lengths of constraints are tested.
 */
TEST_F(ConstraintsTest, ShakeCysValTrpHAngles)
{

    FloatingPointTolerance tolerance            = relativeToleranceAsFloatingPoint(0.1, 0.0001);

    initSystem(54, c_cvwMasses,
               c_cvwHAnglesConstraints, c_cvwHAnglesConstraintsR0, c_cvwHAnglesConstraintsFirstType,
               epbcXYZ, c_cvwBox,
               real(0.0), real(0.001), c_cvwInitialCoordinates);
    applyConstraintsShake(0.0001, true);
    checkConstrained(tolerance);
}

/*
 * All angles and all bonds are constrained by SHAKE. Only final lengths of constraints are tested.
 */
/*TEST_F(ConstraintsTest, ShakeCysValTrpAllAngles)
   {

    FloatingPointTolerance tolerance            = relativeToleranceAsFloatingPoint(0.1, 0.001);

    initSystem(54, c_cvwMasses,
               c_cvwAllAnglesConstraints, c_cvwAllAnglesConstraintsR0, c_cvwAllAnglesConstraintsFirstType,
               epbcXYZ, c_cvwBox,
               real(0.0), real(0.001), c_cvwInitialCoordinates);
    applyConstraintsShake(0.000001, true);
    checkConstrained(tolerance);
   }
 */
} // namespace test
} // namespace gmx
