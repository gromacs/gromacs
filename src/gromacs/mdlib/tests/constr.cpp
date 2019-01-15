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

namespace gmx
{
namespace test
{

//! 0.1/sqrt(2), used to set coordinates
static constexpr real ONE_TENTH_OVER_SQRT_TWO    = 0.0707106781;
//! 0.1/sqrt(3), used to set coordinates
static constexpr real TWO_TENTHS_OVER_SQRT_THREE = 0.1154700538;

/*! \brief
 * Constraints test data structure.
 *
 * Structure to collect all the necessary data, including system coordinates and topology,
 * constraints information, etc. The structure can be reseted and reused.
 */
struct ConstraintsTestData
{
    public:
        std::string           title;                   //!< Human-friendly name for a system
        int                   nAtom;                   //!< Number of atoms
        gmx_mtop_t            mtop;                    //!< Topology
        std::vector<real>     masses;                  //!< Masses
        std::vector<real>     invmass;                 //!< Inverse masses
        t_commrec             cr;                      //!< Communication record
        t_inputrec            ir;                      //!< Input record (info that usually in .mdp file)
        t_idef                idef;                    //!< Local topology
        t_mdatoms             md;                      //!< MD atoms
        gmx_multisim_t        ms;                      //!< Multisim data
        t_nrnb                nrnb;                    //!< Computational time array (normally used for benchmarking)

        real                  invdt      = 1.0/0.001;  //!< Inverse timestep
        int                   nflexcon   = 0;          //!< Number of flexible constraints
        bool                  computeVirial;           //!< Whether the virial should be computed
        tensor                virialScaled;            //!< Scaled virial
        tensor                virialScaledRef;         //!< Scaled virial (reference values)
        bool                  compute_dHdLambda;       //!< If the free energy is computed
        real                  dHdLambda;               //!< For free energy computation
        real                  dHdLambdaRef;            //!< For free energy computation (reference value)

        std::vector<RVec>     x;                       //!< Coordinates before the timestep
        std::vector<RVec>     xPrime;                  //!< Coordinates after timestep, output for the constraints
        std::vector<RVec>     xPrime0;                 //!< Backup for coordinates (for reseting)
        std::vector<RVec>     xPrime2;                 //!< Intermediate set of coordinates used by LINCS and
                                                       //!< SHAKE for different purposes
        std::vector<RVec>     v;                       //!< Velocities
        std::vector<RVec>     v0;                      //!< Backup for velocities (for reseting)

        t_fileio             *log;                     //!< log file (now set to "constraintstest.log")

        // Fields to store constraints data for testing
        std::vector<int>      constraints;             //!< Constraints data (type1-i1-j1-type2-i2-j2-...)
        std::vector<real>     constraintsR0;           //!< Target lengths for al the constrain

        // Tolerances in final values
        real                  shakeExpectedTolerance;  //!< Expected tolerance in values after SHAKE is applied.
        real                  lincsExpectedTolerance;  //!< Expected tolerance in values after LINCS is applied.

        /*! \brief
         * Constructor for the object with all parameters and variables needed to constraints algorithms.
         *
         * This constructor assembles stubs for all the data structures, required to initialize
         * and apply LINCS and SHAKE constraints. The coordinates and velocities before constraining
         * are saved to allow for reseting. The constraints data are stored for testing after constraints
         * were applied.
         *
         * \param[in]  title              Human-friendly name of the system.
         * \param[in]  nAtom              Number of atoms in the system.
         * \param[in]  masses             Atom masses. Size of this vector should be equal to nAtom.
         * \param[in]  constraints        List of constraints, organized in triples of integers.
         *                                First integer is the index of type for a constrain, second
         *                                and third are the indices of constrained atoms. The types
         *                                of constraints should be sequential but not necessarily
         *                                start from zero (which is the way they normally are in
         *                                GROMACS).
         * \param[in]  constraintsR0      Target values for bond lengths for bonds of each type. The
         *                                size of this vector should be equal to the total number of
         *                                unique types in constraints vector.
         * \param[in]  computeVirial      Whether the virial should be computed.
         * \param[in]  virialScaledRef    Reference values for scaled virial tensor.
         * \param[in]  compute_dHdLambda  Whether free energy should be computed.
         * \param[in]  dHdLambdaRef       Reference value for dHdLambda.
         * \param[in]  initialTime        Initial time.
         * \param[in]  timestep           Timestep.
         * \param[in]  x                  Coordinates before integration step.
         * \param[in]  xPrime             Coordinates after integration step, but before constraining.
         * \param[in]  v                  Velocities before constraining.
         * \param[in]  shakeTolerance     Target tolerance for SHAKE.
         * \param[in]  shakeUseSOR        Use successive over-relaxation method for SHAKE iterations.
         *                                The general formula is:
         *                                   x_n+1 = (1-omega)*x_n + omega*f(x_n),
         *                                where omega = 1 if SOR is off and may be < 1 if SOR is on.
         * \param[in]  shakeExpectedTolerance  Tolerance in return values for SHAKE.
         * \param[in]  nLincsIter         Number of iterations used to compute the inverse matrix.
         * \param[in]  nProjOrder         The order for algorithm that adjusts the direction of the
         *                                bond after constraints are applied.
         * \param[in]  lincsWarnAngle     The threshold value for the change in bond angle. When
         *                                exceeded the program will issue a warning.
         * \param[in]  lincsExpectedTolerance  Tolerance in return values for LINCS.
         *
         */
        ConstraintsTestData(const std::string &title,
                            int nAtom, std::vector<real> masses,
                            std::vector<int> constraints, std::vector<real> constraintsR0,
                            bool computeVirial, tensor virialScaledRef,
                            bool compute_dHdLambda, float dHdLambdaRef,
                            real initialTime, real timestep,
                            const std::vector<RVec> &x, const std::vector<RVec> &xPrime, const std::vector<RVec> &v,
                            real shakeTolerance, gmx_bool shakeUseSOR, real shakeExpectedTolerance,
                            int nLincsIter, int nProjOrder, real lincsWarnAngle, real lincsExpectedTolerance)
        {
            this->title = title;
            this->nAtom = nAtom;   // Number of atoms

            this->masses = masses;
            invmass.resize(nAtom); // Vector of inverse masses

            for (int i = 0; i < nAtom; i++)
            {
                invmass.at(i) = 1.0/masses.at(i);
            }

            // Saving constraints to check if they are satisfied after algorithm was applied
            this->constraints   = constraints;   // Constraints indexes (in type-i-j format)
            this->constraintsR0 = constraintsR0; // Equilibrium distances for each type of constraint

            this->invdt  = 1.0/timestep;         // Inverse timestep

            // Communication record
            cr.nnodes = 1;
            cr.dd     = nullptr;

            // Multisim data
            ms.sim  = 0;
            ms.nsim = 1;

            // Input record - data that usually comes from configuration file (.mdp)
            ir.efep           = 0;
            ir.init_t         = initialTime;
            ir.delta_t        = timestep;
            ir.eI             = 0;

            // MD atoms data
            md.nMassPerturbed = 0;
            md.lambda         = 0.0;
            md.invmass        = invmass.data();
            md.nr             = nAtom;
            md.homenr         = nAtom;

            // Virial evaluation
            this->computeVirial = computeVirial;
            if (computeVirial)
            {
                for (int i = 0; i < DIM; i++)
                {
                    for (int j = 0; j < DIM; j++)
                    {
                        this->virialScaled[i][j]    = 0;
                        this->virialScaledRef[i][j] = virialScaledRef[i][j];
                    }
                }
            }


            // Free energy evaluation
            this->compute_dHdLambda = compute_dHdLambda;
            this->dHdLambda         = 0;
            if (compute_dHdLambda)
            {
                ir.efep                    = efepYES;
                this->dHdLambdaRef         = dHdLambdaRef;
            }
            else
            {
                ir.efep                    = efepNO;
                this->dHdLambdaRef         = 0;
            }

            // Constraints and their parameters (local topology)
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
                idef.iparams[constraints.at(3*i)].constr.dA = constraintsR0.at(constraints.at(3*i));
                idef.iparams[constraints.at(3*i)].constr.dB = constraintsR0.at(constraints.at(3*i));
            }

            // Constraints and their parameters (global topology)
            InteractionList interactionList;
            interactionList.iatoms.resize(constraints.size());
            for (unsigned i = 0; i < constraints.size(); i++)
            {
                interactionList.iatoms.at(i) = constraints.at(i);
            }
            InteractionList interactionListEmpty;
            interactionListEmpty.iatoms.resize(0);

            gmx_moltype_t molType;
            molType.atoms.nr             = nAtom;
            molType.ilist.at(F_CONSTR)   = interactionList;
            molType.ilist.at(F_CONSTRNC) = interactionListEmpty;
            mtop.moltype.push_back(molType);

            gmx_molblock_t molBlock;
            molBlock.type = 0;
            molBlock.nmol = 1;
            mtop.molblock.push_back(molBlock);

            mtop.natoms = nAtom;
            mtop.ffparams.iparams.resize(maxType + 1);
            for (int i = 0; i <= maxType; i++)
            {
                mtop.ffparams.iparams.at(i) = idef.iparams[i];
            }
            mtop.bIntermolecularInteractions = false;

            // Log file may be here. Redirect it somewhere useful?
            log = gmx_fio_open("constraintstest.log", "w");

            // Coordinates and velocities
            this->x       = x;
            this->xPrime  = xPrime;
            this->xPrime0 = xPrime;
            this->xPrime2 = xPrime;

            this->v  = v;
            this->v0 = v;

            // SHAKE-specific parameters
            this->ir.shake_tol           = shakeTolerance;
            this->ir.bShakeSOR           = shakeUseSOR;
            this->shakeExpectedTolerance = shakeExpectedTolerance;

            // LINCS-specific parameters
            this->ir.nLincsIter          = nLincsIter;
            this->ir.nProjOrder          = nProjOrder;
            this->ir.LincsWarnAngle      = lincsWarnAngle;
            this->lincsExpectedTolerance = lincsExpectedTolerance;
        }

        /*! \brief
         * Reset the data structure so it can be reused.
         *
         * Set the coordinates and velocities back to their values before
         * constraining. The scaled virial tensor and dHdLambda are zeroed.
         *
         */
        void reset()
        {
            xPrime  = xPrime0;
            xPrime2 = xPrime0;
            v       = v0;

            if (computeVirial)
            {
                for (int i = 0; i < DIM; i++)
                {
                    for (int j = 0; j < DIM; j++)
                    {
                        this->virialScaled[i][j] = 0;
                    }
                }
            }
            this->dHdLambda         = 0;
        }

        /*! \brief
         * Cleaning up the memory.
         */
        ~ConstraintsTestData()
        {
            sfree(idef.il[F_CONSTR].iatoms);
            sfree(idef.iparams);
            gmx_fio_close(log);
        }


};

/*! \brief The three-dimensional parameter space for test.
 *
 * The test will run for all possible combinations of accessible
 * values of the following parameters:
 * 1. PBC setup.
 * 2. System to test on.
 * 3. The algorithm (SHAKE or LINCS).
 */
typedef std::tuple<int, int, int> ConstraintsTestParameters;

/*! \brief Test fixture for SHAKE and LINCS
 *
 * The fixture uses following test systems:
 * 1. Two atoms, connected with one constraint (e.g. NH).
 * 2. Three atoms, connected consequently with two constraints (e.g. CH2).
 * 3. Three atoms, constrained to the fourth atom (e.g. CH3).
 * 4. Four atoms, connected by two independent constraints.
 * 5. Three atoms, connected by three constraints in a triangle
 *      (e.g. H2O with constrained H-O-H angle).
 * 6. Four atoms, connected by three consequential constraints.
 *
 * For all systems, the final lengths of the constraints are tested against the
 * reference values, the direction of each constraint is checked.
 * Test also verifies that the center of mass has not been
 * shifted by the constraints and that its velocity has not changed.
 * For some systems, the value for scaled virial tensor is checked against
 * pre-computed data.
 */
class ConstraintsTest : public ::testing::TestWithParam<ConstraintsTestParameters>
{
    public:

        std::vector<ConstraintsTestData*> systems; //!< A set of system on which constraints will be tested
        std::vector<t_pbc>                pbcs;    //!< PBC setups

        /*! \brief Test setup function.
         *
         * Fills in the systems and pbcs vectors.
         */
        void SetUp() override
        {

            //
            // PBC initialization
            //
            t_pbc pbc;

            // Infinitely small box
            matrix boxNone = { {0, 0, 0}, {0, 0, 0}, {0, 0, 0} };
            //matrix_convert(box, RVec({ 0.0, 0.0, 0.0 }), RVec({ 90.0, 90.0, 90.0 }));
            set_pbc(&pbc, epbcNONE, boxNone);
            pbcs.push_back(pbc);

            // Rectangular box
            matrix boxXyz = { {10.0, 0.0, 0.0}, {0.0, 20.0, 0.0}, {0.0, 0.0, 15.0} };
            //matrix_convert(box, RVec({ 10.0, 20.0, 15.0 }), RVec({ 90.0, 90.0, 90.0 }));
            set_pbc(&pbc, epbcXYZ, boxXyz);
            pbcs.push_back(pbc);

            //
            // Systems
            //
            std::string       title;
            int               nAtom;

            std::vector<real> masses;
            std::vector<int>  constraints;
            std::vector<real> constraintsR0;

            std::vector<RVec> x;
            std::vector<RVec> xPrime;
            std::vector<RVec> v;

            real              defaultShakeTolerance         = 0.0001;
            gmx_bool          defaultShakeUseSOR            = false;
            real              defaultShakeExpectedTolerance = 0.0002;

            int               defaultNLincsIter             = 1;
            int               defaultLincsNProjOrder        = 4;
            real              defaultLincsWarnAngle         = 30.0;
            real              defaultLincsExpectedTolerance = 0.0002;


            tensor emptyTensor = {{0, 0, 0}, {0, 0, 0}, {0, 0, 0}};

            //
            //
            title         = "one constraint (e.g. OH)";
            nAtom         = 2;
            masses        = std::vector<real>( {1.0, 12.0} );
            constraints   = std::vector<int>( {0, 0, 1} );
            constraintsR0 = std::vector<real>( {0.1} );


            x = std::vector<RVec>({{ { 0.0, ONE_TENTH_OVER_SQRT_TWO, 0.0 },
                                     { ONE_TENTH_OVER_SQRT_TWO, 0.0, 0.0 } }});

            xPrime = std::vector<RVec>({{ { 0.01,  0.08,  0.01 },
                                          { 0.06,  0.01, -0.01 } }});

            v = std::vector<RVec>({{ { 1.0, 2.0, 3.0 },
                                     { 3.0, 2.0, 1.0 } }});

            tensor virialScaledRef = { {-5.58e-04,  5.58e-04, 0.00e+00 },
                                       { 5.58e-04, -5.58e-04, 0.00e+00 },
                                       { 0.00e+00,  0.00e+00, 0.00e+00 } };

            systems.push_back(new ConstraintsTestData(
                                      title, nAtom, masses,
                                      constraints, constraintsR0,
                                      true, virialScaledRef,
                                      false, 0,
                                      real(0.0), real(0.001),
                                      x, xPrime, v,
                                      defaultShakeTolerance, defaultShakeUseSOR,
                                      defaultShakeExpectedTolerance,
                                      defaultNLincsIter, defaultLincsNProjOrder, defaultLincsWarnAngle,
                                      defaultLincsExpectedTolerance));

            //
            //
            title         = "two disjoint constraints";
            nAtom         = 4;
            masses        = std::vector<real>( {0.5, 1.0/3.0, 0.25, 1.0} );
            constraints   = std::vector<int>( {0, 0, 1, 1, 2, 3} );
            constraintsR0 = std::vector<real>( {2.0, 1.0} );


            x = std::vector<RVec>({{ {  2.50, -3.10, 15.70 },
                                     {  0.51, -3.02, 15.55 },
                                     { -0.50, -3.00, 15.20 },
                                     { -1.51, -2.95, 15.05 } }});

            xPrime = std::vector<RVec>({{ {  2.50, -3.10, 15.70 },
                                          {  0.51, -3.02, 15.55 },
                                          { -0.50, -3.00, 15.20 },
                                          { -1.51, -2.95, 15.05 } }});

            v = std::vector<RVec>({{ { 0.0, 1.0, 0.0 },
                                     { 1.0, 0.0, 0.0 },
                                     { 0.0, 0.0, 1.0 },
                                     { 0.0, 0.0, 0.0 } }});
            systems.push_back(new ConstraintsTestData(
                                      title, nAtom, masses,
                                      constraints, constraintsR0,
                                      false, emptyTensor,
                                      false, 0,
                                      real(0.0), real(0.001),
                                      x, xPrime, v,
                                      defaultShakeTolerance, defaultShakeUseSOR,
                                      defaultShakeExpectedTolerance,
                                      defaultNLincsIter, defaultLincsNProjOrder, defaultLincsWarnAngle,
                                      defaultLincsExpectedTolerance));

            //
            //
            title         = "three atoms, connected longitudaly (e.g. CH2)";
            nAtom         = 3;
            masses        = std::vector<real>( {1.0, 12.0, 16.0 } );
            constraints   = std::vector<int>( {0, 0, 1, 1, 1, 2} );
            constraintsR0 = std::vector<real>( {0.1, 0.2} );


            x = std::vector<RVec>({{ { ONE_TENTH_OVER_SQRT_TWO, ONE_TENTH_OVER_SQRT_TWO, 0.0 },
                                     { 0.0, 0.0, 0.0 },
                                     { TWO_TENTHS_OVER_SQRT_THREE, TWO_TENTHS_OVER_SQRT_THREE, TWO_TENTHS_OVER_SQRT_THREE } }});

            xPrime = std::vector<RVec>({{ {  0.08, 0.07,  0.01 },
                                          { -0.02, 0.01, -0.02 },
                                          {  0.10, 0.12,  0.11 } }});

            v = std::vector<RVec>({{ { 1.0, 0.0, 0.0 },
                                     { 0.0, 1.0, 0.0 },
                                     { 0.0, 0.0, 1.0 } }});
            systems.push_back(new ConstraintsTestData(
                                      title, nAtom, masses,
                                      constraints, constraintsR0,
                                      false, emptyTensor,
                                      false, 0,
                                      real(0.0), real(0.001),
                                      x, xPrime, v,
                                      defaultShakeTolerance, defaultShakeUseSOR,
                                      defaultShakeExpectedTolerance,
                                      defaultNLincsIter, defaultLincsNProjOrder, defaultLincsWarnAngle,
                                      defaultLincsExpectedTolerance));

            //
            //
            title         = "three atoms, connected to the central atom (e.g. CH3)";
            nAtom         = 4;
            masses        = std::vector<real>( {12.0, 1.0, 1.0, 1.0 } );
            constraints   = std::vector<int>( {0, 0, 1, 0, 0, 2, 0, 0, 3} );
            constraintsR0 = std::vector<real>( {0.1} );


            x = std::vector<RVec>({{  { 0.00,  0.00,  0.00 },
                                      { 0.10,  0.00,  0.00 },
                                      { 0.00, -0.10,  0.00 },
                                      { 0.00,  0.00,  0.10 } }});

            xPrime = std::vector<RVec>({{ { 0.004,  0.009, -0.010 },
                                          { 0.110, -0.006,  0.003 },
                                          {-0.007, -0.102, -0.007 },
                                          {-0.005,  0.011,  0.102 } }});

            v = std::vector<RVec>({{ { 1.0, 0.0, 0.0 },
                                     { 1.0, 0.0, 0.0 },
                                     { 1.0, 0.0, 0.0 },
                                     { 1.0, 0.0, 0.0 } }});
            systems.push_back(new ConstraintsTestData(
                                      title, nAtom, masses,
                                      constraints, constraintsR0,
                                      false, emptyTensor,
                                      false, 0,
                                      real(0.0), real(0.001),
                                      x, xPrime, v,
                                      defaultShakeTolerance, defaultShakeUseSOR,
                                      defaultShakeExpectedTolerance,
                                      defaultNLincsIter, defaultLincsNProjOrder, defaultLincsWarnAngle,
                                      defaultLincsExpectedTolerance));

            //
            //
            title         = "four atoms, connected longitudaly";
            nAtom         = 4;
            masses        = std::vector<real>( {0.5, 1.0/3.0, 0.25, 1.0} );
            constraints   = std::vector<int>( {0, 0, 1, 1, 1, 2, 2, 2, 3} );
            constraintsR0 = std::vector<real>( {2.0, 1.0, 1.0} );


            x = std::vector<RVec>({{ {  2.50, -3.10, 15.70 },
                                     {  0.51, -3.02, 15.55 },
                                     { -0.50, -3.00, 15.20 },
                                     { -1.51, -2.95, 15.05 } }});

            xPrime = std::vector<RVec>({{ {  2.50, -3.10, 15.70 },
                                          {  0.51, -3.02, 15.55 },
                                          { -0.50, -3.00, 15.20 },
                                          { -1.51, -2.95, 15.05 } }});

            v = std::vector<RVec>({{ { 0.0, 0.0,  2.0 },
                                     { 0.0, 0.0,  3.0 },
                                     { 0.0, 0.0, -4.0 },
                                     { 0.0, 0.0, -1.0 } }});
            systems.push_back(new ConstraintsTestData(
                                      title, nAtom, masses,
                                      constraints, constraintsR0,
                                      false, emptyTensor,
                                      false, 0,
                                      real(0.0), real(0.001),
                                      x, xPrime, v,
                                      defaultShakeTolerance, defaultShakeUseSOR,
                                      defaultShakeExpectedTolerance,
                                      2, 8, defaultLincsWarnAngle,
                                      0.002));
            //
            //
            title         = "basic triangle (tree atoms, connected with each other)";
            nAtom         = 3;
            masses        = std::vector<real>( {1.0, 1.0, 1.0} );
            constraints   = std::vector<int>( {0, 0, 1, 2, 0, 2, 1, 1, 2} );
            constraintsR0 = std::vector<real>( {0.1, 0.1, 0.1} );


            x = std::vector<RVec>({{ { ONE_TENTH_OVER_SQRT_TWO, 0.0, 0.0 },
                                     { 0.0, ONE_TENTH_OVER_SQRT_TWO, 0.0 },
                                     { 0.0, 0.0, ONE_TENTH_OVER_SQRT_TWO } }});

            xPrime = std::vector<RVec>({{ {  0.09, -0.02,  0.01 },
                                          { -0.02,  0.10, -0.02 },
                                          {  0.03, -0.01,  0.07 } }});

            v = std::vector<RVec>({{ {  1.0,  1.0,  1.0 },
                                     { -2.0, -2.0, -2.0 },
                                     {  1.0,  1.0,  1.0 } }});
            systems.push_back(new ConstraintsTestData(
                                      title, nAtom, masses,
                                      constraints, constraintsR0,
                                      false, emptyTensor,
                                      false, 0,
                                      real(0.0), real(0.001),
                                      x, xPrime, v,
                                      defaultShakeTolerance, defaultShakeUseSOR,
                                      defaultShakeExpectedTolerance,
                                      defaultNLincsIter, defaultLincsNProjOrder, defaultLincsWarnAngle,
                                      defaultLincsExpectedTolerance));

        }

        /*! \brief
         * Cleaning up the memory.
         */
        void TearDown() override
        {
            for (auto &system : systems)
            {
                delete(system);
            }
        }

        /*! \brief
         * Initialize and apply SHAKE constraints.
         *
         * \param[in] testData        Test data structure.
         * \param[in] testDescription Human-friendly message with test parameters.
         */
        void applyShake(ConstraintsTestData *testData, const std::string &testDescription)
        {
            shakedata* shaked = shake_init();
            make_shake_sblock_serial(shaked, &testData->idef, testData->md);
            bool       success = constrain_shake(
                        gmx_fio_getfp(testData->log),
                        shaked,
                        testData->invmass.data(),
                        testData->idef,
                        testData->ir,
                        as_rvec_array(testData->x.data()),
                        as_rvec_array(testData->xPrime.data()),
                        as_rvec_array(testData->xPrime2.data()),
                        &testData->nrnb,
                        testData->md.lambda,
                        &testData->dHdLambda,
                        testData->invdt,
                        as_rvec_array(testData->v.data()),
                        testData->computeVirial,
                        testData->virialScaled,
                        false,
                        gmx::ConstraintVariable::Positions);
            EXPECT_TRUE(success) << "Something went wrong: SHAKE returned false value. " << testDescription;
            done_shake(shaked);
        }

        /*! \brief
         * Initialize and apply LINCS constraints.
         *
         * \param[in] testData        Test data structure.
         * \param[in] pbc             Periodic boundary data.
         * \param[in] testDescription Human-friendly message with test parameters.
         */
        void applyLincs(ConstraintsTestData *testData, t_pbc pbc, const std::string &testDescription)
        {

            Lincs                *lincsd;
            int                   maxwarn         = 100;
            int                   warncount_lincs = 0;
            gmx_omp_nthreads_set(emntLINCS, 1);

            // Make blocka structure for faster LINCS setup
            std::vector<t_blocka> at2con_mt;
            at2con_mt.reserve(testData->mtop.moltype.size());
            for (const gmx_moltype_t &moltype : testData->mtop.moltype)
            {
                // This function is in constr.cpp
                at2con_mt.push_back(make_at2con(moltype,
                                                testData->mtop.ffparams.iparams,
                                                flexibleConstraintTreatment(EI_DYNAMICS(testData->ir.eI))));
            }
            // Initialize LINCS
            lincsd = init_lincs(gmx_fio_getfp(testData->log), testData->mtop,
                                testData->nflexcon, at2con_mt,
                                false,
                                testData->ir.nLincsIter, testData->ir.nProjOrder);
            set_lincs(testData->idef, testData->md, EI_DYNAMICS(testData->ir.eI), &testData->cr, lincsd);

            // Evaluate constraints
            bool sucess = constrain_lincs(false, testData->ir, 0, lincsd, testData->md,
                                          &testData->cr,
                                          &testData->ms,
                                          as_rvec_array(testData->x.data()),
                                          as_rvec_array(testData->xPrime.data()),
                                          as_rvec_array(testData->xPrime2.data()),
                                          pbc.box, &pbc, testData->md.lambda, &testData->dHdLambda,
                                          testData->invdt,
                                          as_rvec_array(testData->v.data()),
                                          testData->computeVirial, testData->virialScaled,
                                          gmx::ConstraintVariable::Positions, &testData->nrnb,
                                          maxwarn, &warncount_lincs);
            EXPECT_TRUE(sucess) << "Something went wrong: LINCS returned false value. " << testDescription;
            EXPECT_EQ(warncount_lincs, 0) << "There were warnings in LINCS. " << testDescription;
            for (unsigned int i = 0; i < testData->mtop.moltype.size(); i++)
            {
                sfree(at2con_mt.at(i).index);
                sfree(at2con_mt.at(i).a);
            }
            done_lincs(lincsd);
        }

        /*! \brief
         * The test on the final length of constrained bonds.
         *
         * Goes through all the constraints and checks if the final length of all the constraints is equal
         * to the target length with provided tolerance.
         *
         * \param[in] tolerance       Allowed tolerance in final lengths.
         * \param[in] testData        Test data structure.
         * \param[in] pbc             Periodic boundary data.
         * \param[in] testDescription Human-friendly message with test parameters.
         */
        void checkConstrainsLength(FloatingPointTolerance tolerance, ConstraintsTestData *testData, t_pbc pbc,
                                   const std::string &testDescription)
        {

            // Test if all the constraints are satisfied
            for (unsigned c = 0; c < testData->constraints.size()/3; c++)
            {
                real r0 = testData->constraintsR0.at(testData->constraints.at(3*c));
                int  i  = testData->constraints.at(3*c + 1);
                int  j  = testData->constraints.at(3*c + 2);
                RVec xij0, xij1;
                real d0, d1;
                if (pbc.ePBC == epbcXYZ)
                {
                    pbc_dx_aiuc(&pbc, testData->x[i], testData->x[j], xij0);
                    pbc_dx_aiuc(&pbc, testData->xPrime[i], testData->xPrime[j], xij1);
                }
                else
                {
                    rvec_sub(testData->x[i], testData->x[j], xij0);
                    rvec_sub(testData->xPrime[i], testData->xPrime[j], xij1);
                }
                d0 = norm(xij0);
                d1 = norm(xij1);
                EXPECT_REAL_EQ_TOL(r0, d1, tolerance) << "rij = " << d1 << ", which is not equal to r0 = " << r0
                << " for constraint #" << c << ", between atoms " << i << " and " << j
                << " (before constraining rij was " << d0 << "). " << testDescription;
            }
        }

        /*! \brief
         * The test on the final length of constrained bonds.
         *
         * Goes through all the constraints and checks if the direction of constraint has not changed
         * by the algorithm (i.e. the constraints algorithm arrived to the solution that is closest
         * to the initial system conformation).
         *
         * \param[in] testData        Test data structure.
         * \param[in] pbc             Periodic boundary data.
         * \param[in] testDescription Human-friendly message with test parameters.
         */
        void checkConstrainsDirection(ConstraintsTestData *testData, t_pbc pbc,
                                      const std::string &testDescription)
        {

            for (unsigned c = 0; c < testData->constraints.size()/3; c++)
            {
                int  i  = testData->constraints.at(3*c + 1);
                int  j  = testData->constraints.at(3*c + 2);
                RVec xij0, xij1;
                if (pbc.ePBC == epbcXYZ)
                {
                    pbc_dx_aiuc(&pbc, testData->x[i], testData->x[j], xij0);
                    pbc_dx_aiuc(&pbc, testData->xPrime[i], testData->xPrime[j], xij1);
                }
                else
                {
                    rvec_sub(testData->x[i], testData->x[j], xij0);
                    rvec_sub(testData->xPrime[i], testData->xPrime[j], xij1);
                }

                real dot = xij0.dot(xij1);

                EXPECT_GE(dot, 0.0) << "The constraint " << c << " changed direction. Constraining algorithm "
                << "could return the wrong root of the constraints equation. " << testDescription;

            }
        }

        /*! \brief
         * The test on the coordinates of the center of the mass (COM) of the system.
         *
         * Checks if the center of mass has not been shifted by the constraints. Note,
         * that this test does not take into account the periodic boundary conditions.
         * Hence it will not work should the constraints decide to move atoms across
         * PBC borders.
         *
         * \param[in] tolerance       Allowed tolerance in COM coordinates.
         * \param[in] testData        Test data structure.
         * \param[in] testDescription Human-friendly message with test parameters.
         */
        void checkCOMCoordinates(FloatingPointTolerance tolerance, ConstraintsTestData *testData,
                                 const std::string &testDescription)
        {

            RVec comPrime0({0.0, 0.0, 0.0});
            RVec comPrime({0.0, 0.0, 0.0});
            for (int i = 0; i < testData->nAtom; i++)
            {
                comPrime0 += testData->masses[i]*testData->xPrime0[i];
                comPrime  += testData->masses[i]*testData->xPrime[i];
            }

            comPrime0 /= testData->nAtom;
            comPrime  /= testData->nAtom;

            EXPECT_REAL_EQ_TOL(comPrime[XX], comPrime0[XX], tolerance)
            << "Center of mass was shifted by constraints in x-derection. " << testDescription;
            EXPECT_REAL_EQ_TOL(comPrime[YY], comPrime0[YY], tolerance)
            << "Center of mass was shifted by constraints in y-derection. " << testDescription;
            EXPECT_REAL_EQ_TOL(comPrime[ZZ], comPrime0[ZZ], tolerance)
            << "Center of mass was shifted by constraints in z-derection. " << testDescription;

        }

        /*! \brief
         * The test on the velocity of the center of the mass (COM) of the system.
         *
         * Checks if the velocity of the center of mass has not changed.
         *
         * \param[in] tolerance       Allowed tolerance in COM velocity components.
         * \param[in] testData        Test data structure.
         * \param[in] testDescription Human-friendly message with test parameters.
         */
        void checkCOMVelocity(FloatingPointTolerance tolerance, ConstraintsTestData *testData,
                              const std::string &testDescription)
        {

            RVec comV0({0.0, 0.0, 0.0});
            RVec comV({0.0, 0.0, 0.0});
            for (int i = 0; i < testData->nAtom; i++)
            {
                comV0 += testData->masses[i]*testData->v0[i];
                comV  += testData->masses[i]*testData->v[i];
            }
            comV0 /= testData->nAtom;
            comV  /= testData->nAtom;

            EXPECT_REAL_EQ_TOL(comV[XX], comV0[XX], tolerance)
            << "Velocity of the center of mass in x-derection has been changed by constraints. "
            << testDescription;
            EXPECT_REAL_EQ_TOL(comV[YY], comV0[YY], tolerance)
            << "Velocity of the center of mass in y-derection has been changed by constraints. "
            << testDescription;
            EXPECT_REAL_EQ_TOL(comV[ZZ], comV0[ZZ], tolerance)
            << "Velocity of the center of mass in z-derection has been changed by constraints. "
            << testDescription;
        }

        /*! \brief
         * The test on the final coordinates (not used).
         *
         * Goes through all atoms and checks if the final positions correspond to the
         * provided reference set of coordinates.
         *
         * \param[in] xPrimeRef       The reference set of coordinates.
         * \param[in] tolerance       Tolerance for the coordinates test.
         * \param[in] testData        Test data structure.
         * \param[in] testDescription Human-friendly message with test parameters.
         */
        void checkFinalCoordinates(std::vector<RVec> xPrimeRef, FloatingPointTolerance tolerance,
                                   ConstraintsTestData *testData, const std::string &testDescription)
        {
            for (int i = 0; i < testData->nAtom; i++)
            {
                for (int d = 0; d < DIM; d++)
                {
                    EXPECT_REAL_EQ_TOL(xPrimeRef.at(i)[d], testData->xPrime[i][d], tolerance) <<
                    "Coordinates after constrains were applied differ from these in the " <<
                    "reference set for atom #" << i << ". " << testDescription;
                }
            }
        }

        /*! \brief
         * The test of virial tensor.
         *
         * Checks if the values in the scaled virial tensor are equal to pre-computed values.
         *
         * \param[in] tolerance       Tolerance for the tensor values.
         * \param[in] testData        Test data structure.
         * \param[in] testDescription Human-friendly message with test parameters.
         */
        void checkVirialTensor(FloatingPointTolerance tolerance, ConstraintsTestData *testData,
                               const std::string &testDescription)
        {
            for (int i = 0; i < DIM; i++)
            {
                for (int j = 0; j < DIM; j++)
                {
                    EXPECT_REAL_EQ_TOL(testData->virialScaledRef[i][j], testData->virialScaled[i][j],
                                       tolerance) <<
                    "Values in virial tensor at [" << i << "][" << j << "] are not within the " <<
                    "tolerance from reference value. " << testDescription;
                }
            }
        }

        /*! \brief
         * The test for FEP (not used).
         *
         * Checks if the value of dH/dLambda is equal to the reference value.
         *
         * \param[in] dHdLambdaRef    Reference value.
         * \param[in] tolerance       Tolerance.
         * \param[in] testData        Test data structure.
         * \param[in] testDescription Human-friendly message with test parameters.
         */
        void checkFEP(const real dHdLambdaRef, FloatingPointTolerance tolerance,
                      ConstraintsTestData *testData, const std::string &testDescription)
        {
            EXPECT_REAL_EQ_TOL(dHdLambdaRef, testData->dHdLambda, tolerance) <<
            "Computed value for dV/dLambda is not equal to the reference value. " << testDescription;
        }



};

//! Index for SHAKE
static constexpr int TEST_ALGORITHM_SHAKE = 1;
//! Index for LINCS
static constexpr int TEST_ALGORITHM_LINCS = 2;

TEST_P(ConstraintsTest, ShakeLincsTest)
{
    int  pbcId, systemId, algorithmId;
    // Make some symbolic names for the parameter combination under
    // test.
    std::tie(pbcId, systemId, algorithmId) = GetParam();

    t_pbc                  pbc       = pbcs.at(pbcId);
    ConstraintsTestData   *testData  = systems.at(systemId);
    FloatingPointTolerance tolerance = absoluteTolerance(0.0001);

    // Make a string that describes which parameter combination is
    // being tested, to help make failing tests comprehensible.
    std::string algorithmName;
    switch (algorithmId)
    {
        case TEST_ALGORITHM_SHAKE:
            algorithmName = "SHAKE";
            break;
        case TEST_ALGORITHM_LINCS:
            algorithmName = "LINCS";
            break;
        default:
            algorithmName = "unknown";
    }
    std::string testDescription = "(Testing " + algorithmName + " on  " + testData->title;
    if (pbc.ePBC == epbcXYZ)
    {
        testDescription += " with PBC";
    }
    testDescription += ").";

    switch (algorithmId)
    {
        case TEST_ALGORITHM_SHAKE:
            tolerance = absoluteTolerance(testData->shakeExpectedTolerance);
            applyShake(testData, testDescription);
            break;
        case TEST_ALGORITHM_LINCS:
            tolerance = absoluteTolerance(testData->lincsExpectedTolerance);
            applyLincs(testData, pbc, testDescription);
            break;
    }


    checkConstrainsLength(tolerance, testData, pbc, testDescription);
    checkConstrainsDirection(testData, pbc, testDescription);
    checkCOMCoordinates(tolerance, testData, testDescription);
    checkCOMVelocity(tolerance, testData, testDescription);
    if (testData->computeVirial)
    {
        checkVirialTensor(tolerance, testData, testDescription);
    }

    testData->reset();
}

INSTANTIATE_TEST_CASE_P(WithParameters, ConstraintsTest,
                            ::testing::Combine(::testing::Values(0, 1),
                                                   ::testing::Values(0, 1, 2, 3, 4, 5),
                                                   ::testing::Values(TEST_ALGORITHM_SHAKE, TEST_ALGORITHM_LINCS)));

} // namespace test
} // namespace gmx
