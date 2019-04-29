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
 * \todo Better tests for virial are needed.
 * \todo Tests for bigger systems to test threads synchronization,
 *       reduction, etc. on the GPU.
 * \todo Tests for algorithms for derivatives.
 * \todo Free-energy perturbation tests
 *
 * \author Artem Zhmurov <zhmurov@gmail.com>
 * \ingroup module_mdlib
 */

#include "gmxpre.h"

#include "gromacs/mdlib/constr.h"

#include "config.h"

#include <assert.h>

#include <cmath>

#include <algorithm>
#include <unordered_map>
#include <vector>

#include <gtest/gtest.h>

#include "gromacs/fileio/gmxfio.h"
#include "gromacs/gmxlib/nrnb.h"
#include "gromacs/gmxlib/nonbonded/nonbonded.h"
#include "gromacs/gpu_utils/gpu_utils.h"
#include "gromacs/math/paddedvector.h"
#include "gromacs/math/vec.h"
#include "gromacs/math/vectypes.h"
#include "gromacs/mdlib/gmx_omp_nthreads.h"
#include "gromacs/mdlib/lincs.h"
#include "gromacs/mdlib/lincs_cuda.h"
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

/*! \brief
 * Constraints test data structure.
 *
 * Structure to collect all the necessary data, including system coordinates and topology,
 * constraints information, etc. The structure can be reset and reused.
 */
struct ConstraintsTestData
{
    public:
        //! Human-friendly name for a system
        std::string           title_;
        //! Number of atoms
        int                   numAtoms_;
        //! Topology
        gmx_mtop_t            mtop_;
        //! Masses
        std::vector<real>     masses_;
        //! Inverse masses
        std::vector<real>     invmass_;
        //! Communication record
        t_commrec             cr_;
        //! Input record (info that usually in .mdp file)
        t_inputrec            ir_;
        //! Local topology
        t_idef                idef_;
        //! MD atoms
        t_mdatoms             md_;
        //! Multisim data
        gmx_multisim_t        ms_;
        //! Computational time array (normally used to benchmark performance)
        t_nrnb                nrnb_;

        //! Inverse timestep
        real                  invdt_;
        //! Number of flexible constraints
        int                   nflexcon_   = 0;
        //! Whether the virial should be computed
        bool                  computeVirial_;
        //! Scaled virial
        tensor                virialScaled_;
        //! Scaled virial (reference values)
        tensor                virialScaledRef_;
        //! If the free energy is computed
        bool                  compute_dHdLambda_;
        //! For free energy computation
        real                  dHdLambda_;
        //! For free energy computation (reference value)
        real                  dHdLambdaRef_;

        //! Coordinates before the timestep
        PaddedVector<RVec>    x_;
        //! Coordinates after timestep, output for the constraints
        PaddedVector<RVec>    xPrime_;
        //! Backup for coordinates (for reset)
        PaddedVector<RVec>    xPrime0_;
        //! Intermediate set of coordinates (normally used for projection correction)
        PaddedVector<RVec>    xPrime2_;
        //! Velocities
        PaddedVector<RVec>    v_;
        //! Backup for velocities (for reset)
        PaddedVector<RVec>    v0_;

        //! Constraints data (type1-i1-j1-type2-i2-j2-...)
        std::vector<int>      constraints_;
        //! Target lengths for all constraint types
        std::vector<real>     constraintsR0_;

        /*! \brief
         * Constructor for the object with all parameters and variables needed by constraints algorithms.
         *
         * This constructor assembles stubs for all the data structures, required to initialize
         * and apply LINCS and SHAKE constraints. The coordinates and velocities before constraining
         * are saved to allow for reset. The constraints data are stored for testing after constraints
         * were applied.
         *
         * \param[in]  title                Human-friendly name of the system.
         * \param[in]  numAtoms             Number of atoms in the system.
         * \param[in]  masses               Atom masses. Size of this vector should be equal to numAtoms.
         * \param[in]  constraints          List of constraints, organized in triples of integers.
         *                                  First integer is the index of type for a constraint, second
         *                                  and third are the indices of constrained atoms. The types
         *                                  of constraints should be sequential but not necessarily
         *                                  start from zero (which is the way they normally are in
         *                                  GROMACS).
         * \param[in]  constraintsR0        Target values for bond lengths for bonds of each type. The
         *                                  size of this vector should be equal to the total number of
         *                                  unique types in constraints vector.
         * \param[in]  computeVirial        Whether the virial should be computed.
         * \param[in]  virialScaledRef      Reference values for scaled virial tensor.
         * \param[in]  compute_dHdLambda    Whether free energy should be computed.
         * \param[in]  dHdLambdaRef         Reference value for dHdLambda.
         * \param[in]  initialTime          Initial time.
         * \param[in]  timestep             Timestep.
         * \param[in]  x                    Coordinates before integration step.
         * \param[in]  xPrime               Coordinates after integration step, but before constraining.
         * \param[in]  v                    Velocities before constraining.
         * \param[in]  shakeTolerance       Target tolerance for SHAKE.
         * \param[in]  shakeUseSOR          Use successive over-relaxation method for SHAKE iterations.
         *                                  The general formula is:
         *                                     x_n+1 = (1-omega)*x_n + omega*f(x_n),
         *                                  where omega = 1 if SOR is off and may be < 1 if SOR is on.
         * \param[in]  lincsNumIterations   Number of iterations used to compute the inverse matrix.
         * \param[in]  lincsExpansionOrder  The order for algorithm that adjusts the direction of the
         *                                  bond after constraints are applied.
         * \param[in]  lincsWarnAngle       The threshold value for the change in bond angle. When
         *                                  exceeded the program will issue a warning.
         *
         */
        ConstraintsTestData(const std::string &title,
                            int numAtoms, std::vector<real> masses,
                            std::vector<int> constraints, std::vector<real> constraintsR0,
                            bool computeVirial, tensor virialScaledRef,
                            bool compute_dHdLambda, float dHdLambdaRef,
                            real initialTime, real timestep,
                            const std::vector<RVec> &x, const std::vector<RVec> &xPrime, const std::vector<RVec> &v,
                            real shakeTolerance, gmx_bool shakeUseSOR,
                            int lincsNumIterations, int lincsExpansionOrder, real lincsWarnAngle)
        {
            title_    = title;    // Human-friendly name of the system
            numAtoms_ = numAtoms; // Number of atoms

            // Masses of atoms
            masses_ = masses;
            invmass_.resize(numAtoms); // Vector of inverse masses

            for (int i = 0; i < numAtoms; i++)
            {
                invmass_[i] = 1.0/masses.at(i);
            }

            // Saving constraints to check if they are satisfied after algorithm was applied
            constraints_   = constraints;   // Constraints indices (in type-i-j format)
            constraintsR0_ = constraintsR0; // Equilibrium distances for each type of constraint

            invdt_  = 1.0/timestep;         // Inverse timestep

            // Communication record
            cr_.nnodes = 1;
            cr_.dd     = nullptr;

            // Multisim data
            ms_.sim  = 0;
            ms_.nsim = 1;

            // Input record - data that usually comes from configuration file (.mdp)
            ir_.efep           = 0;
            ir_.init_t         = initialTime;
            ir_.delta_t        = timestep;
            ir_.eI             = 0;

            // MD atoms data
            md_.nMassPerturbed = 0;
            md_.lambda         = 0.0;
            md_.invmass        = invmass_.data();
            md_.nr             = numAtoms;
            md_.homenr         = numAtoms;

            // Virial evaluation
            computeVirial_ = computeVirial;
            if (computeVirial)
            {
                for (int i = 0; i < DIM; i++)
                {
                    for (int j = 0; j < DIM; j++)
                    {
                        virialScaled_[i][j]    = 0;
                        virialScaledRef_[i][j] = virialScaledRef[i][j];
                    }
                }
            }


            // Free energy evaluation
            compute_dHdLambda_ = compute_dHdLambda;
            dHdLambda_         = 0;
            if (compute_dHdLambda_)
            {
                ir_.efep                    = efepYES;
                dHdLambdaRef_               = dHdLambdaRef;
            }
            else
            {
                ir_.efep                    = efepNO;
                dHdLambdaRef_               = 0;
            }

            // Constraints and their parameters (local topology)
            for (int i = 0; i < F_NRE; i++)
            {
                idef_.il[i].nr = 0;
            }
            idef_.il[F_CONSTR].nr   = constraints.size();

            snew(idef_.il[F_CONSTR].iatoms, constraints.size());
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
                idef_.il[F_CONSTR].iatoms[i] = constraints.at(i);
            }
            snew(idef_.iparams, maxType + 1);
            for (unsigned i = 0; i < constraints.size()/3; i++)
            {
                idef_.iparams[constraints.at(3*i)].constr.dA = constraintsR0.at(constraints.at(3*i));
                idef_.iparams[constraints.at(3*i)].constr.dB = constraintsR0.at(constraints.at(3*i));
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
            molType.atoms.nr             = numAtoms;
            molType.ilist.at(F_CONSTR)   = interactionList;
            molType.ilist.at(F_CONSTRNC) = interactionListEmpty;
            mtop_.moltype.push_back(molType);

            gmx_molblock_t molBlock;
            molBlock.type = 0;
            molBlock.nmol = 1;
            mtop_.molblock.push_back(molBlock);

            mtop_.natoms = numAtoms;
            mtop_.ffparams.iparams.resize(maxType + 1);
            for (int i = 0; i <= maxType; i++)
            {
                mtop_.ffparams.iparams.at(i) = idef_.iparams[i];
            }
            mtop_.bIntermolecularInteractions = false;

            // Coordinates and velocities
            x_.resizeWithPadding(numAtoms);
            xPrime_.resizeWithPadding(numAtoms);
            xPrime0_.resizeWithPadding(numAtoms);
            xPrime2_.resizeWithPadding(numAtoms);

            v_.resizeWithPadding(numAtoms);
            v0_.resizeWithPadding(numAtoms);

            std::copy(x.begin(), x.end(), x_.begin());
            std::copy(xPrime.begin(), xPrime.end(), xPrime_.begin());
            std::copy(xPrime.begin(), xPrime.end(), xPrime0_.begin());
            std::copy(xPrime.begin(), xPrime.end(), xPrime2_.begin());

            std::copy(v.begin(), v.end(), v_.begin());
            std::copy(v.begin(), v.end(), v0_.begin());

            // SHAKE-specific parameters
            ir_.shake_tol            = shakeTolerance;
            ir_.bShakeSOR            = shakeUseSOR;

            // LINCS-specific parameters
            ir_.nLincsIter     = lincsNumIterations;
            ir_.nProjOrder     = lincsExpansionOrder;
            ir_.LincsWarnAngle = lincsWarnAngle;
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
            xPrime_  = xPrime0_;
            xPrime2_ = xPrime0_;
            v_       = v0_;

            if (computeVirial_)
            {
                for (int i = 0; i < DIM; i++)
                {
                    for (int j = 0; j < DIM; j++)
                    {
                        virialScaled_[i][j] = 0;
                    }
                }
            }
            dHdLambda_         = 0;
        }

        /*! \brief
         * Cleaning up the memory.
         */
        ~ConstraintsTestData()
        {
            sfree(idef_.il[F_CONSTR].iatoms);
            sfree(idef_.iparams);
        }

};

/*! \brief The two-dimensional parameter space for test.
 *
 * The test will run for all possible combinations of accessible
 * values of the:
 * 1. PBC setup ("PBCNONE" or "PBCXYZ")
 * 2. The algorithm ("SHAKE", "LINCS" or "LINCS_GPU").
 */
typedef std::tuple<std::string, std::string> ConstraintsTestParameters;

/*! \brief Names of all availible algorithms
 *
 * Constructed from the algorithms_ field of the test class.
 * Used as the list of values of second parameter in ConstraintsTestParameters.
 */
static std::vector<std::string> algorithmsNames;


/*! \brief Test fixture for constraints.
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
        //! PBC setups
        std::unordered_map <std::string, t_pbc>                                             pbcs_;
        //! Algorithms (SHAKE and LINCS)
        std::unordered_map <std::string, void(*)(ConstraintsTestData *testData, t_pbc pbc)> algorithms_;

        /*! \brief Test setup function.
         *
         * Setting up the pbcs and algorithms. Note, that corresponding string keywords
         * have to be explicitly added at the end of this file when the tests are called.
         *
         */
        void SetUp() override
        {

            //
            // PBC initialization
            //
            t_pbc pbc;

            // Infinitely small box
            matrix boxNone = { {0, 0, 0}, {0, 0, 0}, {0, 0, 0} };
            set_pbc(&pbc, epbcNONE, boxNone);
            pbcs_["PBCNone"] = pbc;

            // Rectangular box
            matrix boxXyz = { {10.0, 0.0, 0.0}, {0.0, 20.0, 0.0}, {0.0, 0.0, 15.0} };
            set_pbc(&pbc, epbcXYZ, boxXyz);
            pbcs_["PBCXYZ"] = pbc;

            //
            // Algorithms
            //
            // SHAKE
            algorithms_["SHAKE"] = applyShake;
            // LINCS
            algorithms_["LINCS"] = applyLincs;
            // LINCS using CUDA (if CUDA is available)
#if GMX_GPU == GMX_GPU_CUDA
            std::string errorMessage;
            if (canDetectGpus(&errorMessage))
            {
                algorithms_["LINCS_CUDA"] = applyLincsCuda;
            }
#endif
            for (const auto &algorithm : algorithms_)
            {
                algorithmsNames.push_back(algorithm.first);
            }
        }

        /*! \brief
         * Initialize and apply SHAKE constraints.
         *
         * \param[in] testData        Test data structure.
         * \param[in] pbc             Periodic boundary data (not used in SHAKE).
         */
        static void applyShake(ConstraintsTestData *testData, t_pbc gmx_unused pbc)
        {
            shakedata* shaked = shake_init();
            make_shake_sblock_serial(shaked, &testData->idef_, testData->md_);
            bool       success = constrain_shake(
                        nullptr,
                        shaked,
                        testData->invmass_.data(),
                        testData->idef_,
                        testData->ir_,
                        as_rvec_array(testData->x_.data()),
                        as_rvec_array(testData->xPrime_.data()),
                        as_rvec_array(testData->xPrime2_.data()),
                        &testData->nrnb_,
                        testData->md_.lambda,
                        &testData->dHdLambda_,
                        testData->invdt_,
                        as_rvec_array(testData->v_.data()),
                        testData->computeVirial_,
                        testData->virialScaled_,
                        false,
                        gmx::ConstraintVariable::Positions);
            EXPECT_TRUE(success) << "Test failed with a false return value in SHAKE.";
            done_shake(shaked);
        }

        /*! \brief
         * Initialize and apply LINCS constraints.
         *
         * \param[in] testData        Test data structure.
         * \param[in] pbc             Periodic boundary data.
         */
        static void applyLincs(ConstraintsTestData *testData, t_pbc pbc)
        {

            Lincs                *lincsd;
            int                   maxwarn         = 100;
            int                   warncount_lincs = 0;
            gmx_omp_nthreads_set(emntLINCS, 1);

            // Make blocka structure for faster LINCS setup
            std::vector<t_blocka> at2con_mt;
            at2con_mt.reserve(testData->mtop_.moltype.size());
            for (const gmx_moltype_t &moltype : testData->mtop_.moltype)
            {
                // This function is in constr.cpp
                at2con_mt.push_back(make_at2con(moltype,
                                                testData->mtop_.ffparams.iparams,
                                                flexibleConstraintTreatment(EI_DYNAMICS(testData->ir_.eI))));
            }
            // Initialize LINCS
            lincsd = init_lincs(nullptr, testData->mtop_,
                                testData->nflexcon_, at2con_mt,
                                false,
                                testData->ir_.nLincsIter, testData->ir_.nProjOrder);
            set_lincs(testData->idef_, testData->md_, EI_DYNAMICS(testData->ir_.eI), &testData->cr_, lincsd);

            // Evaluate constraints
            bool sucess = constrain_lincs(false, testData->ir_, 0, lincsd, testData->md_,
                                          &testData->cr_,
                                          &testData->ms_,
                                          as_rvec_array(testData->x_.data()),
                                          as_rvec_array(testData->xPrime_.data()),
                                          as_rvec_array(testData->xPrime2_.data()),
                                          pbc.box, &pbc,
                                          testData->md_.lambda, &testData->dHdLambda_,
                                          testData->invdt_,
                                          as_rvec_array(testData->v_.data()),
                                          testData->computeVirial_, testData->virialScaled_,
                                          gmx::ConstraintVariable::Positions,
                                          &testData->nrnb_,
                                          maxwarn, &warncount_lincs);
            EXPECT_TRUE(sucess) << "Test failed with a false return value in LINCS.";
            EXPECT_EQ(warncount_lincs, 0) << "There were warnings in LINCS.";
            for (unsigned int i = 0; i < testData->mtop_.moltype.size(); i++)
            {
                sfree(at2con_mt.at(i).index);
                sfree(at2con_mt.at(i).a);
            }
            done_lincs(lincsd);
        }

        /*! \brief
         * Initialize and apply LINCS constraints on CUDA-enabled GPU.
         *
         * \param[in] testData        Test data structure.
         * \param[in] pbc             Periodic boundary data.
         */
        static void applyLincsCuda(ConstraintsTestData *testData, t_pbc pbc)
        {
            auto lincsCuda = std::make_unique<LincsCuda>(testData->numAtoms_,
                                                         testData->ir_.nLincsIter,
                                                         testData->ir_.nProjOrder);
            lincsCuda->set(testData->idef_, testData->md_);
            lincsCuda->setPbc(&pbc);
            lincsCuda->copyCoordinatesToGpu(as_rvec_array(testData->x_.data()),
                                            as_rvec_array(testData->xPrime_.data()));
            lincsCuda->copyVelocitiesToGpu(as_rvec_array(testData->v_.data()));
            lincsCuda->apply(true, testData->invdt_, testData->computeVirial_, testData->virialScaled_);
            lincsCuda->copyCoordinatesFromGpu(as_rvec_array(testData->xPrime_.data()));
            lincsCuda->copyVelocitiesFromGpu(as_rvec_array(testData->v_.data()));
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
         */
        void checkConstrainsLength(FloatingPointTolerance tolerance, const ConstraintsTestData &testData, t_pbc pbc)
        {

            // Test if all the constraints are satisfied
            for (unsigned c = 0; c < testData.constraints_.size()/3; c++)
            {
                real r0 = testData.constraintsR0_.at(testData.constraints_.at(3*c));
                int  i  = testData.constraints_.at(3*c + 1);
                int  j  = testData.constraints_.at(3*c + 2);
                RVec xij0, xij1;
                real d0, d1;
                if (pbc.ePBC == epbcXYZ)
                {
                    pbc_dx_aiuc(&pbc, testData.x_[i], testData.x_[j], xij0);
                    pbc_dx_aiuc(&pbc, testData.xPrime_[i], testData.xPrime_[j], xij1);
                }
                else
                {
                    rvec_sub(testData.x_[i], testData.x_[j], xij0);
                    rvec_sub(testData.xPrime_[i], testData.xPrime_[j], xij1);
                }
                d0 = norm(xij0);
                d1 = norm(xij1);
                EXPECT_REAL_EQ_TOL(r0, d1, tolerance) << gmx::formatString(
                        "rij = %f, which is not equal to r0 = %f for constraint #%u, between atoms %d and %d"
                        " (before constraining rij was %f).", d1, r0, c, i, j, d0);
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
         */
        void checkConstrainsDirection(const ConstraintsTestData &testData, t_pbc pbc)
        {

            for (unsigned c = 0; c < testData.constraints_.size()/3; c++)
            {
                int  i  = testData.constraints_.at(3*c + 1);
                int  j  = testData.constraints_.at(3*c + 2);
                RVec xij0, xij1;
                if (pbc.ePBC == epbcXYZ)
                {
                    pbc_dx_aiuc(&pbc, testData.x_[i], testData.x_[j], xij0);
                    pbc_dx_aiuc(&pbc, testData.xPrime_[i], testData.xPrime_[j], xij1);
                }
                else
                {
                    rvec_sub(testData.x_[i], testData.x_[j], xij0);
                    rvec_sub(testData.xPrime_[i], testData.xPrime_[j], xij1);
                }

                real dot = xij0.dot(xij1);

                EXPECT_GE(dot, 0.0) << gmx::formatString(
                        "The constraint %u changed direction. Constraining algorithm might have returned the wrong root "
                        "of the constraints equation.", c);

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
         */
        void checkCOMCoordinates(FloatingPointTolerance tolerance, const ConstraintsTestData &testData)
        {

            RVec comPrime0({0.0, 0.0, 0.0});
            RVec comPrime({0.0, 0.0, 0.0});
            for (int i = 0; i < testData.numAtoms_; i++)
            {
                comPrime0 += testData.masses_[i]*testData.xPrime0_[i];
                comPrime  += testData.masses_[i]*testData.xPrime_[i];
            }

            comPrime0 /= testData.numAtoms_;
            comPrime  /= testData.numAtoms_;

            EXPECT_REAL_EQ_TOL(comPrime[XX], comPrime0[XX], tolerance)
            << "Center of mass was shifted by constraints in x-direction.";
            EXPECT_REAL_EQ_TOL(comPrime[YY], comPrime0[YY], tolerance)
            << "Center of mass was shifted by constraints in y-direction.";
            EXPECT_REAL_EQ_TOL(comPrime[ZZ], comPrime0[ZZ], tolerance)
            << "Center of mass was shifted by constraints in z-direction.";

        }

        /*! \brief
         * The test on the velocity of the center of the mass (COM) of the system.
         *
         * Checks if the velocity of the center of mass has not changed.
         *
         * \param[in] tolerance       Allowed tolerance in COM velocity components.
         * \param[in] testData        Test data structure.
         */
        void checkCOMVelocity(FloatingPointTolerance tolerance, const ConstraintsTestData &testData)
        {

            RVec comV0({0.0, 0.0, 0.0});
            RVec comV({0.0, 0.0, 0.0});
            for (int i = 0; i < testData.numAtoms_; i++)
            {
                comV0 += testData.masses_[i]*testData.v0_[i];
                comV  += testData.masses_[i]*testData.v_[i];
            }
            comV0 /= testData.numAtoms_;
            comV  /= testData.numAtoms_;

            EXPECT_REAL_EQ_TOL(comV[XX], comV0[XX], tolerance)
            << "Velocity of the center of mass in x-direction has been changed by constraints.";
            EXPECT_REAL_EQ_TOL(comV[YY], comV0[YY], tolerance)
            << "Velocity of the center of mass in y-direction has been changed by constraints.";
            EXPECT_REAL_EQ_TOL(comV[ZZ], comV0[ZZ], tolerance)
            << "Velocity of the center of mass in z-direction has been changed by constraints.";
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
         */
        void checkFinalCoordinates(std::vector<RVec> xPrimeRef, FloatingPointTolerance tolerance,
                                   const ConstraintsTestData &testData)
        {
            for (int i = 0; i < testData.numAtoms_; i++)
            {
                for (int d = 0; d < DIM; d++)
                {
                    EXPECT_REAL_EQ_TOL(xPrimeRef.at(i)[d], testData.xPrime_[i][d], tolerance) << gmx::formatString(
                            "Coordinates after constrains were applied differ from these in the reference set for atom #%d.", i);
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
         */
        void checkVirialTensor(FloatingPointTolerance tolerance, const ConstraintsTestData &testData)
        {
            for (int i = 0; i < DIM; i++)
            {
                for (int j = 0; j < DIM; j++)
                {
                    EXPECT_REAL_EQ_TOL(testData.virialScaledRef_[i][j], testData.virialScaled_[i][j],
                                       tolerance) << gmx::formatString(
                            "Values in virial tensor at [%d][%d] are not within the tolerance from reference value.", i, j);
                }
            }
        }

        /*! \brief
         * The test for FEP (not used).
         *
         * Checks if the value of dH/dLambda is equal to the reference value.
         * \todo Add tests for dHdLambda values.
         *
         * \param[in] dHdLambdaRef    Reference value.
         * \param[in] tolerance       Tolerance.
         * \param[in] testData        Test data structure.
         */
        void checkFEP(const real dHdLambdaRef, FloatingPointTolerance tolerance,
                      const ConstraintsTestData &testData)
        {
            EXPECT_REAL_EQ_TOL(dHdLambdaRef, testData.dHdLambda_, tolerance) <<
            "Computed value for dV/dLambda is not equal to the reference value. ";
        }

};

TEST_P(ConstraintsTest, SingleConstraint){
    std::string       title    = "one constraint (e.g. OH)";
    int               numAtoms = 2;

    std::vector<real> masses        = {1.0, 12.0};
    std::vector<int>  constraints   = {0, 0, 1};
    std::vector<real> constraintsR0 = {0.1};

    real              oneTenthOverSqrtTwo    = 0.1_real / std::sqrt(2.0_real);

    std::vector<RVec> x = {{ 0.0, oneTenthOverSqrtTwo, 0.0 },
                           { oneTenthOverSqrtTwo, 0.0, 0.0 }};
    std::vector<RVec> xPrime = {{ 0.01,  0.08,  0.01 },
                                { 0.06,  0.01, -0.01 }};
    std::vector<RVec> v = {{ 1.0, 2.0, 3.0 },
                           { 3.0, 2.0, 1.0 }};

    tensor            virialScaledRef = {{-5.58e-04,  5.58e-04, 0.00e+00 },
                                         { 5.58e-04, -5.58e-04, 0.00e+00 },
                                         { 0.00e+00,  0.00e+00, 0.00e+00 }};

    real              shakeTolerance         = 0.0001;
    gmx_bool          shakeUseSOR            = false;

    int               lincsNIter                      = 1;
    int               lincslincsExpansionOrder        = 4;
    real              lincsWarnAngle                  = 30.0;

    std::unique_ptr<ConstraintsTestData>   testData = std::make_unique<ConstraintsTestData>
            (title, numAtoms, masses,
            constraints, constraintsR0,
            true, virialScaledRef,
            false, 0,
            real(0.0), real(0.001),
            x, xPrime, v,
            shakeTolerance, shakeUseSOR,
            lincsNIter, lincslincsExpansionOrder, lincsWarnAngle);
    std::string pbcName;
    std::string algorithmName;
    std::tie(pbcName, algorithmName) = GetParam();
    t_pbc                                   pbc       = pbcs_.at(pbcName);

    // Apply constraints
    algorithms_.at(algorithmName)(testData.get(), pbc);

    checkConstrainsLength(absoluteTolerance(0.0002), *testData, pbc);
    checkConstrainsDirection(*testData, pbc);
    checkCOMCoordinates(absoluteTolerance(0.0001), *testData);
    checkCOMVelocity(absoluteTolerance(0.0001), *testData);

    checkVirialTensor(absoluteTolerance(0.0001), *testData);

}

TEST_P(ConstraintsTest, TwoDisjointConstraints){

    std::string       title            = "two disjoint constraints";
    int               numAtoms         = 4;
    std::vector<real> masses           = {0.5, 1.0/3.0, 0.25, 1.0};
    std::vector<int>  constraints      = {0, 0, 1, 1, 2, 3};
    std::vector<real> constraintsR0    = {2.0, 1.0};


    std::vector<RVec> x = {{  2.50, -3.10, 15.70 },
                           {  0.51, -3.02, 15.55 },
                           { -0.50, -3.00, 15.20 },
                           { -1.51, -2.95, 15.05 }};

    std::vector<RVec> xPrime = {{  2.50, -3.10, 15.70 },
                                {  0.51, -3.02, 15.55 },
                                { -0.50, -3.00, 15.20 },
                                { -1.51, -2.95, 15.05 }};

    std::vector<RVec> v = {{ 0.0, 1.0, 0.0 },
                           { 1.0, 0.0, 0.0 },
                           { 0.0, 0.0, 1.0 },
                           { 0.0, 0.0, 0.0 }};

    tensor            virialScaledRef = {{ 3.3e-03, -1.7e-04,  5.6e-04 },
                                         {-1.7e-04,  8.9e-06, -2.8e-05 },
                                         { 5.6e-04, -2.8e-05,  8.9e-05 }};

    real              shakeTolerance         = 0.0001;
    gmx_bool          shakeUseSOR            = false;

    int               lincsNIter                      = 1;
    int               lincslincsExpansionOrder        = 4;
    real              lincsWarnAngle                  = 30.0;

    std::unique_ptr<ConstraintsTestData>   testData = std::make_unique<ConstraintsTestData>
            (title, numAtoms, masses,
            constraints, constraintsR0,
            true, virialScaledRef,
            false, 0,
            real(0.0), real(0.001),
            x, xPrime, v,
            shakeTolerance, shakeUseSOR,
            lincsNIter, lincslincsExpansionOrder, lincsWarnAngle);

    std::string pbcName;
    std::string algorithmName;
    std::tie(pbcName, algorithmName) = GetParam();
    t_pbc                                   pbc       = pbcs_.at(pbcName);

    // Apply constraints
    algorithms_.at(algorithmName)(testData.get(), pbc);

    checkConstrainsLength(absoluteTolerance(0.0002), *testData, pbc);
    checkConstrainsDirection(*testData, pbc);
    checkCOMCoordinates(absoluteTolerance(0.0001), *testData);
    checkCOMVelocity(absoluteTolerance(0.0001), *testData);

    checkVirialTensor(absoluteTolerance(0.0001), *testData);

}

TEST_P(ConstraintsTest, ThreeSequentialConstraints){

    std::string       title            = "three atoms, connected longitudinally (e.g. CH2)";
    int               numAtoms         = 3;
    std::vector<real> masses           = {1.0, 12.0, 16.0 };
    std::vector<int>  constraints      = {0, 0, 1, 1, 1, 2};
    std::vector<real> constraintsR0    = {0.1, 0.2};

    real              oneTenthOverSqrtTwo    = 0.1_real / std::sqrt(2.0_real);
    real              twoTenthsOverSqrtThree = 0.2_real / std::sqrt(3.0_real);

    std::vector<RVec> x = {{ oneTenthOverSqrtTwo, oneTenthOverSqrtTwo, 0.0 },
                           { 0.0, 0.0, 0.0 },
                           { twoTenthsOverSqrtThree, twoTenthsOverSqrtThree, twoTenthsOverSqrtThree }};

    std::vector<RVec> xPrime = {{  0.08, 0.07,  0.01 },
                                { -0.02, 0.01, -0.02 },
                                {  0.10, 0.12,  0.11 }};

    std::vector<RVec> v = {{ 1.0, 0.0, 0.0 },
                           { 0.0, 1.0, 0.0 },
                           { 0.0, 0.0, 1.0 }};

    tensor            virialScaledRef = {{ 4.14e-03, 4.14e-03, 3.31e-03},
                                         { 4.14e-03, 4.14e-03, 3.31e-03},
                                         { 3.31e-03, 3.31e-03, 3.31e-03}};

    real              shakeTolerance         = 0.0001;
    gmx_bool          shakeUseSOR            = false;

    int               lincsNIter                      = 1;
    int               lincslincsExpansionOrder        = 4;
    real              lincsWarnAngle                  = 30.0;

    std::unique_ptr<ConstraintsTestData>   testData = std::make_unique<ConstraintsTestData>
            (title, numAtoms, masses,
            constraints, constraintsR0,
            true, virialScaledRef,
            false, 0,
            real(0.0), real(0.001),
            x, xPrime, v,
            shakeTolerance, shakeUseSOR,
            lincsNIter, lincslincsExpansionOrder, lincsWarnAngle);

    std::string pbcName;
    std::string algorithmName;
    std::tie(pbcName, algorithmName) = GetParam();
    t_pbc                                   pbc       = pbcs_.at(pbcName);

    // Apply constraints
    algorithms_.at(algorithmName)(testData.get(), pbc);

    checkConstrainsLength(absoluteTolerance(0.0002), *testData, pbc);
    checkConstrainsDirection(*testData, pbc);
    checkCOMCoordinates(absoluteTolerance(0.0001), *testData);
    checkCOMVelocity(absoluteTolerance(0.0001), *testData);

    checkVirialTensor(absoluteTolerance(0.0001), *testData);

}

TEST_P(ConstraintsTest, ThreeConstraintsWithCentralAtom){

    std::string       title            = "three atoms, connected to the central atom (e.g. CH3)";
    int               numAtoms         = 4;
    std::vector<real> masses           = {12.0, 1.0, 1.0, 1.0 };
    std::vector<int>  constraints      = {0, 0, 1, 0, 0, 2, 0, 0, 3};
    std::vector<real> constraintsR0    = {0.1};


    std::vector<RVec> x = {{ 0.00,  0.00,  0.00 },
                           { 0.10,  0.00,  0.00 },
                           { 0.00, -0.10,  0.00 },
                           { 0.00,  0.00,  0.10 }};

    std::vector<RVec> xPrime = {{ 0.004,  0.009, -0.010 },
                                { 0.110, -0.006,  0.003 },
                                {-0.007, -0.102, -0.007 },
                                {-0.005,  0.011,  0.102 }};

    std::vector<RVec> v = {{ 1.0, 0.0, 0.0 },
                           { 1.0, 0.0, 0.0 },
                           { 1.0, 0.0, 0.0 },
                           { 1.0, 0.0, 0.0 }};

    tensor            virialScaledRef = {{7.14e-04, 0.00e+00, 0.00e+00},
                                         {0.00e+00, 1.08e-03, 0.00e+00},
                                         {0.00e+00, 0.00e+00, 1.15e-03}};

    real              shakeTolerance         = 0.0001;
    gmx_bool          shakeUseSOR            = false;

    int               lincsNIter                      = 1;
    int               lincslincsExpansionOrder        = 4;
    real              lincsWarnAngle                  = 30.0;

    std::unique_ptr<ConstraintsTestData>   testData = std::make_unique<ConstraintsTestData>
            (title, numAtoms, masses,
            constraints, constraintsR0,
            true, virialScaledRef,
            false, 0,
            real(0.0), real(0.001),
            x, xPrime, v,
            shakeTolerance, shakeUseSOR,
            lincsNIter, lincslincsExpansionOrder, lincsWarnAngle);

    std::string pbcName;
    std::string algorithmName;
    std::tie(pbcName, algorithmName) = GetParam();
    t_pbc                                   pbc       = pbcs_.at(pbcName);

    // Apply constraints
    algorithms_.at(algorithmName)(testData.get(), pbc);

    checkConstrainsLength(absoluteTolerance(0.0002), *testData, pbc);
    checkConstrainsDirection(*testData, pbc);
    checkCOMCoordinates(absoluteTolerance(0.0001), *testData);
    checkCOMVelocity(absoluteTolerance(0.0001), *testData);

    checkVirialTensor(absoluteTolerance(0.0001), *testData);
}

TEST_P(ConstraintsTest, FourSequentialConstraints){

    std::string       title            = "four atoms, connected longitudinally";
    int               numAtoms         = 4;
    std::vector<real> masses           = {0.5, 1.0/3.0, 0.25, 1.0};
    std::vector<int>  constraints      = {0, 0, 1, 1, 1, 2, 2, 2, 3};
    std::vector<real> constraintsR0    = {2.0, 1.0, 1.0};


    std::vector<RVec> x = {{  2.50, -3.10, 15.70 },
                           {  0.51, -3.02, 15.55 },
                           { -0.50, -3.00, 15.20 },
                           { -1.51, -2.95, 15.05 }};

    std::vector<RVec> xPrime = {{  2.50, -3.10, 15.70 },
                                {  0.51, -3.02, 15.55 },
                                { -0.50, -3.00, 15.20 },
                                { -1.51, -2.95, 15.05 }};

    std::vector<RVec> v = {{ 0.0, 0.0,  2.0 },
                           { 0.0, 0.0,  3.0 },
                           { 0.0, 0.0, -4.0 },
                           { 0.0, 0.0, -1.0 }};

    tensor            virialScaledRef = {{ 1.15e-01, -4.20e-03,  2.12e-02},
                                         {-4.20e-03,  1.70e-04, -6.41e-04},
                                         { 2.12e-02, -6.41e-04,  5.45e-03}};

    real              shakeTolerance         = 0.0001;
    gmx_bool          shakeUseSOR            = false;

    int               lincsNIter                      = 4;
    int               lincslincsExpansionOrder        = 8;
    real              lincsWarnAngle                  = 30.0;

    std::unique_ptr<ConstraintsTestData>   testData = std::make_unique<ConstraintsTestData>
            (title, numAtoms, masses,
            constraints, constraintsR0,
            true, virialScaledRef,
            false, 0,
            real(0.0), real(0.001),
            x, xPrime, v,
            shakeTolerance, shakeUseSOR,
            lincsNIter, lincslincsExpansionOrder, lincsWarnAngle);

    std::string pbcName;
    std::string algorithmName;
    std::tie(pbcName, algorithmName) = GetParam();
    t_pbc                                   pbc       = pbcs_.at(pbcName);

    // Apply constraints
    algorithms_.at(algorithmName)(testData.get(), pbc);

    checkConstrainsLength(absoluteTolerance(0.0002), *testData, pbc);
    checkConstrainsDirection(*testData, pbc);
    checkCOMCoordinates(absoluteTolerance(0.0001), *testData);
    checkCOMVelocity(absoluteTolerance(0.0001), *testData);

    checkVirialTensor(absoluteTolerance(0.01), *testData);

}

TEST_P(ConstraintsTest, TriangleOfConstraints){

    std::string       title            = "basic triangle (tree atoms, connected to each other)";
    int               numAtoms         = 3;
    std::vector<real> masses           = {1.0, 1.0, 1.0};
    std::vector<int>  constraints      = {0, 0, 1, 2, 0, 2, 1, 1, 2};
    std::vector<real> constraintsR0    = {0.1, 0.1, 0.1};

    real              oneTenthOverSqrtTwo    = 0.1_real / std::sqrt(2.0_real);

    std::vector<RVec> x = {{ oneTenthOverSqrtTwo, 0.0, 0.0 },
                           { 0.0, oneTenthOverSqrtTwo, 0.0 },
                           { 0.0, 0.0, oneTenthOverSqrtTwo }};

    std::vector<RVec> xPrime = {{  0.09, -0.02,  0.01 },
                                { -0.02,  0.10, -0.02 },
                                {  0.03, -0.01,  0.07 }};

    std::vector<RVec> v = {{  1.0,  1.0,  1.0 },
                           { -2.0, -2.0, -2.0 },
                           {  1.0,  1.0,  1.0 }};

    tensor            virialScaledRef = {{ 6.00e-04, -1.61e-03,  1.01e-03},
                                         {-1.61e-03,  2.53e-03, -9.25e-04},
                                         { 1.01e-03, -9.25e-04, -8.05e-05}};

    real              shakeTolerance         = 0.0001;
    gmx_bool          shakeUseSOR            = false;

    int               lincsNIter                      = 1;
    int               lincslincsExpansionOrder        = 4;
    real              lincsWarnAngle                  = 30.0;

    std::unique_ptr<ConstraintsTestData>   testData = std::make_unique<ConstraintsTestData>
            (title, numAtoms, masses,
            constraints, constraintsR0,
            true, virialScaledRef,
            false, 0,
            real(0.0), real(0.001),
            x, xPrime, v,
            shakeTolerance, shakeUseSOR,
            lincsNIter, lincslincsExpansionOrder, lincsWarnAngle);

    std::string pbcName;
    std::string algorithmName;
    std::tie(pbcName, algorithmName) = GetParam();
    t_pbc                                   pbc       = pbcs_.at(pbcName);

    // Apply constraints
    algorithms_.at(algorithmName)(testData.get(), pbc);

    checkConstrainsLength(absoluteTolerance(0.0002), *testData, pbc);
    checkConstrainsDirection(*testData, pbc);
    checkCOMCoordinates(absoluteTolerance(0.0001), *testData);
    checkCOMVelocity(absoluteTolerance(0.0001), *testData);

    checkVirialTensor(absoluteTolerance(0.00001), *testData);

}


#if GMX_GPU == GMX_GPU_CUDA
INSTANTIATE_TEST_CASE_P(WithParameters, ConstraintsTest,
                            ::testing::Combine(::testing::Values("PBCNone", "PBCXYZ"),
                                                   ::testing::ValuesIn(algorithmsNames)));
#endif

#if GMX_GPU != GMX_GPU_CUDA
INSTANTIATE_TEST_CASE_P(WithParameters, ConstraintsTest,
                            ::testing::Combine(::testing::Values("PBCNone", "PBCXYZ"),
                                                   ::testing::Values("SHAKE", "LINCS")));
#endif

} // namespace test
} // namespace gmx
