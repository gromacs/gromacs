/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2018- The GROMACS Authors
 * and the project initiators Erik Lindahl, Berk Hess and David van der Spoel.
 * Consult the AUTHORS/COPYING files and https://www.gromacs.org for details.
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
 * https://www.gnu.org/licenses, or write to the Free Software Foundation,
 * Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA.
 *
 * If you want to redistribute modifications to GROMACS, please
 * consider that scientific software is very special. Version
 * control is crucial - bugs must be traceable. We will be happy to
 * consider code for inclusion in the official distribution, but
 * derived work must not be called official GROMACS. Details are found
 * in the README & COPYING files - if they are missing, get the
 * official version at https://www.gromacs.org.
 *
 * To help us fund GROMACS development, we humbly ask that you cite
 * the research papers on the package. Check out https://www.gromacs.org.
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

#include <cassert>
#include <cmath>

#include <memory>
#include <ostream>
#include <string>
#include <tuple>
#include <unordered_map>
#include <vector>

#include <gtest/gtest.h>

#include "gromacs/gpu_utils/capabilities.h"
#include "gromacs/math/paddedvector.h"
#include "gromacs/mdlib/gmx_omp_nthreads.h"
#include "gromacs/mdlib/lincs.h"
#include "gromacs/mdlib/shake.h"
#include "gromacs/mdrunutility/multisim.h"
#include "gromacs/mdtypes/inputrec.h"
#include "gromacs/mdtypes/md_enums.h"
#include "gromacs/pbcutil/pbc.h"
#include "gromacs/topology/forcefieldparameters.h"
#include "gromacs/utility/basedefinitions.h"
#include "gromacs/utility/enumerationhelpers.h"
#include "gromacs/utility/listoflists.h"
#include "gromacs/utility/real.h"
#include "gromacs/utility/stringutil.h"
#include "gromacs/utility/vec.h"
#include "gromacs/utility/vectypes.h"

#include "testutils/hardware_test_fixture.h"
#include "testutils/naming.h"
#include "testutils/refdata.h"
#include "testutils/test_device.h"
#include "testutils/test_hardware_environment.h"
#include "testutils/testasserts.h"

#include "constrtestdata.h"

#if GMX_GPU && !GMX_GPU_OPENCL
#    include "gromacs/gpu_utils/device_context.h"
#    include "gromacs/gpu_utils/device_stream.h"
#    include "gromacs/gpu_utils/devicebuffer.h"
#    include "gromacs/gpu_utils/gputraits.h"
#    include "gromacs/mdlib/lincs_gpu.h"
#endif

//! Helper function to convert t_pbc into string and make test failure messages readable
static void PrintTo(const t_pbc& pbc, std::ostream* os)
{
    *os << "PBC: " << c_pbcTypeNames[pbc.pbcType];
}


namespace gmx
{

namespace test
{
namespace
{

void applyShakeCpu(ConstraintsTestData* testData, t_pbc /* pbc */)
{
    shakedata shaked;
    make_shake_sblock_serial(&shaked, testData->idef_.get(), testData->numAtoms_);
    bool success = constrain_shake(nullptr,
                                   &shaked,
                                   testData->invmass_,
                                   *testData->idef_,
                                   testData->ir_,
                                   testData->x_,
                                   testData->xPrime_,
                                   testData->xPrime2_,
                                   nullptr,
                                   &testData->nrnb_,
                                   testData->lambda_,
                                   &testData->dHdLambda_,
                                   testData->invdt_,
                                   testData->v_,
                                   testData->computeVirial_,
                                   testData->virialScaled_,
                                   false,
                                   gmx::ConstraintVariable::Positions);
    EXPECT_TRUE(success) << "Test failed with a false return value in SHAKE.";
}

void applyLincsCpu(ConstraintsTestData* testData, t_pbc pbc)
{
    Lincs* lincsd;
    int    maxwarn         = 100;
    int    warncount_lincs = 0;
    gmx_omp_nthreads_set(ModuleMultiThread::Lincs, 1);

    gmx_domdec_t* dd = nullptr;

    gmx_multisim_t ms{ 1, 0, MPI_COMM_NULL, MPI_COMM_NULL };

    std::vector<ListOfLists<int>> at2con_mt;
    at2con_mt.reserve(testData->mtop_.moltype.size());
    for (const gmx_moltype_t& moltype : testData->mtop_.moltype)
    {
        at2con_mt.push_back(make_at2con(moltype,
                                        testData->mtop_.ffparams.iparams,
                                        flexibleConstraintTreatment(EI_DYNAMICS(testData->ir_.eI))));
    }
    lincsd = init_lincs(nullptr,
                        testData->mtop_,
                        testData->nflexcon_,
                        at2con_mt,
                        false,
                        testData->ir_.nLincsIter,
                        testData->ir_.nProjOrder,
                        nullptr);
    set_lincs(*testData->idef_,
              testData->numAtoms_,
              testData->invmass_,
              testData->lambda_,
              EI_DYNAMICS(testData->ir_.eI),
              dd,
              lincsd);

    bool success = constrain_lincs(false,
                                   testData->ir_,
                                   0,
                                   lincsd,
                                   testData->invmass_,
                                   dd,
                                   &ms,
                                   testData->x_.arrayRefWithPadding(),
                                   testData->xPrime_.arrayRefWithPadding(),
                                   testData->xPrime2_.arrayRefWithPadding().unpaddedArrayRef(),
                                   pbc.box,
                                   &pbc,
                                   testData->hasMassPerturbed_,
                                   testData->lambda_,
                                   &testData->dHdLambda_,
                                   testData->invdt_,
                                   testData->v_.arrayRefWithPadding().unpaddedArrayRef(),
                                   testData->computeVirial_,
                                   testData->virialScaled_,
                                   gmx::ConstraintVariable::Positions,
                                   &testData->nrnb_,
                                   maxwarn,
                                   &warncount_lincs,
                                   nullptr);
    EXPECT_TRUE(success) << "Test failed with a false return value in LINCS.";
    EXPECT_EQ(warncount_lincs, 0) << "There were warnings in LINCS.";
    done_lincs(lincsd);
}

#if GMX_GPU && !GMX_GPU_OPENCL

void applyLincsGpu(const DeviceContext& deviceContext,
                   const DeviceStream&  deviceStream,
                   ConstraintsTestData* testData,
                   t_pbc                pbc)
{
    deviceContext.activate();

    auto lincsGpu = std::make_unique<LincsGpu>(
            testData->ir_.nLincsIter, testData->ir_.nProjOrder, deviceContext, deviceStream);

    bool updateVelocities = true;
    int  numAtoms         = testData->numAtoms_;

    Float3* h_x  = gmx::asGenericFloat3Pointer(testData->x_);
    Float3* h_xp = gmx::asGenericFloat3Pointer(testData->xPrime_);
    Float3* h_v  = gmx::asGenericFloat3Pointer(testData->v_);

    DeviceBuffer<Float3> d_x, d_xp, d_v;

    lincsGpu->set(*testData->idef_, testData->numAtoms_, testData->invmass_);
    PbcAiuc pbcAiuc;
    setPbcAiuc(pbc.ndim_ePBC, pbc.box, &pbcAiuc);

    allocateDeviceBuffer(&d_x, numAtoms, deviceContext);
    allocateDeviceBuffer(&d_xp, numAtoms, deviceContext);
    allocateDeviceBuffer(&d_v, numAtoms, deviceContext);

    copyToDeviceBuffer(&d_x, h_x, 0, numAtoms, deviceStream, GpuApiCallBehavior::Sync, nullptr);
    copyToDeviceBuffer(&d_xp, h_xp, 0, numAtoms, deviceStream, GpuApiCallBehavior::Sync, nullptr);
    if (updateVelocities)
    {
        copyToDeviceBuffer(&d_v, h_v, 0, numAtoms, deviceStream, GpuApiCallBehavior::Sync, nullptr);
    }
    lincsGpu->apply(
            d_x, d_xp, updateVelocities, d_v, testData->invdt_, testData->computeVirial_, testData->virialScaled_, pbcAiuc);

    copyFromDeviceBuffer(h_xp, &d_xp, 0, numAtoms, deviceStream, GpuApiCallBehavior::Sync, nullptr);
    if (updateVelocities)
    {
        copyFromDeviceBuffer(h_v, &d_v, 0, numAtoms, deviceStream, GpuApiCallBehavior::Sync, nullptr);
    }

    freeDeviceBuffer(&d_x);
    freeDeviceBuffer(&d_xp);
    freeDeviceBuffer(&d_v);
}

#endif // GMX_GPU && !GMX_GPU_OPENCL

// Forward declaration for the test system structure
struct ConstraintsTestSystem;

// Define the set of PBCs to run the test for
const std::vector<t_pbc> c_pbcs = []
{
    std::vector<t_pbc> pbcs;
    t_pbc              pbc;

    // Infinitely small box
    matrix boxNone = { { 0, 0, 0 }, { 0, 0, 0 }, { 0, 0, 0 } };
    set_pbc(&pbc, PbcType::No, boxNone);
    pbcs.emplace_back(pbc);

    // Rectangular box
    matrix boxXyz = { { 10.0, 0.0, 0.0 }, { 0.0, 20.0, 0.0 }, { 0.0, 0.0, 15.0 } };
    set_pbc(&pbc, PbcType::Xyz, boxXyz);
    pbcs.emplace_back(pbc);

    return pbcs;
}();

struct ConstraintsTestSystem
{
    //! Human-friendly name of the system.
    std::string title;
    //! Number of atoms in the system.
    int numAtoms;
    //! Atom masses. Size of this vector should be equal to numAtoms.
    std::vector<real> masses;
    /*! \brief List of constraints, organized in triples of integers.
     *
     * First integer is the index of type for a constraint, second
     * and third are the indices of constrained atoms. The types
     * of constraints should be sequential but not necessarily
     * start from zero (which is the way they normally are in
     * GROMACS).
     */
    std::vector<int> constraints;
    /*! \brief Target values for bond lengths for bonds of each type.
     *
     * The size of this vector should be equal to the total number of
     * unique types in constraints vector.
     */
    std::vector<real> constraintsR0;
    //! Whether the constraint system contains a set of constraints in a triangle
    bool hasTriangle = false;
    //! Coordinates before integration step.
    std::vector<RVec> x;
    //! Coordinates after integration step, but before constraining.
    std::vector<RVec> xPrime;
    //! Velocities before constraining.
    std::vector<RVec> v;

    //! Target tolerance for SHAKE.
    real shakeTolerance = 0.00002;
    /*! \brief Use successive over-relaxation method for SHAKE iterations.
     *
     * The general formula is:
     * x_n+1 = (1-omega)*x_n + omega*f(x_n),
     * where omega = 1 if SOR is off and may be < 1 if SOR is on.
     */
    bool shakeUseSOR = false;

    //! Number of iterations used to compute the inverse matrix.
    int lincsNIter = 1;
    //! The order for algorithm that adjusts the direction of the bond after constraints are applied.
    int lincslincsExpansionOrder = 4;
    //! The threshold value for the change in bond angle. When exceeded the program will issue a warning
    real lincsWarnAngle = 30.0;
};

//! Helper function to convert ConstraintsTestSystem into string and make test failure messages readable
void PrintTo(const ConstraintsTestSystem& constraintsTestSystem, std::ostream* os)
{
    *os << constraintsTestSystem.title << " - " << constraintsTestSystem.numAtoms << " atoms";
}

const std::vector<ConstraintsTestSystem> c_constraintsTestSystemList = []
{
    std::vector<ConstraintsTestSystem> constraintsTestSystemList;
    {
        ConstraintsTestSystem constraintsTestSystem;

        constraintsTestSystem.title    = "one constraint (e.g. OH)";
        constraintsTestSystem.numAtoms = 2;

        constraintsTestSystem.masses        = { 1.0, 12.0 };
        constraintsTestSystem.constraints   = { 0, 0, 1 };
        constraintsTestSystem.constraintsR0 = { 0.1 };

        real oneTenthOverSqrtTwo = 0.1_real / std::sqrt(2.0_real);

        constraintsTestSystem.x = { { 0.0, oneTenthOverSqrtTwo, 0.0 }, { oneTenthOverSqrtTwo, 0.0, 0.0 } };
        constraintsTestSystem.xPrime = { { 0.01, 0.08, 0.01 }, { 0.06, 0.01, -0.01 } };
        constraintsTestSystem.v      = { { 1.0, 2.0, 3.0 }, { 3.0, 2.0, 1.0 } };

        constraintsTestSystemList.emplace_back(constraintsTestSystem);
    }
    {
        ConstraintsTestSystem constraintsTestSystem;

        constraintsTestSystem.title         = "two disjoint constraints";
        constraintsTestSystem.numAtoms      = 4;
        constraintsTestSystem.masses        = { 0.5, 1.0 / 3.0, 0.25, 1.0 };
        constraintsTestSystem.constraints   = { 0, 0, 1, 1, 2, 3 };
        constraintsTestSystem.constraintsR0 = { 2.0, 1.0 };


        constraintsTestSystem.x = { { 2.50, -3.10, 15.70 },
                                    { 0.51, -3.02, 15.55 },
                                    { -0.50, -3.00, 15.20 },
                                    { -1.51, -2.95, 15.05 } };

        constraintsTestSystem.xPrime = { { 2.50, -3.10, 15.70 },
                                         { 0.51, -3.02, 15.55 },
                                         { -0.50, -3.00, 15.20 },
                                         { -1.51, -2.95, 15.05 } };

        constraintsTestSystem.v = { { 0.0, 1.0, 0.0 }, { 1.0, 0.0, 0.0 }, { 0.0, 0.0, 1.0 }, { 0.0, 0.0, 0.0 } };

        constraintsTestSystemList.emplace_back(constraintsTestSystem);
    }
    {
        ConstraintsTestSystem constraintsTestSystem;

        constraintsTestSystem.title         = "three atoms, connected longitudinally (e.g. CH2)";
        constraintsTestSystem.numAtoms      = 3;
        constraintsTestSystem.masses        = { 1.0, 12.0, 16.0 };
        constraintsTestSystem.constraints   = { 0, 0, 1, 1, 1, 2 };
        constraintsTestSystem.constraintsR0 = { 0.1, 0.2 };

        real oneTenthOverSqrtTwo    = 0.1_real / std::sqrt(2.0_real);
        real twoTenthsOverSqrtThree = 0.2_real / std::sqrt(3.0_real);

        constraintsTestSystem.x = { { oneTenthOverSqrtTwo, oneTenthOverSqrtTwo, 0.0 },
                                    { 0.0, 0.0, 0.0 },
                                    { twoTenthsOverSqrtThree, twoTenthsOverSqrtThree, twoTenthsOverSqrtThree } };

        constraintsTestSystem.xPrime = { { 0.08, 0.07, 0.01 }, { -0.02, 0.01, -0.02 }, { 0.10, 0.12, 0.11 } };

        constraintsTestSystem.v = { { 1.0, 0.0, 0.0 }, { 0.0, 1.0, 0.0 }, { 0.0, 0.0, 1.0 } };

        constraintsTestSystemList.emplace_back(constraintsTestSystem);
    }
    {
        ConstraintsTestSystem constraintsTestSystem;

        constraintsTestSystem.title         = "four atoms, connected longitudinally";
        constraintsTestSystem.numAtoms      = 4;
        constraintsTestSystem.masses        = { 0.5, 1.0 / 3.0, 0.25, 1.0 };
        constraintsTestSystem.constraints   = { 0, 0, 1, 1, 1, 2, 2, 2, 3 };
        constraintsTestSystem.constraintsR0 = { 2.0, 1.0, 1.0 };


        constraintsTestSystem.x = { { 2.50, -3.10, 15.70 },
                                    { 0.51, -3.02, 15.55 },
                                    { -0.50, -3.00, 15.20 },
                                    { -1.51, -2.95, 15.05 } };

        constraintsTestSystem.xPrime = { { 2.50, -3.10, 15.70 },
                                         { 0.51, -3.02, 15.55 },
                                         { -0.50, -3.00, 15.20 },
                                         { -1.51, -2.95, 15.05 } };

        constraintsTestSystem.v = {
            { 0.0, 0.0, 2.0 }, { 0.0, 0.0, 3.0 }, { 0.0, 0.0, -4.0 }, { 0.0, 0.0, -1.0 }
        };

        // Overriding default values since LINCS converges slowly for this system.
        constraintsTestSystem.lincsNIter               = 4;
        constraintsTestSystem.lincslincsExpansionOrder = 8;

        constraintsTestSystemList.emplace_back(constraintsTestSystem);
    }
    {
        ConstraintsTestSystem constraintsTestSystem;

        constraintsTestSystem.title       = "three atoms, connected to the central atom (e.g. CH3)";
        constraintsTestSystem.numAtoms    = 4;
        constraintsTestSystem.masses      = { 12.0, 1.0, 1.0, 1.0 };
        constraintsTestSystem.constraints = { 0, 0, 1, 0, 0, 2, 0, 0, 3 };
        constraintsTestSystem.constraintsR0 = { 0.1 };


        constraintsTestSystem.x = {
            { 0.00, 0.00, 0.00 }, { 0.10, 0.00, 0.00 }, { 0.00, -0.10, 0.00 }, { 0.00, 0.00, 0.10 }
        };

        constraintsTestSystem.xPrime = { { 0.004, 0.009, -0.010 },
                                         { 0.110, -0.006, 0.003 },
                                         { -0.007, -0.102, -0.007 },
                                         { -0.005, 0.011, 0.102 } };

        constraintsTestSystem.v = { { 1.0, 0.0, 0.0 }, { 1.0, 0.0, 0.0 }, { 1.0, 0.0, 0.0 }, { 1.0, 0.0, 0.0 } };

        constraintsTestSystemList.emplace_back(constraintsTestSystem);
    }
    {
        ConstraintsTestSystem constraintsTestSystem;

        constraintsTestSystem.title       = "basic triangle (three atoms, connected to each other)";
        constraintsTestSystem.numAtoms    = 3;
        constraintsTestSystem.masses      = { 1.0, 1.0, 1.0 };
        constraintsTestSystem.constraints = { 0, 0, 1, 2, 0, 2, 1, 1, 2 };
        constraintsTestSystem.hasTriangle = true;
        constraintsTestSystem.constraintsR0 = { 0.1, 0.1, 0.1 };

        real oneTenthOverSqrtTwo = 0.1_real / std::sqrt(2.0_real);

        constraintsTestSystem.x = { { oneTenthOverSqrtTwo, 0.0, 0.0 },
                                    { 0.0, oneTenthOverSqrtTwo, 0.0 },
                                    { 0.0, 0.0, oneTenthOverSqrtTwo } };

        constraintsTestSystem.xPrime = { { 0.09, -0.02, 0.01 }, { -0.02, 0.10, -0.02 }, { 0.03, -0.01, 0.07 } };

        constraintsTestSystem.v = { { 1.0, 1.0, 1.0 }, { -2.0, -2.0, -2.0 }, { 1.0, 1.0, 1.0 } };

        constraintsTestSystemList.emplace_back(constraintsTestSystem);
    }

    {
        ConstraintsTestSystem singleMolecule;

        singleMolecule.title         = "three atoms, connected longitudinally (e.g. CH2)";
        singleMolecule.numAtoms      = 3;
        singleMolecule.masses        = { 1.0, 12.0, 16.0 };
        singleMolecule.constraints   = { 0, 0, 1, 1, 1, 2 };
        singleMolecule.constraintsR0 = { 0.1, 0.2 };

        real oneTenthOverSqrtTwo    = 0.1_real / std::sqrt(2.0_real);
        real twoTenthsOverSqrtThree = 0.2_real / std::sqrt(3.0_real);

        singleMolecule.x = { { oneTenthOverSqrtTwo, oneTenthOverSqrtTwo, 0.0 },
                             { 0.0, 0.0, 0.0 },
                             { twoTenthsOverSqrtThree, twoTenthsOverSqrtThree, twoTenthsOverSqrtThree } };

        singleMolecule.xPrime = { { 0.08, 0.07, 0.01 }, { -0.02, 0.01, -0.02 }, { 0.10, 0.12, 0.11 } };

        singleMolecule.v = { { 1.0, 0.0, 0.0 }, { 0.0, 1.0, 0.0 }, { 0.0, 0.0, 1.0 } };

        // 150 molecules is 300 constraints, which is larger than the thread block of 256 we are currently using
        const int numMolecules = 150;

        ConstraintsTestSystem constraintsTestSystem;

        constraintsTestSystem.title    = "system of many molecules";
        constraintsTestSystem.numAtoms = numMolecules * singleMolecule.numAtoms;

        constraintsTestSystem.masses.resize(numMolecules * singleMolecule.numAtoms);
        constraintsTestSystem.x.resize(numMolecules * singleMolecule.numAtoms);
        constraintsTestSystem.xPrime.resize(numMolecules * singleMolecule.numAtoms);
        constraintsTestSystem.v.resize(numMolecules * singleMolecule.numAtoms);

        constraintsTestSystem.constraints.resize(numMolecules * singleMolecule.constraints.size());

        for (int mol = 0; mol < numMolecules; mol++)
        {
            int offsetAtoms = mol * singleMolecule.numAtoms;
            for (int i = 0; i < singleMolecule.numAtoms; i++)
            {
                constraintsTestSystem.masses[offsetAtoms + i] = singleMolecule.masses[i];
                for (int d = 0; d < DIM; d++)
                {
                    constraintsTestSystem.x[offsetAtoms + i][d]      = singleMolecule.x[i][d];
                    constraintsTestSystem.xPrime[offsetAtoms + i][d] = singleMolecule.xPrime[i][d];
                    constraintsTestSystem.v[offsetAtoms + i][d]      = singleMolecule.v[i][d];
                }
            }
            int offsetConstraints = mol * singleMolecule.constraints.size();
            constraintsTestSystem.constraints[offsetConstraints + 0] = singleMolecule.constraints[0];
            constraintsTestSystem.constraints[offsetConstraints + 1] =
                    singleMolecule.constraints[1] + offsetAtoms;
            constraintsTestSystem.constraints[offsetConstraints + 2] =
                    singleMolecule.constraints[2] + offsetAtoms;
            constraintsTestSystem.constraints[offsetConstraints + 3] = singleMolecule.constraints[3];
            constraintsTestSystem.constraints[offsetConstraints + 4] =
                    singleMolecule.constraints[4] + offsetAtoms;
            constraintsTestSystem.constraints[offsetConstraints + 5] =
                    singleMolecule.constraints[5] + offsetAtoms;
        }

        constraintsTestSystem.constraintsR0.resize(singleMolecule.constraintsR0.size());
        for (unsigned long i = 0; i < singleMolecule.constraintsR0.size(); i++)
        {
            constraintsTestSystem.constraintsR0[i] = singleMolecule.constraintsR0[i];
        }

        constraintsTestSystemList.emplace_back(constraintsTestSystem);
    }

    return constraintsTestSystemList;
}();

//! Input configuration for Constraints tests (defines the physical scenario: constraint topology, PBC)
using ConstraintsInputConfig = std::tuple<ConstraintsTestSystem, t_pbc>;

/*! \brief Hardware test helper for Constraints
 *
 * \todo There are no execution modes - tests don't vary constraint
 * algorithm, but perhaps they should test LINCS vs SHAKE here rather
 * than reset test data. They should test SIMD vs no SIMD here. */
using ConstraintsTestHelper = HardwareAndExecutionTestHelper<ConstraintsInputConfig, std::tuple<>>;

//! Format PBC type for test names
std::string formatPbcType(const t_pbc& pbc)
{
    return c_pbcTypeNames[pbc.pbcType];
}

//! Format ConstraintsTestSystem for test names
std::string formatConstraintSystem(const ConstraintsTestSystem& system)
{
    // Use a sanitized version of the title
    std::string name = system.title;
    // Replace spaces with underscores and remove special characters
    std::replace(name.begin(), name.end(), ' ', '_');
    std::replace(name.begin(), name.end(), '(', '_');
    std::replace(name.begin(), name.end(), ')', '_');
    std::replace(name.begin(), name.end(), '.', '_');
    std::replace(name.begin(), name.end(), ',', '_');
    // Remove consecutive underscores
    name.erase(
            std::unique(name.begin(), name.end(), [](char a, char b) { return a == '_' && b == '_'; }),
            name.end());
    // Remove trailing underscores
    while (!name.empty() && name.back() == '_')
    {
        name.pop_back();
    }
    return formatString("%datoms_%s", system.numAtoms, name.c_str());
}

//! Formatters for parameters in the config info
static const auto sc_configInfoFormatters = std::make_tuple(formatConstraintSystem, formatPbcType);
//! Formatters for parameters in the execution mode (currently empty)
static const auto sc_executionModeFormatters = std::make_tuple();

//! Helper object to name tests using all parameters
static const NameOfTestFromTuple<ConstraintsTestHelper::DynamicParameters> sc_testNamer =
        ConstraintsTestHelper::testNamer(sc_configInfoFormatters, sc_executionModeFormatters);

//! Helper class for checking constraints against reference data
class ConstraintsVerifier
{
public:
    explicit ConstraintsVerifier(TestReferenceChecker* checker) : checker_(checker) {}

private:
    /*! \brief Test if the final position correspond to the reference data.
     *
     * \param[in] testData        Test data structure.
     */
    void checkFinalPositions(const ConstraintsTestData& testData)
    {
        TestReferenceChecker finalPositionsRef(
                checker_->checkSequenceCompound("FinalPositions", testData.numAtoms_));
        for (int i = 0; i < testData.numAtoms_; i++)
        {
            TestReferenceChecker xPrimeRef(finalPositionsRef.checkCompound("Atom", nullptr));
            const gmx::RVec&     xPrime = testData.xPrime_[i];
            xPrimeRef.checkReal(xPrime[XX], "XX");
            xPrimeRef.checkReal(xPrime[YY], "YY");
            xPrimeRef.checkReal(xPrime[ZZ], "ZZ");
        }
    }

    /*! \brief Test if the final velocities correspond to the reference data.
     *
     * \param[in] testData        Test data structure.
     */
    void checkFinalVelocities(const ConstraintsTestData& testData)
    {
        TestReferenceChecker finalVelocitiesRef(
                checker_->checkSequenceCompound("FinalVelocities", testData.numAtoms_));
        for (int i = 0; i < testData.numAtoms_; i++)
        {
            TestReferenceChecker vRef(finalVelocitiesRef.checkCompound("Atom", nullptr));
            const gmx::RVec&     v = testData.v_[i];
            vRef.checkReal(v[XX], "XX");
            vRef.checkReal(v[YY], "YY");
            vRef.checkReal(v[ZZ], "ZZ");
        }
    }

    /*! \brief
     * The test on the final length of constrained bonds.
     *
     * Goes through all the constraints and checks if the final length of all the constraints is
     * equal to the target length with provided tolerance.
     *
     * \param[in] tolerance       Allowed tolerance in final lengths.
     * \param[in] testData        Test data structure.
     * \param[in] pbc             Periodic boundary data.
     */
    static void checkConstrainsLength(FloatingPointTolerance     tolerance,
                                      const ConstraintsTestData& testData,
                                      t_pbc                      pbc)
    {

        // Test if all the constraints are satisfied
        for (Index c = 0; c < gmx::ssize(testData.constraints_) / 3; c++)
        {
            real r0 = testData.constraintsR0_.at(testData.constraints_.at(3 * c));
            int  i  = testData.constraints_.at(3 * c + 1);
            int  j  = testData.constraints_.at(3 * c + 2);
            RVec xij0, xij1;
            real d0, d1;
            if (pbc.pbcType == PbcType::Xyz)
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
                    "rij = %f, which is not equal to r0 = %f for constraint #%zd, between atoms %d "
                    "and %d"
                    " (before constraining rij was %f).",
                    d1,
                    r0,
                    c,
                    i,
                    j,
                    d0);
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
    static void checkConstrainsDirection(const ConstraintsTestData& testData, t_pbc pbc)
    {

        for (Index c = 0; c < gmx::ssize(testData.constraints_) / 3; c++)
        {
            int  i = testData.constraints_.at(3 * c + 1);
            int  j = testData.constraints_.at(3 * c + 2);
            RVec xij0, xij1;
            if (pbc.pbcType == PbcType::Xyz)
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
                    "The constraint %zd changed direction. Constraining algorithm might have "
                    "returned the wrong root "
                    "of the constraints equation.",
                    c);
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
    static void checkCOMCoordinates(FloatingPointTolerance tolerance, const ConstraintsTestData& testData)
    {

        RVec comPrime0({ 0.0, 0.0, 0.0 });
        RVec comPrime({ 0.0, 0.0, 0.0 });
        for (int i = 0; i < testData.numAtoms_; i++)
        {
            comPrime0 += testData.masses_[i] * testData.xPrime0_[i];
            comPrime += testData.masses_[i] * testData.xPrime_[i];
        }

        comPrime0 /= testData.numAtoms_;
        comPrime /= testData.numAtoms_;

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
    static void checkCOMVelocity(FloatingPointTolerance tolerance, const ConstraintsTestData& testData)
    {

        RVec comV0({ 0.0, 0.0, 0.0 });
        RVec comV({ 0.0, 0.0, 0.0 });
        for (int i = 0; i < testData.numAtoms_; i++)
        {
            comV0 += testData.masses_[i] * testData.v0_[i];
            comV += testData.masses_[i] * testData.v_[i];
        }
        comV0 /= testData.numAtoms_;
        comV /= testData.numAtoms_;

        EXPECT_REAL_EQ_TOL(comV[XX], comV0[XX], tolerance)
                << "Velocity of the center of mass in x-direction has been changed by constraints.";
        EXPECT_REAL_EQ_TOL(comV[YY], comV0[YY], tolerance)
                << "Velocity of the center of mass in y-direction has been changed by constraints.";
        EXPECT_REAL_EQ_TOL(comV[ZZ], comV0[ZZ], tolerance)
                << "Velocity of the center of mass in z-direction has been changed by constraints.";
    }

    /*! \brief
     * The test of virial tensor.
     *
     * Checks if the values in the scaled virial tensor are equal to pre-computed values.
     *
     * \param[in] testData        Test data structure.
     */
    void checkVirialTensor(const ConstraintsTestData& testData)
    {
        const tensor&        virialScaled = testData.virialScaled_;
        TestReferenceChecker virialScaledRef(checker_->checkCompound("VirialScaled", nullptr));

        virialScaledRef.checkReal(virialScaled[XX][XX], "XX");
        virialScaledRef.checkReal(virialScaled[XX][YY], "XY");
        virialScaledRef.checkReal(virialScaled[XX][ZZ], "XZ");
        virialScaledRef.checkReal(virialScaled[YY][XX], "YX");
        virialScaledRef.checkReal(virialScaled[YY][YY], "YY");
        virialScaledRef.checkReal(virialScaled[YY][ZZ], "YZ");
        virialScaledRef.checkReal(virialScaled[ZZ][XX], "ZX");
        virialScaledRef.checkReal(virialScaled[ZZ][YY], "ZY");
        virialScaledRef.checkReal(virialScaled[ZZ][ZZ], "ZZ");
    }

public:
    /*! \brief The workflow for checking whether constraints have been satisfied
     *
     * \param[in] testData   Test data to verify
     * \param[in] algorithm  Constraint algorithm used (SHAKE or LINCS)
     * \param[in] pbc        Periodic boundary conditions
     */
    void check(const ConstraintsTestData& testData, const ConstraintAlgorithm algorithm, t_pbc pbc)
    {
        const FloatingPointTolerance positionsTolerance = absoluteTolerance(0.001F);
        const FloatingPointTolerance velocityTolerance  = absoluteTolerance(0.02F);
        const FloatingPointTolerance lengthTolerance = relativeToleranceAsFloatingPoint(0.1, 0.002F);

        checker_->setDefaultTolerance(positionsTolerance);
        checkFinalPositions(testData);
        checker_->setDefaultTolerance(velocityTolerance);
        checkFinalVelocities(testData);

        checkConstrainsLength(lengthTolerance, testData, pbc);
        checkConstrainsDirection(testData, pbc);
        checkCOMCoordinates(positionsTolerance, testData);
        checkCOMVelocity(velocityTolerance, testData);

        float virialTrace = 0.0F;
        for (int d = 0; d < DIM; d++)
        {
            virialTrace += testData.virialScaled_[d][d];
        }

        // The virial tolerance factor can be:
        // LINCS iter=2, order=4:   0.002
        // LINCS iter=1, order=4:   0.02
        // SHAKE tolerance=0.0001:  0.2
        // SHAKE tolerance=0.00002: 0.1
        const float virialRelativeTolerance = (algorithm == ConstraintAlgorithm::Shake ? 0.1F : 0.02F);
        FloatingPointTolerance virialTolerance =
                absoluteTolerance(fabs(virialTrace) / 3 * virialRelativeTolerance);

        checker_->setDefaultTolerance(virialTolerance);
        checkVirialTensor(testData);
    }

private:
    //! Reference data checker used for all verifications
    TestReferenceChecker* checker_;
};


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
class ConstraintsTest : public HardwareTestFixture<ConstraintsTestHelper>
{
protected:
    ConstraintsTest() : HardwareTestFixture(sc_configInfoFormatters) {}
};

TEST_P(ConstraintsTest, SatisfiesConstraints)
{
    auto [constraintsTestSystem, pbc, _] = GetParam();

    ConstraintsTestData testData(constraintsTestSystem.title,
                                 constraintsTestSystem.numAtoms,
                                 constraintsTestSystem.masses,
                                 constraintsTestSystem.constraints,
                                 constraintsTestSystem.constraintsR0,
                                 true,
                                 false,
                                 real(0.0),
                                 real(0.001),
                                 constraintsTestSystem.x,
                                 constraintsTestSystem.xPrime,
                                 constraintsTestSystem.v,
                                 constraintsTestSystem.shakeTolerance,
                                 constraintsTestSystem.shakeUseSOR,
                                 constraintsTestSystem.lincsNIter,
                                 constraintsTestSystem.lincslincsExpansionOrder,
                                 constraintsTestSystem.lincsWarnAngle);

    testData.reset();

    // Apply constraints based on hardware context
    if (isGpuTest())
    {
#if GMX_GPU && !GMX_GPU_OPENCL
        SCOPED_TRACE(formatString("Testing %s with %s PBC using %s (LINCS).",
                                  testData.title_.c_str(),
                                  c_pbcTypeNames[pbc.pbcType].c_str(),
                                  hardwareContext()->description().c_str()));

        activateHardware();
        applyLincsGpu(*deviceContext(), *deviceStream(), &testData, pbc);

        ConstraintsVerifier verifier(&checker());
        verifier.check(testData, ConstraintAlgorithm::Lincs, pbc);
#else
        GMX_THROW(gmx::InternalError("GPU hardware context with no test code"));
#endif
    }
    else
    {
        // Run both SHAKE and LINCS on CPU
        {
            SCOPED_TRACE(formatString("Testing %s with %s PBC using %s (SHAKE).",
                                      testData.title_.c_str(),
                                      c_pbcTypeNames[pbc.pbcType].c_str(),
                                      hardwareContext()->description().c_str()));
            testData.reset();
            applyShakeCpu(&testData, pbc);

            ConstraintsVerifier verifier(&checker());
            verifier.check(testData, ConstraintAlgorithm::Shake, pbc);
        }
        {
            SCOPED_TRACE(formatString("Testing %s with %s PBC using %s (LINCS).",
                                      testData.title_.c_str(),
                                      c_pbcTypeNames[pbc.pbcType].c_str(),
                                      hardwareContext()->description().c_str()));
            testData.reset();
            applyLincsCpu(&testData, pbc);

            ConstraintsVerifier verifier(&checker());
            verifier.check(testData, ConstraintAlgorithm::Lincs, pbc);
        }
    }
}

INSTANTIATE_TEST_SUITE_P(AllHardware,
                         ConstraintsTest,
                         ::testing::Combine(::testing::ValuesIn(c_constraintsTestSystemList),
                                            ::testing::ValuesIn(c_pbcs),
                                            ::testing::ValuesIn(getHardwareContextsWithCapability(
                                                    GpuConfigurationCapabilities::Update))),
                         sc_testNamer);

/*! \brief Test fixture for topology-based constraint tests.
 *
 * These tests check properties of the constraint topology without needing to
 * run the constraint algorithm on different hardware contexts.
 */
class ConstraintsTopologyTest : public ::testing::TestWithParam<ConstraintsTestSystem>
{
};

TEST_P(ConstraintsTopologyTest, TriangleDetectionWorks)
{
    const ConstraintsTestSystem& constraintsTestSystem = GetParam();
    const ConstraintsTestData    testData(constraintsTestSystem.title,
                                       constraintsTestSystem.numAtoms,
                                       constraintsTestSystem.masses,
                                       constraintsTestSystem.constraints,
                                       constraintsTestSystem.constraintsR0,
                                       true,
                                       false,
                                       real(0.0),
                                       real(0.001),
                                       constraintsTestSystem.x,
                                       constraintsTestSystem.xPrime,
                                       constraintsTestSystem.v,
                                       constraintsTestSystem.shakeTolerance,
                                       constraintsTestSystem.shakeUseSOR,
                                       constraintsTestSystem.lincsNIter,
                                       constraintsTestSystem.lincslincsExpansionOrder,
                                       constraintsTestSystem.lincsWarnAngle);

    EXPECT_EQ(constraintsTestSystem.hasTriangle,
              hasTriangleConstraints(testData.mtop_, FlexibleConstraintTreatment::Include));
}

INSTANTIATE_TEST_SUITE_P(TopologyTests,
                         ConstraintsTopologyTest,
                         ::testing::ValuesIn(c_constraintsTestSystemList),
                         [](const ::testing::TestParamInfo<ConstraintsTestSystem>& i)
                         { return formatConstraintSystem(i.param); });

} // namespace
} // namespace test
} // namespace gmx
