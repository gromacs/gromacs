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
 * \brief
 * Implements test of bonded force routines
 *
 *
 * \todo We should re-work this. For example, a harmonic bond
 * has so few computations that force components should be
 * accurate to a small *and computed* relative error.
 *
 * \author David van der Spoel <david.vanderspoel@icm.uu.se>
 * \ingroup module_listed_forces
 */
#include "gmxpre.h"

#include "gromacs/listed_forces/bonded.h"

#include <cmath>
#include <cstdint>

#include <algorithm>
#include <iterator>
#include <limits>
#include <memory>
#include <ostream>
#include <string>
#include <tuple>
#include <vector>

#include <gtest/gtest.h>

#include "gromacs/gpu_utils/capabilities.h"
#include "gromacs/listed_forces/listed_forces.h"
#include "gromacs/listed_forces/listed_forces_gpu.h"
#include "gromacs/math/paddedvector.h"
#include "gromacs/mdtypes/md_enums.h"
#include "gromacs/pbcutil/ishift.h"
#include "gromacs/pbcutil/pbc.h"
#include "gromacs/topology/idef.h"
#include "gromacs/topology/ifunc.h"
#include "gromacs/utility/arrayref.h"
#include "gromacs/utility/real.h"
#include "gromacs/utility/strconvert.h"
#include "gromacs/utility/vec.h"
#include "gromacs/utility/vectypes.h"

#include "testutils/hardware_test_fixture.h"
#include "testutils/naming.h"
#include "testutils/refdata.h"
#include "testutils/test_device.h"
#include "testutils/test_hardware_environment.h"
#include "testutils/testasserts.h"

#include "bondedtestdata.h"

#if GMX_GPU && !GMX_GPU_OPENCL
#    include "gromacs/gpu_utils/devicebuffer.h"
#    include "gromacs/mdtypes/enerdata.h"
#    include "gromacs/mdtypes/simulation_workload.h"
#    include "gromacs/nbnxm/gpu_types_common.h"
#    include "gromacs/topology/forcefieldparameters.h"
#endif

namespace gmx
{
namespace test
{
namespace
{

//! Input configuration (the physical scenario: bonds, coordinates, PBC, FEP lambda)
using ListedForcesInputConfig = std::tuple<iListInput, PaddedVector<RVec>, PbcType, real>;

/*! \brief Execution modes for listed forces tests
 *
 * These control how values are computed: forces only vs
 * forces+virial+energy, SIMD vs not, etc. */
using ListedForcesExecutionModes = std::tuple<BondedKernelFlavor>;

//! Hardware and execution test helper for listed forces
using ListedForcesTestHelper =
        HardwareAndExecutionTestHelper<ListedForcesInputConfig, ListedForcesExecutionModes>;

void runBondedCpu(const iListInput&           input,
                  const std::vector<t_iatom>& iatoms,
                  const PaddedVector<RVec>&   x,
                  const t_pbc*                pbc,
                  real                        lambda,
                  BondedKernelFlavor          flavor,
                  OutputQuantities*           output)
{
    std::array<int, 4> ddgatindex = { 0, 1, 2, 3 };

    output->energy = calculateSimpleBond(input.ftype.value(),
                                         iatoms.size(),
                                         iatoms.data(),
                                         &input.iparams,
                                         as_rvec_array(x.data()),
                                         output->f,
                                         output->fshift,
                                         pbc,
                                         lambda,
                                         &output->dvdlambda,
                                         c_bondedTestCharges,
                                         nullptr,
                                         nullptr,
                                         nullptr,
                                         ddgatindex.data(),
                                         flavor);
}

#if GMX_GPU && !GMX_GPU_OPENCL

void runBondedGpu(const DeviceContext&        deviceContext,
                  const DeviceStream&         deviceStream,
                  const iListInput&           input,
                  const std::vector<t_iatom>& iatoms,
                  const PaddedVector<RVec>&   x,
                  const t_pbc&                pbc,
                  BondedKernelFlavor          flavor,
                  OutputQuantities*           output)
{
    GMX_RELEASE_ASSERT(flavor == BondedKernelFlavor::ForcesAndVirialAndEnergy,
                       "GPU runner only supports ForcesAndVirialAndEnergy flavor");

    deviceContext.activate();

    const int  ftype    = static_cast<int>(input.ftype.value());
    const int  numAtoms = c_numAtomsBondedTest;
    const real scaleFac = 1.0f;

    // Build minimal force-field parameters (single type)
    gmx_ffparams_t ffparams;
    ffparams.functype = { input.ftype.value() };
    ffparams.iparams  = { input.iparams };

    // Build minimal interaction definitions
    InteractionDefinitions idef(ffparams);
    idef.il[ftype].iatoms                   = iatoms;
    idef.ilsort                             = ilsortNO_FE;
    idef.numNonperturbedInteractions[ftype] = static_cast<int>(iatoms.size());

    // Identity atom order (same as input coordinates)
    std::vector<int> nbnxmAtomOrder(numAtoms);
    for (int i = 0; i < numAtoms; i++)
    {
        nbnxmAtomOrder[i] = i;
    }

    // Host-side xq (x, y, z, q) for upload
    std::vector<real> h_xq;
    h_xq.reserve(numAtoms * 4);
    for (int i = 0; i < numAtoms; i++)
    {
        h_xq.emplace_back(x[i][XX]);
        h_xq.emplace_back(x[i][YY]);
        h_xq.emplace_back(x[i][ZZ]);
        h_xq.emplace_back(c_bondedTestCharges[i]);
    }

    // Allocate device buffers for NBAtomDataGpu (only xq, f, fShift used by bonded)
    NBAtomDataGpu nbAtomDataGpu = {};
    nbAtomDataGpu.numAtoms      = numAtoms;
    nbAtomDataGpu.numAtomsLocal = numAtoms;
    nbAtomDataGpu.numAtomsAlloc = numAtoms;

    allocateDeviceBuffer(&nbAtomDataGpu.xq, numAtoms, deviceContext);
    allocateDeviceBuffer(&nbAtomDataGpu.f, numAtoms, deviceContext);
    allocateDeviceBuffer(&nbAtomDataGpu.fShift, c_numShiftVectors, deviceContext);

    copyToDeviceBuffer(&nbAtomDataGpu.xq,
                       reinterpret_cast<Float4*>(h_xq.data()),
                       0,
                       numAtoms,
                       deviceStream,
                       GpuApiCallBehavior::Sync,
                       nullptr);
    clearDeviceBufferAsync(&nbAtomDataGpu.f, 0, numAtoms, deviceStream);
    clearDeviceBufferAsync(&nbAtomDataGpu.fShift, 0, c_numShiftVectors, deviceStream);

    // ListedForcesGpu (wallcycle can be nullptr – start/stop check for null)
    ListedForcesGpu listedForcesGpu(ffparams, scaleFac, 1, deviceContext, deviceStream, nullptr);

    listedForcesGpu.updateHaveInteractions(idef);
    listedForcesGpu.updateInteractionListsAndDeviceBuffers(nbnxmAtomOrder, idef, &nbAtomDataGpu);

    if (!listedForcesGpu.haveInteractions())
    {
        freeDeviceBuffer(&nbAtomDataGpu.xq);
        freeDeviceBuffer(&nbAtomDataGpu.f);
        freeDeviceBuffer(&nbAtomDataGpu.fShift);
        GMX_THROW(InternalError(formatString("ListedForcesGpu reported no interactions for type %d", ftype)));
    }

    listedForcesGpu.setPbc(pbc.pbcType, pbc.box, false);

    StepWorkload stepWork;
    stepWork.computeVirial = true;
    stepWork.computeEnergy = true;
    listedForcesGpu.launchKernel(stepWork);
    listedForcesGpu.launchEnergyTransfer();

    gmx_enerdata_t enerd(1, nullptr);
    listedForcesGpu.waitAccumulateEnergyTerms(&enerd);

    // Copy forces and shift forces back
    std::vector<Float3> h_f(numAtoms);
    std::vector<Float3> h_fShift(c_numShiftVectors);
    copyFromDeviceBuffer(
            h_f.data(), &nbAtomDataGpu.f, 0, numAtoms, deviceStream, GpuApiCallBehavior::Sync, nullptr);
    copyFromDeviceBuffer(h_fShift.data(),
                         &nbAtomDataGpu.fShift,
                         0,
                         c_numShiftVectors,
                         deviceStream,
                         GpuApiCallBehavior::Sync,
                         nullptr);

    // Fill test output (Float3 is RVec on HIP, uses [0],[1],[2])
    output->energy    = enerd.term[ftype];
    output->dvdlambda = 0.0;
    for (int i = 0; i < numAtoms; i++)
    {
        output->f[i][XX] = h_f[i][XX];
        output->f[i][YY] = h_f[i][YY];
        output->f[i][ZZ] = h_f[i][ZZ];
        output->f[i][3]  = 0.0;
    }
    for (int s = 0; s < c_numShiftVectors; s++)
    {
        output->fshift[s][XX] = h_fShift[s][XX];
        output->fshift[s][YY] = h_fShift[s][YY];
        output->fshift[s][ZZ] = h_fShift[s][ZZ];
    }

    freeDeviceBuffer(&nbAtomDataGpu.xq);
    freeDeviceBuffer(&nbAtomDataGpu.f);
    freeDeviceBuffer(&nbAtomDataGpu.fShift);
}

#endif // GMX_GPU && !GMX_GPU_OPENCL

//! Formatter for iListInput - use function type and FEP status
// NOLINTNEXTLINE(cppcoreguidelines-interfaces-global-init)
auto formatIListInput = [](const iListInput& input)
{
    if (!input.ftype.has_value())
    {
        return std::string("NoFtype");
    }
    // Safe from global-init ordering problems because interaction_function
    // will be defined before this lambda is called.
    std::string name = interaction_function[input.ftype.value()].name;
    if (input.hasFepParameters())
    {
        name.append("_FEP");
    }
    return name;
};

//! Formatter for coordinates - use atom count
auto formatCoordinates = [](const PaddedVector<RVec>& coords)
{ return formatString("%datoms", static_cast<int>(coords.size())); };

//! Formatter for PbcType
// NOLINTNEXTLINE(cppcoreguidelines-interfaces-global-init)
auto formatPbcType = [](PbcType pbc)
{
    // Safe from global-init ordering problems because c_pbcTypeNames
    // will be defined before this lambda is called.
    return std::string(c_pbcTypeNames[pbc]);
};

//! Formatter for BondedKernelFlavor - creates valid GTest parameter names
auto formatFlavor = [](BondedKernelFlavor flavor)
{
    static constexpr gmx::EnumerationArray<BondedKernelFlavor, const char*> flavorNames = {
        "ForcesSimdWhenAvailable", "ForcesNoSimd", "ForcesAndVirialAndEnergy", "ForcesAndEnergy"
    };
    return flavorNames[flavor];
};

//! Formatter for FEP lambda parameter
auto formatFepLambda = [](real lambda) { return formatString("Lambda%s", toString(lambda).c_str()); };

//! Formatters for parameters in the config info
static const auto sc_configInfoFormatters =
        std::make_tuple(formatIListInput, formatCoordinates, formatPbcType, formatFepLambda);
//! Formatters for parameters in the execution mode
static const auto sc_executionModeFormatters = std::make_tuple(formatFlavor);

/*! \brief Helper object to name tests using all parameters
 *
 * Test names include hardware and flavor:
 *   Bond/ListedForcesTest.Ifunc/BONDS_4atoms_no_Lambda0_ForcesAndVirialAndEnergy_GPU1
 */
static const NameOfTestFromTuple<ListedForcesTestHelper::DynamicParameters> sc_testNamer =
        ListedForcesTestHelper::testNamer(sc_configInfoFormatters, sc_executionModeFormatters);

/*! \brief Test fixture for listed forces
 *
 * Uses hardware-independent reference-data naming. The fixture automatically
 * creates the RefDataFilenameMaker from sc_configInfoFormatters, so test code
 * doesn't need to manually construct it.
 */
class ListedForcesTest : public HardwareTestFixture<ListedForcesTestHelper>
{
protected:
    iListInput             input_;
    matrix                 box_;
    t_pbc                  pbc_;
    PaddedVector<RVec>     x_;
    PbcType                pbcType_;
    real                   lambda_;
    BondedKernelFlavor     flavor_;
    FloatingPointTolerance shiftForcesTolerance_ = defaultRealTolerance();

    ListedForcesTest() :
        HardwareTestFixture(sc_configInfoFormatters),
        input_(std::get<0>(GetParam())),
        x_(std::get<1>(GetParam())),
        pbcType_(std::get<2>(GetParam())),
        lambda_(std::get<3>(GetParam())),
        flavor_(std::get<4>(GetParam()))
    {
        clear_mat(box_);
        box_[0][0] = box_[1][1] = box_[2][2] = 1.5;
        set_pbc(&pbc_, pbcType_, box_);
        // The SIMD acos() is only accurate to 2-3 ULP, so the angles
        // computed by it and the non-SIMD code paths (that use
        // std::acos) differ by enough to require quite large
        // tolerances for the shift forces in mixed precision.
        float singleShiftForcesAbsoluteTolerance =
                ((input_.ftype == InteractionFunction::Polarization)
                                 || (input_.ftype == InteractionFunction::AnharmonicPolarization)
                                 || (IS_ANGLE(input_.ftype.value()))
                         ? 5e-3
                         : 5e-5);
        // Note that std::numeric_limits isn't required by the standard to
        // have an implementation for uint64_t(!) but this is likely to
        // work because that type is likely to be a typedef for one of
        // the other numerical types that happens to be 64-bits wide.
        shiftForcesTolerance_ = FloatingPointTolerance(singleShiftForcesAbsoluteTolerance,
                                                       1e-8,
                                                       1e-6,
                                                       1e-12,
                                                       std::numeric_limits<uint64_t>::max(),
                                                       std::numeric_limits<uint64_t>::max(),
                                                       false);
    }

    void testOneIfunc(TestReferenceChecker* checker, const std::vector<t_iatom>& iatoms)
    {
        SCOPED_TRACE(std::string("Testing PBC type: ") + c_pbcTypeNames[pbcType_] + ", on "
                     + hardwareContext()->description() + " with bonded kernel flavor: "
                     + c_bondedKernelFlavorStrings[flavor_] + ", lambda: " + toString(lambda_));

        OutputQuantities output;

        if (isGpuTest())
        {
#if GMX_GPU && !GMX_GPU_OPENCL
            activateHardware();
            runBondedGpu(*deviceContext(), *deviceStream(), input_, iatoms, x_, pbc_, flavor_, &output);
#else
            GMX_THROW(gmx::InternalError("GPU hardware context with no test code"));
#endif
        }
        else
        {
            runBondedCpu(input_, iatoms, x_, &pbc_, lambda_, flavor_, &output);
        }

        // Now check the output
        EXPECT_TRUE((input_.hasFepParameters() || (output.dvdlambda == 0.0)))
                << "dvdlambda was " << output.dvdlambda;
        checkOutput(checker, output, flavor_);
        auto shiftForcesChecker = checker->checkCompound("Shift-Forces", "Shift-forces");
        if (computeVirial(flavor_))
        {
            shiftForcesChecker.setDefaultTolerance(shiftForcesTolerance_);
            shiftForcesChecker.checkVector(output.fshift[c_centralShiftIndex], "Central");
        }
        else
        {
            shiftForcesChecker.disableUnusedEntriesCheck();
        }
    }
    void testIfunc()
    {
        EXPECT_TRUE(input_.ftype.has_value() && input_.ftype < InteractionFunction::Count);

        // Activate device if GPU
        activateHardware();

        // We need quite specific tolerances here since angle functions
        // etc. are not very precise and reproducible.
        test::FloatingPointTolerance tolerance(test::FloatingPointTolerance(
                input_.ftoler, input_.dtoler, 1.0e-6, 1.0e-12, 10000, 100, false));
        checker().setDefaultTolerance(tolerance);

        TestReferenceChecker thisChecker = checker().checkCompound(
                "FunctionType", interaction_function[input_.ftype.value()].name);

        std::vector<t_iatom> iatoms;
        fillIatoms(input_.ftype, &iatoms);

        testOneIfunc(&thisChecker, iatoms);
    }

    //! Override to add parameter-dependent skip logic
    void addCustomSkipReasons(MessageStringCollector& skipReasons) override
    {
        const bool isFep = input_.hasFepParameters();
        const bool isGpu = isGpuTest();

        // CPU + GPU shared skip conditions
        skipReasons.appendIf(!isFep && lambda_ != 0.0, "Non-FEP inputs only support lambda=0.0");
        skipReasons.appendIf(isFep && flavor_ == BondedKernelFlavor::ForcesSimdWhenAvailable,
                             "FEP does not support ForcesSimdWhenAvailable flavor");

        // GPU-specific skip conditions
        if (isGpu)
        {
            // TODO forces-only should also be tested on GPUs
            skipReasons.appendIf(flavor_ != BondedKernelFlavor::ForcesAndVirialAndEnergy,
                                 "GPU bonded only supports ForcesAndVirialAndEnergy flavor");
            skipReasons.appendIf(isFep,
                                 "GPU bonded does not support free-energy perturbation (FEP)");
            skipReasons.appendIf(lambda_ != 0.0, "GPU does not support lambda != 0.0");
            skipReasons.appendIf(!input_.ftype.has_value(), "Interaction type not specified");
            skipReasons.appendIf(
                    input_.ftype.has_value()
                            && std::find(fTypesOnGpu.begin(), fTypesOnGpu.end(), input_.ftype.value())
                                       == fTypesOnGpu.end(),
                    formatString("Interaction type '%s' not implemented on GPU",
                                 interaction_function[input_.ftype.value()].name));
        }
        // CPU-specific skip conditions
        else
        {
            skipReasons.appendIf(flavor_ == BondedKernelFlavor::ForcesSimdWhenAvailable && lambda_ != 0.0,
                                 formatString("CPU does not support %s with lambda=%.1f",
                                              c_bondedKernelFlavorStrings[flavor_].c_str(),
                                              lambda_));
        }
    }
};

TEST_P(ListedForcesTest, Ifunc)
{
    // When updating reference data, skip forces-only flavor tests since they don't compute
    // energy fields. Reference data should be generated by ForcesAndVirialAndEnergy flavor.
    if (!computeEnergy(flavor_) && referenceDataMode() != ReferenceDataMode::Compare)
    {
        GTEST_SKIP() << "Skipping forces-only flavor when updating reference data (use "
                        "ForcesAndVirialAndEnergy for canonical data)";
    }
    testIfunc();
}

//! Function types for testing bonds. Add new terms at the end.
std::vector<iListInput> c_InputBonds = {
    { iListInput(2e-6F, 1e-8).setHarmonic(InteractionFunction::Bonds, 0.15, 500.0) },
    { iListInput(2e-6F, 1e-8).setHarmonic(InteractionFunction::Bonds, 0.15, 500.0, 0.17, 400.0) },
    { iListInput(1e-4F, 1e-8).setHarmonic(InteractionFunction::GROMOS96Bonds, 0.15, 50.0) },
    { iListInput().setHarmonic(InteractionFunction::GROMOS96Bonds, 0.15, 50.0, 0.17, 40.0) },
    { iListInput().setCubic(0.16, 50.0, 2.0) },
    { iListInput(2e-6F, 1e-8).setMorse(0.15, 50.0, 2.0, 0.17, 40.0, 1.6) },
    { iListInput(2e-6F, 1e-8).setMorse(0.15, 30.0, 2.7) },
    { iListInput().setFene(0.4, 5.0) }
};

//! Constants for Quartic Angles
const real cQuarticAngles[5] = { 1.1, 2.3, 4.6, 7.8, 9.2 };

//! Function types for testing angles. Add new terms at the end.
std::vector<iListInput> c_InputAngles = {
    { iListInput(2e-3, 1e-8).setHarmonic(InteractionFunction::Angles, 100.0, 50.0) },
    { iListInput(2e-3, 1e-8).setHarmonic(InteractionFunction::Angles, 100.15, 50.0, 95.0, 30.0) },
    { iListInput(8e-3, 1e-8).setHarmonic(InteractionFunction::GROMOS96Angles, 100.0, 50.0) },
    { iListInput(8e-3, 1e-8).setHarmonic(InteractionFunction::GROMOS96Angles, 100.0, 50.0, 95.0, 30.0) },
    { iListInput().setLinearAngle(50.0, 0.4) },
    { iListInput().setLinearAngle(50.0, 0.4, 40.0, 0.6) },
    { iListInput(2e-6, 1e-8).setCrossBondBonds(0.8, 0.7, 45.0) },
    { iListInput(3e-6, 1e-8).setCrossBondAngles(0.8, 0.7, 0.3, 45.0) },
    { iListInput(2e-2, 1e-8).setUreyBradley(950.0, 46.0, 0.3, 5.0) },
    { iListInput(2e-2, 1e-8).setUreyBradley(100.0, 45.0, 0.3, 5.0, 90.0, 47.0, 0.32, 7.0) },
    { iListInput(2e-3, 1e-8).setQuarticAngles(87.0, cQuarticAngles) }
};

//! Constants for Ryckaert-Bellemans A
const real rbcA[NR_RBDIHS] = { -5.35, 13.6, 8.4, -16.7, 0.3, 12.4 };

//! Constants for Ryckaert-Bellemans B
const real rbcB[NR_RBDIHS] = { -6.35, 12.6, 8.1, -10.7, 0.9, 15.4 };

//! Constants for Ryckaert-Bellemans without FEP
const real rbc[NR_RBDIHS] = { -7.35, 13.6, 8.4, -16.7, 1.3, 12.4 };

//! Function types for testing dihedrals. Add new terms at the end.
std::vector<iListInput> c_InputDihs = {
    { iListInput(5e-4, 1e-8).setPDihedrals(InteractionFunction::ProperDihedrals, -100.0, 10.0, 2, -80.0, 20.0) },
    { iListInput(1e-4, 1e-8).setPDihedrals(InteractionFunction::ProperDihedrals, -105.0, 15.0, 2) },
    { iListInput(2e-4, 1e-8).setHarmonic(InteractionFunction::ImproperDihedrals, 100.0, 50.0) },
    { iListInput(2e-4, 1e-8).setHarmonic(InteractionFunction::ImproperDihedrals, 100.15, 50.0, 95.0, 30.0) },
    { iListInput(4e-4, 1e-8).setRbDihedrals(rbcA, rbcB) },
    { iListInput(4e-4, 1e-8).setRbDihedrals(rbc) }
};

//! Function types for testing polarization. Add new terms at the end.
std::vector<iListInput> c_InputPols = {
    { iListInput(2e-5, 1e-8).setPolarization(0.12) },
    { iListInput(2e-3, 1e-8).setAnharmPolarization(0.0013, 0.02, 1235.6) },
    { iListInput(1.4e-3, 1e-8).setTholePolarization(0.26, 0.07, 0.09) },
    { iListInput(2e-3, 1e-8).setWaterPolarization(0.001, 0.0012, 0.0016, 0.095, 0.15, 0.02) },
};

//! Function types for testing polarization. Add new terms at the end.
std::vector<iListInput> c_InputRestraints = {
    { iListInput(1e-4, 1e-8).setPDihedrals(InteractionFunction::AngleRestraints, -100.0, 10.0, 2, -80.0, 20.0) },
    { iListInput(1e-4, 1e-8).setPDihedrals(InteractionFunction::AngleRestraints, -105.0, 15.0, 2) },
    { iListInput(1e-4, 1e-8).setPDihedrals(InteractionFunction::AngleZAxisRestraints, -100.0, 10.0, 2, -80.0, 20.0) },
    { iListInput(1e-4, 1e-8).setPDihedrals(InteractionFunction::AngleZAxisRestraints, -105.0, 15.0, 2) },
    { iListInput(2e-3, 1e-8).setHarmonic(InteractionFunction::RestrictedBendingPotential, 100.0, 50.0) },
    { iListInput(2e-3, 1e-8).setHarmonic(InteractionFunction::RestrictedBendingPotential, 100.0, 50.0, 110.0, 45.0) }
};

//! Function types for testing bond with zero length, has zero reference length to make physical sense.
std::vector<iListInput> c_InputBondsZeroLength = {
    { iListInput().setHarmonic(InteractionFunction::Bonds, 0.0, 500.0) },
};

//! Function types for testing angles with zero angle, has zero reference angle to make physical sense.
std::vector<iListInput> c_InputAnglesZeroAngle = {
    { iListInput(2e-3, 1e-8).setHarmonic(InteractionFunction::Angles, 0.0, 50.0) },
};

} // namespace
} // namespace test

//! Print an RVec to \c os
static void PrintTo(const RVec& value, std::ostream* os)
{
    *os << value[XX] << " " << value[YY] << " " << value[ZZ] << std::endl;
}

//! Print a padded vector of RVec to \c os
static void PrintTo(const PaddedVector<RVec>& vector, std::ostream* os)
{
    if (vector.empty())
    {
        *os << "Empty vector" << std::endl;
    }
    else
    {
        *os << "Vector of RVec containing:" << std::endl;
        std::for_each(vector.begin(), vector.end(), [os](const RVec& v) { PrintTo(v, os); });
    }
}

namespace test
{
namespace
{

/*! \brief Coordinates for testing
 *
 * Taken from a butane molecule, so we have some
 * normal-sized bonds and angles to test.
 *
 * \todo Test also some weirder values */
std::vector<PaddedVector<RVec>> c_coordinatesForTests = {
    { { 1.382, 1.573, 1.482 }, { 1.281, 1.559, 1.596 }, { 1.292, 1.422, 1.663 }, { 1.189, 1.407, 1.775 } }
};

//! Coordinates for testing bonds with zero length
std::vector<PaddedVector<RVec>> c_coordinatesForTestsZeroBondLength = {
    { { 0.0, 0.0, 0.0 }, { 0.0, 0.0, 0.0 }, { 0.005, 0.0, 0.1 }, { 0.005, 0.0, 0.1 } }
};

//! Coordinates for testing bonds with zero length
std::vector<PaddedVector<RVec>> c_coordinatesForTestsZeroAngle = {
    { { 0.005, 0.0, 0.1 }, { 0.0, 0.0, 0.0 }, { 0.005, 0.0, 0.1 }, { 0.5, 0.18, 0.22 } }
};

//! PBC values for testing
std::vector<PbcType> c_pbcForTests = { PbcType::No, PbcType::XY, PbcType::Xyz };

//! Bonded kernel flavors for testing
std::vector<BondedKernelFlavor> c_flavorsForTests = { BondedKernelFlavor::ForcesAndVirialAndEnergy,
                                                      BondedKernelFlavor::ForcesSimdWhenAvailable };

//! Lambda values for testing
std::vector<real> c_lambdaValuesForTests = { 0.0, 0.5, 1.0 };

INSTANTIATE_TEST_SUITE_P(Bond,
                         ListedForcesTest,
                         ::testing::Combine(::testing::ValuesIn(c_InputBonds),
                                            ::testing::ValuesIn(c_coordinatesForTests),
                                            ::testing::ValuesIn(c_pbcForTests),
                                            ::testing::ValuesIn(c_lambdaValuesForTests),
                                            ::testing::ValuesIn(c_flavorsForTests),
                                            ::testing::ValuesIn(getHardwareContextsWithCapability(
                                                    GpuConfigurationCapabilities::Bonded))),
                         sc_testNamer);

INSTANTIATE_TEST_SUITE_P(Angle,
                         ListedForcesTest,
                         ::testing::Combine(::testing::ValuesIn(c_InputAngles),
                                            ::testing::ValuesIn(c_coordinatesForTests),
                                            ::testing::ValuesIn(c_pbcForTests),
                                            ::testing::ValuesIn(c_lambdaValuesForTests),
                                            ::testing::ValuesIn(c_flavorsForTests),
                                            ::testing::ValuesIn(getHardwareContextsWithCapability(
                                                    GpuConfigurationCapabilities::Bonded))),
                         sc_testNamer);

INSTANTIATE_TEST_SUITE_P(Dihedral,
                         ListedForcesTest,
                         ::testing::Combine(::testing::ValuesIn(c_InputDihs),
                                            ::testing::ValuesIn(c_coordinatesForTests),
                                            ::testing::ValuesIn(c_pbcForTests),
                                            ::testing::ValuesIn(c_lambdaValuesForTests),
                                            ::testing::ValuesIn(c_flavorsForTests),
                                            ::testing::ValuesIn(getHardwareContextsWithCapability(
                                                    GpuConfigurationCapabilities::Bonded))),
                         sc_testNamer);

INSTANTIATE_TEST_SUITE_P(Polarize,
                         ListedForcesTest,
                         ::testing::Combine(::testing::ValuesIn(c_InputPols),
                                            ::testing::ValuesIn(c_coordinatesForTests),
                                            ::testing::ValuesIn(c_pbcForTests),
                                            ::testing::ValuesIn(c_lambdaValuesForTests),
                                            ::testing::ValuesIn(c_flavorsForTests),
                                            ::testing::ValuesIn(getHardwareContextsWithCapability(
                                                    GpuConfigurationCapabilities::Bonded))),
                         sc_testNamer);

INSTANTIATE_TEST_SUITE_P(Restraints,
                         ListedForcesTest,
                         ::testing::Combine(::testing::ValuesIn(c_InputRestraints),
                                            ::testing::ValuesIn(c_coordinatesForTests),
                                            ::testing::ValuesIn(c_pbcForTests),
                                            ::testing::ValuesIn(c_lambdaValuesForTests),
                                            ::testing::ValuesIn(c_flavorsForTests),
                                            ::testing::ValuesIn(getHardwareContextsWithCapability(
                                                    GpuConfigurationCapabilities::Bonded))),
                         sc_testNamer);

INSTANTIATE_TEST_SUITE_P(BondZeroLength,
                         ListedForcesTest,
                         ::testing::Combine(::testing::ValuesIn(c_InputBondsZeroLength),
                                            ::testing::ValuesIn(c_coordinatesForTestsZeroBondLength),
                                            ::testing::ValuesIn(c_pbcForTests),
                                            ::testing::ValuesIn(c_lambdaValuesForTests),
                                            ::testing::ValuesIn(c_flavorsForTests),
                                            ::testing::ValuesIn(getHardwareContextsWithCapability(
                                                    GpuConfigurationCapabilities::Bonded))),
                         sc_testNamer);

INSTANTIATE_TEST_SUITE_P(AngleZero,
                         ListedForcesTest,
                         ::testing::Combine(::testing::ValuesIn(c_InputAnglesZeroAngle),
                                            ::testing::ValuesIn(c_coordinatesForTestsZeroAngle),
                                            ::testing::ValuesIn(c_pbcForTests),
                                            ::testing::ValuesIn(c_lambdaValuesForTests),
                                            ::testing::ValuesIn(c_flavorsForTests),
                                            ::testing::ValuesIn(getHardwareContextsWithCapability(
                                                    GpuConfigurationCapabilities::Bonded))),
                         sc_testNamer);

} // namespace

} // namespace test

} // namespace gmx
