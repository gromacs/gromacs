/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2024- The GROMACS Authors
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
 * Tests for functionality of the NNPot ForceProvider
 *
 * \author Lukas Müllender <lukas.muellender@gmail.com>
 * \ingroup module_applied_forces
 */
#include "gmxpre.h"

#include "gromacs/applied_forces/nnpot/nnpotforceprovider.h"

#include "config.h"

#include <algorithm>
#include <numeric>
#include <string>
#include <tuple>
#include <vector>

#include <gtest/gtest.h>

#include "gromacs/applied_forces/nnpot/nnpotoptions.h"
#include "gromacs/domdec/ga2la.h"
#include "gromacs/domdec/localatomset.h"
#include "gromacs/domdec/localatomsetmanager.h"
#include "gromacs/fileio/confio.h"
#include "gromacs/hardware/device_information.h"
#include "gromacs/mdrunutility/mdmodulesnotifiers.h"
#include "gromacs/mdtypes/enerdata.h"
#include "gromacs/mdtypes/forceoutput.h"
#include "gromacs/pbcutil/ishift.h"
#include "gromacs/pbcutil/pbc.h"
#include "gromacs/topology/mtop_lookup.h"
#include "gromacs/topology/mtop_util.h"
#include "gromacs/topology/topology.h"
#include "gromacs/utility/enumerationhelpers.h"
#include "gromacs/utility/logger.h"
#include "gromacs/utility/message_string_collector.h"
#include "gromacs/utility/mpicomm.h"
#include "gromacs/utility/smalloc.h"
#include "gromacs/utility/stringutil.h"

#include "testutils/hardware_test_fixture.h"
#include "testutils/naming.h"
#include "testutils/refdata.h"
#include "testutils/setenv.h"
#include "testutils/test_hardware_environment.h"
#include "testutils/testasserts.h"
#include "testutils/testfilemanager.h"
#include "testutils/tprfilegenerator.h"

namespace gmx
{

namespace test
{

class NNPotForceProviderTest : public ::testing::Test
{
public:
    /*! \brief Generates tpr file from *.top and *.gro existing in the simulation database directory
     * and loads gmx_mtop_t from it
     */
    void readDataFromFile(const std::string& fileName, const std::string& mdpContent)
    {
        TprAndFileManager tprManager(fileName, mdpContent);
        tprManager.setInputConf(fileName + ".gro");
        bool fullTopology;
        readConfAndTopology(tprManager.tprName(), &fullTopology, &mtop_, &pbcType_, &coords, nullptr, box_);
        x_ = gmx::constArrayRefFromArray(reinterpret_cast<gmx::RVec*>(coords), mtop_.natoms);
    }

    void setDefaultParameters()
    {
        params_.active_        = true;
        params_.modelFileName_ = TestFileManager::getInputFilePath("model.pt").string();
        params_.numAtoms_      = 6;
    }

    void setParametersTopology()
    {
        params_.active_        = true;
        params_.modelFileName_ = TestFileManager::getInputFilePath("model.pt").string();
        params_.atoms_         = gmx_mtop_global_atoms(mtop_);
        params_.numAtoms_      = params_.atoms_.nr;
        params_.modelInput_    = { "atom-positions", "atom-numbers", "box", "pbc" };

        std::vector<gmx::Index> nnpIndices = { 0, 1, 2, 3, 4, 5 };
        LocalAtomSet            set = atomSetManager_.add(ArrayRef<const gmx::Index>(nnpIndices));
        params_.nnpAtoms_           = std::make_unique<LocalAtomSet>(set);
        params_.pbcType_            = std::make_unique<PbcType>(pbcType_);
    }

    void testCalculateForces(TestReferenceData& testData)
    {
        TestReferenceChecker checker(testData.rootChecker());

        std::unique_ptr<NNPotForceProvider> nnpotForceProvider;
        ASSERT_NO_THROW(nnpotForceProvider =
                                std::make_unique<NNPotForceProvider>(params_, logger_, mpiComm_));
        MDModulesAtomsRedistributedSignal signal{ box_, x_, {}, {}, std::nullopt };
        ASSERT_NO_THROW(nnpotForceProvider->gatherAtomNumbersIndices(signal));

        // Prepare input for force provider
        ForceProviderInput fInput(x_, params_.numAtoms_, {}, {}, 0.0, 0, box_, mpiComm_, dd_);

        // Prepare output for force provider
        std::vector<RVec>   forces(params_.numAtoms_, RVec{ 0, 0, 0 });
        ForceWithVirial     forceWithVirial(forces, true);
        gmx_enerdata_t      enerdDummy(1, nullptr);
        ForceProviderOutput forceProviderOutput(&forceWithVirial, &enerdDummy);

        ASSERT_NO_THROW(nnpotForceProvider->calculateForces(fInput, &forceProviderOutput));

        checker.setDefaultTolerance(relativeToleranceAsFloatingPoint(100000.0, 5e-5));
        checker.checkReal(enerdDummy.term[InteractionFunction::NeuralNetworkPotentialEnergy],
                          "Energy");
        checker.setDefaultTolerance(relativeToleranceAsFloatingPoint(100.0, 5e-5));
        checker.checkSequence(forces.begin(), forces.end(), "Forces");
    }

protected:
    NNPotParameters               params_;
    MDLogger                      logger_;
    LocalAtomSetManager           atomSetManager_;
    std::unique_ptr<LocalAtomSet> inpAtomSet_;
    const MpiComm                 mpiComm_ = MpiComm(MpiComm::SingleRank{});
    const gmx_domdec_t*           dd_      = nullptr;

    rvec*                coords = nullptr;
    ArrayRef<const RVec> x_;
    gmx_mtop_t           mtop_;
    PbcType              pbcType_;
    matrix               box_;
};

TEST_F(NNPotForceProviderTest, CanConstruct)
{
    setDefaultParameters();

    // GMX_TORCH is defined by set_source_files_properties() in CMakeLists.txt
    if (GMX_TORCH)
    {
        {
            SCOPED_TRACE("Check construction on CPU");
            EXPECT_NO_THROW(NNPotForceProvider nnpotForceProvider(params_, logger_, mpiComm_));
        }

        {
            SCOPED_TRACE("Check construction with invalid model file name");
            params_.modelFileName_ = "model";
            EXPECT_THROW_GMX(NNPotForceProvider nnpotForceProvider(params_, logger_, mpiComm_), FileIOError);
            params_.modelFileName_ = TestFileManager::getInputFilePath("model.pt").string();
        }

        for (const auto& device : getTestHardwareEnvironment()->getTestDeviceList())
        {
            const DeviceInformation deviceInfo = device->deviceInfo();
            if (deviceInfo.deviceVendor == DeviceVendor::Nvidia
                || (GMX_TORCH_XPU && deviceInfo.deviceVendor == DeviceVendor::Intel))
            {
                SCOPED_TRACE("Check construction on default GPU");
                auto envGuard = GmxEnvGuard("GMX_NN_DEVICE", "gpu");
                EXPECT_NO_THROW(NNPotForceProvider nnpotForceProvider(params_, logger_, mpiComm_));
                break; // Only test one GPU until we have a better way to hande device selection
            }
        }
    }
    else
    {
        EXPECT_ANY_THROW(NNPotForceProvider nnpotForceProvider(params_, logger_, mpiComm_));
    }
}

#if GMX_TORCH
/*! Test if the NNPotForceProvider can calculate forces.
 * This implicitly tests the NNPotModel as well as the input gathering functions.
 */
// Temporary skip while we iron out version incompatibilities in CI
TEST_F(NNPotForceProviderTest, CanCalculateForces)
{
    TestReferenceData data;

    readDataFromFile("spc2", "");
    setParametersTopology();

    {
        SCOPED_TRACE("Check calculate forces on CPU");
        testCalculateForces(data);
    }
    if (!getTestHardwareEnvironment()->hasCompatibleDevices())
    {
        return;
    }
    for (const auto& device : getTestHardwareEnvironment()->getTestDeviceList())
    {
        const DeviceInformation deviceInfo = device->deviceInfo();
        if (deviceInfo.deviceVendor == DeviceVendor::Nvidia
            || (GMX_TORCH_XPU && deviceInfo.deviceVendor == DeviceVendor::Intel))
        {
            SCOPED_TRACE("Check calculate forces on default GPU");
            auto envGuard = GmxEnvGuard("GMX_NN_DEVICE", "gpu");
            testCalculateForces(data);
            break; // Only test one GPU until we have a better way to hande device selection
        }
    }

    done_atom(&params_.atoms_);
}

//! Order of local atom indices relative to global atom indices
enum class LocalAtomOrder : int
{
    Global,   //!< Local order equals global order
    Reversed, //!< Local order differs from global order
    Count
};

//! Names for local atom orders
static const EnumerationArray<LocalAtomOrder, const char*> sc_localAtomOrderNames = {
    { "GlobalOrder", "ReversedOrder" }
};

//! Whether NNP models can run on a GPU in this configuration, cf. getMdrunActiveDevice()
static constexpr bool sc_nnpotSupportsGpu = GMX_GPU_CUDA || GMX_GPU_HIP || GMX_TORCH_XPU;

/*! \brief Hardware test helper for models with pairlist input
 *
 * All test cases compute the same energy and forces, so there is no input
 * configuration, and the local atom order is an execution mode.
 */
using NNPotPairlistTestHelper = HardwareAndExecutionTestHelper<std::tuple<>, std::tuple<LocalAtomOrder>>;

//! Helper object to name tests using all parameters
static const NameOfTestFromTuple<NNPotPairlistTestHelper::DynamicParameters> sc_pairlistTestNamer =
        NNPotPairlistTestHelper::testNamer(
                std::make_tuple(),
                std::make_tuple([](LocalAtomOrder order) { return sc_localAtomOrderNames[order]; }));

/*! \brief Test fixture for models with pairlist input (atom pairs and pair shifts)
 *
 * Runs on CPU and all compatible GPUs, with local atom order equal to the global
 * order and with reordered local atoms, as with domain decomposition (which is
 * also used on a single rank). All cases share the same reference data.
 * The NNP region straddles the periodic boundary, so we also test non-central
 * shift vectors.
 *
 * model_pairs.pt computes E = sum_pairs Z_i Z_j |x_i - x_j + shift|^2. It was generated with:
 *
 *     class PairModel(torch.nn.Module):
 *         def forward(self, positions, atomic_numbers, pairs, shifts):
 *             i = pairs[:, 0].long()
 *             j = pairs[:, 1].long()
 *             d = positions[i] - positions[j] + shifts.to(positions.dtype)
 *             z = atomic_numbers.to(positions.dtype)
 *             return (z[i] * z[j] * (d * d).sum(dim=1)).sum()
 *
 *     torch.jit.script(PairModel()).save("model_pairs.pt")
 */
class NNPotPairlistInputTest : public HardwareTestFixture<NNPotPairlistTestHelper>
{
public:
    NNPotPairlistInputTest() : HardwareTestFixture(std::make_tuple())
    {
        TprAndFileManager tprManager("4water");
        tprManager.setInputConf("4water.gro");
        bool  fullTopology;
        rvec* coords = nullptr;
        readConfAndTopology(tprManager.tprName(), &fullTopology, &mtop_, &pbcType_, &coords, nullptr, box_);
        xGlobal_.assign(reinterpret_cast<RVec*>(coords), reinterpret_cast<RVec*>(coords) + mtop_.natoms);
        sfree(coords);

        params_.active_        = true;
        params_.modelFileName_ = TestFileManager::getInputFilePath("model_pairs.pt").string();
        params_.atoms_         = gmx_mtop_global_atoms(mtop_);
        params_.numAtoms_      = mtop_.natoms;
        params_.modelInput_    = { "atom-positions", "atom-numbers", "atom-pairs", "pair-shifts" };
        params_.pairCutoff_    = sc_pairCutoff;
        params_.pbcType_       = std::make_unique<PbcType>(pbcType_);
        LocalAtomSet nnpSet    = atomSetManager_.add(ArrayRef<const Index>(sc_nnpIndices));
        params_.nnpAtoms_      = std::make_unique<LocalAtomSet>(nnpSet);

        // Translate the system so the NNP waters straddle the periodic boundary,
        // which causes non-central shift vectors in the pairlist
        for (auto& x : xGlobal_)
        {
            x -= RVec{ 0.8, 1.25, 1.2 };
        }
        put_atoms_in_box(pbcType_, box_, xGlobal_);

        calc_shifts(box_, shiftVecs_);
    }

    ~NNPotPairlistInputTest() override { done_atom(&params_.atoms_); }

    void addCustomSkipReasons(MessageStringCollector& skipReasons) override
    {
        // With the XPU backend, the model always runs on the first XPU, since SYCL has no
        // current device that activateHardware() could set, cf. getMdrunActiveDevice()
        skipReasons.appendIf(
                GMX_TORCH_XPU && isGpuTest() && hardwareContext() != &getTestHardwareContexts()[1],
                "With the LibTorch XPU backend, NNP models only run on the first GPU");
    }

    /*! \brief Calculates NNP forces with the given local atom order and checks them
     *
     * \param[in] localToGlobal  localToGlobal[l] is the global index of the atom at local index l
     */
    void testCalculateForces(ArrayRef<const int> localToGlobal)
    {
        const int numAtoms = params_.numAtoms_;

        gmx_ga2la_t       ga2la(numAtoms, numAtoms);
        std::vector<RVec> xLocal(numAtoms);
        for (int i = 0; i < numAtoms; i++)
        {
            ga2la.insert(localToGlobal[i], { i, 0 });
            xLocal[i] = xGlobal_[localToGlobal[i]];
        }
        atomSetManager_.setIndicesInDomainDecomposition(ga2la);

        // Excluded pairlist over all atoms in local indices, as nbnxm provides it,
        // using the convention: distance = x[first] - x[second] + shiftVec[shiftIndex]
        auto isNNPAtom = [](int globalIndex)
        {
            return std::find(sc_nnpIndices.begin(), sc_nnpIndices.end(), globalIndex)
                   != sc_nnpIndices.end();
        };
        std::vector<PairlistEntry> excludedPairs;
        int                        numPeriodicNNPPairs = 0;
        for (int first = 0; first < numAtoms; first++)
        {
            for (int second = first + 1; second < numAtoms; second++)
            {
                for (int shiftIndex = 0; shiftIndex < c_numShiftVectors; shiftIndex++)
                {
                    if (norm(xLocal[first] - xLocal[second] + shiftVecs_[shiftIndex]) < sc_pairCutoff)
                    {
                        excludedPairs.push_back({ { first, second }, shiftIndex });
                        numPeriodicNNPPairs +=
                                (isNNPAtom(localToGlobal[first]) && isNNPAtom(localToGlobal[second])
                                 && shiftIndex != c_centralShiftIndex);
                    }
                }
            }
        }
        ASSERT_GT(numPeriodicNNPPairs, 0)
                << "Test setup should produce NNP pairs across the boundary";

        NNPotForceProvider nnpotForceProvider(params_, logger_, mpiComm_);
        nnpotForceProvider.gatherAtomNumbersIndices(
                MDModulesAtomsRedistributedSignal{ box_, xLocal, {}, {}, std::nullopt });
        nnpotForceProvider.setPairlist(MDModulesPairlistConstructedSignal{ {}, excludedPairs, {} });

        ForceProviderInput  fInput(xLocal, numAtoms, {}, {}, 0.0, 0, box_, mpiComm_, nullptr);
        std::vector<RVec>   forces(numAtoms, RVec{ 0, 0, 0 });
        ForceWithVirial     forceWithVirial(forces, true);
        gmx_enerdata_t      enerd(1, nullptr);
        ForceProviderOutput forceProviderOutput(&forceWithVirial, &enerd);
        nnpotForceProvider.calculateForces(fInput, &forceProviderOutput);

        // Check forces in global order, so all local orders share the same reference
        std::vector<RVec> forcesGlobal(numAtoms);
        for (int i = 0; i < numAtoms; i++)
        {
            forcesGlobal[localToGlobal[i]] = forces[i];
        }

        checker().setDefaultTolerance(relativeToleranceAsFloatingPoint(10.0, 5e-5));
        checker().checkReal(enerd.term[InteractionFunction::NeuralNetworkPotentialEnergy], "Energy");
        checker().setDefaultTolerance(relativeToleranceAsFloatingPoint(100.0, 5e-5));
        checker().checkSequence(forcesGlobal.begin(), forcesGlobal.end(), "Forces");
    }

protected:
    //! NNP region: first two of four waters
    static inline const std::vector<Index> sc_nnpIndices = { 0, 1, 2, 3, 4, 5 };
    //! Cutoff for NNP pairs
    static constexpr real sc_pairCutoff = 0.5;

    NNPotParameters                     params_;
    MDLogger                            logger_;
    LocalAtomSetManager                 atomSetManager_;
    const MpiComm                       mpiComm_ = MpiComm(MpiComm::SingleRank{});
    gmx_mtop_t                          mtop_;
    PbcType                             pbcType_;
    matrix                              box_;
    std::vector<RVec>                   xGlobal_;
    std::array<RVec, c_numShiftVectors> shiftVecs_;
};

TEST_P(NNPotPairlistInputTest, CanCalculateForces)
{
    auto [localAtomOrder, _] = GetParam();

    SCOPED_TRACE(formatString("Testing on %s with %s",
                              hardwareContext()->description().c_str(),
                              sc_localAtomOrderNames[localAtomOrder]));

    // The model runs on the active GPU device if GMX_NN_DEVICE=gpu
    activateHardware();
    auto envGuard = GmxEnvGuard("GMX_NN_DEVICE", isGpuTest() ? "gpu" : "cpu");

    std::vector<int> localToGlobal(params_.numAtoms_);
    std::iota(localToGlobal.begin(), localToGlobal.end(), 0);
    if (localAtomOrder == LocalAtomOrder::Reversed)
    {
        std::reverse(localToGlobal.begin(), localToGlobal.end());
    }

    testCalculateForces(localToGlobal);
}

INSTANTIATE_TEST_SUITE_P(
        AllHardware,
        NNPotPairlistInputTest,
        ::testing::Combine(::testing::Values(LocalAtomOrder::Global, LocalAtomOrder::Reversed),
                           ::testing::ValuesIn(getHardwareContextsWithCapability(sc_nnpotSupportsGpu))),
        sc_pairlistTestNamer);

#endif // GMX_TORCH

} // namespace test

} // namespace gmx
