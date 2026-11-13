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

#include <gmock/gmock.h>
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
#include "testutils/testmatchers.h"
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

//! Returns the global index of each local atom for the given local atom order
std::vector<int> localToGlobalIndices(LocalAtomOrder localAtomOrder, int numAtoms)
{
    std::vector<int> localToGlobal(numAtoms);
    std::iota(localToGlobal.begin(), localToGlobal.end(), 0);
    if (localAtomOrder == LocalAtomOrder::Reversed)
    {
        std::reverse(localToGlobal.begin(), localToGlobal.end());
    }
    return localToGlobal;
}

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
 *         def forward(self, input_data: Dict[str, torch.Tensor]) -> Dict[str, torch.Tensor]:
 *             positions = input_data["positions"]
 *             atomic_numbers = input_data["atomic-numbers"]
 *             pairs = input_data["pairs"]
 *             shifts = input_data["shifts"]
 *             i = pairs[:, 0].long()
 *             j = pairs[:, 1].long()
 *             d = positions[i] - positions[j] + shifts.to(positions.dtype)
 *             z = atomic_numbers.to(positions.dtype)
 *             return {"energy": (z[i] * z[j] * (d * d).sum(dim=1)).sum()}
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

    //! Returns whether the atom with global index \p globalIndex is part of the NNP region
    static bool isNNPAtom(int globalIndex)
    {
        return std::find(sc_nnpIndices.begin(), sc_nnpIndices.end(), globalIndex) != sc_nnpIndices.end();
    }

    /*! \brief Sets up the local atom order and returns the local coordinates
     *
     * \param[in] localToGlobal  localToGlobal[l] is the global index of the atom at local index l
     */
    std::vector<RVec> setLocalAtomOrder(ArrayRef<const int> localToGlobal)
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

        return xLocal;
    }

    /*! \brief Returns all pairs within \p range in local indices, as nbnxm provides them
     *
     * Uses the convention: distance = x[first] - x[second] + shiftVec[shiftIndex]
     *
     * \param[in] x         Local coordinates
     * \param[in] box       Box for computing the shift vectors
     * \param[in] range     Pair search range
     * \param[in] nnpOnly   When not empty, only return pairs of NNP atoms; holds the local to
     *                      global index mapping
     */
    static std::vector<PairlistEntry> pairsWithinRange(ArrayRef<const RVec> x,
                                                       const matrix         box,
                                                       real                 range,
                                                       ArrayRef<const int>  nnpOnly = {})
    {
        std::array<RVec, c_numShiftVectors> shiftVecs;
        calc_shifts(box, shiftVecs);

        std::vector<PairlistEntry> pairs;
        for (int first = 0; first < gmx::ssize(x); first++)
        {
            for (int second = first + 1; second < gmx::ssize(x); second++)
            {
                if (!nnpOnly.empty() && !(isNNPAtom(nnpOnly[first]) && isNNPAtom(nnpOnly[second])))
                {
                    continue;
                }
                for (int shiftIndex = 0; shiftIndex < c_numShiftVectors; shiftIndex++)
                {
                    if (norm(x[first] - x[second] + shiftVecs[shiftIndex]) <= range)
                    {
                        pairs.push_back({ { first, second }, shiftIndex });
                    }
                }
            }
        }
        return pairs;
    }

    //! \brief Helper struct to hold a pair of energy and forces
    struct EnergyAndForces
    {
        double            energy;
        std::vector<RVec> forces;
    };

    /*! \brief Computes energy and forces of model_pairs.pt for all NNP pairs within the cutoff
     *
     * This is the expected output of the force provider, independent of the pairlist buffer.
     */
    EnergyAndForces referenceEnergyAndForces(ArrayRef<const RVec> x,
                                             ArrayRef<const int>  localToGlobal,
                                             const matrix         box) const
    {
        std::array<RVec, c_numShiftVectors> shiftVecs;
        calc_shifts(box, shiftVecs);

        double            energy = 0;
        std::vector<RVec> forces(x.size(), RVec{ 0, 0, 0 });
        for (const auto& [atomPair, shiftIndex] : pairsWithinRange(x, box, sc_pairCutoff, localToGlobal))
        {
            const auto [first, second] = atomPair;
            const RVec d               = x[first] - x[second] + shiftVecs[shiftIndex];
            const real zz              = params_.atoms_.atom[localToGlobal[first]].atomnumber
                            * params_.atoms_.atom[localToGlobal[second]].atomnumber;
            energy += zz * norm2(d);
            forces[first] -= 2 * zz * d;
            forces[second] += 2 * zz * d;
        }
        return { energy, forces };
    }

    //! Runs the force provider on local coordinates \p x and returns energy and forces
    EnergyAndForces calculateEnergyAndForces(NNPotForceProvider*  nnpotForceProvider,
                                             ArrayRef<const RVec> x,
                                             const matrix         box)
    {
        ForceProviderInput  fInput(x, params_.numAtoms_, {}, {}, 0.0, 0, box, mpiComm_, nullptr);
        std::vector<RVec>   forces(params_.numAtoms_, RVec{ 0, 0, 0 });
        ForceWithVirial     forceWithVirial(forces, true);
        gmx_enerdata_t      enerd(1, nullptr);
        ForceProviderOutput forceProviderOutput(&forceWithVirial, &enerd);
        nnpotForceProvider->calculateForces(fInput, &forceProviderOutput);

        return { enerd.term[InteractionFunction::NeuralNetworkPotentialEnergy], forces };
    }

    //! Checks that the force provider output matches the brute-force reference
    void compareEnergyAndForces(const EnergyAndForces& first, const EnergyAndForces& second)
    {
        EXPECT_REAL_EQ_TOL(
                first.energy, second.energy, relativeToleranceAsFloatingPoint(second.energy, 5e-5));
        EXPECT_THAT(first.forces,
                    ::testing::Pointwise(RVecEq(relativeToleranceAsFloatingPoint(100.0, 5e-5)),
                                         second.forces));
    }

    /*! \brief Calculates NNP forces with the given local atom order and checks them
     *
     * \param[in] localToGlobal  localToGlobal[l] is the global index of the atom at local index l
     */
    void testCalculateForces(ArrayRef<const int> localToGlobal)
    {
        const int               numAtoms = params_.numAtoms_;
        const std::vector<RVec> xLocal   = setLocalAtomOrder(localToGlobal);

        // Excluded pairlist over all atoms in local indices
        const std::vector<PairlistEntry> excludedPairs = pairsWithinRange(xLocal, box_, sc_pairCutoff);
        const auto nnpPairs = pairsWithinRange(xLocal, box_, sc_pairCutoff, localToGlobal);
        ASSERT_TRUE(std::any_of(nnpPairs.begin(),
                                nnpPairs.end(),
                                [](const PairlistEntry& p)
                                { return p.second != c_centralShiftIndex; }))
                << "Test setup should produce NNP pairs across the boundary";

        NNPotForceProvider nnpotForceProvider(params_, logger_, mpiComm_);
        nnpotForceProvider.gatherAtomNumbersIndices(
                MDModulesAtomsRedistributedSignal{ box_, xLocal, {}, {}, std::nullopt });
        nnpotForceProvider.setPairlist(MDModulesPairlistConstructedSignal{ {}, excludedPairs, {} });

        const auto result = calculateEnergyAndForces(&nnpotForceProvider, xLocal, box_);

        // Check forces in global order, so all local orders share the same reference
        std::vector<RVec> forcesGlobal(numAtoms);
        for (int i = 0; i < numAtoms; i++)
        {
            forcesGlobal[localToGlobal[i]] = result.forces[i];
        }

        checker().setDefaultTolerance(relativeToleranceAsFloatingPoint(10.0, 5e-5));
        checker().checkReal(result.energy, "Energy");
        checker().setDefaultTolerance(relativeToleranceAsFloatingPoint(100.0, 5e-5));
        checker().checkSequence(forcesGlobal.begin(), forcesGlobal.end(), "Forces");
    }

protected:
    //! NNP region: first two of four waters
    static inline const std::vector<Index> sc_nnpIndices = { 0, 1, 2, 3, 4, 5 };
    //! Cutoff for NNP pairs
    static constexpr real sc_pairCutoff = 0.5;
    //! Pairlist buffer, as computed from the pairlist lifetime in mdrun
    static constexpr real sc_pairlistBuffer = 0.1;

    NNPotParameters     params_;
    MDLogger            logger_;
    LocalAtomSetManager atomSetManager_;
    const MpiComm       mpiComm_ = MpiComm(MpiComm::SingleRank{});
    gmx_mtop_t          mtop_;
    PbcType             pbcType_;
    matrix              box_;
    std::vector<RVec>   xGlobal_;
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

    testCalculateForces(localToGlobalIndices(localAtomOrder, params_.numAtoms_));
}

/*! \brief Test that pairs from the buffer region are pruned at every step
 *
 * The pairlist is constructed once with a buffer, as mdrun does. Then the system is rescaled
 * without constructing a new pairlist, so that pairs from the buffer region move inside the
 * cutoff and pairs inside the cutoff move out, and the shift vectors change.
 * This check makes sure that the pairs within the buffer region do not contribute to
 * the energy and force calculation by the model, by checking against a reference implementation
 * (in referenceEnergyAndForces()). That is, the model should always see exactly the NNP pairs within the cutoff.
 */
TEST_P(NNPotPairlistInputTest, PrunesBufferedPairlistEveryStep)
{
    auto [localAtomOrder, _] = GetParam();

    SCOPED_TRACE(formatString("Testing on %s with %s",
                              hardwareContext()->description().c_str(),
                              sc_localAtomOrderNames[localAtomOrder]));

    activateHardware();
    auto envGuard = GmxEnvGuard("GMX_NN_DEVICE", isGpuTest() ? "gpu" : "cpu");

    const std::vector<int>  localToGlobal = localToGlobalIndices(localAtomOrder, params_.numAtoms_);
    const std::vector<RVec> xLocal        = setLocalAtomOrder(localToGlobal);

    const real pairlistRange = sc_pairCutoff + sc_pairlistBuffer;
    const auto numNNPPairs   = [&](ArrayRef<const RVec> x, const matrix box, real range)
    { return pairsWithinRange(x, box, range, localToGlobal).size(); };
    ASSERT_GT(numNNPPairs(xLocal, box_, pairlistRange), numNNPPairs(xLocal, box_, sc_pairCutoff))
            << "Test setup should produce NNP pairs in the buffer region";

    NNPotForceProvider nnpotForceProvider(params_, logger_, mpiComm_);
    nnpotForceProvider.gatherAtomNumbersIndices(
            MDModulesAtomsRedistributedSignal{ box_, xLocal, {}, {}, std::nullopt });
    nnpotForceProvider.setPairlist(MDModulesPairlistConstructedSignal{
            {}, pairsWithinRange(xLocal, box_, pairlistRange), {} });

    {
        SCOPED_TRACE("Step with pairlist construction");
        const auto energyAndForcesFromModel = calculateEnergyAndForces(&nnpotForceProvider, xLocal, box_);
        const auto energyAndForcesFromReference = referenceEnergyAndForces(xLocal, localToGlobal, box_);
        compareEnergyAndForces(energyAndForcesFromModel, energyAndForcesFromReference);
    }

    // Scale factors are chosen such that all pairs within the scaled cutoff are still within
    // the pairlist range, i.e. the buffer is sufficient
    for (const real scale : { real(0.9), real(1.1) })
    {
        SCOPED_TRACE(formatString("Step without pairlist construction, scaled by %g", scale));

        matrix box;
        for (int d = 0; d < DIM; d++)
        {
            for (int e = 0; e < DIM; e++)
            {
                box[d][e] = scale * box_[d][e];
            }
        }
        std::vector<RVec> x(xLocal);
        for (auto& xi : x)
        {
            xi *= scale;
        }
        ASSERT_NE(numNNPPairs(x, box, sc_pairCutoff), numNNPPairs(xLocal, box_, sc_pairCutoff))
                << "Scaling should change the set of NNP pairs within the cutoff";
        ASSERT_LE(sc_pairCutoff / scale, pairlistRange)
                << "Test setup should have all pairs within the cutoff in the pairlist";

        const auto energyAndForcesFromModel = calculateEnergyAndForces(&nnpotForceProvider, x, box);
        const auto energyAndForcesFromReference = referenceEnergyAndForces(x, localToGlobal, box);
        compareEnergyAndForces(energyAndForcesFromModel, energyAndForcesFromReference);
    }
}

//! Test that a newly constructed pairlist replaces the previous one
TEST_P(NNPotPairlistInputTest, NewPairlistReplacesPrevious)
{
    auto [localAtomOrder, _] = GetParam();

    SCOPED_TRACE(formatString("Testing on %s with %s",
                              hardwareContext()->description().c_str(),
                              sc_localAtomOrderNames[localAtomOrder]));

    activateHardware();
    auto envGuard = GmxEnvGuard("GMX_NN_DEVICE", isGpuTest() ? "gpu" : "cpu");

    const std::vector<int>  localToGlobal = localToGlobalIndices(localAtomOrder, params_.numAtoms_);
    const std::vector<RVec> xLocal        = setLocalAtomOrder(localToGlobal);
    const std::vector<PairlistEntry> excludedPairs =
            pairsWithinRange(xLocal, box_, sc_pairCutoff + sc_pairlistBuffer);

    NNPotForceProvider nnpotForceProvider(params_, logger_, mpiComm_);
    nnpotForceProvider.gatherAtomNumbersIndices(
            MDModulesAtomsRedistributedSignal{ box_, xLocal, {}, {}, std::nullopt });
    // Repeatedly set a (new) pairlist and check against reference, to make sure that the
    // pairlist is reset upon receiving a new signal, instead of appended to
    for (int construction = 0; construction < 3; construction++)
    {
        SCOPED_TRACE(formatString("Pairlist construction %d", construction));
        nnpotForceProvider.setPairlist(MDModulesPairlistConstructedSignal{ {}, excludedPairs, {} });

        const auto energyAndForcesFromModel = calculateEnergyAndForces(&nnpotForceProvider, xLocal, box_);
        const auto energyAndForcesFromReference = referenceEnergyAndForces(xLocal, localToGlobal, box_);
        compareEnergyAndForces(energyAndForcesFromModel, energyAndForcesFromReference);
    }
}

INSTANTIATE_TEST_SUITE_P(
        AllHardware,
        NNPotPairlistInputTest,
        ::testing::Combine(::testing::Values(LocalAtomOrder::Global, LocalAtomOrder::Reversed),
                           ::testing::ValuesIn(getHardwareContextsWithCapability(sc_nnpotSupportsGpu))),
        sc_pairlistTestNamer);

/*! \brief Test fixture for the virial contribution of NNP models
 *
 * model_virial.pt computes a Gaussian pair repulsion with minimum image convention,
 * E = 10 sum_{i<j} exp(-|x_ij|^2 / 0.15^2), and returns forces and the virial
 * -1/2 sum_i x_i (x) f_i = 1/2 (dE/dstrain)^T, obtained by automatic differentiation
 * with respect to a homogeneous strain of positions and box. It was generated with:
 *
 *     class VirialModel(torch.nn.Module):
 *         def forward(self, input_data: Dict[str, torch.Tensor]) -> Dict[str, torch.Tensor]:
 *             positions = input_data["positions"].to(torch.float64)
 *             box = input_data["box"].to(torch.float64)
 *             strain = torch.zeros(3, 3, dtype=torch.float64,
 *                                  device=positions.device).requires_grad_()
 *             x = positions + positions @ strain.t()
 *             b = box + box @ strain.t()
 *             s = x @ torch.inverse(b)
 *             ds = s.unsqueeze(0) - s.unsqueeze(1)
 *             dx = (ds - torch.round(ds)) @ b
 *             r2 = (dx * dx).sum(dim=-1)
 *             mask = torch.triu(torch.ones_like(r2, dtype=torch.bool), diagonal=1)
 *             energy = 10.0 * torch.exp(-r2[mask] / 0.15**2).sum()
 *             gradX, gradStrain = torch.autograd.grad([energy], [positions, strain])
 *             assert gradX is not None and gradStrain is not None
 *             dtype = input_data["positions"].dtype
 *             return {"energy": energy.to(dtype),
 *                     "forces": (-gradX).to(dtype),
 *                     "virial": (0.5 * gradStrain.t()).to(dtype)}
 *
 *     torch.jit.script(VirialModel()).save("model_virial.pt")
 */
class NNPotVirialTest : public ::testing::Test
{
public:
    NNPotVirialTest()
    {
        TprAndFileManager tprManager("4water");
        tprManager.setInputConf("4water.gro");
        bool  fullTopology;
        rvec* coords = nullptr;
        readConfAndTopology(tprManager.tprName(), &fullTopology, &mtop_, &pbcType_, &coords, nullptr, box_);
        positions_.assign(reinterpret_cast<RVec*>(coords), reinterpret_cast<RVec*>(coords) + mtop_.natoms);
        sfree(coords);

        std::vector<Index> nnpIndices(mtop_.natoms);
        std::iota(nnpIndices.begin(), nnpIndices.end(), 0);
        params_.active_     = true;
        params_.atoms_      = gmx_mtop_global_atoms(mtop_);
        params_.numAtoms_   = mtop_.natoms;
        params_.pbcType_    = std::make_unique<PbcType>(pbcType_);
        LocalAtomSet nnpSet = atomSetManager_.add(ArrayRef<const Index>(nnpIndices));
        params_.nnpAtoms_   = std::make_unique<LocalAtomSet>(nnpSet);
    }

    ~NNPotVirialTest() override { done_atom(&params_.atoms_); }

    //! Energy and virial returned by the force provider
    struct EnergyAndVirial
    {
        real   energy;
        matrix virial;
    };

    //! Evaluates the model in \p modelFile on the given configuration
    EnergyAndVirial evaluate(const std::string& modelFile, ArrayRef<const RVec> x, const matrix box)
    {
        params_.modelFileName_ = TestFileManager::getInputFilePath(modelFile).string();
        NNPotForceProvider nnpotForceProvider(params_, logger_, mpiComm_);
        nnpotForceProvider.gatherAtomNumbersIndices(
                MDModulesAtomsRedistributedSignal{ box, x, {}, {}, std::nullopt });

        ForceProviderInput  fInput(x, params_.numAtoms_, {}, {}, 0.0, 0, box, mpiComm_, nullptr);
        std::vector<RVec>   forces(params_.numAtoms_, RVec{ 0, 0, 0 });
        ForceWithVirial     forceWithVirial(forces, true);
        gmx_enerdata_t      enerd(1, nullptr);
        ForceProviderOutput forceProviderOutput(&forceWithVirial, &enerd);
        nnpotForceProvider.calculateForces(fInput, &forceProviderOutput);

        EnergyAndVirial result;
        result.energy = enerd.term[InteractionFunction::NeuralNetworkPotentialEnergy];
        copy_mat(forceWithVirial.getVirial(), result.virial);
        return result;
    }

protected:
    NNPotParameters     params_;
    MDLogger            logger_;
    LocalAtomSetManager atomSetManager_;
    const MpiComm       mpiComm_ = MpiComm(MpiComm::SingleRank{});
    gmx_mtop_t          mtop_;
    PbcType             pbcType_;
    matrix              box_;
    std::vector<RVec>   positions_;
};

/*! \brief Checks the virial from the model against the energy change under box deformation
 *
 * Under a homogeneous strain x' = (1 + strain) x of positions and box, the virial
 * Xi = -1/2 sum_i x_i (x) f_i satisfies dE/dstrain_ij = 2 Xi_ji. We use central
 * differences for each strain component that keeps the box lower triangular.
 */
TEST_F(NNPotVirialTest, VirialMatchesBoxDeformation)
{
    auto envGuard = GmxEnvGuard("GMX_NN_DEVICE", "cpu");

    matrix virialFromModel;
    copy_mat(evaluate("model_virial.pt", positions_, box_).virial, virialFromModel);

    constexpr real c_strain = 5e-3;
    for (int i = 0; i < DIM; i++)
    {
        for (int j = i; j < DIM; j++)
        {
            real energyPosStrain, energyNegStrain;
            for (int sign : { 1, -1 })
            {
                std::vector<RVec> strainedPositions = positions_;
                matrix            strainedBox;
                copy_mat(box_, strainedBox);
                for (auto& pos : strainedPositions)
                {
                    pos[i] += sign * c_strain * pos[j];
                }
                for (int d = 0; d < DIM; d++)
                {
                    strainedBox[d][i] += sign * c_strain * box_[d][j];
                }
                if (sign > 0)
                {
                    energyPosStrain = evaluate("model_virial.pt", strainedPositions, strainedBox).energy;
                }
                else
                {
                    energyNegStrain = evaluate("model_virial.pt", strainedPositions, strainedBox).energy;
                }
            }
            const real virialFromEnergy = (energyPosStrain - energyNegStrain) / (4 * c_strain);

            SCOPED_TRACE(formatString("Strain component %d %d", i, j));
            EXPECT_REAL_EQ_TOL(virialFromEnergy, virialFromModel[j][i], absoluteTolerance(3e-3));
        }
    }
}

//! Checks that no virial is added when the model does not return one
TEST_F(NNPotVirialTest, NoVirialWithoutModelOutput)
{
    auto envGuard = GmxEnvGuard("GMX_NN_DEVICE", "cpu");

    const EnergyAndVirial result = evaluate("model.pt", positions_, box_);
    for (int i = 0; i < DIM; i++)
    {
        for (int j = 0; j < DIM; j++)
        {
            EXPECT_EQ(result.virial[i][j], 0.0) << "virial component " << i << " " << j;
        }
    }
}

#endif // GMX_TORCH

} // namespace test

} // namespace gmx
