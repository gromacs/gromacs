/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2026- The GROMACS Authors
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
 * \brief Hardware-test fixture and utilities for CPU and GPU
 * parameterized tests
 *
 * Provides hardware abstraction for parameterized tests that run on
 * both CPU and all detected compatible GPUs. Enables hardware
 * selection via GoogleTest parameterization with capability-based
 * filtering and shared reference data between CPU and GPU runs.
 *
 * \author Mark Abraham
 * \ingroup module_testutils
 *
 * \section overview Overview
 *
 * This fixture enables writing parameterized tests that run on both
 * CPU and GPU hardware with minimal boilerplate. Tests have
 * descriptive names, automatically run on all available hardware,
 * share reference data between CPU/GPU runs, and support
 * capability-based filtering.
 *
 * \section three_layers The Three-Layer Parameter Model
 *
 * Tests are parameterized across three conceptually distinct layers:
 *
 * 1. **Input Configuration** (required): What to compute
 *    - Examples: system size, timestep, PBC type, physical parameters
 *    - Affects BOTH test names AND reference data names
 *    - Use when: its value changes the expected result values
 *
 * 2. **Execution Modes** (optional): How to compute results
 *    - Examples: forces-only vs all, SIMD vs not, update velocities or not
 *    - Affects test names but NOT reference data names or contents
 *    - Use when: different implementations compute the same values or subsets thereof
 *
 * 3. **Hardware Context** (automatic): CPU vs GPU execution
 *    - Distinct named GoogleTest on each piece of compatible hardware detected
 *    - Added via testing::Combine with getHardwareContextsWithCapability()
 *    - Skipped from reference data names (CPU and GPU share same expected results)
 *
 * \subsection decision_tree When to Use Execution Modes vs Input Configuration
 *
 * Ask: "If I change this parameter, do the **physics results** change?"
 * - YES → Input configuration (affects reference data)
 * - NO → Execution mode (same result values, different computation path)
 *
 * Examples:
 * - Number of atoms: Input (different physics)
 * - SHAKE vs LINCS: Execution mode (same constraint satisfaction, different algorithm)
 * - Compute virial: Execution mode (computes additional outputs, doesn't change forces)
 * - Update velocities: Execution mode (doesn't change constraint forces)
 * - Lambda for FEP: Input (changes force values)
 *
 * \section complete_example Complete Example with Execution Modes
 *
 * This synthetic example demonstrates all features including input configuration, execution modes,
 * custom skip logic, and hardware-capability filtering:
 *
 * \code{.cpp}
 * // Step 1: Define your test's input configuration (what physics to test - affects reference data)
 * struct BondedConfig {
 *     int numAtoms;
 *     PbcType pbcType;
 * };
 *
 * // Step 2: Define execution modes (how to compute - affects test names only)
 * enum class BondedKernelFlavor {
 *     ForcesOnly,                // Compute forces only
 *     ForcesAndVirialAndEnergy   // Also compute virial and energy (same forces, more outputs)
 * };
 *
 * struct BondedExecutionMode {
 *     BondedKernelFlavor flavor;
 *     bool useSIMD;  // Use SIMD vectorization (doesn't change results)
 * };
 *
 * // Step 3: Create test configurations
 * const std::vector<BondedConfig> sc_bondedConfigs = {
 *     {10, PbcType::Xyz},
 *     {100, PbcType::No}
 * };
 *
 * const std::vector<BondedExecutionMode> sc_bondedModes = {
 *     {BondedKernelFlavor::ForcesOnly, false},
 *     {BondedKernelFlavor::ForcesOnly, true},
 *     {BondedKernelFlavor::ForcesAndVirialAndEnergy, false}
 * };
 *
 * // Step 4: Define parameter tuples
 * using BondedInputConfig = std::tuple<int, PbcType>;
 * using BondedExecutionModes = std::tuple<BondedKernelFlavor, bool>;
 * using BondedTestHelper = HardwareAndExecutionTestHelper<BondedInputConfig, BondedExecutionModes>;
 *
 * // Step 5: Create formatters for input config and execution modes
 * static const auto sc_configFormatters = std::make_tuple(
 *     [](int n) { return formatString("_%dAtoms", n); },
 *     [](PbcType pbc) { return formatString("_%s", c_pbcTypeNames[pbc].c_str()); });
 *
 * static const auto sc_executionFormatters = std::make_tuple(
 *     [](BondedKernelFlavor f) {
 *         return f == BondedKernelFlavor::ForcesOnly ? "_ForcesOnly" : "_All";
 *     },
 *     [](bool simd) { return simd ? "_SIMD" : ""; });
 *
 * // Step 6: Create test namer (used at instantiation)
 * static const auto sc_testNamer =
 *     BondedTestHelper::testNamer(sc_configFormatters, sc_executionFormatters);
 *
 * // Step 7: Define test fixture with custom skip logic
 * class BondedTest : public HardwareTestFixture<BondedTestHelper>
 * {
 * protected:
 *     BondedKernelFlavor flavor_;
 *     bool useSIMD_;
 *
 *     // Pass config formatters directly - RefDataFilenameMaker created automatically
 *     BondedTest() : HardwareTestFixture(sc_configFormatters)
 *     {
 *         // Extract execution modes from parameters
 *         auto [_, numAtoms, pbcType, flavor, useSIMD, __] = GetParam();
 *         flavor_ = flavor;
 *         useSIMD_ = useSIMD;
 *     }
 *
 *     // Override to add parameter-dependent skip logic
 *     void addCustomSkipReasons(MessageStringCollector& skipReasons) override
 *     {
 *         // Example: GPU only supports ForcesOnly flavor
 *         skipReasons.appendIf(
 *             isGpuTest() && flavor_ != BondedKernelFlavor::ForcesOnly,
 *             "GPU does not support ForcesAndVirialAndEnergy mode");
 *
 *         // Example: GPU doesn't use CPU SIMD paths
 *         skipReasons.appendIf(
 *             isGpuTest() && useSIMD_,
 *             "SIMD execution mode only applies to CPU");
 *     }
 * };
 *
 * // Step 8: Write test body
 * TEST_P(BondedTest, Works)
 * {
 *     auto [_, numAtoms, pbcType, flavor, useSIMD, __] = GetParam();
 *     activateHardware();  // Activate GPU device if needed
 *
 *     // IMPORTANT: Skip ForcesOnly when generating reference data
 *     // Reference data should be complete (Forces + Virial + Energy).
 *     // Then ForcesOnly execution can check against a subset of that complete data.
 *     if (flavor == BondedKernelFlavor::ForcesOnly && referenceDataMode() != ReferenceDataMode::Compare)
 *     {
 *         GTEST_SKIP() << "Skipping ForcesOnly mode in reference data generation - "
 *                      << "use ForcesAndVirialAndEnergy to generate complete reference data";
 *     }
 *
 *     // Your test logic here - framework handles CPU/GPU differences
 *     if (isGpuTest()) {
 *         // GPU-specific code
 *         computeBondedGpu(numAtoms, pbcType, flavor, deviceContext(), deviceStream());
 *     } else {
 *         // CPU-specific code
 *         computeBondedCpu(numAtoms, pbcType, flavor, useSIMD);
 *     }
 *
 *     // Always check forces (present in all execution modes)
 *     checker().checkSequence(forces.begin(), forces.end(), "Forces");
 *
 *     if (flavor == BondedKernelFlavor::ForcesAndVirialAndEnergy) {
 *         checker().checkReal(virial, "Virial");
 *         checker().checkReal(energy, "Energy");
 *     }
 *     else
 *     {
 *         // ForcesOnly doesn't compute or check Virial/Energy - tell checker
 *         // to expect unused entries This prevents "unused reference data"
 *         // warnings when ForcesOnly runs
 *         checker().disableUnusedEntriesCheck();
 *     }
 *
 *     // Pattern: Skip incomplete execution modes during ref data generation,
 *     // then check subsets conditionally when those modes run against complete ref data.
 *     // Use disableUnusedEntriesCheck() when execution modes intentionally skip data.
 * }
 *
 * // Step 8: Instantiate test with all combinations
 * INSTANTIATE_TEST_SUITE_P(
 *     AllHardware, BondedTest,
 *     ::testing::Combine(
 *         ::testing::ValuesIn(sc_bondedConfigs),
 *         ::testing::ValuesIn(sc_bondedModes),
 *         ::testing::ValuesIn(getHardwareContextsWithCapability(
 *             GpuConfigurationCapabilities::Bonded))),
 *     sc_testNamer);
 * \endcode
 *
 * This creates test names like:
 * - "AllHardware/BondedTest.Works/10Atoms_Xyz_ForcesOnly_CPU"
 * - "AllHardware/BondedTest.Works/10Atoms_Xyz_ForcesOnly_SIMD_CPU"
 * - "AllHardware/BondedTest.Works/100Atoms_No_All_CPU"
 * - "AllHardware/BondedTest.Works/10Atoms_Xyz_ForcesOnly_GPU0"
 * - (10Atoms_Xyz_ForcesOnly_SIMD_GPU0 is skipped - SIMD only applies to CPU)
 * - (10Atoms_Xyz_All_GPU0 is skipped - GPU doesn't support All flavor)
 *
 * And reference data files like (execution modes and hardware omitted):
 * - "BondedTest_Works_10Atoms_Xyz.xml"
 * - "BondedTest_Works_100Atoms_No.xml"
 *
 * Key points:
 * - CPU and GPU share the same reference data files
 * - Execution modes don't affect reference data filenames
 * - Custom skip logic can filter parameter combinations at runtime
 * - Capability filtering prevents instantiation of unsupported hardware combinations
 *
 */
#ifndef GMX_TESTUTILS_HARDWARE_TEST_FIXTURE_H
#define GMX_TESTUTILS_HARDWARE_TEST_FIXTURE_H

#include "config.h"

#include <optional>
#include <string>

#include <gtest/gtest.h>

#include "gromacs/gpu_utils/capabilities.h"
#include "gromacs/utility/arrayref.h"
#include "gromacs/utility/basedefinitions.h"
#include "gromacs/utility/gmxassert.h"
#include "gromacs/utility/message_string_collector.h"
#include "gromacs/utility/stringutil.h"

#include "testutils/naming.h"
#include "testutils/refdata.h"

class DeviceContext;
class DeviceStream;

namespace gmx
{

namespace test
{

class TestDevice;
class TestReferenceData;
class TestReferenceChecker;

/*! \libinternal
 * \brief Hardware context for parameterized tests
 *
 * Represents either CPU or GPU execution context. Encapsulates device
 * selection, activation, and context/stream access. Supports test
 * naming and tracing with string representations. See file-level
 * documentation for usage examples.
 */
class TestHardwareContext
{
public:
    //! Constructor (nullptr for CPU context, valid pointer for GPU)
    explicit TestHardwareContext(TestDevice* testDevice = nullptr) : testDevice_(testDevice) {}

    //! Returns human-readable context description
    std::string description() const;

    //! Returns sanitized name suitable for GoogleTest parameter names
    std::string testName() const;

    //! Returns GPU ID if this is a GPU context, std::nullopt otherwise
    std::optional<int> gpuId() const;

    //! Returns true if this is a GPU context
    bool isGpuTest() const { return testDevice_ != nullptr; }

    //! Get device context (nullptr for CPU)
    const DeviceContext* deviceContext() const;

    //! Get device stream (nullptr for CPU)
    const DeviceStream* deviceStream() const;

    //! Activate the context (set the device for GPU)
    void activate() const;

private:
    //! Pointer to the test device (nullptr for CPU)
    TestDevice* testDevice_ = nullptr;
};

/*! \brief Get available hardware contexts (CPU + detected GPUs)
 *
 * Returns a view of statically-allocated hardware contexts. The contexts
 * are initialized lazily on first call. Storage lifetime exceeds test
 * execution, making it safe to store pointers in test parameters.
 *
 * The first context is always CPU. Remaining contexts (if any) are GPUs
 * detected by the test hardware environment.
 */
ArrayRef<const TestHardwareContext> getTestHardwareContexts();

/*! \brief Get hardware contexts filtered by GPU capability requirement
 *
 * Returns pointers to hardware contexts where GPU contexts are included
 * only if the specified capability is supported. CPU context is always included.
 *
 * \param[in] hardwareHasCapability  The GPU capability required (e.g., GpuConfigurationCapabilities::Update)
 * \return Vector of pointers to matching hardware contexts (backed by static storage)
 *
 * See file-level documentation for complete usage examples.
 */
std::vector<const TestHardwareContext*> getHardwareContextsWithCapability(bool hardwareHasCapability);


/*! \libinternal
 * \brief Helper for hardware-parameterized tests with input config and execution modes
 *
 * Parameterizes tests across three layers:
 * - **Input configuration**: What to compute (affects reference data)
 * - **Execution modes**: How to compute (affects test names only)
 * - **Hardware**: Where to compute (CPU/GPU - automatic)
 *
 * Key property: Execution modes change how results are computed, but when the same
 * physical quantity is computed, results match regardless of mode. Therefore execution
 * modes affect test names but not reference data names.
 *
 * See file-level documentation for complete usage examples.
 *
 * \tparam InputConfig      Input configuration tuple (affects reference data)
 * \tparam ExecutionModes   Execution mode tuple (skipped from reference data)
 */
template<typename InputConfig, typename ExecutionModes>
class HardwareAndExecutionTestHelper
{
private:
    /*! \brief Create tuple of toEmptyString<T> formatters for each type in a tuple
     *
     * Given a tuple type like std::tuple<int, bool, float>, generates
     * std::make_tuple(toEmptyString<int>, toEmptyString<bool>, toEmptyString<float>).
     *
     * This is used to automatically skip execution mode parameters from reference data
     * filenames by generating formatters that return empty strings.
     *
     * C++17 implementation: Requires two-function pattern with std::index_sequence
     * to generate compile-time indices (0, 1, 2...), then use those indices to
     * extract each tuple element type via std::tuple_element_t.
     *
     * C++20 simplification: Template lambdas eliminate the helper function:
     * \code
     *   return []<std::size_t... Is>(std::index_sequence<Is...>) {
     *       return std::make_tuple(toEmptyString<std::tuple_element_t<Is, Tuple>>...);
     *   }(std::make_index_sequence<std::tuple_size_v<Tuple>>{});
     * \endcode
     *
     * \tparam Tuple      The tuple type to generate formatters for
     * \tparam Is         Index pack (0, 1, 2, ..., N-1) for tuple element access
     */
    template<typename Tuple, std::size_t... Is>
    static auto makeEmptyStringTupleImpl(std::index_sequence<Is...> /*unused*/)
    {
        // For each index I, get the I-th type from Tuple and create toEmptyString<T>
        return std::make_tuple(toEmptyString<std::tuple_element_t<Is, Tuple>>...);
    }

    //! Create tuple of empty-string formatters (see makeEmptyStringTupleImpl for details)
    template<typename Tuple>
    static auto makeEmptyStringTuple()
    {
        // Generate index sequence [0, 1, 2, ..., tuple_size-1] and forward to impl
        return makeEmptyStringTupleImpl<Tuple>(std::make_index_sequence<std::tuple_size_v<Tuple>>{});
    }

public:
    //! Full parameter tuple: input config + execution modes + hardware context
    using DynamicParameters =
            decltype(std::tuple_cat(std::declval<InputConfig>(),
                                    std::declval<ExecutionModes>(),
                                    std::declval<std::tuple<const TestHardwareContext*>>()));

    //! Type for input configuration formatters tuple
    using InputConfigFormatters = typename detail::ParamsToFormatterVariants<InputConfig>::type;

    //! Type for execution mode formatters tuple
    using ExecutionModeFormatters = typename detail::ParamsToFormatterVariants<ExecutionModes>::type;

    //! Create test namer, using all parameters
    static auto testNamer(InputConfigFormatters inputFormatters, ExecutionModeFormatters executionFormatters)
    {
        return NameOfTestFromTuple<DynamicParameters>{ std::tuple_cat(
                inputFormatters,
                executionFormatters,
                std::make_tuple([](const TestHardwareContext* ctx) { return ctx->testName(); })) };
    }
    //! Create reference data filename maker, using only input-config paramters
    static auto refDataFilenameMaker(InputConfigFormatters inputFormatters)
    {
        auto executionModeSkippers = makeEmptyStringTuple<ExecutionModes>();
        return RefDataFilenameMaker<DynamicParameters>{ std::tuple_cat(
                inputFormatters,
                executionModeSkippers,
                std::make_tuple(toEmptyString<const TestHardwareContext*>)) };
    }
};

/*! \brief Flatten nested tuples from Combine(ValuesIn(tuples), ValuesIn(...))
 *
 * Sometimes we want to paramterize tests on all-vs-all across
 * parameter space via testing::Combine(...), and other times we want
 * to be more selective and specify e.g. a list of TestConfig
 * tuples. This helper is used with testing::ConvertGenerator(...) to
 * flatten the inner tuple in the latter case.
 *
 * \tparam ConfigTuple  The input configuration tuple type
 * \returns Lambda that flattens tuple<ConfigTuple, const TestHardwareContext*>
 */
template<typename ConfigTuple>
auto flattenTupleWithHardwareContext()
{
    return [](const std::tuple<ConfigTuple, const TestHardwareContext*>& nested)
    {
        // Unpack the config tuple and append the hardware context
        return std::tuple_cat(std::get<0>(nested), std::make_tuple(std::get<1>(nested)));
    };
}

/*! \libinternal
 * \brief Base class for hardware-parameterized tests
 *
 * Provides infrastructure for tests using HardwareAndExecutionTestHelper:
 * - Safe reference data and checker lifecycle management
 * - Hardware activation convenience method (activateHardware())
 * - Extensible custom skip logic via addCustomSkipReasons() virtual hook
 *
 * See file-level documentation for complete usage examples.
 *
 * \tparam TestHelper  A HardwareAndExecutionTestHelper instantiation that provides
 *                     the parameterization over input configs and execution modes.
 *
 * \warning Do NOT call virtual functions from constructors! Use SetUp() for skip logic.
 *
 * \see HardwareAndExecutionTestHelper, TestHardwareContext
 */
template<typename TestHelper>
class HardwareTestFixture : public ::testing::TestWithParam<typename TestHelper::DynamicParameters>
{
public:
    //! Parameter tuple type, including input config, execution mode, and hardware contexts
    using ParametersTuple = typename TestHelper::DynamicParameters;

protected:
    /*! \brief Constructor
     *
     * Takes input config formatters and automatically creates the RefDataFilenameMaker
     * internally. See file-level documentation for complete usage examples.
     *
     * \param[in] inputFormatters  Tuple of formatter functions for input config parameters
     *
     * \note High-level hardware capability filtering should be done at test
     *        instantiation via testing::Combine with getHardwareContextsWithCapability().
     *       Custom parameter-dependent skip logic can be added via the virtual
     *       addCustomSkipReasons() hook, which runs in SetUp(). The reference data
     *       checker is created in SetUp() after all skip checks complete.
     */
    explicit HardwareTestFixture(typename TestHelper::InputConfigFormatters inputFormatters) :
        refData_(TestHelper::refDataFilenameMaker(inputFormatters)(this->GetParam()))
    {
        hardwareContext_ = std::get<std::tuple_size_v<ParametersTuple> - 1>(this->GetParam());
    }

    /*! \brief SetUp override - performs custom skip checks and creates checker
     *
     * Called by GoogleTest after construction is complete. This is where
     * custom skip logic runs via the virtual addCustomSkipReasons() hook.
     *
     * If derived classes override SetUp(), they MUST call base SetUp() first:
     *
     * \code{.cpp}
     * void SetUp() override {
     *     HardwareTestFixture::SetUp();
     *     // Your setup here...
     * }
     * \endcode
     */
    void SetUp() override
    {
        // Call parent SetUp() first (GoogleTest best practice)
        ::testing::TestWithParam<ParametersTuple>::SetUp();

        MessageStringCollector skipReasons;

        // Call virtual hook for custom skip logic - dynamic dispatch works correctly here!
        addCustomSkipReasons(skipReasons);

        // Skip if any reasons accumulated
        if (!skipReasons.isEmpty())
        {
            GTEST_SKIP() << skipReasons.toString();
        }

        // NOW safe to create checker (after all skip checks)
        checker_ = refData_.rootChecker();
    }

    /*! \brief Virtual hook for derived classes to add custom skip conditions
     *
     * Called from SetUp() (NOT constructor) to add test-specific skip logic.
     * Derived classes can override this to add parameter-dependent skip conditions.
     *
     * \param[in,out] skipReasons  Collector to append skip reasons to
     *
     * Default implementation does nothing (no custom skip logic).
     */
    virtual void addCustomSkipReasons(MessageStringCollector& skipReasons)
    {
        // Default: no custom skip logic
        // Derived classes override to add parameter-dependent skip conditions
        GMX_UNUSED_VALUE(skipReasons);
    }

    //! Get hardware context (extracted from parameter tuple)
    const TestHardwareContext* hardwareContext() const { return hardwareContext_; }

    /*! \brief Get reference checker
     *
     * Safe to call after SetUp() completes. Asserts if called before
     * checker is initialized (which should be impossible if used correctly).
     */
    TestReferenceChecker& checker()
    {
        GMX_RELEASE_ASSERT(checker_.has_value(), "Checker not initialized - was SetUp() called?");
        return *checker_;
    }

    //! Activate hardware (sets device for GPU, no-op for CPU)
    void activateHardware() const { hardwareContext_->activate(); }

    //! Check if this is a GPU test
    bool isGpuTest() const { return hardwareContext_->isGpuTest(); }

    //! Get GPU ID if GPU test, std::nullopt otherwise
    std::optional<int> gpuId() const { return hardwareContext_->gpuId(); }

    //! Get device context (nullptr for CPU, convenience wrapper for hardwareContext()->deviceContext())
    const DeviceContext* deviceContext() const { return hardwareContext_->deviceContext(); }

    //! Get device stream (nullptr for CPU, convenience wrapper for hardwareContext()->deviceStream())
    const DeviceStream* deviceStream() const { return hardwareContext_->deviceStream(); }

private:
    //! Reference data (created in constructor, before skip checks)
    TestReferenceData refData_;

    //! Reference checker (created in SetUp() after skip checks, hence optional)
    std::optional<TestReferenceChecker> checker_;

    //! Hardware context (extracted from parameters for convenience)
    const TestHardwareContext* hardwareContext_;
};


} // namespace test
} // namespace gmx

#endif // GMX_TESTUTILS_HARDWARE_TEST_FIXTURE_H
