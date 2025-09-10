/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2012- The GROMACS Authors
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
 * Functions for initialing \Gromacs unit test executables.
 *
 * \author Teemu Murtola <teemu.murtola@gmail.com>
 * \ingroup module_testutils
 */
#ifndef GMX_TESTUTILS_TESTINIT_H
#define GMX_TESTUTILS_TESTINIT_H

#include <filesystem>
#include <string>

#include <gtest/gtest.h>

#include "naming.h"

namespace gmx
{

namespace test
{

//! \cond internal
/*! \internal
 * \brief
 * Initializes the test utilities library.
 *
 * Does not throw.  Terminates the program with a non-zero error code if an
 * error occurs.
 *
 * This function is automatically called by unittest_main.cpp.
 *
 * \param[in] dataPath Filepath to input files.
 * \param[in] tempPath Filepath to temporary files.
 * \param[in] usesMpi  If the test is run with MPI or not.
 * \param[in] usesHardwareDetection If hardwaredetection is enabled.
 * \param[in] registersDynamically  Whether dynamical GoogleTest registration
 *                                  is used by this test binary
 * \param[in] argc Number of cmdline options
 * \param[in] argv Cmdline options.
 *
 * \ingroup module_testutils
 */
void initTestUtils(const std::filesystem::path& dataPath,
                   const std::filesystem::path& tempPath,
                   bool                         usesMpi,
                   bool                         usesHardwareDetection,
                   bool                         registersDynamically,
                   int*                         argc,
                   char***                      argv);

/*! \internal
 * \brief
 * Finalizes the test utilities library.
 *
 * \param[in] usesHardwareDetection If hardwaredetection is enabled.
 * \param[in] registersDynamically  Whether dynamical GoogleTest registration
 *                                  is used by this test binary
 *
 * Does not throw.  Terminates the program with a non-zero error code if an
 * error occurs.
 *
 * This function is automatically called by unittest_main.cpp.
 *
 * \ingroup module_testutils
 */
void finalizeTestUtils(bool usesHardwareDetection, bool registersDynamically);
//! \endcond

/*! \brief Declaration of function used to dynamically register
 * GoogleTest tests.
 *
 * When a test binary uses gmx_add_gtest_executable(exename
 * DYNAMIC_REGISTRATION ...) it must provide an implementation of this
 * method. The method is called before RUN_ALL_TESTS() and is expected
 * to call testing::RegisterTest to register tests dynamically. This
 * approach might be necessary to run the tests in a stable way over
 * whichever hardware is detected at run time.
 *
 * Normal test binaries do not need to implement this function.
 *
 * ::gmx::test::TestHardwareEnvironment::gmxSetUp() should be called before
 * this method, in case the test hardware environment is needed to help
 * decide which tests to register. */
void registerTestsDynamically();

/*! \brief Register tests dynamically based on the execution context
 *
 * This template reduces code duplication across the dynamically
 * registered tests, letting them register their tests more tersely.
 *
 * \tparam TestFixture  The type of the test fixture
 * \tparam TestCase     The type of the test case (derived from \c TestFixture)
 * \tparam Combinations The interal GoogleTest type describing the
 *                        return from \c testing::Combine() intended
 *                        to generate test parameters of type \c
 *                        TestFixture::ParamType (which is typically a
 *                        tuple).
 *
 * \param[in] testSuiteName           The name of the test suite that shares the \c TestFixture
 * \param[in] testNamer               Functor to make the full name of the test case
 * \param[in] combinations            A generator of test values produced with
 *                                      \c testing::Combine()
 *
 * Note that \c Combinations is actually well defined relative to \c
 * TestFixture::ParamType, but its concrete type is an internal
 * GoogleTest type, so we do not want to express it in code. In C++20,
 * it would be better to declare the function parameter like `const
 * auto& combinations` to achieve the same effect of hiding the
 * concrete type. */
template<typename TestFixture, typename TestCase, typename Combinations>
void registerTests(const std::string&                                          testSuiteName,
                   const NameOfTestFromTuple<typename TestFixture::ParamType>& testNamer,
                   const Combinations&                                         combinations)
{
    // It is not good practice to use GoogleTest's internal type here,
    // but it's a practical alternative that lets us use GoogleTest's Combine
    // in the code that declares the values passed to registered tests.
    //
    // Normally the use of this type is hidden behind a call to a
    // INSTANTIATE_TEST_SUITE_P macro. If GoogleTest do change things
    // such that this breaks, it may be simple to fix it. Or if not,
    // we can always manually build or enumerate the Cartesian product
    // of values that it generates and use that in the loop below.
    const auto testParamGenerator =
            ::testing::internal::ParamGenerator<typename TestFixture::ParamType>(combinations);
    for (const auto& parameters : testParamGenerator)
    {
        ::testing::TestParamInfo<typename TestFixture::ParamType> testParamInfo(parameters, 0);
        const std::string testName = testNamer(testParamInfo);
        // This returns a testing::TestInfo object that leaks. That's currently not
        // a problem for a short-lived test binary. But if it did become one, then we
        // should fill and return a std::vector<std::unique_ptr<TestInfo>> whose lifetime
        // is managed in testinit.cpp.
        ::testing::RegisterTest(testSuiteName.c_str(),
                                testName.c_str(),
                                nullptr,
                                testName.c_str(),
                                __FILE__,
                                __LINE__,
                                // Important to use the fixture type as the return type here, even
                                // though we construct an object of the derived type.
                                [=]() -> TestFixture* { return new TestCase(parameters); });
    }
}

/*! \brief DEPRECATED Register tests dynamically based on the execution context
 *
 * New tests should use the above registerTests method. This method
 * will be removed as part of #5027.
 *
 * This template reduces code duplication across the dynamically
 * registered tests, letting them register their tests more tersely.
 *
 * \tparam TestFixture  The type of the test fixture
 * \tparam TestCase     The type of the test case (derived from \c TestFixture)
 * \tparam Combinations The interal GoogleTest type describing the
 *                        return from \c testing::Combine() intended
 *                        to generate test parameters of type \c
 *                        TestFixture::ParamType (which is typically a
 *                        tuple).
 *
 * \param[in] testSuiteName           The name of the test suite that shares the \c TestFixture
 * \param[in] makeBriefNameOfTestCase Function that will make the brief name of the test case,
 *                                      used for naming the refdata file
 * \param[in] makeFullNameOfTestCase  Function that will make the full name of the test case
 * \param[in] combinations            A generator of test values produced with
 *                                      \c testing::Combine()
 *
 * Note that \c Combinations is actually well defined relative to \c
 * TestFixture::ParamType, but its concrete type is an internal
 * GoogleTest type, so we do not want to express it in code. In C++20,
 * it would be better to declare the function parameter like `const
 * auto& combinations` to achieve the same effect of hiding the
 * concrete type. */
template<typename TestFixture, typename TestCase, typename Combinations>
void registerTests(const std::string& testSuiteName,
                   std::string (*makeBriefNameOfTestCase)(
                           const typename ::testing::TestParamInfo<typename TestFixture::ParamType>&),
                   std::string (*makeFullNameOfTestCase)(
                           const typename ::testing::TestParamInfo<typename TestFixture::ParamType>&,
                           const std::string&),
                   const Combinations& combinations)
{
    // It is not good practice to use GoogleTest's internal type here,
    // but it's a practical alternative that lets us use GoogleTest's Combine
    // in the code that declares the values passed to registered tests.
    //
    // Normally the use of this type is hidden behind a call to a
    // INSTANTIATE_TEST_SUITE_P macro. If GoogleTest do change things
    // such that this breaks, it may be simple to fix it. Or if not,
    // we can always manually build or enumerate the Cartesian product
    // of values that it generates and use that in the loop below.
    const auto testParamGenerator =
            ::testing::internal::ParamGenerator<typename TestFixture::ParamType>(combinations);
    for (const auto& parameters : testParamGenerator)
    {
        ::testing::TestParamInfo<typename TestFixture::ParamType> testParamInfo(parameters, 0);
        std::string testName = makeBriefNameOfTestCase(testParamInfo);
        // This returns a testing::TestInfo object that leaks. That's currently not
        // a problem for a short-lived test binary. But if it did become one, then we
        // should fill and return a std::vector<std::unique_ptr<TestInfo>> whose lifetime
        // is managed in testinit.cpp.
        ::testing::RegisterTest(testSuiteName.c_str(),
                                makeFullNameOfTestCase(testParamInfo, testName).c_str(),
                                nullptr,
                                testName.c_str(),
                                __FILE__,
                                __LINE__,
                                // Important to use the fixture type as the return type here, even
                                // though we construct an object of the derived type.
                                [=]() -> TestFixture* { return new TestCase(parameters); });
    }
}

} // namespace test
} // namespace gmx

#endif
