/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2025- The GROMACS Authors
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
 * Tests various text-dumping routines
 *
 * \ingroup module_math
 */
#include "gmxpre.h"

#include "gromacs/math/vecdump.h"

#include <iostream>
#include <vector>

#include <gtest/gtest.h>

#include "gromacs/utility/arrayref.h"
#include "gromacs/utility/stringstream.h"
#include "gromacs/utility/textwriter.h"

#include "testutils/posixmemstream.h"
#include "testutils/refdata.h"
#include "testutils/setenv.h"
#include "testutils/stringtest.h"

namespace gmx
{
namespace test
{
namespace
{

struct TestParameters
{
    bool valuesAvailable;
    bool showIndices;
};

//! GoogleTest helper method for pretty-printing during failure
void PrintTo(const TestParameters& parameters, std::ostream* os)
{
    *os << "Values available: " << parameters.valuesAvailable << "\n";
    *os << "Show indices: " << parameters.showIndices << "\n";
}

class DumpingVectorsTest : public StringTestBase, public testing::WithParamInterface<TestParameters>
{
public:
    template<typename T>
    using LegacyTestFunctionType = void (*)(FILE*, int, const char*, const T[], int, gmx_bool);
    template<typename T>
    void testLegacyFunction(LegacyTestFunctionType<T> legacyFunctionToTest,
                            const T*                  v,
                            int                       size,
                            const char*               id,
                            const bool                testingLongFormat = false)
    {
        PosixMemstream       stream;
        const TestParameters parameters = GetParam();
        legacyFunctionToTest(stream.stream(), amountToIndent_, title_, v, size, parameters.showIndices);
        // If running on a system that supports memory streams, check
        // the buffer contents, except that the long format for RVec
        // output only serializes reliably when the values were in
        // double precision, so we only check its output in that case.
        if (stream.canCheckBufferContents() and (GMX_DOUBLE or !testingLongFormat))
        {
            checkText(stream.toString(), id);
        }
        else
        {
            checker().disableUnusedEntriesCheck();
        }
    }
    template<typename T>
    using TestFunctionType = void (*)(TextWriter*, const char*, ArrayRef<const T>, bool);
    template<typename T>
    void testModernFunction(TestFunctionType<T> functionToTest,
                            ArrayRef<T>         testArrayRef,
                            const char*         id,
                            const bool          testingLongFormat = false)
    {
        StringOutputStream stringStream;
        {
            TextWriter           writer(&stringStream);
            ScopedIndenter       indenter   = writer.addScopedIndentation(amountToIndent_);
            const TestParameters parameters = GetParam();
            functionToTest(&writer, title_, testArrayRef, parameters.showIndices);
        }
        // Check the buffer contents, except that the long format for
        // RVec output only serializes reliably when the values were
        // in double precision, so we only check its output in that
        // case.
        if (GMX_DOUBLE or !testingLongFormat)
        {
            checkText(stringStream.toString(), id);
        }
        else
        {
            checker().disableUnusedEntriesCheck();
        }
    }

private:
    const int   amountToIndent_ = 2;
    const char* title_          = "The title";
};

TEST_P(DumpingVectorsTest, CanDumpIntIdentically)
{
    const TestParameters parameters = GetParam();
    using TestType                  = int;
    std::vector<TestType> vector    = { 1, 2, 3 };
    const auto            testArrayRef =
            parameters.valuesAvailable ? ArrayRef<TestType>(vector) : ArrayRef<TestType>{};
    testLegacyFunction(pr_ivec, testArrayRef.data(), testArrayRef.size(), "output");
    testModernFunction(dumpIntArrayRef, testArrayRef, "output");
}

TEST_P(DumpingVectorsTest, CanDumpIntBlocks)
{
    const TestParameters parameters = GetParam();
    using TestType                  = int;
    std::vector<TestType> vector    = { 1, 2, 3, 4, 5, 10, 11, 12, 13, 99, 127 };
    const auto            testArrayRef =
            parameters.valuesAvailable ? ArrayRef<TestType>(vector) : ArrayRef<TestType>{};
    testLegacyFunction(pr_ivec_block, testArrayRef.data(), testArrayRef.size(), "output");
}

TEST_P(DumpingVectorsTest, CanDumpFloatIdentically)
{
    const TestParameters parameters = GetParam();
    using TestType                  = float;
    std::vector<TestType> vector    = { 1.1, 2.2, 3.3 };
    const auto            testArrayRef =
            parameters.valuesAvailable ? ArrayRef<TestType>(vector) : ArrayRef<TestType>{};
    testLegacyFunction(pr_fvec, testArrayRef.data(), testArrayRef.size(), "output");
    testModernFunction(dumpFloatArrayRef, testArrayRef, "output");
}

TEST_P(DumpingVectorsTest, CanDumpDoubleIdentically)
{
    const TestParameters parameters = GetParam();
    using TestType                  = double;
    std::vector<TestType> vector    = { 1.1, 2.2, 3.3 };
    const auto            testArrayRef =
            parameters.valuesAvailable ? ArrayRef<TestType>(vector) : ArrayRef<TestType>{};
    testLegacyFunction(pr_dvec, testArrayRef.data(), testArrayRef.size(), "output");
    testModernFunction(dumpDoubleArrayRef, testArrayRef, "output");
}

TEST_P(DumpingVectorsTest, CanDumpRealIdentically)
{
    const TestParameters parameters = GetParam();
    using TestType                  = real;
    std::vector<TestType> vector    = { 1.1, 2.2, 3.3 };
    const auto            testArrayRef =
            parameters.valuesAvailable ? ArrayRef<TestType>(vector) : ArrayRef<TestType>{};
    testLegacyFunction(pr_rvec, testArrayRef.data(), testArrayRef.size(), "output");
    testModernFunction(dumpRealArrayRef, testArrayRef, "output");
}

//! Wrapper that ignores the "show indices" boolean
void pr_rvecs_wrapper(FILE* fp, int indent, const char* title, const rvec vec[], int n, gmx_bool /* unused */)
{
    pr_rvecs(fp, indent, title, vec, n);
}

//! Wrapper that ignores the "show indices" boolean
void dumpRvecArrayRefWrapper(TextWriter* writer, const char* description, ArrayRef<const RVec> values, bool /* unused */)
{
    dumpRvecArrayRef(writer, description, values);
}

TEST_P(DumpingVectorsTest, CanDumpRvecIdentically)
{
    const TestParameters parameters = GetParam();
    if (!parameters.showIndices)
    {
        // This facility is not supported for these methods.  We don't
        // use GTEST_SKIP because we don't want to imply to the user
        // of the test binary that something that *might* be tested
        // was not tested.
        return;
    }
    using TestType               = RVec;
    std::vector<TestType> vector = { { 1.1, 2.2, 3.3 }, { 4.4, 5.5, 6.6 }, { 7.7, 8.8, 9.9 } };
    const auto            testArrayRef =
            parameters.valuesAvailable ? ArrayRef<TestType>(vector) : ArrayRef<TestType>{};
    const rvec* rvecPtr = reinterpret_cast<rvec*>(testArrayRef.data());
    testLegacyFunction(pr_rvecs_wrapper, rvecPtr, testArrayRef.size(), "output");
    testModernFunction(dumpRvecArrayRefWrapper, testArrayRef, "output");
    {
        SCOPED_TRACE("In long format");
        const bool overWriteEnvironmentVariable = true;
        gmxSetenv("GMX_PRINT_LONGFORMAT", "1", overWriteEnvironmentVariable);
        const bool testingLongFormat = true;
        testLegacyFunction(pr_rvecs_wrapper, rvecPtr, testArrayRef.size(), "long format", testingLongFormat);
        testModernFunction(dumpRvecArrayRefWrapper, testArrayRef, "long format", testingLongFormat);
        gmxUnsetenv("GMX_PRINT_LONGFORMAT");
    }
}

//! Define the space for testing
TestParameters g_testParameters[] = { { true, true }, { true, false }, { false, true }, { false, false } };
INSTANTIATE_TEST_SUITE_P(Works, DumpingVectorsTest, ::testing::ValuesIn(g_testParameters));

} // namespace
} // namespace test
} // namespace gmx
