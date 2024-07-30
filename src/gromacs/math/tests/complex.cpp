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
 * Tests for math operations on complex numbers.
 *
 * \author Paul Bauer <paul.bauer.q@gmail.com>
 * \ingroup module_math
 */
#include "gmxpre.h"

#include <string>
#include <vector>

#include <gtest/gtest.h>

#include "gromacs/math/gmxcomplex.h"
#include "gromacs/utility/real.h"

#include "testutils/refdata.h"

namespace gmx
{
namespace test
{
namespace
{

TEST(ComplexNumberTest, RealComplexMultiply)
{
    gmx::test::TestReferenceData    data;
    gmx::test::TestReferenceChecker checker(data.rootChecker());

    std::vector<t_complex> result;

    result.emplace_back(rcmul(2.5, { 13.37, 42.23 }));
    result.emplace_back(rcmul(0, { 13.37, 42.23 }));
    result.emplace_back(rcmul(2.5, { 0, 42.23 }));
    result.emplace_back(rcmul(2.5, { 13.37, 0 }));

    checker.checkSequence(result.begin(), result.end(), "RealComplexMultiply");
}

TEST(ComplexNumberTest, RealComplexExp)
{
    gmx::test::TestReferenceData    data;
    gmx::test::TestReferenceChecker checker(data.rootChecker());

    std::vector<t_complex> result;

    result.emplace_back(rcexp(2.5));
    result.emplace_back(rcexp(13.37));
    result.emplace_back(rcexp(0.0));
    result.emplace_back(rcexp(-2.3));

    checker.checkSequence(result.begin(), result.end(), "RealComplexExp");
}

TEST(ComplexNumberTest, ComplexAdd)
{
    gmx::test::TestReferenceData    data;
    gmx::test::TestReferenceChecker checker(data.rootChecker());

    std::vector<t_complex> result;

    result.emplace_back(cadd({ 13.37, 42.23 }, { 2.5, 1.2 }));
    result.emplace_back(cadd({ 13.37, 0 }, { 2.5, 1.2 }));
    result.emplace_back(cadd({ 0, 42.23 }, { 2.5, 1.2 }));
    result.emplace_back(cadd({ 13.37, 42.23 }, { 0, 1.2 }));
    result.emplace_back(cadd({ 13.37, 42.23 }, { 2.5, 0 }));
    result.emplace_back(cadd({ 13.37, 42.23 }, { 0, 0 }));
    result.emplace_back(cadd({ 0, 0 }, { 2.5, 1.2 }));

    checker.checkSequence(result.begin(), result.end(), "ComplexAdd");
}

TEST(ComplexNumberTest, ComplexSubtract)
{
    gmx::test::TestReferenceData    data;
    gmx::test::TestReferenceChecker checker(data.rootChecker());

    std::vector<t_complex> result;

    result.emplace_back(csub({ 13.37, 42.23 }, { 2.5, 1.2 }));
    result.emplace_back(csub({ 13.37, 0 }, { 2.5, 1.2 }));
    result.emplace_back(csub({ 0, 42.23 }, { 2.5, 1.2 }));
    result.emplace_back(csub({ 13.37, 42.23 }, { 0, 1.2 }));
    result.emplace_back(csub({ 13.37, 42.23 }, { 2.5, 0 }));
    result.emplace_back(csub({ 13.37, 42.23 }, { 0, 0 }));
    result.emplace_back(csub({ 0, 0 }, { 2.5, 1.2 }));

    checker.checkSequence(result.begin(), result.end(), "ComplexSubtract");
}

TEST(ComplexNumberTest, ComplexMultiply)
{
    gmx::test::TestReferenceData    data;
    gmx::test::TestReferenceChecker checker(data.rootChecker());

    std::vector<t_complex> result;

    result.emplace_back(cmul({ 13.37, 42.23 }, { 2.5, 1.2 }));
    result.emplace_back(cmul({ 13.37, 0 }, { 2.5, 1.2 }));
    result.emplace_back(cmul({ 0, 42.23 }, { 2.5, 1.2 }));
    result.emplace_back(cmul({ 13.37, 42.23 }, { 0, 1.2 }));
    result.emplace_back(cmul({ 13.37, 42.23 }, { 2.5, 0 }));
    result.emplace_back(cmul({ 13.37, 42.23 }, { 0, 0 }));
    result.emplace_back(cmul({ 0, 0 }, { 2.5, 1.2 }));

    checker.checkSequence(result.begin(), result.end(), "ComplexMultiply");
}

TEST(ComplexNumberTest, ComplexDivision)
{
    gmx::test::TestReferenceData    data;
    gmx::test::TestReferenceChecker checker(data.rootChecker());

    std::vector<t_complex> result;

    result.emplace_back(cdiv({ 13.37, 42.23 }, { 2.5, 1.2 }));
    result.emplace_back(cdiv({ 13.37, 0 }, { 2.5, 1.2 }));
    result.emplace_back(cdiv({ 0, 42.23 }, { 2.5, 1.2 }));
    result.emplace_back(cdiv({ 13.37, 42.23 }, { 2.5, 0 }));
    result.emplace_back(cdiv({ 0, 0 }, { 2.5, 1.2 }));

    checker.checkSequence(result.begin(), result.end(), "ComplexDivision");
}

TEST(ComplexNumberTest, ComplexConjugate)
{
    gmx::test::TestReferenceData    data;
    gmx::test::TestReferenceChecker checker(data.rootChecker());

    std::vector<t_complex> result;

    result.emplace_back(conjugate({ 13.37, 42.23 }));
    result.emplace_back(conjugate({ 13.37, 0 }));
    result.emplace_back(conjugate({ 0, 42.23 }));

    checker.checkSequence(result.begin(), result.end(), "ComplexConjugate");
}

TEST(ComplexNumberTest, ComplexAbs2)
{
    gmx::test::TestReferenceData    data;
    gmx::test::TestReferenceChecker checker(data.rootChecker());

    std::vector<real> result;

    result.emplace_back(cabs2({ 13.37, 42.23 }));
    result.emplace_back(cabs2({ 13.37, 0 }));
    result.emplace_back(cabs2({ 0, 42.23 }));

    checker.checkSequence(result.begin(), result.end(), "ComplexAbs2");
}

} // namespace
} // namespace test
} // namespace gmx
