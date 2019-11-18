/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2019, by the GROMACS development team, led by
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
 * \brief
 * Tests PBC code
 *
 * \author Mark Abraham <mark.j.abraham@gmail.com>
 * \ingroup module_pbcutil
 */
#include "gmxpre.h"

#include "gromacs/pbcutil/pbc.h"

#include <gtest/gtest.h>

#include "gromacs/math/vectypes.h"
#include "gromacs/pbcutil/ishift.h"

#include "testutils/refdata.h"

namespace gmx
{

namespace test
{

TEST(PbcTest, CalcShiftsWorks)
{
    // Choose box vector entries whose magnitudes will lead to unique
    // shift vector values when the largest box shift in any dimension
    // is two.
    const matrix box = { { 0.01, 1, -100 }, { 300, -0.03, 3 }, { -6, -600, 0.06 } };
    rvec         shiftVectors[SHIFTS];

    calc_shifts(box, shiftVectors);

    gmx::test::TestReferenceData    data;
    gmx::test::TestReferenceChecker checker(data.rootChecker());
    checker.checkSequence(std::begin(shiftVectors), std::end(shiftVectors), "ShiftVectors");
}

} // namespace test

} // namespace gmx
