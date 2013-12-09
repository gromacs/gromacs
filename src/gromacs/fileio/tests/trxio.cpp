/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2013, by the GROMACS development team, led by
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
 * Tests utilities for trxio routines
 *
 * \author Mark Abraham <mark.j.abraham@gmail.com>
 * \ingroup module_fileio
 */
#include "../trxio.h"

#include <gtest/gtest.h>

namespace
{

//! Helper function for EXPECT_*
::testing::AssertionResult RunningBRealModulo(double a, double b, double c)
{
    if (bRealModulo(a, b, c))
    {
        return ::testing::AssertionSuccess();
    }
    else
    {
        return ::testing::AssertionFailure();
    }
}

//! Helper typedef
typedef ::testing::Test bRealModuloFunction;

TEST_F(bRealModuloFunction, HandlesInputsCorrectly)
{
    EXPECT_TRUE(RunningBRealModulo(2, 0, 1));
    /* There used to be a bug (#1397) such that
     *
     * gmx eneconv -offset 1.998 -dt 2
     *
     * would work correctly up to t=16777, and then fail by writing
     * extra frames because of some floating-point catastrophe. */
    EXPECT_FALSE(RunningBRealModulo(16777.996, 1.998, 2));
    EXPECT_TRUE(RunningBRealModulo(16777.998, 1.998, 2));
    EXPECT_FALSE(RunningBRealModulo(16778.000, 1.998, 2));
    /* Also, check that some differences between large times work
     * correctly. */
    EXPECT_FALSE(RunningBRealModulo(16779.998, 16779.996, 2));
    EXPECT_TRUE(RunningBRealModulo(16779.998, 16777.998, 2));
    EXPECT_TRUE(RunningBRealModulo(16780.000, 16778.000, 2));
}

} // namespace
