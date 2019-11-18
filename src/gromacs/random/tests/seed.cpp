/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2015,2016,2018,2019, by the GROMACS development team, led by
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
 * \brief Tests for random seed functions
 *
 * \author Erik Lindahl <erik.lindahl@gmail.com>
 * \ingroup module_random
 */
#include "gmxpre.h"

#include "gromacs/random/seed.h"

#include <gtest/gtest.h>

namespace gmx
{

namespace
{

// Test the random device call
TEST(SeedTest, makeRandomSeed)
{
    // Unlike Sony, we do not use "4" as a constant random value, so the only
    // thing we can check for the random device is that multiple calls to
    // it produce different results.
    // We choose to ignore the 2^-64 probability this will happen by chance;
    // if you execute the unit tests once per second you might have to run them
    // an extra time rougly once per 300 billion years - apologies in advance!

    uint64_t i0 = makeRandomSeed();
    uint64_t i1 = makeRandomSeed();

    EXPECT_NE(i0, i1);
}

} // namespace

} // namespace gmx
