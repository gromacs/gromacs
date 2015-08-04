/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2014,2015, by the GROMACS development team, led by
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
 * Tests utilities for random number generation.
 *
 * \author Roland Schulz <roland@utk.edu>
 * \ingroup module_random
 */
#include "gmxpre.h"

#include <vector>

#include <gtest/gtest.h>

#include "external/Random123-1.08/include/Random123/threefry.h"

#include "testutils/refdata.h"

namespace
{

class Threefry : public ::testing::TestWithParam<std::pair<threefry2x64_ctr_t,
                                                           threefry2x64_key_t> >
{
};

TEST_P(Threefry, 2x64)
{
    gmx::test::TestReferenceData    data;

    gmx::test::TestReferenceChecker checker(data.rootChecker());

    const std::pair<threefry2x64_ctr_t, threefry2x64_key_t> input = GetParam();

    threefry2x64_ctr_t rand = threefry2x64(input.first, input.second);

    checker.checkSequenceArray(2, rand.v, "Threefry2x64");
}

//The input values are the same as the ones used by Known Answer Tests (kat) in
//Random123. Reference values agree with those in kat_vectors
/** input value: zero */
const threefry2x64_ctr_t tf_zero = {{0, 0}};
/** input value: max unit64 */
const threefry2x64_ctr_t tf_max  = {{std::numeric_limits<gmx_uint64_t>::max(),
                                     std::numeric_limits<gmx_uint64_t>::max()}};
/** input value: Pi */
const threefry2x64_ctr_t tf_pi1  = {{0x243f6a8885a308d3ULL, 0x13198a2e03707344ULL}};
/** input value: More Pi */
const threefry2x64_ctr_t tf_pi2  = {{0xa4093822299f31d0ULL, 0x082efa98ec4e6c89ULL}};

INSTANTIATE_TEST_CASE_P(0_ff_pi, Threefry,
                            ::testing::Values(std::make_pair(tf_zero, tf_zero),
                                              std::make_pair(tf_max, tf_max),
                                              std::make_pair(tf_pi1, tf_pi2)));

} // namespace
