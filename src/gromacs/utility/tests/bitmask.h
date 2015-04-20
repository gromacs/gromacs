/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2014, by the GROMACS development team, led by
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
 * Tests for bitmask functionality.
 *
 * These tests check the functionality of bitmask.h

 * \author Roland Schulz <roland@rschulz.eu>
 * \ingroup module_utility
 */
#include <gtest/gtest.h>

#include "gromacs/utility/bitmask.h"

//! Implemenation of BITMASK_CLASSNAME
#define BITMASK_CLASSNAME_(S) BitmaskTest ## S
//! Returns name of Bitmask test fixture class
#define BITMASK_CLASSNAME(S) BITMASK_CLASSNAME_(S)
//! Implementation of BITMASK_TEST_P
#define BITMASK_TEST_P_(C, T) TEST_P(C, T)
//! Defines a parameterized bitmask test
#define BITMASK_TEST_P(T) BITMASK_TEST_P_(BITMASK_CLASSNAME(BITMASK_SIZE), T)

class BITMASK_CLASSNAME(BITMASK_SIZE) : public ::testing::TestWithParam<int>
{
};

BITMASK_TEST_P(SetAndClear)
{
    gmx_bitmask_t m;
    int           i = GetParam();
    bitmask_clear(&m);
    EXPECT_TRUE(bitmask_is_zero(m));
    EXPECT_FALSE(bitmask_is_set(m, i));
    bitmask_set_bit(&m, i);
    for (int j = 0; j < BITMASK_SIZE; j++)
    {
        EXPECT_EQ(bitmask_is_set(m, j), j == i);
    }
    bitmask_clear(&m);
    EXPECT_TRUE(bitmask_is_zero(m));
}

BITMASK_TEST_P(InitBit)
{
    gmx_bitmask_t m1, m2;
    int           i = GetParam();
    bitmask_init_bit(&m1, i);
    bitmask_clear(&m2);
    EXPECT_FALSE(bitmask_is_equal(m1, m2));
    bitmask_set_bit(&m2, i);
    EXPECT_TRUE(bitmask_is_equal(m1, m2));
}

BITMASK_TEST_P(InitLowBits)
{
    gmx_bitmask_t m;
    int           i = GetParam();
    bitmask_init_low_bits(&m, i);
    for (int j = 0; j < BITMASK_SIZE; j++)
    {
        EXPECT_EQ(bitmask_is_set(m, j), j < i);
    }
}

BITMASK_TEST_P(Disjoint)
{
    gmx_bitmask_t m1, m2;
    int           i = GetParam();
    bitmask_init_bit(&m1, i);
    bitmask_init_bit(&m2, i);
    EXPECT_FALSE(bitmask_is_disjoint(m1, m2));
    bitmask_init_low_bits(&m2, i);
    EXPECT_TRUE(bitmask_is_disjoint(m1, m2));
}

BITMASK_TEST_P(Union)
{
    gmx_bitmask_t m1, m2;
    int           i = GetParam();
    int           j = (i + BITMASK_SIZE/2)%BITMASK_SIZE;
    bitmask_init_bit(&m1, i);
    bitmask_init_bit(&m2, j);
    bitmask_union(&m1, m2);
    for (int k = 0; k < BITMASK_SIZE; k++)
    {
        EXPECT_EQ(bitmask_is_set(m1, k), k == i || k == j);
    }

    bitmask_init_bit(&m1, i);
    bitmask_clear(&m2);
    bitmask_union(&m1, m2);
    bitmask_init_bit(&m2, i);
    EXPECT_TRUE(bitmask_is_equal(m1, m2));

    bitmask_clear(&m1);
    bitmask_init_bit(&m2, i);
    bitmask_union(&m1, m2);
    EXPECT_TRUE(bitmask_is_equal(m1, m2));
}
