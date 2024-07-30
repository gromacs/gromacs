/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2018- The GROMACS Authors
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
 * Tests for the HashedMap class.
 *
 * \author berk Hess <hess@kth.se>
 * \ingroup module_domdec
 */
#include "gmxpre.h"

#include "gromacs/domdec/hashedmap.h"

#include <string>

#include <gtest/gtest.h>

#include "testutils/testasserts.h"

namespace gmx
{
namespace test
{
namespace
{

/*! \brief Checks that the key is found and if so also checks the value */
void checkFinds(const gmx::HashedMap<char>& map, int key, char value)
{
    const char* pointer = map.find(key);
    EXPECT_FALSE(pointer == nullptr);
    if (pointer)
    {
        EXPECT_EQ(*pointer, value);
    }
}

/*! \brief Checks that the key is not found */
void checkDoesNotFind(const gmx::HashedMap<char>& map, int key)
{
    const char* pointer = map.find(key);
    EXPECT_TRUE(pointer == nullptr);
}

TEST(HashedMap, InsertsFinds)
{
    gmx::HashedMap<char> map(2);

    map.insert(10, 'a');
    map.insert(5, 'b');
    map.insert(7, 'c');

    checkFinds(map, 10, 'a');
    checkFinds(map, 5, 'b');
    checkFinds(map, 7, 'c');
    checkDoesNotFind(map, 4);
}

TEST(HashedMap, NegativeKeysWork)
{
    gmx::HashedMap<char> map(5);

    map.insert(-1, 'a');
    map.insert(1, 'b');
    map.insert(-3, 'c');

    checkFinds(map, -1, 'a');
    checkFinds(map, 1, 'b');
    checkFinds(map, -3, 'c');
}

TEST(HashedMap, InsertsErases)
{
    gmx::HashedMap<char> map(3);

    map.insert(10, 'a');
    map.insert(5, 'b');
    map.insert(7, 'c');

    checkFinds(map, 10, 'a');
    map.erase(10);
    checkDoesNotFind(map, 10);
}

TEST(HashedMap, InsertsOrAssigns)
{
    gmx::HashedMap<char> map(3);

    map.insert(10, 'a');
    map.insert(5, 'b');

    map.insert_or_assign(7, 'c');
    checkFinds(map, 7, 'c');

    checkFinds(map, 10, 'a');
    map.insert_or_assign(10, 'd');
    checkFinds(map, 10, 'd');
}

TEST(HashedMap, Clears)
{
    gmx::HashedMap<char> map(3);

    map.insert(10, 'a');
    map.insert(5, 'b');
    map.insert(7, 'c');

    map.clear();
    checkDoesNotFind(map, 10);
    checkDoesNotFind(map, 5);
    checkDoesNotFind(map, 7);
}

// Check that entries with the same hash are handled correctly
TEST(HashedMap, LinkedEntries)
{
    // HashedMap uses bit masking, so keys that differ by exactly
    // a power of 2 larger than the table size will have the same hash

    gmx::HashedMap<char> map(20);

    const int largePowerOf2 = 2048;

    map.insert(3 + 0 * largePowerOf2, 'a');
    map.insert(3 + 1 * largePowerOf2, 'b');
    map.insert(3 + 2 * largePowerOf2, 'c');

    checkFinds(map, 3 + 0 * largePowerOf2, 'a');
    checkFinds(map, 3 + 1 * largePowerOf2, 'b');
    checkFinds(map, 3 + 2 * largePowerOf2, 'c');

    // Erase the middle entry in the linked list
    map.erase(3 + 1 * largePowerOf2);

    checkFinds(map, 3 + 0 * largePowerOf2, 'a');
    checkDoesNotFind(map, 3 + 1 * largePowerOf2);
    checkFinds(map, 3 + 2 * largePowerOf2, 'c');
}

// HashedMap only throws in debug mode, so only test in debug mode
#ifndef NDEBUG

TEST(HashedMap, CatchesDuplicateKey)
{
    gmx::HashedMap<char> map(15);

    map.insert(10, 'a');
    map.insert(5, 'b');
    EXPECT_THROW_GMX(map.insert(10, 'c'), gmx::InvalidInputError);
}

#endif // NDEBUG

// Check the table is resized after clear()
TEST(HashedMap, ResizesTable)
{
    gmx::HashedMap<char> map(1);

    // This test assumes the minimum bucket count is 64 or less
    EXPECT_LT(map.bucket_count(), 128);

    for (int i = 0; i < 60; i++)
    {
        map.insert(2 * i + 3, 'a');
    }
    EXPECT_LT(map.bucket_count(), 128);

    // Check that the table size is double #elements after clear()
    map.clearAndResizeHashTable();
    EXPECT_EQ(map.bucket_count(), 128);

    // Check that calling clear() a second time does not resize
    map.clearAndResizeHashTable();
    EXPECT_EQ(map.bucket_count(), 128);

    map.insert(2, 'b');
    EXPECT_EQ(map.bucket_count(), 128);

    // Check that calling clear with 1 elements sizes down
    map.clearAndResizeHashTable();
    EXPECT_LT(map.bucket_count(), 128);
}

} // namespace
} // namespace test
} // namespace gmx
