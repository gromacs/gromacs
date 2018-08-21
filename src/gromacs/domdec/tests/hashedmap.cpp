/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2018, by the GROMACS development team, led by
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
 * Tests for the HashedMap class.
 *
 * \author berk Hess <hess@kth.se>
 * \ingroup module_domdec
 */
#include "gmxpre.h"

#include "gromacs/domdec/hashedmap.h"

#include <gtest/gtest.h>

#include "testutils/testasserts.h"

namespace
{

static void checkFinds(const gmx::HashedMap<char> &map,
                       int                         key,
                       char                        value)
{
    const char *pointer = map.find(key);
    EXPECT_FALSE(pointer == nullptr);
    if (pointer)
    {
        EXPECT_EQ(*pointer, value);
    }
}

static void checkDoesNotFind(const gmx::HashedMap<char> &map,
                             int                         key)
{
    const char *pointer = map.find(key);
    EXPECT_TRUE(pointer == nullptr);
}

TEST(HashedMap, InsertsFinds)
{
    gmx::HashedMap<char> map(2);

    map.insert(10, 'a');
    map.insert(5,  'b');
    map.insert(7,  'c');

    checkFinds(map, 10, 'a');
    checkFinds(map, 5,  'b');
    checkFinds(map, 7,  'c');
    checkDoesNotFind(map, 4);
}

TEST(HashedMap, InsertsErases)
{
    gmx::HashedMap<char> map(3);

    map.insert(10, 'a');
    map.insert(5,  'b');
    map.insert(7,  'c');

    checkFinds(map, 10, 'a');
    map.erase(10);
    checkDoesNotFind(map, 10);
}

TEST(HashedMap, Clears)
{
    gmx::HashedMap<char> map(3);

    map.insert(10, 'a');
    map.insert(5,  'b');
    map.insert(7,  'c');

    map.clear();
    checkDoesNotFind(map, 10);
    checkDoesNotFind(map, 5);
    checkDoesNotFind(map, 7);
}

// HashedMap only throws in debug mode, so only test in debug mode
#ifndef NDEBUG

TEST(HashedMap, CatchesInvalidKey)
{
    gmx::HashedMap<char> map(101);

    EXPECT_THROW_GMX(map.insert(-1, 'a'), gmx::InvalidInputError);
}

TEST(HashedMap, CatchesDuplicateKey)
{
    gmx::HashedMap<char> map(15);

    map.insert(10, 'a');
    map.insert(5,  'b');
    EXPECT_THROW_GMX(map.insert(10, 'c'), gmx::InvalidInputError);
}

#endif // NDEBUG

}      // namespace
