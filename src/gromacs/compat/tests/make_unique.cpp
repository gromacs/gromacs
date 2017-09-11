/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2017, by the GROMACS development team, led by
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
 * \brief Tests for gmx::compat::make_unique
 *
 * \author M. Eric Irrgang <ericirrgang@gmail.com>
 * \ingroup module_compat
 */
#include "gmxpre.h"

#include "gromacs/compat/make_unique.h"

#include <gtest/gtest.h>

namespace gmx
{

namespace
{

struct dummy
{
    char foo;
    char bar;
    dummy() :
        foo {0},
    bar {0}
    {};
    dummy(const char a, const char b) :
        foo {a},
    bar {b}
    {};
};

TEST(CompatibilityHelper, MakeUniqueCompiles)
{
    // Check template parameters
    auto ptr = gmx::compat::make_unique<dummy>();
    ASSERT_NE(ptr, nullptr);
    ASSERT_NE(ptr.get(), nullptr);
    constexpr bool is_dummy = std::is_same < decltype(ptr), std::unique_ptr < dummy>>::value;
    ASSERT_TRUE(is_dummy);

    // Check template and function parameters
    ptr = gmx::compat::make_unique<dummy>('a', 'b');
    ASSERT_EQ(ptr->foo, 'a');
}


} // anonymous namespace
} // namespace gmx
