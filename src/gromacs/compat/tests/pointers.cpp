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
 * \brief Tests for pointers.h, e.g. gmx::compat::not_null
 *
 * \author Mark Abraham <mark.j.abraham@gmail.com>
 * \ingroup module_compat
 */
#include "gmxpre.h"

#include "gromacs/compat/pointers.h"

#include <memory>
#include <string>
#include <vector>

#include <gtest/gtest.h>

#include "gromacs/utility/basedefinitions.h"

#include "testutils/testasserts.h"

namespace gmx
{
namespace compat
{
namespace
{

TEST(NotNullConstruction, Works)
{
    // shared_ptr<int> is nullptr assignable
    not_null<std::shared_ptr<int>> sharedPointer(std::make_shared<int>(10));

#ifndef NDEBUG
    int* nullPointer = nullptr;
    GMX_EXPECT_DEATH_IF_SUPPORTED(not_null<int*> invalidNullPointer(nullPointer), "");
#endif

    int  value        = 20;
    int* validPointer = &value;
    {
        not_null<int*> validNotNullPointer(validPointer);
        GMX_UNUSED_VALUE(validNotNullPointer);
    }
    {
        not_null<int*> validNotNullPointer = not_null<int*>(validPointer);
        GMX_UNUSED_VALUE(validNotNullPointer);
    }
}

TEST(NotNullCasting, Works)
{
    struct MyBase
    {
    };
    struct MyDerived : public MyBase
    {
    };
    struct Unrelated
    {
    };

    MyBase    base;
    MyDerived derived;
    Unrelated unrelated;

    not_null<Unrelated*> u{ &unrelated };
    (void)u;
    not_null<MyDerived*> p{ &derived };
    not_null<MyBase*>    q(&base);
    // Allowed with heterogeneous copy constructor
    q = p;

    not_null<Unrelated*> t(reinterpret_cast<Unrelated*>(p.get()));
    EXPECT_EQ(reinterpret_cast<void*>(p.get()), reinterpret_cast<void*>(t.get()));
}

TEST(NotNullAssignment, Works)
{
    int            i = 12;
    not_null<int*> p(&i);
    EXPECT_EQ(*p, 12);
}

TEST(MakeNotNull, Works)
{
    {
        int i = 42;

        const not_null<int*> x = make_not_null(&i);
        EXPECT_EQ(*x, 42);
        not_null<int*> y = make_not_null(&i);
        EXPECT_EQ(*y, 42);
        not_null<const int*> z = make_not_null(&i);
        EXPECT_EQ(*z, 42);
    }

    {
        // TODO These should work, but the GSL version of
        // make_not_null has an auto return type that we can't use
        // here, so maybe the issue is there.
        /*
           int i = 42;
           int* p = &i;

           not_null<int *> x = make_not_null(p);
           EXPECT_EQ(*x, 42);
           not_null<int *> y = make_not_null(p);
           EXPECT_EQ(*y, 42);
           not_null<const int *> z = make_not_null(p);
           EXPECT_EQ(*z, 42);
         */
    }

    {
        std::unique_ptr<int> i = std::make_unique<int>(42);

        const not_null<int*> x = make_not_null(i);
        EXPECT_EQ(*x, 42);
        not_null<int*> y = make_not_null(i);
        EXPECT_EQ(*y, 42);
        not_null<const int*> z = make_not_null(i);
        EXPECT_EQ(*z, 42);
    }

    {
        std::unique_ptr<const int> i = std::make_unique<int>(42);

        // not_null<int *> does not compile, as expected
        not_null<const int*> z = make_not_null(i);
        EXPECT_EQ(*z, 42);
    }
}

TEST(NotNull, WorksInContainers)
{
    int            i = 12;
    not_null<int*> p(&i);

    std::vector<not_null<int*>> v;
    v.push_back(p);
    EXPECT_EQ(*v.back(), 12);
}

// TODO We currently don't have infrastructure for checking that e.g.
// expected static assertions fire and calls to deleted functions do
// not compile. When we do, there are more tests that should be found
// here.

} // anonymous namespace
} // namespace compat
} // namespace gmx
