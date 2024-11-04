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
 * Tests for scope guard functionality.
 *
 * These tests check the functionality of scope_guard.h
 * This header has its own test suite which we trust, but we still want
 * to ensure that some functionality works the way we need it to.
 *
 * \author Petter Johansson <pettjoha@kth.se>
 * \ingroup module_utility
 */

#include "external/scope_guard/scope_guard.h"

#include <utility>

#include <gtest/gtest.h>

namespace gmx
{
namespace test
{
namespace
{

// \brief Function signature for a deleter.
using DeleterType = int (*)(int);

/*! \brief Global value which is set by the below function.
 *
 *  A global value is used here because we want to ensure that the scope
 *  guards execute a callback function of form int(*)(int) (similar to
 *  closing functions used by the HDF5 library, e.g. H5Gclose) when going
 *  out of scope.
 *
 *  We need to inspect the value after the guard exits and the simplest way is to
 *  use a global variable.
 */
static int s_deleterValue = 0;

/* \brief Callback function which sets a value to s_deleterValue.
 *
 * \param[in] object Value to set.
 */
int deleterSetter(int object)
{
    s_deleterValue = object;

    return 0;
}

/* \brief Class which binds a deleter function and can be returned as a concrete type.
 */
class BoundDeleter
{
public:
    BoundDeleter(int object, DeleterType deleter) : object_{ object }, deleter_{ deleter } {}
    void operator()() { deleter_(object_); }

private:
    int         object_;
    DeleterType deleter_;
};

/* \brief Create a scope guard for a given object, using a given function.
 *
 * \param[in] object Object to guard with the function
 * \param[in] deleter Deleting function, using object as an argument
 * \returns An (object, guard) pair.
 */
std::pair<int, sg::detail::scope_guard<BoundDeleter>> makeGuard(int object, DeleterType deleter)
{
    return std::make_pair(object, sg::make_scope_guard(BoundDeleter(object, deleter)));
}

TEST(ScopeGuardTest, ScopeGuardExecutesCallbackOnExit)
{
    int i = 1;

    {
        const auto guard = sg::make_scope_guard([&i]() { i++; });
        EXPECT_EQ(i, 1);
    }

    EXPECT_EQ(i, 2);
}

TEST(ScopeGuardTest, ScopeGuardCanFreePointers)
{
    constexpr int value = 5;
    int*          ptr   = new int{ value };

    {
        auto guard = sg::make_scope_guard(
                [&ptr]()
                {
                    delete ptr;
                    ptr = nullptr;
                });

        EXPECT_EQ(*ptr, value);
    }

    EXPECT_EQ(ptr, nullptr);
}

TEST(ScopeGuardTest, ScopeGuardsCanBeCreatedByHelperFunctions)
{
    constexpr int originalValue = 0;
    constexpr int setValue      = 5;

    s_deleterValue = originalValue;

    {
        const auto [object, guard] = makeGuard(setValue, deleterSetter);

        EXPECT_EQ(object, setValue);
        EXPECT_EQ(s_deleterValue, originalValue);
    }

    EXPECT_EQ(s_deleterValue, setValue);
}

} // namespace
} // namespace test
} // namespace gmx
