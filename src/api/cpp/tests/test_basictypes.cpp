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
// #include "gmxapi/basictypes.h"

#include <array>

#include "gmxapi/basictypes.h"
#include <gtest/gtest.h>

#define real float

using namespace gmxapi;

template<typename Scalar>
struct Vector3D
{
    Scalar x;
    Scalar y;
    Scalar z;

    operator std::array<Scalar, 3>() { return std::array<Scalar, 3>{{x, y, z}}; }
};

class foo
{
    public:
        explicit foo(int in) : a(in){};
        static int get(const foo &f) {return f.a; };
        struct optionA{};
        struct optionB{};
    private:
        int a;
};
template<class T> bool check();
template<> bool check<typename foo::optionA>()
{
    return true;
};
template<> bool check<typename foo::optionB>()
{
    return false;
};

TEST(Eric, Concepts)
{
    ASSERT_EQ(check<foo::optionA>(), true);
    ASSERT_EQ(check<foo::optionB>(), false);
    foo bar {
        3
    };
    ASSERT_EQ(foo::get(bar), 3);
}

TEST(VectorTypes, ImplicitConversion)
{
    Vector3D<real> v;
    v.x = 1;
    v.y = 2;
    v.z = 3;
    auto a = std::array<real, 3>(v);
    ASSERT_EQ(a[0], 1);
    ASSERT_EQ(a[1], 2);
    ASSERT_EQ(a[2], 3);
}

TEST(Managed, BasicSemantics)
{
    // using ManagedScalar = Managed<float>;
    // // Default constructor acquires and initializes resouces
    // ManagedScalar s;
    // {
    // auto sh = s.handle<read_only>();
    // ASSERT_EQ(sh, 0);
    // }
    using ManagedVectorList = Managed < std::vector < float>>;
    ManagedVectorList r {{
                             1., 2., 3.
                         }};
    // Test initialization and access modes
    {
        // Test initialization and first read.
        auto rh = r.handle();
        ASSERT_EQ(rh[2], 3);
    }
    // using ManagedVectorList = Managed<std::vector<std::array<float, 3>>>;
    // ManagedVectorList r{{{1.,2.,3.},{0.,0.,0.},{4.,5.,6.}}};
    // // Test initialization and access modes
    // {
    //     // Test initialization and first read.
    //     auto rh = r.handle();
    //     ASSERT_EQ(rh[2][2], 6);
    // }
    // {
    //     // Test write and read from write handle.
    //     auto rh = r.handle<Access::LOCKED_WRITE>();
    //     rh[1][1] = float(2);
    //     ASSERT_EQ(rh[1][1], float(2));
    // }
    // {
    //     // Test that the write stuck.
    //     auto rh = r.handle<Access::LOCKED_READ>();
    //     ASSERT_EQ(rh[1][1], float(2));
    // }
    // {
    //     // Test that COPY_ON_WRITE does what it is supposed to (i.e. writes _don't_ stick)
    //     {
    //         auto rh = r.handle<Access::COPY_ON_WRITE>();
    //         rh[0][0] = float(42);
    //         ASSERT_EQ(rh[0][0], 42);
    //     }
    //     ASSERT_EQ(r.handle<Access::READ_ONLY>()[0][0], float(1));
    //     ASSERT_EQ(r.handle<Access::LOCKED_READ>()[0][0], float(1));
    //     ASSERT_EQ(r.handle<Access::LOCKED_WRITE>()[0][0], float(1));
    //     ASSERT_EQ(r.handle<Access::COPY_ON_WRITE>()[0][0], float(1));
    // }
}
