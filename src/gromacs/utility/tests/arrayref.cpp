/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2015, by the GROMACS development team, led by
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
 * \brief Tests for gmx::ArrayRef and gmx::ConstArrayRef.
 *
 * \author Mark Abraham <mark.j.abraham@gmail.com>
 * \ingroup module_utility
 */
#include "gmxpre.h"

#include "config.h"

#include "gromacs/utility/arrayref.h"

#include <vector>

#include <gtest/gtest.h>

#include "gromacs/utility/basedefinitions.h"
#include "gromacs/utility/real.h"

namespace gmx
{

namespace
{

TEST(EmptyArrayRefTest, IsEmpty)
{
    EmptyArrayRef  emptyArray = {};
    ArrayRef<real> empty(emptyArray);
    size_t         zero = 0; // stop warnings about int vs unsigned comparison

    EXPECT_EQ(zero, empty.size());
    EXPECT_TRUE(empty.empty());
}

TEST(EmptyConstArrayRefTest, IsEmpty)
{
    EmptyArrayRef       emptyArray = {};
    ConstArrayRef<real> empty(emptyArray);
    size_t              zero = 0; // stop warnings about int vs unsigned comparison

    EXPECT_EQ(zero, empty.size());
    EXPECT_TRUE(empty.empty());
}

//! Define the types that end up being available as TypeParam in the test cases for both kinds of ArrayRef
typedef ::testing::Types<
        char,
        unsigned char,
        int,
        unsigned int,
        long,
        unsigned long,
        gmx_int64_t,
        gmx_uint64_t,
        float,
        double
        > ArrayRefTypes;

// === Tests for gmx::ConstArrayRef

//! Permit all the tests to run on all kinds of ArrayRefs, that have been constructed differently
template <typename TypeParam>
class ConstArrayRefTest : public ::testing::Test
{
    public:
        typedef ConstArrayRef<TypeParam> ArrayRefType;

        /*! \brief Run the same tests all the time
         *
         * Note that test cases must call this->runTests(), because
         * that's how the derived class templates that implement
         * type-parameterized tests actually work. */
        void runTests(TypeParam     *a,
                      size_t         aSize,
                      TypeParam     *aData,
                      ArrayRefType  &arrayRef)
        {
            EXPECT_EQ(aSize, arrayRef.size());
            EXPECT_FALSE(arrayRef.empty());
            EXPECT_EQ(aData, arrayRef.data());
            EXPECT_EQ(a[0], arrayRef.front());
            EXPECT_EQ(a[aSize-1], arrayRef.back());
            for (size_t i = 0; i != aSize; ++i)
            {
                EXPECT_EQ(a[i], arrayRef[i]);
            }
        }
};

TYPED_TEST_CASE_P(ConstArrayRefTest);

/* Welcome back to the past. While you can declare a static array[] of
   templated type in a class, in C++98, you have to define it outside
   the class, and when you do, the compiler knows the declaration is
   incomplete and can't match the types to actual functions. So,
   declaring locals is the only choice available, so we need macros to
   avoid duplication. Lovely. */
#define DEFINE_ARRAY(a, aSize)                          \
    TypeParam (a)[] = { static_cast<TypeParam>(1.2),    \
                        static_cast<TypeParam>(2.4),    \
                        static_cast<TypeParam>(-3.1)    \
    };                                                  \
    size_t (aSize) = sizeof((a)) / sizeof(TypeParam);

#define DEFINE_VECTOR_FROM_ARRAY(v, a, aSize)   \
    std::vector<TypeParam>(v);                  \
    for (size_t i = 0; i != (aSize); ++i)       \
    {                                           \
        (v).push_back((a)[i]);                  \
    }

TYPED_TEST_P(ConstArrayRefTest, MakeWithAssignmentWorks)
{
    DEFINE_ARRAY(a, aSize);
    ConstArrayRef<TypeParam> arrayRef = a;
    this->runTests(a, aSize, a, arrayRef);
}

TYPED_TEST_P(ConstArrayRefTest, ConstructWithTemplateConstructorWorks)
{
    DEFINE_ARRAY(a, aSize);
    ConstArrayRef<TypeParam> arrayRef(a);
    this->runTests(a, aSize, a, arrayRef);
}

TYPED_TEST_P(ConstArrayRefTest, ConstructFromPointersWorks)
{
    DEFINE_ARRAY(a, aSize);
    ConstArrayRef<TypeParam> arrayRef(a, a + aSize);
    this->runTests(a, aSize, a, arrayRef);
}

TYPED_TEST_P(ConstArrayRefTest, MakeFromPointersWorks)
{
    DEFINE_ARRAY(a, aSize);
    ConstArrayRef<TypeParam> arrayRef = gmx::constArrayRefFromPointers<TypeParam>(a, a + aSize);;
    this->runTests(a, aSize, a, arrayRef);
}

TYPED_TEST_P(ConstArrayRefTest, MakeFromArrayWorks)
{
    DEFINE_ARRAY(a, aSize);
    ConstArrayRef<TypeParam> arrayRef = gmx::constArrayRefFromArray<TypeParam>(a, aSize);
    this->runTests(a, aSize, a, arrayRef);
}

TYPED_TEST_P(ConstArrayRefTest, ConstructFromVectorWorks)
{
    DEFINE_ARRAY(a, aSize);
    DEFINE_VECTOR_FROM_ARRAY(v, a, aSize);
    ConstArrayRef<TypeParam> arrayRef(v);
    this->runTests(a, aSize, v.data(), arrayRef);
}

TYPED_TEST_P(ConstArrayRefTest, MakeFromVectorWorks)
{
    DEFINE_ARRAY(a, aSize);
    DEFINE_VECTOR_FROM_ARRAY(v, a, aSize);
    ConstArrayRef<TypeParam> arrayRef = gmx::constArrayRefFromVector<TypeParam>(v.begin(), v.end());
    this->runTests(a, aSize, v.data(), arrayRef);
}

REGISTER_TYPED_TEST_CASE_P(ConstArrayRefTest,
                           MakeWithAssignmentWorks,
                           ConstructWithTemplateConstructorWorks,
                           ConstructFromPointersWorks,
                           MakeFromPointersWorks,
                           MakeFromArrayWorks,
                           ConstructFromVectorWorks,
                           MakeFromVectorWorks);

INSTANTIATE_TYPED_TEST_CASE_P(My, ConstArrayRefTest, ArrayRefTypes);

/* === Tests for gmx::ArrayRef

   These follow the same approach as used for gmx::ConstArrayRef above,
   but there's no good way to avoid repetition (in C++98, at least).

   There is also machinery to test the usage in mdrun
   signalling. There, we take a non-const struct-field array of static
   length and make an ArrayRef to it using the template constructor
   that is supposed to infer the length from the static size. But on
   xlc on BlueGene/Q, if the base type is not char (or unsigned char),
   the generated code ends up with an ArrayRef of zero size, so
   everything breaks. Presumably the default code path accidentally
   works for char.

   Fortunately, all current uses of that constructor have a base type
   of char, so there's no big problem. Using a pointer-based
   constructor does work, if there's ever a problem (and that is
   tested above). */

//! Permit all the tests to run on all kinds of ArrayRefs, that have been constructed differently
template <typename TypeParam>
class NonConstArrayRefTest : public ::testing::Test
{
    public:
        typedef ArrayRef<TypeParam> ArrayRefType;

        /*! \brief Run the same tests all the time
         *
         * Note that test cases must call this->runTests(), because
         * that's how the derived class templates that implement
         * type-parameterized tests actually work. */
        void runTests(TypeParam           *a,
                      size_t               aSize,
                      TypeParam           *aData,
                      ArrayRef<TypeParam> &arrayRef,
                      bool                 bTestingProblematicConstructor)
        {
#if defined(GMX_TARGET_BGQ) && defined(__xlC__)
            if (bTestingProblematicConstructor &&
                sizeof(TypeParam) != sizeof(char))
            {
                /* This is the buggy case, failure is expected */
                EXPECT_NE(aSize, arrayRef.size());
                /* Because back() refers to the wrong memory, it might
                   sometimes have the right value, so we cannot assert
                   on that value! */
                //EXPECT_NE(a[aSize-1], arrayRef.back());
            }
            else
#endif
            {
                EXPECT_EQ(aSize, arrayRef.size());
                EXPECT_EQ(a[aSize-1], arrayRef.back());
            }
            EXPECT_FALSE(arrayRef.empty());
            EXPECT_EQ(aData, arrayRef.data());
            EXPECT_EQ(a[0], arrayRef.front());
            for (size_t i = 0; i != aSize; ++i)
            {
                EXPECT_EQ(a[i], arrayRef[i]);
            }
        }
};

TYPED_TEST_CASE_P(NonConstArrayRefTest);

//! Helper struct for the case actually used in mdrun signalling
template <typename T>
struct Helper
{
    public:
        T   a[3];
        int size;
};

//! Test of the case actually used in mdrun signalling
TYPED_TEST_P(NonConstArrayRefTest, ConstructFromStructFieldWithTemplateConstructorWorks)
{
    DEFINE_ARRAY(a, aSize);
    Helper<TypeParam> h;
    h.size = aSize;
    for (size_t i = 0; i != h.size; ++i)
    {
        h.a[i] = a[i];
    }
    ArrayRef<TypeParam> arrayRef(h.a);
    this->runTests(h.a, h.size, h.a, arrayRef, true);
}

// Subsequent test cases are the same as for gmx::ConstArrayRef

TYPED_TEST_P(NonConstArrayRefTest, MakeWithAssignmentWorks)
{
    DEFINE_ARRAY(a, aSize);
    ArrayRef<TypeParam> arrayRef = a;
    this->runTests(a, aSize, a, arrayRef, false);
}

TYPED_TEST_P(NonConstArrayRefTest, ConstructWithTemplateConstructorWorks)
{
    DEFINE_ARRAY(a, aSize);
    ArrayRef<TypeParam> arrayRef(a);
    this->runTests(a, aSize, a, arrayRef, false);
}

TYPED_TEST_P(NonConstArrayRefTest, ConstructFromPointersWorks)
{
    DEFINE_ARRAY(a, aSize);
    ArrayRef<TypeParam> arrayRef(a, a + aSize);
    this->runTests(a, aSize, a, arrayRef, false);
}

TYPED_TEST_P(NonConstArrayRefTest, MakeFromPointersWorks)
{
    DEFINE_ARRAY(a, aSize);
    ArrayRef<TypeParam> arrayRef = gmx::arrayRefFromPointers<TypeParam>(a, a + aSize);;
    this->runTests(a, aSize, a, arrayRef, false);
}

TYPED_TEST_P(NonConstArrayRefTest, MakeFromArrayWorks)
{
    DEFINE_ARRAY(a, aSize);
    ArrayRef<TypeParam> arrayRef = gmx::arrayRefFromArray<TypeParam>(a, aSize);
    this->runTests(a, aSize, a, arrayRef, false);
}

TYPED_TEST_P(NonConstArrayRefTest, ConstructFromVectorWorks)
{
    DEFINE_ARRAY(a, aSize);
    DEFINE_VECTOR_FROM_ARRAY(v, a, aSize);
    ArrayRef<TypeParam> arrayRef(v);
    this->runTests(a, aSize, v.data(), arrayRef, false);
}

TYPED_TEST_P(NonConstArrayRefTest, MakeFromVectorWorks)
{
    DEFINE_ARRAY(a, aSize);
    DEFINE_VECTOR_FROM_ARRAY(v, a, aSize);
    ArrayRef<TypeParam> arrayRef = gmx::arrayRefFromVector<TypeParam>(v.begin(), v.end());
    this->runTests(a, aSize, v.data(), arrayRef, false);
}

REGISTER_TYPED_TEST_CASE_P(NonConstArrayRefTest,
                           MakeWithAssignmentWorks,
                           ConstructWithTemplateConstructorWorks,
                           ConstructFromPointersWorks,
                           MakeFromPointersWorks,
                           MakeFromArrayWorks,
                           ConstructFromVectorWorks,
                           MakeFromVectorWorks,
                           ConstructFromStructFieldWithTemplateConstructorWorks);

INSTANTIATE_TYPED_TEST_CASE_P(My, NonConstArrayRefTest, ArrayRefTypes);

} // namespace

} // namespace
