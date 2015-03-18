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

/* Tests for gmx::ArrayRef and gmx::ConstArrayRef follow. There's a
   bunch of templating that supports testing const and non-const array
   references of a wide range of integral types that have been created
   through various mechanisms. Probably that's better than a forest of
   cut-and-paste.

   Ideally, we would not need both of ArrayRefTypeWrapper and
   ConstArrayRefTypeWrapper, but I haven't thought of a way that it
   can work in C++98, particularly since the free functions have
   different names for const and non-const case. */

//! Helper template so we can parameterize typed tests both on kind of ArrayRef and base type
template <typename T>
class ArrayRefTypeWrapper
{
    public:
        typedef ArrayRef<T> ArrayRefType;
        typedef T BaseType;
        static ArrayRefType makeFromPointers(BaseType *begin, BaseType *end);
        static ArrayRefType makeFromArray(BaseType *begin, size_t size);
        static ArrayRefType makeFromVector(typename std::vector<BaseType>::iterator begin,
                                           typename std::vector<BaseType>::iterator end);
};

template <typename T>
typename ArrayRefTypeWrapper<T>::ArrayRefType ArrayRefTypeWrapper<T>::makeFromPointers(BaseType * begin, BaseType *end)
{
    return gmx::arrayRefFromPointers<BaseType>(begin, end);
}

template <typename T>
typename ArrayRefTypeWrapper<T>::ArrayRefType ArrayRefTypeWrapper<T>::makeFromArray(BaseType * begin, size_t size)
{
    return gmx::arrayRefFromArray<BaseType>(begin, size);
}

template <typename T>
typename ArrayRefTypeWrapper<T>::ArrayRefType ArrayRefTypeWrapper<T>::makeFromVector(typename std::vector<BaseType>::iterator begin,
                                                                                     typename std::vector<BaseType>::iterator end)
{
    return gmx::arrayRefFromVector<BaseType>(begin, end);
}

//! Helper template so we can parameterize typed tests both on kind of ArrayRef and base type
template <typename T>
class ConstArrayRefTypeWrapper
{
    public:
        typedef ConstArrayRef<T> ArrayRefType;
        typedef T BaseType;
        static ArrayRefType makeFromPointers(BaseType *begin, BaseType *end);
        static ArrayRefType makeFromArray(BaseType *begin, size_t size);
        static ArrayRefType makeFromVector(typename std::vector<BaseType>::const_iterator begin,
                                           typename std::vector<BaseType>::const_iterator end);
};

template <typename T>
typename ConstArrayRefTypeWrapper<T>::ArrayRefType ConstArrayRefTypeWrapper<T>::makeFromPointers(BaseType * begin, BaseType *end)
{
    return gmx::constArrayRefFromPointers<BaseType>(begin, end);
}

template <typename T>
typename ConstArrayRefTypeWrapper<T>::ArrayRefType ConstArrayRefTypeWrapper<T>::makeFromArray(BaseType * begin, size_t size)
{
    return gmx::constArrayRefFromArray<BaseType>(begin, size);
}

template <typename T>
typename ConstArrayRefTypeWrapper<T>::ArrayRefType ConstArrayRefTypeWrapper<T>::makeFromVector(typename std::vector<BaseType>::const_iterator begin,
                                                                                               typename std::vector<BaseType>::const_iterator end)
{
    return gmx::constArrayRefFromVector<BaseType>(begin, end);
}

//! Permit all the tests to run on all kinds of ArrayRefs, that have been constructed differently
template <typename T>
class ArrayRefTest : public ::testing::Test
{
    public:
        typedef typename T::ArrayRefType ArrayRefType;
        typedef typename T::BaseType BaseType;

        /*! \brief Run the same tests all the time
         *
         * Note that test cases must call this->runTests(), because
         * that's how the derived class templates that implement
         * type-parameterized tests actually work. */
        void runTests(BaseType     *a,
                      size_t        aSize,
                      BaseType     *aData,
                      ArrayRefType &arrayRef)
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

TYPED_TEST_CASE_P(ArrayRefTest);

/* Welcome back to the past. While you can declare a static array[] of
   templated type in a class, in C++98, you have to define it outside
   the class, and when you do, the compiler knows the declaration is
   incomplete and can't match the types to actual functions. So,
   declaring locals is the only choice available, so we need macros to
   avoid duplication. Lovely.

   And while we're at it, some typedefs will make reading easier. */
#define DEFINE_ARRAY(a, aSize)                              \
    typedef typename TypeParam::BaseType BaseType;          \
    typedef typename TypeParam::ArrayRefType ArrayRefType;  \
    BaseType (a)[] = { static_cast<BaseType>(1.2),          \
                       static_cast<BaseType>(2.4),          \
                       static_cast<BaseType>(-3.1) \
    };       \
    size_t (aSize) = sizeof((a)) / sizeof(BaseType);

#define DEFINE_VECTOR_FROM_ARRAY(v, a, aSize)   \
    std::vector<BaseType>(v);                   \
    for (size_t i = 0; i != (aSize); ++i)       \
    {                                           \
        (v).push_back((a)[i]);                  \
    }

TYPED_TEST_P(ArrayRefTest, MakeWithAssignmentWorks)
{
    DEFINE_ARRAY(a, aSize);
    ArrayRefType arrayRef = a;
    this->runTests(a, aSize, a, arrayRef);
}

TYPED_TEST_P(ArrayRefTest, ConstructWithTemplateConstructorWorks)
{
    DEFINE_ARRAY(a, aSize);
    ArrayRefType arrayRef(a);
    this->runTests(a, aSize, a, arrayRef);
}

TYPED_TEST_P(ArrayRefTest, ConstructFromPointersWorks)
{
    DEFINE_ARRAY(a, aSize);
    ArrayRefType arrayRef(a, a + aSize);
    this->runTests(a, aSize, a, arrayRef);
}

TYPED_TEST_P(ArrayRefTest, MakeFromPointersWorks)
{
    DEFINE_ARRAY(a, aSize);
    ArrayRefType arrayRef = TypeParam::makeFromPointers(a, a + aSize);
    this->runTests(a, aSize, a, arrayRef);
}

TYPED_TEST_P(ArrayRefTest, MakeFromArrayWorks)
{
    DEFINE_ARRAY(a, aSize);
    ArrayRefType arrayRef = TypeParam::makeFromArray(a, aSize);
    this->runTests(a, aSize, a, arrayRef);
}

TYPED_TEST_P(ArrayRefTest, ConstructFromVectorWorks)
{
    DEFINE_ARRAY(a, aSize);
    DEFINE_VECTOR_FROM_ARRAY(v, a, aSize);
    ArrayRefType arrayRef(v);
    this->runTests(a, aSize, v.data(), arrayRef);
}

TYPED_TEST_P(ArrayRefTest, MakeFromVectorWorks)
{
    DEFINE_ARRAY(a, aSize);
    DEFINE_VECTOR_FROM_ARRAY(v, a, aSize);
    ArrayRefType arrayRef = TypeParam::makeFromVector(v.begin(), v.end());
    this->runTests(a, aSize, v.data(), arrayRef);
}

REGISTER_TYPED_TEST_CASE_P(ArrayRefTest,
                           MakeWithAssignmentWorks,
                           ConstructWithTemplateConstructorWorks,
                           ConstructFromPointersWorks,
                           MakeFromPointersWorks,
                           MakeFromArrayWorks,
                           ConstructFromVectorWorks,
                           MakeFromVectorWorks);

//! Define the kinds of types that end up being available as TypeParam in the test cases
typedef ::testing::Types<
        ArrayRefTypeWrapper<char>,
        ConstArrayRefTypeWrapper<char>,
        ArrayRefTypeWrapper<unsigned char>,
        ConstArrayRefTypeWrapper<unsigned char>,
        ArrayRefTypeWrapper<int>,
        ConstArrayRefTypeWrapper<unsigned int>,
        ArrayRefTypeWrapper<long>,
        ConstArrayRefTypeWrapper<unsigned long>,
        ArrayRefTypeWrapper<gmx_int64_t>,
        ConstArrayRefTypeWrapper<gmx_int64_t>,
        ArrayRefTypeWrapper<gmx_uint64_t>,
        ConstArrayRefTypeWrapper<gmx_uint64_t>,
        ArrayRefTypeWrapper<float>,
        ConstArrayRefTypeWrapper<float>,
        ArrayRefTypeWrapper<double>,
        ConstArrayRefTypeWrapper<double>
        > ArrayRefTypes;

INSTANTIATE_TYPED_TEST_CASE_P(My, ArrayRefTest, ArrayRefTypes);

} // namespace

} // namespace
