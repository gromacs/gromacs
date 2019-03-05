/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2019, by the GROMACS development team, led by
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
/*! \file
    \ingroup module_utility
    \inpublicapi

    \brief declares gmx::IrregArray4D, a 4D-array that can have irregular shape

    \author R. Thomas Ullmann <tullman@gwdg.de>
 */
#ifndef GMX_UTILITY_DATASTRUCTURES_IRREG_ARRAY_4D_H
#define GMX_UTILITY_DATASTRUCTURES_IRREG_ARRAY_4D_H

#include "gromacs/utility/data_structures/irreg_array_3d.h"

namespace gmx
{


/*! \class IrregArray4D

    \brief irregularly shaped 4-dimensional array

     This array can have nonuniform second, third and fourth dimensions.
     That is, instead of constant dimensions M x N x O x P, this array
     has dimensions M x N(M) x O(M, N) x P(M, N, O). isIrreg() returns
     true if the array has irregular shape and false for regular shape.

     example application:
     interaction energy matrix indexed by site1, form1, site2, form2

    \ingroup module_utility
    \inpublicapi

    \tparam   T           data type to be stored
    \tparam   Allocator   allocator to be used in creating arrays storing T
                          allocators for storing additional book-keeping data
                          are derived from this allocator via rebind if necessary
    \tparam   isFlat      underlying array is 1-dimensional if true and multi-dimensional if false
 */
template<class T = size_t, class Allocator = std::allocator<T>, bool IsFlat = false>
class IrregArray4D
{
    public:
        //! the data type stored in the array
        typedef T value_type;
        //! the allocator type to allocate value_type***
        typedef typename std::allocator_traits<Allocator>::template rebind_alloc<value_type***> p3_allocator_type;
        //! the allocator type to allocate value_type**
        typedef typename std::allocator_traits<Allocator>::template rebind_alloc<value_type**> p2_allocator_type;
        //! the allocator type to allocate value_type*
        typedef typename std::allocator_traits<Allocator>::template rebind_alloc<value_type*> p_allocator_type;
        //! the allocator type to allocate value_type
        typedef typename std::allocator_traits<Allocator>::template rebind_alloc<value_type> allocator_type;
        //! the pointer type returned by the allocator
        typedef typename std::allocator_traits<Allocator>::pointer allocator_return_ptr_type;
        //! the (integer) datatype used to represent sizes, e.g., array lengths
        typedef  size_t  size_type;
        //! the (integer) data type used to represent array indices
        typedef ssize_t index_type;
        //! the allocator type to allocate size_type
        typedef typename std::allocator_traits<Allocator>::template rebind_alloc<size_type> allocator_size_type;
        //! type of array used for storing the lengths of array stripes
        typedef std::vector<size_type, allocator_size_type>  size_1d_array_type;
        //! type of array used for storing the lengths of array stripes
        typedef IrregArray2D<size_type, allocator_size_type, IsFlat> size_2d_array_type;
        //! type of array used for storing the lengths of array stripes
        typedef IrregArray3D<size_type, allocator_size_type, IsFlat> size_3d_array_type;
    private:
        //! array length in dimension 1
        size_type               length1_   = 0;
        //! array length in dimension 2
        size_type               length2_   = 0;
        //! array length in dimension 3
        size_type               length3_   = 0;
        //! array length in dimension 4
        size_type               length4_   = 0;
        //! total number of array elements
        size_type               nelements_ = 0;
        //! dimension1,2,3-specific number of elements of dimension 4
        IrregArray3D<size_type, allocator_size_type, IsFlat> irreg_;
        //! flag says that the array length in dimensions 2,3,4 is not a constant
        bool                    isIrreg_   = false;
        //! array storing the the actual data, 1D if IsFlat, 4D if !isFlat
        typename std::conditional<IsFlat, value_type*, value_type****>::type arr_ = nullptr;
        //! the allocator_ used to allocate arr_
        p3_allocator_type       p3_allocator_;
        //! the allocator_ used to allocate arr_[i]
        p2_allocator_type       p2_allocator_;
        //! the allocator_ used to allocate arr_[i][j]
        p_allocator_type        p_allocator_;
        //! the allocator_ used to allocate arr_[i][j][k]
        allocator_type          allocator_;
    public:
        /*! \brief default constructor without memory allocation

            The initArray member functions allow for memory allocation
            and content initialization.
         */
        IrregArray4D()
            : length1_(0), length2_(0), length3_(0), length4_(0), nelements_(0), isIrreg_(false), arr_(nullptr)
        {
// The properties of std::initializer list are tracked incorrectly by the
// clang static analyzer. In effect, the array allocation is not registered
// and a nullptr dereference reported when the array is accessed.
// https://bugs.llvm.org/show_bug.cgi?id=39042
// Scan-build reports a false positive nullptr dereference upon array access because
// it misses calls to template member functions within a constructor, but registers
// the array initialization to nullptr in the standard constructor that precedes allocation.
// https://bugs.llvm.org/show_bug.cgi?id=39028
// Trick the clang_static_analyzer by a dummy call to an init function that it understands
#ifdef __clang_analyzer__
            initArray(1, 1, 1, 1);
#endif      // end scan-build bug work-around
        }

        /*! \brief  copy constructor

            \param[in]   ini   array to be copied
         */
        IrregArray4D(const IrregArray4D &ini)
            : IrregArray4D()
        {
            initArray(ini);
        }

        /*! \brief  copy constructor

            \tparam      AForeign   allocator used by the input data structure
            \tparam      FForeign   true if the input data structure uses a flat data array

            \param[in]   ini   array to be copied
         */
        template<class AForeign = allocator_type, bool FForeign>
        IrregArray4D(const IrregArray4D<value_type, AForeign, FForeign> &ini)
            : IrregArray4D()
        {
            initArray(ini);
        }

        //! move constructor
        //! \param[in]   ini   source to be moved
        IrregArray4D(IrregArray4D &&ini) noexcept
            : length1_(ini.length1_), length2_(ini.length2_), length3_(ini.length3_), length4_(ini.length4_),
              nelements_(ini.nelements_), irreg_(std::move(ini.irreg_)), isIrreg_(ini.isIrreg_), arr_(ini.arr_),
              p3_allocator_(std::move(ini.p3_allocator_)), p2_allocator_(std::move(ini.p2_allocator_)),
              p_allocator_(std::move(ini.p_allocator_)), allocator_(std::move(ini.allocator_))
        {
            // set ini to a valid default state
            ini.length1_   = 0;
            ini.length2_   = 0;
            ini.length3_   = 0;
            ini.length4_   = 0;
            ini.nelements_ = 0;
            ini.isIrreg_   = false;
            ini.arr_       = nullptr;
        }

        /*! \brief  move assignment operator this = rhs
            \param[in]   rhs   right-hand side of the assignment
         */
        IrregArray4D &operator=(IrregArray4D &&rhs) noexcept
        {
            if (this != &rhs)
            {
                // free any previously allocated memory
                deallocateArray();
                // move the data from rhs
                length1_      = rhs.length1_;
                length2_      = rhs.length2_;
                length3_      = rhs.length3_;
                length4_      = rhs.length4_;
                nelements_    = rhs.nelements_;
                isIrreg_      = rhs.isIrreg_;
                arr_          = rhs.arr_;
                irreg_        = std::move(rhs.irreg_);
                allocator_    = std::move(rhs.allocator_);
                p_allocator_  = std::move(rhs.p_allocator_);
                p2_allocator_ = std::move(rhs.p2_allocator_);
                p3_allocator_ = std::move(rhs.p3_allocator_);
                // set rhs to a valid default state
                rhs.length1_   = 0;
                rhs.length2_   = 0;
                rhs.length3_   = 0;
                rhs.length4_   = 0;
                rhs.nelements_ = 0;
                rhs.isIrreg_   = false;
                rhs.arr_       = nullptr;
            }
            return *this;
        }

        /*! \brief   assignment operator *this = rhs
            \param[in]   rhs   right-hand side object of the assignment
         */
        IrregArray4D &operator=(const IrregArray4D &rhs)
        {
            if (this != &rhs)
            {
                deallocateArray();
                initArray(rhs);
            }
            return *this;
        }

        /*! \brief   assignment operator *this = rhs

            \tparam      AForeign   allocator used by the input data structure
            \tparam      FForeign   true if the input data structure uses a flat data array

            \param[in]   rhs   right-hand side object of the assignment
         */
        template<class AForeign = allocator_size_type, bool FForeign>
        IrregArray4D &operator=(const IrregArray4D<value_type, AForeign, FForeign> &rhs)
        {
            if (this != &rhs)
            {
                deallocateArray();
                initArray(rhs);
            }
            return *this;
        }

        /*! \brief  conversion operator (explicit to avoid unintended implicit conversion)

            usable via

            \code
                IrregArray4D<TForeign, AForeign> newArr =
                    static_cast<IrregArray4D<TForeign, AForeign> >
                    (
                        IrregArray4D<TCurrent, ACurrent> oldArr
                    );
            \endcode

            \tparam      TForeign   data type stored by the input data structure
            \tparam      AForeign   allocator used by the input data structure
            \tparam      FForeign   true if the output data structure uses a flat data array
         */
        template<typename TForeign, class AForeign = std::allocator<TForeign>, bool FForeign>
        explicit
        operator IrregArray4D<TForeign, AForeign, FForeign>() const
        {
            IrregArray4D<TForeign, AForeign, FForeign> result;
            result.initArray(*this);
            return result;
        }

        /*! \brief   constructor for a regularly shaped array

            \param[in]    l1      array length in dimension 1
            \param[in]    l2      array length in dimension 2
            \param[in]    l3      array length in dimension 3
            \param[in]    l4      array length in dimension 4
         */
        IrregArray4D(size_type l1, size_type l2, size_type l3, size_type l4)
            : IrregArray4D()
        {
            initArray(l1, l2, l3, l4);
        }
        /*! \brief   constructor for a regularly shaped array with initialization

            \param[in]    l1      array length in dimension 1
            \param[in]    l2      array length in dimension 2
            \param[in]    l3      array length in dimension 3
            \param[in]    l4      array length in dimension 4
            \param[in]    ini     initialization value for the data array elements
         */
        IrregArray4D(size_type l1, size_type l2, size_type l3, size_type l4, const value_type &ini)
            : IrregArray4D()
        {
            initArray(l1, l2, l3, l4, ini);
        }

        /*! \brief   general constructor helper function for irregularly shaped matrices and all special cases

            \tparam      AForeign   allocator used by the input data structure
            \tparam      FForeign   true if the input data structure uses a flat data array

            \param[in]    l4      number of elements in dimension 4 as function of the indices in the lower dimensions
                                  number of elements in dimensions 1,2,3 also deducible from l4
         */
        template <class AForeign = allocator_size_type, bool FForeign>
        IrregArray4D(const IrregArray3D<size_type, AForeign, FForeign> &l4)
            : IrregArray4D()
        {
            initArray(l4);
        }
        /*! \brief   general constructor helper function for irregularly shaped matrices and all special cases

            \tparam      AForeign   allocator used by the input data structure
            \tparam      FForeign   true if the input data structure uses a flat data array

            \param[in]    l4       number of elements in dimension 4 as function of the indices in the lower dimensions
                                   number of elements in dimensions 1,2,3 also deducible from l4
            \param[in]    ini      initial value for the array elements
         */
        template <class AForeign = allocator_size_type, bool FForeign>
        IrregArray4D(const IrregArray3D<size_type, AForeign, FForeign> &l4, const value_type &ini)
            : IrregArray4D()
        {
            initArray(l4, ini);
        }

        /*! \brief  initialize from a std::initializer list
                    \code
                        arr = {
                                  {
                                      {
                                          {value0000, value0001, ...},
                                          {value0010, value0011, ...}
                                          ...
                                      },
                                      {
                                          {value0100, value0101, ...},
                                          {value0110, value0111, ...}
                                          ...
                                      }
                                      ...
                                  },
                                  {
                                      {
                                          {value1000, value1001, ...},
                                          {value1010, value1011, ...}
                                          ...
                                      },
                                      {
                                          {value1100, value0101, ...},
                                          {value1110, value1111, ...}
                                          ...
                                      }
                                      ...
                                  }
                                  ...
                              };
                    \endcode

            \param[in]   ini  nested list of initialization values
         */
        IrregArray4D(const std::initializer_list<std::initializer_list<std::initializer_list<std::initializer_list<value_type> > > > &ini)
            : IrregArray4D()
        {
            typedef std::initializer_list<value_type>         l1_ilist_type;
            typedef std::initializer_list<l1_ilist_type>      l2_ilist_type;
            typedef std::initializer_list<l2_ilist_type>      l3_ilist_type;
            typedef const value_type*                     l1_ilist_ptr_type;
            typedef const l1_ilist_type*                  l2_ilist_ptr_type;
            typedef const l2_ilist_type*                  l3_ilist_ptr_type;
            typedef const l3_ilist_type*                  l4_ilist_ptr_type;

            // check whether the list contains empty elements in any but the innermost nesting level
            if (ini.size() == 0)
            {
                return;
            }
            else
            {
                for (l4_ilist_ptr_type lPtr = ini.begin(); lPtr != ini.end(); ++lPtr)
                {
                    if (lPtr->size() == 0)
                    {
                        return;
                    }
                }

                for (l4_ilist_ptr_type lPtr = ini.begin(); lPtr != ini.end(); ++lPtr)
                {
                    for (l3_ilist_ptr_type kPtr = lPtr->begin(); kPtr != lPtr->end(); ++kPtr)
                    {
                        if (kPtr->size() == 0)
                        {
                            return;
                        }
                    }
                }

                for (l4_ilist_ptr_type lPtr = ini.begin(); lPtr != ini.end(); ++lPtr)
                {
                    for (l3_ilist_ptr_type kPtr = lPtr->begin(); kPtr != lPtr->end(); ++kPtr)
                    {
                        for (l2_ilist_ptr_type jPtr = kPtr->begin(); jPtr != kPtr->end(); ++jPtr)
                        {
                            if (jPtr->size() == 0)
                            {
                                return;
                            }
                        }
                    }
                }
            }

            // check whether the array has a simple n x m x k x l shape or an irregular shape, where m, k or l differ between array stripes
            const size_type iniLength2 = ini.begin()->size();
            for (l4_ilist_ptr_type iPtrL4 = ini.begin(); iPtrL4 != ini.end(); ++iPtrL4)
            {
                if (iniLength2 != iPtrL4->size())
                {
                    isIrreg_ = true;
                    break;
                }
            }
            if (!isIrreg_)
            {
                const size_type iniLength3 = ini.begin()->begin()->size();
                for (l4_ilist_ptr_type iPtrL4 = ini.begin(); iPtrL4 != ini.end(); ++iPtrL4)
                {
                    for (l3_ilist_ptr_type iPtrL3 = iPtrL4->begin(); iPtrL3 != iPtrL4->end(); ++iPtrL3)
                    {
                        if (iniLength3 != iPtrL3->size())
                        {
                            isIrreg_ = true;
                            break;
                        }
                    }
                    if (isIrreg_)
                    {
                        break;
                    }
                }
            }
            if (!isIrreg_)
            {
                const size_type iniLength4 = ini.begin()->begin()->begin()->size();
                for (l4_ilist_ptr_type iPtrL4 = ini.begin(); iPtrL4 != ini.end(); ++iPtrL4)
                {
                    for (l3_ilist_ptr_type iPtrL3 = iPtrL4->begin(); iPtrL3 != iPtrL4->end(); ++iPtrL3)
                    {
                        for (l2_ilist_ptr_type iPtrL2 = iPtrL3->begin(); iPtrL2 != iPtrL3->end(); ++iPtrL2)
                        {
                            if (iniLength4 != iPtrL2->size())
                            {
                                isIrreg_ = true;
                                break;
                            }
                        }
                        if (isIrreg_)
                        {
                            break;
                        }
                    }
                    if (isIrreg_)
                    {
                        break;
                    }
                }
            }

            // allocate memory depending on the array shape
            if (isIrreg_)
            {
                size_1d_array_type l2(ini.size(), 0);
                size_t             i = 0;
                for (l4_ilist_ptr_type iPtrL4 = ini.begin(); iPtrL4 != ini.end(); ++iPtrL4)
                {
                    l2[i] = iPtrL4->size();
                    i++;
                }

                size_2d_array_type l3(l2, 0);
                i = 0;
                for (l4_ilist_ptr_type iPtrL4 = ini.begin(); iPtrL4 != ini.end(); ++iPtrL4)
                {
                    size_t j = 0;
                    for (l3_ilist_ptr_type iPtrL3 = iPtrL4->begin(); iPtrL3 != iPtrL4->end(); ++iPtrL3)
                    {
                        l3(i, j) = iPtrL3->size();
                        j++;
                    }
                    i++;
                }

                size_3d_array_type l4(l3, 0);
                i = 0;
                for (l4_ilist_ptr_type iPtrL4 = ini.begin(); iPtrL4 != ini.end(); ++iPtrL4)
                {
                    size_t j = 0;
                    for (l3_ilist_ptr_type iPtrL3 = iPtrL4->begin(); iPtrL3 != iPtrL4->end(); ++iPtrL3)
                    {
                        size_t k = 0;
                        for (l2_ilist_ptr_type iPtrL2 = iPtrL3->begin(); iPtrL2 != iPtrL3->end(); ++iPtrL2)
                        {
                            l4(i, j, k) = iPtrL2->size();
                            k++;
                        }
                        j++;
                    }
                    i++;
                }

                initArray(l4);
            }
            else
            {
                initArray(ini.size(),
                          ini.begin()->size(),
                          ini.begin()->begin()->size(),
                          ini.begin()->begin()->begin()->size());
            }

            // assign the values from the input initializer_list
            for (size_type l = 0; l < length1(); ++l)
            {
                l4_ilist_ptr_type iPtr4 = (ini.begin() + static_cast<std::ptrdiff_t>(l));
                for (size_type k = 0; k < length2(l); ++k)
                {
                    l3_ilist_ptr_type iPtr3 = (iPtr4->begin() + static_cast<std::ptrdiff_t>(k));
                    for (size_type j = 0; j < length3(l, k); ++j)
                    {
                        l2_ilist_ptr_type iPtr2 = (iPtr3->begin() + static_cast<std::ptrdiff_t>(j));
                        for (size_type i = 0; i < length4(l, k, j); ++i)
                        {
                            l1_ilist_ptr_type iPtr1 = (iPtr2->begin() + static_cast<std::ptrdiff_t>(i));
                            operator()(l, k, j, i) = *iPtr1;
                        }
                    }
                }
            }
        }

        /*! \brief initialize all array elements with a user-defined value
            \param[in]    ini     value for all array elements
         */
        void initArrayElements(const value_type &ini)
        {
            for (size_type i = 0; i < length1(); i++)
            {
                for (size_type j = 0; j < length2(i); j++)
                {
                    for (size_type k = 0; k < length3(i, j); k++)
                    {
                        for (size_type l = 0; l < length4(i, j, k); l++)
                        {
                            operator()(i, j, k, l) = ini;
                        }
                    }
                }
            }
        }

        /*! \brief  initialize IrregArray4D<value_type> from IrregArray4D<TForeign>,
                    used by the copy constructor and by the explicit conversion operator,
                    TForeign must be convertible to value_type via static_cast

            \tparam      TForeign   data type stored by the input data structure
            \tparam      AForeign   allocator used by the input data structure
            \tparam      FForeign   true if the input data structure uses a flat data array

            \param[in]    ini    array to be copied
         */
        template <typename TForeign, class AForeign = std::allocator<TForeign>, bool FForeign>
        void initArray(const IrregArray4D<TForeign, AForeign, FForeign> &ini)
        {
            if (ini.data() != nullptr)
            {
                if (ini.isIrreg())
                {
                    size_1d_array_type l2(ini.length1(), 0);
                    for (size_type i = 0; i < ini.length1(); ++i)
                    {
                        l2[i] = ini.length2(i);
                    }
                    size_2d_array_type l3(l2, 0);
                    for (size_type i = 0; i < ini.length1(); ++i)
                    {
                        for (size_type j = 0; j < ini.length2(i); ++j)
                        {
                            l3(i, j) = ini.length3(i, j);
                        }
                    }
                    size_3d_array_type l4(l3, 0);
                    for (size_type i = 0; i < ini.length1(); ++i)
                    {
                        for (size_type j = 0; j < ini.length2(i); ++j)
                        {
                            for (size_type k = 0; k < ini.length3(i, j); ++k)
                            {
                                l4(i, j, k) = ini.length4(i, j, k);
                            }
                        }
                    }
                    initArray(l4);
                }
                else
                {
                    initArray(ini.length1(), ini.length2(), ini.length3(), ini.length4());
                }
            }
            if (ini.data() != nullptr && arr_ != nullptr)
            {
                for (size_type i = 0; i < length1(); i++)
                {
                    for (size_type j = 0; j < length2(i); j++)
                    {
                        for (size_type k = 0; k < length3(i, j); k++)
                        {
                            for (size_type l = 0; l < length4(i, j, k); l++)
                            {
// the static analyzer of Clang 7.0.1-4 reports a false-positive garbage value when copying the array elements
// version 7.0.0 (trunk 338398) (llvm/trunk 338397) does not
#ifndef __clang_analyzer__
                                operator()(i, j, k, l) = static_cast<value_type>(ini(i, j, k, l));
#endif                          // __clang_analyzer__
                            }
                        }
                    }
                }
            }
        }

        /*! \brief specialized constructor helper function for regularly shaped matrices

            \param[in]    l1      array length in dimension 1
            \param[in]    l2      array length in dimension 2
            \param[in]    l3      array length in dimension 3
            \param[in]    l4      array length in dimension 4
         */
        void initArray(size_type l1, size_type l2, size_type l3, size_type l4)
        {
            deallocateArray();
            // otherwise, just keep the empty array
            if (l1 > 0 && l2 > 0 && l3 > 0 && l4 > 0)
            {
                length1_   = l1;
                length2_   = l2;
                length3_   = l3;
                length4_   = l4;
                nelements_ = l1 * l2 * l3 * l4;
                allocateArray();
            }
        }

        /*! \brief specialized constructor helper function for regularly shaped matrices

            \param[in]    l1      array length in dimension 1
            \param[in]    l2      array length in dimension 2
            \param[in]    l3      array length in dimension 3
            \param[in]    l4      array length in dimension 4
            \param[in]    ini     initialization value for the data array elements
         */
        void initArray(size_type l1, size_type l2, size_type l3, size_type l4, const value_type &ini)
        {
            initArray(l1, l2, l3, l4);
            initArrayElements(ini);
        }

        /*! \brief constructor helper function for irregularly shaped matrices and all special cases (IsFlat == false)

            \tparam      AForeign   allocator used by the input data structure
            \tparam      FForeign      true if the input data structure uses a flat data array
            \tparam      DummyIsFlat   dummy template parameter for SFINAE selection according to IsFlat

            \param[in]   l4   number of elements in dimension 4 as function of the indices in the lower dimensions
                              number of elements in dimensions 1,2,3 also deducible from l4
         */
        template<class AForeign = std::allocator<value_type>, bool FForeign, bool DummyIsFlat = IsFlat>
        typename std::enable_if<!DummyIsFlat, void>::type
        initArray(const IrregArray3D<size_type, AForeign, FForeign> &l4)
        {
            GMX_ASSERT(IsFlat == DummyIsFlat, "Explicit instantiation with DummyIsFlat != IsFlat prohibited.");
            deallocateArray();
            bool length2AlwaysGreaterZero  = true;
            bool length3AlwaysGreaterZero  = true;
            for (size_type i = 0; i < l4.length1(); ++i)
            {
                if (l4.length2(i) <= 0)
                {
                    length2AlwaysGreaterZero = false;
                }
                for (size_type j = 0; j < l4.length2(i); ++j)
                {
                    if (l4.length3(i, j) <= 0)
                    {
                        length3AlwaysGreaterZero = false;
                    }
                }
            }

            // otherwise, just keep the empty array
            if (l4.length1() > 0 && length2AlwaysGreaterZero && length3AlwaysGreaterZero)
            {
                isIrreg_   = true;
                length1_   = l4.length1();
                irreg_.initArray(l4);
                computeNelements();
                allocateArray();
            }
        }
        /*! \brief constructor helper function for irregularly shaped matrices and all special cases (IsFlat == true)

            \tparam      AForeign   allocator used by the input data structure
            \tparam      FForeign      true if the input data structure uses a flat data array
            \tparam      DummyIsFlat   dummy template parameter for SFINAE selection according to IsFlat

            \param[in]   l4   number of elements in dimension 4 as function of the indices in the lower dimensions
                              number of elements in dimensions 1,2,3 also deducible from l4
         */
        template<class AForeign = std::allocator<value_type>, bool FForeign, bool DummyIsFlat = IsFlat>
        typename std::enable_if<DummyIsFlat, void>::type
        initArray(const IrregArray3D<size_type, AForeign, FForeign> &l4)
        {
            GMX_ASSERT(IsFlat == DummyIsFlat, "Explicit instantiation with DummyIsFlat != IsFlat prohibited.");
            deallocateArray();
            bool length2AlwaysGreaterZero  = true;
            bool length3AlwaysGreaterZero  = true;
            for (size_type i = 0; i < l4.length1(); ++i)
            {
                if (l4.length2(i) <= 0)
                {
                    length2AlwaysGreaterZero = false;
                }
                for (size_type j = 0; j < l4.length2(i); ++j)
                {
                    if (l4.length3(i, j) <= 0)
                    {
                        length3AlwaysGreaterZero = false;
                    }
                }
            }

            // otherwise, just keep the empty array
            if (l4.length1() > 0 && length2AlwaysGreaterZero && length3AlwaysGreaterZero)
            {
                isIrreg_ = true;
                length1_ = l4.length1();
                length2_ = 0;
                length3_ = 0;
                length4_ = 0;
                size_1d_array_type l2(l4.length1(), 0);
                for (size_type i = 0; i < l4.length1(); i++)
                {
                    l2[i] = l4.length2(i);
                }
                size_2d_array_type l3(l2);
                for (size_type i = 0; i < l4.length1(); i++)
                {
                    for (size_type j = 0; j < l4.length2(i); j++)
                    {
                        l3(i, j) = l4.length3(i, j);
                    }
                }

                irreg_.initArray(l3);
                size_type  last_length = 0;
                index_type last_index  = 0;
                for (size_type i = 0; i < l4.length1(); i++)
                {
                    for (size_type j = 0; j < l4.length2(i); j++)
                    {
                        for (size_type k = 0; k < l4.length3(i, j); k++)
                        {
                            const size_type thisLength4 = l4(i, j, k);
                            if (i == 0 && j == 0 && k == 0)
                            {
                                irreg_(i, j, k) = 0;
                                last_index      = irreg_(i, j, k);
                                last_length     = thisLength4;
                            }
                            else
                            {
                                irreg_(i, j, k) = last_index + last_length;
                                last_index      = irreg_(i, j, k);
                                last_length     = thisLength4;
                            }
                        }
                    }
                }
                // the number of elements is given by the starting index of the last array stripe
                // FlatIrregArray2D[i][:] in dimension 4 within arr + the length of this stripe
                nelements_ = last_index + last_length;
                allocateArray();
            }
        }
        /*! \brief general constructor helper function for irregularly shaped matrices and all special cases

            \tparam      AForeign   allocator used by the input data structure
            \tparam      FForeign   true if the input data structure uses a flat data array

            \param[in]   l4    number of elements in dimension 4 as function of the indices in the lower dimensions
                               number of elements in dimensions 1,2,3 also deducible from l4
            \param[in]   ini   initialization value for the data array elements
         */
        template<class AForeign = allocator_size_type, bool FForeign>
        void initArray(const IrregArray3D<size_type, AForeign, FForeign> &l4, const value_type &ini)
        {
            initArray(l4);
            initArrayElements(ini);
        }

        //! returns true if the array has irregular shape
        inline bool isIrreg() const noexcept { return isIrreg_; }

        //! returns true if the underlying data array is 1-D (flat), false if it is 4D
        constexpr inline bool isFlat() const noexcept { return IsFlat; }

        //! number of elements in dimension 1
        inline size_type length1() const noexcept { return length1_; }
        //! number of elements in dimension 2 (regular array)
        inline size_type length2() const
        {
            if (isIrreg_)
            {
// clang-tidy diagnoses that an exception may be thrown in cases there isIrreg_ is guaranteed to be false
#ifndef __clang_analyzer__
                std::string errorMessage(formatString("Error in %s::length2(): called without argument for an irregular array.", typeid(*this).name()));
                GMX_THROW(APIError(errorMessage));
#endif          // __clang_analyzer__
            }
            return length2_;
        }
        /*! \brief number of elements in dimension 2 at lower dimension index x

            \param[in]    x       index in dimension 1 for which the array length in dimension 2 is requested */
        inline size_type length2(index_type x) const noexcept
        {
            return (isIrreg_ ? irreg_.length2(x) : length2_);
        }
        //! number of elements in dimension 3 (regular array)
        inline size_type length3() const
        {
            if (isIrreg_)
            {
// clang-tidy diagnoses that an exception may be thrown in cases there isIrreg_ is guaranteed to be false
#ifndef __clang_analyzer__
                std::string errorMessage(formatString("Error in %s::length3(): called without argument for an irregular array.", typeid(*this).name()));
                GMX_THROW(APIError(errorMessage));
#endif          // __clang_analyzer__
            }
            return length3_;
        }
        /*! \brief number of elements in dimension 3 at lower dimension indices x, y
            \param[in]    x       index in dimension 1
            \param[in]    y       index in dimension 2  */
        inline size_type length3(index_type x, index_type y) const noexcept
        {
            return (isIrreg_ ? irreg_.length3(x, y) : length3_);
        }
        //! number of elements in dimension 4 (regular array)
        inline size_type length4() const
        {
            if (isIrreg_)
            {
// clang-tidy diagnoses that an exception may be thrown in cases there isIrreg_ is guaranteed to be false
#ifndef __clang_analyzer__
                std::string errorMessage(formatString("Error in %s::length4(): called without argument for an irregular array.", typeid(*this).name()));
                GMX_THROW(APIError(errorMessage));
#endif          // __clang_analyzer__
            }
            return length4_;
        }
        /*! \brief number of elements in dimension 4 at lower dimension indices x, y
            \tparam       DummyIsFlat   dummy template parameter for SFINAE selection according to IsFlat
            \param[in]    x       index in dimension 1
            \param[in]    y       index in dimension 2
            \param[in]    z       index in dimension 3  */
        template<bool DummyIsFlat = IsFlat>
        inline typename std::enable_if<!DummyIsFlat, size_type>::type
        length4(index_type x, index_type y, index_type z) const noexcept
        {
            GMX_ASSERT(IsFlat == DummyIsFlat, "Explicit instantiation with DummyIsFlat != IsFlat prohibited.");
            return (isIrreg_ ? irreg_(x, y, z) : length4_);
        }
        /*! \brief number of elements in dimension 4 at lower dimension indices x, y
            \tparam       DummyIsFlat   dummy template parameter for SFINAE selection according to IsFlat
            \param[in]    x       index in dimension 1
            \param[in]    y       index in dimension 2
            \param[in]    z       index in dimension 3  */
        template<bool DummyIsFlat = IsFlat>
        inline typename std::enable_if<DummyIsFlat, size_type>::type
        length4(index_type x, index_type y, index_type z) const noexcept
        {
            GMX_ASSERT(IsFlat == DummyIsFlat, "Explicit instantiation with DummyIsFlat != IsFlat prohibited.");
            if (!isIrreg_)
            {
                return length4_;
            }
            else
            {
                const index_type last1 = length1_ - 1;
                const index_type last2 = length2(last1) - 1;
                const index_type last3 = length3(last1, last2) - 1;
                if (x == last1 && y == last2 && z == last3)
                {
                    return nelements_ - irreg_(x, y, z);
                }
                else
                {
                    // get first element of next stripe in the underlying data array of irreg_mat to avoid further conditionals
                    const size_type        i = irreg_.arrayIndex(x, y, z);
                    const size_type* const a = irreg_.data();
                    return a[i + 1] - a[i];
                }
            }
        }

        /*! \brief determine the index of the array element IrregArray2D(i, j) in the underlying 1D array arr_ (IsFlat = true)
            This function is public because it may be useful in client code, too.

            size_type instead of index_type used intentionally here to use the maximum possible length

            \tparam      DummyIsFlat   dummy template parameter for SFINAE selection according to IsFlat
            \param[in]   x    index in dimension 1
            \param[in]   y    index in dimension 2
            \param[in]   z    index in dimension 3
            \param[in]   k    index in dimension 4
         */
        template<bool DummyIsFlat = IsFlat>
        inline typename std::enable_if<DummyIsFlat, size_type>::type
        arrayIndex(index_type x, index_type y, index_type z, index_type k) const noexcept
        {
            if (isIrreg_)
            {
                return irreg_(x, y, z) + k;
            }
            else
            {
                return length2() * length3() * length4() * x
                       + length3() * length4() * y
                       + length4() * z
                       + k;
            }
        }

        //! returns the number of data array elements
        size_type nElements() const noexcept { return nelements_; }
        //! compute the number of data array elements
        size_type computeNelements() noexcept
        {
            nelements_ = 0;
            for (size_type i = 0; i < length1(); i++)
            {
                for (size_type j = 0; j < length2(i); j++)
                {
                    for (size_type k = 0; k < length3(i, j); k++)
                    {
                        nelements_ += length4(i, j, k);
                    }
                }
            }
            return nelements_;
        }

        //! deallocate the memory of all constituent arrays
        void deallocateArray() noexcept
        {
            if (arr_ != nullptr)
            {
                destroyArrayElements();
                deallocateMemory();
            }
            // deallocate book-keeping arrays
            if (isIrreg_)
            {
                irreg_.deallocateArray();
            }
            isIrreg_   = false;
            length1_   = 0;
            length2_   = 0;
            length3_   = 0;
            length4_   = 0;
            nelements_ = 0;
        }
        //! destructor
        ~IrregArray4D() noexcept
        {
            deallocateArray();
        }

        //! doxygen can not handle nested classes
        /// \cond DEV
        /*! \class Proxy

            \brief helper functor used to implement index operator *this[i][j][k][l]
            The trick here is to chain operator[] of the nested subclasses.
         */
        class Proxy
        {
            private:
                IrregArray4D    * optr_;
                index_type        x_;
                class Proxy2
                {
                    private:
                        IrregArray4D    * optr_;
                        index_type        x_;
                        index_type        y_;
                        class Proxy3
                        {
                            private:
                                IrregArray4D    * optr_;
                                index_type        x_;
                                index_type        y_;
                                index_type        z_;
                            public:
                                Proxy3(IrregArray4D* optr, index_type x, index_type y, index_type z) noexcept
                                    : optr_(optr), x_(x), y_(y), z_(z) { }

                                value_type &operator[](index_type k) noexcept
                                {
                                    return optr_->operator()(x_, y_, z_, k);
                                }
                                const value_type &operator[](index_type k) const noexcept
                                {
                                    return optr_->operator()(x_, y_, z_, k);
                                }
                        };
                    public:
                        Proxy2(IrregArray4D* optr, index_type x, index_type y) noexcept
                            : optr_(optr), x_(x), y_(y) { }

                        Proxy3 operator[](index_type z) noexcept
                        {
                            return Proxy3(optr_, x_, y_, z);
                        }
                        Proxy3 operator[](index_type z) const noexcept
                        {
                            return Proxy3(optr_, x_, y_, z);
                        }
                };
            public:
                Proxy(IrregArray4D* optr, index_type x) noexcept
                    : optr_(optr), x_(x) { }

                Proxy2 operator[](index_type y) noexcept
                {
                    return Proxy2(optr_, x_, y);
                }
                Proxy2 operator[](index_type y) const noexcept
                {
                    return Proxy2(optr_, x_, y);
                }
        };
        /// \endcond DEV
        /*! \brief   index operator[x][y][z][k] for array-like usage realized with chained functors "Proxy*"
            \param[in]   x   index in dimension 1
            \cond
            \param[in]   y   index in dimension 2
            \param[in]   z   index in dimension 3
            \param[in]   k   index in dimension 4
            \endcond
         */
        Proxy operator[](index_type x) noexcept
        {
            return Proxy(this, x);
        }
        /*! \brief   const index operator[x][y][z][k] for array-like usage realized with chained functors "Proxy*"
            \param[in]   x   index in dimension 1
            \cond
            \param[in]   y   index in dimension 2
            \param[in]   z   index in dimension 3
            \param[in]   k   index in dimension 4
            \endcond
         */
        Proxy operator[](index_type x) const noexcept
        {
            return Proxy(this, x);
        }

        /*! \brief   index operator(x,y,z) for functor-like usage
            \tparam      DummyIsFlat   dummy template parameter for SFINAE selection according to IsFlat
            \param[in]   x   index in dimension 1
            \param[in]   y   index in dimension 2
            \param[in]   z   index in dimension 3
            \param[in]   k   index in dimension 4
         */
        template<bool DummyIsFlat = IsFlat>
        typename std::enable_if<!DummyIsFlat, value_type &>::type
        operator()(index_type x, index_type y, index_type z, index_type k) noexcept
        {
            GMX_ASSERT(IsFlat == DummyIsFlat, "Explicit instantiation with DummyIsFlat != IsFlat prohibited.");
            return arr_[x][y][z][k];
        }
        /*! \brief   const index operator(x,y,z) for functor-like usage
            \tparam      DummyIsFlat   dummy template parameter for SFINAE selection according to IsFlat
            \param[in]   x   index in dimension 1
            \param[in]   y   index in dimension 2
            \param[in]   z   index in dimension 3
            \param[in]   k   index in dimension 4
         */
        template<bool DummyIsFlat = IsFlat>
        typename std::enable_if<!DummyIsFlat, const value_type &>::type
        operator()(index_type x, index_type y, index_type z, index_type k) const noexcept
        {
            GMX_ASSERT(IsFlat == DummyIsFlat, "Explicit instantiation with DummyIsFlat != IsFlat prohibited.");
            return arr_[x][y][z][k];
        }
        /*! \brief   index operator(x,y,z) for functor-like usage
            \tparam      DummyIsFlat   dummy template parameter for SFINAE selection according to IsFlat
            \param[in]   x   index in dimension 1
            \param[in]   y   index in dimension 2
            \param[in]   z   index in dimension 3
         */
        template<bool DummyIsFlat = IsFlat>
        inline typename std::enable_if<!DummyIsFlat, value_type*>::type
        operator()(index_type x, index_type y, index_type z) noexcept
        {
            GMX_ASSERT(IsFlat == DummyIsFlat, "Explicit instantiation with DummyIsFlat != IsFlat prohibited.");
            return arr_[x][y][z];
        }
        /*! \brief   const index operator(x,y,z) for functor-like usage
            \tparam      DummyIsFlat   dummy template parameter for SFINAE selection according to IsFlat
            \param[in]   x   index in dimension 1
            \param[in]   y   index in dimension 2
            \param[in]   z   index in dimension 3
         */
        template<bool DummyIsFlat = IsFlat>
        inline typename std::enable_if<!DummyIsFlat, const value_type*>::type
        operator()(index_type x, index_type y, index_type z) const noexcept
        {
            GMX_ASSERT(IsFlat == DummyIsFlat, "Explicit instantiation with DummyIsFlat != IsFlat prohibited.");
            return arr_[x][y][z];
        }

        //! get a pointer to the data array
        //! \tparam   DummyIsFlat   dummy template parameter for SFINAE selection according to IsFlat
        template<bool DummyIsFlat = IsFlat>
        typename std::enable_if<!DummyIsFlat, value_type****>::type
        data() noexcept
        {
            GMX_ASSERT(IsFlat == DummyIsFlat, "Explicit instantiation with DummyIsFlat != IsFlat prohibited.");
            return arr_;
        }
        //! get a pointer to the data array
        //! \tparam   DummyIsFlat   dummy template parameter for SFINAE selection according to IsFlat
        template<bool DummyIsFlat = IsFlat>
        typename std::enable_if<!DummyIsFlat, const value_type****>::type
        data() const noexcept
        {
            GMX_ASSERT(IsFlat == DummyIsFlat, "Explicit instantiation with DummyIsFlat != IsFlat prohibited.");
            return const_cast<const value_type****>(arr_);
        }

        /*! \brief   index operator(x,y,z) for functor-like usage
            \tparam      DummyIsFlat   dummy template parameter for SFINAE selection according to IsFlat
            \param[in]   x   index in dimension 1
            \param[in]   y   index in dimension 2
            \param[in]   z   index in dimension 3
            \param[in]   k   index in dimension 4
         */
        template<bool DummyIsFlat = IsFlat>
        inline typename std::enable_if<DummyIsFlat, value_type &>::type
        operator()(index_type x, index_type y, index_type z, const index_type k) noexcept
        {
            GMX_ASSERT(IsFlat == DummyIsFlat, "Explicit instantiation with DummyIsFlat != IsFlat prohibited.");
            return arr_[arrayIndex(x, y, z, k)];
        }
        /*! \brief   const index operator(x,y,z) for functor-like usage
            \tparam      DummyIsFlat   dummy template parameter for SFINAE selection according to IsFlat
            \param[in]   x   index in dimension 1
            \param[in]   y   index in dimension 2
            \param[in]   z   index in dimension 3
            \param[in]   k   index in dimension 4
         */
        template<bool DummyIsFlat = IsFlat>
        inline typename std::enable_if<DummyIsFlat, const value_type &>::type
        operator()(index_type x, index_type y, index_type z, index_type k) const noexcept
        {
            GMX_ASSERT(IsFlat == DummyIsFlat, "Explicit instantiation with DummyIsFlat != IsFlat prohibited.");
            return arr_[arrayIndex(x, y, z, k)];
        }

        /*! \brief   index operator(x,y,z) for functor-like usage
            \tparam      DummyIsFlat   dummy template parameter for SFINAE selection according to IsFlat
            \param[in]   x   index in dimension 1
            \param[in]   y   index in dimension 2
            \param[in]   z   index in dimension 3
         */
        template<bool DummyIsFlat = IsFlat>
        inline typename std::enable_if<DummyIsFlat, value_type*>::type
        operator()(index_type x, index_type y, index_type z) noexcept
        {
            GMX_ASSERT(IsFlat == DummyIsFlat, "Explicit instantiation with DummyIsFlat != IsFlat prohibited.");
            return &arr_[arrayIndex(x, y, z, 0)];
        }
        /*! \brief   const index operator(x,y,z) for functor-like usage
            \tparam      DummyIsFlat   dummy template parameter for SFINAE selection according to IsFlat
            \param[in]   x   index in dimension 1
            \param[in]   y   index in dimension 2
            \param[in]   z   index in dimension 3
         */
        template<bool DummyIsFlat = IsFlat>
        inline typename std::enable_if<DummyIsFlat, const value_type*>::type
        operator()(index_type x, index_type y, index_type z) const noexcept
        {
            GMX_ASSERT(IsFlat == DummyIsFlat, "Explicit instantiation with DummyIsFlat != IsFlat prohibited.");
            return &arr_[arrayIndex(x, y, z, 0)];
        }

        //! get a pointer to the data array
        //! \tparam   DummyIsFlat   dummy template parameter for SFINAE selection according to IsFlat
        template<bool DummyIsFlat = IsFlat>
        typename std::enable_if<DummyIsFlat, value_type*>::type
        data() noexcept
        {
            GMX_ASSERT(IsFlat == DummyIsFlat, "Explicit instantiation with DummyIsFlat != IsFlat prohibited.");
            return arr_;
        }
        //! get a pointer to the data array
        //! \tparam   DummyIsFlat   dummy template parameter for SFINAE selection according to IsFlat
        template<bool DummyIsFlat = IsFlat>
        typename std::enable_if<DummyIsFlat, const value_type*>::type
        data() const noexcept
        {
            GMX_ASSERT(IsFlat == DummyIsFlat, "Explicit instantiation with DummyIsFlat != IsFlat prohibited.");
            return const_cast<const value_type*>(arr_);
        }

        /*! \brief   stream operator for convenient printing of the array to an ostream
            \param[in]   output   output stream for the array contents
            \param[in]   ten      the array to be printed
         */
        friend std::ostream &operator<<(std::ostream &output, const IrregArray4D &ten)
        {
            // safe the current ostream format for restoring it after output insertion
            std::ios  state(nullptr);
            state.copyfmt(output);

            if (ten.arr_ != nullptr)
            {
                constexpr size_t floatWidth = 10;
                constexpr size_t floatPrec  =  2;
                for (size_type i = 0; i < ten.length1(); i++)
                {
                    for (size_type j = 0; j < ten.length2(i); j++)
                    {
                        output << "------------------------------" << std::endl << i << ", " << j << std::endl;
                        for (size_type k = 0; k < ten.length3(i, j); k++)
                        {
                            for (size_type l = 0; l < ten.length4(i, j, k); l++)
                            {
                                if (std::is_floating_point<value_type>::value)
                                {
                                    output << std::setw(floatWidth) << std::fixed << std::setprecision(floatPrec) << ten(i, j, k, l);
                                }
                                else
                                {
                                    output << ten(i, j, k, l);
                                }
                                if (k < (ten.length4(i, j, k) - 1))
                                {
                                    output << " ";
                                }
                            }
                            output << "\n";
                        }
                        output << "\n";
                    }
                    output << "\n";
                }
                output << std::endl;
            }

            // restore the original ostream format
            output.copyfmt(state);

            return output;
        }

    private:
        //! allocate memory for the array
        void allocateArray()
        {
            if (arr_ != nullptr)
            {
                std::string errorMessage(formatString("Error in %s::allocateArray(): called allocate without deallocating previously allocated memory.", typeid(*this).name()));
                GMX_THROW(InternalError(errorMessage));
            }
            if (length1() <= 0)
            {
                std::string errorMessage(formatString("Error in %s::allocateArray(): called with illegal bounds length1 <= 0.", typeid(*this).name()));
                GMX_THROW(InvalidInputError(errorMessage));
            }
            for (size_type i = 0; i < length1(); ++i)
            {
                if (length2(i) <= 0)
                {
                    std::string errorMessage(formatString("Error in %s::allocateArray(): called with illegal bounds length2 <= 0 at i = %lu.", typeid(*this).name(), i));
                    GMX_THROW(InvalidInputError(errorMessage));
                }
                for (size_type j = 0; j < length2(i); ++j)
                {
                    if (length3(i, j) <= 0)
                    {
                        std::string errorMessage(formatString("Error in %s::allocateArray(): called with illegal bounds length3 <= 0 at i, j = %lu, %lu.", typeid(*this).name(), i, j));
                        GMX_THROW(InvalidInputError(errorMessage));
                    }
                    for (size_type k = 0; k < length3(i, j); ++k)
                    {
                        if (static_cast<std::make_signed<size_type>::type>(length4(i, j, k)) < 0)
                        {
                            std::string errorMessage(formatString("Error in %s::allocateArray(): called with illegal bounds length4 < 0 at i, j, k = %lu, %lu, %lu.", typeid(*this).name(), i, j, k));
                            GMX_THROW(InvalidInputError(errorMessage));
                        }
                    }
                }
            }
            allocateMemory();
            constructArrayElements();
        }

        //! allocate memory for the array (IsFlat == false), only to be called by allocateArray()
        template<bool DummyIsFlat = IsFlat>
        typename std::enable_if<!DummyIsFlat, void>::type
        allocateMemory()
        {
            try
            {
                arr_  = std::allocator_traits<p3_allocator_type>::allocate(p3_allocator_, length1());
                for (size_type i = 0; i < length1(); ++i)
                {
                    arr_[i] = std::allocator_traits<p2_allocator_type>::allocate(p2_allocator_, length2(i));
                    for (size_type j = 0; j < length2(i); ++j)
                    {
                        arr_[i][j] = std::allocator_traits<p_allocator_type>::allocate(p_allocator_, length3(i, j));
                        for (size_type k = 0; k < length3(i, j); ++k)
                        {
                            arr_[i][j][k] = std::allocator_traits<allocator_type>::allocate(allocator_, length4(i, j, k));
                        }
                    }
                }
            }
            catch (const std::bad_alloc &)
            {
                std::fprintf(stderr, "Error in %s::allocateMemory(): could not allocate memory for the data array.", typeid(*this).name());
                throw;
            }
        }
        //! allocate memory for the array (IsFlat == true), only to be called by allocateArray()
        template<bool DummyIsFlat = IsFlat>
        typename std::enable_if<DummyIsFlat, void>::type
        allocateMemory()
        {
            try
            {
                arr_ = std::allocator_traits<allocator_type>::allocate(allocator_, nelements_);
            }
            catch (const std::bad_alloc &)
            {
                std::fprintf(stderr, "Error in %s::allocateMemory(): could not allocate memory for the data array.", typeid(*this).name());
                throw;
            }
        }

        //! default construct each array element, only to be called by allocateArray()
        void constructArrayElements()
        {
            if (arr_ != nullptr)
            {
                // destroy the array elements
                for (size_type i = 0; i < length1(); ++i)
                {
                    for (size_type j = 0; j < length2(i); ++j)
                    {
                        for (size_type k = 0; k < length3(i, j); ++k)
                        {
                            for (size_type l = 0; l < length4(i, j, k); ++l)
                            {
                                std::allocator_traits<allocator_type>::construct(allocator_, &this->operator()(i, j, k, l), value_type());
                            }
                        }
                    }
                }
            }
        }

        //! destroy each array element, only to be called by deallocateArray()
        void destroyArrayElements() noexcept
        {
            if (arr_ != nullptr)
            {
                // destroy the array elements
                for (size_type i = 0; i < length1(); ++i)
                {
                    for (size_type j = 0; j < length2(i); ++j)
                    {
                        for (size_type k = 0; k < length3(i, j); ++k)
                        {
                            for (size_type l = 0; l < length4(i, j, k); ++l)
                            {
                                std::allocator_traits<allocator_type>::destroy(allocator_, &this->operator()(i, j, k, l));
                            }
                        }
                    }
                }
            }
        }

        //! deallocate array memory (IsFlat == false), only to be called by deallocateArray()
        template<bool DummyIsFlat = IsFlat>
        typename std::enable_if<!DummyIsFlat, void>::type
        deallocateMemory() noexcept
        {
            if (arr_ != nullptr)
            {
                for (size_type i = 0; i < length1(); ++i)
                {
                    for (size_type j = 0; j < length2(i); ++j)
                    {
                        for (size_type k = 0; k < length3(i, j); ++k)
                        {
                            std::allocator_traits<allocator_type>::deallocate(allocator_, arr_[i][j][k], length4(i, j, k));
                        }
                        std::allocator_traits<p_allocator_type>::deallocate(p_allocator_, arr_[i][j], length3(i, j));
                    }
                    std::allocator_traits<p2_allocator_type>::deallocate(p2_allocator_, arr_[i], length2(i));
                }
                std::allocator_traits<p3_allocator_type>::deallocate(p3_allocator_, arr_, length1());
                arr_ = nullptr;
            }
        }
        //! deallocate array memory (IsFlat == true), only to be called by deallocateArray()
        template<bool DummyIsFlat = IsFlat>
        typename std::enable_if<DummyIsFlat, void>::type
        deallocateMemory() noexcept
        {
            if (arr_ != nullptr)
            {
                std::allocator_traits<allocator_type>::deallocate(allocator_, arr_, nelements_);
                arr_ = nullptr;
            }
        }
};

} // end namespace gmx

#endif
