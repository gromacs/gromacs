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
/*! \file
    \ingroup module_utility
    \inpublicapi

    \brief IrregArray4D is an irregularly shaped 4-dimensional array

    \author R. Thomas Ullmann <tullman@gwdg.de>
 */
#ifndef GMX_UTILITY_DATASTRUCTURES_IRREG_ARRAY_4D_H
#define GMX_UTILITY_DATASTRUCTURES_IRREG_ARRAY_4D_H

#include "gromacs/utility/data_structures/irreg_array_3d.h"

namespace gmx
{


/*! \class IrregArray4D

    \brief irregularly shaped 4-dimensional array

     this array has nonuniform second, third and fourth dimensions,
     that is, instead of a regular array with constant dimensions
     M x N x O x P, an irregular array with dimensions
     M x N(M) x O(M, N) x P(M, N, O).

     example application:
     interaction matrix indexed by site1, form1, site2, form2

    \ingroup module_utility
    \inpublicapi

    \tparam   T           data type to be stored
    \tparam   Allocator   allocator to be used in creating arrays storing T
                          allocators for storing additional book-keeping data
                          are derived from this allocator via rebind if necessary
 */
template<class T = size_t, class Allocator = std::allocator<T> >
class IrregArray4D
{
    public:
        /* adopted these typedefs from AnBe's TriangularArray, for consistency with the other data
            structures, can be retrieved from outside to make sure that the intended types are used
            when accessing or manipulating the array */
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
        //! the allocator type to allocate value_type*
        typedef typename std::allocator_traits<Allocator>::pointer allocator_return_ptr_type;
        //! the (integer) datatype used to represent sizes, e.g., array lengths
        typedef  size_t  size_type;
        //! the (integer) datatype used to represent array indices
        typedef ssize_t index_type;
        //! the allocator type to allocate size_type
        typedef typename std::allocator_traits<Allocator>::template rebind_alloc<size_type> allocator_size_type;
        //! type of array used for storing the lengths of array stripes
        typedef IrregArray1D<size_type, allocator_size_type> size_1d_array_type;
        //! type of array used for storing the lengths of array stripes
        typedef IrregArray2D<size_type, allocator_size_type> size_2d_array_type;
        //! type of array used for storing the lengths of array stripes
        typedef IrregArray3D<size_type, allocator_size_type> size_3d_array_type;
    private:
        //! first index of dimension 1
        index_type              first1_;
        //! first index of dimension 2
        index_type              first2_;
        //! first index of dimension 3
        index_type              first3_;
        //! first index of dimension 4
        index_type              first4_;
        //! ending index of dimension 1
        index_type              last1_;
        //! ending index of dimension 2
        index_type              last2_;
        //! ending index of dimension 3
        index_type              last3_;
        //! ending index of dimension 4
        index_type              last4_;
        //! dimension1,2,3-specific number of elements of dimension 4
        IrregArray3D<size_type, allocator_size_type> irreg_;
        //! flag says that dimX(i) = last_X(i) - first_X(i) + 1 is not a constant
        bool                    isIrreg_;
        //! array storing the the actual data
        value_type          ****arr_;
        //! the allocator_ used to allocate arr_
        p3_allocator_type       p3_allocator_;
        //! the allocator_ used to allocate arr_[i]
        p2_allocator_type       p2_allocator_;
        //! the allocator_ used to allocate arr_[i][j]
        p_allocator_type        p_allocator_;
        //! the allocator_ used to allocate arr_[i][j][k]
        allocator_type          allocator_;
        //! let other IrregArray4D instantiations access private members (for the explicit conversion operator)
        template<class W, typename AW> friend class IrregArray4D;
    public:
        /*! \brief   default constructor without memory allocation, allocateArray has to be called for
                     allocating memory initArray for memory allocation and content initialization
         */
        IrregArray4D()
            : first1_(0), first2_(0), first3_(0), first4_(0), last1_(-1), last2_(-1), last3_(-1), last4_(-1), isIrreg_(false), arr_(nullptr)
        {
        }

        /*! \brief  copy constructor

            \param[in]   ini   array to be copied
         */
        IrregArray4D(const IrregArray4D &ini)
            : first1_(0), first2_(0), first3_(0), first4_(0), last1_(-1), last2_(-1), last3_(-1), last4_(-1), isIrreg_(false), arr_(nullptr)
        {
            initArray(ini);
        }

        /*! \brief  copy constructor

            \tparam      AForeign   allocator used by the input data structure

            \param[in]   ini   array to be copied
         */
        template<class AForeign = allocator_type>
        IrregArray4D(const IrregArray4D<value_type, AForeign> &ini)
            : first1_(0), first2_(0), first3_(0), first4_(0), last1_(-1), last2_(-1), last3_(-1), last4_(-1), isIrreg_(false), arr_(nullptr)
        {
            initArray(ini);
        }

        //! move constructor
        //! \param[in]    ini   source to be moved
        IrregArray4D(IrregArray4D &&ini) noexcept
            : first1_(ini.first1_), first2_(ini.first2_), first3_(ini.first3_), first4_(ini.first4_),
              last1_(ini.last1_), last2_(ini.last2_), last3_(ini.last3_), last4_(ini.last4_),
              irreg_(std::move(ini.irreg_)), isIrreg_(ini.isIrreg_), arr_(ini.arr_),
              p3_allocator_(std::move(ini.p3_allocator_)), p2_allocator_(std::move(ini.p2_allocator_)),
              p_allocator_(std::move(ini.p_allocator_)), allocator_(std::move(ini.allocator_))
        {
            // set ini to a valid default state
            ini.first1_   =  0;
            ini.first2_   =  0;
            ini.first3_   =  0;
            ini.first4_   =  0;
            ini.last1_    = -1;
            ini.last2_    = -1;
            ini.last3_    = -1;
            ini.last4_    = -1;
            ini.isIrreg_  = false;
            ini.arr_      = nullptr;
        }

        /*! \brief   constructor for a regularly shaped array

            \param[in]    x1       starting index in dimension 1
            \param[in]    x2       last index in dimension 1
            \param[in]    y1       starting index in dimension 2
            \param[in]    y2       last index in dimension 2
            \param[in]    z1       starting index in dimension 3
            \param[in]    z2       last index in dimension 3
            \param[in]    k1       starting index in dimension 4
            \param[in]    k2       last index in dimension 4
         */
        IrregArray4D(const index_type x1, const index_type x2, const index_type y1, const index_type y2, const index_type z1, const index_type z2, const index_type k1, const index_type k2)
            : first1_(0), first2_(0), first3_(0), first4_(0), last1_(-1), last2_(-1), last3_(-1), last4_(-1), isIrreg_(false), arr_(nullptr)
        {
            initArray(x1, x2, y1, y2, z1, z2, k1, k2);
        }
        /*! \brief   constructor for a regularly shaped data array with initialization

            \param[in]    x1       starting index in dimension 1
            \param[in]    x2       last index in dimension 1
            \param[in]    y1       starting index in dimension 2
            \param[in]    y2       last index in dimension 2
            \param[in]    z1       starting index in dimension 3
            \param[in]    z2       last index in dimension 3
            \param[in]    k1       starting index in dimension 4
            \param[in]    k2       last index in dimension 4
            \param[in]    ini      initialization value for the data array elements
         */
        IrregArray4D(const index_type x1, const index_type x2, const index_type y1, const index_type y2, const index_type z1, const index_type z2, const index_type k1, const index_type k2, const value_type &ini)
            : first1_(0), first2_(0), first3_(0), first4_(0), last1_(-1), last2_(-1), last3_(-1), last4_(-1), isIrreg_(false), arr_(nullptr)
        {
            initArray(x1, x2, y1, y2, z1, z2, k1, k2, ini);
        }

        /*! \brief   general constructor helper function for irregularly shaped matrices and all special cases

            \tparam      AForeign   allocator used by the input data structure

            \param[in]    x1       starting index in dimension 1
            \param[in]    x2       last index in dimension 1
            \param[in]    sizes1   number of elements in dimension 2 as function of the index in dimension 1
            \param[in]    sizes2   number of elements in dimension 3 as function of the indices in the lower dimensions
            \param[in]    sizes3   number of elements in dimension 4 as function of the indices in the lower dimensions
         */
        template <class AForeign = allocator_size_type>
        IrregArray4D(const index_type x1, const index_type x2, const IrregArray1D<size_type, AForeign> &sizes1, const IrregArray2D<size_type, AForeign> &sizes2, const IrregArray3D<size_type, AForeign> &sizes3)
            : first1_(0), first2_(0), first3_(0), first4_(0), last1_(-1), last2_(-1), last3_(-1), last4_(-1), isIrreg_(true), arr_(nullptr)
        {
            initArray(x1, x2, sizes1, sizes2, sizes3);
        }
        /*! \brief   general constructor helper function for irregularly shaped matrices and all special cases

            \tparam      AForeign   allocator used by the input data structure

            \param[in]    x1       starting index in dimension 1
            \param[in]    x2       last index in dimension 1
            \param[in]    sizes1   number of elements in dimension 2 as function of the index in dimension 1
            \param[in]    sizes2   number of elements in dimension 3 as function of the indices in the lower dimensions
            \param[in]    sizes3   number of elements in dimension 4 as function of the indices in the lower dimensions
            \param[in]    ini      initial value for the array elements
         */
        template <class AForeign = allocator_size_type>
        IrregArray4D(const index_type x1, const index_type x2, const IrregArray1D<size_type, AForeign> &sizes1, const IrregArray2D<size_type, AForeign> &sizes2, const IrregArray3D<size_type, AForeign> &sizes3, const value_type &ini)
            : first1_(0), first2_(0), first3_(0), first4_(0), last1_(-1), last2_(-1), last3_(-1), last4_(-1), isIrreg_(true), arr_(nullptr)
        {
            initArray(x1, x2, sizes1, sizes2, sizes3, ini);
        }

        /*! \brief   general constructor helper function for irregularly shaped matrices and all special cases

            \tparam      AForeign   allocator used by the input data structure

            \param[in]    sizes3   number of elements in dimension 4 as function of the indices in the lower dimensions
         */
        template <class AForeign = allocator_size_type>
        IrregArray4D(const IrregArray3D<size_type, AForeign> &sizes3)
            : first1_(0), first2_(0), first3_(0), first4_(0), last1_(-1), last2_(-1), last3_(-1), last4_(-1), isIrreg_(true), arr_(nullptr)
        {
            initArray(sizes3);
        }
        /*! \brief   general constructor helper function for irregularly shaped matrices and all special cases

            \tparam      AForeign   allocator used by the input data structure

            \param[in]    sizes3   number of elements in dimension 4 as function of the indices in the lower dimensions
            \param[in]    ini      initial value for the array elements
         */
        template <class AForeign = allocator_size_type>
        IrregArray4D(const IrregArray3D<size_type, AForeign> &sizes3, const value_type &ini)
            : first1_(0), first2_(0), first3_(0), first4_(0), last1_(-1), last2_(-1), last3_(-1), last4_(-1), isIrreg_(true), arr_(nullptr)
        {
            initArray(sizes3, ini);
        }

        /*! \brief   specialized constructor for irregular, symmetrically shaped matrices
                     application: interaction energy matrix site1 form1(site1) site2 form2(site2)

            \tparam      AForeign   allocator used by the input data structure

            \param[in]    x1       starting index in dimensions 1 and 2
            \param[in]    x2       last index in dimension 1
            \param[in]    sizes1   number of elements in dimension 2 as function of the index in dimension 1
            \param[in]    z1       starting index in dimension 3
            \param[in]    z2       last index in dimension 2
            \param[in]    sizes2   number of elements in dimension 3 as function of the indices in the lower dimensions
            \param[in]    ini      initial value for the array elements
         */
        template <class AForeign = allocator_size_type>
        IrregArray4D(const index_type x1, const index_type x2, const IrregArray1D<size_type, AForeign> &sizes1, const index_type z1, const index_type z2, const IrregArray1D<size_type, AForeign> &sizes2, const value_type &ini)
            : first1_(x1), first2_(x1), first3_(z1), first4_(z1), last1_(x2), last2_(0), last3_(z2), last4_(0), isIrreg_(true), arr_(nullptr)
        {
            initArray(x1, x2, sizes1, z1, z2, sizes2, ini);
        }
        /*! \brief   specialized constructor for irregular, symmetrically shaped matrices
                     application: interaction energy matrix site1 form1(site1) site2 form2(site2)

            \tparam      AForeign   allocator used by the input data structure

            \param[in]    x1       starting index in dimensions 1 and 2
            \param[in]    x2       last index in dimension 1
            \param[in]    sizes1   number of elements in dimension 2 as function of the index in dimension 1
            \param[in]    z1       starting index in dimension 3
            \param[in]    z2       last index in dimension 2
            \param[in]    sizes2   number of elements in dimension 3 as function of the indices in the lower dimensions
         */
        template <class AForeign = allocator_size_type>
        IrregArray4D(const index_type x1, const index_type x2, const IrregArray1D<size_type, AForeign> &sizes1, const index_type z1, const index_type z2, const IrregArray1D<size_type, AForeign> &sizes2)
            : first1_(x1), first2_(x1), first3_(z1), first4_(z1), last1_(x2), last2_(0), last3_(z2), last4_(0), isIrreg_(true), arr_(nullptr)
        {
            initArray(x1, x2, sizes1, z1, z2, sizes2);
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

            \param[in]   ini  initialization value
         */
        IrregArray4D(const std::initializer_list<std::initializer_list<std::initializer_list<std::initializer_list<value_type> > > > &ini)
            : first1_(0), first2_(0), first3_(0), first4_(0), last1_(-1), last2_(-1), last3_(-1), last4_(-1), isIrreg_(true), arr_(nullptr)
        {
            typedef std::initializer_list<value_type>         l1_ilist_type;
            typedef std::initializer_list<l1_ilist_type>      l2_ilist_type;
            typedef std::initializer_list<l2_ilist_type>      l3_ilist_type;
            typedef const value_type*                     l1_ilist_ptr_type;
            typedef const l1_ilist_type*                  l2_ilist_ptr_type;
            typedef const l2_ilist_type*                  l3_ilist_ptr_type;
            typedef const l3_ilist_type*                  l4_ilist_ptr_type;

            // check whether the list constains empty elements in any but the innermost nesting level
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
            size_type prev_length = (ini.begin())->size();
            for (l4_ilist_ptr_type iPtrL4 = ini.begin(); iPtrL4 != ini.end(); ++iPtrL4)
            {
                if (prev_length != iPtrL4->size())
                {
                    isIrreg_ = true;
                    break;
                }
                prev_length = iPtrL4->size();
            }
            if (!isIrreg_)
            {
                size_type prev_length = ini.begin()->begin()->size();
                for (l4_ilist_ptr_type iPtrL4 = ini.begin(); iPtrL4 != ini.end(); ++iPtrL4)
                {
                    for (l3_ilist_ptr_type iPtrL3 = iPtrL4->begin(); iPtrL3 != iPtrL4->end(); ++iPtrL3)
                    {
                        if (prev_length != iPtrL3->size())
                        {
                            isIrreg_ = true;
                            break;
                        }
                        prev_length = iPtrL3->size();
                    }
                    if (isIrreg_)
                    {
                        break;
                    }
                }
            }
            if (!isIrreg_)
            {
                size_type prev_length = ini.begin()->begin()->begin()->size();
                for (l4_ilist_ptr_type iPtrL4 = ini.begin(); iPtrL4 != ini.end(); ++iPtrL4)
                {
                    for (l3_ilist_ptr_type iPtrL3 = iPtrL4->begin(); iPtrL3 != iPtrL4->end(); ++iPtrL3)
                    {
                        for (l2_ilist_ptr_type iPtrL2 = iPtrL3->begin(); iPtrL2 != iPtrL3->end(); ++iPtrL2)
                        {
                            if (prev_length != iPtrL2->size())
                            {
                                isIrreg_ = true;
                                break;
                            }
                            prev_length = iPtrL2->size();
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
                IrregArray1D<size_type, allocator_size_type> irregArr1((index_type)0, (index_type)ini.size() - 1, (size_t)0);
                size_t i = 0;
                for (l4_ilist_ptr_type iPtrL4 = ini.begin(); iPtrL4 != ini.end(); ++iPtrL4)
                {
                    irregArr1(i) = iPtrL4->size();
                    i++;
                }

                IrregArray2D<size_type, allocator_size_type> irregArr2((index_type)0, (index_type)ini.size() - 1, irregArr1, (size_t)0);
                i = 0;
                for (l4_ilist_ptr_type iPtrL4 = ini.begin(); iPtrL4 != ini.end(); ++iPtrL4)
                {
                    size_t j = 0;
                    for (l3_ilist_ptr_type iPtrL3 = iPtrL4->begin(); iPtrL3 != iPtrL4->end(); ++iPtrL3)
                    {
                        irregArr2(i, j) = iPtrL3->size();
                        j++;
                    }
                    i++;
                }

                IrregArray3D<size_type, allocator_size_type> irregArr3((index_type)0, (index_type)ini.size() - 1, irregArr1, irregArr2, (size_t)0);
                i = 0;
                for (l4_ilist_ptr_type iPtrL4 = ini.begin(); iPtrL4 != ini.end(); ++iPtrL4)
                {
                    size_t j = 0;
                    for (l3_ilist_ptr_type iPtrL3 = iPtrL4->begin(); iPtrL3 != iPtrL4->end(); ++iPtrL3)
                    {
                        size_t k = 0;
                        for (l2_ilist_ptr_type iPtrL2 = iPtrL3->begin(); iPtrL2 != iPtrL3->end(); ++iPtrL2)
                        {
                            irregArr3(i, j, k) = iPtrL2->size();
                            k++;
                        }
                        j++;
                    }
                    i++;
                }

                initArray(0, ini.size() - 1, irregArr1, irregArr2, irregArr3);
            }
            else
            {
                initArray(0, static_cast<index_type>(ini.size()) - 1,
                          0, static_cast<index_type>((ini.begin())->size()) - 1,
                          0, static_cast<index_type>((ini.begin())->begin()->size()) - 1,
                          0, static_cast<index_type>((ini.begin())->begin()->begin()->size()) - 1);
            }

            // assign the values from the input initializer_list
            for (size_type l = 0; l < getLength1(); ++l)
            {
                l4_ilist_ptr_type iPtr4 = (ini.begin() + static_cast<std::ptrdiff_t>(l));
                for (size_type k = 0; k < getLength2(l); ++k)
                {
                    l3_ilist_ptr_type iPtr3 = (iPtr4->begin() + static_cast<std::ptrdiff_t>(k));
                    for (size_type j = 0; j < getLength3(l, k); ++j)
                    {
                        l2_ilist_ptr_type iPtr2 = (iPtr3->begin() + static_cast<std::ptrdiff_t>(j));
                        for (size_type i = 0; i < getLength4(l, k, j); ++i)
                        {
                            l1_ilist_ptr_type iPtr1 = (iPtr2->begin() + static_cast<std::ptrdiff_t>(i));
                            operator()(l, k, j, i) = *iPtr1;
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

            \param[in]    ini    array to be copied
         */
        template <typename TForeign, class AForeign = std::allocator<TForeign> >
        void initArray(const IrregArray4D<TForeign, AForeign> &ini)
        {
            if (ini.arr_ != nullptr)
            {
                if (ini.isIrreg_)
                {
                    initArray(ini.irreg_);
                }
                else
                {
                    initArray(ini.first1_, ini.last1_, ini.first2_, ini.last2_, ini.first3_, ini.last3_, ini.first4_, ini.last4_);
                }
            }

            if (arr_ != nullptr)
            {
                for (index_type h1 = first1_; h1 <= last1_; h1++)
                {
                    for (index_type h2 = first2_; h2 <= getLast2(h1); h2++)
                    {
                        for (index_type h3 = first3_; h3 <= getLast3(h1, h2); h3++)
                        {
                            for (index_type h4 = first4_; h4 <= getLast4(h1, h2, h3); h4++)
                            {
                                arr_[h1][h2][h3][h4] = static_cast<value_type>(ini(h1, h2, h3, h4));
                            }
                        }
                    }
                }
            }
        }
        /*! \brief specialized constructor helper function for regularly shaped matrices

            \param[in]    x1       starting index in dimension 1
            \param[in]    x2       last index in dimension 1
            \param[in]    y1       starting index in dimension 2
            \param[in]    y2       last index in dimension 2
            \param[in]    z1       starting index in dimension 3
            \param[in]    z2       last index in dimension 3
            \param[in]    k1       starting index in dimension 4
            \param[in]    k2       last index in dimension 4
         */
        void initArray(const index_type x1, const index_type x2, const index_type y1, const index_type y2, const index_type z1, const index_type z2, const index_type k1, const index_type k2)
        {
            deallocateArray();
            isIrreg_   = false;
            first1_    = x1;
            first2_    = y1;
            first3_    = z1;
            first4_    = k1;
            last1_     = x2;
            last2_     = y2;
            last3_     = z2;
            last4_     = k2;
            allocateArray();
        }
        /*! \brief specialized constructor helper function for regularly shaped matrices

            \param[in]    x1       starting index in dimension 1
            \param[in]    x2       last index in dimension 1
            \param[in]    y1       starting index in dimension 2
            \param[in]    y2       last index in dimension 2
            \param[in]    z1       starting index in dimension 3
            \param[in]    z2       last index in dimension 3
            \param[in]    k1       starting index in dimension 4
            \param[in]    k2       last index in dimension 4
            \param[in]    ini      initialization value for the data array elements
         */
        void initArray(const index_type x1, const index_type x2, const index_type y1, const index_type y2, const index_type z1, const index_type z2, const index_type k1, const index_type k2, const value_type &ini)
        {
            initArray(x1, x2, y1, y2, z1, z2, k1, k2);
            for (index_type h1 = x1; h1 <= x2; h1++)
            {
                for (index_type h2 = y1; h2 <= y2; h2++)
                {
                    for (index_type h3 = z1; h3 <= z2; h3++)
                    {
                        for (index_type h4 = k1; h4 <= k2; h4++)
                        {
                            arr_[h1][h2][h3][h4] = ini;
                        }
                    }
                }
            }
        }
        /*! \brief general constructor helper function for irregularly shaped matrices and all special cases

            \tparam      AForeign   allocator used by the input data structure

            \param[in]    x1       starting index in dimension 1
            \param[in]    y1       starting index in dimension 2
            \param[in]    z1       starting index in dimension 3
            \param[in]    k1       starting index in dimension 4
            \param[in]    sizes1   number of elements in dimension 2 as function of the index in dimension 1
            \param[in]    sizes2   number of elements in dimension 3 as function of the indices in the lower dimensions
            \param[in]    sizes3   number of elements in dimension 4 as function of the indices in the lower dimensions
         */
        template<class AForeign = std::allocator<value_type> >
        void initArray(const index_type x1, const index_type y1, const index_type z1, const index_type k1, const IrregArray1D<size_type, AForeign> &sizes1, const IrregArray2D<size_type, AForeign> &sizes2, const IrregArray3D<size_type, AForeign> &sizes3)
        {
            deallocateArray();
            isIrreg_   = true;
            first1_    = x1;
            first2_    = y1;
            first3_    = z1;
            first4_    = k1;
            last1_     = first1_ + sizes1.getLength1() - 1;
            last2_     = 0;
            last3_     = 0;
            last4_     = 0;
            irreg_.initArray(first1_, last1_, sizes1, sizes2);
            for (index_type i = getFirst1(); i <= getLast1(); i++)
            {
                const index_type end2 = sizes1(i) + getFirst2();
                for (index_type j = getFirst2(); j < end2; j++)
                {
                    const index_type end3 = getFirst3() + sizes2(i, j);
                    for (index_type k = getFirst3(); k < end3; k++)
                    {
                        irreg_(i, j, k) = sizes3(i, j, k);
                    }
                }
            }
            allocateArray();
        }
        /*! \brief general constructor helper function for irregularly shaped matrices and all special cases

            \tparam      AForeign   allocator used by the input data structure

            \param[in]    sizes3   number of elements in dimension 4 as function of the indices in the lower dimensions
         */
        template<class AForeign = std::allocator<value_type> >
        void initArray(const IrregArray3D<size_type, AForeign> &sizes3)
        {
            deallocateArray();
            isIrreg_   = true;
            first1_    = sizes3.getFirst1();
            first2_    = sizes3.getFirst1();
            first3_    = sizes3.getFirst1();
            first4_    = sizes3.getFirst1();
            last1_     = sizes3.getLast1();
            last2_     = 0;
            last3_     = 0;
            last4_     = 0;
            irreg_.initArray(sizes3.irreg_);
            for (index_type i = getFirst1(); i <= getLast1(); i++)
            {
                const index_type end2 = getFirst2() + sizes3.getLength2(i);
                for (index_type j = getFirst2(); j < end2; j++)
                {
                    const index_type end3 = getFirst3() + sizes3.getLength3(i, j);
                    for (index_type k = getFirst3(); k < end3; k++)
                    {
                        irreg_(i, j, k) = sizes3(i, j, k);
                    }
                }
            }
            allocateArray();
        }

        /*! \brief general constructor helper function for irregularly shaped matrices and all special cases

            \tparam      AForeign   allocator used by the input data structure

            \param[in]    sizes3   number of elements in dimension 4 as function of the indices in the lower dimensions
            \param[in]    ini      initialization value for the data array elements
         */
        template<class AForeign = allocator_size_type>
        void initArray(const IrregArray3D<size_type, AForeign> &sizes3, const value_type &ini)
        {
            initArray(sizes3);
            for (index_type h1 = getFirst1(); h1 <= getLast1(); h1++)
            {
                for (index_type h2 = getFirst2(); h2 <= getLast2(h1); h2++)
                {
                    for (index_type h3 = getFirst3(); h3 <= getLast3(h1, h2); h3++)
                    {
                        for (index_type h4 = getFirst4(); h4 <= getLast4(h1, h2, h3); h4++)
                        {
                            arr_[h1][h2][h3][h4] = ini;
                        }
                    }
                }
            }
        }

        /*! \brief general constructor helper function for irregularly shaped matrices and all special cases

            \tparam      AForeign   allocator used by the input data structure

            \param[in]    x1       starting index in dimension 1
            \param[in]    y1       starting index in dimension 2
            \param[in]    z1       starting index in dimension 3
            \param[in]    k1       starting index in dimension 4
            \param[in]    sizes1   number of elements in dimension 2 as function of the index in dimension 1
            \param[in]    sizes2   number of elements in dimension 3 as function of the indices in the lower dimensions
            \param[in]    sizes3   number of elements in dimension 4 as function of the indices in the lower dimensions
            \param[in]    ini      initialization value for the data array elements
         */
        template<class AForeign = allocator_size_type>
        void initArray(const index_type x1, const index_type y1, const index_type z1, const index_type k1, const IrregArray1D<size_type, AForeign> &sizes1, const IrregArray2D<size_type, AForeign> &sizes2, const IrregArray3D<size_type, AForeign> &sizes3, const value_type &ini)
        {
            initArray(x1, y1, z1, k1, sizes1, sizes2, sizes3);
            for (index_type h1 = x1; h1 <= getLast1(); h1++)
            {
                for (index_type h2 = y1; h2 <= getLast2(h1); h2++)
                {
                    for (index_type h3 = z1; h3 <= getLast3(h1, h2); h3++)
                    {
                        for (index_type h4 = k1; h4 <= getLast4(h1, h2, h3); h4++)
                        {
                            arr_[h1][h2][h3][h4] = ini;
                        }
                    }
                }
            }
        }

        /*! \brief general constructor helper function for irregularly shaped matrices and all special cases

            \tparam      AForeign   allocator used by the input data structure

            \param[in]    x1       starting index in dimension 1
            \param[in]    x2       last     index in dimension 1
            \param[in]    sizes1   number of elements in dimension 2 as function of the index in dimension 1
            \param[in]    sizes2   number of elements in dimension 3 as function of the indices in the lower dimensions
            \param[in]    sizes3   number of elements in dimension 4 as function of the indices in the lower dimensions
         */
        template<class AForeign = allocator_size_type>
        void initArray(const index_type x1, const index_type x2, const IrregArray1D<size_type, AForeign> &sizes1, const IrregArray2D<size_type, AForeign> &sizes2, const IrregArray3D<size_type, AForeign> &sizes3)
        {
            deallocateArray();
            isIrreg_   = true;
            first1_    = x1;
            first2_    = x1;
            first3_    = x1;
            first4_    = x1;
            last1_     = x2;
            last2_     = 0;
            last3_     = 0;
            last4_     = 0;
            // make sure that irreg_ matches the range prescribed by x1 and x2
            IrregArray1D<size_type, allocator_size_type> tmpIrreg(first1_, last1_, (size_type)0);
            for (index_type i = first1_; i <= last1_; ++i)
            {
                tmpIrreg(i) = sizes1(i);
            }
            initArray(first1_, first2_, first3_, first4_, tmpIrreg, sizes2, sizes3);
        }
        /*! \brief general constructor helper function for irregularly shaped matrices and all special cases

            \tparam      AForeign   allocator used by the input data structure

            \param[in]    x1       starting index in dimension 1
            \param[in]    x2       last index in dimension 1
            \param[in]    sizes1   number of elements in dimension 2 as function of the index in dimension 1
            \param[in]    sizes2   number of elements in dimension 3 as function of the indices in the lower dimensions
            \param[in]    sizes3   number of elements in dimension 4 as function of the indices in the lower dimensions
            \param[in]    ini      initialization value for the data array elements
         */
        template<class AForeign = allocator_size_type>
        void initArray(const index_type x1, const index_type x2, const IrregArray1D<size_type, AForeign> &sizes1, const IrregArray2D<size_type, AForeign> &sizes2, const IrregArray3D<size_type, AForeign> &sizes3, const value_type &ini)
        {
            initArray(x1, x2, sizes1, sizes2, sizes3);
            for (index_type h1 = getFirst1(); h1 <= getLast1(); h1++)
            {
                for (index_type h2 = getFirst2(); h2 <= getLast2(h1); h2++)
                {
                    for (index_type h3 = getFirst3(); h3 <= getLast3(h1, h2); h3++)
                    {
                        for (index_type h4 = getFirst4(); h4 <= getLast4(h1, h2, h3); h4++)
                        {
                            arr_[h1][h2][h3][h4] = ini;
                        }
                    }
                }
            }
        }

        /*! \brief general constructor helper function for irregularly shaped matrices and all special cases

            \param[in]    x1       starting index in dimension 1
            \param[in]    x2       last index in dimension 1
            \param[in]    sizes1   number of elements in dimension 2 as function of the index in dimension 1
            \param[in]    sizes2   number of elements in dimension 3 as function of the indices in the lower dimensions
            \param[in]    sizes3   number of elements in dimension 4 as function of the indices in the lower dimensions
         */
        template<class AForeign = allocator_size_type>
        void initArray(const index_type x1, const index_type x2, const size_type *sizes1, const IrregArray2D<size_type, AForeign> &sizes2, const IrregArray3D<size_type, AForeign> &sizes3)
        {
            initArray(x1, x2, IrregArray1D<size_type, allocator_size_type>(x1, x2, sizes1), sizes2, sizes3);
        }
        /*! \brief specialized constructor helper function for irregular, symmetrically shaped matrices
                   application: interaction energy matrix site1 form1(site1) site2 form2(site2)

            \tparam      AForeign   allocator used by the input data structure

            \param[in]    x1       starting index in dimension 1
            \param[in]    x2       last index in dimension 1
            \param[in]    sizes1   number of index values in dimension 4
                                   as function of the indices in dimension 3
                                   (in this special case, independent on the
                                   index in dimensions 1 and 2)
            \param[in]    z1       starting index in dimension 3
            \param[in]    z2       last index in dimension 3
            \param[in]    sizes2   number of index values in dimension 4
                                   as function of the indices in dimension 3
                                   (in this special case, independent on the
                                   index in dimensions 1 and 2)
         */
        template<class AForeign = allocator_size_type>
        void initArray(const index_type x1, const index_type x2, const IrregArray1D<size_type, AForeign> &sizes1, const index_type z1, const index_type z2, const IrregArray1D<size_type, AForeign> &sizes2)
        {
            deallocateArray();
            isIrreg_   = true;
            first1_    = x1;
            first2_    = first1_;
            first3_    = z1;
            first4_    = first3_;
            last1_     = x2;
            last2_     = 0;
            last3_     = z2;
            last4_     = 0;
            IrregArray1D<size_type, allocator_size_type> irregVec(first1_, last1_);
            for (index_type i = first1_; i <= last1_; i++)
            {
                irregVec[i] = sizes1[i];
            }
            IrregArray2D<size_type, allocator_size_type> irregMat(first1_, last1_, sizes1, (size_type)0);
            const index_type        dim3 = z2 - z1 + 1;
            for (index_type i = first1_; i <= last1_; i++)
            {
                for (index_type j = first2_; j <= getLast2(i); j++)
                {
                    irregMat(i, j) = dim3;
                }
            }
            irreg_.initArray(first1_, last1_, irregVec, irregMat, (size_type)0);
            for (index_type i = first1_; i <= last1_; i++)
            {
                for (index_type j = first2_; j <= getLast2(i); j++)
                {
                    for (index_type k = first3_; k <= getLast3(i, j); k++)
                    {
                        irreg_(i, j, k) = sizes2[k];
                    }
                }
            }
            allocateArray();
        }
        /*! \brief specialized constructor helper function for irregular, symmetrically shaped matrices
                   application: interaction energy matrix site1 form1(site1) site2 form2(site2)

            \tparam      AForeign   allocator used by the input data structure

            \param[in]    x1       starting index in dimension 1
            \param[in]    x2       last index in dimension 1
            \param[in]    sizes1   number of index values in dimension 4
                                   as function of the indices in dimension 3
                                   (in this special case, independent on the
                                   index in dimensions 1 and 2)
            \param[in]    z1       starting index in dimension 3
            \param[in]    z2       last index in dimension 3
            \param[in]    sizes2   number of index values in dimension 4
                                   as function of the indices in dimension 3
                                   (in this special case, independent on the
                                   index in dimensions 1 and 2)
            \param[in]    ini      initialization value for the data array elements
         */
        template<class AForeign = allocator_size_type>
        void initArray(const index_type x1, const index_type x2, const IrregArray1D<size_type, AForeign> &sizes1, const index_type z1, const index_type z2, const IrregArray1D<size_type, AForeign> &sizes2, const value_type &ini)
        {
            initArray(x1, x2, sizes1, z1, z2, sizes2);
            for (index_type h1 = getFirst1(); h1 <= getLast1(); h1++)
            {
                for (index_type h2 = getFirst2(); h2 <= getLast2(h1); h2++)
                {
                    for (index_type h3 = getFirst3(); h3 <= getLast3(h1, h2); h3++)
                    {
                        for (index_type h4 = getFirst4(); h4 <= getLast4(h1, h2, h3); h4++)
                        {
                            arr_[h1][h2][h3][h4] = ini;
                        }
                    }
                }
            }
        }

        //! returns true if the array is irregular
        inline bool isIrreg() const { return isIrreg_; }

        //! number of elements in dimension 1
        size_type getLength1() const { return (last1_ - first1_ + 1); }
        //! number of elements in dimension 2 (regular array)
        size_type getLength2() const
        {
            if (!isIrreg_)
            {
                return (last2_ - first2_ + 1);
            }
            else
            {
                std::string errorMessage(formatString("Error in %s::getLength2(): called without argument for an irregular array.", typeid(*this).name()));
                GMX_THROW(APIError(errorMessage));
            }
        }
        /*! \brief number of elements in dimension 2 at lower dimension index x

            \param[in]    x       index in dimension 1 for which the array length in dimension 2 is requested */
        size_type getLength2(const index_type x) const
        {
            return (isIrreg_ ? irreg_.getLength2(x) : (last2_ - first2_ + 1));
        }
        //! number of elements in dimension 3 (regular array)
        size_type getLength3() const
        {
            if (!isIrreg_)
            {
                return (last3_ - first3_ + 1);
            }
            else
            {
                std::string errorMessage(formatString("Error in %s::getLength3(): called without argument for an irregular array.", typeid(*this).name()));
                GMX_THROW(APIError(errorMessage));
            }
        }
        /*! \brief number of elements in dimension 3 at lower dimension indices x, y

            \param[in]    x       index in dimension 1
            \param[in]    y       index in dimension 2  */
        size_type getLength3(const index_type x, const index_type y) const
        {
            return (isIrreg_ ? irreg_.getLength3(x, y) : (last3_ - first3_ + 1));
        }
        //! number of elements in dimension 4 (regular array)
        size_type getLength4() const
        {
            if (!isIrreg_)
            {
                return (last4_ - first4_ + 1);
            }
            else
            {
                std::string errorMessage(formatString("Error in %s::getLength4(): called without argument for an irregular array.", typeid(*this).name()));
                GMX_THROW(APIError(errorMessage));
            }
        }
        /*! \brief number of elements in dimension 4 at lower dimension indices x, y
            \param[in]    x       index in dimension 1
            \param[in]    y       index in dimension 2
            \param[in]    z       index in dimension 3  */
        size_type getLength4(const index_type x, const index_type y, const index_type z) const
        {
            return (isIrreg_ ? irreg_(x, y, z) : (last4_ - first4_ + 1));
        }

        //! get the first index of dimension 1
        index_type getBegin1() const { return first1_; }
        //! get the first index of dimension 2
        index_type getBegin2() const { return first2_; }
        //! get the first index of dimension 3
        index_type getBegin3() const { return first3_; }
        //! get the first index of dimension 4
        index_type getBegin4() const { return first4_; }

        //! get the first index of dimension 1
        index_type getFirst1() const { return first1_; }
        //! get the first index of dimension 2
        index_type getFirst2() const { return first2_; }
        //! get the first index of dimension 3
        index_type getFirst3() const { return first3_; }
        //! get the first index of dimension 4
        index_type getFirst4() const { return first4_; }

        //! get the last index of dimension 1
        index_type getLast1() const { return last1_; }
        //! get the last index of dimension 2
        index_type getLast2() const
        {
            if (!isIrreg_)
            {
                return last2_;
            }
            else
            {
                std::string errorMessage(formatString("Error in %s::getLast2(): called without argument for an irregular array.", typeid(*this).name()));
                GMX_THROW(APIError(errorMessage));
            }
        }
        /*! \brief get the last index of dimension 2 at lower dimension index x

            \param[in]    x       index in dimension 1 for which the last index in dimension 2 is requested */
        index_type getLast2(const index_type x) const
        {
            return (isIrreg_ ? irreg_.getLast2(x) : last2_);
        }
        //! get the last index of dimension 3
        index_type getLast3() const
        {
            if (!isIrreg_)
            {
                return last3_;
            }
            else
            {
                std::string errorMessage(formatString("Error in %s::getLast3(): called without argument for an irregular array.", typeid(*this).name()));
                GMX_THROW(APIError(errorMessage));
            }
        }
        /*! \brief get the last index of dimension 3 for lower dimension indices x, y

            \param[in]    x       index in dimension 1
            \param[in]    y       index in dimension 2                       */
        index_type getLast3(const index_type x, const index_type y) const
        {
            return (isIrreg_ ? irreg_.getLast3(x, y) : last3_);
        }
        //! get the last index of dimension 4
        index_type getLast4() const
        {
            if (!isIrreg_)
            {
                return last4_;
            }
            else
            {
                std::string errorMessage(formatString("Error in %s::getLast4(): called without argument for an irregular array.", typeid(*this).name()));
                GMX_THROW(APIError(errorMessage));
            }
        }
        /*! \brief get the last index of dimension 3 for lower dimension indices x, y, z

            \param[in]    x       index in dimension 1
            \param[in]    y       index in dimension 2
            \param[in]    z       index in dimension 3                       */
        index_type getLast4(const index_type x, const index_type y, const index_type z) const
        {
            return (isIrreg_ ? (irreg_(x, y, z) + first4_ - 1) : last4_);
        }

        //! get the index one past the last valid index of dimension 1
        index_type getEnd1() const { return last1_; }
        //! get the index one past the last valid index of dimension 2
        index_type getEnd2() const
        {
            if (!isIrreg_)
            {
                return last2_ + 1;
            }
            else
            {
                std::string errorMessage(formatString("Error in %s::getEnd2(): called without argument for an irregular array.", typeid(*this).name()));
                GMX_THROW(APIError(errorMessage));
            }
        }
        /*! \brief get the index one past the last valid index of dimension 2 at lower dimension index x

            \param[in]    x       index in dimension 1 for which the last index in dimension 2 is requested */
        index_type getEnd2(const index_type x) const
        {
            return (isIrreg_ ? irreg_.getEnd2(x) : (last2_ + 1));
        }
        //! get the index one past the last valid index of dimension 3
        index_type getEnd3() const
        {
            if (!isIrreg_)
            {
                return last3_ + 1;
            }
            else
            {
                std::string errorMessage(formatString("Error in %s::getEnd3(): called without argument for an irregular array.", typeid(*this).name()));
                GMX_THROW(APIError(errorMessage));
            }
        }
        /*! \brief get the index one past the last valid index of dimension 3 for lower dimension indices x, y

            \param[in]    x       index in dimension 1
            \param[in]    y       index in dimension 2                       */
        index_type getEnd3(const index_type x, const index_type y) const
        {
            return (isIrreg_ ? irreg_.getEnd3(x, y) : (last3_ + 1));
        }
        //! get the index one past the last valid index of dimension 4
        index_type getEnd4() const
        {
            if (!isIrreg_)
            {
                return last4_ + 1;
            }
            else
            {
                std::string errorMessage(formatString("Error in %s::getEnd4(): called without argument for an irregular array.", typeid(*this).name()));
                GMX_THROW(APIError(errorMessage));
            }
        }
        /*! \brief get the index one past the last valid index of dimension 3 for lower dimension indices x, y, z

            \param[in]    x       index in dimension 1
            \param[in]    y       index in dimension 2
            \param[in]    z       index in dimension 3                       */
        index_type getEnd4(const index_type x, const index_type y, const index_type z) const
        {
            return (isIrreg_ ? (irreg_(x, y, z) + first4_) : (last4_ + 1));
        }

        //! returns the amount of memory occupied by the object
        size_type getSize() const
        {
            // the book-keeping data
            size_type size = sizeof(*this);
            size += irreg_.getSize() - sizeof(size_3d_array_type);
            // the storage data
            size += getLength1() * sizeof(value_type***);
            for (index_type i = getFirst1(); i <= getLast1(); i++)
            {
                size += getLength2(i) * sizeof(value_type**);
                for (index_type j = getFirst2(); j <= getLast2(i); j++)
                {
                    size += getLength3(i, j) * sizeof(value_type*);
                    for (index_type k = getFirst3(); k <= getLast3(i, j); k++)
                    {
                        size += getLength4(i, j, k) * sizeof(value_type);
                    }
                }
            }
            return size;
        }

        //! determine the number of data array elements
        size_type getNelements() const
        {
            // the storage data
            size_type nelements_ = 0;
            for (index_type i = getFirst1(); i <= getLast1(); i++)
            {
                for (index_type j = getFirst2(); j <= getLast2(i); j++)
                {
                    for (index_type k = getFirst3(); k <= getLast3(i, j); k++)
                    {
                        nelements_ += getLength4(i, j, k);
                    }
                }
            }
            return nelements_;
        }

        //! deallocate the memory of all constituent arrays
        void deallocateArray()
        {
            if (arr_ != nullptr)
            {
                // destroy the array elements
                for (index_type i = getFirst1(); i <= getLast1(); ++i)
                {
                    for (index_type j = getFirst2(); j <= getLast2(i); ++j)
                    {
                        for (index_type k = getFirst3(); k <= getLast3(i, j); ++k)
                        {
                            for (index_type l = getFirst4(); l <= getLast4(i, j, k); ++l)
                            {
                                std::allocator_traits<allocator_type>::destroy(allocator_, &this->operator()(i, j, k, l));
                            }
                        }
                    }
                }

                // deallocate data array
                for (index_type i = getFirst1(); i <= getLast1(); i++)
                {
                    for (index_type j = getFirst2(); j <= getLast2(i); j++)
                    {
                        for (index_type k = getFirst3(); k <= getLast3(i, j); k++)
                        {
                            arr_[i][j][k] += first4_;
                            std::allocator_traits<allocator_type>::deallocate(allocator_, arr_[i][j][k], getLength4(i, j, k));
                        }
                        arr_[i][j] += first3_;
                        std::allocator_traits<p_allocator_type>::deallocate(p_allocator_, arr_[i][j], getLength3(i, j));
                    }
                    arr_[i] += first2_;
                    std::allocator_traits<p2_allocator_type>::deallocate(p2_allocator_, arr_[i], getLength2(i));
                }
                arr_ += first1_;
                std::allocator_traits<p3_allocator_type>::deallocate(p3_allocator_, arr_, getLength1());
                arr_ = nullptr;

/* if block allocation was used
                else
                {
                    const size_type d1 = getLength1();
                    const size_type d2 = getLength2();
                    const size_type d3 = getLength3();
                    const size_type d4 = getLength4();
                    arr_[first1_][first2_][first3_] += first4_;
                    std::allocator_traits<allocator_type>::deallocate(allocator_, arr_[first1_][first2_][first3_], d1*d2*d3*d4);
                    arr_[first1_][first2_] += first3_;
                    std::allocator_traits<p_allocator_type>::deallocate(p_allocator_, arr_[first1_][first2_], d1*d2*d3);
                    arr_[first1_] += first2_;
                    std::allocator_traits<p2_allocator_type>::deallocate(p2_allocator_, arr_[first1_], d1*d2);
                    arr_ += first1_;
                    std::allocator_traits<p3_allocator_type>::deallocate(p3_allocator_, arr_, d1);
                }
 */

            }
            // deallocate book-keeping arrays
            if (isIrreg_)
            {
                irreg_.deallocateArray();
                isIrreg_ = false;
            }
            first1_ = 0;
            first2_ = 0;
            first3_ = 0;
            first4_ = 0;
            last1_  = -1;
            last2_  = -1;
            last3_  = -1;
            last4_  = -1;
        }

        //! destructor
        ~IrregArray4D()
        {
            // deallocate arrays
            deallocateArray();
        }

        // index operators

        /*! index operator for bare array-like usage
            IrregArray2d<value_type> a;
            value_type  val   = a[i][j][k]
         */

        //! doxygen can not handle nested classes
        /// \cond DEV
        /*! \class Proxy

            \brief helper functor used to implement index operator *this[i][j][k][l]
         */
        class Proxy
        {
            private:
                value_type*** array_;
                class Proxy2
                {
                    private:
                        value_type** array_;
                        class Proxy3
                        {
                            private:
                                value_type* array_;
                            public:
                                Proxy3(value_type* array) : array_(array) { }

                                value_type &operator[](index_type k)
                                {
                                    return array_[k];
                                }
                                const value_type &operator[](index_type k) const
                                {
                                    return array_[k];
                                }
                        };
                    public:
                        Proxy2(value_type** array) : array_(array) { }

                        Proxy3 operator[](index_type z)
                        {
                            return Proxy3(array_[z]);
                        }
                        const Proxy3 operator[](index_type z) const
                        {
                            return array_[z];
                        }
                };
                //! let IrregArray4D access private members
                template<class W, typename AW> friend class IrregArray4D;
            public:
                Proxy(value_type*** array) : array_(array) { }

                Proxy2 operator[](index_type y)
                {
                    return Proxy2(array_[y]);
                }
                const Proxy2 operator[](index_type y) const
                {
                    return Proxy2(array_[y]);
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
        Proxy operator[](index_type x)
        {
            return Proxy(arr_[x]);
        }
        /*! \brief   const index operator[x][y][z][k] for array-like usage realized with chained functors "Proxy*"
            \param[in]   x   index in dimension 1
            \cond
            \param[in]   y   index in dimension 2
            \param[in]   z   index in dimension 3
            \param[in]   k   index in dimension 4
            \endcond
         */
        const Proxy operator[](index_type x) const
        {
            return Proxy(arr_[x]);
        }

        /*! \brief   index operator(x,y,z) for functor-like usage
            \param[in]   x   index in dimension 1
            \param[in]   y   index in dimension 2
            \param[in]   z   index in dimension 3
            \param[in]   k   index in dimension 4
         */
        inline value_type &operator()(const index_type x, const index_type y, const index_type z, const index_type k)
        {
            return arr_[x][y][z][k];
        }
        /*! \brief   const index operator(x,y,z) for functor-like usage
            \param[in]   x   index in dimension 1
            \param[in]   y   index in dimension 2
            \param[in]   z   index in dimension 3
            \param[in]   k   index in dimension 4
         */
        inline const value_type &operator()(const index_type x, const index_type y, const index_type z, const index_type k) const
        {
            return arr_[x][y][z][k];
        }
        /*! \brief   index operator(x,y,z) for functor-like usage
            \param[in]   x   index in dimension 1
            \param[in]   y   index in dimension 2
            \param[in]   z   index in dimension 3
         */
        inline value_type *operator()(const index_type x, const index_type y, const index_type z)
        {
            return arr_[x][y][z];
        }
        /*! \brief   const index operator(x,y,z) for functor-like usage
            \param[in]   x   index in dimension 1
            \param[in]   y   index in dimension 2
            \param[in]   z   index in dimension 3
         */
        inline const value_type *operator()(const index_type x, const index_type y, const index_type z) const
        {
            return arr_[x][y][z];
        }
        /*! \brief   index operator(x,y,z) for functor-like usage
            \param[in]   x   index in dimension 1
            \param[in]   y   index in dimension 2
         */
        inline value_type **operator()(const index_type x, const index_type y)
        {
            return arr_[x][y];
        }
        /*! \brief   const index operator(x,y,z) for functor-like usage
            \param[in]   x   index in dimension 1
            \param[in]   y   index in dimension 2
         */
        inline const value_type **operator()(const index_type x, const index_type y) const
        {
            return arr_[x][y];
        }
        /*! \brief   index operator(x,y,z) for functor-like usage
            \param[in]   x   index in dimension 1
         */
        inline value_type ***operator()(const index_type x)
        {
            return arr_[x];
        }
        /*! \brief   const index operator(x,y,z) for functor-like usage
            \param[in]   x   index in dimension 1
         */
        inline const value_type ***operator()(const index_type x) const
        {
            return arr_[x];
        }

        //! get a pointer to the data array
        const value_type**** getArray() const
        {
            return const_cast<const value_type****>(arr_);
        }
        //! get a pointer to the data array
        value_type**** getArray()
        {
            return arr_;
        }

        //! get a pointer to the data array
        const value_type**** data() const
        {
            return const_cast<const value_type****>(arr_);
        }
        //! get a pointer to the data array
        value_type**** data()
        {
            return arr_;
        }

        /*! \brief   stream operator for convenient printing of the array to an ostream
            \param[in]   output   output stream for the array contents
            \param[in]   ten      the array to be printed
         */
        friend std::ostream &operator<<(std::ostream &output, const IrregArray4D &ten)
        {
            // safe the current ostream format for restoring it after output insertion
            std::ios  state(NULL);
            state.copyfmt(output);

            if (ten.arr_ != nullptr)
            {
                for (index_type i = ten.first1_; i <= ten.last1_; i++)
                {
                    for (index_type j = ten.first2_; j <= ten.getLast2(i); j++)
                    {
                        output << "------------------------------" << std::endl << i << ", " << j << std::endl;
                        for (index_type k = ten.first3_; k <= ten.getLast3(i, j); k++)
                        {
                            for (index_type l = ten.first4_; l <= ten.getLast4(i, j, k); l++)
                            {
                                if (std::is_floating_point<value_type>::value)
                                {
                                    output << std::setw(10) << std::fixed << std::setprecision(2) << ten(i, j, k, l);
                                }
                                else
                                {
                                    output << ten(i, j, k, l);
                                }
                                if (k != ten.getLast4(i, j, k))
                                {
                                    output << " ";
                                }
                            }
                            output << std::endl;
                        }
                        output << std::endl;
                    }
                    output << std::endl;
                }
                output << std::endl;
            }

            // restore the original ostream format
            output.copyfmt(state);

            return output;
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
                first1_       = rhs.first1_;
                first2_       = rhs.first2_;
                first3_       = rhs.first3_;
                first4_       = rhs.first4_;
                last1_        = rhs.last1_;
                last2_        = rhs.last2_;
                last3_        = rhs.last3_;
                last4_        = rhs.last4_;
                isIrreg_      = rhs.isIrreg_;
                arr_          = rhs.arr_;
                irreg_        = std::move(rhs.irreg_);
                allocator_    = std::move(rhs.allocator_);
                p_allocator_  = std::move(rhs.p_allocator_);
                p2_allocator_ = std::move(rhs.p2_allocator_);
                p3_allocator_ = std::move(rhs.p3_allocator_);
                // set rhs to a valid default state
                rhs.first1_   =  0;
                rhs.first2_   =  0;
                rhs.first3_   =  0;
                rhs.first4_   =  0;
                rhs.last1_    = -1;
                rhs.last2_    = -1;
                rhs.last3_    = -1;
                rhs.last4_    = -1;
                rhs.isIrreg_  = false;
                rhs.arr_      = nullptr;
            }
            return *this;
        }

        /*! \brief   assignment operator *this = rhs
            \param[in]   rhs   right-hand side object of the assignment
         */
        IrregArray4D &operator=(const IrregArray4D &rhs)
        {
            deallocateArray();
            initArray(rhs);
            return *this;
        }

        /*! \brief   assignment operator *this = rhs

            \tparam      AForeign   allocator used by the input data structure

            \param[in]   rhs   right-hand side object of the assignment
         */
        template<class AForeign = allocator_size_type>
        IrregArray4D &operator=(const IrregArray4D<value_type, AForeign> &rhs)
        {
            deallocateArray();
            initArray(rhs);
            return *this;
        }

        /*! \brief  conversion operator (explicit to avoid unintended implicit conversion), usable via
                    IrregArray4D<TForeign, AForeign> newArr = static_cast<IrregArray4D<TForeign, AForeign> >(IrregArray4D<TCurrent, ACurrent> oldArr);

            \tparam      TForeign   data type stored by the input data structure
            \tparam      AForeign   allocator used by the input data structure
         */
        template<typename TForeign, class AForeign = std::allocator<TForeign> >
        explicit
        operator IrregArray4D<TForeign, AForeign>() const
        {
            IrregArray4D<TForeign, AForeign> result;
            result.initArray(*this);
            return result;
        }

    private:

        //! \brief allocate memory for the array using the prreset array geometry
        void allocateArray()
        {
            if (!isIrreg_)
            {
                if (first1_ > last1_)
                {
                    std::string errorMessage("Error in %s::allocateArray(): called with illegal bounds first1 > last1.", typeid(*this).name());
                    GMX_THROW(InvalidInputError(errorMessage));
                }
                if (first2_ > last2_)
                {
                    std::string errorMessage("Error in %s::allocateArray(): called with illegal bounds first2 > last2.", typeid(*this).name());
                    GMX_THROW(InvalidInputError(errorMessage));
                }
                if (first3_ > last3_)
                {
                    std::string errorMessage("Error in %s::allocateArray(): called with illegal bounds first3 > last3.", typeid(*this).name());
                    GMX_THROW(InvalidInputError(errorMessage));
                }
                if (first4_ > last4_)
                {
                    std::string errorMessage("Error in %s::allocateArray(): called with illegal bounds first4 > last4.", typeid(*this).name());
                    GMX_THROW(InvalidInputError(errorMessage));
                }
            }

            try
            {
                arr_  = std::allocator_traits<p3_allocator_type>::allocate(p3_allocator_, getLength1());
                arr_ -= getFirst1();
                for (index_type i = getFirst1(); i <= getLast1(); ++i)
                {
                    arr_[i]  = std::allocator_traits<p2_allocator_type>::allocate(p2_allocator_, getLength2(i));
                    arr_[i] -= getFirst2();
                }
                for (index_type i = getFirst1(); i <= getLast1(); i++)
                {
                    for (index_type j = getFirst2(); j <= getLast2(i); ++j)
                    {
                        arr_[i][j]  = std::allocator_traits<p_allocator_type>::allocate(p_allocator_, getLength3(i, j));
                        arr_[i][j] -= getFirst3();
                    }
                }
                for (index_type i = getFirst1(); i <= getLast1(); i++)
                {
                    for (index_type j = getFirst2(); j <= getLast2(i); ++j)
                    {
                        for (index_type k = getFirst3(); k <= getLast3(i, j); ++k)
                        {
                            arr_[i][j][k]  = std::allocator_traits<allocator_type>::allocate(allocator_, getLength4(i, j, k));
                            arr_[i][j][k] -= getFirst4();
                        }
                    }
                }
            }
            catch (const std::bad_alloc &)
            {
                std::fprintf(stderr, "Error in %s::allocateArray(): could not allocate memory for the data array.", typeid(*this).name());
                throw;
            }

            // initialize the array elements
            for (index_type i = getFirst1(); i <= getLast1(); ++i)
            {
                for (index_type j = getFirst2(); j <= getLast2(i); ++j)
                {
                    for (index_type k = getFirst3(); k <= getLast3(i, j); ++k)
                    {
                        for (index_type l = getFirst4(); l <= getLast4(i, j, k); ++l)
                        {
                            std::allocator_traits<allocator_type>::construct(allocator_, &this->operator()(i, j, k, l), value_type());
                        }
                    }
                }
            }
        }

        /*! \brief allocate memory for a regular special case of the array

            \note  currently unused because block allocation is not easily combinable with alignment requirements

            \param[in]    x1       starting index in dimension 1
            \param[in]    x2       last index in dimension 1
            \param[in]    y1       starting index in dimension 2
            \param[in]    y2       last index in dimension 2
            \param[in]    z1       starting index in dimension 3
            \param[in]    z2       last index in dimension 3
            \param[in]    k1       starting index in dimension 4
            \param[in]    k2       last index in dimension 4
         */
        void allocateArray(const index_type x1, const index_type x2, const index_type y1, const index_type y2, const index_type z1, const index_type z2, const index_type k1, const index_type k2)
        {
            isIrreg_ = false;
            if (first1_ > last1_)
            {
                std::string errorMessage("Error in %s::allocateArray(): called with illegal bounds first1 > last1.", typeid(*this).name());
                GMX_THROW(InvalidInputError(errorMessage));
            }
            if (first2_ > last2_)
            {
                std::string errorMessage("Error in %s::allocateArray(): called with illegal bounds first2 > last2.", typeid(*this).name());
                GMX_THROW(InvalidInputError(errorMessage));
            }
            if (first3_ > last3_)
            {
                std::string errorMessage("Error in %s::allocateArray(): called with illegal bounds first3 > last3.", typeid(*this).name());
                GMX_THROW(InvalidInputError(errorMessage));
            }
            if (first4_ > last4_)
            {
                std::string errorMessage("Error in %s::allocateArray(): called with illegal bounds first4 > last4.", typeid(*this).name());
                GMX_THROW(InvalidInputError(errorMessage));
            }
            const size_type d1 = x2 - x1 + 1;
            const size_type d2 = y2 - y1 + 1;
            const size_type d3 = z2 - z1 + 1;
            const size_type d4 = k2 - k1 + 1;
            try
            {
                arr_              = std::allocator_traits<p3_allocator_type>::allocate(p3_allocator_, d1);
                arr_              = arr_ - first1_;
                arr_[x1]          = std::allocator_traits<p2_allocator_type>::allocate(p2_allocator_, d1*d2);
                arr_[x1]         -= first2_;
                arr_[x1][y1]      = std::allocator_traits<p_allocator_type>::allocate(p_allocator_, d1*d2*d3);
                arr_[x1][y1]     -= first3_;
                arr_[x1][y1][z1]  = std::allocator_traits<allocator_type>::allocate(allocator_, d1*d2*d3*d4);
                arr_[x1][y1][z1] -= first4_;
                // set second dimension
                for (index_type i = y1 + 1; i <= y2; i++)
                {
                    arr_[x1][i] = arr_[x1][i-1] + d3;
                }
                // set third dimension
                for (index_type i = x1 + 1; i <= x2; i++)
                {
                    arr_[i]         = arr_[i-1] + d2;
                    arr_[i][y1]     = arr_[i-1][y1] + (d2 * d3);
                    arr_[i][y1][z1] = arr_[i-1][y1][z1] + (d2 * d3 * d4);
                    for (index_type j = y1 + 1; j <= y2; j++)
                    {
                        arr_[i][j] = arr_[i][j-1]+d3;
                    }
                }
                for (index_type i = x1; i <= x2; i++)
                {
                    for (index_type j = y1 + 1; j <= y2; j++)
                    {
                        arr_[i][j][z1] = arr_[i][j-1][z1] + (d3 * d4);
                    }
                }
                for (index_type i = x1; i <= x2; i++)
                {
                    for (index_type j = y1; j <= y2; j++)
                    {
                        for (index_type h = z1 + 1; h <= z2; h++)
                        {
                            arr_[i][j][h] = arr_[i][j][h-1] + d4;
                        }
                    }
                }
            }
            catch (const std::bad_alloc &)
            {
                std::fprintf(stderr, "Error in %s::allocateArray(): could not allocate memory for the data array.", typeid(*this).name());
                throw;
            }

            // initialize the array elements
            for (index_type i = getFirst1(); i <= getLast1(); ++i)
            {
                for (index_type j = getFirst2(); j <= getLast2(); ++j)
                {
                    for (index_type k = getFirst3(); k <= getLast3(); ++k)
                    {
                        for (index_type l = getFirst4(); l <= getLast4(); ++l)
                        {
                            std::allocator_traits<allocator_type>::construct(allocator_, &this->operator()(i, j, k, l), value_type());
                        }
                    }
                }
            }

        }
};

} // end namespace gmx

#endif
