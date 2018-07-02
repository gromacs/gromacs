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

    \brief FlatIrregArray4D irregularly shaped 4-dimensional array
           internally stored in an 1D-array

    \author R. Thomas Ullmann <tullman@gwdg.de>
 */
#ifndef GMX_UTILITY_DATASTRUCTURES_FLAT_IRREG_ARRAY_4D_H
#define GMX_UTILITY_DATASTRUCTURES_FLAT_IRREG_ARRAY_4D_H

#include "gromacs/utility/data_structures/flat_irreg_array_3d.h"

namespace gmx
{

/*! \class FlatIrregArray4D

    \brief irregularly shaped 4-dimensional array
           internally stored in an 1D-array

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
class FlatIrregArray4D
{
    public:
        /*  adopted these typedefs from AnBe's TriangularArray, for consistency with the other data
            structures, can be retrieved from outside to make sure that the intended types are used
            when accessing or manipulating the array */
        //! the data type stored in the array
        typedef T value_type;
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
        typedef IrregArray1D<size_type, allocator_size_type>     size_1d_array_type;
        //! type of array used for storing the lengths of array stripes
        typedef FlatIrregArray2D<size_type, allocator_size_type> size_2d_array_type;
        //! type of array used for storing the lengths of array stripes
        typedef FlatIrregArray3D<size_type, allocator_size_type> size_3d_array_type;
    private:
        //! first index of dimension 1
        index_type                  first1_;
        //! first index of dimension 2
        index_type                  first2_;
        //! first index of dimension 3
        index_type                  first3_;
        //! first index of dimension 4
        index_type                  first4_;
        //! number of elements of dimension 1
        index_type                  last1_;
        //! number of elements of dimension 2
        index_type                  last2_;
        //! number of elements of dimension 3
        index_type                  last3_;
        //! number of elements of dimension 4
        index_type                  last4_;
        //! total number of array elements arr_[:][:][:][:]
        size_type                   nelements_;
        //! irreg_[i][j][k] = x is the index of the first actual array element arr_[x] of the
        //! stripe FlatIrregArray3D[i][j][k][:] in dimension 4 of the of the represented 4D array
        size_3d_array_type          irreg_;
        //! flag says that dimX(i) = endX(i) - beginX(i) + 1 is not a constant
        bool                        isIrreg_;
        //! array storing the the actual data
        value_type                 *arr_;
        //! the allocator used to allocate arr_
        allocator_type              allocator_;
        //! let other FlatIrregArray4D instantations access private members (for the explicit conversion operator)
        template<class W, typename AW> friend class FlatIrregArray4D;
    public:
        /*! \brief   default constructor without memory allocation, allocateArray has to be called for
                     allocating memory initArray for memory allocation and content initialization
         */
        FlatIrregArray4D()
            : first1_(0), first2_(0), first3_(0), first4_(0), last1_(-1), last2_(-1), last3_(-1), last4_(-1), nelements_(0), isIrreg_(false), arr_(nullptr)
        {
        }

        /*! \brief  copy constructor

            \param[in]   ini   array to be copied
         */
        FlatIrregArray4D(const FlatIrregArray4D &ini)
            : first1_(0), first2_(0), first3_(0), first4_(0), last1_(-1), last2_(-1), last3_(-1), last4_(-1), nelements_(0), isIrreg_(false), arr_(nullptr)
        {
            initArray(ini);
        }

        /*! \brief  copy constructor

            \tparam      AForeign   allocator used by the input data structure

            \param[in]   ini   array to be copied
         */
        template<class AForeign = allocator_type>
        FlatIrregArray4D(const FlatIrregArray4D<value_type, AForeign> &ini)
            : first1_(0), first2_(0), first3_(0), first4_(0), last1_(-1), last2_(-1), last3_(-1), last4_(-1), nelements_(0), isIrreg_(false), arr_(nullptr)
        {
            initArray(ini);
        }

        //! move constructor
        //! \param[in]    ini   source to be moved
        FlatIrregArray4D(FlatIrregArray4D &&ini) noexcept
            : first1_(ini.first1_), first2_(ini.first2_), first3_(ini.first3_), first4_(ini.first4_),
              last1_(ini.last1_), last2_(ini.last2_), last3_(ini.last3_), last4_(ini.last4_),
              nelements_(ini.nelements_), irreg_(std::move(ini.irreg_)), isIrreg_(ini.isIrreg_),
              arr_(ini.arr_), allocator_(std::move(ini.allocator_))
        {
            // set ini to a valid default state
            ini.first1_    =  0;
            ini.first2_    =  0;
            ini.first3_    =  0;
            ini.first4_    =  0;
            ini.last1_     = -1;
            ini.last2_     = -1;
            ini.last3_     = -1;
            ini.last4_     = -1;
            ini.nelements_ =  0;
            ini.isIrreg_   = false;
            ini.arr_       = nullptr;
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
        FlatIrregArray4D(const index_type x1, const index_type x2, const index_type y1, const index_type y2, const index_type z1, const index_type z2, const index_type k1, const index_type k2)
            : first1_(0), first2_(0), first3_(0), first4_(0), last1_(-1), last2_(-1), last3_(-1), last4_(-1), nelements_(0), isIrreg_(false), arr_(nullptr)
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
        FlatIrregArray4D(const index_type x1, const index_type x2, const index_type y1, const index_type y2, const index_type z1, const index_type z2, const index_type k1, const index_type k2, const value_type &ini)
            : first1_(0), first2_(0), first3_(0), first4_(0), last1_(-1), last2_(-1), last3_(-1), last4_(-1), nelements_(0), isIrreg_(false), arr_(nullptr)
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
        FlatIrregArray4D(const index_type x1, const index_type x2, const IrregArray1D<size_type, AForeign> &sizes1, const FlatIrregArray2D<size_type, AForeign> &sizes2, const FlatIrregArray3D<size_type, AForeign> &sizes3)
            : first1_(0), first2_(0), first3_(0), first4_(0), last1_(-1), last2_(-1), last3_(-1), last4_(-1), nelements_(0), isIrreg_(true), arr_(nullptr)
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
        FlatIrregArray4D(const index_type x1, const index_type x2, const IrregArray1D<size_type, AForeign> &sizes1, const FlatIrregArray2D<size_type, AForeign> &sizes2, const FlatIrregArray3D<size_type, AForeign> &sizes3, const value_type &ini)
            : first1_(0), first2_(0), first3_(0), first4_(0), last1_(-1), last2_(-1), last3_(-1), last4_(-1), nelements_(0), isIrreg_(true), arr_(nullptr)
        {
            initArray(x1, x2, sizes1, sizes2, sizes3, ini);
        }

        /*! \brief   general constructor helper function for irregularly shaped matrices and all special cases

            \tparam      AForeign   allocator used by the input data structure

            \param[in]    sizes3   number of elements in dimension 4 as function of the indices in the lower dimensions
         */
        template <class AForeign = allocator_size_type>
        FlatIrregArray4D(const FlatIrregArray3D<size_type, AForeign> &sizes3)
            : first1_(0), first2_(0), first3_(0), first4_(0), last1_(-1), last2_(-1), last3_(-1), last4_(-1), nelements_(0), isIrreg_(true), arr_(nullptr)
        {
            initArray(sizes3);
        }

        /*! \brief   general constructor helper function for irregularly shaped matrices and all special cases

            \tparam      AForeign   allocator used by the input data structure

            \param[in]    sizes3   number of elements in dimension 4 as function of the indices in the lower dimensions
            \param[in]    ini      initial value for the array elements
         */
        template <class AForeign = allocator_size_type>
        FlatIrregArray4D(const FlatIrregArray3D<size_type, AForeign> &sizes3, const value_type &ini)
            : first1_(0), first2_(0), first3_(0), first4_(0), last1_(-1), last2_(-1), last3_(-1), last4_(-1), nelements_(0), isIrreg_(true), arr_(nullptr)
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
        FlatIrregArray4D(const index_type x1, const index_type x2, const IrregArray1D<size_type, AForeign> &sizes1, const index_type z1, const index_type z2, const IrregArray1D<size_type, AForeign> &sizes2, const value_type &ini)
            : first1_(x1), first2_(x1), first3_(z1), first4_(z1), last1_(x2), last2_(0), last3_(z2), last4_(0), nelements_(0), isIrreg_(true), arr_(nullptr)
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
        FlatIrregArray4D(const index_type x1, const index_type x2, const IrregArray1D<size_type, AForeign> &sizes1, const index_type z1, const index_type z2, const IrregArray1D<size_type, AForeign> &sizes2)
            : first1_(x1), first2_(x1), first3_(z1), first4_(z1), last1_(x2), last2_(0), last3_(z2), last4_(0), nelements_(0), isIrreg_(true), arr_(nullptr)
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
        FlatIrregArray4D(const std::initializer_list<std::initializer_list<std::initializer_list<std::initializer_list<value_type> > > > &ini)
            : first1_(0), first2_(0), first3_(0), first4_(0), last1_(-1), last2_(-1), last3_(-1), last4_(-1), nelements_(0), isIrreg_(true), arr_(nullptr)
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

                FlatIrregArray2D<size_type, allocator_size_type> irregArr2((index_type)0, (index_type)ini.size() - 1, irregArr1, (size_t)0);
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

                FlatIrregArray3D<size_type, allocator_size_type> irregArr3((index_type)0, (index_type)ini.size() - 1, irregArr1, irregArr2, (size_t)0);
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

        /*! \brief  initialize FlatIrregArray4D<value_type> from FlatIrregArray4D<TForeign>,
                    used by the copy constructor and by the explicit conversion operator,
                    TForeign must be convertible to value_type via static_cast

            \tparam      TForeign   data type stored by the input data structure
            \tparam      AForeign   allocator used by the input data structure

            \param[in]   ini   the template to be copied
         */
        template <typename TForeign, typename Allocator_foreign = std::allocator<TForeign> >
        void initArray(const FlatIrregArray4D<TForeign, Allocator_foreign> &ini)
        {
            irreg_     = ini.irreg_;
            isIrreg_   = ini.isIrreg_;
            first1_    = ini.first1_;
            last1_     = ini.last1_;
            first2_    = ini.first2_;
            last2_     = ini.last2_;
            first3_    = ini.first3_;
            last3_     = ini.last3_;
            first4_    = ini.first4_;
            last4_     = ini.last4_;
            nelements_ = ini.nelements_;
            allocateArray(nelements_);
            for (size_type i = 0; i < nelements_; ++i)
            {
                arr_[i] = static_cast<value_type>(ini.arr_[i]);
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
            allocateArray(first1_, last1_, first2_, last2_, first3_, last3_, first4_, last4_);
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
            for (size_type i = 0; i < nelements_; ++i)
            {
                arr_[i] = ini;
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
        template <class AForeign = allocator_size_type>
        void initArray(const index_type x1, const index_type y1, const index_type z1, const index_type k1, const IrregArray1D<size_type, AForeign> &sizes1, const FlatIrregArray2D<size_type, AForeign> &sizes2, const FlatIrregArray3D<size_type, AForeign> &sizes3)
        {
            deallocateArray();
            isIrreg_   = true;
            first1_    = x1;
            first2_    = y1;
            first3_    = z1;
            first4_    = k1;
            last1_     = x1 + sizes1.getLength1() - 1;
            last2_     = 0;
            last3_     = 0;
            last4_     = 0;
            allocateArray(first1_, first2_, first3_, first4_, sizes1, sizes2, sizes3);
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
        template <class AForeign = allocator_size_type>
        void initArray(const index_type x1, const index_type y1, const index_type z1, const index_type k1, const IrregArray1D<size_type, AForeign> &sizes1, const FlatIrregArray2D<size_type, AForeign> &sizes2, const FlatIrregArray3D<size_type, AForeign> &sizes3, const value_type &ini)
        {
            initArray(x1, y1, z1, k1, sizes1, sizes2, sizes3);
            for (size_type i = 0; i < nelements_; ++i)
            {
                arr_[i] = ini;
            }
        }
        /*! \brief general constructor helper function for irregularly shaped matrices and all special cases

            \tparam      AForeign   allocator used by the input data structure

            \param[in]    x1       starting index in dimension 1
            \param[in]    x2       last index in dimension 1
            \param[in]    sizes1   number of elements in dimension 2 as function of the index in dimension 1
            \param[in]    sizes2   number of elements in dimension 3 as function of the indices in the lower dimensions
            \param[in]    sizes3   number of elements in dimension 4 as function of the indices in the lower dimensions
         */
        template <class AForeign = allocator_size_type>
        void initArray(const index_type x1, const index_type x2, const IrregArray1D<size_type, AForeign> &sizes1, const FlatIrregArray2D<size_type, AForeign> &sizes2, const FlatIrregArray3D<size_type, AForeign> &sizes3)
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
            allocateArray(first1_, first2_, first3_, first4_, sizes1, sizes2, sizes3);
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
        template <class AForeign = allocator_size_type>
        void initArray(const index_type x1, const index_type x2, const IrregArray1D<size_type, AForeign> &sizes1, const FlatIrregArray2D<size_type, AForeign> &sizes2, const FlatIrregArray3D<size_type, AForeign> &sizes3, const value_type &ini)
        {
            initArray(x1, x2, sizes1, sizes2, sizes3);
            for (size_type i = 0; i < nelements_; ++i)
            {
                arr_[i] = ini;
            }
        }

        /*! \brief general constructor helper function for irregularly shaped matrices and all special cases

            \tparam      AForeign   allocator used by the input data structure

            \param[in]    sizes3   number of elements in dimension 4 as function of the indices in the lower dimensions
            \param[in]    ini      initialization value for the data array elements
         */
        template <class AForeign = allocator_size_type>
        void initArray(const FlatIrregArray3D<size_type, AForeign> &sizes3, const value_type &ini)
        {
            initArray(sizes3);
            for (size_type i = 0; i < nelements_; ++i)
            {
                arr_[i] = ini;
            }
        }
        /*! \brief general constructor helper function for irregularly shaped matrices and all special cases

            \tparam      AForeign   allocator used by the input data structure

            \param[in]    sizes3   number of elements in dimension 4 as function of the indices in the lower dimensions
         */
        template <class AForeign = allocator_size_type>
        void initArray(const FlatIrregArray3D<size_type, AForeign> &sizes3)
        {
            IrregArray1D<size_type, AForeign> sizes1(sizes3.getFirst1(), sizes3.getLast1(), (size_type)0);
            for (index_type i = sizes3.getFirst1(); i <= sizes3.getLast1(); ++i)
            {
                sizes1(i) = sizes3.getLength2(i);
            }
            FlatIrregArray2D<size_type, AForeign> sizes2(sizes3.getFirst1(), sizes3.getLast1(), sizes1, (size_type)0);
            for (index_type i = sizes3.getFirst1(); i <= sizes3.getLast1(); ++i)
            {
                for (index_type j = sizes3.getFirst2(); j <= sizes3.getLast2(i); ++j)
                {
                    sizes2(i, j) = sizes3.getLength3(i, j);
                }
            }
            initArray(sizes3.getFirst1(), sizes3.getLast1(), sizes1, sizes2, sizes3);
        }

        /*! \brief specialized constructor helper function for irregular, symmetrically shaped matrices
                   application: interaction energy matrix site1 form1(site1) site2 form2(site2)

            \tparam      AForeign   allocator used by the input data structure

            \param[in]    x1       starting index in dimension 1
            \param[in]    x2       last index in dimension 1
            \param[in]    sizes1   number of index values in dimension 2
                                   as function of the indices in dimension 1
                                   (in this special case, independent on the
                                   index in dimensions 1 and 2)
            \param[in]    z1       starting index in dimension 3
            \param[in]    z2       last index in dimension 3
            \param[in]    sizes2   number of index values in dimension 4
                                   as function of the indices in dimension 3
                                   (in this special case, independent on the
                                   index in dimensions 1 and 2)
         */
        template <class AForeign = allocator_size_type>
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
            const index_type         dim3 = z2 - z1 + 1;
            FlatIrregArray2D<size_t, allocator_size_type> tmp_irreg_mat(first1_, last1_, sizes1, dim3);
            for (index_type i = first1_; i <= last1_; i++)
            {
                for (index_type j = first2_; j <= first2_ + static_cast<index_type>(sizes1(i)) - 1; j++)
                {
                    tmp_irreg_mat(i, j) = dim3;
                }
            }
            FlatIrregArray3D<size_t, allocator_size_type> tmp_irreg_ten(first1_, last1_, sizes1, tmp_irreg_mat, (size_type)0);
            for (index_type i = first1_; i <= last1_; i++)
            {
                for (index_type j = first2_; j <= first2_ + static_cast<index_type>(sizes1(i)) - 1; j++)
                {
                    for (index_type k = first3_; k <= last3_; k++)
                    {
                        tmp_irreg_ten(i, j, k) = sizes2(k);
                    }
                }
            }
            allocateArray(first1_, first2_, first3_, first4_, sizes1, tmp_irreg_mat, tmp_irreg_ten);
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
        template <class AForeign = allocator_size_type>
        void initArray(const index_type x1, const index_type x2, const IrregArray1D<size_type, AForeign> &sizes1, const index_type z1, const index_type z2, const IrregArray1D<size_type, AForeign> &sizes2, const value_type &ini)
        {
            initArray(x1, x2, sizes1, z1, z2, sizes2);
            for (size_type i = 0; i < nelements_; ++i)
            {
                arr_[i] = ini;
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
                std::fprintf(stderr, "Error in %s::getLength2(): called without argument for an irregular array.", typeid(*this).name());
                throw InternalError("Extent of a variable dimension in an irregular array requested without specifying the index/indices in the lower dimension(s).");
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
                std::fprintf(stderr, "Error in %s::getLength3(): called without arguments for an irregular array.", typeid(*this).name());
                throw InternalError("Extent of a variable dimension in an irregular array requested without specifying the index/indices in the lower dimension(s).");
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
                std::fprintf(stderr, "Error in %s::getLength4(): called without arguments for an irregular array.", typeid(*this).name());
                throw InternalError("Extent of a variable dimension in an irregular array requested without specifying the index/indices in the lower dimension(s).");
            }
        }
        /*! \brief number of elements in dimension 4 at lower dimension indices x, y
            \param[in]    x       index in dimension 1
            \param[in]    y       index in dimension 2
            \param[in]    z       index in dimension 3  */
        size_type getLength4(const index_type x, const index_type y, const index_type z) const
        {
            if (isIrreg_)
            {
                if (x == last1_ && y == getLast2(last1_) && z == getLast3(last1_, getLast2(last1_)))
                {
                    return nelements_ - irreg_(x, y, z);
                }
                else
                {
                    // get first element of next stripe in the underlying data array of irreg_mat to avoid further conditionals
                    const size_type        i       = irreg_.getArrayIndex(x, y, z);
                    const size_type* const a       = irreg_.getArray();
                    return a[i + 1] - a[i];
                }
            }
            else
            {
                return last4_ - first4_ + 1;
            }
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
                std::fprintf(stderr, "Error in %s::getLast2(): called without argument for an irregular array.", typeid(*this).name());
                throw InternalError("Extent of a variable dimension in an irregular array requested without specifying the index/indices in the lower dimension(s).");
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
                std::fprintf(stderr, "Error in %s::getLast3(): called without arguments for an irregular array.", typeid(*this).name());
                throw InternalError("Extent of a variable dimension in an irregular array requested without specifying the index/indices in the lower dimension(s).");
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
                std::fprintf(stderr, "Error in %s::getLast4(): called without arguments for an irregular array.", typeid(*this).name());
                throw InternalError("Extent of a variable dimension in an irregular array requested without specifying the index/indices in the lower dimension(s).");
            }
        }
        /*! \brief get the last index of dimension 3 for lower dimension indices x, y, z

            \param[in]    x       index in dimension 1
            \param[in]    y       index in dimension 2
            \param[in]    z       index in dimension 3                       */
        index_type getLast4(const index_type x, const index_type y, const index_type z) const
        {
            if (isIrreg_)
            {
                if (x == last1_ && y == getLast2(last1_) && z == getLast3(last1_, getLast2(last1_)))
                {
                    return nelements_ - 1 - irreg_(x, y, z) + first4_;
                }
                else
                {
                    // get first element of next stripe in the underlying data array of irreg_mat to avoid further conditionals
                    const size_type        i       = irreg_.getArrayIndex(x, y, z);
                    const size_type* const a       = irreg_.getArray();
                    return a[i + 1] - 1 - a[i] + first4_;
                }
            }
            else
            {
                return last4_;
            }
        }

        //! get the index one past the last valid index of dimension 1
        index_type getEnd1() const { return last1_ + 1; }

        //! get the index one past the last valid index of dimension 2
        index_type getEnd2() const
        {
            if (!isIrreg_)
            {
                return last2_ + 1;
            }
            else
            {
                std::fprintf(stderr, "Error in %s::getEnd2(): called without argument for an irregular array.", typeid(*this).name());
                throw InternalError("Extent of a variable dimension in an irregular array requested without specifying the index/indices in the lower dimension(s).");
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
                std::fprintf(stderr, "Error in %s::getEnd3(): called without arguments for an irregular array.", typeid(*this).name());
                throw InternalError("Extent of a variable dimension in an irregular array requested without specifying the index/indices in the lower dimension(s).");
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
                std::fprintf(stderr, "Error in %s::getEnd4(): called without arguments for an irregular array.", typeid(*this).name());
                throw InternalError("Extent of a variable dimension in an irregular array requested without specifying the index/indices in the lower dimension(s).");
            }
        }
        /*! \brief get the index one past the last valid index of dimension 4 for lower dimension indices x, y, z

            \param[in]    x       index in dimension 1
            \param[in]    y       index in dimension 2
            \param[in]    z       index in dimension 3                       */
        index_type getEnd4(const index_type x, const index_type y, const index_type z) const
        {
            if (isIrreg_)
            {
                if (x == last1_ && y == getLast2(last1_) && z == getLast3(last1_, getLast2(last1_)))
                {
                    return nelements_ - irreg_(x, y, z) + first4_;
                }
                else
                {
                    // get first element of next stripe in the underlying data array of irreg_mat to avoid further conditionals
                    const size_type        i       = irreg_.getArrayIndex(x, y, z);
                    const size_type* const a       = irreg_.getArray();
                    return a[i + 1] - a[i] + first4_;
                }
            }
            else
            {
                return last4_ + 1;
            }
        }

        /*! \brief determine the index of the array element IrregArray2D(i, j) in the underlying 1D array arr
            this function needs to be public to be compatible to C++ standards < C++11, where nested classes
            can not access private members of the enclosing class. A friend declaration can also not be used
            because the friends are defined to be entities that are not members.

            length_type instead of index_type used intentionally here to use the maximum possible length

            \param[in]   x    index in dimension 1
            \param[in]   y    index in dimension 2
            \param[in]   z    index in dimension 3
            \param[in]   k    index in dimension 4
         */
        size_type getArrayIndex(const index_type x, const index_type y, const index_type z, const index_type k) const
        {
            if (isIrreg_)
            {
                return irreg_(x, y, z) + (k - first4_);
            }
            else
            {
                return (getLength2() * getLength3() * getLength4() * (x - first1_))
                       + (getLength3() * getLength4() * (y - first2_))
                       + (getLength4() * (z - first3_))
                       + (k - first4_);
            }
        }

        //! returns the amount of memory occupied by the object
        size_type getSize() const
        {
            // the book-keeping data
            size_type size = sizeof(*this);
            size += irreg_.getSize() - sizeof(FlatIrregArray3D<size_type>);
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
        size_type getNelements() const { return nelements_; }

        //! deallocate the memory of all constituent arrays
        void deallocateArray()
        {
            if (arr_ != nullptr)
            {

                // destroy the array elements
                for (size_type i = 0; i < nelements_; ++i)
                {
                    std::allocator_traits<allocator_type>::destroy(allocator_, &arr_[i]);
                }

                std::allocator_traits<allocator_type>::deallocate(allocator_, arr_, nelements_);
                arr_       = nullptr;
                nelements_ = 0;
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
        ~FlatIrregArray4D()
        {
            // deallocate arrays
            deallocateArray();
        }

        // index operators

        /*! index operator for bare array-like usage a[i][j][k][l] */

        //! doxygen can not handle nested classes
        /// \cond DEV
        /*! \class Proxy

            \brief helper functor used to implement index operator *this[i][j][k][l]
         */
        class Proxy
        {
            private:
                FlatIrregArray4D* optr_;
                index_type        x_;
                class Proxy2
                {
                    private:
                        FlatIrregArray4D* optr_;
                        index_type        x_;
                        index_type        y_;
                        class Proxy3
                        {
                            private:
                                FlatIrregArray4D* optr_;
                                index_type        x_;
                                index_type        y_;
                                index_type        z_;
                            public:
                                Proxy3(FlatIrregArray4D* optr, index_type x, index_type y, index_type z)
                                    : optr_(optr), x_(x), y_(y), z_(z) { }

                                value_type &operator[](index_type k)
                                {
                                    return optr_->operator()(x_, y_, z_, k);
                                }
                                const value_type &operator[](index_type k) const
                                {
                                    return optr_->operator()(x_, y_, z_, k);
                                }
                        };
                    public:
                        Proxy2(FlatIrregArray4D* optr, index_type x, index_type y)
                            : optr_(optr), x_(x), y_(y) { }

                        Proxy3 operator[](index_type z)
                        {
                            return Proxy3(optr_, x_, y_, z);
                        }
                        const Proxy3 operator[](index_type z) const
                        {
                            return Proxy3(optr_, x_, y_, z);
                        }
                };
            public:
                Proxy(FlatIrregArray4D* optr, index_type x)
                    : optr_(optr), x_(x) { }

                Proxy2 operator[](index_type y)
                {
                    return Proxy2(optr_, x_, y);
                }
                const Proxy2 operator[](index_type y) const
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
        Proxy operator[](index_type x)
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
        const Proxy operator[](index_type x) const
        {
            return Proxy(this, x);
        }

        /*! \brief   index operator(x,y,z) for functor-like usage
            \param[in]   x   index in dimension 1
            \param[in]   y   index in dimension 2
            \param[in]   z   index in dimension 3
            \param[in]   k   index in dimension 4
         */
        inline value_type &operator()(const index_type x, const index_type y, const index_type z, const index_type k)
        {
            return arr_[getArrayIndex(x, y, z, k)];
        }
        /*! \brief   const index operator(x,y,z) for functor-like usage
            \param[in]   x   index in dimension 1
            \param[in]   y   index in dimension 2
            \param[in]   z   index in dimension 3
            \param[in]   k   index in dimension 4
         */
        inline const value_type &operator()(const index_type x, const index_type y, const index_type z, const index_type k) const
        {
            return arr_[getArrayIndex(x, y, z, k)];
        }

        /*! \brief   index operator(x,y,z) for functor-like usage
            \param[in]   x   index in dimension 1
            \param[in]   y   index in dimension 2
            \param[in]   z   index in dimension 3
         */
        inline value_type* operator()(const index_type x, const index_type y, const index_type z)
        {
            return &arr_[getArrayIndex(x, y, z, getBegin4())];
        }
        /*! \brief   const index operator(x,y,z) for functor-like usage
            \param[in]   x   index in dimension 1
            \param[in]   y   index in dimension 2
            \param[in]   z   index in dimension 3
         */
        inline const value_type* operator()(const index_type x, const index_type y, const index_type z) const
        {
            return &arr_[getArrayIndex(x, y, z, getBegin4())];
        }

        //! get a pointer to the data array
        const value_type* getArray() const
        {
            return const_cast<const value_type*>(arr_);
        }
        //! get a pointer to the data array
        value_type* getArray()
        {
            return arr_;
        }

        //! get a pointer to the data array
        const value_type* data() const
        {
            return const_cast<const value_type*>(arr_);
        }
        //! get a pointer to the data array
        value_type* data()
        {
            return arr_;
        }

        /*! \brief   stream operator for convenient printing of the array to an ostream

            \param[in]   output   output stream for the array contents
            \param[in]   ten      the array to be printed
         */
        friend std::ostream &operator<<(std::ostream &output, const FlatIrregArray4D &ten)
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
        FlatIrregArray4D &operator=(FlatIrregArray4D &&rhs) noexcept
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
                nelements_    = rhs.nelements_;
                isIrreg_      = rhs.isIrreg_;
                arr_          = rhs.arr_;
                irreg_        = std::move(rhs.irreg_);
                allocator_    = std::move(rhs.allocator_);
                // set rhs to a valid default state
                rhs.first1_    =  0;
                rhs.first2_    =  0;
                rhs.first3_    =  0;
                rhs.first4_    =  0;
                rhs.last1_     = -1;
                rhs.last2_     = -1;
                rhs.last3_     = -1;
                rhs.last4_     = -1;
                rhs.nelements_ =  0;
                rhs.isIrreg_   = false;
                rhs.arr_       = nullptr;
            }
            return *this;
        }

        /*! \brief   assignment operator *this = rhs

            \tparam      AForeign   allocator used by the input data structure

            \param[in]   rhs   right-hand side object of the assignment
         */
        FlatIrregArray4D &operator=(const FlatIrregArray4D &rhs)
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
        FlatIrregArray4D &operator=(const FlatIrregArray4D<value_type, AForeign> &rhs)
        {
            deallocateArray();
            initArray(rhs);
            return *this;
        }

        /*! \brief  conversion operator (explicit to avoid unintended implicit conversion), usable via
                    FlatIrregArray4D<TForeign, AForeign> newArr = static_cast<FlatIrregArray4D<TForeign, AForeign> >(FlatIrregArray4D<TCurrent, ACurrent> oldArr);

            \tparam      TForeign   data type stored by the input data structure
            \tparam      AForeign   allocator used by the input data structure
         */
        template<typename TForeign, class AForeign = std::allocator<TForeign> >
        explicit
        operator FlatIrregArray4D<TForeign, AForeign>() const
        {
            FlatIrregArray4D<TForeign, AForeign> result;
            result.initArray(*this);
            return result;
        }

    private:

        /*! \brief allocate memory for a regular special case of the array

            \param[in]    x1       starting index in dimension 1
            \param[in]    x2       last index in dimension 1
            \param[in]    y1       starting index in dimension 2
            \param[in]    y2       last index in dimension 2
            \param[in]    z1       starting index in dimension 3
            \param[in]    z2       last index in dimension 3
            \param[in]    k1       starting index in dimension 4
            \param[in]    k2       last index in dimension 4              */
        void allocateArray(const index_type x1, const index_type x2, const index_type y1, const index_type y2, const index_type z1, const index_type z2, const index_type k1, const index_type k2)
        {
            deallocateArray(); // prevent memory leaks
            isIrreg_   = false;
            first1_    = x1;
            last1_     = x2;
            first2_    = y1;
            last2_     = y2;
            first3_    = z1;
            last3_     = z2;
            first4_    = k1;
            last4_     = k2;
            if (first1_ > last1_)
            {
                std::fprintf(stderr, "Error in %s::allocateArray(): called with illegal bounds begin1 > end1.", typeid(*this).name());
                throw InternalError("Trying to create an irregular array with illegal bounds begin1 > end1");
            }
            if (first2_ > last2_)
            {
                std::fprintf(stderr, "Error in %s::allocateArray(): called with illegal bounds begin2 > end2.", typeid(*this).name());
                throw InternalError("Trying to create an irregular array with illegal bounds begin2 > end2");
            }
            if (first3_ > last3_)
            {
                std::fprintf(stderr, "Error in %s::allocateArray(): called with illegal bounds begin3 > end3.", typeid(*this).name());
                throw InternalError("Trying to create an irregular array with illegal bounds begin3 > end3");
            }
            if (first4_ > last4_)
            {
                std::fprintf(stderr, "Error in %s::allocateArray(): called with illegal bounds begin4 > end4.", typeid(*this).name());
                throw InternalError("Trying to create an irregular array with illegal bounds begin4 > end4");
            }
            const size_type dim1 = getLength1();
            const size_type dim2 = getLength2();
            const size_type dim3 = getLength3();
            const size_type dim4 = getLength4();
            nelements_ = dim1 * dim2 * dim3 * dim4;
            allocateArray(nelements_);
        }

        /*! \brief allocate memory for the array

            \tparam      AForeign   allocator used by the input data structure

            \param[in]    s1       starting index in dimension 1
            \param[in]    s2       starting index in dimension 2
            \param[in]    s3       starting index in dimension 3
            \param[in]    s4       starting index in dimension 4
            \param[in]    d1       number of elements in dimension 1
            \param[in]    sizes1   number of elements in dimension 2 as function of the index in dimension 1
            \param[in]    sizes2   number of elements in dimension 3 as function of the indices in the lower dimensions
            \param[in]    sizes3   number of elements in dimension 4 as function of the indices in the lower dimensions
         */
        template<class AForeign = allocator_size_type>
        void allocateArray(const index_type s1, const index_type s2, const index_type s3, const index_type s4, const index_type d1, size_type *sizes1, const FlatIrregArray2D<size_type, AForeign> &sizes2, const FlatIrregArray3D<size_type, AForeign> &sizes3)
        {
            allocateArray(s1, s2, s3, s4, IrregArray1D<size_type, allocator_size_type>(s1, (s1 + d1 - 1), sizes1), sizes2, sizes3);
        }
        /*! \brief allocate memory for the array

            \tparam      AForeign   allocator used by the input data structure

            \param[in]    s1       starting index in dimension 1
            \param[in]    s2       starting index in dimension 2
            \param[in]    s3       starting index in dimension 3
            \param[in]    s4       starting index in dimension 4
            \param[in]    sizes1   number of elements in dimension 2 as function of the index in dimension 1
            \param[in]    sizes2   number of elements in dimension 3 as function of the indices in the lower dimensions
            \param[in]    sizes3   number of elements in dimension 4 as function of the indices in the lower dimensions
         */
        template<class AForeign = allocator_size_type>
        void allocateArray(const index_type s1, const index_type s2, const index_type s3, const index_type s4, const IrregArray1D<size_type, AForeign> &sizes1, const FlatIrregArray2D<size_type, AForeign> &sizes2, const FlatIrregArray3D<size_type, AForeign> &sizes3)
        {
            deallocateArray(); // prevent memory leaks
            isIrreg_   = true;
            first1_    = s1;
            first2_    = s2;
            first3_    = s3;
            first4_    = s4;
            last1_     = sizes1.getLast();
            last2_     = 0;
            last3_     = 0;
            last4_     = 0;

            const size_type dim1 = getLength1();
            // zero length seems to be OK with the standard allocator, but size_type may be signed
            if (static_cast<long int>(dim1) < 0)
            {
                std::fprintf(stderr, "Error in %s::allocateArray(): called with illegal bounds begin1 >= end1.", typeid(*this).name());
                throw InternalError("Trying to create an irregular array with illegal bounds begin1 >= end1");
            }

            IrregArray1D<size_type, allocator_size_type> irregVec(first1_, last1_);
            for (index_type i = first1_; i <= last1_; i++)
            {
                if (static_cast<long int>(sizes1(i)) <= 0)
                {
                    std::fprintf(stderr, "Error in %s::allocateArray(): called with illegal length %li <= 0 for element %li\n",
                                 typeid(*this).name(), static_cast<long int>(sizes1(i)), static_cast<long int>(i));
                    throw InternalError("Trying to create an irregular array with illegal bounds begin2 >= end2");
                }
                irregVec(i) = sizes1(i);
            }

            FlatIrregArray2D<size_type, allocator_size_type> irregMat(first1_, last1_, irregVec);
            for (index_type i = first1_; i <= last1_; i++)
            {
                for (index_type j = first2_; j <= first2_ + static_cast<index_type>(irregVec(i)) - 1; j++)
                {
                    if (static_cast<long int>(sizes2(i, j)) < 1)
                    {
                        std::fprintf(stderr, "Error in %s::allocateArray(): called with illegal length %ld < 1 for element %li, %li\n",
                                     typeid(*this).name(), static_cast<long int>(sizes2(i, j)),
                                     static_cast<long int>(i), static_cast<long int>(j));
                        throw InternalError("Trying to create an irregular array with illegal bounds begin3 >= end3");
                    }
                    irregMat(i, j) = sizes2(i, j);
                }
            }

            irreg_.initArray(first1_, last1_, irregVec, irregMat);
            size_type  last_length = 0;
            index_type last_index  = 0;
            for (index_type i = first1_; i <= last1_; i++)
            {
                for (index_type j = first2_; j <= first2_ + static_cast<index_type>(irregVec(i)) - 1; j++)
                {
                    for (index_type k = first3_; k <= first3_ + static_cast<index_type>(irregMat(i, j)) - 1; k++)
                    {
                        // zero length seems to be OK with the standard allocator, but size_type may be signed
                        if (static_cast<long int>(sizes3(i, j, k)) < 0)
                        {
                            std::fprintf(stderr, "Error in %s::allocateArray(): called with illegal length %ld < 0 for element %ld, %ld, %ld\n",
                                         typeid(*this).name(), static_cast<long int>(sizes3(i, j, k)),
                                         static_cast<long int>(i), static_cast<long int>(j), static_cast<long int>(k));
                            throw InternalError("Trying to create an irregular array with illegal bounds begin4 > end4");
                        }
                        if (i == first1_ && j == first2_ && k == first3_)
                        {
                            irreg_(i, j, k)  = 0;
                            last_index       = irreg_(i, j, k);
                            last_length      = sizes3(i, j, k);
                        }
                        else
                        {
                            irreg_(i, j, k) = last_index + last_length;
                            last_index      = irreg_(i, j, k);
                            last_length     = sizes3(i, j, k);
                        }
                    }
                }
            }
            // the number of elements is given by the starting index of the last array stripe
            // FlatIrregArray2D[i][:] in dimension 2 within arr + the length of this stripe + 1
            const size_type n = last_index + sizes3(last1_, getLast2(last1_), getLast3(last1_, getLast2(last1_)));
            allocateArray(n);
        }
        /*! \brief  allocate memory for the array

            \param[in]   n  total number of array elements
         */
        void allocateArray(const size_type n)
        {
            nelements_ = n;
            try
            {
                arr_ = std::allocator_traits<allocator_type>::allocate(allocator_, nelements_);
            }
            catch (const std::bad_alloc &)
            {
                std::fprintf(stderr, "Error in %s::allocateArray(): could not allocate memory for the data array.", typeid(*this).name());
                throw;
            }

            // initialize the array elements
            for (size_type i = 0; i < nelements_; ++i)
            {
                std::allocator_traits<allocator_type>::construct(allocator_, &arr_[i], value_type());
            }
        }
};

} // end namespace gmx

#endif
