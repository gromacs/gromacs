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

    \brief IrregArray3D is an irregularly shaped 3-dimensional array

    \author R. Thomas Ullmann <tullman@gwdg.de>
 */
#ifndef GMX_UTILITY_DATASTRUCTURES_IRREG_ARRAY_3D_H
#define GMX_UTILITY_DATASTRUCTURES_IRREG_ARRAY_3D_H

#include "gromacs/utility/data_structures/irreg_array_2d.h"

namespace gmx
{


/*! \class IrregArray3D

    \brief irregularly shaped 3 dimensional array

    this array has nonuniform second and third dimensions,
    that is, instead of a regular array with constant dimensions
    M x N x O, an irregular array with dimensions M x N(M) x O(M, N).

    example application:
    number of ligands of each ligand typebound by a site form
    indexed by site, form, ligand

    \author R. Thomas Ullmann <tullman@gwdg.de>

    \copyright GROMACS license

    \date Feb 2015

    \ingroup module_utility
    \inpublicapi

    \tparam   T           data type to be stored
    \tparam   Allocator   allocator to be used in creating arrays storing T
                          allocators for storing additional book-keeping data
                          are derived from this allocator via rebind if necessary
 */
template<class T = size_t, class Allocator = std::allocator<T> >
class IrregArray3D
{
    public:
        /*  adopted these typedefs from AnBe's TriangularArray, for consistency with the other data
            structures, can be retrieved from outside to make sure that the intended types are used
            when accessing or manipulating the array */
        //! the data type stored in the array
        typedef T value_type;
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
    private:
        //! first index of dimension 1
        index_type              first1_;
        //! first index of dimension 2
        index_type              first2_;
        //! first index of dimension 3
        index_type              first3_;
        //! ending index of dimension 1
        index_type              last1_;
        //! ending index of dimension 2
        index_type              last2_;
        //! ending index of dimension 3
        index_type              last3_;
        //! row,column-specific number of elements of dimension 3
        size_2d_array_type      irreg_;
        //! flag says that dimX(i) = last_X(i) - first_X(i) + 1 is not a constant
        bool                    isIrreg_;
        //! array storing the the actual data
        value_type           ***arr_;
        //! the allocator_ used to allocate arr_
        p2_allocator_type       p2_allocator_;
        //! the allocator_ used to allocate arr_[i]
        p_allocator_type        p_allocator_;
        //! the allocator_ used to allocate arr_[i][j]
        allocator_type          allocator_;
        //! let other IrregArray3D instantations access private members (for the explicit conversion operator)
        template<class V, typename AV> friend class IrregArray3D;
        //! let higher-dimensional IrregArrayNDs access private members of their book-keeping arrays
        template<class W, typename AW> friend class IrregArray4D;
    public:
        /*! \brief default constructor without memory allocation, allocateArray has to be called for
                   allocating memory initArray for memory allocation and content initialization
         */
        IrregArray3D() : first1_(0), first2_(0), first3_(0), last1_(-1), last2_(-1), last3_(-1), isIrreg_(false), arr_(nullptr)
        {
        }

        //! copy constructor
        //! \param[in]   ini
        IrregArray3D(const IrregArray3D &ini) : first1_(0), first2_(0), first3_(0), isIrreg_(false), arr_(nullptr)
        {
            initArray(ini);
        }

        /*! \brief templated copy constructor allowing for a different allocator of the input array

            \tparam      AForeign   allocator used by the input data structure
         */
        template<class AForeign = allocator_type>
        IrregArray3D(const IrregArray3D<value_type, AForeign> &ini) : first1_(0), first2_(0), first3_(0), isIrreg_(false), arr_(nullptr)
        {
            initArray(ini);
        }

        //! move constructor
        //! \param[in]    ini   source to be moved
        IrregArray3D(IrregArray3D &&ini) noexcept
            : first1_(ini.first1_), first2_(ini.first2_), first3_(ini.first3_), last1_(ini.last1_), last2_(ini.last2_), last3_(ini.last3_),
              irreg_(std::move(ini.irreg_)), isIrreg_(ini.isIrreg_), arr_(ini.arr_),
              p2_allocator_(std::move(ini.p2_allocator_)), p_allocator_(std::move(ini.p_allocator_)), allocator_(std::move(ini.allocator_))
        {
            // set ini to a valid default state
            ini.first1_   =  0;
            ini.first2_   =  0;
            ini.first3_   =  0;
            ini.last1_    = -1;
            ini.last2_    = -1;
            ini.last3_    = -1;
            ini.isIrreg_  = false;
            ini.arr_      = nullptr;
        }

        /*! \brief general constructor for irregularly shaped matrices and all special cases

            \tparam      AForeign   allocator used by the input data structure

            \param[in]    sizes2   number of elements in dimension 3 as function of the indices in the lower dimensions
         */
        template<class AForeign = allocator_size_type>
        IrregArray3D(const IrregArray2D<size_type, AForeign> &sizes2)
            : first1_(0), first2_(0), first3_(0), last1_(-1), last2_(-1), last3_(-1), isIrreg_(false), arr_(nullptr)
        {
            initArray(sizes2);
        }

        /*! \brief general constructor for irregularly shaped matrices and all special cases

            \tparam      AForeign   allocator used by the input data structure

            \param[in]    sizes2   number of elements in dimension 3 as function of the indices in the lower dimensions
            \param[in]    ini      initialization value for the data array elements
         */
        template<class AForeign = allocator_size_type>
        IrregArray3D(const IrregArray2D<size_type, AForeign> &sizes2, const value_type &ini)
            : first1_(0), first2_(0), first3_(0), last1_(-1), last2_(-1), last3_(-1), isIrreg_(false), arr_(nullptr)
        {
            initArray(sizes2, ini);
        }

        /*! \brief   constructor for a regularly shaped array

            \param[in]    x1       starting index in dimension 1
            \param[in]    x2       last index in dimension 1
            \param[in]    y1       starting index in dimension 2
            \param[in]    y2       last index in dimension 2
            \param[in]    z1       starting index in dimension 3
            \param[in]    z2       last index in dimension 3
         */
        IrregArray3D(const index_type x1, const index_type x2, const index_type y1, const index_type y2, const index_type z1, const index_type z2)
            : first1_(0), first2_(0), first3_(0), last1_(-1), last2_(-1), last3_(-1), isIrreg_(false), arr_(nullptr)
        {
            initArray(x1, x2, y1, y2, z1, z2);
        }
        /*! \brief   constructor for a regularly shaped data data with initialization

            \param[in]    x1       starting index in dimension 1
            \param[in]    x2       last index in dimension 1
            \param[in]    y1       starting index in dimension 2
            \param[in]    y2       last index in dimension 2
            \param[in]    z1       starting index in dimension 3
            \param[in]    z2       last index in dimension 3
            \param[in]    ini      initialization value for the data array elements
         */
        IrregArray3D(const index_type x1, const index_type x2, const index_type y1, const index_type y2, const index_type z1, const index_type z2, const value_type &ini)
            : first1_(0), first2_(0), first3_(0), last1_(-1), last2_(-1), last3_(-1), isIrreg_(false), arr_(nullptr)

        {
            initArray(x1, x2, y1, y2, z1, z2, ini);
        }

        /*! \brief general constructor for irregularly shaped matrices and all special cases

            \tparam      AForeign   allocator used by the input data structure

            \param[in]    x1       starting index in dimension 1
            \param[in]    x2       last index in dimension 1
            \param[in]    sizes    number of elements in dimension 2 as function of the index in dimension 1
            \param[in]    z1       starting index in dimension 2
            \param[in]    z2       last index in dimension 2
         */
        template<class AForeign = allocator_size_type>
        IrregArray3D(const index_type x1, const index_type x2, const IrregArray1D<size_type, AForeign> &sizes, const index_type z1, const index_type z2)
            : first1_(0), first2_(0), first3_(0), last1_(-1), last2_(-1), last3_(-1), isIrreg_(false), arr_(nullptr)
        {
            initArray(x1, x2, sizes, z1, z2);
        }
        /*! \brief general constructor for irregularly shaped matrices and all special cases

            \tparam      AForeign   allocator used by the input data structure

            \param[in]    x1       starting index in dimension 1
            \param[in]    x2       last index in dimension 1
            \param[in]    sizes    number of elements in dimension 2 as function of the index in dimension 1
            \param[in]    z1       starting index in dimension 2
            \param[in]    z2       last index in dimension 2
            \param[in]    ini      initialization value for the data array elements
         */
        template<class AForeign = allocator_size_type>
        IrregArray3D(const index_type x1, const index_type x2, const IrregArray1D<size_type, AForeign> &sizes, const index_type z1, const index_type z2, const value_type &ini)
            : first1_(0), first2_(0), first3_(0), last1_(-1), last2_(-1), last3_(-1), isIrreg_(false), arr_(nullptr)
        {
            initArray(x1, x2, sizes, z1, z2, ini);
        }

        /*! \brief general constructor for irregularly shaped matrices and all special cases

            \tparam      AForeign   allocator used by the input data structure

            \param[in]    x1       starting index in dimension 1
            \param[in]    x2       last index in dimension 1
            \param[in]    sizes1   number of elements in dimension 2 as function of the index in dimension 1
            \param[in]    sizes2   number of elements in dimension 3 as function of the indices in the lower dimensions
         */
        template<class AForeign = allocator_size_type>
        IrregArray3D(const index_type x1, const index_type x2, const IrregArray1D<size_type, AForeign> &sizes1, const IrregArray2D<size_type, AForeign> &sizes2)
            : first1_(0), first2_(0), first3_(0), last1_(-1), last2_(-1), last3_(-1), isIrreg_(false), arr_(nullptr)
        {
            initArray(x1, x2, sizes1, sizes2);
        }
        /*! \brief general constructor for irregularly shaped matrices and all special cases

            \tparam      AForeign   allocator used by the input data structure

            \param[in]    x1       starting index in dimension 1
            \param[in]    x2       last index in dimension 1
            \param[in]    sizes1   number of elements in dimension 2 as function of the index in dimension 1
            \param[in]    sizes2   number of elements in dimension 3 as function of the indices in the lower dimensions
            \param[in]    ini      initialization value for the data array elements
         */
        template<class AForeign = allocator_size_type>
        IrregArray3D(const index_type x1, const index_type x2, const IrregArray1D<size_type, AForeign> &sizes1, const IrregArray2D<size_type, AForeign> &sizes2, const value_type &ini)
            : first1_(0), first2_(0), first3_(0), last1_(-1), last2_(-1), last3_(-1), isIrreg_(false), arr_(nullptr)
        {
            initArray(x1, x2, sizes1, sizes2, ini);
        }

        /*! \brief  initialize from a std::initializer list
                    \code
                        arr = {
                                  {
                                      {value000, value001, ...},
                                      {value010, value011, ...}
                                      ...
                                  },
                                  {
                                      {value100, value101, ...},
                                      {value110, value111, ...}
                                      ...
                                  }
                                  ...
                              };
                    \endcode

            \param[in]   ini  initialization value
         */
        IrregArray3D(const std::initializer_list<std::initializer_list<std::initializer_list<value_type> > > &ini)
            : first1_(0), first2_(0), first3_(0), last1_(-1), last2_(-1), last3_(-1), isIrreg_(false), arr_(nullptr)
        {
            typedef std::initializer_list<value_type>         l1_ilist_type;
            typedef std::initializer_list<l1_ilist_type>      l2_ilist_type;
            typedef const value_type*                     l1_ilist_ptr_type;
            typedef const l1_ilist_type*                  l2_ilist_ptr_type;
            typedef const l2_ilist_type*                  l3_ilist_ptr_type;

            // check whether the list constains empty elements in any but the innermost nesting level
            if (ini.size() == 0)
            {
                return;
            }
            else
            {
                for (l3_ilist_ptr_type kPtr = ini.begin(); kPtr != ini.end(); ++kPtr)
                {
                    if (kPtr->size() == 0)
                    {
                        return;
                    }
                }
            }

            // check whether the array has a simple n x m x k shape or an irregular shape, where m or k differ between array stripes
            size_type prev_length = (ini.begin())->size();
            for (l3_ilist_ptr_type iPtrL3 = ini.begin(); iPtrL3 != ini.end(); ++iPtrL3)
            {
                if (prev_length != iPtrL3->size())
                {
                    isIrreg_ = true;
                    break;
                }
                prev_length = iPtrL3->size();
            }
            if (!isIrreg_)
            {
                size_type prev_length = ini.begin()->begin()->size();
                for (l3_ilist_ptr_type iPtrL3 = ini.begin(); iPtrL3 != ini.end(); ++iPtrL3)
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
            }

            // allocate memory depending on the array shape
            if (isIrreg_)
            {
                IrregArray1D<size_type, allocator_size_type> irregArr1((index_type)0, (index_type)ini.size() - 1, (size_t)0);
                size_t i = 0;
                for (l3_ilist_ptr_type iPtrL3 = ini.begin(); iPtrL3 != ini.end(); ++iPtrL3)
                {
                    irregArr1(i) = iPtrL3->size();
                    i++;
                }

                IrregArray2D<size_type, allocator_size_type> irregArr2((index_type)0, (index_type)ini.size() - 1, irregArr1, (size_t)0);
                i = 0;
                for (l3_ilist_ptr_type iPtrL3 = ini.begin(); iPtrL3 != ini.end(); ++iPtrL3)
                {
                    size_t j = 0;
                    for (l2_ilist_ptr_type iPtrL2 = iPtrL3->begin(); iPtrL2 != iPtrL3->end(); ++iPtrL2)
                    {
                        irregArr2(i, j) = iPtrL2->size();
                        j++;
                    }
                    i++;
                }

                initArray(0, ini.size() - 1, irregArr1, irregArr2);
            }
            else
            {
                initArray(0, static_cast<index_type>(ini.size()) - 1,
                          0, static_cast<index_type>((ini.begin())->size()) - 1,
                          0, static_cast<index_type>((ini.begin())->begin()->size()) - 1);
            }

            // assign the values from the input initializer_list
            for (size_type k = 0; k < getLength1(); ++k)
            {
                l3_ilist_ptr_type iPtr3 = (ini.begin() + static_cast<std::ptrdiff_t>(k));
                for (size_type j = 0; j < getLength2(k); ++j)
                {
                    l2_ilist_ptr_type iPtr2 = (iPtr3->begin() + static_cast<std::ptrdiff_t>(j));
                    for (size_type i = 0; i < getLength3(k, j); ++i)
                    {
                        l1_ilist_ptr_type iPtr1 = (iPtr2->begin() + static_cast<std::ptrdiff_t>(i));
                        operator()(k, j, i) = *iPtr1;
                    }
                }
            }
        }

        /*! \brief  initialize IrregArray3D<value_type> from IrregArray3D<TForeign>,
                    used by the copy constructor and by the explicit conversion operator,
                    TForeign must be convertible to value_type via static_cast

            \tparam      TForeign   data type stored by the input data structure
            \tparam      AForeign   allocator used by the input data structure

            \param[in]    ini    array to be copied
         */
        template <typename TForeign, class AForeign = std::allocator<TForeign> >
        void initArray(const IrregArray3D<TForeign, AForeign> &ini)
        {
            if (ini.arr_ != nullptr)
            {
                if (ini.isIrreg_)
                {
                    initArray(ini.first1_, ini.last1_, ini.irreg_.irreg_, ini.irreg_);
                }
                else
                {
                    initArray(ini.first1_, ini.last1_, ini.first2_, ini.last2_, ini.first3_, ini.last3_);
                }
            }

            if (arr_ != nullptr)
            {
                for (index_type i = first1_; i <= last1_; i++)
                {
                    for (index_type j = first2_; j <= getLast2(i); j++)
                    {
                        for (index_type k = first3_; k <= getLast3(i, j); k++)
                        {
                            arr_[i][j][k] = static_cast<value_type>(ini(i, j, k));
                        }
                    }
                }
            }
        }
        /*! \brief   assign the book keeping data and allocate memory for the
                     book keeping and data arrays

            \param[in]    x1      first index in dimension 1
            \param[in]    x2      last index in dimension 1
            \param[in]    y1      first index in dimension 2
            \param[in]    y2      last index in dimension 2
            \param[in]    z1      first index in dimension 3
            \param[in]    z2      last index in dimension 3
         */
        void initArray(const index_type x1, const index_type x2, const index_type y1, const index_type y2, const index_type z1, const index_type z2)
        {
            isIrreg_   = false;
            first1_    = x1;
            first2_    = y1;
            first3_    = z1;
            last1_     = x2;
            last2_     = y2;
            last3_     = z2;
            allocateArray();
        }
        /*! \brief   assign the book keeping data and allocate memory for the book keeping
                     and data arrays, assign an initial value to each array element

            \param[in]    x1      first index in dimension 1
            \param[in]    x2      last index in dimension 1
            \param[in]    y1      first index in dimension 2
            \param[in]    y2      last index in dimension 2
            \param[in]    z1      first index in dimension 3
            \param[in]    z2      last index in dimension 3
            \param[in]    ini     initialization value for the data array elements
         */
        void initArray(const index_type x1, const index_type x2, const index_type y1, const index_type y2, const index_type z1, const index_type z2, const value_type &ini)
        {
            isIrreg_ = false;
            initArray(x1, x2, y1, y2, z1, z2);
            for (index_type i = x1; i <= x2; i++)
            {
                for (index_type j = y1; j <= y2; j++)
                {
                    for (index_type k = z1; k <= z2; k++)
                    {
                        arr_[i][j][k] = ini;
                    }
                }
            }
        }
        /*! \brief   constructor helper function, assign the book keeping data and allocate memory for
                     the book keeping and data arrays

            \tparam      AForeign   allocator used by the input data structure

            \param[in]    x1      first index in dimension 1
            \param[in]    x2      last index in dimension 1
            \param[in]    sizes1  number of array elements in dimension 2 as a
                                  function of the index in dimension 1
            \param[in]    sizes2  number of array elements in dimension 3 as a
                                  function of the indices in the lower dimensions
         */
        template<class AForeign = allocator_size_type>
        void initArray(const index_type x1, const index_type x2, const IrregArray1D<size_type, AForeign> &sizes1, const IrregArray2D<size_type, AForeign> &sizes2)
        {
            isIrreg_   = true;
            first1_    = x1;
            first2_    = x1;
            first3_    = x1;
            last1_     = x2;
            last2_     = 0;
            last3_     = 0;
            IrregArray1D<size_type, allocator_size_type> irregVec(x1, x2);
            for (index_type i = x1; i <= x2; i++)
            {
                irregVec(i) = sizes1(i);
            }
            irreg_.initArray(getBegin1(), getLast1(), irregVec, 0);
            for (index_type i = getFirst1(); i <= getLast1(); i++)
            {
                for (index_type j = getFirst2(); j <= getLast2(i); j++)
                {
                    irreg_(i, j) = sizes2(i, j);
                }
            }
            allocateArray();
        }
        /*! \brief   constructor helper function, assign the book keeping data and allocate memory for
                     the book keeping and data arrays

            \tparam      AForeign   allocator used by the input data structure

            \param[in]    sizes2  number of array elements in dimensions, 1, 2, 3 as a
                                  function of the indices in the lower dimensions
            \param[in]    ini     initialization value for the array elements
         */
        template<class AForeign = allocator_size_type>
        void initArray(const IrregArray2D<size_type, AForeign> &sizes2, const value_type &ini)
        {
            initArray(sizes2);
            for (index_type i = getFirst1(); i <= getLast1(); i++)
            {
                for (index_type j = getFirst2(); j <= getLast2(i); j++)
                {
                    for (index_type k = getFirst3(); k <= getLast3(i, j); k++)
                    {
                        arr_[i][j][k] = ini;
                    }
                }
            }
        }
        /*! \brief   constructor helper function, assign the book keeping data and allocate memory for
                     the book keeping and data arrays

            \tparam      AForeign   allocator used by the input data structure

            \param[in]    sizes2  number of array elements in dimensions, 1, 2, 3 as a
                                  function of the indices in the lower dimensions
         */
        template<class AForeign = allocator_size_type>
        void initArray(const IrregArray2D<size_type, AForeign> &sizes2)
        {
            isIrreg_   = true;
            first1_    = sizes2.getFirst1();
            first2_    = sizes2.getFirst1();
            first3_    = sizes2.getFirst1();
            last1_     = sizes2.getLast1();
            last2_     = 0;
            last3_     = 0;
            IrregArray1D<size_type, allocator_size_type> irregVec(first1_, last1_);
            for (index_type i = first1_; i <= last1_; i++)
            {
                irregVec(i) = sizes2.getLength2(i);
            }
            irreg_.initArray(getBegin1(), getLast1(), irregVec, 0);
            for (index_type i = getFirst1(); i <= getLast1(); i++)
            {
                for (index_type j = getFirst2(); j <= getLast2(i); j++)
                {
                    irreg_(i, j) = sizes2(i, j);
                }
            }
            allocateArray();
        }
        /*! \brief   constructor helper function, assign the book keeping data and allocate memory for
                     the book keeping and data arrays, assign an initial value to each array element

            \tparam      AForeign   allocator used by the input data structure

            \param[in]    x1      first index in dimension 1
            \param[in]    x2      last index in dimension 1
            \param[in]    sizes1  number of array elements in dimension 2 as a
                                  function of the index in dimension 1
            \param[in]    sizes2  number of array elements in dimension 3 as a
                                  function of the indices in the lower dimensions
            \param[in]    ini     initialization value for the data array elements
         */
        template<class AForeign = allocator_size_type>
        void initArray(const index_type x1, const index_type x2, const IrregArray1D<size_type, AForeign> &sizes1, const IrregArray2D<size_type, AForeign> &sizes2, const value_type &ini)
        {
            initArray(x1, x2, sizes1, sizes2);
            for (index_type i = x1; i <= x2; i++)
            {
                for (index_type j = x1; j <= getLast2(i); j++)
                {
                    for (index_type k = x1; k <= getLast3(i, j); k++)
                    {
                        arr_[i][j][k] = ini;
                    }
                }
            }
        }
        /*! \brief   constructor helper function, assign the book keeping data and allocate memory for
                     the book keeping and data arrays, assign an initial value to each array element

            \tparam      AForeign   allocator used by the input data structure

            \param[in]    x1      first index in dimension 1
            \param[in]    x2      last index in dimension 1
            \param[in]    sizes2  number of array elements in dimension 3 as a
                                  function of the indices in the lower dimensions,
                                  and in dimension2 accordingly from its irreg member
            \param[in]    ini     initialization value for the data array elements
         */
        template<class AForeign = allocator_size_type>
        void initArray(const index_type x1, const index_type x2, const IrregArray2D<size_type, AForeign> &sizes2, const value_type &ini)
        {
            initArray(x1, x2, sizes2.irreg_, sizes2, ini);
        }
        /*! \brief   constructor helper function assign the book keeping data and allocate memory for
                     the book keeping and data arrays, assign an initial value to each array element

            \tparam      AForeign   allocator used by the input data structure

            \param[in]    x1      first index in dimension 1
            \param[in]    x2      last index in dimension 1
            \param[in]    sizes   number of array elements in dimension 2 as a
                                  function of the index in dimension 1
            \param[in]    z1      first index in dimension 3
            \param[in]    z2      last index in dimension 3
            \param[in]    ini     initialization value for the data array elements
         */
        template<class AForeign = allocator_size_type>
        void initArray(const index_type x1, const index_type x2, const IrregArray1D<size_type, AForeign> &sizes, const index_type z1, const index_type z2, const value_type &ini)
        {
            initArray(x1, x2, sizes, z1, z2);
            for (index_type i = getFirst1(); i <= getLast1(); i++)
            {
                for (index_type j = getFirst2(); j <= getLast2(i); j++)
                {
                    for (index_type k = getFirst3(); k <= getLast3(i, j); k++)
                    {
                        arr_[i][j][k] = ini;
                    }
                }
            }
        }
        /*! \brief   constructor helper function assign the book keeping data and allocate memory for
                     the book keeping and data arrays, assign an initial value to each array element

            \tparam      AForeign   allocator used by the input data structure

            \param[in]    x1      first index in dimension 1
            \param[in]    x2      last index in dimension 1
            \param[in]    y1      first index in dimension 2
            \param[in]    y2      last index in dimension 2
            \param[in]    sizes   number of array elements in dimension 2 as a
                                  function of the indices in the lower dimensions
            \param[in]    ini     initialization value for the data array elements
         */
        template<class AForeign = allocator_size_type>
        void initArray(const index_type x1, const index_type x2, const index_type y1, const index_type y2, const IrregArray1D<size_type, AForeign> &sizes, const value_type &ini)
        {
            initArray(x1, x2, y1, y2, sizes);
            for (index_type i = getFirst1(); i <= getLast1(); i++)
            {
                for (index_type j = getFirst2(); j <= getLast2(i); j++)
                {
                    for (index_type k = getFirst3(); k <= getLast3(i, j); k++)
                    {
                        arr_[i][j][k] = ini;
                    }
                }
            }
        }
        /*! \brief   assign the book keeping data and allocate memory for the
                     book keeping and data arrays

            \tparam      AForeign   allocator used by the input data structure

            \param[in]    x1      first index in dimension 1
            \param[in]    x2      last index in dimension 1
            \param[in]    y1      first index in dimension 2
            \param[in]    y2      last index in dimension 2
            \param[in]    sizes   number of array elements in dimension 2 as a
                                  function of the indices in the lower dimensions
         */
        template<class AForeign = allocator_size_type>
        void initArray(const index_type x1, const index_type x2, const index_type y1, const index_type y2, const IrregArray1D<size_type, AForeign> &sizes)
        {
            isIrreg_   = true;
            first1_    = x1;
            first2_    = y1;
            first3_    = first1_;
            last1_     = x2;
            last2_     = y2;
            last3_     = 0;
            irreg_.initArray(x1, x2, y1, y2, 0);
            for (index_type i = x1; i <= x2; i++)
            {
                for (index_type j = y1; j <= y2; j++)
                {
                    irreg_(i, j) = sizes(j);
                }
            }
            allocateArray();
        }
        /*! \brief   assign the book keeping data and allocate memory for the
                     book keeping and data arrays

            \tparam      AForeign   allocator used by the input data structure

            \param[in]    x1      first index in dimension 1
            \param[in]    x2      last index in dimension 1
            \param[in]    sizes   number of array elements in dimension 2
                                  as a function of the index in dimension 1
            \param[in]    z1      first index in dimension 3
            \param[in]    z2      last index in dimension 3
         */
        template<class AForeign = allocator_size_type>
        void initArray(index_type x1, index_type x2, const IrregArray1D<size_type, AForeign> &sizes, index_type z1, index_type z2)
        {
            isIrreg_ = true;
            const index_type dim3 = z2 - z1 + 1;
            first1_  = x1;
            first2_  = first1_;
            first3_  = z1;
            last1_   = x2;
            last2_   = 0;
            last3_   = z2;
            irreg_.initArray(x1, x2, sizes, dim3);
            for (index_type i = x1; i <= x2; i++)
            {
                for (index_type j = first2_; j < first2_ + sizes[i]; j++)
                {
                    irreg_(i, j) = dim3;
                }
            }
            allocateArray();
        }

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
        //! number of elements in dimension 2 at lower dimension index x
        //! \param[in]    x       index in dimension 1 for which the array length in dimension 2 is requested
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
                std::fprintf(stderr, "Error in %s::getLength3(): called without argument for an irregular array.", typeid(*this).name());
                throw InternalError("Extent of a variable dimension in an irregular array requested without specifying the index/indices in the lower dimension(s).");
            }
        }
        //! number of elements in dimension 3 at lower dimension indices x, y
        //! \param[in]    x       index in dimension 1
        //! \param[in]    y       index in dimension 2
        size_type getLength3(const index_type x, const index_type y) const
        {
            return (isIrreg_ ? irreg_(x, y) : (last3_ - first3_ + 1));
        }

        //! get the first index of dimension 1
        index_type getBegin1() const { return first1_; }
        //! get the first index of dimension 2
        index_type getBegin2() const { return first2_; }
        //! get the first index of dimension 3
        index_type getBegin3() const { return first3_; }

        //! get the first index of dimension 1
        index_type getFirst1() const { return first1_; }
        //! get the first index of dimension 2
        index_type getFirst2() const { return first2_; }
        //! get the first index of dimension 3
        index_type getFirst3() const { return first3_; }

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
        //! get the last index of dimension 2 at lower dimension index x
        //! \param[in]    x       index in dimension 1 for which the last index in dimension 2 is requested
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
        //! get the last index of dimension 3 for lower dimension indices x, y
        //! \param[in]    x       index in dimension 1
        //! \param[in]    y       index in dimension 2
        index_type getLast3(const index_type x, const index_type y) const
        {
            return (isIrreg_ ? (irreg_(x, y) + first3_ - 1) : last3_);
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
        //! get the index one past the last valid index of dimension 2 at lower dimension index x
        //! \param[in]    x       index in dimension 1 for which the last index in dimension 2 is requested
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
        //! get the index one past the last valid index of dimension 3 for lower dimension indices x, y
        //! \param[in]    x       index in dimension 1
        //! \param[in]    y       index in dimension 2
        index_type getEnd3(const index_type x, const index_type y) const
        {
            return (isIrreg_ ? (irreg_(x, y) + first3_) : last3_ + 1);
        }

        //! returns the amount of memory occupied by the object
        size_type getSize() const
        {
            // the book-keeping data
            size_type size = sizeof(*this);
            size += irreg_.getSize() - sizeof(IrregArray2D<size_type>);
            // the storage data
            size += getLength1() * sizeof(value_type**);
            for (index_type i = getFirst1(); i <= getLast1(); i++)
            {
                size += getLength2(i) * sizeof(value_type*);
                for (index_type j = getFirst2(); j <= getLast2(i); j++)
                {
                    size += getLength3(i, j) * sizeof(value_type);
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
                    nelements_ += getLength3(i, j);
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
                            std::allocator_traits<allocator_type>::destroy(allocator_, &this->operator()(i, j, k));
                        }
                    }
                }

                for (index_type i = getFirst1(); i <= getLast1(); i++)
                {
                    for (index_type j = getFirst2(); j <= getLast2(i); ++j)
                    {
                        arr_[i][j] += first3_;
                        std::allocator_traits<allocator_type>::deallocate(allocator_, arr_[i][j], getLength3(i, j));
                    }
                    arr_[i] += first2_;
                    std::allocator_traits<p_allocator_type>::deallocate(p_allocator_, arr_[i], getLength2(i));
                }
                arr_ += first1_;
                std::allocator_traits<p2_allocator_type>::deallocate(p2_allocator_, arr_, getLength1());

                arr_ = nullptr;
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
            last1_  = -1;
            last2_  = -1;
            last3_  = -1;
        }
        //! destructor
        ~IrregArray3D()
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

            \brief helper functor used to implement index operator *this[i][j][k]
         */
        class Proxy
        {
            private:
                value_type** array_;
                class Proxy2
                {
                    private:
                        value_type* array_;
                    public:
                        Proxy2(value_type* array) : array_(array) { }

                        value_type &operator[](index_type z)
                        {
                            return array_[z];
                        }
                        const value_type &operator[](index_type z) const
                        {
                            return array_[z];
                        }
                };
            public:
                Proxy(value_type** array) : array_(array) { }

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
        /*! \brief   index operator[x][y][z] for array-like usage realized with chained functors "Proxy*"
            \param[in]   x   index in dimension 1
            \cond
            \param[in]   y   index in dimension 2
            \param[in]   z   index in dimension 3
            \endcond
         */
        Proxy operator[](index_type x)
        {
            return Proxy(arr_[x]);
        }
        /*! \brief   index operator[x][y][z] for array-like usage realized with chained functors "Proxy*"
            \param[in]   x   index in dimension 1
            \cond
            \param[in]   y   index in dimension 2
            \param[in]   z   index in dimension 3
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
         */
        value_type &operator()(const index_type x, const index_type y, const index_type z)
        {
            return arr_[x][y][z];
        }
        /*! \brief   index operator(x,y,z) for functor-like usage
            \param[in]   x   index in dimension 1
            \param[in]   y   index in dimension 2
            \param[in]   z   index in dimension 3
         */
        value_type &operator()(const index_type x, const index_type y, const index_type z) const
        {
            return arr_[x][y][z];
        }
        /*! \brief index operator(x,y) for functor-like usage
            \param[in]   x   index in dimension 1
            \param[in]   y   index in dimension 2
         */
        value_type *operator()(const index_type x, const index_type y)
        {
            return arr_[x][y];
        }
        /*! \brief   const index operator(x,y) for functor-like usage
            \param[in]   x   index in dimension 1
            \param[in]   y   index in dimension 2
         */
        const value_type *operator()(const index_type x, const index_type y) const
        {
            return arr_[x][y];
        }
        /*! \brief   index operator(x) for functor-like usage
            \param[in]   x   index in dimension 1
         */
        value_type **operator()(const index_type x)
        {
            return arr_[x];
        }
        /*! \brief   const index operator(x) for functor-like usage
            \param[in]   x   index in dimension 1
         */
        const value_type **operator()(const index_type x) const
        {
            return arr_[x];
        }

        //! get a pointer to the data array
        const value_type*** getArray() const
        {
            return const_cast<const value_type***>(arr_);
        }
        //! get a pointer to the data array
        value_type*** getArray()
        {
            return arr_;
        }

        //! get a pointer to the data array
        const value_type*** data() const
        {
            return const_cast<const value_type***>(arr_);
        }
        //! get a pointer to the data array
        value_type*** data()
        {
            return arr_;
        }

        /*! \brief   assign data from a vector to an array stripe in dimension 3
            \param[in]   x      index in dimension 1
            \param[in]   y      index in dimension 2
            \param[in]   data   data array to be assigned to the array elements
                                arr[x][x][:]
         */
        void assign_vector(const index_type x, const index_type y, const value_type *data)
        {
            for (index_type z = getFirst3(); z <= getLast3(x, y); z++)
            {
                arr_[x][y][z] = data[z];
            }
        }

        /*! \brief   stream operator for convenient printing of the array to an ostream
            \param[in]   output   output stream for the array contents
            \param[in]   ten      the array to be printed
         */
        friend std::ostream &operator<<(std::ostream &output, const IrregArray3D &ten)
        {
            // safe the current ostream format for restoring it after output insertion
            std::ios  state(NULL);
            state.copyfmt(output);

            if (ten.arr_ != nullptr)
            {
                for (index_type i = ten.first1_; i <= ten.last1_; i++)
                {
                    output << "------------------------------" << std::endl << i << std::endl;
                    for (index_type j = ten.first2_; j <= ten.getLast2(i); j++)
                    {
                        for (index_type k = ten.first3_; k <= ten.getLast3(i, j); k++)
                        {
                            if (std::is_floating_point<value_type>::value)
                            {
                                output << std::setw(10) << std::fixed << std::setprecision(2) << ten(i, j, k);
                            }
                            else
                            {
                                output << ten(i, j, k);
                            }
                            if (k != ten.getLast3(i, j))
                            {
                                output << " ";
                            }
                        }
                        output << std::endl;
                    }
                    output << std::endl;
                }
            }

            // restore the original ostream format
            output.copyfmt(state);

            return output;
        }

        /*! \brief  move assignment operator this = rhs

            \param[in]   rhs   right-hand side of the assignment
         */
        IrregArray3D &operator=(IrregArray3D &&rhs) noexcept
        {
            if (this != &rhs)
            {
                // free any previously allocated memory
                deallocateArray();
                // move the data from rhs
                first1_       = rhs.first1_;
                first2_       = rhs.first2_;
                first3_       = rhs.first3_;
                last1_        = rhs.last1_;
                last2_        = rhs.last2_;
                last3_        = rhs.last3_;
                isIrreg_      = rhs.isIrreg_;
                arr_          = rhs.arr_;
                irreg_        = std::move(rhs.irreg_);
                allocator_    = std::move(rhs.allocator_);
                p_allocator_  = std::move(rhs.p_allocator_);
                p2_allocator_ = std::move(rhs.p2_allocator_);
                // set rhs to a valid default state
                rhs.first1_   =  0;
                rhs.first2_   =  0;
                rhs.first3_   =  0;
                rhs.last1_    = -1;
                rhs.last2_    = -1;
                rhs.last3_    = -1;
                rhs.isIrreg_  = false;
                rhs.arr_      = nullptr;
            }
            return *this;
        }

        /*! \brief   assignment operator *this = rhs

            \param[in]   rhs   right-hand side object of the assignment
         */
        IrregArray3D &operator=(const IrregArray3D &rhs)
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
        IrregArray3D &operator=(const IrregArray3D<value_type, AForeign> &rhs)
        {
            deallocateArray();
            initArray(rhs);
            return *this;
        }

        /*! \brief  conversion operator (explicit to avoid unintended implicit conversion), usable via
                    IrregArray3D<TForeign, AForeign> newArr = static_cast<IrregArray3D<TForeign, AForeign> >(IrregArray3D<TCurrent, ACurrent> oldArr);

            \tparam      TForeign   data type stored by the input data structure
            \tparam      AForeign   allocator used by the input data structure
         */
        template<typename TForeign, class AForeign = std::allocator<TForeign> >
        explicit
        operator IrregArray3D<TForeign, AForeign>() const
        {
            IrregArray3D<TForeign, AForeign> result;
            result.initArray(*this);
            return result;
        }

        //! returns true if the array is irregular
        inline bool isIrreg() const { return isIrreg_; }

    private:

        //! \brief allocate memory for the array using the preset array geometry
        void allocateArray()
        {
            if (last1_ < first1_)
            {
                std::fprintf(stderr, "Error in %s::allocateArray(): called with illegal bounds begin1 > end1.", typeid(*this).name());
                throw InternalError("Trying to create an irregular array with illegal bounds begin1 > end1");
            }
            if (!isIrreg_)
            {
                if (last2_ < first2_)
                {
                    std::fprintf(stderr, "Error in %s::allocateArray(): called with illegal bounds begin2 > end2.", typeid(*this).name());
                    throw InternalError("Trying to create an irregular array with illegal bounds begin2 > end2");
                }
                if (last3_ < first3_)
                {
                    std::fprintf(stderr, "Error in %s::allocateArray(): called with illegal bounds begin3 > end3.", typeid(*this).name());
                    throw InternalError("Trying to create an irregular array with illegal bounds begin3 > end3");
                }
            }

            try
            {
                arr_ = std::allocator_traits<p2_allocator_type>::allocate(p2_allocator_, getLength1());
                arr_ = arr_ - getBegin1();
                for (index_type i = getFirst1(); i <= getLast1(); ++i)
                {
                    arr_[i]  = std::allocator_traits<p_allocator_type>::allocate(p_allocator_, getLength2(i));
                    arr_[i] -= getFirst2();
                }
                for (index_type i = getFirst1(); i <= getLast1(); i++)
                {
                    for (index_type j = getFirst2(); j <= getLast2(i); ++j)
                    {
                        arr_[i][j]  = std::allocator_traits<allocator_type>::allocate(allocator_, getLength3(i, j));
                        arr_[i][j] -= getFirst3();
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
                        std::allocator_traits<allocator_type>::construct(allocator_, &this->operator()(i, j, k), value_type());
                    }
                }
            }
        }
};

} // end namespace gmx

#endif
