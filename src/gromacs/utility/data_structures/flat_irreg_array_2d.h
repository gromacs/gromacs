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

    \brief FlatIrregArray2D irregularly shaped 2-dimensional array
           internally stored in an 1D-array

    \author R. Thomas Ullmann <tullman@gwdg.de>
 */
#ifndef GMX_UTILITY_DATASTRUCTURES_FLAT_IRREG_ARRAY_2D_H
#define GMX_UTILITY_DATASTRUCTURES_FLAT_IRREG_ARRAY_2D_H

#include "gromacs/utility/data_structures/irreg_array_1d.h"

namespace gmx
{


/*! \class FlatIrregArray2D

     \brief irregularly shaped 2D array

     this array has a nonuniform second dimension, i.e. instead
     of a regular array with constant dimensions M x N, an
     irregular array with dimensions M x N(M).

     \verbatim

              dim2
              ----------------->
            / xxxxxxxxxxxxxxxx
     dim1  / x
          / xxx
         / xxxxxx
        / xxx
       V xxxxxxxx

     \endverbatim

     example application intrinsic energy matrix indexed by
     site, form, where the number of forms differs between sites

    \ingroup module_utility
    \inpublicapi

    \tparam   T           data type to be stored
    \tparam   Allocator   allocator to be used in creating arrays storing T
                          allocators for storing additional book-keeping data
                          are derived from this allocator via rebind if necessary
 */
template<class T = size_t, class Allocator = std::allocator<T> >
class FlatIrregArray2D
{
    public:
        /*  adopted these typedefs from AnBe's TriangularArray, for consistency with the other data
            structures, can be retrieved from outside to make sure that the intended types are used
            when accessing or manipulating the array */
        //! the data type stored in the array
        typedef T value_type; //! the data type stored in the array
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
    private:
        //! first index of dimension 1
        index_type          first1_;
        //! first index of dimension 2
        index_type          first2_;
        //! ending index of dimension 1
        index_type          last1_;
        //! ending index of dimension 2
        index_type          last2_;
        //! total number of array elements arr[:][:]
        size_type           nelements_;
        //! irreg_[i] = x is the index of the first actual array element arr_[x] of the
        //! stripe FlatIrregArray2D[i][:] in dimension 2 of the of the represented 2D array */
        size_1d_array_type   irreg_;
        //! flag says that dim2(i) = last2_(i) - first2_(i) + 1 is not a constant
        bool                 isIrreg_;
        //! array storing the the actual data
        value_type          *arr_;
        //! the allocator used to allocate arr_
        allocator_type       allocator_;
        //! let other FlatIrregArray2D instantiations access private members (for the explicit conversion operator)
        template<class U, typename AU> friend class FlatIrregArray2D;
        //! let higher-dimensional FlatIrregArrayNDs access private members of their book-keeping arrays
        template<class V, typename AV> friend class FlatIrregArray3D;
        //! let higher-dimensional FlatIrregArrayNDs access private members of their book-keeping arrays
        template<class W, typename AW> friend class FlatIrregArray4D;
    public:

        /*! \brief default constructor without memory allocation

            The initArray member functions allow for memory allocation
            and content initialization.
         */
        FlatIrregArray2D()
            : first1_(0), first2_(0), last1_(-1), last2_(-1),
              nelements_(0), isIrreg_(false), arr_(nullptr)
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

        /*! \brief  copy constructor with allocation and initialization
            \param[in]   ini   the template to be copied
         */
        FlatIrregArray2D(const FlatIrregArray2D &ini)
            : FlatIrregArray2D()
        {
            initArray(ini);
        }

        /*! \brief  copy constructor with allocation and initialization

            \tparam      AForeign   allocator used by the input data structure

            \param[in]   ini   the template to be copied
         */
        template<class AForeign = allocator_type>
        FlatIrregArray2D(const FlatIrregArray2D<value_type, AForeign> &ini)
            : FlatIrregArray2D()
        {
            initArray(ini);
        }

        //! move constructor
        //! \param[in]    ini   source to be moved
        FlatIrregArray2D(FlatIrregArray2D &&ini) noexcept
            : first1_(ini.first1_), first2_(ini.first2_),
              last1_(ini.last1_), last2_(ini.last2_),
              nelements_(ini.nelements_),
              irreg_(std::move(ini.irreg_)), isIrreg_(ini.isIrreg_),
              arr_(ini.arr_), allocator_(std::move(ini.allocator_))
        {
            // set ini to a valid default state
            ini.first1_    =  0;
            ini.first2_    =  0;
            ini.last1_     = -1;
            ini.last2_     = -1;
            ini.nelements_ =  0;
            ini.isIrreg_   = false;
            ini.arr_       = nullptr;
        }

        /*! \brief   constructor with memory allocation but without initialization

            \param[in]   x1   first index for dimension 1
            \param[in]   x2   last index for dimension 1
            \param[in]   y1   first index for dimension 2
            \param[in]   y2   last index for dimension 2
         */
        FlatIrregArray2D(const index_type x1, const index_type x2, const index_type y1, const index_type y2)
            : FlatIrregArray2D()
        {
            allocateArray(x1, x2, y1, y2);
        }
        /*! \brief  constructor for a regularly shaped array with allocation and initialization

            \tparam      AForeign   allocator used by the input data structure

            \param[in]   x1   first index for dimension 1
            \param[in]   x2   last index for dimension 1
            \param[in]   y1   first index for dimension 2
            \param[in]   y2   last index for dimension 2
            \param[in]   ini  initial value for all array elements
         */
        FlatIrregArray2D(const index_type x1, const index_type x2, const index_type y1, const index_type y2, const value_type &ini)
            : FlatIrregArray2D()
        {
            initArray(x1, x2, y1, y2, ini);
        }

        /*! \brief  constructor with memory allocation

            \tparam      AForeign   allocator used by the input data structure

            \param[in]   x1      first index for dimension 1
            \param[in]   x2      last index for dimension 1
            \param[in]   sizes   sizes[i] specifies the number of array elements [j] in arr[i]
         */
        template<class AForeign = allocator_size_type>
        FlatIrregArray2D(const index_type x1, const index_type x2, const IrregArray1D<size_type, AForeign> &sizes)
            : FlatIrregArray2D()
        {
            allocateArray(x1, x2, sizes);
        }
        /*! \brief  constructor with memory allocation and content initialization

            \tparam      AForeign   allocator used by the input data structure

            \param[in]   x1      first index for dimension 1
            \param[in]   x2      last index for dimension 1
            \param[in]   sizes   sizes[i] specifies the number of array elements [j] in arr[i]
            \param[in]   ini     initial value for all array elements
         */
        template<class AForeign = allocator_size_type>
        FlatIrregArray2D(const index_type x1, const index_type x2, const IrregArray1D<size_type, AForeign> &sizes, const value_type &ini)
            : FlatIrregArray2D()
        {
            allocateArray(x1, x2, sizes);
            for (size_type i = 0; i < nelements_; ++i)
            {
                arr_[i] = ini;
            }
        }

        /*! \brief  constructor with memory allocation

            \tparam      AForeign   allocator used by the input data structure

            \param[in]   sizes   sizes[i] specifies the number of array elements [j] in arr[i]
         */
        template<class AForeign = allocator_size_type>
        FlatIrregArray2D(const IrregArray1D<size_type, AForeign> &sizes)
            : FlatIrregArray2D()
        {
            allocateArray(sizes.getFirst1(), sizes.getLast1(), sizes);
        }
        /*! \brief  constructor with memory allocation and content initialization

            \tparam      AForeign   allocator used by the input data structure

            \param[in]   sizes   sizes[i] specifies the number of array elements [j] in arr[i]
            \param[in]   ini     initial value for all array elements
         */
        template<class AForeign = allocator_size_type>
        FlatIrregArray2D(const IrregArray1D<size_type, AForeign> &sizes, const value_type &ini)
            : FlatIrregArray2D()
        {
            allocateArray(sizes.getFirst1(), sizes.getLast1(), sizes);
            for (size_type i = 0; i < nelements_; ++i)
            {
                arr_[i] = ini;
            }
        }

        /*! \brief  constructor with memory allocation

            \tparam      AForeign   allocator used by the input data structure

            \param[in]   sizes   sizes[i] specifies the number of array elements [j] in arr[i]
         */
        template<class AForeign = allocator_size_type>
        FlatIrregArray2D(const std::vector<size_type, AForeign> &sizes)
            : FlatIrregArray2D()
        {
            allocateArray(0, static_cast<index_type>(sizes.size()) - 1, sizes);
        }
        /*! \brief  constructor with memory allocation and content initialization

            \tparam      AForeign   allocator used by the input data structure

            \param[in]   sizes   sizes[i] specifies the number of array elements [j] in arr[i]
            \param[in]   ini     initial value for all array elements
         */
        template<class AForeign = allocator_size_type>
        FlatIrregArray2D(const std::vector<size_type, AForeign> &sizes, const value_type &ini)
            : FlatIrregArray2D()
        {
            allocateArray(0, static_cast<index_type>(sizes.size()) - 1, sizes);
            for (size_type i = 0; i < nelements_; ++i)
            {
                arr_[i] = ini;
            }
        }

        /*! \brief  constructor with memory allocation and content initialization

            \tparam      AForeign   allocator used by the input data structure

            \param[in]   x1      first index for dimension 1
            \param[in]   x2      last index for dimension 1
            \param[in]   sizes   sizes[i] specifies the number of array elements arr[i][j] in arr[i]
            \param[in]   ini     initial value for all array elements
         */
        template<class AForeign = allocator_size_type>
        FlatIrregArray2D(const index_type x1, const index_type x2, const std::vector<size_type, AForeign> &sizes, const value_type &ini)
            : FlatIrregArray2D()
        {
            initArray(x1, x2, IrregArray1D<size_type, allocator_size_type>(x1, x2, sizes), ini);
        }

        /*! \brief  constructor with memory allocation but without content initialization

            \tparam      AForeign   allocator used by the input data structure

            \param[in]   x1      first index for dimension 1
            \param[in]   x2      last index for dimension 1
            \param[in]   sizes   sizes[i] specifies the number of array elements arr[i][j] in arr[i]
         */
        template<class AForeign = allocator_size_type>
        FlatIrregArray2D(const index_type x1, const index_type x2, const std::vector<size_type, AForeign> &sizes)
            : FlatIrregArray2D()
        {
            initArray(x1, x2, IrregArray1D<size_type, allocator_size_type>(x1, x2, sizes));
        }

        /*! \brief  initialize from a std::initializer list
                    \code
                        arr = { {value00, value01, ...}, {value10, value11, ...} ...};
                    \endcode

            \param[in]   ini  initialization value
         */
        FlatIrregArray2D(const std::initializer_list<std::initializer_list<value_type> > &ini)
            : FlatIrregArray2D()
        {
            typedef std::initializer_list<value_type>         l1_ilist_type;
            typedef const value_type*                     l1_ilist_ptr_type;
            typedef const l1_ilist_type*                  l2_ilist_ptr_type;

            // check whether the list contains empty elements in any but the innermost nesting level
            if (ini.size() == 0)
            {
                return;
            }

            // check whether the array has a simple n x m shape or an irregular shape, where m differs between array stripes
            size_type prev_length = (ini.begin())->size();
            for (l2_ilist_ptr_type iPtrL2 = ini.begin(); iPtrL2 != ini.end(); ++iPtrL2)
            {
                if (prev_length != iPtrL2->size())
                {
                    isIrreg_ = true;
                }
                prev_length = iPtrL2->size();
            }

// The properties of the initializer list are not tracked correctly by the clang static analyzer
// the condition in the loop above is recognized as false by clang for a non-empty input list
// below the loop length of the list is also incorrectly recorded as zero so resulting in the
// incorrect finding that no array is allocated. Strangely, this doesn't happen for all array classes.
// https://bugs.llvm.org/show_bug.cgi?id=39042
#ifdef __clang_analyzer__
            isIrreg_ = true;
#endif

            // allocate memory depending on the array shape
            if (isIrreg_)
            {
                IrregArray1D<size_type, allocator_size_type> irregArr1(0, static_cast<index_type>(ini.size()) - 1, 0);
                size_t i = 0;
                for (l2_ilist_ptr_type iPtrL2 = ini.begin(); iPtrL2 != ini.end(); ++iPtrL2)
                {
                    irregArr1[i] = iPtrL2->size();
                    i++;
                }

                initArray(0, static_cast<index_type>(ini.size()) - 1, irregArr1);
            }
            else
            {
                initArray(0, static_cast<index_type>(ini.size()) - 1, 0, ini.begin()->size() - 1);
            }

            // assign the values from the input initializer_list
            for (size_type j = 0; j < getLength1(); ++j)
            {
                l2_ilist_ptr_type iPtr2 = (ini.begin() + static_cast<std::ptrdiff_t>(j));
                for (size_type i = 0; i < getLength2(j); ++i)
                {
                    l1_ilist_ptr_type iPtr1 = (iPtr2->begin() + static_cast<std::ptrdiff_t>(i));
                    operator()(j, i) = *iPtr1;
                }
            }
        }

        /*! \brief  initialize FlatIrregArray2D<value_type> from FlatIrregArray2D<TForeign>,
                    used by the copy constructor and by the explicit conversion operator,
                    TForeign must be convertible to value_type via static_cast

            \tparam      TForeign   data type stored by the input data structure
            \tparam      AForeign   allocator used by the input data structure

            \param[in]   ini   the template to be copied
         */
        template <typename TForeign, class AForeign = std::allocator<TForeign> >
        void initArray(const FlatIrregArray2D<TForeign, AForeign> &ini)
        {
            if (ini.arr_ != nullptr)
            {
                deallocateArray();
                irreg_      = ini.irreg_;
                isIrreg_    = ini.isIrreg_;
                first1_     = ini.first1_;
                last1_      = ini.last1_;
                first2_     = ini.first2_;
                last2_      = ini.last2_;
                nelements_  = ini.nelements_;
                allocateArray(nelements_);
                for (size_type i = 0; i < nelements_; ++i)
                {
                    arr_[i] = static_cast<value_type>(ini.arr_[i]);
                }
            }
        }

        /*! \brief  initialize a regular special case of FlatIrregArray2D
            \param[in]   x1    first index for dimension 1
            \param[in]   x2    last index for dimension 1
            \param[in]   y1    first index for dimension 2
            \param[in]   y2    last index for dimension 2
         */
        void initArray(const index_type x1, const index_type x2,
                       const index_type y1, const index_type y2)
        {
            allocateArray(x1, x2, y1, y2);
        }

        /*! \brief  initialize a regular special case of FlatIrregArray2D
            \param[in]   x1    first index for dimension 1
            \param[in]   x2    last index for dimension 1
            \param[in]   y1    first index for dimension 2
            \param[in]   y2    last index for dimension 2
            \param[in]   ini   initialization value
         */
        void initArray(const index_type x1, const index_type x2,
                       const index_type y1, const index_type y2,
                       const value_type &ini)
        {
            allocateArray(x1, x2, y1, y2);
            for (size_type i = 0; i < nelements_; ++i)
            {
                arr_[i] = ini;
            }
        }

        /*! \brief  allocate memory for an FlatIrregArray2D and initialize its values

            \tparam      AForeign   allocator used by the input data structure

            \param[in]   x1       first index for dimension 1 and dimension 2
            \param[in]   x2       last index for dimension 1
            \param[in]   sizes    vector storing the array lengths in dimension 2
            \param[in]   ini      initialization value in second variant
         */
        template<class AForeign = allocator_size_type>
        void initArray(const index_type x1, const index_type x2,
                       const IrregArray1D<size_type, AForeign> &sizes,
                       const value_type &ini)
        {
            allocateArray(x1, x2, sizes);
            for (size_type i = 0; i < nelements_; ++i)
            {
                arr_[i] = ini;
            }
        }

        /*! \brief  allocate memory for an FlatIrregArray2D

            \tparam      AForeign   allocator used by the input data structure

            \param[in]   x1       first index for dimension 1 and dimension 2
            \param[in]   x2       last index for dimension 1
            \param[in]   sizes    vector storing the array lengths in dimension 2
         */
        template<class AForeign = allocator_size_type>
        void initArray(const index_type x1, const index_type x2,
                       const std::vector<size_type, AForeign> &sizes)
        {
            initArray(x1, x2, IrregArray1D<size_type, allocator_size_type>(sizes));
        }
        /*! \brief  allocate memory for an FlatIrregArray2D

            \tparam      AForeign   allocator used by the input data structure

            \param[in]   x1       first index for dimension 1 and dimension 2
            \param[in]   x2       last index for dimension 1
            \param[in]   sizes    vector storing the array lengths in dimension 2
            \param[in]   ini      initialization value in second variant
         */
        template<class AForeign = allocator_size_type>
        void initArray(const index_type x1, const index_type x2,
                       const std::vector<size_type, AForeign> &sizes,
                       const value_type &ini)
        {
            initArray(x1, x2, IrregArray1D<size_type, allocator_size_type>(sizes));
            for (size_type i = 0; i < nelements_; ++i)
            {
                arr_[i] = ini;
            }
        }

        /*! \brief  allocate memory for an FlatIrregArray2D

            \tparam      AForeign   allocator used by the input data structure

            \param[in]   sizes    vector storing the array lengths in dimension 2
         */
        template<class AForeign = allocator_size_type>
        void initArray(const std::vector<size_type, AForeign> &sizes)
        {
            initArray(0, static_cast<index_type>(sizes.size()) - 1,
                      IrregArray1D<size_type, allocator_size_type>(sizes));
        }
        /*! \brief  allocate memory for an FlatIrregArray2D

            \tparam      AForeign   allocator used by the input data structure

            \param[in]   sizes    vector storing the array lengths in dimension 2
            \param[in]   ini      initialization value in second variant
         */
        template<class AForeign = allocator_size_type>
        void initArray(const std::vector<size_type, AForeign> &sizes, const value_type &ini)
        {
            initArray(0, static_cast<index_type>(sizes.size()) - 1,
                      IrregArray1D<size_type, allocator_size_type>(sizes));
            for (size_type i = 0; i < nelements_; ++i)
            {
                arr_[i] = ini;
            }
        }

        /*! \brief  allocate memory for an FlatIrregArray2D

            \tparam      AForeign   allocator used by the input data structure

            \param[in]   sizes    vector storing the array lengths in dimension 2
         */
        template<class AForeign = allocator_size_type>
        void initArray(const IrregArray1D<size_type, AForeign> &sizes)
        {
            initArray(sizes.getFirst1(), sizes.getLast1(),
                      IrregArray1D<size_type, allocator_size_type>(sizes));
        }
        /*! \brief  allocate memory for an FlatIrregArray2D

            \tparam      AForeign   allocator used by the input data structure

            \param[in]   sizes    vector storing the array lengths in dimension 2
            \param[in]   ini      initialization value in second variant
         */
        template<class AForeign = allocator_size_type>
        void initArray(const IrregArray1D<size_type, AForeign> &sizes, const value_type &ini)
        {
            initArray(sizes.getFirst1(), sizes.getLast1(),
                      IrregArray1D<size_type, allocator_size_type>(sizes));
            for (size_type i = 0; i < nelements_; ++i)
            {
                arr_[i] = ini;
            }
        }

        /*! \brief  allocate memory for an FlatIrregArray2D

            \tparam      AForeign   allocator used by the input data structure

            \param[in]   x1       first index for dimension 1 and dimension 2
            \param[in]   x2       last index for dimension 1
            \param[in]   &sizes   vector storing the array lengths in dimension 2
         */
        template<class AForeign = allocator_size_type>
        void initArray(const index_type x1, const index_type x2,
                       const IrregArray1D<size_type, AForeign> &sizes)
        {
            allocateArray(x1, x2, sizes);
        }

        //! returns true if the array is irregular
        inline bool isIrreg() const { return isIrreg_; }

        //! number of elements in dimension 1
        inline size_type getLength1() const { return (last1_ - first1_ + 1); }
        //! number of elements in dimension 2 (regular array)
        inline size_type getLength2() const
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
        /*! \brief  number of elements in dimension 2
            \param[in]    x       index in dimension 1 for which the array length in dimension 2 is requested */
        inline size_type getLength2(const index_type x) const
        {
            if (isIrreg_)
            {
                // the last element of the array stripe of interest is located
                // one element before the first element of the next array stripe
                // the length of the stripe is then given by the difference between
                // the indices of the first elements of the current and next stripes
                if (x == last1_)
                {
                    return nelements_ - irreg_(last1_);
                }
                else
                {
                    return irreg_(x + 1) - irreg_(x);
                }
            }
            else
            {
                return last2_ - first2_ + 1;
            }
        }

        //! get the first index of dimension 1
        inline index_type getBegin1() const { return first1_; }
        //! get the first index of dimension 2
        inline index_type getBegin2() const { return first2_; }

        //! get the first index of dimension 1
        inline index_type getFirst1() const { return first1_; }
        //! get the first index of dimension 2
        inline index_type getFirst2() const { return first2_; }

        //! get the last index of dimension 1
        inline index_type getLast1() const { return last1_; }

        //! get the end as in std containers, one past the last valid index of dimension 1
        inline index_type getEnd1() const { return last1_ + 1; }

        //! get the last index of dimension 2
        inline index_type getLast2() const
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
        /*! \brief  get the last index of dimension 2 at lower dimension index x

            \param[in]    x       index in dimension 1 for which the last index in dimension 2 is requested
         */
        inline index_type getLast2(const index_type x) const
        {
            if (isIrreg_)
            {
                // the last element of the array stripe of interest is located
                // one element before the first element of the next array stripe
                // the end of the stripe is then given by the difference between
                // the indices of the first and last elements of the current stripe
                // plus the starting index of the array stripe
                if (x == last1_)
                {
                    return (nelements_ - 1) - irreg_(last1_) + first2_;
                }
                else
                {
                    return (irreg_(x + 1) - 1) - irreg_(x) + first2_;
                }
            }
            else
            {
                return last2_;
            }
        }

        //! get the end as in std containers, one past the last valid index of dimension 2
        inline index_type getEnd2() const
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
        /*! \brief get the end as in std containers, one past the last valid index of dimension 2

            \param[in]    x       index in dimension 1 for which the last index in dimension 2 is requested
         */
        inline index_type getEnd2(const index_type x) const
        {
            if (isIrreg_)
            {
                // the last element of the array stripe of interest is located
                // one element before the first element of the next array stripe
                // the end of the stripe is then given by the difference between
                // the indices of the first and last elements of the current stripe
                // plus the starting index of the array stripe
                if (x == last1_)
                {
                    return nelements_ - irreg_(last1_) + first2_;
                }
                else
                {
                    return irreg_(x + 1) - irreg_(x) + first2_;
                }
            }
            else
            {
                return last2_ + 1;
            }
        }

        /*! \brief determine the index of the array element FlatIrregArray2D(i, j) in the underlying 1D array arr
            this function needs to be public to be compatible to C++ standards < C++11, where nested classes
            can not access private members of the enclosing class. A friend declaration can also not be used
            because the friends are defined to be entities that are not members.

            size_type instead of index_type used intentionally here to use the maximum possible length

            \param[in]   x    index in dimension 1
            \param[in]   y    index in dimension 2
         */
        inline size_type getArrayIndex(const index_type x, const index_type y) const
        {
            if (isIrreg_)
            {
                return irreg_(x) + (y - first2_);
            }
            else
            {
                return (getLength2() * (x - first1_)) + (y - first2_);
            }
        }

        //! determine the amount of memory occupied by the object, usually not needed, therefore computed on the fly
        size_type getSize() const
        {
            size_type size  = sizeof(*this);
            size           += nelements_ * sizeof(value_type);
            size           += irreg_.getSize() - sizeof(size_1d_array_type);
            return size;
        }

        //! determine the number of data array elements
        size_type getNelements() const { return nelements_; }

        //! deallocate memory
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
            irreg_.deallocateArray();
            isIrreg_ = false;
            first1_  = 0;
            first2_  = 0;
            last1_   = -1;
            last2_   = -1;
        }
        //! destructor
        ~FlatIrregArray2D()
        {
            // deallocate arrays
            deallocateArray();
            // deallocate book-keeping array
            if (isIrreg_)
            {
                irreg_.deallocateArray();
            }
            first1_ = first2_ = last1_ = last2_ = 0;
        }

        // index operators

        //! doxygen can not handle nested classes
        /// \cond DEV
        /*!  \class Proxy

             \brief  proxy functor for index operator for bare array-like usage a[i][j]
                   as described by Mike Seymour on stackoverflow
         */
        class Proxy
        {
            public:
                /*!
                   \brief constructor

                   \param[in]   optr   the array containing the element to be accessed
                   \param[in]   x   index in dimension 1
                 */
                Proxy(FlatIrregArray2D* optr, index_type x) : optr_(optr), x_(x) { }
                /*!
                   \brief square bracket operator for bare array-like access

                   \param[in]   y   index in dimension 2
                 */
                value_type &operator[](index_type y)
                {
                    return optr_->operator()(x_, y);
                }
                /*!
                   \brief const square bracket operator for bare array-like access

                   \param[in]   y   index in dimension 2
                 */
                const value_type &operator[](index_type y) const
                {
                    return optr_->operator()(x_, y);
                }
            private:
                //! pointer to the array to be accessed
                FlatIrregArray2D            * optr_;
                //! index of in dimension 1
                const index_type              x_;
        };
        /// \endcond DEV
        /*! index operator[x] [y] for array-like usage
            \param[in]   x   index in dimension 1
            \cond
            \param[in]   y   index in dimension 2
            \endcond
         */
        Proxy operator[](index_type x)
        {
            return Proxy(this, x);
        }
        /*! \brief  index operator[x] [y] for array-like usage
            \param[in]   x   index in dimension 1
            \cond
            \param[in]   y   index in dimension 2
            \endcond
         */
        const Proxy operator[](index_type x) const
        {
            return Proxy(this, x);
        }

        /*! \brief  index operator(x,y) for functor-like usage
            \param[in]   x   index in dimension 1
            \param[in]   y   index in dimension 2
         */
        const value_type &operator()(const index_type x, const index_type y) const
        {
            return arr_[getArrayIndex(x, y)];
        }
        /*! \brief  index operator(x,y) for functor-like usage
            \param[in]   x   index in dimension 1
            \param[in]   y   index in dimension 2
         */
        value_type &operator()(const index_type x, const index_type y)
        {
            return arr_[getArrayIndex(x, y)];
        }

        /*! \brief  index operator(x) for functor-like usage
            \param[in]   x   index in dimension 1
         */
        const value_type* operator()(const index_type x) const
        {
            return &arr_[getArrayIndex(x, getBegin2())];
        }
        /*! \brief  index operator(x) for functor-like usage
            \param[in]   x   index in dimension 1
         */
        value_type* operator()(const index_type x)
        {
            return &arr_[getArrayIndex(x, getBegin2())];
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

        /*! \brief  output stream operator
            \param[in]   output   ostream in which the content is to be inserted
            \param[in]   mat      the object whose contents are to be inserted
         */
        friend std::ostream &operator<<(std::ostream &output, const FlatIrregArray2D &mat)
        {
            // safe the current ostream format for restoring it after output insertion
            std::ios  state(nullptr);
            state.copyfmt(output);

            if (mat.arr_ != nullptr)
            {
                constexpr size_t floatWidth = 10;
                constexpr size_t floatPrec  =  2;
                for (index_type i = mat.first1_; i <= mat.last1_; i++)
                {
                    for (index_type j = mat.first2_; j <= mat.getLast2(i); j++)
                    {
                        if (std::is_floating_point<value_type>::value)
                        {
                            output << std::setw(floatWidth) << std::fixed << std::setprecision(floatPrec) << mat(i, j);
                        }
                        else
                        {
                            output << mat(i, j);
                        }
                        if (j != mat.getLast2(i))
                        {
                            output << " ";
                        }
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
        FlatIrregArray2D &operator=(FlatIrregArray2D &&rhs) noexcept
        {
            if (this != &rhs)
            {
                // free any previously allocated memory
                deallocateArray();
                // move the data from rhs
                first1_    = rhs.first1_;
                first2_    = rhs.first2_;
                last1_     = rhs.last1_;
                last2_     = rhs.last2_;
                nelements_ = rhs.nelements_;
                isIrreg_   = rhs.isIrreg_;
                arr_       = rhs.arr_;
                irreg_     = std::move(rhs.irreg_);
                allocator_ = std::move(rhs.allocator_);
                // set rhs to a valid default state
                rhs.first1_    =  0;
                rhs.first2_    =  0;
                rhs.last1_     = -1;
                rhs.last2_     = -1;
                rhs.nelements_ =  0;
                rhs.isIrreg_   = false;
                rhs.arr_       = nullptr;
            }
            return *this;
        }

        /*! \brief  assignment operator this = rhs
            \param[in]   rhs   right-hand side of the assignment
         */
        FlatIrregArray2D &operator=(const FlatIrregArray2D &rhs)
        {
            if (this != &rhs)
            {
                deallocateArray();
                initArray(rhs);
            }
            return *this;
        }

        /*! \brief  assignment operator this = rhs

            \tparam      AForeign   allocator used by the input data structure

            \param[in]   rhs   right-hand side of the assignment
         */
        template<class AForeign = allocator_size_type>
        FlatIrregArray2D &operator=(const FlatIrregArray2D<value_type, AForeign> &rhs)
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
                FlatIrregArray2D<TForeign, AForeign> newArr =
                    static_cast<FlatIrregArray2D<TForeign, AForeign> >
                    (
                        FlatIrregArray2D<TCurrent, ACurrent> oldArr
                    );
            \endcode

            \tparam      TForeign   data type stored by the input data structure
            \tparam      AForeign   allocator used by the input data structure
         */
        template<typename TForeign, class AForeign = std::allocator<TForeign> >
        explicit
        operator FlatIrregArray2D<TForeign, AForeign>() const
        {
            FlatIrregArray2D<TForeign, AForeign> result;
            result.initArray(*this);
            return result;
        }

    private:

        /*! \brief  allocate memory for the array (regular shape)
            \param[in]   x1   first index for dimension 1
            \param[in]   x2   last index for dimension 1
            \param[in]   y1   first index for dimension 2
            \param[in]   y2   last index for dimension 2
         */
        void allocateArray(const index_type x1, const index_type x2,
                           const index_type y1, const index_type y2)
        {
            deallocateArray(); // prevent memory leaks
            if (x2 >= x1 && y2 >= y1)
            {
                isIrreg_   = false;
                first1_    = x1;
                last1_     = x2;
                first2_    = y1;
                last2_     = y2;
                const size_type dim1 = getLength1();
                const size_type dim2 = getLength2();
                nelements_ = dim1 * dim2;
                allocateArray(nelements_);
            }
            else
            {
                if (x1 > x2)
                {
                    std::string errorMessage(formatString("Error in %s::allocateArray(): called with illegal bounds first1 > last1.", typeid(*this).name()));
                    GMX_THROW(InvalidInputError(errorMessage));
                }
                if (y1 > y2)
                {
                    std::string errorMessage(formatString("Error in %s::allocateArray(): called with illegal bounds first2 > last2.", typeid(*this).name()));
                    GMX_THROW(InvalidInputError(errorMessage));
                }
            }
        }
        /*! \brief  allocate memory for the array (irregular shape)

            \tparam      AForeign   allocator used by the input data structure

            \param[in]   x1       first index for dimension 1 and dimension 2
            \param[in]   x2       last index for dimension 1
            \param[in]   &sizes   vector storing the array lengths in dimension 2
         */
        template<class AForeign = allocator_size_type>
        void allocateArray(const index_type x1, const index_type x2,
                           const IrregArray1D<size_type, AForeign> &sizes)
        {
            deallocateArray(); // prevent memory leaks
            bool length2AlwaysNonNegative = true;
            for (index_type i = x1; i <= x2; ++i)
            {
                // size_type may be signed or unsigned, quite compiler warnings
                // for the comparison < 0 if size_type is unsigned
                if (static_cast<std::make_signed<size_type>::type>(sizes(i)) < 0)
                {
                    length2AlwaysNonNegative = false;
                    std::string errorMessage(formatString("Error in %s::allocateArray(): called with illegal length %li < 0 for element %li\n",
                                                          typeid(*this).name(), static_cast<long int>(sizes(i)), static_cast<long int>(i)));
                    GMX_THROW(InvalidInputError(errorMessage));
                }
            }

            if (x2 >= x1 && length2AlwaysNonNegative)
            {
                isIrreg_   = true;
                first1_    = x1;
                first2_    = first1_;
                last1_     = x2;
                last2_     = 0;

                if (first1_ > last1_)
                {
                    std::string errorMessage(formatString("Error in %s::allocateArray(): called with illegal bounds first1 > last1.", typeid(*this).name()));
                    GMX_THROW(InvalidInputError(errorMessage));
                }

                irreg_.initArray(x1, x2);
                for (index_type i = x1; i <= x2; i++)
                {
                    irreg_(i) = (i == x1 ? 0 : irreg_(i-1) + sizes(i-1));
                }
                // the number of elements is given by the starting index of the last array stripe
                // FlatIrregArray2D[i][:] in dimension 2 within arr_ + the length of this stripe + 1
                const size_type n = irreg_(last1_) + sizes(last1_);

                allocateArray(n);
            }
            else
            {
                if (x2 < x1)
                {
                    std::string errorMessage("Error in %s::allocateArray(): called with illegal bounds first1 > last1.", typeid(*this).name());
                    GMX_THROW(InvalidInputError(errorMessage));
                }
                if (!length2AlwaysNonNegative)
                {
                    std::string errorMessage("Error in %s::allocateArray(): called with illegal bounds first2 > last2.", typeid(*this).name());
                    GMX_THROW(InvalidInputError(errorMessage));
                }
            }
        }
        /*! \brief  allocate memory for the array
            \param[in]   n  total number of array elements
         */
        void allocateArray(const size_type n)
        {
            nelements_ = n;
            if (nelements_ > 0)
            {
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
        }
};

} // end namespace gmx

#endif
