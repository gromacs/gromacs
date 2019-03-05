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

    \brief declares gmx::IrregArray2D, a 2D-array that can have irregular shape

    \author R. Thomas Ullmann <tullman@gwdg.de>
 */
#ifndef GMX_UTILITY_DATASTRUCTURES_IRREG_ARRAY_2D_H
#define GMX_UTILITY_DATASTRUCTURES_IRREG_ARRAY_2D_H

#include <iomanip>
#include <iostream>
#include <memory>
#include <stdexcept>
#include <typeinfo>
#include <vector>

#include "gromacs/utility/exceptions.h"
#include "gromacs/utility/stringutil.h"

// older(?) MS Visual C++ compilers don' know ssize_t
#if defined(_MSC_VER)
#include <BaseTsd.h>
typedef SSIZE_T ssize_t;
#endif

namespace gmx
{

/*! \class IrregArray2D

    \brief irregularly shaped 2-dimensional array

     This array can have a nonuniform second dimension, i.e.,
     instead of a regular array with constant dimensions M x N,
     an irregular array with dimensions M x N(M). The regular
     special case is also possible. isIrreg() returns true if
     the array has irregular shape and false if it has regular
     rectangular shape.

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

     example application: energy matrix indexed by site, form,
     where the number of forms differs between sites

    \ingroup module_utility
    \inpublicapi

    \tparam   T           data type to be stored
    \tparam   Allocator   allocator to be used in creating arrays storing T
                          allocators for storing additional book-keeping data
                          are derived from this allocator via rebind if necessary
    \tparam   isFlat      underlying array is 1-dimensional if true and multi-dimensional if false
 */
template<class T = size_t, class Allocator = std::allocator<T>, bool IsFlat = false>
class IrregArray2D
{
    public:
        //! the data type stored in the array
        typedef T value_type;
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
        typedef std::vector<size_type, allocator_size_type> size_1d_array_type;
    private:
        //! array length in dimension 1
        size_type               length1_   = 0;
        //! array length in dimension 2
        size_type               length2_   = 0;
        //! total number of array elements
        size_type               nelements_ = 0;
        //! number of array elements in dimension 2 as a function of the index in dimension 1
        size_1d_array_type      irreg_;
        //! flag says that length2(i) is not a constant
        bool                    isIrreg_   = false;
        //! array storing the the actual data, 1D if IsFlat, 2D if !isFlat
        typename std::conditional<IsFlat, value_type*, value_type**>::type arr_ = nullptr;
        //! the allocator_ used to allocate arr_
        p_allocator_type        p_allocator_;
        //! the allocator_ used to allocate arr_[i]
        allocator_type          allocator_;
    public:
        /*! \brief default constructor without memory allocation

            The initArray member functions allow for memory allocation
            and content initialization.
         */
        IrregArray2D()
            : length1_(0), length2_(0), nelements_(0), isIrreg_(false), arr_(nullptr)
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
            initArray(1, 1);
#endif      // end scan-build bug work-around
        }

        /*! \brief  copy constructor with allocation and initialization
            \param[in]   ini   the template to be copied
         */
        IrregArray2D(const IrregArray2D &ini)
            : IrregArray2D()
        {
            initArray(ini);
        }

        /*! \brief  templated copy constructor
            \tparam      AForeign   allocator used by the input data structure
            \tparam      FForeign   true if the input data structure uses a flat data array
            \param[in]   ini        the template to be copied
         */
        template<class AForeign, bool FForeign>
        IrregArray2D(const IrregArray2D<value_type, AForeign, FForeign> &ini)
            : IrregArray2D()
        {
            initArray(ini);
        }

        //! move constructor
        //! \param[in]    ini   source to be moved
        IrregArray2D(IrregArray2D &&ini) noexcept
            : length1_(ini.length1_), length2_(ini.length2_),
              nelements_(ini.nelements_),
              irreg_(std::move(ini.irreg_)), isIrreg_(ini.isIrreg_),
              arr_(ini.arr_), p_allocator_(std::move(ini.p_allocator_)),
              allocator_(std::move(ini.allocator_))
        {
            // set ini to a valid default state
            ini.length1_   = 0;
            ini.length2_   = 0;
            ini.nelements_ = 0;
            ini.isIrreg_   = false;
            ini.arr_       = nullptr;
        }

        /*! \brief  move assignment operator this = rhs
            \param[in]   rhs   right-hand side of the assignment
         */
        IrregArray2D &operator=(IrregArray2D &&rhs) noexcept
        {
            if (this != &rhs)
            {
                // free any previously allocated memory
                deallocateArray();
                // move the data from rhs
                length1_     = rhs.length1_;
                length2_     = rhs.length2_;
                nelements_   = rhs.nelements_;
                isIrreg_     = rhs.isIrreg_;
                arr_         = rhs.arr_;
                irreg_       = std::move(rhs.irreg_);
                allocator_   = std::move(rhs.allocator_);
                p_allocator_ = std::move(rhs.p_allocator_);
                // set rhs to a valid default state
                rhs.length1_   = 0;
                rhs.length2_   = 0;
                rhs.nelements_ = 0;
                rhs.isIrreg_   = false;
                rhs.arr_       = nullptr;
            }
            return *this;
        }

        /*! \brief  copy assignment operator this = rhs
            \param[in]   rhs   right-hand side of the assignment
         */
        IrregArray2D &operator=(const IrregArray2D &rhs)
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
            \tparam      FForeign   true if the input data structure uses a flat data array
            \param[in]   rhs   right-hand side of the assignment
         */
        template<class AForeign = allocator_size_type, bool FForeign>
        IrregArray2D &operator=(const IrregArray2D<value_type, AForeign, FForeign> &rhs)
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
                IrregArray2D<TForeign, AForeign> newArr =
                    static_cast<IrregArray2D<TForeign, AForeign, FForeign> >
                    (
                        IrregArray2D<TCurrent, ACurrent, FCurrent> oldArr
                    );
            \endcode

            \tparam      TForeign   data type stored by the input data structure
            \tparam      AForeign   allocator used by the input data structure
            \tparam      FForeign   true if the output data structure shall use a flat data array
         */
        template<typename TForeign, class AForeign = std::allocator<TForeign>, bool FForeign = false>
        explicit
        operator IrregArray2D<TForeign, AForeign, FForeign>() const
        {
            IrregArray2D<TForeign, AForeign, FForeign> result;
            result.initArray(*this);
            return result;
        }

        /*! \brief   construct a regularly shaped array
            \param[in]   l1   array length in dimension 1
            \param[in]   l2   array length in dimension 2
         */
        IrregArray2D(size_type l1, size_type l2)
            : IrregArray2D()
        {
            initArray(l1, l2);
        }
        /*! \brief  construct a regularly shaped array and initialize the array elements to ini
            \param[in]   l1   array length in dimension 1
            \param[in]   l2   array length in dimension 2
            \param[in]   ini  initial value for all array elements
         */
        IrregArray2D(size_type l1, size_type l2, const value_type &ini)
            : IrregArray2D()
        {
            initArray(l1, l2, ini);
        }

        /*! \brief  construct an irregularly shaped array

            \tparam      AForeign   allocator used by the input data structure

            \param[in]   l2    l2[i] specifies the number of array elements [j] in arr[i],
                               l2.size() is the array length in dimension 1
         */
        template<class AForeign = allocator_size_type>
        IrregArray2D(const std::vector<size_type, AForeign> &l2)
            : IrregArray2D()
        {
            initArray(l2);
        }
        /*! \brief  construct an irregularly shaped array and initialize the array elements to ini

            \tparam      AForeign   allocator used by the input data structure

            \param[in]   l2    l2[i] specifies the number of array elements [j] in arr[i],
                               l2.size() is the array length in dimension 1
            \param[in]   ini   initial value for all array elements
         */
        template<class AForeign = allocator_size_type>
        IrregArray2D(const std::vector<size_type, AForeign> &l2, const value_type &ini)
            : IrregArray2D()
        {
            initArray(l2, ini);
        }

        /*! \brief  initialize from a std::initializer list
                    \code
                        arr = { {value00, value01, ...}, {value10, value11, ...} ...};
                    \endcode

            \param[in]   ini  nested list of initialization values
         */
        IrregArray2D(const std::initializer_list<std::initializer_list<value_type> > &ini)
            : IrregArray2D()
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
            const size_type ini_length2 = ini.begin()->size();
            for (l2_ilist_ptr_type iPtrL2 = ini.begin(); iPtrL2 != ini.end(); ++iPtrL2)
            {
                if (ini_length2 != iPtrL2->size())
                {
                    isIrreg_ = true;
                    break;
                }
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
                size_1d_array_type l2(ini.size(), 0);
                size_t             i = 0;
                for (l2_ilist_ptr_type iPtrL2 = ini.begin(); iPtrL2 != ini.end(); ++iPtrL2)
                {
                    l2[i] = iPtrL2->size();
                    i++;
                }
                initArray(l2);
            }
            else
            {
                initArray(ini.size(), ini.begin()->size());
            }

            // assign the values from the input initializer_list
            for (size_type j = 0; j < length1(); ++j)
            {
                l2_ilist_ptr_type iPtr2 = ini.begin() + static_cast<std::ptrdiff_t>(j);
                for (size_type i = 0; i < length2(j); ++i)
                {
                    l1_ilist_ptr_type iPtr1 = iPtr2->begin() + static_cast<std::ptrdiff_t>(i);
                    operator()(j, i) = *iPtr1;
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
                    operator()(i, j) = ini;
                }
            }
        }

        /*! \brief  initialize IrregArray2D<value_type> from IrregArray2D<TForeign>,
                    used by the copy constructor and by the explicit conversion operator,
                    TForeign must be convertible to value_type via static_cast

            \tparam      TForeign   data type stored by the input data structure
            \tparam      AForeign   allocator used by the input data structure
            \tparam      FForeign   true if the input data structure uses a flat data array

            \param[in]   ini   array to be copied
         */
        template <typename TForeign, class AForeign, bool FForeign>
        void initArray(const IrregArray2D<TForeign, AForeign, FForeign> &ini)
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
                    initArray(l2);
                }
                else
                {
                    initArray(ini.length1(), ini.length2());
                }
            }
            if (ini.data() != nullptr && arr_ != nullptr)
            {
                for (size_type i = 0; i < ini.length1(); ++i)
                {
                    for (size_type j = 0; j < ini.length2(i); ++j)
                    {
                        operator()(i, j) = static_cast<value_type>(ini(i, j));
                    }
                }
            }
        }

        /*! \brief  initialize a regularly shaped array
            \param[in]   l1    array length in dimension 1
            \param[in]   l2    array length in dimension 2
         */
        void initArray(size_type l1, size_type l2)
        {
            deallocateArray();
            // otherwise, we simply keep the empty array
            if (l1 > 0 && l2 > 0)
            {
                isIrreg_   = false;
                length1_   = l1;
                length2_   = l2;
                nelements_ = l1 * l2;
                allocateArray();
            }
        }
        /*! \brief  initialize a regularly shaped array and assign ini to all array elements
            \param[in]   l1    array length in dimension 1
            \param[in]   l2    array length in dimension 2
            \param[in]   ini   initial value for the array elements
         */
        void initArray(size_type l1, size_type l2, const value_type &ini)
        {
            initArray(l1, l2);
            initArrayElements(ini);
        }

        /*! \brief  create an irregularly shaped array with underlying 2D data array

            \tparam      AForeign      allocator used by the input data structure
            \tparam      DummyIsFlat   dummy template parameter for SFINAE selection according to IsFlat

            \param[in]   l2   array length in dimension 2 as function of the index in dimension 1,
                              l2.size() is the array length in dimension 1
         */
        template<class AForeign = allocator_size_type, bool DummyIsFlat = IsFlat>
        typename std::enable_if<!DummyIsFlat, void>::type
        initArray(const std::vector<size_type, AForeign> &l2)
        {
            GMX_ASSERT(IsFlat == DummyIsFlat, "Explicit instantiation with DummyIsFlat != IsFlat prohibited.");
            deallocateArray();
            // simply keep an empty array otherwise
            if (!l2.empty())
            {
                isIrreg_ = true;
                length1_ = l2.size();
                length2_ = 0;
                irreg_.resize(l2.size());
                for (size_type i = 0; i < length1(); i++)
                {
                    irreg_[i] = l2[i];
                }
                computeNelements();
                allocateArray();
            }
        }
        /*! \brief  create an irregularly shaped array with underlying 1D data array

            \tparam      AForeign      allocator used by the input data structure
            \tparam      DummyIsFlat   dummy template parameter for SFINAE selection according to IsFlat

            \param[in]   l2   array length in dimension 2 as function of the index in dimension 1
                              l2.size() is the array length in dimension 1
         */
        template<class AForeign = allocator_size_type, bool DummyIsFlat = IsFlat>
        typename std::enable_if<DummyIsFlat, void>::type
        initArray(const std::vector<size_type, AForeign> &l2)
        {
            GMX_ASSERT(IsFlat == DummyIsFlat, "Explicit instantiation with DummyIsFlat != IsFlat prohibited.");
            deallocateArray();

            // otherwise, just keep the empty array
            if (!l2.empty())
            {
                isIrreg_ = true;
                length1_ = l2.size();
                length2_ = 0;
                irreg_.resize(l2.size());
                for (size_type i = 0; i < l2.size(); i++)
                {
                    irreg_[i] = (i == 0 ? 0 : irreg_[i-1] + l2[i-1]);
                }
                // the number of elements is given by the starting index of the last array stripe
                // FlatIrregArray2D[i][:] in dimension 2 within arr_ + the length of this stripe
                nelements_ = irreg_[l2.size() - 1] + l2[l2.size() - 1];
                allocateArray();
            }
        }
        /*! \brief  create an irregularly shaped array and initialize its elements to ini (IsFlat == false)

            \tparam      AForeign   allocator used by the input data structure

            \param[in]   l2       array lengths in dimension 2 as function of the index in dimension 1,
                                  l2.size() is the array length in dimension 1
            \param[in]   ini      initialization value in second variant
         */
        template<class AForeign = allocator_size_type>
        void initArray(const std::vector<size_type, AForeign> &l2, const value_type &ini)
        {
            initArray(l2);
            initArrayElements(ini);
        }

        //! returns true if the array has irregular shape
        inline bool isIrreg() const noexcept { return isIrreg_; }

        //! returns true if the underlying data array is 1-D (flat), false if it is 2D
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
        //! \brief  number of elements in dimension 2 (IsFlat = false)
        //! \tparam      DummyIsFlat   dummy template parameter for SFINAE selection according to IsFlat
        //! \param[in]   x             index in dimension 1 for which the array length in dimension 2 is requested
        template<bool DummyIsFlat = IsFlat>
        inline typename std::enable_if<!DummyIsFlat, size_type>::type
        length2(index_type x) const noexcept
        {
            GMX_ASSERT(IsFlat == DummyIsFlat, "Explicit instantiation with DummyIsFlat != IsFlat prohibited.");
            return (isIrreg_ ? irreg_[x] : length2_);
        }
        //! \brief  number of elements in dimension 2 (IsFlat = true)
        //! \tparam      DummyIsFlat   dummy template parameter for SFINAE selection according to IsFlat
        //! \param[in]   x             index in dimension 1 for which the array length in dimension 2 is requested
        template<bool DummyIsFlat = IsFlat>
        inline typename std::enable_if<DummyIsFlat, size_type>::type
        length2(index_type x) const noexcept
        {
            GMX_ASSERT(IsFlat == DummyIsFlat, "Explicit instantiation with DummyIsFlat != IsFlat prohibited.");
            if (!isIrreg_)
            {
                return length2_;
            }
            else
            {
                // the last element of the array stripe of interest is located
                // one element before the first element of the next array stripe
                // the length of the stripe is then given by the difference between
                // the indices of the first elements of the current and next stripes
                const index_type last1 = length1_ - 1;
                if (x == last1)
                {
                    return nelements_ - irreg_[last1];
                }
                else
                {
                    return irreg_[x + 1] - irreg_[x];
                }
            }
        }

        /*! \brief determine the index of an array element in the underlying 1D array arr_ (IsFlat = true)
            This function is public because it may be useful in client code, too.

            size_type instead of index_type returned intentionally to use the maximum possible length

            \tparam      DummyIsFlat   dummy template parameter for SFINAE selection according to IsFlat

            \param[in]   x    index in dimension 1
            \param[in]   y    index in dimension 2
         */
        template<bool DummyIsFlat = IsFlat>
        inline typename std::enable_if<DummyIsFlat, size_type>::type
        arrayIndex(index_type x, index_type y) const noexcept
        {
            GMX_ASSERT(IsFlat == DummyIsFlat, "Explicit instantiation with DummyIsFlat != IsFlat prohibited.");
            if (isIrreg_)
            {
                return irreg_[x] + y;
            }
            else
            {
                return length2() * x + y;
            }
        }

        //! returns the number of data array elements
        size_type nElements() const noexcept { return nelements_; }
        //! compute the number of data array elements
        size_type computeNelements() noexcept
        {
            // the storage data
            nelements_ = 0;
            for (size_type i = 0; i < length1(); ++i)
            {
                nelements_ += length2(i);
            }
            return nelements_;
        }

        //! deallocate memory
        void deallocateArray() noexcept
        {
            if (arr_ != nullptr)
            {
                destroyArrayElements();
                deallocateMemory();
            }
            irreg_.clear();
            isIrreg_   = false;
            length1_   = 0;
            length2_   = 0;
            nelements_ = 0;
        }
        //! destructor
        ~IrregArray2D() noexcept
        {
            deallocateArray();
        }

        //! doxygen can not handle nested classes
        /// \cond DEV
        /*!  \class Proxy

             \brief  proxy functor for index operator for bare array-like usage a[i][j]

             The trick here is to chain operator[] of the nested subclasses.
             The original idea was provided by Mike Seymor on Stack Overflow.
             I modified it here, so that the same implementation can be used
             for flat and non-flat data arrays.
         */
        class Proxy
        {
            public:
                /*!
                   \brief constructor

                   \param[in]   optr   the array containing the element to be accessed
                   \param[in]   x   index in dimension 1
                 */
                Proxy(IrregArray2D* optr, index_type x) noexcept : optr_(optr), x_(x) { }
                /*!
                   \brief square bracket operator for bare array-like access

                   \param[in]   y   index in dimension 2
                 */
                value_type &operator[](index_type y) noexcept
                {
                    return optr_->operator()(x_, y);
                }
                /*!
                   \brief const square bracket operator for bare array-like access

                   \param[in]   y   index in dimension 2
                 */
                const value_type &operator[](index_type y) const noexcept
                {
                    return optr_->operator()(x_, y);
                }
            private:
                //! pointer to the array to be accessed
                IrregArray2D                * optr_;
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
        Proxy operator[](index_type x) noexcept
        {
            return Proxy(this, x);
        }
        /*! \brief  index operator[x] [y] for array-like usage
            \param[in]   x   index in dimension 1
            \cond
            \param[in]   y   index in dimension 2
            \endcond
         */
        Proxy operator[](index_type x) const noexcept
        {
            return Proxy(this, x);
        }

        /*! \brief  index operator(x,y) for functor-like usage (IsFlat = false)
            \tparam      DummyIsFlat   dummy template parameter for SFINAE selection according to IsFlat
            \param[in]   x   index in dimension 1
            \param[in]   y   index in dimension 2
         */
        template<bool DummyIsFlat = IsFlat>
        typename std::enable_if<!DummyIsFlat, value_type &>::type
        operator()(index_type x, index_type y) noexcept
        {
            GMX_ASSERT(IsFlat == DummyIsFlat, "Explicit instantiation with DummyIsFlat != IsFlat prohibited.");
            return arr_[x][y];
        }
        /*! \brief  index operator(x,y) for functor-like usage (IsFlat = false)
            \tparam      DummyIsFlat   dummy template parameter for SFINAE selection according to IsFlat
            \param[in]   x   index in dimension 1
            \param[in]   y   index in dimension 2
         */
        template<bool DummyIsFlat = IsFlat>
        typename std::enable_if<!DummyIsFlat, const value_type &>::type
        operator()(index_type x, index_type y) const noexcept
        {
            GMX_ASSERT(IsFlat == DummyIsFlat, "Explicit instantiation with DummyIsFlat != IsFlat prohibited.");
            return arr_[x][y];
        }
        //! index operator(x) for functor-like usage (IsFlat = false)
        //! \tparam      DummyIsFlat   dummy template parameter for SFINAE selection according to IsFlat
        //! \param[in]   x   index in dimension 1
        template<bool DummyIsFlat = IsFlat>
        typename std::enable_if<!DummyIsFlat, value_type*>::type
        operator()(index_type x) noexcept
        {
            GMX_ASSERT(IsFlat == DummyIsFlat, "Explicit instantiation with DummyIsFlat != IsFlat prohibited.");
            return arr_[x];
        }
        //! \brief  index operator(x) for functor-like usage (IsFlat = false)
        //! \tparam      DummyIsFlat   dummy template parameter for SFINAE selection according to IsFlat
        //! \param[in]   x   index in dimension 1
        template<bool DummyIsFlat = IsFlat>
        typename std::enable_if<!DummyIsFlat, const value_type*>::type
        operator()(index_type x) const noexcept
        {
            GMX_ASSERT(IsFlat == DummyIsFlat, "Explicit instantiation with DummyIsFlat != IsFlat prohibited.");
            return arr_[x];
        }

        //! get a pointer to the data array (IsFlat = false)
        //! \tparam   DummyIsFlat   dummy template parameter for SFINAE selection according to IsFlat
        template<bool DummyIsFlat = IsFlat>
        typename std::enable_if<!DummyIsFlat, value_type**>::type
        data() noexcept
        {
            GMX_ASSERT(IsFlat == DummyIsFlat, "Explicit instantiation with DummyIsFlat != IsFlat prohibited.");
            return arr_;
        }
        //! get a pointer to the data array (IsFlat = false)
        //! \tparam   DummyIsFlat   dummy template parameter for SFINAE selection according to IsFlat
        template<bool DummyIsFlat = IsFlat>
        typename std::enable_if<!DummyIsFlat, const value_type**>::type
        data() const noexcept
        {
            GMX_ASSERT(IsFlat == DummyIsFlat, "Explicit instantiation with DummyIsFlat != IsFlat prohibited.");
            return const_cast<const value_type**>(arr_);
        }

        //! \brief  index operator(x,y) for functor-like usage (IsFlat = true)
        //! \tparam   DummyIsFlat   dummy template parameter for SFINAE selection according to IsFlat
        //! \param[in]   x   index in dimension 1
        //! \param[in]   y   index in dimension 2
        template<bool DummyIsFlat = IsFlat>
        typename std::enable_if<DummyIsFlat, value_type &>::type
        operator()(index_type x, index_type y) noexcept
        {
            GMX_ASSERT(IsFlat == DummyIsFlat, "Explicit instantiation with DummyIsFlat != IsFlat prohibited.");
            return arr_[arrayIndex(x, y)];
        }
        //! \brief  index operator(x,y) for functor-like usage (IsFlat = true)
        //! \tparam   DummyIsFlat   dummy template parameter for SFINAE selection according to IsFlat
        //! \param[in]   x   index in dimension 1
        //! \param[in]   y   index in dimension 2
        template<bool DummyIsFlat = IsFlat>
        typename std::enable_if<DummyIsFlat, const value_type &>::type
        operator()(index_type x, index_type y) const noexcept
        {
            GMX_ASSERT(IsFlat == DummyIsFlat, "Explicit instantiation with DummyIsFlat != IsFlat prohibited.");
            return arr_[arrayIndex(x, y)];
        }
        //! \brief  index operator(x) for functor-like usage (IsFlat = true)
        //! \tparam   DummyIsFlat   dummy template parameter for SFINAE selection according to IsFlat
        //! \param[in]   x   index in dimension 1
        template<bool DummyIsFlat = IsFlat>
        typename std::enable_if<DummyIsFlat, value_type*>::type
        operator()(index_type x) noexcept
        {
            GMX_ASSERT(IsFlat == DummyIsFlat, "Explicit instantiation with DummyIsFlat != IsFlat prohibited.");
            return &arr_[arrayIndex(x, 0)];
        }
        //! \brief  index operator(x) for functor-like usage (IsFlat = true)
        //! \tparam   DummyIsFlat   dummy template parameter for SFINAE selection according to IsFlat
        //! \param[in]   x   index in dimension 1
        template<bool DummyIsFlat = IsFlat>
        typename std::enable_if<DummyIsFlat, const value_type*>::type
        operator()(index_type x) const noexcept
        {
            GMX_ASSERT(IsFlat == DummyIsFlat, "Explicit instantiation with DummyIsFlat != IsFlat prohibited.");
            return &arr_[arrayIndex(x, 0)];
        }

        //! get a pointer to the data array (IsFlat = true)
        //! \tparam   DummyIsFlat   dummy template parameter for SFINAE selection according to IsFlat
        template<bool DummyIsFlat = IsFlat>
        typename std::enable_if<DummyIsFlat, value_type*>::type
        data() noexcept
        {
            GMX_ASSERT(IsFlat == DummyIsFlat, "Explicit instantiation with DummyIsFlat != IsFlat prohibited.");
            return arr_;
        }
        //! get a pointer to the data array (IsFlat = true)
        //! \tparam   DummyIsFlat   dummy template parameter for SFINAE selection according to IsFlat
        template<bool DummyIsFlat = IsFlat>
        typename std::enable_if<DummyIsFlat, const value_type*>::type
        data() const noexcept
        {
            GMX_ASSERT(IsFlat == DummyIsFlat, "Explicit instantiation with DummyIsFlat != IsFlat prohibited.");
            return const_cast<const value_type*>(arr_);
        }

        //! output stream operator
        //! \param[in]   output   ostream in which the content is to be inserted
        //! \param[in]   mat      the object whose contents are to be inserted
        friend std::ostream &operator<<(std::ostream &output, const IrregArray2D &mat)
        {
            // safe the current ostream format for restoring it after output insertion
            std::ios  state(nullptr);
            state.copyfmt(output);

            if (mat.arr_ != nullptr)
            {
                constexpr size_t floatWidth = 10;
                constexpr size_t floatPrec  =  2;
                for (size_type i = 0; i < mat.length1(); i++)
                {
                    for (size_type j = 0; j < mat.length2(i); j++)
                    {
                        if (std::is_floating_point<value_type>::value)
                        {
                            output << std::setw(floatWidth) << std::fixed << std::setprecision(floatPrec) << mat(i, j);
                        }
                        else
                        {
                            output << mat(i, j);
                        }
                        if (j < (mat.length2(i) - 1))
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

    private:
        //! allocate memory for the array given the array geometry preset by the public initArray functions that accept arguments
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
                if (static_cast<std::make_signed<size_type>::type>(length2(i)) < 0)
                {
                    std::string errorMessage(formatString("Error in %s::allocateArray(): called with illegal bounds length2 < 0 at i1 = %lu.", typeid(*this).name(), i));
                    GMX_THROW(InvalidInputError(errorMessage));
                }
            }
            allocateMemory();
            constructArrayElements();
        }

        //! allocate memory for the underlying 2D array if !IsFlat, only to be called by allocateArray()
        //! \tparam   DummyIsFlat   dummy template parameter for SFINAE selection according to IsFlat
        template<bool DummyIsFlat = IsFlat>
        typename std::enable_if<!DummyIsFlat, void>::type
        allocateMemory()
        {
            try
            {
                arr_  = std::allocator_traits<p_allocator_type>::allocate(p_allocator_, length1());
                for (size_type i = 0; i < length1(); ++i)
                {
                    arr_[i] = std::allocator_traits<allocator_type>::allocate(allocator_, length2(i));
                }
            }
            catch (const std::bad_alloc &)
            {
                std::fprintf(stderr, "Error in %s::allocateMemory(): could not allocate memory for the data array.", typeid(*this).name());
                throw;
            }
        }

        //! allocate memory for the underlying 1D array if IsFlat, only to be called by allocateArray()
        //! \tparam   DummyIsFlat   dummy template parameter for SFINAE selection according to IsFlat
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
                        std::allocator_traits<allocator_type>::construct(allocator_, &this->operator()(i, j), value_type());
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
                        std::allocator_traits<allocator_type>::destroy(allocator_, &this->operator()(i, j));
                    }
                }
            }
        }
        //! deallocate memory for the underlying 2D array if !IsFlat, only to be called by deallocateArray()
        //! \tparam   DummyIsFlat   dummy template parameter for SFINAE selection according to IsFlat
        template<bool DummyIsFlat = IsFlat>
        typename std::enable_if<!DummyIsFlat, void>::type
        deallocateMemory() noexcept
        {
            if (arr_ != nullptr)
            {
                for (size_type i = 0; i < length1(); ++i)
                {
                    std::allocator_traits<allocator_type>::deallocate(allocator_, arr_[i], length2(i));
                }
                std::allocator_traits<p_allocator_type>::deallocate(p_allocator_, arr_, length1());
                arr_ = nullptr;
            }
        }

        //! deallocate memory for the underlying 1D array if IsFlat, only to be called by deallocateArray()
        //! \tparam   DummyIsFlat   dummy template parameter for SFINAE selection according to IsFlat
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
