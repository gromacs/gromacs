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

    \brief IrregArray1D is a 1-dimensional array
           with custom begin and end index

    \author R. Thomas Ullmann <tullman@gwdg.de>
 */
#ifndef GMX_UTILITY_DATASTRUCTURES_IRREG_ARRAY_1D_H
#define GMX_UTILITY_DATASTRUCTURES_IRREG_ARRAY_1D_H

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


/*! \class IrregArray1D

    \brief 1-dimensional array with custom begin and end index

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
class IrregArray1D
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
    private:
        //! the allocator type to allocate size_type
        typedef typename std::allocator_traits<Allocator>::template rebind_alloc<size_type> allocator_size_type;
        //! first index of dimension 1
        index_type     first1_;
        //! ending index of dimension 1
        index_type     last1_;
        //! array storing the the actual data
        value_type    *arr_;
        //! the allocator used to allocate arr
        allocator_type allocator_;
        //! let other IrregArray1D instantiations access private members (for the explicit conversion operator)
        template<class S, typename AS> friend class IrregArray1D;
        //! let higher-dimensional irregular arrays access private members of their book-keeping arrays
        template<class U, typename AU> friend class IrregArray2D;
        //! let higher-dimensional irregular arrays access private members of their book-keeping arrays
        template<class V, typename AV> friend class IrregArray3D;
        //! let higher-dimensional irregular arrays access private members of their book-keeping arrays
        template<class W, typename AW> friend class IrregArray4D;
        //! let higher-dimensional irregular arrays access private members of their book-keeping arrays
        template<class X, typename AX> friend class FlatIrregArray2D;
        //! let higher-dimensional irregular arrays access private members of their book-keeping arrays
        template<class Y, typename AY> friend class FlatIrregArray3D;
        //! let higher-dimensional irregular arrays access private members of their book-keeping arrays
        template<class Z, typename AZ> friend class FlatIrregArray4D;
    public:
        /*! \brief default constructor without memory allocation, allocateArray has to be called for
                   allocating memory initArray for memory allocation and content initialization
         */
        IrregArray1D()
            : first1_(0), last1_(-1), arr_(nullptr) {}

        //! copy constructor
        //! \param[in]    ini   template to be copied
        IrregArray1D(const IrregArray1D &ini)
            : first1_(0), last1_(-1), arr_(nullptr)
        {
            initArray(ini);
        }

        /*! \brief copy constructor, allows for a different allocator of the input array

            \tparam      AForeign   allocator used by the input data structure

            \param[in]    ini   template to be copied
         */
        template<class AForeign = allocator_type>
        IrregArray1D(const IrregArray1D<value_type, AForeign> &ini)
            : first1_(0), last1_(-1), arr_(nullptr)
        {
            initArray(ini);
        }

        //! move constructor
        //! \param[in]    ini   source to be moved
        IrregArray1D(IrregArray1D &&ini) noexcept
            : first1_(ini.first1_), last1_(ini.last1_), arr_(ini.arr_), allocator_(std::move(ini.allocator_))
        {
            // set ini to a valid default state
            ini.first1_ =  0;
            ini.last1_  = -1;
            ini.arr_    = nullptr;
        }

        /*! \brief  initialize IrregArray1D<value_type> from IrregArray1D<TForeign>,
                    used by the copy constructor and by the explicit conversion operator,
                    TForeign must be convertible to value_type via static_cast

            \tparam      AForeign   allocator used by the input data structure

            \param[in]   ini   the template to be copied
         */
        template <typename TForeign, class AForeign = std::allocator<TForeign> >
        void initArray(const IrregArray1D<TForeign, AForeign> &ini)
        {
            deallocateArray();
            if (ini.getLast1() >= ini.getFirst1())
            {
                allocateArray(ini.getFirst1(), ini.getLast1());
            }
            if (arr_ != nullptr && ini.arr_ != nullptr)
            {
                for (index_type i = ini.getFirst1(); i <= ini.getLast1(); i++)
                {
                    arr_[i] = static_cast<value_type>(ini(i));
                }
            }
        }

        /*! \brief constructor with memory allocation but without initialization

            \param[in]   x1   first index for dimension 1
            \param[in]   x2   last index for dimension 1
         */
        IrregArray1D(const index_type x1, const index_type x2)
            : first1_(0), last1_(-1), arr_(nullptr)
        {
            allocateArray(x1, x2);
        }

        /*! \brief constructor with allocation and initialization
            \param[in]   x1   first index for dimension 1
            \param[in]   x2   last index for dimension 1
            \param[in]   ini  initial value for all array elements
         */
        IrregArray1D(const index_type x1, const index_type x2, const value_type &ini)
            : first1_(0), last1_(-1), arr_(nullptr)
        {
            initArray(x1, x2, ini);
        }

        /*! \brief constructor with allocation but without initialization

            \param[in]   x1        first index for dimension 1
            \param[in]   x2        last index for dimension 1
            \param[in]   ini_vec   initial values for all array elements
         */
        IrregArray1D(const index_type x1, const index_type x2, const value_type* const ini_vec)
            : first1_(0), last1_(-1), arr_(nullptr)
        {
            allocateArray(x1, x2);
            if (arr_ != nullptr && ini_vec != nullptr)
            {
                for (index_type i = getFirst1(); i <= getLast1(); i++)
                {
                    arr_[i] = ini_vec[i];
                }
            }
        }

        /*! \brief constructor with initilization values taken form an IrregArray1D with a possibly differing allocator

            \tparam      AForeign   allocator used by the input data structure

            \param[in]   x1        first index for dimension 1
            \param[in]   x2        last index for dimension 1
            \param[in]   ini_vec   vector to be copied to a IrregArray1D
         */
        template<class AForeign = allocator_type>
        IrregArray1D(const index_type x1, const index_type x2, const std::vector<value_type, AForeign> &ini_vec)
            : first1_(0), last1_(-1), arr_(nullptr)
        {
            if (ini_vec.size() != 0)
            {
                allocateArray(x1, x2);
                if (arr_ != nullptr)
                {
                    for (index_type i = x1; i <= x2; i++)
                    {
                        arr_[i] = ini_vec[i];
                    }
                }
            }
        }

        /*! \brief copy constructor from a std::vector

            \tparam      AForeign   allocator used by the input data structure

            \param[in]   ini_vec   vector to be copied to a IrregArray1D
         */
        template<class AForeign = allocator_type>
        IrregArray1D(const std::vector<value_type, AForeign> &ini_vec)
            : first1_(0), last1_(-1), arr_(nullptr)
        {
            if (ini_vec.size() != 0)
            {
                const index_type x2 = ini_vec.size() - 1;
                allocateArray(0, x2);
                for (index_type i = 0; i <= x2; ++i)
                {
                    arr_[i] = ini_vec[i];
                }
            }
        }

        /*! \brief  constructor with allocation but without initialization

            \tparam      AForeign   allocator used by the input data structure

            \param[in]   x1        first index for dimension 1
            \param[in]   x2        last index for dimension 1
            \param[in]   ini_vec   IrregArray1D to be copied
         */
        template<class AForeign = allocator_type>
        IrregArray1D(const index_type x1, const index_type x2, const IrregArray1D<value_type, AForeign> &ini_vec)
            : first1_(0), last1_(-1), arr_(nullptr)
        {
            if (ini_vec.arr_ != nullptr)
            {
                allocateArray(x1, x2);
                for (index_type i = x1; i <= x2; i++)
                {
                    arr_[i] = ini_vec[i];
                }
            }
        }

        /*! \brief  initialize from a std::initializer list \code arr = {value0, value1, ...}; \endcode

            \param[in]   ini  initialization value
         */
        IrregArray1D(const std::initializer_list<value_type> &ini)
            : first1_(0), last1_(-1), arr_(nullptr)
        {
            if (ini.size() != 0)
            {
                initArray(0, ini.size() - 1);
                for (size_t i = 0; i < ini.size(); ++i)
                {
                    arr_[i] = *(ini.begin() + static_cast<std::ptrdiff_t>(i));
                }
            }
        }

        /*! \brief  initialize a regular special case of IrregArray1D

            \param[in]   x1   first index for dimension 1
            \param[in]   x2   last index for dimension 1
         */
        void initArray(const index_type x1, const index_type x2)
        {
            allocateArray(x1, x2);
        }
        /*! \brief  initialize a regular special case of IrregArray1D

            \param[in]   x1   first index for dimension 1
            \param[in]   x2   last index for dimension 1
            \param[in]   ini  initialization value
         */
        void initArray(const index_type x1, const index_type x2, const value_type &ini)
        {
            initArray(x1, x2);
            for (index_type i = x1; i <= x2; i++)
            {
                arr_[i] = ini;
            }
        }

        //! for consistency with the higher dimensional arrays, 1D array is always regular
        bool isIrreg() const { return false; }

        //! number of elements in dimension 1
        size_type getLength1() const { return (last1_ - first1_ + 1); }
        //! number of elements in dimension 1
        size_type getLength() const { return getLength1(); }

        //! get the first index of dimension 1
        index_type getBegin1() const { return first1_; }
        //! get the first index of dimension 1
        index_type getBegin() const { return first1_; }

        //! get the first index of dimension 1
        index_type getFirst1() const { return first1_; }
        //! get the first index of dimension 1
        index_type getFirst() const { return first1_; }

        //! get the last index of dimension 1
        index_type getLast1() const { return last1_; }
        //! get the last index of dimension 1
        index_type getLast() const { return last1_; }

        //! get the end as in std containers, one past the last valid index of dimension 1
        index_type getEnd1() const { return last1_ + 1; }
        //! get the end as in std containers, one past the last valid index of dimension 1
        index_type getEnd() const { return last1_ + 1; }

        //! determine the amount of memory occupied by the object, usually not needed, therefore computed on the fly
        size_type getSize() const
        {
            size_type size = sizeof(*this)+ getLength1() * sizeof(value_type);
            return size;
        }

        //! determine the number of data array elements
        size_type getNelements() const { return getLength(); }
        //! for compatibility with std::vector, determine the number of data array elements
        size_type size() const { return getLength(); }

        //! deallocate memory
        void deallocateArray()
        {
            if (arr_ != nullptr)
            {
                for (index_type i = getFirst1(); i < getLast1(); ++i)
                {
                    std::allocator_traits<allocator_type>::destroy(allocator_, &arr_[i]);
                }
                arr_ += first1_;
                std::allocator_traits<allocator_type>::deallocate(allocator_, arr_, getLength1());
                arr_ = nullptr;
            }
            first1_ =  0;
            last1_  = -1;
        }

        //! destructor
        ~IrregArray1D()
        {
            // deallocate arrays
            deallocateArray();
        }

        // index operators

        /*! \brief  index operator[x] for array-like usage
            \param[in]   x   index in dimension 1  */
        value_type &operator[](const index_type x)
        {
            return arr_[x];
        }
        /*! \brief  index operator[x] for array-like usage
            \param[in]   x   index in dimension 1  */
        const value_type &operator[](const index_type x) const
        {
            return arr_[x];
        }

        /*! \brief  index operator(x) for functor-like usage
            \param[in]   x   index in dimension 1  */
        value_type &operator()(const index_type x)
        {
            return arr_[x];
        }
        /*! \brief index operator(x) for functor-like usage
            \param[in]   x   index in dimension 1  */
        const value_type &operator()(const index_type x) const
        {
            return arr_[x];
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

        /*! \brief  output stream operator to insert a string representation of the
                    object into an output stream, e.g., cout

            \param[in]   output   ostream in which the content is to be inserted
            \param[in]   vec      the object whose contents are to be inserted   */
        friend std::ostream &operator<<(std::ostream &output, const IrregArray1D &vec)
        {
            // safe the current ostream format for restoring it after output insertion
            std::ios  state(NULL);
            state.copyfmt(output);

            if (vec.arr_ != nullptr)
            {
                for (index_type i = vec.first1_; i <= vec.last1_; i++)
                {
                    if (std::is_floating_point<value_type>::value)
                    {
                        output << std::setw(10) << std::fixed << std::setprecision(2) << vec(i);
                    }
                    else
                    {
                        output << vec(i);
                    }
                    if (i != vec.last1_)
                    {
                        output << " ";
                    }
                }
            }

            // restore the original ostream format
            output.copyfmt(state);

            return output;
        }

        /*! \brief  move assignment operator this = rhs

            \param[in]   rhs   right-hand side of the assignment
         */
        IrregArray1D &operator=(IrregArray1D &&rhs) noexcept
        {
            if (this != &rhs)
            {
                // free any previously allocated memory
                deallocateArray();
                // move the data from rhs
                first1_       = rhs.first1_;
                last1_        = rhs.last1_;
                arr_          = rhs.arr_;
                allocator_    = std::move(rhs.allocator_);
                // set rhs to a valid default state
                rhs.first1_   =  0;
                rhs.last1_    = -1;
                rhs.arr_      = nullptr;
            }
            return *this;
        }

        /*! \brief  assignment operator this = rhs

            \param[in]   rhs   right-hand side of the assignment
         */
        IrregArray1D &operator=(const IrregArray1D &rhs)
        {
            deallocateArray();
            initArray(rhs);
            return *this;
        }

        /*! \brief  assignment operator this = rhs

            \tparam      AForeign   allocator used by the input data structure

            \param[in]   rhs   right-hand side of the assignment
         */
        template<class AForeign = allocator_size_type>
        IrregArray1D &operator=(const IrregArray1D<value_type, AForeign> &rhs)
        {
            deallocateArray();
            initArray(rhs);
            return *this;
        }

        /*! \brief  conversion operator (explicit to avoid unintended implicit conversion), usable via
                    IrregArray1D<TForeign, AForeign> newArr = static_cast<IrregArray1D<TForeign, AForeign> >(IrregArray1D<TCurrent, ACurrent> oldArr);

            \tparam      TForeign   data type stored by the input data structure
            \tparam      AForeign   allocator used by the input data structure
         */
        template<typename TForeign, class AForeign = std::allocator<TForeign> >
        explicit
        operator IrregArray1D<TForeign, AForeign>() const
        {
            IrregArray1D<TForeign, AForeign> result;
            result.initArray(*this);
            return result;
        }

    private:

        /*! \brief  allocate memory for the array (regular shape)

            \param[in]   x1   first index for dimension 1
            \param[in]   x2   last index for dimension 1
         */
        void allocateArray(const index_type x1, const index_type x2)
        {
            deallocateArray(); // prevent memory leaks
            first1_    = x1;
            last1_     = x2;
            if (last1_ < first1_)
            {
                std::string errorMessage("Error in %s::allocateArray(): called with illegal bounds first1 > last1.", typeid(*this).name());
                GMX_THROW(InvalidInputError(errorMessage));
            }
            size_type dim1 = getLength1();
            try
            {
                arr_  = std::allocator_traits<allocator_type>::allocate(allocator_, dim1);
                arr_ -= first1_;
                for (index_type i = getFirst1(); i <= getLast1(); ++i)
                {
                    std::allocator_traits<allocator_type>::construct(allocator_, &arr_[i], value_type());
                }
            }
            catch (const std::bad_alloc &)
            {
                std::fprintf(stderr, "Error in %s::allocateArray(): could not allocate memory for the data array.", typeid(*this).name());
                throw;
            }
        }
};

} // end namespace gmx

#endif
