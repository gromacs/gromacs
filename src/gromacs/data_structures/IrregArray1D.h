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
/*! \file
    \ingroup module_data_structures
    \inpublicapi

    \brief IrregArray1D is a 1-dimensional array
           with custom begin and end index

    \author R. Thomas Ullmann <tullman@gwdg.de>

    \copyright GROMACS license

    \date Feb 2015
 */
#ifndef GMX_IrregArray1D_H
#define GMX_IrregArray1D_H

#include <vector>
#include <iostream>
#include <iomanip>
#include <typeinfo>
#include <stdexcept>

#include "gromacs/utility/exceptions.h"

namespace gmx
{


/*! \class IrregArray1D IrregArray1D.h "gromacs/data_structures/IrregArray1D.h"

    \brief 1-dimensional array with custom begin and end index

    \author R. Thomas Ullmann <tullman@gwdg.de>

    \copyright GROMACS license

    \date Feb 2015

    \ingroup module_data_structures
    \inpublicapi

    \tparam   T           data type to be stored
    \tparam   Allocator   allocator to be used in creating arrays storing T
                          allocators for storing additional book-keeping data
                          are derived from this allocator via rebind if necessary
 */
template<class T = size_t, typename Allocator = std::allocator<T> >
class IrregArray1D
{
    public:
        /*! adopted these typedefs from AnBe's TriangularArray, for consistency with the other data
            structures, can be retrieved from outside to make sure that the intended types are used
            when accessing or manipulating the array */
        //! the data type stored in the array
        typedef T value_type;
        //! the allocator type to allocate value_type
        typedef typename Allocator::template rebind<value_type>::other allocator_type;
        //! the allocator type to allocate value_type*
        typedef typename Allocator::template rebind<value_type>::other::pointer allocator_return_ptr_type;
        //! the (integer) datatype used to represent sizes, e.g., array lengths
        typedef  size_t  size_type;
        //! the (integer) datatype used to represent array indices
        typedef ssize_t index_type;
    private:
        //! first index of dimension 1
        index_type     begin1;
        //! ending index of dimension 1
        index_type     end1;
        //! array storing the the actual data
        value_type    *arr;
        //! the allocator used to allocate arr
        allocator_type allocator;
        //! let other IrregArray1D instantations access private members (for the explicit conversion operator)
        template<class W, typename AW> friend class IrregArray1D;
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
            : begin1(0), end1(0), arr(NULL) {}
        /*! \brief constructor with memory allocation but without initialization
            \param[in]   x1   first index for dimension 1
            \param[in]   x2   last index for dimension 1
         */
        IrregArray1D(const index_type x1, const index_type x2)
            : begin1(0), end1(0), arr(NULL)
        {
            allocateArray(x1, x2);
        }
        /*! \brief constructor with allocation and initialization
            \param[in]   x1   first index for dimension 1
            \param[in]   x2   last index for dimension 1
            \param[in]   ini  initial value for all array elements
         */
        IrregArray1D(const index_type x1, const index_type x2, const value_type &ini)
            : begin1(0), end1(0), arr(NULL)
        {
            initArray(x1, x2, ini);
        }
        /*! \brief constructor with allocation but without initialization
            \param[in]   x1        first index for dimension 1
            \param[in]   x2        last index for dimension 1
            \param[in]   ini_vec   initial values for all array elements
         */
        IrregArray1D(const index_type x1, const index_type x2, const value_type* const ini_vec)
            : begin1(0), end1(0), arr(NULL)
        {
            allocateArray(x1, x2);
            for (index_type i = getBegin1(); i <= getEnd1(); i++)
            {
                arr[i] = ini_vec[i];
            }

        }
        /*! \brief constructor with allocation but without initialization
            \param[in]   x1        first index for dimension 1
            \param[in]   x2        last index for dimension 1
            \param[in]   ini_vec   vector to be copied to a IrregArray1D
         */
        IrregArray1D(const index_type x1, const index_type x2, const std::vector<value_type> &ini_vec)
            : begin1(0), end1(0), arr(NULL)
        {
            allocateArray(x1, x2);
            for (index_type i = x1; i <= x2; i++)
            {
                arr[i] = ini_vec[i];
            }

        }
        /*! \brief copy constructor from a std::vector
            \param[in]   ini_vec   vector to be copied to a IrregArray1D
         */
        IrregArray1D(const std::vector<value_type> &ini_vec)
            : begin1(0), end1(0), arr(NULL)
        {
            const index_type x2 = ini_vec.size() - 1;
            allocateArray(0, x2);
            for (index_type i = 0; i <= x2; i++)
            {
                arr[i] = ini_vec[i];
            }

        }
        /*! \brief  constructor with allocation but without initialization

            \param[in]   x1        first index for dimension 1
            \param[in]   x2        last index for dimension 1
            \param[in]   ini_vec   IrregArray1D to be copied
         */
        IrregArray1D(const index_type x1, const index_type x2, const IrregArray1D<value_type> &ini_vec)
            : begin1(0), end1(0), arr(NULL)
        {
            allocateArray(x1, x2);
            for (index_type i = x1; i <= x2; i++)
            {
                arr[i] = ini_vec[i];
            }

        }
        //! copy constructor
        IrregArray1D(const IrregArray1D<value_type> &ini_vec)
            : begin1(0), end1(0), arr(NULL)
        {
            initArray(ini_vec);
        }
        /*! \brief  initialize IrregArray1D<value_type> from IrregArray1D<T_foreign>,
                    used by the copy constructor and by the explicit conversion operator,
                    T_foreign must be convertible to value_type via static_cast

            \param[in]   ini   the template to be copied
         */
        template <typename T_foreign, typename A_foreign = std::allocator<T_foreign> >
        void initArray(const IrregArray1D<T_foreign, A_foreign> &ini)
        {
            allocateArray(ini.getBegin1(), ini.getEnd1());
            for (index_type i = ini.getBegin1(); i <= ini.getEnd1(); i++)
            {
                arr[i] = static_cast<value_type>(ini(i));
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
                arr[i] = ini;
            }
        }
        /*! \brief  allocate memory for the array (regular shape)

            \param[in]   x1   first index for dimension 1
            \param[in]   x2   last index for dimension 1
         */
        void allocateArray(const index_type x1, const index_type x2)
        {
            deallocateArray(); // prevent memory leaks
            begin1   = x1;
            end1     = x2;
            if (end1 < begin1)
            {
                std::fprintf(stderr, "Error in %s::allocateArray(): called with illegal bounds begin1 > end1.", typeid(*this).name());
                throw InternalError("Trying to create an irregular array with illegal bounds begin1 > end1");
            }
            size_type dim1 = getLength1();
            try
            {
                arr  = allocator.allocate(dim1);
                arr -= begin1;
                for (index_type i = getBegin1(); i <= getEnd1(); ++i)
                {
                    allocator.construct(&arr[i], value_type());
                }
            }
            catch (std::bad_alloc)
            {
                std::fprintf(stderr, "Error in %s::allocateArray(): could not allocate memory for the data array.", typeid(*this).name());
                throw;
            }
            ;
        }
        //! number of elements in dimension 1
        size_type getLength1() const { return (end1 - begin1 + 1); }
        //! number of elements in dimension 1
        size_type getLength() const { return getLength1(); }
        //! get the first index of dimension 1
        index_type getBegin1() const { return begin1; }
        //! get the first index of dimension 1
        index_type getBegin() const { return begin1; }
        //! get the last index of dimension 1
        index_type getEnd1() const { return end1; }
        //! get the last index of dimension 1
        index_type getEnd() const { return end1; }

        //! determine the amount of memory occupied by the object, usually not needed, therefore computed on the fly
        size_type getSize() const
        {
            size_type size = sizeof(*this)+ getLength1() * sizeof(value_type);
            return size;
        }

        //! determine the number of data array elements
        size_type getNelements() const { return getLength(); }

        //! sum all elements
        value_type sum() const
        {
            value_type sum = 0;
            for (index_type i = begin1; i <= end1; i++)
            {
                sum += arr[i];
            }
            return sum;
        }
        //! deallocate memory
        void deallocateArray()
        {
            if (arr != NULL)
            {
                for (index_type i = getBegin1(); i < getEnd1(); ++i)
                {
                    allocator.destroy(&arr[i]);
                }
                arr += begin1;
                allocator.deallocate(arr, getLength1());
                arr = NULL;
            }
        }
        //! destructor
        ~IrregArray1D()
        {
            // deallocate arrays
            deallocateArray();
            begin1 = end1 = 0;
        }

        // index operators

        /*! \brief  index operator[x] for array-like usage
            \param[in]   x   index in dimension 1  */
        value_type &operator[](const index_type x)
        {
            return arr[x];
        }
        /*! \brief  index operator[x] for array-like usage
            \param[in]   x   index in dimension 1  */
        const value_type &operator[](const index_type x) const
        {
            return arr[x];
        }

        /*! \brief  index operator(x) for functor-like usage
            \param[in]   x   index in dimension 1  */
        value_type &operator()(const index_type x)
        {
            return arr[x];
        }
        /*! \brief index operator(x) for functor-like usage
            \param[in]   x   index in dimension 1  */
        const value_type &operator()(const index_type x) const
        {
            return arr[x];
        }

        //! get a pointer to the data array
        const value_type* getArray() const
        {
            return arr;
        }
        //! get a pointer to the data array
        value_type* getArray()
        {
            return arr;

        }
        /*! \brief  output stream operator to insert a string representation of the
                    object into an output stream, e.g., cout

            \param[in]   output   ostream in which the content is to be inserted
            \param[in]   vec      the object whose contents are to be inserted   */
        friend std::ostream &operator<<(std::ostream &output, const IrregArray1D &vec)
        {
            if (vec.arr == NULL)
            {
                return output;
            }
            for (index_type i = vec.begin1; i <= vec.end1; i++)
            {
                output << std::setw(10) << std::fixed << std::setprecision(2) << vec(i) << " ";
            }
            return output;
        }
        /*! \brief  assignment operator this = rhs
            \param[in]   rhs   right-hand side of the equation to be assigned to this object
         */
        IrregArray1D<value_type> &operator=(const IrregArray1D<value_type> &rhs)
        {
            deallocateArray();
            initArray(rhs);
            return *this;
        }
        /*! \brief  conversion operator (explicit to avoid unintended implicit conversion), usable via
                    IrregArray1D<T_foreign> newArr = static_cast<IrregArray1D<T_foreign> >(FlatIrregArray<T> oldArr);
         */
        template<typename T_foreign, typename A_foreign = std::allocator<T_foreign> >
#if __cplusplus < 201103L
        explicit
#endif
        operator IrregArray1D<T_foreign, A_foreign>() const
        {
            IrregArray1D<T_foreign, A_foreign> result;
            result.initArray(*this);
            return result;
        }
};

} // end namespace gmx

#endif
