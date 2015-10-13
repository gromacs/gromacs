/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2015 by the GROMACS development team, led by
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

    \brief IrregArray2D is an irregularly shaped 2-dimensional array

    \author R. Thomas Ullmann <tullman@gwdg.de>

    \copyright GROMACS license

    \date Feb 2015
*/
#ifndef GMX_IrregArray2D_H
#define GMX_IrregArray2D_H

#include "IrregArray1D.h"

namespace gmx
{


/*! \class IrregArray2D IrregArray2D.h "gromacs/data_structures/IrregArray2D.h"

    \brief irregularly shaped 2-dimensional array

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
class IrregArray2D
{
    public:
        /*! adopted these typedefs from AnBe's TriangularArray, for consistency with the other data
            structures, can be retrieved from outside to make sure that the intended types are used
            when accessing or manipulating the array */
        //! the data type stored in the array
        typedef T value_type;
        //! the allocator type to allocate value_type*
        typedef typename Allocator::template rebind<value_type*>::other p_allocator_type;
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
        index_type begin1;
        //! first index of dimension 2
        index_type begin2;
        //! ending index of dimension 1
        index_type end1;
        //! ending index of dimension 2
        index_type end2;
        //! number of array elements in dimension 2 as a function of the index in dimension 1
        IrregArray1D<size_type> irreg;
        //! flag says that dim2(i) = end2(i) - begin2(i) + 1 is not a constant
        bool is_irreg;
        //! array storing the the actual data
        value_type **arr;
        //! the allocator used to allocate arr
        p_allocator_type p_allocator;
        //! the allocator used to allocate arr[i]
        allocator_type     allocator;
        //! let other IrregArray2D instantations access private members (for the explicit conversion operator)
        template<class W, typename AW> friend class IrregArray2D;
        //! let higher-dimensional irregular arrays access private members of their book-keeping arrays
        template<class V, typename AV> friend class IrregArray3D;
        //! let higher-dimensional irregular arrays access private members of their book-keeping arrays
        template<class W, typename AW> friend class IrregArray4D;
    public:
        /*! \brief default constructor without memory allocation, allocateArray has to be called for
                   allocating memory initArray for memory allocation and content initialization
        */
        IrregArray2D()
            : begin1(0), begin2(0), end1(0), end2(0), is_irreg(false), arr(NULL) {}
        /*! \brief   constructor with memory allocation but without initialization
            \param[in]   x1   first index for dimension 1
            \param[in]   x2   last index for dimension 1
            \param[in]   y1   first index for dimension 2
            \param[in]   y2   last index for dimension 2
        */
        IrregArray2D(const index_type x1, const index_type x2, const index_type y1, const index_type y2)
            : begin1(0), begin2(0), end1(0), end2(0), is_irreg(false), arr(NULL)
        {
            allocateArray(x1, x2, y1, y2);
        }
        /*! \brief  constructor for a regularly shaped array with allocation and initialization
            \param[in]   x1   first index for dimension 1
            \param[in]   x2   last index for dimension 1
            \param[in]   y1   first index for dimension 2
            \param[in]   y2   last index for dimension 2
            \param[in]   ini  initial value for all array elements
        */
        IrregArray2D(const index_type x1, const index_type x2, const index_type y1, const index_type y2, const value_type& ini)
            : begin1(0), begin2(0), end1(0), end2(0), is_irreg(false), arr(NULL)
        {
            initArray(x1, x2, y1, y2, ini);
        }
        /*! \brief  constructor with allocation but without initialization
            \param[in]   x1      first index for dimension 1
            \param[in]   x2      last index for dimension 1
            \param[in]   sizes   sizes[i] specifies the number of array elements [j] in arr[i]
        */
        IrregArray2D(const index_type x1, const index_type x2, const size_type * const sizes)
            : begin1(0), begin2(0), end1(0), end2(0), is_irreg(false), arr(NULL)
        {
            allocateArray(x1, x2, IrregArray1D<size_type>(x1, x2, sizes));
        }
        /*! \brief  constructor with memory allocation
            \param[in]   x1      first index for dimension 1
            \param[in]   x2      last index for dimension 1
            \param[in]   sizes   sizes[i] specifies the number of array elements [j] in arr[i]
        */
        IrregArray2D(const index_type x1, const index_type x2, const IrregArray1D<size_type> &sizes)
            : begin1(0), begin2(0), end1(0), end2(0), is_irreg(false), arr(NULL)
        {
            allocateArray(x1, x2, sizes);
        }
        /*! \brief  constructor with memory allocation and content initialization
            \param[in]   x1      first index for dimension 1
            \param[in]   x2      last index for dimension 1
            \param[in]   sizes   sizes[i] specifies the number of array elements [j] in arr[i]
            \param[in]   ini     initial value for all array elements
        */
        IrregArray2D(const index_type x1, const index_type x2, const IrregArray1D<size_type> &sizes, const value_type& ini)
            : begin1(0), begin2(0), end1(0), end2(0), is_irreg(false), arr(NULL)
        {
            allocateArray(x1, x2, sizes);
            for (index_type i = begin1; i <= end1; ++i)
            {
                for (index_type j = begin2; j <= getEnd2(i); ++j)
                {
                    arr[i][j] = ini;
                }
            }
        }
        /*! \brief  allocate memory for an IrregArray2D and initialize its values
            \param[in]   x1       first index for dimension 1 and dimension 2
            \param[in]   x2       last index for dimension 1
            \param[in]   sizes    vector storing the array lengths in dimension 2
            \param[in]   ini      initialization value in second variant
        */
        void initArray(const index_type x1, const index_type x2, const IrregArray1D<size_type> &sizes, const value_type& ini)
        {
            allocateArray(x1, x2, sizes);
            for(index_type i=getBegin1(); i<=getEnd1(); i++)
            {
                for(index_type j=getBegin2(); j <= getEnd2(i); j++)
                {
                    arr[i][j] = ini;
                }
            }
        }
        /*! \brief  allocate memory for an IrregArray2D and initialize its values
            \param[in]   x1       first index for dimension 1 and dimension 2
            \param[in]   x2       last index for dimension 1
            \param[in]   sizes    vector storing the array lengths in dimension 2
            \param[in]   ini      initialization value in second variant
        */
        void initArray(const index_type x1, const index_type x2, const std::vector<size_type> &sizes, const value_type& ini)
        {
            initArray(x1, x2, IrregArray1D<size_type>(x1, x2, sizes), ini);
        }
        /*! \brief  constructor with allocation but without initialization
            \param[in]   x1      first index for dimension 1
            \param[in]   x2      last index for dimension 1
            \param[in]   sizes   sizes[i] specifies the number of array elements arr[i][j] in arr[i]
            \param[in]   ini     initial value for all array elements
        */
        IrregArray2D(const index_type x1, const index_type x2, const index_type * const sizes, const value_type& ini)
            : begin1(0), begin2(0), end1(0), end2(0), is_irreg(false), arr(NULL)
        {
            initArray(x1, x2, sizes, ini);
        }
        /*! \brief  constructor with memory allocation and content initialization
            \param[in]   x1      first index for dimension 1
            \param[in]   x2      last index for dimension 1
            \param[in]   sizes   sizes[i] specifies the number of array elements arr[i][j] in arr[i]
            \param[in]   ini     initial value for all array elements
        */
        IrregArray2D(const index_type x1, const index_type x2, const std::vector<size_type> &sizes, const value_type& ini)
            : begin1(0), begin2(0), end1(0), end2(0), is_irreg(false), arr(NULL)
        {
            initArray(x1, x2, IrregArray1D<size_type>(x1, x2, sizes), ini);
        }
        /*! \brief  constructor with memory allocation but without content initialization
            \param[in]   x1      first index for dimension 1
            \param[in]   x2      last index for dimension 1
            \param[in]   sizes   sizes[i] specifies the number of array elements arr[i][j] in arr[i]
        */
        IrregArray2D(const index_type x1, const index_type x2, const std::vector<size_type> &sizes)
            : begin1(0), begin2(0), end1(0), end2(0), is_irreg(false), arr(NULL)
        {
            initArray(x1, x2, IrregArray1D<size_type>(x1, x2, sizes));
        }
        /*! \brief  copy constructor with allocation and initialization
            \param[in]   ini_mat   the template to be copied
        */
        IrregArray2D(const IrregArray2D<value_type> &ini_mat)
            : begin1(0), begin2(0), end1(0), end2(0), is_irreg(false), arr(NULL)
        {
            initArray(ini_mat);
        }
        /*! \brief  initialize IrregArray2D<value_type> from IrregArray2D<T_foreign>,
                    used by the copy constructor and by the explicit conversion operator,
                    T_foreign must be convertible to value_type via static_cast

            \param[in]    ini    array to be copied
         */
        template <typename T_foreign, typename A_foreign = std::allocator<T_foreign> >
        void initArray(const IrregArray2D<T_foreign, A_foreign> &ini)
        {
            if (ini.is_irreg)
            {
                allocateArray(ini.getBegin1(), ini.getEnd1(), ini.irreg);
            }
            else
            {
                allocateArray(ini.getBegin1(), ini.getEnd1(), ini.getBegin2(), ini.getEnd2());
            }

            for(index_type i = ini.getBegin1(); i <= ini.getEnd1(); i++)
            {
                for(index_type j=ini.getBegin2(); j<=ini.getEnd2(i); j++)
                {
                    arr[i][j] = static_cast<value_type>(ini(i, j));
                }
            }
        }
        /*! \brief  initialize a regular special case of IrregArray2D
            \param[in]   x1    first index for dimension 1
            \param[in]   x2    last index for dimension 1
            \param[in]   y1    first index for dimension 2
            \param[in]   y2    last index for dimension 2
            \param[in]   ini   initialization value
        */
        void initArray(const index_type x1, const index_type x2, const index_type y1, const index_type y2, const value_type& ini)
        {
            allocateArray(x1, x2, y1, y2);
            for(index_type i=x1; i<= x2; i++)
            {
                for(index_type j=y1; j<= y2; j++)
                {
                    arr[i][j] = ini;
                }
            }
        }
        /*! \brief  allocate memory for the array (regular shape)
            \param[in]   x1   first index for dimension 1
            \param[in]   x2   last index for dimension 1
            \param[in]   y1   first index for dimension 2
            \param[in]   y2   last index for dimension 2
        */
        void allocateArray(const index_type x1, const index_type x2, const index_type y1, const index_type y2)
        {
            deallocateArray(); // prevent memory leaks
            is_irreg = false;
            begin1   = x1;
            end1     = x2;
            begin2   = y1;
            end2     = y2;
            index_type dim1 = getLength1();
            index_type dim2 = getLength2();
            if(dim1 < 0)
            {
                std::fprintf(stderr, "Error in %s::allocateArray(): called with illegal bounds begin1 > end1.", typeid(*this).name());
                throw InternalError("Trying to create an irregular array with illegal bounds begin1 > end1");
            }
            if(dim2 < 0)
            {
                std::fprintf(stderr, "Error in %s::allocateArray(): called with illegal bounds begin2 > end2.", typeid(*this).name());
                throw InternalError("Trying to create an irregular array with illegal bounds begin2 > end2");
            }
            try
            {
                arr = p_allocator.allocate(dim1);
                arr -= begin1;
		arr[begin1]  = allocator.allocate(dim1*dim2);
                arr[begin1] -= y1;
                for (index_type i = begin1 + 1; i <= x2; ++i)
                {
                    arr[i] = arr[i-1] + dim2;
                }
            }
            catch(std::bad_alloc)
            {
                std::fprintf(stderr, "Error in %s::allocateArray(): could not allocate memory for the data array.", typeid(*this).name());
                throw;
            };

            //! initialize the array elements
            for (index_type i = getBegin1(); i <= getEnd1(); ++i)
            {
                for (index_type j = getBegin2(); j <= getEnd2(); ++j)
                {
                    allocator.construct(&this->operator()(i, j), value_type());
                }
            }
        }
        /*! \brief  allocate memory for an IrregArray2D
            \param[in]   x1       first index for dimension 1 and dimension 2
            \param[in]   x2       last index for dimension 1
            \param[in]   &sizes   vector storing the array lengths in dimension 2
        */
        void initArray(const index_type x1, const index_type x2, const std::vector<size_type> &sizes)
        {
            initArray(x1, x2, IrregArray1D<size_type>(sizes));
        }
        /*! \brief  allocate memory for an IrregArray2D
            \param[in]   x1       first index for dimension 1 and dimension 2
            \param[in]   x2       last index for dimension 1
            \param[in]   &sizes   vector storing the array lengths in dimension 2
        */
        void initArray(const index_type x1, const index_type x2, const IrregArray1D<size_type> &sizes)
        {
            allocateArray(x1, x2, sizes);
        }
        /*! \brief  allocate memory for the array (irregular shape)
            \param[in]   x1       first index for dimension 1 and dimension 2
            \param[in]   x2       last index for dimension 1
            \param[in]   &sizes   vector storing the array lengths in dimension 2
        */
        void allocateArray(const index_type x1, const index_type x2, const IrregArray1D<size_type> &sizes)
        {
            deallocateArray(); // prevent memory leaks
            is_irreg = true;
            begin1   = x1;
            begin2   = begin1;
            end1     = x2;
            end2     = 0;

            const index_type dim1 = getLength1();
            if(dim1 <= 0){
                std::fprintf(stderr, "Error in %s::allocateArray(): called with illegal bounds begin1 > end1.", typeid(*this).name());
                throw InternalError("Trying to create an irregular array with illegal bounds begin1 > end1");
            }

            irreg.initArray(x1, x2);
            for(index_type i = x1; i <= x2; i++)
            {
                irreg[i] = sizes[i];
            }

            try {
                arr = p_allocator.allocate(dim1);
                arr -= begin1;
                for (index_type i = begin1; i <= end1; ++i){
                    arr[i] = allocator.allocate(sizes[i]);
                    arr[i] -= begin2;
                }
            } catch(std::bad_alloc) {
                std::fprintf(stderr, "Error in %s::allocateArray(): could not allocate memory for the data array.", typeid(*this).name());
                throw;
            };

            //! initialize the array elements
            for (index_type i = getBegin1(); i <= getEnd1(); ++i)
            {
                for (index_type j = getBegin2(); j <= getEnd2(i); ++j)
                {
                    allocator.construct(&this->operator()(i, j), value_type());
                }
            }

        }
        /*! \brief  initialize an IrregArray2D and values in it.
            \param[in]    x1      first index for dimension 1
            \param[in]    x2      last index for dimension 1
            \param[in]    sizes   vector storing the array lengths in dimension 2
        */
        void initArray(const index_type x1, const index_type x2, const size_type *sizes)
        {
            IrregArray1D<size_type> tmpvec(x1, x2, sizes);
            initArray(x1, x2, tmpvec);
        }
        /*! \brief  initialize an IrregArray2D and values in it.
            \param[in]    x1      first index for dimension 1 and dimension 2
            \param[in]    x2      last index for dimension 1
            \param[in]    sizes   vector storing the array lengths in dimension 2
            \param[in]    ini     value used to initialize the array elements
        */
        void initArray(const index_type x1, const index_type x2, const size_type *sizes, const value_type& ini)
        {
            IrregArray1D<size_type> tmpvec(x1, x2, sizes);
            initArray(x1, x2, tmpvec, ini);
        }
        //! number of elements in dimension 1
        size_type getLength1() const { return (end1 - begin1 + 1); }
        //! number of elements in dimension 2 (regular array)
        size_type getLength2() const
        {
            if (!is_irreg) {
                return (end2 - begin2 + 1);
            }
            else
            {
                std::fprintf(stderr, "Error in %s::getLength2(): called without argument for an irregular array.", typeid(*this).name());
                throw InternalError("Extent of a variable dimension in an irregular array requested without specifying the index/indices in the lower dimension(s).");
            }
        }
        /*! \brief  number of elements in dimension 2
            \param[in]    x       index in dimension 1 for which the array length in dimension 2 is requested */
        size_type getLength2(const index_type x) const
        {
            return (is_irreg ? irreg[x] : (end2 - begin2 + 1));
        }
        //! get the first index of dimension 1
        index_type getBegin1() const { return begin1; }
        //! get the first index of dimension 2
        index_type getBegin2() const { return begin2; }
        //! get the last index of dimension 1
        index_type getEnd1() const { return end1; }
        //! get the last index of dimension 2
        index_type getEnd2() const
        {
            if (!is_irreg)
            {
                return end2;
            }
            else
            {
                std::fprintf(stderr, "Error in %s::getEnd2(): called without argument for an irregular array.", typeid(*this).name());
                throw InternalError("Extent of a variable dimension in an irregular array requested without specifying the index/indices in the lower dimension(s).");
            }
        }
        /*! \brief  initialize an IrregArray2D and values in it
            \param[in]    x       index in dimension 1 for which the last index in dimension 2 is requested */
        index_type getEnd2(const index_type x) const
        {
            return (is_irreg ? (irreg[x] + begin2 - 1) : end2);
        }

        //! determine the amount of memory occupied by the object, usually not needed, therefore computed on the fly
        size_type getSize() const
        {
            size_type size  = sizeof(*this);
            size           += irreg.getSize() - sizeof(IrregArray1D<size_type>);
            if (arr != NULL)
            {
                if(is_irreg)
                {
                    for (index_type i = begin1; i <= end1; ++i)
                    {
                        size += getLength2(i) * sizeof(value_type);
                    }
                    size += getLength1() * sizeof(value_type*);
                }
                else
                {
                    size += getLength1() * getLength2() * sizeof(value_type);
                    size += getLength1() * sizeof(value_type*);
                }
            }
            return size;
        }

        //! determine the number of data array elements
        size_type getNelements() const
        {
            // the storage data
            size_type nelements = 0;
            for(index_type i = getBegin1(); i <= getEnd1(); i++)
            {
                nelements += getLength2(i);
            }
            return nelements;
        }

        /*! \brief  sum all elements in a row (elements along the second dimension at fixed index in the first dimension)
            \param[in]   row    index in dimension 1  */
	value_type sumRow(const index_type row) const
        {
            value_type sum = 0;
            for(index_type i = begin2; i <= getEnd2(row); i++)
            {
                sum += this->operator()(row,i);
            }
	    return sum;
        }

        //! deallocate memory
        void deallocateArray()
        {
            if(arr != NULL)
            {

                //! destroy the array elements
                for (index_type i = getBegin1(); i <= getEnd1(); ++i)
                {
                    for (index_type j = getBegin2(); j <= getEnd2(i); ++j)
                    {
                        allocator.destroy(&this->operator()(i, j));
                    }
                }

                if(is_irreg)
                {
                    for (index_type i = begin1 ; i <= end1; ++i)
                    {
                        arr[i] += begin2;
                        allocator.deallocate(arr[i], getLength2(i));
                    }
                    arr += begin1;
                    p_allocator.deallocate(arr, getLength1());
                    arr = NULL;
                    irreg.deallocateArray();
                }
                else
                {
                    arr[begin1] += begin2;
                    allocator.deallocate(arr[begin1], getLength1() * getLength2());
                    arr += begin1;
                    p_allocator.deallocate(arr, getLength1());
                    arr = NULL;
                }
            }
        }
        //! destructor
        ~IrregArray2D()
        {
            // deallocate arrays
            deallocateArray();
            begin1 = begin2 = end1 = end2 = 0;
        }

        // index operators

        /*! \brief  index operator for bare array-like usage
            IrregArray2d<value_type> a;
            value_type  val   = a[i][j]
            FIXME: I did not get
                   value_type *val   = a[i]
                   to work, because then there would
                   be different return types for different
                   instantiations of operator[]
                   recursive templates might be a way, but this variant
                   seems cumbersome, at best, to be married to the irregular shape
        */
        //! doxygen can not handle nested classes
        /// \cond DEV
        /*! \class Proxy

            \brief used to implement index operator *this[i][j]
                   as described by Mike Seymour on stackoverflow

            \param[in]   x   index in dimension 1
            \param[in]   y   index in dimension 2
        */
        class Proxy {
            public:
                Proxy(value_type* array) : _array(array) { }

                value_type& operator[](index_type y)
                {
                    return _array[y];
                }
                const value_type& operator[](index_type y) const
                {
                    return _array[y];
                }
            private:
                value_type* _array;
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
           return Proxy(arr[x]);
        }
        /*! \brief  index operator[x] [y] for array-like usage
            \param[in]   x   index in dimension 1
            \cond
            \param[in]   y   index in dimension 2
            \endcond
        */
        const Proxy operator[](index_type x) const
        {
           return Proxy(arr[x]);
        }

        /*! \brief  index operator(x,y) for functor-like usage
            \param[in]   x   index in dimension 1
            \param[in]   y   index in dimension 2
        */
        const value_type& operator()(const index_type x, const index_type y) const
        {
                return arr[x][y];
        }
        /*! \brief  index operator(x,y) for functor-like usage
            \param[in]   x   index in dimension 1
            \param[in]   y   index in dimension 2
        */
        value_type& operator()(const index_type x, const index_type y)
        {
                return arr[x][y];
        }
        /*! \brief  index operator(x) for functor-like usage
            \param[in]   x   index in dimension 1
        */
        value_type* operator()(const index_type x)
        {
            return arr[x];
        }
        /*! \brief  index operator(x) for functor-like usage
            \param[in]   x   index in dimension 1
        */
        const value_type* operator()(const index_type x) const
        {
            return arr[x];
        }

        //! get a pointer to the data array
        const value_type** getArray() const
        {
            return arr;
        }
        //! get a pointer to the data array
        value_type** getArray()
        {
            return arr;
        }

        /*! \brief  output stream operator
            \param[in]   output   ostream in which the content is to be inserted
            \param[in]   mat      the object whose contents are to be inserted
        */
        friend std::ostream& operator<<(std::ostream &output, const IrregArray2D &mat)
        {
            if(!mat.is_irreg)
            {
                for(index_type i = mat.begin1; i <= mat.end1; i++)
                {
                    for(index_type j = mat.begin2;j <= mat.end2; j++)
                    {
                        output << std::setw(10) << std::fixed << std::setprecision(2) << mat(i, j) << " ";
                    }
                    output << std::endl;
                }
            }
            else
            {
                for(index_type i = mat.begin1; i <= mat.end1; i++)
                {
                    for(index_type j = mat.begin2; j <= mat.getEnd2(i); j++)
                    {
                        output << std::fixed << std::setprecision(2) << mat(i, j) << " ";
                    }
                    output << std::endl;
                }
            }
            return output;
        }
        /*! \brief  assignment operator this = rhs
            \param[in]   rhs   right-hand side of the equation to be assigned to this object
        */
        IrregArray2D<value_type>& operator=(const IrregArray2D<value_type> &rhs)
        {
            deallocateArray();
            initArray(rhs);
            return *this;
        }
        /*! \brief  conversion operator (explicit to avoid unintended implicit conversion), usable via
                    IrregArray2D<T_foreign> newArr = static_cast<IrregArray2D<T_foreign> >(FlatIrregArray<T> oldArr);
        */
        template<typename T_foreign, typename A_foreign = std::allocator<T_foreign> >
#if __cplusplus < 201103L
            explicit
#endif
            operator IrregArray2D<T_foreign, A_foreign>() const
        {
            IrregArray2D<T_foreign, A_foreign> result;
            result.initArray(*this);
            return result;
        }
};

} // end namespace gmx

#endif
