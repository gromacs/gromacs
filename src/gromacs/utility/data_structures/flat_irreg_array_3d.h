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
    \ingroup module_utility
    \inpublicapi

    \brief FlatIrregArray3D irregularly shaped 3-dimensional
          array internally stored in an 1D-array

    \author R. Thomas Ullmann <tullman@gwdg.de>

    \copyright GROMACS license

    \date Feb 2015
 */
#ifndef GMX_UTILITY_FLAT_IRREG_ARRAY_3D_H
#define GMX_UTILITY_FLAT_IRREG_ARRAY_3D_H

#include "gromacs/utility/data_structures/flat_irreg_array_2d.h"

namespace gmx
{

/*! \class FlatIrregArray3D

    \brief irregularly shaped 3 dimensional array
           internally stored in an 1D-array

    this array has nonuniform second and third dimensions,
    that is, instead of a regular array with constant dimensions
    M x N x O, an irregular array with dimensions M x N(M) x O(M, N).

    example application:
    number of ligands of each ligand typebound by a site form
    indexed by site, form, ligand

    \ingroup module_utility
    \inpublicapi

    \tparam   T           data type to be stored
    \tparam   Allocator   allocator to be used in creating arrays storing T
                          allocators for storing additional book-keeping data
                          are derived from this allocator via rebind if necessary
 */
template<class T = size_t, typename Allocator = std::allocator<T> >
class FlatIrregArray3D
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
        index_type              begin1_;
        //! first index of dimension 2
        index_type              begin2_;
        //! first index of dimension 3
        index_type              begin3_;
        //! ending index of dimension 1
        index_type              end1_;
        //! ending index of dimension 2
        index_type              end2_;
        //! ending index of dimension 3
        index_type              end3_;
        //! total number of array elements arr_[:][:][:]
        size_type               nelements_;
        //! number of array elements in dimension 2 in dependece on the index in dimension 1
        IrregArray1D<size_type> irreg_;
        //! row,column-specific ending index of dimension 3
        /*! irreg_mat_[i][j] = x  is the index of the first actual array element arr_[x] of the
            stripe FlatIrregArray3D[i][j][:] in dimension 3 of the of the represented 3D array */
        FlatIrregArray2D<size_type> irreg_mat_;
        //! flag says that dimX(i) = endX(i) - beginX(i) + 1 is not a constant
        bool                        is_irreg_;
        //! array storing the the actual data
        value_type                 *arr_;
        //! the allocator used to allocate arr_
        allocator_type              allocator_;
        //! let other FlatIrregArray3D instantations access private members (for the explicit conversion operator)
        template<class W, typename AW> friend class FlatIrregArray3D;
        //! let higher-dimensional FlatIrregArrayNDs access private members of their book-keeping arrays
        template<class W, typename AW> friend class FlatIrregArray4D;
    public:
        /*! \brief default constructor without memory allocation, allocateArray has to be called for
                   allocating memory initArray for memory allocation and content initialization
         */
        FlatIrregArray3D()
            : begin1_(0), begin2_(0), begin3_(0), end1_(0), end2_(0), end3_(0), nelements_(0), is_irreg_(false), arr_(nullptr)
        {
        }
        //! copy constructor
        FlatIrregArray3D(const FlatIrregArray3D<value_type> &ini_t)
            : begin1_(0), begin2_(0), begin3_(0), end1_(0), end2_(0), end3_(0), nelements_(0), is_irreg_(false), arr_(nullptr)
        {
            initArray(ini_t);
        }
        /*! \brief   constructor for a regularly shaped array

            \param[in]    x1       starting index in dimension 1
            \param[in]    x2       last index in dimension 1
            \param[in]    y1       starting index in dimension 2
            \param[in]    y2       last index in dimension 2
            \param[in]    z1       starting index in dimension 3
            \param[in]    z2       last index in dimension 3
         */
        FlatIrregArray3D(const index_type x1, const index_type x2, const index_type y1, const index_type y2, const index_type z1, const index_type z2)
            : begin1_(0), begin2_(0), begin3_(0), end1_(0), end2_(0), end3_(0), nelements_(0), is_irreg_(false), arr_(nullptr)

        {
            initArray(x1, x2, y1, y2, z1, z2);
        }
        /*! \brief   constructor for a regularly shaped array data with initialization

            \param[in]    x1       starting index in dimension 1
            \param[in]    x2       last index in dimension 1
            \param[in]    y1       starting index in dimension 2
            \param[in]    y2       last index in dimension 2
            \param[in]    z1       starting index in dimension 3
            \param[in]    z2       last index in dimension 3
            \param[in]    ini      initialization value for the data array elements
         */
        FlatIrregArray3D(const index_type x1, const index_type x2, const index_type y1, const index_type y2, const index_type z1, const index_type z2, const value_type &ini)
            : begin1_(0), begin2_(0), begin3_(0), end1_(0), end2_(0), end3_(0), nelements_(0), is_irreg_(false), arr_(nullptr)

        {
            initArray(x1, x2, y1, y2, z1, z2, ini);
        }
        /*! \brief specialized constructor for irregular, symmetrically shaped matrices
                   application: interaction energy matrix site1 form1(site1) site2 form2(site2)

            \param[in]    x1       starting index in dimension 1
            \param[in]    x2       last index in dimension 1
            \param[in]    sizes    number of index values in dimension 4
                                   as function of the indices in dimension 3
                                   (in this special case, independent on the
                                   index in dimensions 1 and 2)
            \param[in]    z1       starting index in dimension 3
            \param[in]    z2       last index in dimension 3
         */
        FlatIrregArray3D(const index_type x1, const index_type x2, const size_type *sizes, const index_type z1, const index_type z2)
            : begin1_(0), begin2_(0), begin3_(0), end1_(0), end2_(0), end3_(0), nelements_(0), is_irreg_(false), arr_(nullptr)

        {
            initArray(x1, x2, sizes, z1, z2);
        }
        /*! \brief specialized constructor for irregular, symmetrically shaped matrices
                   application: interaction energy matrix site1 form1(site1) site2 form2(site2)

            \param[in]    x1       starting index in dimension 1
            \param[in]    x2       last index in dimension 1
            \param[in]    sizes    number of index values in dimension 4
                                   as function of the indices in dimension 3
                                   (in this special case, independent on the
                                   index in dimensions 1 and 2)
            \param[in]    z1       starting index in dimension 3
            \param[in]    z2       last index in dimension 3
            \param[in]    ini      initialization value for the data array elements
         */
        FlatIrregArray3D(const index_type x1, const index_type x2, const size_type *sizes, const index_type z1, const index_type z2, const value_type &ini)
            : begin1_(0), begin2_(0), begin3_(0), end1_(0), end2_(0), end3_(0), nelements_(0), is_irreg_(false), arr_(nullptr)

        {
            initArray(x1, x2, sizes, z1, z2, ini);
        }
        /*! \brief general constructor for irregularly shaped matrices and all special cases

            \param[in]    x1       starting index in dimension 1
            \param[in]    x2       last index in dimension 1
            \param[in]    sizes1   number of elements in dimension 2 as function of the index in dimension 1
            \param[in]    sizes2   number of elements in dimension 3 as function of the indices in the lower dimensions
         */
        FlatIrregArray3D(const index_type x1, const index_type x2, const size_type *sizes1, const FlatIrregArray2D<size_type> &sizes2)
            : begin1_(0), begin2_(0), begin3_(0), end1_(0), end2_(0), end3_(0), nelements_(0), is_irreg_(false), arr_(nullptr)

        {
            initArray(x1, x2, sizes1, sizes2);
        }
        /*! \brief general constructor for irregularly shaped matrices and all special cases

            \param[in]    x1       starting index in dimension 1
            \param[in]    x2       last index in dimension 1
            \param[in]    sizes1   number of elements in dimension 2 as function of the index in dimension 1
            \param[in]    sizes2   number of elements in dimension 3 as function of the indices in the lower dimensions
            \param[in]    ini      initialization value for the data array elements
         */
        FlatIrregArray3D(const index_type x1, const index_type x2, const size_type *sizes1, const FlatIrregArray2D<size_type> &sizes2, const value_type &ini)
            : begin1_(0), begin2_(0), begin3_(0), end1_(0), end2_(0), end3_(0), nelements_(0), is_irreg_(false), arr_(nullptr)

        {
            initArray(x1, x2, sizes1, sizes2, ini);
        }
        /*! \brief general constructor for irregularly shaped matrices and all special cases

            \param[in]    x1       starting index in dimension 1
            \param[in]    x2       last index in dimension 1
            \param[in]    sizes    number of elements in dimension 2 as function of the index in dimension 1
            \param[in]    z1       starting index in dimension 2
            \param[in]    z2       last index in dimension 2
         */
        FlatIrregArray3D(const index_type x1, const index_type x2, const IrregArray1D<size_type> &sizes, const index_type z1, const index_type z2)
            : begin1_(0), begin2_(0), begin3_(0), end1_(0), end2_(0), end3_(0), nelements_(0), is_irreg_(false), arr_(nullptr)

        {
            initArray(x1, x2, sizes, z1, z2);
        }
        /*! \brief general constructor for irregularly shaped matrices and all special cases

            \param[in]    x1       starting index in dimension 1
            \param[in]    x2       last index in dimension 1
            \param[in]    sizes    number of elements in dimension 2 as function of the index in dimension 1
            \param[in]    z1       starting index in dimension 2
            \param[in]    z2       last index in dimension 2
            \param[in]    ini      initialization value for the data array elements
         */
        FlatIrregArray3D(const index_type x1, const index_type x2, const IrregArray1D<size_type> &sizes, const index_type z1, const index_type z2, const value_type &ini)
            : begin1_(0), begin2_(0), begin3_(0), end1_(0), end2_(0), end3_(0), nelements_(0), is_irreg_(false), arr_(nullptr)

        {
            initArray(x1, x2, sizes, z1, z2, ini);
        }
        /*! \brief general constructor for irregularly shaped matrices and all special cases

            \param[in]    x1       starting index in dimension 1
            \param[in]    x2       last index in dimension 1
            \param[in]    sizes1   number of elements in dimension 2 as function of the index in dimension 1
            \param[in]    sizes2   number of elements in dimension 3 as function of the indices in the lower dimensions
         */
        FlatIrregArray3D(const index_type x1, const index_type x2, const IrregArray1D<size_type> &sizes1, const FlatIrregArray2D<size_type> &sizes2)
            : begin1_(0), begin2_(0), begin3_(0), end1_(0), end2_(0), end3_(0), nelements_(0), is_irreg_(false), arr_(nullptr)

        {
            initArray(x1, x2, sizes1, sizes2);
        }
        /*! \brief general constructor for irregularly shaped matrices and all special cases

            \param[in]    x1       starting index in dimension 1
            \param[in]    x2       last index in dimension 1
            \param[in]    sizes1   number of elements in dimension 2 as function of the index in dimension 1
            \param[in]    sizes2   number of elements in dimension 3 as function of the indices in the lower dimensions
            \param[in]    ini      initialization value for the data array elements
         */
        FlatIrregArray3D(const index_type x1, const index_type x2, const IrregArray1D<size_type> &sizes1, const FlatIrregArray2D<size_type> &sizes2, const value_type &ini)
            : begin1_(0), begin2_(0), begin3_(0), end1_(0), end2_(0), end3_(0), nelements_(0), is_irreg_(false), arr_(nullptr)

        {
            initArray(x1, x2, sizes1, sizes2, ini);
        }
        /*! \brief general constructor for irregularly shaped matrices and all special cases

            \param[in]    x1       starting index in dimension 1
            \param[in]    x2       last index in dimension 1
            \param[in]    y1       starting index in dimension 2
            \param[in]    y2       last index in dimension 2
            \param[in]    sizes    number of elements in dimension 3 as function of the index in dimension 2
         */
        FlatIrregArray3D(const index_type x1, const index_type x2, const index_type y1, const index_type y2, const size_type *sizes)
            : begin1_(0), begin2_(0), begin3_(0), end1_(0), end2_(0), end3_(0), nelements_(0), is_irreg_(false), arr_(nullptr)

        {
            initArray(x1, x2, y1, y2, sizes);
        }
        /*! \brief general constructor for irregularly shaped matrices and all special cases

            \param[in]    x1       starting index in dimension 1
            \param[in]    x2       last index in dimension 1
            \param[in]    y1       starting index in dimension 2
            \param[in]    y2       last index in dimension 2
            \param[in]    sizes    number of elements in dimension 3 as function of the index in dimension 2
            \param[in]    ini      initialization value for the data array elements
         */
        FlatIrregArray3D(const index_type x1, index_type x2, const index_type y1, const index_type y2, const size_type *sizes, const value_type &ini)
            : begin1_(0), begin2_(0), begin3_(0), end1_(0), end2_(0), end3_(0), nelements_(0), is_irreg_(false), arr_(nullptr)

        {
            initArray(x1, x2, y1, y2, sizes, ini);
        }
        /*! \brief  initialize FlatIrregArray3D<value_type> from FlatIrregArray3D<T_foreign>,
                    used by the copy constructor and by the explicit conversion operator,
                    T_foreign must be convertible to value_type via static_cast

            \param[in]   ini   the template to be copied
         */
        template <typename T_foreign, typename A_foreign = std::allocator<T_foreign> >
        void initArray(const FlatIrregArray3D<T_foreign, A_foreign> &ini)
        {
            irreg_     = ini.irreg_;
            irreg_mat_ = ini.irreg_mat_;
            is_irreg_  = ini.is_irreg_;
            begin1_    = ini.begin1_;
            end1_      = ini.end1_;
            begin2_    = ini.begin2_;
            end2_      = ini.end2_;
            begin3_    = ini.begin3_;
            end3_      = ini.end3_;
            nelements_ = ini.nelements_;
            allocateArray(nelements_);
            for (size_type i = 0; i < nelements_; ++i)
            {
                arr_[i] = static_cast<value_type>(ini.arr_[i]);
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
            is_irreg_ = false;
            begin1_   = x1;
            begin2_   = y1;
            begin3_   = z1;
            end1_     = x2;
            end2_     = y2;
            end3_     = z2;
            allocateArray(x1, x2, y1, y2, z1, z2);
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
            is_irreg_ = false;
            initArray(x1, x2, y1, y2, z1, z2);
            for (size_type i = 0; i < nelements_; ++i)
            {
                arr_[i] = ini;
            }
        }
        /*! \brief   constructor helper function, assign the book keeping data and allocate memory for
                     the book keeping and data arrays

            \param[in]    x1      first index in dimension 1
            \param[in]    x2      last index in dimension 1
            \param[in]    sizes1  number of array elements in dimension 2 as a
                                  function of the index in dimension 1
            \param[in]    sizes2  number of array elements in dimension 3 as a
                                  function of the indices in the lower dimensions
         */
        void initArray(const index_type x1, const index_type x2, const size_type *sizes1, const FlatIrregArray2D<size_type> &sizes2)
        {
            is_irreg_ = true;
            IrregArray1D<size_type> tmp_vec(x1, x2);
            for (index_type i = x1; i <= x2; i++)
            {
                tmp_vec(i) = sizes1[i];
            }
            initArray(x1, x2, tmp_vec, sizes2);
        }
        /*! \brief   constructor helper function, assign the book keeping data and allocate memory for
                     the book keeping and data arrays

            \param[in]    x1      first index in dimension 1
            \param[in]    x2      last index in dimension 1
            \param[in]    sizes1  number of array elements in dimension 2 as a
                                  function of the index in dimension 1
            \param[in]    sizes2  number of array elements in dimension 3 as a
                                  function of the indices in the lower dimensions
         */
        void initArray(const index_type x1, const index_type x2, const IrregArray1D<size_type> &sizes1, const FlatIrregArray2D<size_type> &sizes2)
        {
            is_irreg_ = true;
            const index_type dim1 = x2 - x1 + 1;
            begin1_ = x1;
            begin2_ = x1;
            begin3_ = x1;
            end1_   = x2;
            end2_   = 0;
            end3_   = 0;
            allocateArray(getBegin1(), getBegin2(), getBegin3(), dim1, sizes1, sizes2);
        }
        /*! \brief   constructor helper function, assign the book keeping data and allocate memory for
                     the book keeping and data arrays, assign an initial value to each array element

            \param[in]    x1      first index in dimension 1
            \param[in]    x2      last index in dimension 1
            \param[in]    sizes1  number of array elements in dimension 2 as a
                                  function of the index in dimension 1
            \param[in]    sizes2  number of array elements in dimension 3 as a
                                  function of the indices in the lower dimensions
            \param[in]    ini     initialization value for the data array elements
         */
        void initArray(const index_type x1, const index_type x2, const IrregArray1D<size_type> &sizes1, const FlatIrregArray2D<size_type> &sizes2, const value_type &ini)
        {
            initArray(x1, x2, sizes1, sizes2);
            for (size_type i = 0; i < nelements_; ++i)
            {
                arr_[i] = ini;
            }
        }
        /*! \brief   constructor helper function, assign the book keeping data and allocate memory for
                     the book keeping and data arrays, assign an initial value to each array element

            \param[in]    x1      first index in dimension 1
            \param[in]    x2      last index in dimension 1
            \param[in]    sizes2  number of array elements in dimension 3 as a
                                  function of the indices in the lower dimensions,
                                  and in dimension2 accordingly from its irreg member
            \param[in]    ini     initialization value for the data array elements
         */
        void initArray(const index_type x1, const index_type x2, const FlatIrregArray2D<size_type> &sizes2, const value_type &ini)
        {
            initArray(x1, x2, sizes2.irreg_, sizes2, ini);
        }
        /*! \brief   constructor helper function, assign the book keeping data and allocate memory for
                     the book keeping and data arrays

            \param[in]    x1      first index in dimension 1
            \param[in]    x2      last index in dimension 1
            \param[in]    y1      first index in dimension 2
            \param[in]    y2      last index in dimension 2
            \param[in]    sizes   number of array elements in dimension 3 as a
                                  function of the indices in the lower dimensions
         */
        void initArray(const index_type x1, const index_type x2, const index_type y1, const index_type y2, const size_type *sizes)
        {
            IrregArray1D<size_type> tmpvec(y1, y2, (size_type)0);
            for (index_type i = y1; i <= y2; i++)
            {
                tmpvec(i) = sizes[i];
            }
            initArray(x1, x2, y1, y2, tmpvec);
        }
        /*! \brief   constructor helper function assign the book keeping data and allocate memory for
                     the book keeping and data arrays, assign an initial value to each array element

            \param[in]    x1      first index in dimension 1
            \param[in]    x2      last index in dimension 1
            \param[in]    y1      first index in dimension 2
            \param[in]    y2      last index in dimension 2
            \param[in]    sizes   number of array elements in dimension 3 as a
                                  function of the indices in the lower dimensions
            \param[in]    ini     initialization value for the data array elements
         */
        void initArray(const index_type x1, const index_type x2, const index_type y1, const index_type y2, const size_type *sizes, const value_type &ini)
        {
            IrregArray1D<size_type> tmpvec(y1, y2, (size_type)0);
            for (index_type i = y1; i <= y2; i++)
            {
                tmpvec(i) = sizes[i];
            }
            initArray(x1, x2, y1, y2, tmpvec, ini);
        }
        /*! \brief   constructor helper function assign the book keeping data and allocate memory for
                     the book keeping and data arrays, assign an initial value to each array element

            \param[in]    x1      first index in dimension 1
            \param[in]    x2      last index in dimension 1
            \param[in]    sizes   number of array elements in dimension 2 as a
                                  function of the index in dimension 1
            \param[in]    z1      first index in dimension 3
            \param[in]    z2      last index in dimension 3
         */
        void initArray(const index_type x1, const index_type x2, const size_type *sizes, const index_type z1, const index_type z2)
        {
            IrregArray1D<size_type> tmpvec(x1, x2, (size_type)0);
            for (index_type i = x1; i <= x2; i++)
            {
                tmpvec(i) = sizes[i];
            }
            initArray(x1, x2, tmpvec, z1, z2);
        }
        /*! \brief   constructor helper function assign the book keeping data and allocate memory for
                     the book keeping and data arrays, assign an initial value to each array element

            \param[in]    x1      first index in dimension 1
            \param[in]    x2      last index in dimension 1
            \param[in]    sizes   number of array elements in dimension 2 as a
                                  function of the index in dimension 1
            \param[in]    z1      first index in dimension 3
            \param[in]    z2      last index in dimension 3
            \param[in]    ini     initialization value for the data array elements
         */
        void initArray(const index_type x1, const index_type x2, const size_type *sizes, const index_type z1, const index_type z2, const value_type &ini)
        {
            IrregArray1D<size_type> tmpvec(x1, x2, (size_type)0);
            for (index_type i = x1; i <= x2; i++)
            {
                tmpvec(i) = sizes[i];
            }
            initArray(x1, x2, tmpvec, z1, z2);
            for (size_type i = 0; i < nelements_; ++i)
            {
                arr_[i] = ini;
            }
        }
        /*! \brief   constructor helper function assign the book keeping data and allocate memory for
                     the book keeping and data arrays, assign an initial value to each array element

            \param[in]    x1      first index in dimension 1
            \param[in]    x2      last index in dimension 1
            \param[in]    sizes   number of array elements in dimension 2 as a
                                  function of the index in dimension 1
            \param[in]    z1      first index in dimension 3
            \param[in]    z2      last index in dimension 3
            \param[in]    ini     initialization value for the data array elements
         */
        void initArray(const index_type x1, const index_type x2, const IrregArray1D<size_type> &sizes, const index_type z1, const index_type z2, const value_type &ini)
        {
            initArray(x1, x2, sizes, z1, z2);
            for (size_type i = 0; i < nelements_; ++i)
            {
                arr_[i] = ini;
            }
        }
        /*! \brief   constructor helper function assign the book keeping data and allocate memory for
                     the book keeping and data arrays, assign an initial value to each array element

            \param[in]    x1      first index in dimension 1
            \param[in]    x2      last index in dimension 1
            \param[in]    y1      first index in dimension 2
            \param[in]    y2      last index in dimension 2
            \param[in]    sizes   number of array elements in dimension 2 as a
                                  function of the indices in the lower dimensions
            \param[in]    ini     initialization value for the data array elements
         */
        void initArray(const index_type x1, const index_type x2, const index_type y1, const index_type y2, const IrregArray1D<size_type> &sizes, const value_type &ini)
        {
            initArray(x1, x2, y1, y2, sizes);
            for (size_type i = 0; i < nelements_; ++i)
            {
                arr_[i] = ini;
            }
        }
        /*! \brief   assign the book keeping data and allocate memory for the
                     book keeping and data arrays

            \param[in]    x1      first index in dimension 1
            \param[in]    x2      last index in dimension 1
            \param[in]    y1      first index in dimension 2
            \param[in]    y2      last index in dimension 2
            \param[in]    sizes   number of array elements in dimension 2 as a
                                  function of the indices in dimension 2
         */
        void initArray(const index_type x1, const index_type x2, const index_type y1, const index_type y2, const IrregArray1D<size_type> &sizes)
        {
            is_irreg_ = true;
            const index_type dim1 = x2 - x1 + 1;
            const index_type dim2 = y2 - y1 + 1;
            begin1_ = x1;
            begin2_ = y1;
            begin3_ = begin1_;
            end1_   = x2;
            end2_   = y2;
            end3_   = 0;
            IrregArray1D<size_type>     tmp_irreg(begin1_, end1_, dim2);
            FlatIrregArray2D<size_type> tmp_irreg_mat(begin1_, end1_, begin2_, end2_, (size_type)0);
            for (index_type i = begin1_; i <= end1_; ++i)
            {
                for (index_type j = begin2_; j <= end2_; ++j)
                {
                    tmp_irreg_mat(i, j) = sizes(j);
                }
            }
            allocateArray(begin1_, begin2_, begin3_, dim1, tmp_irreg, tmp_irreg_mat);
        }
        /*! \brief   assign the book keeping data and allocate memory for the
                     book keeping and data arrays

            \param[in]    x1      first index in dimension 1
            \param[in]    x2      last index in dimension 1
            \param[in]    sizes   number of array elements in dimension 2
                                  as a function of the index in dimension 1
            \param[in]    z1      first index in dimension 3
            \param[in]    z2      last index in dimension 3
         */
        void initArray(index_type x1, index_type x2, const IrregArray1D<size_type> &sizes, index_type z1, index_type z2)
        {
            is_irreg_ = true;
            const index_type dim1 = x2 - x1 + 1;
            const index_type dim3 = z2 - z1 + 1;
            begin1_ = x1;
            begin2_ = begin1_;
            begin3_ = z1;
            end1_   = x2;
            end2_   = 0;
            end3_   = z2;
            FlatIrregArray2D<size_type> tmp_irreg_mat(begin1_, end1_, begin2_, end2_, (size_type)0);
            for (index_type i = begin1_; i <= end1_; ++i)
            {
                for (index_type j = begin2_; j <= end2_; ++j)
                {
                    tmp_irreg_mat(i, j) = dim3;
                }
            }
            allocateArray(begin1_, begin2_, begin3_, dim1, sizes, tmp_irreg_mat);
        }
        //! number of elements in dimension 1
        size_type getLength1() const { return (end1_ - begin1_ + 1); }
        //! number of elements in dimension 2 (regular array)
        size_type getLength2() const
        {
            if (!is_irreg_)
            {
                return (end2_ - begin2_ + 1);
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
            return (is_irreg_ ? irreg_[x] : (end2_ - begin2_ + 1));
        }
        //! number of elements in dimension 3 (regular array)
        size_type getLength3() const
        {
            if (!is_irreg_)
            {
                return (end3_ - begin3_ + 1);
            }
            else
            {
                std::fprintf(stderr, "Error in %s::getLength3(): called without argument for an irregular array.", typeid(*this).name());
                throw InternalError("Extent of a variable dimension in an irregular array requested without specifying the index/indices in the lower dimension(s).");
            }
        }
        /*! \brief number of elements in dimension 3 at lower dimension indices x, y
            \param[in]    x       index in dimension 1
            \param[in]    y       index in dimension 2  */
        size_type getLength3(const index_type x, const index_type y) const
        {
            if (is_irreg_)
            {
                if (x == end1_ && y == getEnd2(end1_))
                {
                    return nelements_ - irreg_mat_(x, y);
                }
                else
                {
                    // get first element of next stripe in the underlying data array of irreg_mat_ to avoid further conditionals
                    const size_type        i       = irreg_mat_.getArrayIndex(x, y);
                    const size_type* const a       = irreg_mat_.getArray();
                    return a[i + 1] - a[i];
                }
            }
            else
            {
                return end3_ - begin3_ + 1;
            }
        }
        //! get the first index of dimension 1
        index_type getBegin1() const { return begin1_; }
        //! get the first index of dimension 2
        index_type getBegin2() const { return begin2_; }
        //! get the first index of dimension 3
        index_type getBegin3() const { return begin3_; }
        //! get the last index of dimension 1
        index_type getEnd1() const { return end1_; }
        //! get the last index of dimension 2
        index_type getEnd2() const
        {
            if (!is_irreg_)
            {
                return end2_;
            }
            else
            {
                std::fprintf(stderr, "Error in %s::getEnd2(): called without argument for an irregular array.", typeid(*this).name());
                throw InternalError("Extent of a variable dimension in an irregular array requested without specifying the index/indices in the lower dimension(s).");
            }
        }
        //! get the last index of dimension 2 at lower dimension index x
        //! \param[in]    x       index in dimension 1 for which the last index in dimension 2 is requested
        index_type getEnd2(const index_type x) const
        {
            return (is_irreg_ ? (irreg_[x] + begin2_ - 1) : end2_);
        }
        //! get the last index of dimension 3
        index_type getEnd3() const
        {
            if (!is_irreg_)
            {
                return end3_;
            }
            else
            {
                std::fprintf(stderr, "Error in %s::getEnd3(): called without arguments for an irregular array.", typeid(*this).name());
                throw InternalError("Extent of a variable dimension in an irregular array requested without specifying the index/indices in the lower dimension(s).");
            }
        }
        /*! \brief get the last index of dimension 3 for lower dimension indices x, y
            \param[in]    x       index in dimension 1
            \param[in]    y       index in dimension 2                       */
        index_type getEnd3(const index_type x, const index_type y) const
        {
            if (is_irreg_)
            {
                if (x == end1_ && y == getEnd2(end1_))
                {
                    return nelements_ - 1 - irreg_mat_(x, y) + begin3_;
                }
                else
                {
                    // get first element of next stripe in the underlying data array of irreg_mat_ to avoid further conditionals
                    const size_type        i       = irreg_mat_.getArrayIndex(x, y);
                    const size_type* const a       = irreg_mat_.getArray();
                    return a[i + 1] - 1 - a[i] + begin3_;
                }
            }
            else
            {
                return end3_;
            }
        }

        /*! \brief determine the index of the array element FlatIrregArray2D(i, j) in the underlying 1D array arr
            this function needs to be public to be compatible to C++ standards < C++11, where nested classes
            can not access private members of the enclosing class. A friend declaration can also not be used
            because the friends are defined to be entities that are not members.

            length_type instead of index_type used intentionally here to use the maximum possible length

            \param[in]   x    index in dimension 1
            \param[in]   y    index in dimension 2
            \param[in]   z    index in dimension 3
         */
        size_type getArrayIndex(const index_type x, const index_type y, const index_type z) const
        {
            if (is_irreg_)
            {
                return irreg_mat_(x, y) + (z - begin3_);
            }
            else
            {
                return (getLength2() * getLength3() * (x - begin1_))
                       + (getLength3() * (y - begin2_))
                       + (z - begin3_);
            }
        }

        //! returns the amount of memory occupied by the object
        size_type getSize() const
        {
            // the book-keeping data
            size_type size = sizeof(*this);
            size += irreg_.getSize()     - sizeof(IrregArray1D<size_type>);
            size += irreg_mat_.getSize() - sizeof(FlatIrregArray2D<size_type>);
            // the storage data
            size += getLength1() * sizeof(value_type**);
            for (index_type i = getBegin1(); i <= getEnd1(); i++)
            {
                size += getLength2(i) * sizeof(value_type*);
                for (index_type j = getBegin2(); j <= getEnd2(i); j++)
                {
                    size += getLength3(i, j) * sizeof(value_type);
                }
            }
            return size;
        }

        //! determine the number of data array elements
        size_type getNelements() const { return nelements_; }

        /*! \brief allocate memory for a regular special case of the array

            \param[in]    x1       starting index in dimension 1
            \param[in]    x2       last index in dimension 1
            \param[in]    y1       starting index in dimension 2
            \param[in]    y2       last index in dimension 2
            \param[in]    z1       starting index in dimension 3
            \param[in]    z2       last index in dimension 3              */
        void allocateArray(const index_type x1, const index_type x2, const index_type y1, const index_type y2, const index_type z1, const index_type z2)
        {
            is_irreg_ = false;
            begin1_   = x1;
            end1_     = x2;
            begin2_   = y1;
            end2_     = y2;
            begin3_   = z1;
            end3_     = z2;

            if (begin1_ > end1_)
            {
                std::fprintf(stderr, "Error in %s::allocateArray(): called with illegal bounds begin1 >= end1.", typeid(*this).name());
                throw InternalError("Trying to create an irregular array with illegal bounds begin1 >= end1");
            }
            if (begin2_ > end2_)
            {
                std::fprintf(stderr, "Error in %s::allocateArray(): called with illegal bounds begin2 >= end2.", typeid(*this).name());
                throw InternalError("Trying to create an irregular array with illegal bounds begin2 >= end2");
            }
            if (begin3_ > end3_)
            {
                std::fprintf(stderr, "Error in %s::allocateArray(): called with illegal bounds begin3 >= end3.", typeid(*this).name());
                throw InternalError("Trying to create an irregular array with illegal bounds begin3 >= end3");
            }

            const size_type dim1 = getLength1();
            const size_type dim2 = getLength2();
            const size_type dim3 = getLength3();
            nelements_ = dim1 * dim2 * dim3;
            allocateArray(nelements_);
        }
        /*! \brief allocate memory for the array

            \param[in]    s1       starting index in dimension 1
            \param[in]    s2       starting index in dimension 2
            \param[in]    s3       starting index in dimension 3
            \param[in]    d1       number of elements in dimension 1
            \param[in]    irr      number of elements in dimension 2 as function of the index in dimension 1
            \param[in]    irr_mat  number of elements in dimension 3 as function of the indices in the lower dimensions
         */
        void allocateArray(const index_type s1, const index_type s2, const index_type s3, const index_type d1, const IrregArray1D<size_type> &irr, const FlatIrregArray2D<size_type> &irr_mat)
        {
            is_irreg_ = true;
            begin1_   = s1;
            begin2_   = s2;
            begin3_   = s3;
            end1_     = s1 + d1 - 1;
            end2_     = 0;
            end3_     = 0;

            const size_type dim1 = d1;
            // zero length seems to be OK with the standard allocator, but size_type may be signed
            if (static_cast<long int>(dim1) < 0)
            {
                std::fprintf(stderr, "Error in %s::allocateArray(): called with illegal bounds begin1 >= end1.", typeid(*this).name());
                throw InternalError("Trying to create an irregular array with illegal bounds begin1 >= end1");
            }

            irreg_.initArray(begin1_, end1_);
            for (index_type i = begin1_; i <= end1_; i++)
            {
                if (static_cast<long int>(irr(i)) < 0)
                {
                    std::fprintf(stderr, "Error in %s::allocateArray(): called with illegal length %li <= 0 for element %li\n", typeid(*this).name(), static_cast<long int>(irr(i)), i);
                    throw InternalError("Trying to create an irregular array with illegal bounds begin2 >= end2");
                }
                irreg_(i) = irr(i);
            }

            irreg_mat_.initArray(begin1_, end1_, irreg_);
            size_type last_index = 0;
            for (index_type i = begin1_; i <= end1_; i++)
            {
                for (index_type j = begin2_; j <= begin2_ + static_cast<index_type>(irreg_(i)) - 1; j++)
                {
                    // zero length seems to be OK with the standard allocator, but size_type may be signed
                    if (static_cast<long int>(irr_mat(i, j)) < 0)
                    {
                        std::fprintf(stderr, "Error in %s::allocateArray(): called with illegal length %li < 0 for element %li, %li\n", typeid(*this).name(), static_cast<long int>(irr_mat(i, j)), static_cast<long int>(i), static_cast<long int>(j));
                        throw InternalError("Trying to create an irregular array with illegal bounds begin2 > end2");
                    }
                    if (i == begin1_ && j == begin2_)
                    {
                        irreg_mat_(i, j) = 0;
                    }
                    else
                    {
                        irreg_mat_(i, j) = irr_mat(i, j) + last_index;
                        last_index       = irreg_mat_(i, j);
                    }
                }
            }
            // the number of elements is given by the starting index of the last array stripe
            // FlatIrregArray2D[i][:] in dimension 2 within arr + the length of this stripe + 1
            const size_type n = last_index + irr_mat(end1_, begin2_ + irreg_[end1_] - 1);

            allocateArray(n);
        }
        /*! \brief  allocate memory for the array
            \param[in]   n  total number of array elements
         */
        void allocateArray(const size_type n)
        {
            deallocateArray(); // prevent memory leaks
            nelements_ = n;
            try
            {
                arr_ = allocator_.allocate(nelements_);
            }
            catch (std::bad_alloc)
            {
                std::fprintf(stderr, "Error in %s::allocateArray(): could not allocate memory for the data array.", typeid(*this).name());
                throw;
            }

            // initialize the array elements
            for (size_type i = 0; i < nelements_; ++i)
            {
                allocator_.construct(&arr_[i], value_type());
            }
        }
        //! deallocate the memory of all constituent arrays
        void deallocateArray()
        {
            if (arr_ != nullptr)
            {

                // destroy the array elements
                for (size_type i = 0; i < nelements_; ++i)
                {
                    allocator_.destroy(&arr_[i]);
                }

                allocator_.deallocate(arr_, nelements_);
                arr_       = nullptr;
                nelements_ = 0;
            }
        }

        //! destructor
        ~FlatIrregArray3D()
        {
            // deallocate arrays
            deallocateArray();

            // deallocate book-keeping arrays
            if (is_irreg_)
            {
                irreg_mat_.deallocateArray();
                irreg_.deallocateArray();
                is_irreg_ = false;
            }
        }

        // index operators

        /*! index operator for bare array-like usage a[i][j][k] */

        //! doxygen can not handle nested classes
        /// \cond DEV
        /*! \class Proxy

            \brief helper functor used to implement index operator *this[i][j][k]
         */
        class Proxy
        {
            private:
                FlatIrregArray3D<value_type>* optr_;
                index_type                    x_;
                class Proxy2
                {
                    private:
                        FlatIrregArray3D<value_type>* optr_;
                        index_type                    x_;
                        index_type                    y_;
                    public:
                        Proxy2(FlatIrregArray3D<value_type>* optr, index_type x, index_type y)
                            : optr_(optr), x_(x), y_(y)  { }

                        value_type &operator[](index_type z)
                        {
                            return optr_->operator()(x_, y_, z);
                        }
                        const value_type &operator[](index_type z) const
                        {
                            return optr_->operator()(x_, y_, z);
                        }
                };
            public:
                Proxy(FlatIrregArray3D<value_type>* optr, index_type x) : optr_(optr), x_(x) { }

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
        /*! \brief   index operator[x][y][z] for array-like usage realized with chained functors "Proxy*"
            \param[in]   x   index in dimension 1
            \cond
            \param[in]   y   index in dimension 2
            \param[in]   z   index in dimension 3
            \endcond
         */
        Proxy operator[](index_type x)
        {
            return Proxy(this, x);
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
            return Proxy(this, x);
        }

        /*! \brief   index operator(x,y,z) for functor-like usage
            \param[in]   x   index in dimension 1
            \param[in]   y   index in dimension 2
            \param[in]   z   index in dimension 3
         */
        value_type &operator()(const index_type x, const index_type y, const index_type z)
        {
            return arr_[getArrayIndex(x, y, z)];
        }
        /*! \brief   index operator(x,y,z) for functor-like usage
            \param[in]   x   index in dimension 1
            \param[in]   y   index in dimension 2
            \param[in]   z   index in dimension 3
         */
        value_type &operator()(const index_type x, const index_type y, const index_type z) const
        {
            return arr_[getArrayIndex(x, y, z)];
        }

        //! get a pointer to the data array
        const value_type* getArray() const
        {
            return arr_;
        }
        //! get a pointer to the data array
        value_type* getArray()
        {
            return arr_;
        }

        /*! \brief   assign data from a vector to an array stripe in dimension 3
            \param[in]   x      index in dimension 1
            \param[in]   y      index in dimension 2
            \param[in]   data   data array to be assigned to the array elements
                                arr_[x][x][:]
         */
        void assign_vector(const index_type x, const index_type y, const value_type *data)
        {
            for (index_type z = getBegin3(); z <= getEnd3(x, y); z++)
            {
                this->operator()(x, y, z) = data[z];
            }
        }
        /*! \brief   stream operator for convenient printing of the array to an ostream
            \param[in]   output   output stream for the array contents
            \param[in]   ten      the array to be printed
         */
        friend std::ostream &operator<<(std::ostream &output, const FlatIrregArray3D &ten)
        {
            for (index_type i = ten.begin1_; i <= ten.end1_; i++)
            {
                output << "------------------------------" << std::endl << i << std::endl;
                for (index_type j = ten.begin2_; j <= ten.getEnd2(i); j++)
                {
                    for (index_type k = ten.begin3_; k <= ten.getEnd3(i, j); k++)
                    {
                        output << std::setw(10) << std::fixed << std::setprecision(2) << ten(i, j, k) << " ";
                    }
                    output << std::endl;
                }
                output << std::endl;
            }
            return output;
        }
        /*! \brief   assignment operator *this = rhs
            \param[in]   rhs   right-hand side object of the assignment
         */
        FlatIrregArray3D<value_type> &operator=(const FlatIrregArray3D<value_type> &rhs)
        {
            deallocateArray();
            initArray(rhs);
            return *this;
        }
        /*! \brief  conversion operator (explicit to avoid unintended implicit conversion), usable via
                    FlatIrregArray3D<T_foreign> newArr = static_cast<FlatIrregArray3D<T_foreign> >(FlatIrregArray<T> oldArr);
         */
        template<typename T_foreign, typename A_foreign = std::allocator<T_foreign> >
        explicit
        operator FlatIrregArray3D<T_foreign, A_foreign>() const
        {
            FlatIrregArray3D<T_foreign, A_foreign> result;
            result.initArray(*this);
            return result;
        }
};

} // end namespace gmx

#endif
