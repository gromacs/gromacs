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

    \brief IrregArray4D is an irregularly shaped 4-dimensional array

    \author R. Thomas Ullmann <tullman@gwdg.de>

    \copyright GROMACS license

    \date Feb 2015
*/
#ifndef GMX_IrregArray4D_H
#define GMX_IrregArray4D_H

#include "IrregArray3D.h"

namespace gmx
{


/*! \class IrregArray4D

    \brief irregularly shaped 4-dimensional array

     this array has nonuniform second and third dimensions,
     that is, instead of a regular array with constant dimensions
     M x N x O x P, an irregular array with dimensions
     M x N(M) x O(M, N) x P(M, N, O).

     example application:
     interaction matrix indexed by site1, form1, site2, form2

    \ingroup module_data_structures
    \inpublicapi

    \tparam   T           data type to be stored
    \tparam   Allocator   allocator to be used in creating arrays storing T
                          allocators for storing additional book-keeping data
                          are derived from this allocator via rebind if necessary
*/
template<class T = size_t, typename Allocator = std::allocator<T> >
class IrregArray4D
{
    public:
        /*! adopted these typedefs from AnBe's TriangularArray, for consistency with the other data
            structures, can be retrieved from outside to make sure that the intended types are used
            when accessing or manipulating the array */
        //! the data type stored in the array
        typedef T value_type;
        //! the allocator type to allocate value_type***
        typedef typename Allocator::template rebind<value_type***>::other p3_allocator_type;
        //! the allocator type to allocate value_type**
        typedef typename Allocator::template rebind<value_type**>::other p2_allocator_type;
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
        //! first index of dimension 3
        index_type begin3;
        //! first index of dimension 4
        index_type begin4;
        //! ending index of dimension 1
        index_type end1;
        //! ending index of dimension 2
        index_type end2;
        //! ending index of dimension 3
        index_type end3;
        //! ending index of dimension 4
        index_type end4;
        //! row-specific ending index of dimension 2
        IrregArray1D<size_type> irreg;
        //! row,column-specific ending index of dimension 3
        IrregArray2D<size_type> irreg_mat;
        //! dimension1,2,3-specific ending index of dimension 4
        IrregArray3D<size_type> irreg_ten;
        //! flag says that dimX(i) = endX(i) - beginX(i) + 1 is not a constant
        bool is_irreg;
        //! array storing the the actual data
        value_type ****arr;
        //! the allocator used to allocate arr
        p3_allocator_type p3_allocator;
        //! the allocator used to allocate arr[i]
        p2_allocator_type p2_allocator;
        //! the allocator used to allocate arr[i][j]
        p_allocator_type   p_allocator;
        //! the allocator used to allocate arr[i][j][k]
        allocator_type       allocator;
        //! let other IrregArray4D instantations access private members (for the explicit conversion operator)
        template<class W, typename AW> friend class IrregArray4D;
    public:
        /*! \brief   default constructor without memory allocation, allocateArray has to be called for
                     allocating memory initArray for memory allocation and content initialization
        */
        IrregArray4D()
            : begin1(0), begin2(0), begin3(0), begin4(0), end1(0), end2(0), end3(0), end4(0), is_irreg(false), arr(NULL)
        {
        }
        /*! \brief  copy constructor

            \param[in]   ini_t   array to be copied
         */
        IrregArray4D(const IrregArray4D<value_type> &ini_t)
            : begin1(0), begin2(0), begin3(0), begin4(0), end1(0), end2(0), end3(0), end4(0), is_irreg(false), arr(NULL)
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
            \param[in]    k1       starting index in dimension 4
            \param[in]    k2       last index in dimension 4
         */
        IrregArray4D(const index_type x1, const index_type x2, const index_type y1, const index_type y2, const index_type z1, const index_type z2, const index_type k1, const index_type k2)
            : begin1(0), begin2(0), begin3(0), begin4(0), end1(0), end2(0), end3(0), end4(0), is_irreg(false), arr(NULL)
        {
            initArray(x1, x2, y1, y2, z1, z2, k1, k2);
        }
        /*! \brief   constructor for a regularly shaped array array data with initialization

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
        IrregArray4D(const index_type x1, const index_type x2, const index_type y1, const index_type y2, const index_type z1, const index_type z2, const index_type k1, const index_type k2, const value_type& ini)
            : begin1(0), begin2(0), begin3(0), begin4(0), end1(0), end2(0), end3(0), end4(0), is_irreg(false), arr(NULL)
        {
            initArray(x1, x2, y1, y2, z1, z2, k1, k2, ini);
        }
        /*! \brief   general constructor helper function for irregularly shaped matrices and all special cases

            \param[in]    x1       starting index in dimension 1
            \param[in]    x2       last index in dimension 1
            \param[in]    sizes1   number of elements in dimension 2 as function of the index in dimension 1
            \param[in]    sizes2   number of elements in dimension 3 as function of the indices in the lower dimensions
            \param[in]    sizes3   number of elements in dimension 4 as function of the indices in the lower dimensions
         */
        IrregArray4D(const index_type x1, const index_type x2, const size_type *sizes1, const IrregArray2D<size_type> &sizes2, const IrregArray3D<size_type> &sizes3)
            : begin1(0), begin2(0), begin3(0), begin4(0), end1(0), end2(0), end3(0), end4(0), is_irreg(true), arr(NULL)
        {
            initArray(x1, x2, IrregArray1D<size_type>(x1, x2, sizes1), sizes2, sizes3);
        }
        /*! \brief   general constructor helper function for irregularly shaped matrices and all special cases

            \param[in]    x1       starting index in dimension 1
            \param[in]    x2       last index in dimension 1
            \param[in]    sizes1   number of elements in dimension 2 as function of the index in dimension 1
            \param[in]    sizes2   number of elements in dimension 3 as function of the indices in the lower dimensions
            \param[in]    sizes3   number of elements in dimension 4 as function of the indices in the lower dimensions
            \param[in]    ini      initial value for the array elements
         */
        IrregArray4D(const index_type x1, const index_type x2, const size_type *sizes1, const IrregArray2D<size_type> &sizes2, const IrregArray3D<size_type> &sizes3, const value_type& ini)
            : begin1(0), begin2(0), begin3(0), begin4(0), end1(0), end2(0), end3(0), end4(0), is_irreg(true), arr(NULL)
        {
            initArray(x1, x2, IrregArray1D<size_type>(x1, x2, sizes1), sizes2, sizes3, ini);
        }
        /*! \brief   general constructor helper function for irregularly shaped matrices and all special cases

            \param[in]    x1       starting index in dimension 1
            \param[in]    x2       last index in dimension 1
            \param[in]    sizes1   number of elements in dimension 2 as function of the index in dimension 1
            \param[in]    sizes2   number of elements in dimension 3 as function of the indices in the lower dimensions
            \param[in]    sizes3   number of elements in dimension 4 as function of the indices in the lower dimensions
         */
        IrregArray4D(const index_type x1, const index_type x2, const IrregArray1D<size_type> &sizes1, const IrregArray2D<size_type> &sizes2, const IrregArray3D<size_type> &sizes3)
            : begin1(0), begin2(0), begin3(0), begin4(0), end1(0), end2(0), end3(0), end4(0), is_irreg(true), arr(NULL)
        {
            initArray(x1, x2, sizes1, sizes2, sizes3);
        }
        /*! \brief   general constructor helper function for irregularly shaped matrices and all special cases

            \param[in]    x1       starting index in dimension 1
            \param[in]    x2       last index in dimension 1
            \param[in]    sizes1   number of elements in dimension 2 as function of the index in dimension 1
            \param[in]    sizes2   number of elements in dimension 3 as function of the indices in the lower dimensions
            \param[in]    sizes3   number of elements in dimension 4 as function of the indices in the lower dimensions
            \param[in]    ini      initial value for the array elements
         */
        IrregArray4D(const index_type x1, const index_type x2, const IrregArray1D<size_type> &sizes1, const IrregArray2D<size_type> &sizes2, const IrregArray3D<size_type> &sizes3, const value_type& ini)
            : begin1(0), begin2(0), begin3(0), begin4(0), end1(0), end2(0), end3(0), end4(0), is_irreg(true), arr(NULL)
        {
            initArray(x1, x2, sizes1, sizes2, sizes3, ini);
        }
        /*! \brief   specialized constructor for irregular, symmetrically shaped matrices
                     application: interaction energy matrix site1 form1(site1) site2 form2(site2)

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
        IrregArray4D(const index_type x1, const index_type x2, const size_type *sizes1, const index_type z1, const index_type z2, const size_type *sizes2)
            : begin1(x1), begin2(1), begin3(z1), begin4(1), end1(0), end2(0), end3(0), end4(0), is_irreg(true), arr(NULL)
        {
            initArray(x1, x2, IrregArray1D<size_type>(x1, x2, sizes1), z1, z2, sizes2);
        }
        /*! \brief   specialized constructor for irregular, symmetrically shaped matrices
                     application: interaction energy matrix site1 form1(site1) site2 form2(site2)

            \param[in]    x1       starting index in dimensions 1 and 2
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
            \param[in]    ini      initial value for the array elements
         */
        IrregArray4D(const index_type x1, const index_type x2, const size_type *sizes1, const index_type z1, const index_type z2, const size_type *sizes2, const value_type& ini)
            : begin1(x1), begin2(1), begin3(z1), begin4(1), end1(0), end2(0), end3(0), end4(0), is_irreg(true), arr(NULL)
        {
            initArray(x1, x2, IrregArray1D<size_type>(x1, x2, sizes1), z1, z2, IrregArray1D<size_type>(z1, z2, sizes2), ini);
        }
        /*! \brief   specialized constructor for irregular, symmetrically shaped matrices
                     application: interaction energy matrix site1 form1(site1) site2 form2(site2)

            \param[in]    s1       starting index in dimension 1
            \param[in]    s2       starting index in dimension 2
            \param[in]    s3       starting index in dimension 3
            \param[in]    s4       starting index in dimension 4
            \param[in]    d1       number of elements in dimension 1
            \param[in]    sizes1   number of index values in dimension 4
            \param[in]    sizes2   number of index values in dimension 4
                                   as function of the indices in dimension 3
                                   (in this special case, independent on the
                                   index in dimensions 1 and 2)
            \param[in]    ini      initial value for the array elements
         */
        IrregArray4D(const index_type s1, const index_type s2, const index_type s3, const index_type s4, const index_type d1, const size_type *sizes1, const IrregArray2D<size_type> &sizes2, const value_type& ini)
            : begin1(s1), begin2(s2), begin3(s3), begin4(s4), end1(s1+d1-1), end2(0), end3(0), end4(0), is_irreg(true), arr(NULL)
        {
            initArray(s1, s2, s3, IrregArray1D<size_type>(s1, s1 + d1 - 1, sizes1), sizes2, ini);
        }
        /*! \brief   constructor for irregularly shaped matrices

            \param[in]    s1       starting index in dimension 1
            \param[in]    s2       starting index in dimension 2
            \param[in]    s3       starting index in dimension 3
            \param[in]    s4       starting index in dimension 4
            \param[in]    d1       number of elements in dimension 1
            \param[in]    sizes1   number of elements in dimension 2 as function of the index in dimension 1
            \param[in]    sizes2   number of elements in dimension 3 as function of the indices in the lower dimensions
         */
        IrregArray4D(const index_type s1, const index_type s2, const index_type s3, const index_type s4, const index_type d1, const size_type *sizes1, const IrregArray2D<size_type> &sizes2)
            : begin1(s1), begin2(s2), begin3(s3), begin4(s4), end1(s1 + d1 - 1), end2(0), end3(0), end4(0), is_irreg(true), arr(NULL)
        {
            initArray(s1, s2, s3, IrregArray1D<size_type>(s1, (s1 + d1 - 1), sizes1), sizes2);
        }
        /*! \brief   specialized constructor for irregular, symmetrically shaped matrices
                     application: interaction energy matrix site1 form1(site1) site2 form2(site2)

            \param[in]    x1       starting index in dimensions 1 and 2
            \param[in]    x2       last index in dimension 1
            \param[in]    sizes1   number of elements in dimension 2 as function of the index in dimension 1
            \param[in]    z1       starting index in dimension 3
            \param[in]    z2       last index in dimension 2
            \param[in]    sizes2   number of elements in dimension 3 as function of the indices in the lower dimensions
            \param[in]    ini      initial value for the array elements
         */
        IrregArray4D(const index_type x1, const index_type x2, const IrregArray1D<size_type> &sizes1, const index_type z1, const index_type z2, const IrregArray1D<size_type> &sizes2, const value_type& ini)
            : begin1(x1), begin2(x1), begin3(z1), begin4(z1), end1(x2), end2(0), end3(z2), end4(0), is_irreg(true), arr(NULL)
        {
            initArray(x1, x2, sizes1, z1, z2, sizes2, ini);
        }
        /*! \brief   specialized constructor for irregular, symmetrically shaped matrices
                     application: interaction energy matrix site1 form1(site1) site2 form2(site2)

            \param[in]    x1       starting index in dimensions 1 and 2
            \param[in]    x2       last index in dimension 1
            \param[in]    sizes1   number of elements in dimension 2 as function of the index in dimension 1
            \param[in]    z1       starting index in dimension 3
            \param[in]    z2       last index in dimension 2
            \param[in]    sizes2   number of elements in dimension 3 as function of the indices in the lower dimensions
         */
        IrregArray4D(const index_type x1, const index_type x2, const IrregArray1D<size_type> &sizes1, const index_type z1, const index_type z2, const IrregArray1D<size_type> &sizes2)
            : begin1(x1), begin2(x1), begin3(z1), begin4(z1), end1(x2), end2(0), end3(z2), end4(0), is_irreg(true), arr(NULL)
        {
            initArray(x1, x2, sizes1, z1, z2, sizes2);
        }
        /*! \brief  initialize IrregArray4D<value_type> from IrregArray4D<T_foreign>,
                    used by the copy constructor and by the explicit conversion operator,
                    T_foreign must be convertible to value_type via static_cast

            \param[in]    ini    array to be copied
         */
        template <typename T_foreign, typename A_foreign = std::allocator<T_foreign> >
        void initArray(const IrregArray4D<T_foreign, A_foreign> &ini)
        {
            if (ini.is_irreg)
            {
                initArray(ini.begin1, ini.end1, ini.irreg, ini.irreg_mat, ini.irreg_ten);
            }
            else
            {
                initArray(ini.begin1, ini.end1, ini.begin2, ini.end2, ini.begin3, ini.end3, ini.begin4, ini.end4);
            }
            for(index_type h1=begin1; h1<=end1; h1++)
            {
                for(index_type h2=begin2; h2<=getEnd2(h1); h2++)
                {
                    for(index_type h3=begin3; h3<=getEnd3(h1, h2); h3++)
                    {
                        for(index_type h4=begin4; h4<=getEnd4(h1, h2, h3); h4++)
                        {
                            arr[h1][h2][h3][h4] = static_cast<value_type>(ini(h1, h2, h3, h4));
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
            is_irreg = false;
            begin1 = x1;
            begin2 = y1;
            begin3 = z1;
            begin4 = k1;
            end1 = x2;
            end2 = y2;
            end3 = z2;
            end4 = k2;
            allocateArray(begin1, end1, begin2, end2, begin3, end3, begin4, end4);
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
        void initArray(const index_type x1, const index_type x2, const index_type y1, const index_type y2, const index_type z1, const index_type z2, const index_type k1, const index_type k2, const value_type& ini)
        {
            initArray(x1, x2, y1, y2, z1, z2, k1, k2);
            for(index_type h1 = x1; h1 <= x2; h1++)
            {
                for(index_type h2 = y1; h2 <= y2; h2++)
                {
                    for(index_type h3 = z1; h3 <= z2; h3++)
                    {
                        for(index_type h4 = k1; h4 <= k2; h4++)
                        {
                            arr[h1][h2][h3][h4] = ini;
                        }
                    }
                }
            }
        }
        /*! \brief general constructor helper function for irregularly shaped matrices and all special cases

            \param[in]    x1       starting index in dimension 1
            \param[in]    y1       starting index in dimension 2
            \param[in]    z1       starting index in dimension 3
            \param[in]    k1       starting index in dimension 4
            \param[in]    sizes1   number of elements in dimension 2 as function of the index in dimension 1
            \param[in]    sizes2   number of elements in dimension 3 as function of the indices in the lower dimensions
            \param[in]    sizes3   number of elements in dimension 4 as function of the indices in the lower dimensions
         */
        void initArray(const index_type x1, const index_type y1, const index_type z1, const index_type k1, const IrregArray1D<size_type> &sizes1, const IrregArray2D<size_type> &sizes2, const IrregArray3D<size_type> &sizes3)
        {
         deallocateArray();
         is_irreg= true;
         begin1 = x1;
         begin2 = y1;
         begin3 = z1;
         begin4 = k1;
         end1 = x1 + sizes1.getLength1() - 1;
         end2 = 0;
         end3 = 0;
         end4 = 0;
         allocateArray(begin1, begin2, begin3, begin4, sizes1, sizes2, sizes3);
        }
        /*! \brief general constructor helper function for irregularly shaped matrices and all special cases

            \param[in]    x1       starting index in dimension 1
            \param[in]    y1       starting index in dimension 2
            \param[in]    z1       starting index in dimension 3
            \param[in]    k1       starting index in dimension 4
            \param[in]    sizes1   number of elements in dimension 2 as function of the index in dimension 1
            \param[in]    sizes2   number of elements in dimension 3 as function of the indices in the lower dimensions
            \param[in]    sizes3   number of elements in dimension 4 as function of the indices in the lower dimensions
            \param[in]    ini      initialization value for the data array elements
         */
        void initArray(const index_type x1, const index_type y1, const index_type z1, const index_type k1, const IrregArray1D<size_type> &sizes1, const IrregArray2D<size_type> &sizes2, const IrregArray3D<size_type> &sizes3, const value_type& ini)
        {
            initArray(x1, y1, z1, k1, sizes1, sizes2, sizes3);
            for(index_type h1 = x1; h1 <= getEnd1(); h1++)
            {
                for(index_type h2 = y1; h2 <= getEnd2(h1); h2++)
                {
                    for(index_type h3 = z1; h3 <= getEnd3(h1,h2); h3++)
                    {
                        for(index_type h4 = k1; h4 <= getEnd4(h1,h2,h3); h4++)
                        {
                            arr[h1][h2][h3][h4] = ini;
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
        void initArray(const index_type x1, const index_type x2, const IrregArray1D<size_type> &sizes1, const IrregArray2D<size_type> &sizes2, const IrregArray3D<size_type> &sizes3)
        {
            deallocateArray();
            is_irreg = true;
            begin1 = x1;
            begin2 = x1;
            begin3 = x1;
            begin4 = x1;
            end1 = x2;
            end2 = 0;
            end3 = 0;
            end4 = 0;
            irreg.initArray(x1, x2);
            for(index_type i = getBegin1(); i <= x2; i++)
            {
                irreg[i] = sizes1(i);
            }
            irreg_mat.initArray(x1, x2, irreg);
            for(index_type i = getBegin1(); i <= x2; i++)
            {
                for(index_type j = getBegin2(); j <= getEnd2(i); j++)
                {
                    irreg_mat(i, j) = sizes2(i, j);
                }
            }
            irreg_ten.initArray(x1, x2, irreg, irreg_mat);
            for(index_type i = getBegin1(); i <= x2; i++)
            {
                for(index_type j = getBegin2(); j <= getEnd2(i); j++)
                {
                    for(index_type k = getBegin3(); k <= getEnd3(i, j); k++)
                    {
                        irreg_ten(i, j, k) = sizes3(i, j, k);
                    }
                }
            }
            allocateArray(begin1, begin2, begin3, begin4, irreg, irreg_mat, irreg_ten);
        }
        /*! \brief general constructor helper function for irregularly shaped matrices and all special cases

            \param[in]    x1       starting index in dimension 1
            \param[in]    x2       last index in dimension 1
            \param[in]    sizes1   number of elements in dimension 2 as function of the index in dimension 1
            \param[in]    sizes2   number of elements in dimension 3 as function of the indices in the lower dimensions
            \param[in]    sizes3   number of elements in dimension 4 as function of the indices in the lower dimensions
            \param[in]    ini      initialization value for the data array elements
         */
        void initArray(const index_type x1, const index_type x2, const IrregArray1D<size_type> &sizes1, const IrregArray2D<size_type> &sizes2, const IrregArray3D<size_type> &sizes3, const value_type& ini)
        {
            initArray(x1, x2, sizes1, sizes2, sizes3);
            for(index_type h1 = getBegin1(); h1 <= getEnd1(); h1++)
            {
                for(index_type h2 = getBegin2(); h2 <= getEnd2(h1); h2++)
                {
                    for(index_type h3 = getBegin3(); h3 <= getEnd3(h1,h2); h3++)
                    {
                        for(index_type h4 = getBegin4(); h4 <= getEnd4(h1,h2,h3); h4++)
                        {
                            arr[h1][h2][h3][h4] = ini;
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
        void initArray(const index_type x1, const index_type x2, const size_type *sizes1, const IrregArray2D<size_type> &sizes2, const IrregArray3D<size_type> &sizes3)
        {
            initArray(x1, x2, IrregArray1D<size_type>(x1, x2, sizes1), sizes2, sizes3);
        }
        /*! \brief specialized constructor helper function for irregular, symmetrically shaped matrices
                   application: interaction energy matrix site1 form1(site1) site2 form2(site2)

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
        void initArray(const index_type x1, const index_type x2, const IrregArray1D<size_type> &sizes1, const index_type z1, const index_type z2, const IrregArray1D<size_type> &sizes2)
        {
            deallocateArray();
            is_irreg = true;
            begin1 = x1;
            begin2 = begin1;
            begin3 = z1;
            begin4 = begin3;
            end1 = x2;
            end2 = 0;
            end3 = z2;
            end4 = 0;
            irreg.initArray(begin1, end1);
            for(index_type i = begin1; i <= end1; i++)
            {
                irreg[i] = sizes1[i];
            }
            irreg_mat.initArray(begin1, end1, sizes1, (size_type)0);
            const index_type dim3 = z2 - z1 + 1;
            for(index_type i = begin1; i <= end1; i++)
            {
                for(index_type j = begin2; j <= getEnd2(i); j++)
                {
                    irreg_mat(i, j) = dim3;
                }
            }
            irreg_ten.initArray(begin1, end1, irreg, irreg_mat, (size_type)0);
            for(index_type i = begin1; i <= end1; i++)
            {
                for(index_type j = begin2; j <= getEnd2(i); j++)
                {
                    for(index_type k = begin3; k <= getEnd3(i, j); k++)
                    {
                        irreg_ten(i, j, k) = sizes2[k];
                    }
                }
            }
            allocateArray(begin1, begin2, begin3, begin4, irreg, irreg_mat, irreg_ten);
        }
        /*! \brief specialized constructor helper function for irregular, symmetrically shaped matrices
                   application: interaction energy matrix site1 form1(site1) site2 form2(site2)


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
        void initArray(const index_type x1, const index_type x2, const IrregArray1D<size_type> &sizes1, const index_type z1, const index_type z2, const IrregArray1D<size_type> &sizes2, const value_type& ini)
        {
            initArray(x1, x2, sizes1, z1, z2, sizes2);
            for(index_type h1 = getBegin1(); h1 <= getEnd1(); h1++)
            {
                for(index_type h2 = getBegin2(); h2 <= getEnd2(h1); h2++)
                {
                    for(index_type h3 = getBegin3(); h3 <= getEnd3(h1, h2); h3++)
                    {
                        for(index_type h4 = getBegin4(); h4 <= getEnd4(h1, h2, h3); h4++)
                        {
                            arr[h1][h2][h3][h4] = ini;
                        }
                    }
                }
            }
        }
        /*! \brief specialized constructor helper function for irregular, symmetrically shaped matrices
                   application: interaction energy matrix site1 form1(site1) site2 form2(site2)

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
        void initArray(const index_type x1, const index_type x2, const size_type *sizes1, const index_type z1, const index_type z2, const size_type *sizes2, const value_type& ini)
        {
            initArray(x1, x2, IrregArray1D<size_type>(x1, x2, sizes1), z1, z2, IrregArray1D<size_type>(z1, z2, sizes2), ini);
        }
        /*! \brief specialized constructor helper function for irregular, symmetrically shaped matrices
                   application: interaction energy matrix site1 form1(site1) site2 form2(site2)

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
        void initArray(const index_type x1, const index_type x2, const size_type *sizes1, const index_type z1, const index_type z2, const size_type *sizes2)
        {
            initArray(x1, x2, IrregArray1D<size_type>(x1, x2, sizes1), z1, z2, IrregArray1D<size_type>(z1, z2, sizes2));
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
        /*! \brief number of elements in dimension 2 at lower dimension index x

            \param[in]    x       index in dimension 1 for which the array length in dimension 2 is requested */
        size_type getLength2(const index_type x) const
        {
            return (is_irreg ? irreg[x] : (end2 - begin2 + 1));
        }
        //! number of elements in dimension 3 (regular array)
        size_type getLength3() const
        {
            if (!is_irreg) {
                return (end3 - begin3 + 1);
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
            return (is_irreg ? irreg_mat(x,y) : (end3 - begin3 + 1));
        }
        //! number of elements in dimension 4 (regular array)
        size_type getLength4() const
        {
            if (!is_irreg) {
                return (end4 - begin4 + 1);
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
            return (is_irreg ? irreg_ten(x,y,z) : (end4 - begin4 + 1));
        }
        //! get the first index of dimension 1
        index_type getBegin1() const { return begin1; }
        //! get the first index of dimension 2
        index_type getBegin2() const { return begin2; }
        //! get the first index of dimension 3
        index_type getBegin3() const { return begin3; }
        //! get the first index of dimension 4
        index_type getBegin4() const { return begin4; }
        //! get the last index of dimension 1
        index_type getEnd1() const { return end1; }
        //! get the last index of dimension 2
        index_type getEnd2() const
        {
            if (!is_irreg) {
                return end2;
            }
            else
            {
                std::fprintf(stderr, "Error in %s::getEnd2(): called without argument for an irregular array.", typeid(*this).name());
                throw InternalError("Extent of a variable dimension in an irregular array requested without specifying the index/indices in the lower dimension(s).");
            }
        }
        /*! \brief get the last index of dimension 2 at lower dimension index x

            \param[in]    x       index in dimension 1 for which the last index in dimension 2 is requested */
        index_type getEnd2(const index_type x) const
        {
            return (is_irreg ? (irreg[x] + begin2 - 1) : end2);
        }
        //! get the last index of dimension 3
        index_type getEnd3() const
        {
            if (!is_irreg) {
                return end3;
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
            return (is_irreg ? (irreg_mat(x, y) + begin3 - 1) : end3);
        }
        //! get the last index of dimension 4
        index_type getEnd4() const
        {
            if (!is_irreg) {
                return end4;
            }
            else
            {
                std::fprintf(stderr, "Error in %s::getEnd4(): called without arguments for an irregular array.", typeid(*this).name());
                throw InternalError("Extent of a variable dimension in an irregular array requested without specifying the index/indices in the lower dimension(s).");
            }
        }
        /*! \brief get the last index of dimension 3 for lower dimension indices x, y, z

            \param[in]    x       index in dimension 1
            \param[in]    y       index in dimension 2
            \param[in]    z       index in dimension 3                       */
        index_type getEnd4(const index_type x, const index_type y, const index_type z) const
        {
            return (is_irreg ? (irreg_ten(x, y, z) + begin4 - 1) : end4);
        }

        //! returns the amount of memory occupied by the object
        size_type getSize() const
        {
            // the book-keeping data
            size_type size = sizeof(*this);
            size += irreg.getSize()     - sizeof(IrregArray1D<size_type>);
            size += irreg_mat.getSize() - sizeof(IrregArray2D<size_type>);
            size += irreg_ten.getSize() - sizeof(IrregArray3D<size_type>);
            // the storage data
            size += getLength1() * sizeof(value_type***);
            for(index_type i = getBegin1(); i <= getEnd1(); i++)
            {
                size += getLength2(i) * sizeof(value_type**);
                for(index_type j = getBegin2(); j <= getEnd2(i); j++)
                {
                    size += getLength3(i,j) * sizeof(value_type*);
                    for(index_type k = getBegin3(); k <= getEnd3(i, j); k++)
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
            size_type nelements = 0;
            for(index_type i = getBegin1(); i <= getEnd1(); i++)
            {
                for(index_type j = getBegin2(); j <= getEnd2(i); j++)
                {
                    for(index_type k = getBegin3(); k <= getEnd3(i, j); k++)
                    {
                        nelements += getLength4(i, j, k);
                    }
                }
            }
            return nelements;
        }

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
            const index_type d1 = x2 - x1 + 1;
            const index_type d2 = y2 - y1 + 1;
            const index_type d3 = z2 - z1 + 1;
            const index_type d4 = k2 - k1 + 1;
            if(d1 < 0)
            {
                std::fprintf(stderr, "Error in %s::allocateArray(): called with illegal bounds begin1 > end1.", typeid(*this).name());
                throw InternalError("Trying to create an irregular array with illegal bounds begin1 > end1");
            }
            if(d2 < 0)
            {
                std::fprintf(stderr, "Error in %s::allocateArray(): called with illegal bounds begin2 > end2.", typeid(*this).name());
                throw InternalError("Trying to create an irregular array with illegal bounds begin2 > end2");
            }
            if(d3 < 0)
            {
                std::fprintf(stderr, "Error in %s::allocateArray(): called with illegal bounds begin3 > end3.", typeid(*this).name());
                throw InternalError("Trying to create an irregular array with illegal bounds begin3 > end3");
            }
            if(d4 < 0)
            {
                std::fprintf(stderr, "Error in %s::allocateArray(): called with illegal bounds begin4 > end4.", typeid(*this).name());
                throw InternalError("Trying to create an irregular array with illegal bounds begin4 > end4");
            }
            try
            {
                arr = p3_allocator.allocate(d1);
                arr = arr - begin1;
                arr[x1] = p2_allocator.allocate(d1*d2);
                arr[x1] -= begin2;
                arr[x1][y1] = p_allocator.allocate(d1*d2*d3);
                arr[x1][y1] -= begin3;
                arr[x1][y1][z1] = allocator.allocate(d1*d2*d3*d4);
                arr[x1][y1][z1] -= begin4;
                // set second dimension
                for(index_type i = y1 + 1; i <= y2; i++)
                {
                    arr[x1][i] = arr[x1][i-1] + d3;
                }
                // set third dimension
                for(index_type i = x1 + 1; i <= x2; i++)
                {
                    arr[i] = arr[i-1] + d2;
                    arr[i][y1] = arr[i-1][y1] + (d2 * d3);
                    arr[i][y1][z1] = arr[i-1][y1][z1] + (d2 * d3 * d4);
                    for(index_type j = y1 + 1; j <= y2; j++)
                    {
                        arr[i][j] = arr[i][j-1]+d3;
                    }
                }
                for(index_type i = x1; i <= x2; i++)
                {
                    for(index_type j = y1 + 1;j <= y2; j++)
                    {
                        arr[i][j][z1] = arr[i][j-1][z1] + (d3 * d4);
                    }
                }
                for(index_type i = x1; i <= x2; i++)
                {
                    for(index_type j = y1; j <= y2; j++)
                    {
                        for(index_type h = z1 + 1; h <= z2; h++)
                        {
                            arr[i][j][h] = arr[i][j][h-1] + d4;
                        }
                    }
                }
            }
            catch(std::bad_alloc)
            {
                std::fprintf(stderr, "Error in %s::init_Arr(): could not allocate memory for the data array.", typeid(*this).name());
                throw;
            };

            //! initialize the array elements
            for (index_type i = getBegin1(); i <= getEnd1(); ++i)
            {
                for (index_type j = getBegin2(); j <= getEnd2(); ++j)
                {
                    for (index_type k = getBegin3(); k <= getEnd3(); ++k)
                    {
                        for (index_type l = getBegin4(); l <= getEnd4(); ++l)
                        {
                            allocator.construct(&this->operator()(i, j, k, l), value_type());
                        }
                    }
                }
            }

        }
        /*! \brief allocate memory for the array

            \param[in]    s1       starting index in dimension 1
            \param[in]    s2       starting index in dimension 2
            \param[in]    s3       starting index in dimension 3
            \param[in]    s4       starting index in dimension 4
            \param[in]    d1       number of elements in dimension 1
            \param[in]    sizes1   number of elements in dimension 2 as function of the index in dimension 1
            \param[in]    sizes2   number of elements in dimension 3 as function of the indices in the lower dimensions
            \param[in]    sizes3   number of elements in dimension 4 as function of the indices in the lower dimensions
        */
        void allocateArray(const index_type s1, const index_type s2, const index_type s3, const index_type s4, const index_type d1, size_type *sizes1, const IrregArray2D<size_type> &sizes2, const IrregArray3D<size_type> &sizes3)
        {
            allocateArray(s1, s2, s3, s4, IrregArray1D<size_type>(s1, (s1 + d1 - 1), sizes1), sizes2, sizes3);
        }
        /*! \brief allocate memory for the array

            \param[in]    s1       starting index in dimension 1
            \param[in]    s2       starting index in dimension 2
            \param[in]    s3       starting index in dimension 3
            \param[in]    s4       starting index in dimension 4
            \param[in]    sizes1   number of elements in dimension 2 as function of the index in dimension 1
            \param[in]    sizes2   number of elements in dimension 3 as function of the indices in the lower dimensions
            \param[in]    sizes3   number of elements in dimension 4 as function of the indices in the lower dimensions
        */
        void allocateArray(const index_type s1, const index_type s2, const index_type s3, const index_type s4, const IrregArray1D<size_type> &sizes1, const IrregArray2D<size_type> &sizes2, const IrregArray3D<size_type> &sizes3)
        {
            try
            {
                arr = p3_allocator.allocate(sizes1.getLength1());
                arr -= s1;
                for (index_type i = s1; i <= sizes1.getEnd1(); ++i)
                {
                    arr[i] = p2_allocator.allocate(sizes1[i]);
                    arr[i] -= s2;
                }
                for (index_type i = s1; i <= sizes1.getEnd1(); i++)
                {
                    for (index_type j = s2; j <= (s2 + static_cast<index_type>(sizes1[i]) - 1); ++j)
                    {
                        arr[i][j] = p_allocator.allocate(sizes2(i, j));
                        arr[i][j] -= s3;
                    }
                }
                for (index_type i = s1; i <= sizes1.getEnd1(); i++)
                {
                    for (index_type j = s2; j <= (s2 + static_cast<index_type>(sizes1[i]) - 1); ++j)
                    {
                        for(index_type t = s3; t <= (s3 + static_cast<index_type>(sizes2(i, j)) - 1); t++)
                        {
                            arr[i][j][t] = allocator.allocate(sizes3(i, j, t));
                            arr[i][j][t] -= s4;
                        }
                    }
                }
            }
            catch(std::bad_alloc)
            {
                std::fprintf(stderr, "Error in %s::init_Arr(): could not allocate memory for the data array.", typeid(*this).name());
                throw;
            };

            //! initialize the array elements
            for (index_type i = getBegin1(); i <= getEnd1(); ++i)
            {
                for (index_type j = getBegin2(); j <= getEnd2(i); ++j)
                {
                    for (index_type k = getBegin3(); k <= getEnd3(i, j); ++k)
                    {
                        for (index_type l = getBegin4(); l <= getEnd4(i, j, k); ++l)
                        {
                            allocator.construct(&this->operator()(i, j, k, l), value_type());
                        }
                    }
                }
            }
        }
        //! deallocate the memory of all constituent arrays
        void deallocateArray()
        {
            if (arr != NULL)
            {
                //! destroy the array elements
                for (index_type i = getBegin1(); i <= getEnd1(); ++i)
                {
                    for (index_type j = getBegin2(); j <= getEnd2(i); ++j)
                    {
                        for (index_type k = getBegin3(); k <= getEnd3(i, j); ++k)
                        {
                            for (index_type l = getBegin4(); l <= getEnd4(i, j, k); ++l)
                            {
                                allocator.destroy(&this->operator()(i, j, k, l));
                            }
                        }
                    }
                }

                // deallocate data array
                if (is_irreg)
                {
                    for (index_type i = getBegin1(); i <= getEnd1(); i++)
                    {
                            for (index_type j = getBegin2(); j <= getEnd2(i); j++)
                            {
                                for(index_type t = getBegin3(); t <= getEnd3(i, j); t++)
                                {
                                    arr[i][j][t] += begin4;
                                    allocator.deallocate(arr[i][j][t], getLength4(i, j, t));
                                }
                                arr[i][j] += begin3;
                                p_allocator.deallocate(arr[i][j], getLength3(i, j));
                            }
                            arr[i] += begin2;
                            p2_allocator.deallocate(arr[i], getLength2(i));
                    }
                    arr += begin1;
                    p3_allocator.deallocate(arr, getLength1());
                }
                else
                {
                    const size_type d1 = getLength1();
                    const size_type d2 = getLength2();
                    const size_type d3 = getLength3();
                    const size_type d4 = getLength4();
                    arr[begin1][begin2][begin3] += begin4;
                    allocator.deallocate(arr[begin1][begin2][begin3], d1*d2*d3*d4);
                    arr[begin1][begin2] += begin3;
                    p_allocator.deallocate(arr[begin1][begin2], d1*d2*d3);
                    arr[begin1] += begin2;
                    p2_allocator.deallocate(arr[begin1], d1*d2);
                    arr += begin1;
                    p3_allocator.deallocate(arr, d1);
                }
            }
            arr = NULL;
            // deallocate book-keeping arrays
            if (is_irreg)
            {
                irreg_ten.deallocateArray();
                irreg_mat.deallocateArray();
                irreg.deallocateArray();
                is_irreg = false;
            }
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
        class Proxy {
            private:
                value_type*** _array;
                class Proxy2 {
                    private:
                        value_type** _array;
                        class Proxy3 {
                            private:
                                value_type* _array;
                            public:
                                Proxy3(value_type* array) : _array(array) { }

                                value_type& operator[](index_type k)
                                {
                                    return _array[k];
                                }
                                const value_type& operator[](index_type k) const
                                {
                                    return _array[k];
                                }
                        };
                    public:
                        Proxy2(value_type** array) : _array(array) { }

                        Proxy3 operator[](index_type z)
                        {
                            return Proxy3(_array[z]);
                        }
                        const Proxy3 operator[](index_type z) const
                        {
                            return _array[z];
                        }
                };
                //! let IrregArray4D access private members
                template<class W, typename AW> friend class IrregArray4D;
            public:
                Proxy(value_type*** array) : _array(array) { }

                Proxy2 operator[](index_type y)
                {
                    return Proxy2(_array[y]);
                }
                const Proxy2 operator[](index_type y) const
                {
                    return Proxy2(_array[y]);
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
        Proxy operator[](index_type x) {
           return Proxy(arr[x]);
        }
        /*! \brief   const index operator[x][y][z][k] for array-like usage realized with chained functors "Proxy*"
            \param[in]   x   index in dimension 1
            \cond
            \param[in]   y   index in dimension 2
            \param[in]   z   index in dimension 3
            \param[in]   k   index in dimension 4
            \endcond
        */
        const Proxy operator[](index_type x) const {
           return Proxy(arr[x]);
        }

        /*! \brief   index operator(x,y,z) for functor-like usage
            \param[in]   x   index in dimension 1
            \param[in]   y   index in dimension 2
            \param[in]   z   index in dimension 3
            \param[in]   k   index in dimension 4
        */
        inline value_type &operator()(const index_type x, const index_type y, const index_type z, const index_type k)
        {
            return arr[x][y][z][k];
        }
        /*! \brief   const index operator(x,y,z) for functor-like usage
            \param[in]   x   index in dimension 1
            \param[in]   y   index in dimension 2
            \param[in]   z   index in dimension 3
            \param[in]   k   index in dimension 4
        */
        inline const value_type &operator()(const index_type x, const index_type y, const index_type z, const index_type k) const
        {
            return arr[x][y][z][k];
        }
        /*! \brief   index operator(x,y,z) for functor-like usage
            \param[in]   x   index in dimension 1
            \param[in]   y   index in dimension 2
            \param[in]   z   index in dimension 3
        */
        inline value_type *operator()(const index_type x, const index_type y, const index_type z)
        {
            return arr[x][y][z];
        }
        /*! \brief   const index operator(x,y,z) for functor-like usage
            \param[in]   x   index in dimension 1
            \param[in]   y   index in dimension 2
            \param[in]   z   index in dimension 3
        */
        inline const value_type *operator()(const index_type x, const index_type y, const index_type z) const
        {
            return arr[x][y][z];
        }
        /*! \brief   index operator(x,y,z) for functor-like usage
            \param[in]   x   index in dimension 1
            \param[in]   y   index in dimension 2
        */
        inline value_type **operator()(const index_type x, const index_type y)
        {
            return arr[x][y];
        }
        /*! \brief   const index operator(x,y,z) for functor-like usage
            \param[in]   x   index in dimension 1
            \param[in]   y   index in dimension 2
        */
        inline const value_type **operator()(const index_type x, const index_type y) const
        {
            return arr[x][y];
        }
        /*! \brief   index operator(x,y,z) for functor-like usage
            \param[in]   x   index in dimension 1
        */
        inline value_type ***operator()(const index_type x)
        {
            return arr[x];
        }
        /*! \brief   const index operator(x,y,z) for functor-like usage
            \param[in]   x   index in dimension 1
        */
        inline const value_type ***operator()(const index_type x) const
        {
            return arr[x];
        }

        //! get a pointer to the data array
        const value_type**** getArray() const
        {
            return arr;
        }
        //! get a pointer to the data array
        value_type**** getArray()
        {
            return arr;
        }

        /*! \brief   stream operator for convenient printing of the array to an ostream
            \param[in]   output   output stream for the array contents
            \param[in]   ten      the array to be printed
        */
        friend std::ostream& operator<<(std::ostream &output, const IrregArray4D &ten)
        {
            for(index_type i = ten.begin1; i<=ten.end1; i++)
            {
               for(index_type j = ten.begin2; j<= ten.getEnd2(i); j++)
               {
                  output << "------------------------------" << std::endl << i << ", " << j << std::endl;
                  for(index_type k = ten.begin3; k <= ten.getEnd3(i, j); k++)
                  {
                     for(index_type l = ten.begin4; l <= ten.getEnd4(i, j, k); l++)
                     {
                         output << std::setw(10) << std::fixed << std::setprecision(2) << ten(i, j, k, l) << " ";
                     }
                     output << std::endl;
                  }
                  output << std::endl;
               }
               output << std::endl;
            }
            output << std::endl;
            return output;
        }
        /*! \brief   assignment operator *this = rhs
            \param[in]   rhs   right-hand side object of the assignment
        */
        IrregArray4D<value_type>& operator=(const IrregArray4D<value_type> &rhs)
        {
            deallocateArray();
            initArray(rhs);
            return *this;
        }
        /*! \brief  conversion operator (explicit to avoid unintended implicit conversion), usable via
                    IrregArray4D<T_foreign> newArr = static_cast<IrregArray4D<T_foreign> >(FlatIrregArray<T> oldArr);
        */
        template<typename T_foreign, typename A_foreign = std::allocator<T_foreign> >
#if __cplusplus < 201103L
            explicit
#endif
            operator IrregArray4D<T_foreign, A_foreign>() const
        {
            IrregArray4D<T_foreign, A_foreign> result;
            result.initArray(*this);
            return result;
        }
};

} // end namespace gmx

#endif
