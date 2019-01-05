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
/*! \libinternal
 * \file
 * \brief Declares MultiDimArray.
 *
 * \author Christian Blau <cblau@gwdg.de>
 * \ingroup module_math
 * \ingroup module_mdspan
 */

#ifndef GMX_MATH_MULTIDIMARRAY_H_
#define GMX_MATH_MULTIDIMARRAY_H_

#include "gromacs/mdspan/mdspan.h"
#include "gromacs/utility/arrayref.h"

namespace gmx
{

/*! \libinternal
 * \brief
 * Multidimensional array that manages its own memory.
 *
 * Holds as many elements as required by a multidimensional view.
 * \tparam TContainer Data type container for the data to be stored
 *                    as MultiDimArray with random element access and
 *                    value_type, refence and const_reference exposed
 * \tparam Extents    An extents class describing the array dimensions
 *                    as used in module_mdspan
 */
template <class TContainer, class Extents>
class MultiDimArray
{
    public:
        //! the type of values that are stored
        using value_type = typename TContainer::value_type;
        //! reference type to the stored values
        using reference = typename TContainer::reference;
        //! const reference type to the stored values
        using const_reference = typename TContainer::const_reference;

        /*! \libinternal \brief
         * Set the dynamic extents of the array if any and set container size accordingly.
         *
         * \param[in] index An parameter pack that describes the dynamic size of the array.
         *
         * \tparam TContainer The container of the multidimensional data
         * \tparam Extents The extents of the multidimensional array
         * \tparam IndexType Parameter pack holding the dynamic extents of the multidimensional array
         *
         * \returns A MultiDimArray
         */
        template <class ... IndexType>
        static MultiDimArray<TContainer, Extents>
        allocate(IndexType && ... index)
        {
            MultiDimArray<TContainer, Extents> result;
            static_assert( sizeof ... (index) == Extents::rank_dynamic(),
                           "Number of indices must match the dynamic rank." );
            // use a mapping object to determine the required span size;
            layout_right::mapping<Extents>     map {
                Extents {
                    std::forward<IndexType>(index) ...
                }
            };
            result.data_.resize(map.required_span_size());
            // to construct a valid view on the data, the container has to be resized before
            // the assignment, so that data_.data() is valid
            result.view_ = basic_mdspan<value_type, Extents>(result.data_.data(), index ...);
            return result;
        }

        //! Copy constructor
        constexpr MultiDimArray(const MultiDimArray &o) :
            view_ {o.view_},
        data_ {o.data_}
        {}
        //! Move constructor
        MultiDimArray(const MultiDimArray &&o) :
            view_(std::move(o.view_)),
            data_(std::move(o.data_))
        {}
        //! Move assignment operator
        MultiDimArray &operator=(const MultiDimArray &&o)
        {
            if (&o != this)
            {
                view_ = std::move(o.view_);
                data_ = std::move(o.data_);
            }
            return *this;
        }

        //! Copy assignment operator
        MultiDimArray &operator=(const MultiDimArray &o)
        {
            if (&o != this)
            {
                view_ = o.view_;
                data_ = o.data_;
            }
            return *this;
        }
        /*! \brief Data access via multidimensional indices.
         * \param[in] index multidimensional indices as parameter
         * pack; the number of parameters matching the rank of the array.
         * This allows referencing rank R array elements as
         * array(x_0,x_1,x_2, .., x_R)
         * \returns reference to array element
         */
        template< class... IndexType>
        reference operator()(IndexType && ... index)
        {
            static_assert( sizeof ... (index) == Extents::rank(),
                           "Number of indices must match rank" );
            return view_(std::forward<IndexType>(index) ...);
        }
        /*! \brief Const data access via multidimensional indices.
         * \param[in] index multidimensional indices as parameter
         * pack; the number of parameters matching the rank of the array.
         * This allows referencing rank R array elements as
         * array(x_0,x_1,x_2, .., x_R)
         * \returns const reference to array element
         */
        template< class... IndexType>
        const_reference operator()(IndexType && ... index) const
        {
            static_assert( sizeof ... (index) == Extents::rank(),
                           "Number of indices must match rank" );
            return view_(std::forward<IndexType>(index) ...);
        }
        /*! \brief Linearised access to the data.
         * \returns ArrayRef to stored data.
         */
        ArrayRef<value_type> data()
        {
            return {data_.data(), data_.data()+data_.size()};
        }
        /*! \brief Linearised const access to the data.
         * \returns ArrayRef to stored data.
         */
        ArrayRef<const value_type> data() const
        {
            return {data_.data(), data_.data()+data_.size()};
        }
    private:
        //! Private default constructor for use in builder
        MultiDimArray() = default;
        /*! \brief  Declaring a friend builder class.
         * Allows canonical copy constructor while setting dynamic ranks from parameter packs.
         * If this were a constructor, this would always match better than the
         * copy constructor, thus rendering the copy constructor unusable.
         */
        template <class T, class E, class ... IndexType>
        MultiDimArray<T, E>
        friend makeMultiDimArray(IndexType && ... index);
        //! Multidimensional view on to the stored data.
        basic_mdspan<value_type, Extents> view_;
        //! The linear data that is equipped with multidimensional indexing in this class
        TContainer  data_;
};

}      // namespace gmx

#endif // GMX_MATH_MULTIDIMARRAY_H_
