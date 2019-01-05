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
        //! the view used to access the data
        using view_type = basic_mdspan<value_type, Extents>;
        //! iterator type for contiguous iteration over the stored data, used, e.g., in free begin and end functions
        using iterator = typename TContainer::value_type *;
        //! const iterator type for contiguous iteration over the stored data, used, e.g., in free begin and end functions
        using const_iterator = const typename TContainer::value_type *;

        /*! \libinternal \brief
         * Set the dynamic extents of the array if any and set container size accordingly.
         *
         * Allows canonical copy constructor while setting dynamic ranks from parameter packs.
         * If this were a constructor, this would always match better than the
         * copy constructor, thus rendering the copy constructor unusable.

         * \param[in] dynamicExtent A parameter pack that describes the dynamic size of the array.
         *                          Empty if purely static.
         *
         * \tparam TContainer       The container of the multidimensional data
         * \tparam Extents          The extents of the multidimensional array
         * \tparam IndexType        Parameter pack type holding the dynamic extents of the multidimensional array
         *
         * \returns A MultiDimArray with memory for its dynamic extents allocated
         */
        template <class ... IndexType>
        static MultiDimArray<TContainer, Extents>
        allocate(IndexType && ... dynamicExtent)
        {
            MultiDimArray<TContainer, Extents> result;
            // use a mapping object to determine the required span size;
            layout_right::mapping<Extents>     map {
                Extents {
                    std::forward<IndexType>(dynamicExtent) ...
                }
            };
            result.data_.resize(map.required_span_size());
            // to construct a valid view on the data, the container has to be resized before
            // the assignment, so that data_.data() is valid
            result.view_ = view_type(result.data_.data(), dynamicExtent ...);
            return result;
        }

        //! Copy constructor
        constexpr MultiDimArray(const MultiDimArray &o) : data_ {o.data_}, //calls the copy constructor of the data, then sets a new matching view to the copied data
        view_ {view_type(data_.data(), o.view_.extents())}
        {}
        //! Move constructor
        MultiDimArray(MultiDimArray &&o) noexcept :
            data_ {},
        view_ {}
        {
            swap(o);
        }
        //! Move assignment operator
        MultiDimArray &operator=(MultiDimArray &&o) noexcept
        {
            data_ = {};
            view_ = {};
            swap(o);
            return *this;
        }
        //! Copy assignment operator
        MultiDimArray &operator=(const MultiDimArray &o) noexcept
        {
            MultiDimArray<TContainer, Extents> tmp(o);
            tmp.swap(*this);
            return *this;
        }
        //! Swaps content with other
        void swap(MultiDimArray &o)
        {
            using std::swap;
            swap(data_, o.data_);
            swap(view_, o.view_);
        }
        /*! \brief Data access via multidimensional indices.
         * This allows referencing rank R array elements as array(x_0,x_1,x_2, .., x_R)
         *
         * \param[in] index multidimensional indices as parameter pack
         *                  the number of parameters must match the rank of the array.
         *
         * \returns reference to array element
         */
        template< class... IndexType>
        reference operator()(IndexType && ... index)
        {
            return view_(std::forward<IndexType>(index) ...);
        }
        /*! \brief Const data access via multidimensional indices.
         * This allows referencing rank R array elements as array(x_0,x_1,x_2, .., x_R)
         *
         * \param[in] index multidimensional indices as parameter pack
         *                  the number of parameters must match the rank of the array.
         *
         * \returns const reference to array element
         */
        template< class... IndexType>
        constexpr const_reference operator()(IndexType && ... index) const
        {
            return view_(std::forward<IndexType>(index) ...);
        }
        /*! \brief Contiguous access to the data.
         * \returns ArrayRef to stored data.
         */
        ArrayRef<value_type> toArrayRef()
        {
            return {data_.data(), data_.data()+data_.size()};
        }
        /*! \brief Contiguous const access to the data.
         * \returns ArrayRef to stored data.
         */
        constexpr ArrayRef<const value_type> toArrayRef() const
        {
            return {data_.data(), data_.data()+data_.size()};
        }
        /*! \brief Return the extent.
         * \param[in] k dimension to query for extent
         * \returns extent along specified dimension
         */
        constexpr typename view_type::index_type extent(int k) const noexcept
        {
            return view_.extent(k);
        }

    private:
        //! Hidden default constructor to enforce use of allocate function to generate this object
        MultiDimArray() = default;
        //! The contiguous data that is equipped with multidimensional indexing in this class
        TContainer data_;
        //! Multidimensional view into data_.
        view_type  view_;
};

//! Free MultiDimArray begin function addressing its contiguous memory.
template <class TContainer, class Extents>
typename MultiDimArray<TContainer, Extents>::const_iterator
begin(const MultiDimArray<TContainer, Extents> &b)
{
    return b.toArrayRef().begin();
}

//! Free MultiDimArray begin function addressing its contiguous memory.
template <class TContainer, class Extents>
typename MultiDimArray<TContainer, Extents>::iterator
begin(MultiDimArray<TContainer, Extents> &b)
{
    return b.toArrayRef().begin();
}

//! Free MultiDimArray end function addressing its contiguous memory.
template <class TContainer, class Extents>
typename MultiDimArray<TContainer, Extents>::const_iterator
end(const MultiDimArray<TContainer, Extents> &b)
{
    return b.toArrayRef().end();
}

//! Free MultiDimArray end function addressing its contiguous memory.
template <class TContainer, class Extents>
typename MultiDimArray<TContainer, Extents>::iterator
end(MultiDimArray<TContainer, Extents> &b)
{
    return b.toArrayRef().end();
}

//! Swap function
template <class TContainer, class Extents>
void swap(MultiDimArray<TContainer, Extents> &a, MultiDimArray<TContainer, Extents > &b) noexcept
{
    a.swap(b);
}

}      // namespace gmx

#endif // GMX_MATH_MULTIDIMARRAY_H_
