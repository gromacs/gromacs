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

namespace detail
{
/*! \internal \brief
 * Determine static array size at compile time.
 * Statically evaluates the product of static array extents.
 * \tparam Extents Extents of an multidimensional array as in module_mdspan
 * \returns the product of the static extents
 */
template <typename Extents>
constexpr typename Extents::index_type staticExtentsProduct(size_t const i = 0)
{
    return (i < Extents::rank()) ? Extents::static_extent(i) * staticExtentsProduct<Extents>(i + 1) : 1;
}
} // namespace detail
/*! \libinternal
 * \brief
 * Multidimensional array that manages its own memory.
 *
 * Holds as many elements as required by a multidimensional view.
 * \tparam TContainer   Data type container for the data to be stored
 *                      as MultiDimArray with random element access and
 *                      value_type, refence and const_reference exposed
 * \tparam Extents      An extents class describing the array dimensions
 *                      as used in module_mdspan
 * \tparam LayoutPolicy The data layout as in module_mdspan describes
 *                      translation of indices to memory offset.
 *                      Defaults to right aligned, so that the right-most index
 *                      is contiguous in memory.
 */
template <class TContainer, class Extents, class LayoutPolicy = layout_right>
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
        using view_type = basic_mdspan<value_type, Extents, LayoutPolicy>;
        /*! \brief Iterator type for contiguous iteration over the stored data.
         * Used, e.g., in free begin and end functions
         */
        using iterator = typename ArrayRef<value_type>::iterator;
        /*! \brief Const iterator type for contiguous iteration over the stored data.
         *  used, e.g., in free begin and end functions
         */
        using const_iterator = const typename ArrayRef<const value_type>::const_iterator;

        /*! \brief
         * Allocate dynamic array data and set view with the dynamic extents.
         *
         * \param[in] dynamicExtent A parameter pack that describes the dynamic
         *                          size of the array. Empty if purely static.
         *
         * \tparam IndexType        Parameter pack type holding the dynamic
         *                          extents of the multidimensional array
         */
        template <class ... IndexType>
        MultiDimArray(IndexType... dynamicExtent)
        {
            resize(dynamicExtent ...);
        }
        /*! \brief
         * Construction from fixed sized arrays if the array size is static and
         * layout policy allows compile time determination of the container size.
         *
         * Enables ,e.g., MultiDimArray<int[8], extents<2,2,2>> or
         * MultiDimArray<std::array<float, 9>, extents<3,3>>
         * \tparam StaticSizeAndLayout Template parameter for SNFINAE activation.
         */
        template <bool StaticSizeAndLayoutRight = Extents::rank_dynamic() == 0 &&
                      std::is_same<LayoutPolicy, layout_right>::value,
                  typename std::enable_if<StaticSizeAndLayoutRight, int>::type = 0>
        constexpr MultiDimArray() : view_(data_.data())
        {
            // \todo replace staticExtentsProduct by the required_span_size
            // of the LayoutPolicy, once required_span_size is a constexpr
            // with C++14
            static_assert(std::tuple_size<TContainer>() ==
                          detail::staticExtentsProduct<Extents>(),
                          "Non-resizable container type size must match static MultiDimArray size.");
        }

        /*! \libinternal \brief
         * Resize the dynamic extents of the array if any and set container size
         * accordingly.
         *
         * Invalidates data and views of this array.
         *
         * \param[in] dynamicExtent A parameter pack that describes the dynamic
         *                          size of the array. Empty if purely static.
         * \tparam IndexType        Parameter pack type holding the dynamic
         *                          extents of the multidimensional array
         */
        template <class ... IndexType,
                  typename std::enable_if<sizeof ... (IndexType) != 0, int>::type = 0>
        void resize(IndexType... dynamicExtent)
        {
            // use a mapping object to determine the required span size;
            layout_right::mapping<Extents> map { Extents { dynamicExtent ...}};
            data_.resize(map.required_span_size());
            // to construct a valid view on the data, the container has to be resized before
            // the assignment, so that data_.data() is valid
            view_ = view_type(data_.data(), dynamicExtent ...);
        }
        //! Copy constructor
        constexpr MultiDimArray(const MultiDimArray &o) : data_ {o.data_},
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
        reference operator()(IndexType ... index)
        {
            return view_(index ...);
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
        constexpr const_reference operator()(IndexType ... index) const
        {
            return view_(index ...);
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
        /*! \brief
         * Hide default constructor when size is dynamic to enforce RAII with proper
         * resizing and view setting.
         */
        template <bool DynamicInSize = Extents::rank_dynamic() != 0,
                  typename std::enable_if<DynamicInSize, int>::type = 0, typename T = TContainer>
        MultiDimArray()
        {}
        //! The contiguous data that is equipped with multidimensional indexing in this class
        TContainer data_;
        //! Multidimensional view into data_.
        view_type  view_;
};

//! Free MultiDimArray begin function addressing its contiguous memory.
template <class TContainer, class Extents>
constexpr typename MultiDimArray<TContainer, Extents>::const_iterator
begin(const MultiDimArray<TContainer, Extents> &multiDimArray)
{
    return multiDimArray.toArrayRef().begin();
}

//! Free MultiDimArray begin function addressing its contiguous memory.
template <class TContainer, class Extents>
constexpr typename MultiDimArray<TContainer, Extents>::iterator
begin(MultiDimArray<TContainer, Extents> &multiDimArray)
{
    return multiDimArray.toArrayRef().begin();
}

//! Free MultiDimArray end function addressing its contiguous memory.
template <class TContainer, class Extents>
constexpr typename MultiDimArray<TContainer, Extents>::const_iterator
end(const MultiDimArray<TContainer, Extents> &multiDimArray)
{
    return multiDimArray.toArrayRef().end();
}

//! Free MultiDimArray end function addressing its contiguous memory.
template <class TContainer, class Extents>
constexpr typename MultiDimArray<TContainer, Extents>::iterator
end(MultiDimArray<TContainer, Extents> &multiDimArray)
{
    return multiDimArray.toArrayRef().end();
}

//! Swap function
template <class TContainer, class Extents>
void swap(MultiDimArray<TContainer, Extents> &a, MultiDimArray<TContainer, Extents > &b) noexcept
{
    a.swap(b);
}

}      // namespace gmx

#endif // GMX_MATH_MULTIDIMARRAY_H_
