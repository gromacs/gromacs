/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2019- The GROMACS Authors
 * and the project initiators Erik Lindahl, Berk Hess and David van der Spoel.
 * Consult the AUTHORS/COPYING files and https://www.gromacs.org for details.
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
 * https://www.gnu.org/licenses, or write to the Free Software Foundation,
 * Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA.
 *
 * If you want to redistribute modifications to GROMACS, please
 * consider that scientific software is very special. Version
 * control is crucial - bugs must be traceable. We will be happy to
 * consider code for inclusion in the official distribution, but
 * derived work must not be called official GROMACS. Details are found
 * in the README & COPYING files - if they are missing, get the
 * official version at https://www.gromacs.org.
 *
 * To help us fund GROMACS development, we humbly ask that you cite
 * the research papers on the package. Check out https://www.gromacs.org.
 */
/*! \libinternal
 * \file
 * \brief Declares MultiDimArray.
 *
 * \author Christian Blau <cblau@gwdg.de>
 * \ingroup module_math
 * \ingroup module_mdspan
 * \inlibraryapi
 */

#ifndef GMX_MATH_MULTIDIMARRAY_H_
#define GMX_MATH_MULTIDIMARRAY_H_

#include <type_traits>
#include <utility>

#include "gromacs/mdspan/mdspan.h"
#include "gromacs/utility/arrayref.h"

namespace gmx
{

namespace detail
{
//! Same as std::void_t from C++17
template<class...>
using void_t = void;

template<typename T, typename = void>
struct is_resizable : std::false_type
{
};

template<typename T>
struct is_resizable<T, void_t<decltype(std::declval<T>().resize(size_t(0)))>> : std::true_type
{
};

//! Type has a resize member function callable with size_t argument
template<typename T>
// NOLINTNEXTLINE misc-definitions-in-headers
constexpr bool is_resizable_v = is_resizable<T>::value;
} // namespace detail

/*! \libinternal \brief
 * Multidimensional array that manages its own memory.
 *
 * \note No bounds checking when accessing memory
 *
 * \note That the view holds a valid pointer to the data is a class invariant.
 *
 * The Container type that stores the data may be resizable (std::vector or similar)
 * or static (std::array or similar). Copy and move assignment routines as well as
 * swapping are designed to yield good performances in both cases, notably
 * foregoing the copy-swap idiom due to the bad performance in swapping std::array.
 *
 * This class avoids throwing exceptions, apart from the ones that might be thrown
 * from the containers during resizing an allocation. (bad_malloc from std::vector
 * is a likely candidate)
 *
 *
 * Holds as many elements as required by a multidimensional view.
 * \tparam TContainer   Data type container for the data to be stored
 *                      as MultiDimArray with random element access and
 *                      value_type, reference and const_reference exposed
 * \tparam Extents      An extents class describing the array dimensions
 *                      as used in module_mdspan
 * \tparam LayoutPolicy The data layout as in module_mdspan describes
 *                      translation of indices to memory offset.
 *                      Defaults to right aligned, so that the right-most index
 *                      is contiguous in memory.
 */
template<class TContainer, class Extents, class LayoutPolicy = layout_right>
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
    //! const view on the data
    using const_view_type = basic_mdspan<const value_type, Extents, LayoutPolicy>;
    /*! \brief Iterator type for contiguous iteration over the stored data.
     * Used, e.g., in free begin and end functions
     */
    using iterator = typename ArrayRef<value_type>::iterator;
    /*! \brief Const iterator type for contiguous iteration over the stored data.
     *  used, e.g., in free begin and end functions
     */
    using const_iterator = const typename ArrayRef<const value_type>::const_iterator;

    static_assert(detail::is_resizable_v<TContainer> == (Extents::rank_dynamic() > 0),
                  "Resizable container (e.g. std::vector) requires at least one dynamic rank. "
                  "Non-resizable container (e.g. std::array) requires zero dynamic ranks.");

    /*! \brief
     * Allocate dynamic array data and set view with the dynamic extents.
     *
     * \param[in] dynamicExtent A parameter pack that describes the dynamic
     *                          size of the array. Empty if purely static.
     *
     * \tparam IndexType        Parameter pack type holding the dynamic
     *                          extents of the multidimensional array
     */
    template<class... IndexType, typename T = TContainer, typename = typename std::enable_if_t<detail::is_resizable_v<T>>>
    MultiDimArray(IndexType... dynamicExtent)
    {
        resize(dynamicExtent...);
    }
    /*! \brief
     * Construction from fixed sized arrays if the array size is static and
     * layout policy allows compile time determination of the container size.
     *
     * Enables the expected initialization
     * MultiDimArray<std::array<float, 9>, extents<3,3>> arr = {{1,2...}}
     * \tparam T Template parameter for activation via SFINAE.
     */
    // SFINAE required because std::vector::size isn't constexpr and is_constexpr doesn't exist.
    template<typename T = TContainer, typename = typename std::enable_if_t<!detail::is_resizable_v<T>>>
    constexpr MultiDimArray(const TContainer& data = {}) noexcept : data_(data), view_(data_.data())
    {
        static_assert(TContainer().size() == typename view_type::mapping_type().required_span_size(),
                      "Non-resizable container type size must match static MultiDimArray size.");
    }
    //! Copy constructor
    constexpr MultiDimArray(const MultiDimArray& o) :
        data_(o.data_), view_(data_.data(), o.view_.extents())
    {
    }
    //! Move constructor
    MultiDimArray(MultiDimArray&& o) noexcept :
        data_(std::move(o.data_)), view_(data_.data(), o.view_.extents())
    {
    }
    //! Copy assignment
    MultiDimArray& operator=(const MultiDimArray& o)
    {
        data_ = o.data_;
        view_ = view_type(data_.data(), o.view_.extents());
        return *this;
    }
    //! Move assignment
    MultiDimArray& operator=(MultiDimArray&& o) noexcept
    {
        data_ = std::move(o.data_);
        view_ = view_type(data_.data(), o.view_.extents());
        return *this;
    }
    //! Swaps content with other
    void swap(MultiDimArray& o) noexcept
    {
        using std::swap;
        swap(data_, o.data_);
        // swap(view_, o.view_) also swaps the pointer to the data and thus does not work
        // instead, restore the view as class invariant after the data swapping operation
        o.view_ = view_type(o.data_.data(), view_.extents());
        view_   = view_type(data_.data(), o.view_.extents());
    }
    /*! \brief
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
    template<class... IndexType>
    void resize(IndexType... dynamicExtent)
    {
        // use a mapping object to determine the required span size;
        layout_right::mapping<Extents> map{ Extents{ dynamicExtent... } };
        data_.resize(map.required_span_size());
        // to construct a valid view on the data, the container has to be resized before
        // the assignment, so that data_.data() is valid
        view_ = view_type(data_.data(), dynamicExtent...);
    }
    /*! \brief Data access via multidimensional indices.
     * This allows referencing rank R array elements as array(x_0,x_1,x_2, .., x_R)
     *
     * \param[in] index multidimensional indices as parameter pack
     *                  the number of parameters must match the rank of the array.
     *
     * \returns reference to array element
     */
    template<class... IndexType>
    reference operator()(IndexType... index) noexcept
    {
        return view_(index...);
    }
    /*! \brief Const data access via multidimensional indices.
     * This allows referencing rank R array elements as array(x_0,x_1,x_2, .., x_R)
     *
     * \param[in] index multidimensional indices as parameter pack
     *                  the number of parameters must match the rank of the array.
     *
     * \returns const reference to array element
     */
    template<class... IndexType>
    constexpr const_reference operator()(IndexType... index) const noexcept
    {
        return view_(index...);
    }
    /*! \brief Contiguous access to the data.
     * \returns ArrayRef to stored data.
     */
    ArrayRef<value_type> toArrayRef() { return { data_.data(), data_.data() + data_.size() }; }
    /*! \brief Contiguous const access to the data.
     * \returns ArrayRef to stored data.
     */
    constexpr ArrayRef<const value_type> toArrayRef() const
    {
        return { data_.data(), data_.data() + data_.size() };
    }
    /*! \brief Return the extent.
     * \param[in] k dimension to query for extent
     * \returns extent along specified dimension
     */
    constexpr typename view_type::index_type extent(int k) const noexcept
    {
        return view_.extent(k);
    }
    //! Conversion to multidimensional view on the data
    constexpr view_type asView() noexcept { return view_; }
    //! Conversion to const multidimensional view on the data
    constexpr const_view_type asConstView() const noexcept
    {
        return { data_.data(), view_.mapping() };
    }

private:
    //! The contiguous data that is equipped with multidimensional indexing in this class
    TContainer data_;
    //! Multidimensional view into data_.
    view_type view_;
};

//! Free MultiDimArray begin function addressing its contiguous memory.
template<class TContainer, class Extents>
constexpr typename MultiDimArray<TContainer, Extents>::const_iterator
begin(const MultiDimArray<TContainer, Extents>& multiDimArray)
{
    return multiDimArray.toArrayRef().begin();
}

//! Free MultiDimArray begin function addressing its contiguous memory.
template<class TContainer, class Extents>
constexpr typename MultiDimArray<TContainer, Extents>::iterator begin(MultiDimArray<TContainer, Extents>& multiDimArray)
{
    return multiDimArray.toArrayRef().begin();
}

//! Free MultiDimArray end function addressing its contiguous memory.
template<class TContainer, class Extents>
constexpr typename MultiDimArray<TContainer, Extents>::const_iterator
end(const MultiDimArray<TContainer, Extents>& multiDimArray)
{
    return multiDimArray.toArrayRef().end();
}

//! Free MultiDimArray end function addressing its contiguous memory.
template<class TContainer, class Extents>
constexpr typename MultiDimArray<TContainer, Extents>::iterator end(MultiDimArray<TContainer, Extents>& multiDimArray)
{
    return multiDimArray.toArrayRef().end();
}

//! Swap function
template<class TContainer, class Extents>
void swap(MultiDimArray<TContainer, Extents>& a, MultiDimArray<TContainer, Extents>& b) noexcept
{
    a.swap(b);
}

} // namespace gmx

#endif // GMX_MATH_MULTIDIMARRAY_H_
