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
/*! \file
 * \brief Defines helper types for class enumerations.
 *
 * These helper types facilitate iterating over class enums, and
 * maintaining a type-safe and value-safe matching list of names. The
 * code is closely based on the public-domain code by Guilherme
 * R. Lampert, found in commit c94c18a of
 * https://github.com/glampert/enum_helpers/blob/master/enum_helpers.hpp
 * Thanks Guilherme!
 *
 * NOTE This functionality only works for enumerations of monotonically
 * increasing values, starting with the value zero.
 *
 * Usage examples:
 *
 *  enum class Foo : int
 *  {
 *      Bar,
 *      Baz,
 *      Fooz,
 *      Count
 *  };
 *
 *  for (Foo c : EnumerationWrapper<Foo>{})
 *  {
 *      // 'c' is a constant from Foo
 *  }
 *
 *
 *  const EnumerationArray<Foo, std::string> fooStrings = { { "Bar", "Baz", "Fooz" } };
 *  std::cout << fooStrings[Foo::Baz];
 *  std::cout << fooStrings[Foo::Count]; // Triggers an assertion
 *
 *  for (Foo c : keysOf(fooStrings))
 *  {
 *      print(fooStrings[c]);
 *  }
 *
 *  ArrayRef<const std::string> namesRef(fooStrings);
 *
 * \author Mark Abraham <mark.j.abraham@gmail.com>
 * \inlibraryapi
 * \ingroup module_utility
 */
#ifndef GMX_UTILITY_ENUMHELPERS_H
#define GMX_UTILITY_ENUMHELPERS_H

#include <cstddef>

#include <iterator>
#include <type_traits>

#if __has_include("external/boost/stl_interfaces/iterator_interface.hpp")
#    include "external/boost/stl_interfaces/iterator_interface.hpp"
#else // fallback for installed headers
#    include <gromacs/external/boost/stl_interfaces/iterator_interface.hpp>
#endif

#include "gromacs/utility/gmxassert.h"

namespace gmx
{

/*! \libinternal
 * \brief Allows iterating sequential enumerators.
 *
 * You can also provide an increment step > 1 if each constant is
 * spaced by a larger value.  Terminating constant is assumed to be a
 * 'Count' member, which is never iterated. A different name for the
 * terminating constant can also be specified on declaration.
 *
 * NOTE This functionality only works for enumerations of monotonically
 * increasing values, starting with the value zero.
 *
 * See file documentation for usage example.
 *
 * \tparam  EnumType   The enum (class) type.
 * \tparam  Last       Last constant or number thereof (assumes a default 'Count' member).
 * \tparam  Step       Step increment.
 */
template<typename EnumType, EnumType Last = EnumType::Count, std::ptrdiff_t Step = 1>
class EnumerationIterator final :
    public gmx::boost::stl_interfaces::iterator_interface<EnumerationIterator<EnumType, Last, Step>, std::random_access_iterator_tag, EnumType>
{
public:
    // TODO: Use std::is_enum_v when CUDA 11 is a requirement.
    static_assert(std::is_enum<EnumType>::value, "Enumeration iterator must be over an enum type.");
    //! Convenience alias
    using IntegerType = std::underlying_type_t<EnumType>;

    constexpr EnumerationIterator() noexcept : m_current{ 0 } // Assumes 0 is the first constant
    {
    }
    //! Conversion constructor
    explicit constexpr EnumerationIterator(const EnumType enumIndex) noexcept :
        m_current(static_cast<IntegerType>(enumIndex))
    {
    }
    //! Addition-assignment operator
    constexpr EnumerationIterator& operator+=(std::ptrdiff_t i) noexcept
    {
        m_current += Step * i;
        return *this;
    }
    //! Dereference operator
    constexpr EnumType operator*() const noexcept
    {
        GMX_ASSERT(m_current < static_cast<IntegerType>(Last), "dereferencing out of range");
        return static_cast<EnumType>(m_current);
    }
    //! Difference operator
    constexpr std::ptrdiff_t operator-(const EnumerationIterator other) const noexcept
    {
        return (static_cast<std::ptrdiff_t>(m_current) - static_cast<std::ptrdiff_t>(other.m_current)) / Step;
    }

private:
    IntegerType m_current;
};

/*! \libinternal
 * \brief Allows constructing iterators for looping over sequential enumerators.
 *
 * These are particularly useful for range-based for statements.
 *
 * You can also provide an increment step > 1 if each constant is
 * spaced by a larger value.  Terminating constant is assumed to be a
 * 'Count' member, which is never iterated. A different name for the
 * terminating constant can also be specified on declaration.
 *
 * See file documentation for usage example.
 *
 * \tparam  EnumType   The enum (class) type.
 * \tparam  Last       Last constant or number thereof (assumes a default 'Count' member).
 * \tparam  Step       Step increment.
 */
template<typename EnumType, EnumType Last = EnumType::Count, unsigned int Step = 1>
class EnumerationWrapper final
{
public:
    //! Convenience alias.
    using IteratorType = EnumerationIterator<EnumType, Last, Step>;
    //! Convenience alias
    using value_type = EnumType;

    //! Functions required for range-based for statements to work.
    /*!@{*/
    IteratorType begin() const { return IteratorType{}; }
    IteratorType end() const { return IteratorType{ Last }; }
    /*!@}*/
};

/*! \libinternal
 * \brief Wrapper for a C-style array with size and indexing defined
 * by an enum. Useful for declaring arrays of enum names for debug
 * or other printing. An ArrayRef<DataType> may be constructed from
 * an object of this type.
 *
 * See file documentation for usage example.
 *
 * Note that if clang-tidy gives strange errors referring to the line
 * number of the struct declaration, these likely refer to the
 * compiler-generated constructors. Simplification of the calling code
 * might eliminate that call and thus the clang-tidy error.
 *
 * \tparam  EnumType   The enum (class) type.
 * \tparam  DataType   Type of the data stored in the array.
 * \tparam  ArraySize  Size in entries of the array.
 */
template<typename EnumType, typename DataType, EnumType ArraySize = EnumType::Count>
struct EnumerationArray final // NOLINT
{
    //! Convenience alias
    using EnumerationWrapperType = EnumerationWrapper<EnumType, ArraySize>;

    //! Convenience alias
    using value_type = DataType;

    /*! \brief Data for names.
     *
     * Data is kept public so we can use direct aggregate
     * initialization just like in a plain C-style array. */
    DataType m_elements[std::size_t(ArraySize)];

    //! Returns an object that provides iterators over the keys.
    static constexpr EnumerationWrapperType keys() { return EnumerationWrapperType{}; }
    //! Returns the size of the enumeration.
    constexpr std::size_t size() const { return std::size_t(ArraySize); }

    /*!@{*/
    //! Array access with asserts:
    DataType& operator[](const std::size_t arrayIndex)
    {
        GMX_ASSERT(arrayIndex < size(), "index out of range");
        return m_elements[arrayIndex];
    }
    const DataType& operator[](const std::size_t arrayIndex) const
    {
        GMX_ASSERT(arrayIndex < size(), "index out of range");
        return m_elements[arrayIndex];
    }

    DataType& operator[](const EnumType enumIndex)
    {
        GMX_ASSERT(std::size_t(enumIndex) < size(), "index out of range");
        return m_elements[std::size_t(enumIndex)];
    }
    const DataType& operator[](const EnumType enumIndex) const
    {
        GMX_ASSERT(std::size_t(enumIndex) < size(), "index out of range");
        return m_elements[std::size_t(enumIndex)];
    }
    /*!@}*/

    /*!@{*/
    //! Range iterators (unchecked)
    using iterator               = DataType*;
    using const_iterator         = const DataType*;
    using reverse_iterator       = std::reverse_iterator<iterator>;
    using const_reverse_iterator = std::reverse_iterator<const_iterator>;
    /*!@}*/

    /*!@{*/
    //! Getters for forward iterators for ranges
    iterator                 begin() { return &m_elements[0]; }
    iterator                 end() { return &m_elements[size()]; }
    constexpr const_iterator begin() const { return &m_elements[0]; }
    constexpr const_iterator end() const { return &m_elements[size()]; }
    /*!@}*/

    /*!@{*/
    //! Getters for reverse iterators for ranges
    reverse_iterator                 rbegin() { return reverse_iterator{ end() }; }
    reverse_iterator                 rend() { return reverse_iterator{ begin() }; }
    constexpr const_reverse_iterator rbegin() const { return const_reverse_iterator{ end() }; }
    constexpr const_reverse_iterator rend() const { return const_reverse_iterator{ begin() }; }
    /*!@}*/

    /*!@{*/
    //! Pointers (unchecked)
    using pointer       = DataType*;
    using const_pointer = const DataType*;
    /*!@}*/

    //! Returns a const raw pointer to the contents of the array.
    const_pointer data() const { return &m_elements[0]; }
    //! Returns a raw pointer to the contents of the array.
    pointer data() { return &m_elements[0]; }
};

/*! \brief Returns an object that provides iterators over the keys
 * associated with \c EnumerationArrayType.
 *
 * This helper function is useful in contexts where there is an object
 * of an EnumerationArray, and we want to use a range-based for loop
 * over the keys associated with it, and it would be inconvenient to
 * use the very word EnumerationArray<...> type, nor introduce a using
 * statement for this purpose. It is legal in C++ to call a static
 * member function (such as keys()) via an object rather than the
 * type, but clang-tidy warns about that. So instead we make available
 * a free function that calls that static method. */
template<typename EnumerationArrayType>
typename EnumerationArrayType::EnumerationWrapperType keysOf(const EnumerationArrayType& /* arrayObject */)
{
    return EnumerationArrayType::keys();
}

} // namespace gmx

#endif
