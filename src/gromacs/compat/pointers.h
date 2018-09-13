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

/*! \libinternal \file
 * \brief Provides ported functions/classes from gsl/pointers
 *
 * Adapted from the Guidelines Support Library at
 * https://github.com/Microsoft/GSL
 *
 * \author Mark Abraham <mark.j.abraham@gmail.com>
 * \ingroup module_compat
 * \inlibraryapi
 */
 #ifndef GMX_COMPAT_POINTERS_H
 #define GMX_COMPAT_POINTERS_H

#include <type_traits>
#include <utility>

#include "gromacs/utility/gmxassert.h"

namespace gmx
{
namespace compat
{

namespace detail
{

/*! \brief Helpers found in C++14
 *
 * \todo Replace these with std flavours when we support C++14 */
//! \{
template <bool b, class T = void>
using enable_if_t         = typename std::enable_if<b, T>::type;

template<class T>
using remove_cv_t = typename std::remove_cv<T>::type;

template<class T>
using remove_reference_t = typename std::remove_reference<T>::type;
//! \}

}       // namespace detail

/*! \def GMX_RELEASE_CONSTEXPR
 * \brief
 * Macro for debug assertions in constexpr functions.
 *
 * If NDEBUG is defined, this macro expands to constexpr, which in C++11
 * cannot contain the kind of logic used in GMX_ASSERT. If NDEBUG is not
 * defined, then this is a debug build, which should assert, and thus
 * cannot yet be constexpr.
 *
 * \todo Eliminate this when we required C++14.
 */
#ifdef NDEBUG
#define GMX_RELEASE_CONSTEXPR constexpr
#else
#define GMX_RELEASE_CONSTEXPR
#endif

/*! \libinternal
 * \brief Restricts a pointer or smart pointer to only hold non-null values.
 *
 * Has zero size overhead over T.
 *
 * If T is a pointer (i.e. T == U*) then
 * - allow construction from U*
 * - disallow construction from nullptr_t
 * - disallow default construction
 * - ensure construction from null U* fails
 * - allow implicit conversion to U*
 *
 * \todo Eliminate this when we require a version of C++ that supports
 * std::not_null.
 */
template <class T>
class not_null
{
    public:
        static_assert(std::is_assignable<T &, std::nullptr_t>::value, "T cannot be assigned nullptr.");

        //! Move constructor. Asserts in debug mode if \c is nullptr.
        template <typename U, typename = detail::enable_if_t<std::is_convertible<U, T>::value> >
        GMX_RELEASE_CONSTEXPR explicit not_null(U &&u) : ptr_(std::forward<U>(u))
        {
#ifndef NDEBUG
            GMX_ASSERT(ptr_ != nullptr, "Precondition violation: null pointers cannot be used to construct not_null pointers");
#endif
        }

        //! Simple constructor. Asserts in debug mode if \c u is nullptr.
        template <typename = detail::enable_if_t<!std::is_same<std::nullptr_t, T>::value> >
        GMX_RELEASE_CONSTEXPR explicit not_null(T u) : ptr_(u)
        {
#ifndef NDEBUG
            GMX_ASSERT(ptr_ != nullptr, "Precondition violation: null pointers cannot be used to construct not_null pointers");
#endif
        }

        //! Copy constructor.
        template <typename U, typename = detail::enable_if_t<std::is_convertible<U, T>::value> >
        constexpr not_null(const not_null<U> &other) : not_null(other.get())
        {
        }

        //! Default constructors and assignment.
        //! \{
        not_null(not_null &&other) noexcept        = default;
        not_null(const not_null &other)            = default;
        not_null &operator=(const not_null &other) = default;
        //! \}

        //! Getters
        //! \{
        GMX_RELEASE_CONSTEXPR T get() const
        {
#ifndef NDEBUG
            GMX_ASSERT(ptr_ != nullptr, "Postcondition violation: Null pointers cannot be returned from not_null pointers get() method");
#endif
            return ptr_;
        }

        constexpr operator T() const { return get(); }
        constexpr T operator->() const { return get(); }
        //! \}

        //! Deleted to prevent compilation when someone attempts to assign a null pointer constant.
        //! \{
        not_null(std::nullptr_t)            = delete;
        not_null &operator=(std::nullptr_t) = delete;
        //! \}

        //! Deleted unwanted operators because pointers only point to single objects.
        //! \{
        not_null &operator++()                = delete;
        not_null &operator--()                = delete;
        not_null operator++(int)              = delete;
        not_null operator--(int)              = delete;
        not_null &operator+=(std::ptrdiff_t)  = delete;
        not_null &operator-=(std::ptrdiff_t)  = delete;
        void operator[](std::ptrdiff_t) const = delete;
        //! \}

    private:
        T ptr_;
};

//! Convenience function for making not_null pointers.
template <class T>
not_null<T> make_not_null(T &&t)
{
    return not_null < detail::remove_cv_t < detail::remove_reference_t<T>>>{
               std::forward<T>(t)
    };
}

//! Operators to compare not_null pointers.
//! \{
template <class T, class U>
auto operator==(const not_null<T> &lhs, const not_null<U> &rhs)->decltype(lhs.get() == rhs.get())
{
    return lhs.get() == rhs.get();
}

template <class T, class U>
auto operator!=(const not_null<T> &lhs, const not_null<U> &rhs)->decltype(lhs.get() != rhs.get())
{
    return lhs.get() != rhs.get();
}

template <class T, class U>
auto operator<(const not_null<T> &lhs, const not_null<U> &rhs)->decltype(lhs.get() < rhs.get())
{
    return lhs.get() < rhs.get();
}

template <class T, class U>
auto operator<=(const not_null<T> &lhs, const not_null<U> &rhs)->decltype(lhs.get() <= rhs.get())
{
    return lhs.get() <= rhs.get();
}

template <class T, class U>
auto operator>(const not_null<T> &lhs, const not_null<U> &rhs)->decltype(lhs.get() > rhs.get())
{
    return lhs.get() > rhs.get();
}

template <class T, class U>
auto operator>=(const not_null<T> &lhs, const not_null<U> &rhs)->decltype(lhs.get() >= rhs.get())
{
    return lhs.get() >= rhs.get();
}
//! \}

//! Deleted unwanted arithmetic operators.
//! \{
template <class T, class U>
std::ptrdiff_t operator-(const not_null<T> &, const not_null<U> &) = delete;
template <class T>
not_null<T> operator-(const not_null<T> &, std::ptrdiff_t) = delete;
template <class T>
not_null<T> operator+(const not_null<T> &, std::ptrdiff_t) = delete;
template <class T>
not_null<T> operator+(std::ptrdiff_t, const not_null<T> &) = delete;
//! \}

} // namespace compat
} // namespace gmx

#endif
