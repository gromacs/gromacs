/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2020, by the GROMACS development team, led by
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
/*! \inpublicapi \file
 * \brief
 * Implements nblib utilities
 *
 * \author Victor Holanda <victor.holanda@cscs.ch>
 * \author Joe Jordan <ejjordan@kth.se>
 * \author Prashanth Kanduri <kanduri@cscs.ch>
 * \author Sebastian Keller <keller@cscs.ch>
 * \author Artem Zhmurov <zhmurov@gmail.com>
 */

#ifndef NBLIB_UTIL_USER_H
#define NBLIB_UTIL_USER_H

#include <functional>
#include <iostream>
#include <sstream>
#include <string>
#include <tuple>
#include <type_traits>
#include <vector>

#include "nblib/basicdefinitions.h"
#include "nblib/vector.h"

namespace gmx
{
template<typename T>
class ArrayRef;
} // namespace gmx

namespace nblib
{

//! Generate velocities from a Maxwell Boltzmann distribution, masses should be the
//! same as the ones specified for the Topology object
std::vector<Vec3> generateVelocity(real Temperature, unsigned int seed, std::vector<real> const& masses);

//! Check within the container of gmx::RVecs for a NaN or inf
bool isRealValued(gmx::ArrayRef<const Vec3> values);

//! Zero a cartesian buffer
void zeroCartesianArray(gmx::ArrayRef<Vec3> cartesianArray);

//! Used to ignore unused arguments of a lambda functions
inline void ignore_unused() {}

//! Variadic argument version of the ignore_unused function
template<class T, class... Ts>
inline void ignore_unused(T& x, Ts&... xs)
{
    static_cast<void>(x);
    ignore_unused(xs...);
}

/*! \brief A template to create structs as a type-safe alternative to using declarations
 *
 * \inpublicapi
 * \ingroup nblib
 *
 * Used in public API functions where a distinction between different
 * arguments of the same underlying type is desired. This provides a type-safe
 * version to using declarations. Instead of naming a type alias, the name
 * is used to define a struct that inherits from StrongType<T>, where T is
 * the underlying type. For example:
 *
 * struct C6 : StrongType<real>
 * {
 *     using StrongType::StrongType;
 *     using StrongType::operator=;
 * };
 *
 * Due to the T() conversion and assignment from T,
 * an instance of the resulting C6 struct behaves essentially like a real, while construction
 * from real is disabled. This makes it impossible to pass a real as a function parameter
 * of type C6.
 */
template<class T, class Phantom>
struct StrongType
{
    //! default ctor
    StrongType() : value_{} {}
    //! construction from the underlying type T, implicit conversions disabled
    explicit StrongType(T v) : value_(std::move(v)) {}

    //! assignment from T
    StrongType& operator=(T v)
    {
        value_ = std::move(v);
        return *this;
    }

    //! conversion to T
    operator T() const { return value_; }

    //! access the underlying value
    T value() const { return value_; }

private:
    T value_;
};

//! Equality comparison. For the case where a comparison between StrongTypes with matching T, but differing Phantom
//! parameters is desired, the underlying value attribute should be compared instead
template<class T, class Phantom>
[[maybe_unused]] inline bool operator==(const StrongType<T, Phantom>& lhs, const StrongType<T, Phantom>& rhs)
{
    return lhs.value() == rhs.value();
}

//! comparison function <
template<class T, class Phantom>
inline bool operator<(const StrongType<T, Phantom>& lhs, const StrongType<T, Phantom>& rhs)
{
    return lhs.value() < rhs.value();
}

//! comparison function >
template<class T, class Phantom>
inline bool operator>(const StrongType<T, Phantom>& lhs, const StrongType<T, Phantom>& rhs)
{
    return lhs.value() > rhs.value();
}


//! Base template for a holder of entries of different data types
template<class... Ts>
struct TypeList
{
};

//! Base template for mapping between a datatype P templated separately with instances of type list L
template<template<class...> class P, class L>
struct Map_
{
};

//! this is a specialization of the Map_ base template
//! for the case that the L template parameter itself has template parameters
//! in this case, the template parameters of L are caught in Ts...
template<template<class...> class P, template<class...> class L, class... Ts>
struct Map_<P, L<Ts...>>
{
    //! resulting type is a TypeList of the P-template instantiated
    //! with all template parameters of L
    typedef TypeList<P<Ts>...> type;
};

//! Maps a datatype P to create instances where each is templated with entries of type list L
template<template<class...> class P, class L>
using Map = typename Map_<P, L>::type;

//! Base template for expressing a datatype P templated with all the entries in type list L
template<template<class...> class P, class L>
struct Reduce_
{
};

//! Specialization of the Reduce_ base template
template<template<class...> class P, template<class...> class L, class... Ts>
struct Reduce_<P, L<Ts...>>
{
    //! resulting type is P instantiated
    //! with all template parameters of L
    typedef P<Ts...> type;
};

//! Expresses a data type P instantiated with all the entries in list L as template arguments
template<template<class...> class P, class L>
using Reduce = typename Reduce_<P, L>::type;

} // namespace nblib

#endif // NBLIB_UTIL_USER_H
