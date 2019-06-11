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
/*! \file
 * \brief Tools for managing mappings of gmxapi data types.
 *
 * \author M. Eric Irrgang <ericirrgang@gmail.com>
 */

#ifndef GMXPY_TYPETEMPLATES_H
#define GMXPY_TYPETEMPLATES_H


#include <type_traits>

#include "mdparams.h"

namespace gmxapicompat
{

namespace traits
{

// These can be more than traits. We might as well make them named types.
struct gmxNull {
    static const GmxapiType value = GmxapiType::gmxNull;
};
struct gmxMap {
    static const GmxapiType value = GmxapiType::gmxMap;
};
struct gmxInt32 {
    static const GmxapiType value = GmxapiType::gmxInt32;
};
struct gmxInt64 {
    static const GmxapiType value = GmxapiType::gmxInt64;
};
struct gmxFloat32 {
    static const GmxapiType value = GmxapiType::gmxFloat32;
};
struct gmxFloat64 {
    static const GmxapiType value = GmxapiType::gmxFloat64;
};
struct gmxBool {
    static const GmxapiType value = GmxapiType::gmxBool;
};
struct gmxString {
    static const GmxapiType value = GmxapiType::gmxString;
};
struct gmxMDArray {
    static const GmxapiType value = GmxapiType::gmxMDArray;
};
//struct gmxFloat32Vector3 {
//    static const GmxapiType value = GmxapiType::gmxFloat32Vector3;
//};
//struct gmxFloat32SquareMatrix3 {
//    static const GmxapiType value = GmxapiType::gmxFloat32SquareMatrix3;
//};

}   // end namespace traits

// Use an anonymous namespace to restrict these template definitions to file scope.
namespace
{
// Partial specialization of functions is not allowed, which makes the following tedious.
// To-do: switch to type-based logic, struct templates, etc.
template<typename T, size_t s>
GmxapiType mapCppType()
{
    return GmxapiType::gmxNull;
}

template<typename T>
GmxapiType mapCppType()
{
    return mapCppType<T, sizeof(T)>();
};

template<>
GmxapiType mapCppType<bool>()
{
    return GmxapiType::gmxBool;
}

template<>
GmxapiType mapCppType<int, 4>()
{
    return GmxapiType::gmxInt32;
}

template<>
GmxapiType mapCppType<int, 8>()
{
    return GmxapiType::gmxInt64;
};


template<>
GmxapiType mapCppType<float, 4>()
{
    return GmxapiType::gmxFloat32;
}

template<>
GmxapiType mapCppType<double, 8>()
{
    return GmxapiType::gmxFloat64;
};

}      // end anonymous namespace

}      // end namespace gmxapicompat

#endif //GMXPY_TYPETEMPLATES_H
