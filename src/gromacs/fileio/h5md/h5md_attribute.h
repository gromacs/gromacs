/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2024- The GROMACS Authors
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

/*! \brief Declarations of utility functions for working with H5MD attributes.
 *
 * \author Magnus Lundborg <lundborg.magnus@gmail.com>
 * \author Petter Johansson <pettjoha@kth.se>
 * \author Yang Zhang <yang.zhang@scilifelab.se>
 */

#ifndef GMX_FILEIO_H5MD_ATTRIBUTE_H
#define GMX_FILEIO_H5MD_ATTRIBUTE_H

#include <hdf5.h>

#include <cstring>

#include <optional>
#include <string>
#include <vector>

#include "gromacs/utility/arrayref.h"

namespace gmx
{

constexpr char c_unitAttributeKey[] = "unit";

/*! \brief Read an attribute of a given data type
 *
 * \tparam ValueType The type of the attribute to read
 * \param[in] container The HDF5 container (group or dataset) to read the attribute from.
 * \param[in] attributeName The name of the attribute to read.
 * \returns The value of the desired type, or std::nullopt if the attribute does not exist.
 */
template<typename ValueType>
std::optional<ValueType> getAttribute(const hid_t container, const char* attributeName);

/*! \copydoc getAttribute()
 * \brief Specialization of getAttribute() for reading attributes of string type.
 */
template<>
std::optional<std::string> getAttribute<std::string>(const hid_t container, const char* attributeName);

/*! \brief Read a vector-like attribute of a given data type
 *
 * \tparam ValueType The type of the attribute to read
 * \param[in] container The HDF5 container (group or dataset) to read the attribute from.
 * \param[in] attributeName The name of the attribute to read.
 * \returns The 1D vector of the desired data, or std::nullopt if the attribute does not exist.
 */
template<typename ValueType>
std::optional<std::vector<ValueType>> getAttributeVector(const hid_t container, const char* attributeName);

/*! \copydoc getAttributeVector()
 * \brief Specialization of getAttributeVector() for reading attributes of string type.
 */
template<>
std::optional<std::vector<std::string>> getAttributeVector<std::string>(const hid_t container,
                                                                        const char* attributeName);

/*! \brief Write a scalar attribute of a given data type
 *
 * \tparam ValueType The type of the attribute to write
 * \param[in] container The HDF5 container (group or dataset) to write the attribute to.
 * \param[in] attributeName The name of the attribute to write.
 * \param[in] value The scalar to write.
 */
template<typename ValueType>
void setAttribute(const hid_t container, const char* attributeName, const ValueType& value);

/*! \copydoc setAttribute()
 * \brief Specialization of setAttribute() for writing attributes of char* type.
 */
void setAttribute(const hid_t container, const char* attributeName, const char* value);

/*! \copydoc setAttribute()
 * \brief Specialization of setAttribute() for writing attributes of std::string type.
 */
void setAttribute(const hid_t container, const char* attributeName, const std::string& value);

/*! \brief Write a vector-like attribute of a given data type
 *
 * \tparam    ValueType The type of the attribute to write.
 * \param[in] container The path of the group or dataset to write the attribute to.
 * \param[in] attributeName The name of the attribute to write.
 * \param[in] value The 1D-vector to write, provided as an ArrayRef.
 */
template<typename ValueType>
void setAttributeVector(const hid_t container, const char* attributeName, ArrayRef<const ValueType> value);

/*! \brief Write a vector-like string attribute via iterators
 * String data type is treated specially because of the HDF5 requirements for
 * contiguous memory storage for fixed-size strings. A reusable character buffer
 * is needed to avoid frequent memory allocations for multiple writing of different
 * string attributes.
 *
 * \tparam    Iterator      The type of the iterator to write.
 * \param[in] container     The HDF5 container (group or dataset) to write the attribute to.
 * \param[in] attributeName The name of the attribute to write.
 * \param[in] buffer        The buffer to write the strings into, could leave it empty {}.
 * \param[in] begin         The beginning iterator of the data to write.
 * \param[in] end           The ending iterator of the data to write.
 */
template<typename Iterator>
std::vector<char> setAttributeStringVector(const hid_t         container,
                                           const char*         attributeName,
                                           std::vector<char>&& buffer,
                                           Iterator            begin,
                                           Iterator            end);

/*! \brief Low-level utility function to write a string-vector attribute via a character buffer.
 *
 * The size of the character buffer must be greater than or equal to
 * numberOfStrings * (maxStrLength + 1).
 *
 * \param[in] container       The HDF5 container to write the attribute to.
 * \param[in] attributeName   The name of the attribute to write.
 * \param[in] numberOfStrings The number of strings in the vector.
 * \param[in] maxStrLength    The maximum length of the strings in the vector.
 * \param[in] buffer          The character buffer to write.
 */
void setStringAttributeByBuffer(const hid_t          container,
                                const char*          attributeName,
                                const size_t         numberOfStrings,
                                const int            maxStrLength,
                                ArrayRef<const char> buffer);

/// @cond DO_NOT_DOCUMENT
extern template std::optional<int32_t>  getAttribute(const hid_t, const char*);
extern template std::optional<int64_t>  getAttribute(const hid_t, const char*);
extern template std::optional<uint32_t> getAttribute(const hid_t, const char*);
extern template std::optional<uint64_t> getAttribute(const hid_t, const char*);
extern template std::optional<float>    getAttribute(const hid_t, const char*);
extern template std::optional<double>   getAttribute(const hid_t, const char*);

extern template void setAttribute(const hid_t, const char*, const int32_t&);
extern template void setAttribute(const hid_t, const char*, const int64_t&);
extern template void setAttribute(const hid_t, const char*, const uint32_t&);
extern template void setAttribute(const hid_t, const char*, const uint64_t&);
extern template void setAttribute(const hid_t, const char*, const float&);
extern template void setAttribute(const hid_t, const char*, const double&);

extern template std::optional<std::vector<int32_t>>  getAttributeVector(const hid_t, const char*);
extern template std::optional<std::vector<int64_t>>  getAttributeVector(const hid_t, const char*);
extern template std::optional<std::vector<uint32_t>> getAttributeVector(const hid_t, const char*);
extern template std::optional<std::vector<uint64_t>> getAttributeVector(const hid_t, const char*);
extern template std::optional<std::vector<float>>    getAttributeVector(const hid_t, const char*);
extern template std::optional<std::vector<double>>   getAttributeVector(const hid_t, const char*);

extern template void setAttributeVector(const hid_t, const char*, ArrayRef<const int32_t>);
extern template void setAttributeVector(const hid_t, const char*, ArrayRef<const int64_t>);
extern template void setAttributeVector(const hid_t, const char*, ArrayRef<const uint32_t>);
extern template void setAttributeVector(const hid_t, const char*, ArrayRef<const uint64_t>);
extern template void setAttributeVector(const hid_t, const char*, ArrayRef<const float>);
extern template void setAttributeVector(const hid_t, const char*, ArrayRef<const double>);

extern template std::vector<char> setAttributeStringVector(const hid_t,
                                                           const char*,
                                                           std::vector<char>&&,
                                                           ArrayRef<const std::string>::const_iterator,
                                                           ArrayRef<const std::string>::const_iterator);
extern template std::vector<char> setAttributeStringVector(const hid_t,
                                                           const char*,
                                                           std::vector<char>&&,
                                                           ArrayRef<const char*>::const_iterator,
                                                           ArrayRef<const char*>::const_iterator);
/// @endcond

} // namespace gmx

#endif // GMX_FILEIO_H5MD_ATTRIBUTE_H
