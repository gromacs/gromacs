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

#include <optional>
#include <string>

#include "h5md_error.h"
#include "h5md_guard.h"
#include "h5md_type.h"

namespace gmx
{
/*! \brief Read an attribute of a given data type
 *
 * \tparam ValueType The type of the attribute to read
 * \param[in] container The path of the group or dataset to read the attribute from.
 * \param[in] attributeName The name of the attribute to read.
 * \returns The value of the desired type, or std::nullopt if the attribute does not exist.
 */
template<typename ValueType>
std::optional<ValueType> getAttribute(const hid_t container, const std::string& attributeName);

/*! \copydoc getAttribute()
 * \brief Specialization of getAttribute() for reading attributes of string type.
 */
template<>
std::optional<std::string> getAttribute<std::string>(const hid_t container, const std::string& attributeName);

/*! \brief Read a vector-like attribute of a given data type
 *
 * \tparam ValueType The type of the attribute to read
 * \param[in] container The path of the group or dataset to read the attribute from.
 * \param[in] attributeName The name of the attribute to read.
 * \returns The 1D vector of the desired data, or std::nullopt if the attribute does not exist.
 */
template<typename ValueType>
std::optional<std::vector<ValueType>> getAttributeVector(const hid_t        container,
                                                         const std::string& attributeName);

/*! \copydoc getAttributeVector()
 * \brief Specialization of getAttributeVector() for reading attributes of string type.
 */
template<>
std::optional<std::vector<std::string>> getAttributeVector<std::string>(const hid_t container,
                                                                        const std::string& attributeName);

/*! \brief Write a scalar attribute of a given data type
 *
 * \tparam ValueType The type of the attribute to write
 * \param[in] container The path of the group or dataset to write the attribute to.
 * \param[in] attributeName The name of the attribute to write.
 * \param[in] value The scalar to write.
 */
template<typename ValueType>
void setAttribute(const hid_t container, const std::string& attributeName, const ValueType& value);

/*! \copydoc setAttribute()
 * \brief Specialization of setAttribute() for writing attributes of char* type.
 */
template<>
void setAttribute<const char*>(const hid_t        container,
                               const std::string& attributeName,
                               const char* const& value);

/*! \copydoc setAttribute()
 * \brief Specialization of setAttribute() for writing attributes of std::string type.
 */
template<>
void setAttribute<std::string>(const hid_t        container,
                               const std::string& attributeName,
                               const std::string& value);

/*! \brief Write a vector-like attribute of a given data type
 *
 * \tparam ValueType The type of the attribute to write.
 * \param[in] container The path of the group or dataset to write the attribute to.
 * \param[in] attributeName The name of the attribute to write.
 * \param[in] value The 1D-vector to write.
 */
template<typename ValueType>
void setAttributeVector(const hid_t                   container,
                        const std::string&            attributeName,
                        const std::vector<ValueType>& value);

/*! \copydoc setAttributeVector()
 * \brief Specialization of setAttributeVector() for writing attributes of char* type.
 */
template<>
void setAttributeVector<const char*>(const hid_t                     container,
                                     const std::string&              attributeName,
                                     const std::vector<const char*>& value);

/*! \copydoc setAttributeVector()
 * \brief Specialization of setAttributeVector() for writing attributes of std::string type.
 */
template<>
void setAttributeVector<std::string>(const hid_t                     container,
                                     const std::string&              attributeName,
                                     const std::vector<std::string>& value);


extern template std::optional<int32_t>  getAttribute<int32_t>(const hid_t        container,
                                                             const std::string& attributeName);
extern template std::optional<int64_t>  getAttribute<int64_t>(const hid_t        container,
                                                             const std::string& attributeName);
extern template std::optional<uint32_t> getAttribute<uint32_t>(const hid_t        container,
                                                               const std::string& attributeName);
extern template std::optional<uint64_t> getAttribute<uint64_t>(const hid_t        container,
                                                               const std::string& attributeName);
extern template std::optional<float>    getAttribute<float>(const hid_t        container,
                                                         const std::string& attributeName);
extern template std::optional<double>   getAttribute<double>(const hid_t        container,
                                                           const std::string& attributeName);


extern template std::optional<std::vector<int32_t>>
getAttributeVector<int32_t>(const hid_t container, const std::string& attributeName);
extern template std::optional<std::vector<int64_t>>
getAttributeVector<int64_t>(const hid_t container, const std::string& attributeName);
extern template std::optional<std::vector<uint32_t>>
getAttributeVector<uint32_t>(const hid_t container, const std::string& attributeName);
extern template std::optional<std::vector<uint64_t>>
getAttributeVector<uint64_t>(const hid_t container, const std::string& attributeName);
extern template std::optional<std::vector<float>>  getAttributeVector<float>(const hid_t container,
                                                                            const std::string& attributeName);
extern template std::optional<std::vector<double>> getAttributeVector<double>(const hid_t container,
                                                                              const std::string& attributeName);


extern template void setAttribute<int32_t>(const hid_t        container,
                                           const std::string& attributeName,
                                           const int32_t&     value);
extern template void setAttribute<int64_t>(const hid_t        container,
                                           const std::string& attributeName,
                                           const int64_t&     value);
extern template void setAttribute<uint32_t>(const hid_t        container,
                                            const std::string& attributeName,
                                            const uint32_t&    value);
extern template void setAttribute<uint64_t>(const hid_t        container,
                                            const std::string& attributeName,
                                            const uint64_t&    value);
extern template void setAttribute<float>(const hid_t        container,
                                         const std::string& attributeName,
                                         const float&       value);
extern template void setAttribute<double>(const hid_t        container,
                                          const std::string& attributeName,
                                          const double&      value);


extern template void setAttributeVector<int32_t>(const hid_t                 container,
                                                 const std::string&          attributeName,
                                                 const std::vector<int32_t>& value);
extern template void setAttributeVector<int64_t>(const hid_t                 container,
                                                 const std::string&          attributeName,
                                                 const std::vector<int64_t>& value);
extern template void setAttributeVector<uint32_t>(const hid_t                  container,
                                                  const std::string&           attributeName,
                                                  const std::vector<uint32_t>& value);
extern template void setAttributeVector<uint64_t>(const hid_t                  container,
                                                  const std::string&           attributeName,
                                                  const std::vector<uint64_t>& value);
extern template void setAttributeVector<float>(const hid_t               container,
                                               const std::string&        attributeName,
                                               const std::vector<float>& value);
extern template void setAttributeVector<double>(const hid_t                container,
                                                const std::string&         attributeName,
                                                const std::vector<double>& value);


} // namespace gmx

#endif // GMX_FILEIO_H5MD_ATTRIBUTE_H
