/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2023- The GROMACS Authors
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

#ifndef GMX_FILEIO_H5MD_UTIL_H
#define GMX_FILEIO_H5MD_UTIL_H

#include <string>

#include "gromacs/utility/real.h"

enum class PbcType : int;
namespace gmx
{
namespace h5mdio
{
typedef int64_t            hid_t;
typedef unsigned long long hsize_t;

/*! \brief An enumeration of compression options */
enum class CompressionAlgorithm
{
    None,                //!< No compression
    LosslessNoShuffle,   //!< Lossless, only gzip (deflate)
    LosslessWithShuffle, //!< Lossless, byte shuffle followed by gzip (deflate)
    LossySz3,            //!< Lossy SZ3 compression.
    Count
};

/*! \brief Open an existing HDF5 group or create it if it did not exist already.
 *
 * \param[in] container  The container where the group is located, or should be created.
 * \param[in] name       The name of the group.
 * \returns the ID of the group.
 */
hid_t openOrCreateGroup(hid_t container, const char* name);

/*! \brief Registers the SZ3 filter by using the automatic registration mechanism by H5Pset_filter().
 * Must be done before appending (e.g. when restarting from a checkpoint) to a compressed dataset. */
void registerSz3FilterImplicitly();

/*! \brief Open an existing dataset (called name, in container). If it does not exist create a new dataset.
 *
 * \tparam numDims The number of dimensions.
 * \param[in] container The ID of the container of the data. This can be a group in the HDF5 or the HDF5 file itself.
 * \param[in] name The name of the data set.
 * \param[in] unit The unit of the data. See de Buyl et al., 2014 (https://www.sciencedirect.com/science/article/pii/S0010465514000447) for more information.
 * \param[in] datatype The HDF5 data type of the data.
 * \param[in] chunkDims The dimensions of each chunk. The size should match numDims.
 * \param[in] compression The compression algorithm to use.
 * \param[in] compressionError The required precision of lossy compression.
 * \returns The ID of the dataset.
 */

template<int numDims>
hid_t openOrCreateDataSet(hid_t                container,
                          const char*          name,
                          const char*          unit,
                          hid_t                dataType,
                          const hsize_t*       chunkDims,
                          CompressionAlgorithm compression,
                          double               compressionError);

/*! \brief Writes data to an HDF5 data set, labelled by name.
 *
 * \tparam numDims The number of dimensions of the data.
 * \tparam writeFullDataSet Whether to write the whole data set at once. If true, frameToWrite must be 0.
 * \param[in] dataSet The ID of the dataset to write to.
 * \param[in] data The data to write.
 * \param[in] frameToWrite The frame number to write (starting from 0).
 */
template<int numDims, bool writeFullDataSet>
void writeData(hid_t dataSet, const void* data, hsize_t frameToWrite);

/*! \brief Returns the size of the data type used in a dataset.
 * \param[in] The ID of the dataset.
 * \returns the size of the data type used in the dataset.
 */
size_t getDataTypeSize(hid_t dataSet);

/*! \brief Reads data from an HDF5 data set, labelled by name.
 *
 * \tparam numDims The number of dimensions of the data.
 * \tparam readFullDataSet Whether to read the whole data set at once.
 *                         If true, frameToRead must be 0.
 * \param[in] dataSet The ID of the dataset to read from.
 * \param[in] frameToRead The frame number to read (starting from 0).
 * \param[in] dataTypeSize The size (in bytes) per value.
 * \param[out] buffer The buffer to fill with the read data. Memory will be allocated.
 *                    Must be freed by the caller.
 * \param[out] totalNumElements The number of data values read.
 */
template<int numDims, bool readFullDataSet>
void readData(hid_t dataSet, hsize_t frameToRead, size_t dataTypeSize, void** buffer, size_t* totalNumElements);

/*! Set the periodic box attribute of a box group.
 * \param[in] boxGroup The ID of the box group.
 * \param[in] pbcType The periodic box type.
 */
void setBoxGroupAttributes(hid_t boxGroup, PbcType pbcType);

/*! Set an attribute value in a data set.
 * \tparam T The type of the data to write.
 * \param[in] dataSet The ID of the HDF5 data set.
 * \param[in] name The name of the attribute.
 * \param[in] value The new value of the attribute.
 * \param[in] dataType The HDF5 type of the output.
 */
template<typename T>
void setAttribute(hid_t dataSet, const char* name, const T value, hid_t dataType);

/*! Set a string attribute value in a data set.
 * \param[in] dataSet The ID of the HDF5 data set.
 * \param[in] name The name of the attribute.
 * \param[in] value The string to set as attribute value.
 */
void setAttribute(hid_t dataSet, const char* name, const char* value);

/*! Get an attribute value from a data set.
 * \tparam T The type of the data to read.
 * \param[in] dataSet The ID of the HDF5 data set.
 * \param[in] name The name of the attribute.
 * \param[out] value The returned value of the attribute.
 * \param[in] dataType The HDF5 type of the output.
 */
template<typename T>
bool getAttribute(hid_t dataSet, const char* name, T* value, hid_t dataType);

/*! Get a string attribute value from a data set.
 * \param[in] dataSet The ID of the HDF5 data set.
 * \param[in] name The name of the attribute. Memory must be allocated by the caller.
 * \param[out] value The returned string value of the attribute.
 */
bool getAttribute(hid_t dataSet, const char* name, char** value);

/*! Set a list of attribute strings in a data set.
 * \tparam numEntries The number of strings to set.
 * \tparam stringLength The length of the strings. Must be the same for all.
 * \param[in] dataSet The ID of the HDF5 data set.
 * \param[in] name The name of the attribute.
 * \param[in] value The list of string to set as attribute value.
 */
template<hid_t numEntries, size_t stringLength>
void setAttributeStringList(hid_t dataSet, const char* name, const char value[numEntries][stringLength]);

/*! Get the SZ3 lossy compression error setting (absolute or relative) from a data set.
 * \param[in] dataSet The ID of the HDF5 data set.
 * \returns The compression error setting of the data set or -1 if it is not SZ3 compressed or if the error is not absolute or relative.
 */
real getDataSetSz3CompressionError(hid_t dataSet);

/*! Check if an object exists in a container
 * \param[in] container The ID of the HDF5 container.
 * \param[in] name The name of the HDF5 object to look for.
 */
bool objectExists(hid_t container, const char* name);

} // namespace h5mdio
} // namespace gmx
#endif
