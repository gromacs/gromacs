/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2025- The GROMACS Authors
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

/*! \brief Definitions of H5md data set manipulation routines.
 *
 * \author Magnus Lundborg <lundborg.magnus@gmail.com>
 * \author Petter Johansson <pettjoha@kth.se>
 */

#include "gmxpre.h"

#include "h5md_dataset.h"

#include "gromacs/utility/basedefinitions.h"

#include "h5md_error.h"
#include "h5md_guard.h"
#include "h5md_type.h"

// HDF5 constants use old style casts.
CLANG_DIAGNOSTIC_IGNORE("-Wold-style-cast")

namespace
{

using namespace gmx;

/*! \brief Set the fill value to -1 in a data set creation property list.
 * \tparam ValueType The compiled type of data to set fill value for.
 * \param[in] createPropertyList The ID of the propery list to update.
 * \param[in] dataType The ID of the HDF5 data type of the data set.
 */
template<typename ValueType>
static void setNumericFillValue(const hid_t createPropertyList, const hid_t dataType)
{
    throwUponH5mdError(!valueTypeIsDataType<ValueType>(dataType),
                       "ValueType != dataType mismatch when setting fill value for data set");

    constexpr ValueType fillValue = -1;
    H5Pset_fill_value(createPropertyList, dataType, &fillValue);
}

/*! \brief Return the dimensions of a data set.
 *
 * \tparam numDims Number of dimensions of data set.
 * \param[in] dataSet Handle to data set.
 * \returns Fixed-size array with the dimensions.
 *
 * \throws gmx::FileIOError if the dimensions cannot be read.
 */
template<int numDims>
static DataSetDims<numDims> getDataSetDims(const hid_t dataSet)
{
    DataSetDims<numDims> dataSetDims;

    const auto [dataSpace, dataSpaceGuard] = makeH5mdDataSpaceGuard(H5Dget_space(dataSet));
    throwUponH5mdError(H5Sget_simple_extent_dims(dataSpace, dataSetDims.data(), nullptr) < 0,
                       "Could not get dimensions from data set.");

    return dataSetDims;
}

/*! \brief Return the per-dimension offset to the data corresponding to a given frame.
 *
 * Since frames are always the major axis the first value of the returned offset is
 * the given frame index. We are only considering the case when all data for that frame
 * is selected, so the offset for all other axis is 0. For example, to reach frame 15
 * of a 3d set we use the offset [15, 0, 0] to access the memory.
 *
 * \tparam numDims Number of dimensions of data set.
 * \param[in] frameIndex Index of frame.
 * \returns Fixed-size array with the offset.
 */
template<int numDims>
static DataSetDims<numDims> getFrameChunkOffset(const hsize_t frameIndex)
{
    return DataSetDims<numDims>{ frameIndex };
}

/* \brief Return the dimensions of the memory chunk for a single frame.
 *
 * Since frames are always the major axis the first value off the chunk dimensions is 1.
 * We only consider the case of selecting all data corresponding to a single frame so
 * all other chunk dimensions are equal to the data set dimensions. For example, the
 * chunk which contains a single frame in a data set with dimensions [30, 150, 3] is
 * [1, 150, 3].
 *
 * \tparam numDims Number of dimensions of data set.
 * \param[in] dataSetDims Dimensions of data set.
 * \returns Fixed-size array with the chunk dimensions.
 */
template<int numDims>
static DataSetDims<numDims> getFrameChunkDims(const DataSetDims<numDims>& dataSetDims)
{
    DataSetDims<numDims> chunkDims = dataSetDims;
    chunkDims[0]                   = 1;

    return chunkDims;
}

/*! \brief Create a data set and return its handle.
 *
 * The returned handle must be closed with H5Dclose to avoid resource leaks.
 *
 * \tparam ValueType The compiled type of data to create data set for.
 * \tparam numDims Number of dimensions of created data set.
 * \param[in] container Container in which to create the data set.
 * \param[in] dataSetName Name of data set.
 * \param[in] dataType Handle to HDF5 data type for the data set values.
 * \param[in] dataSetDims Dimensions to create the data set with.
 * \param[in] dataSetMaxDims Maximum dimensions of the data set.
 * \returns Handle to created data set.
 */
template<typename ValueType, int numDims>
static hid_t createDataSet(const hid_t                 container,
                           const std::string&          dataSetName,
                           const hid_t                 dataType,
                           const DataSetDims<numDims>& dataSetDims,
                           const DataSetDims<numDims>& dataSetMaxDims)
{
    // We need to set some options for the new data set, using this property list
    const auto [creationPropertyList, creationPropertyListGuard] =
            makeH5mdPropertyListGuard(H5Pcreate(H5P_DATASET_CREATE));

    // Memory chunk sizes must be > 0 along every dimension when reading or writing
    const DataSetDims<numDims> chunkDims = [&]()
    {
        DataSetDims<numDims> chunkDims = dataSetDims;
        for (auto& v : chunkDims)
        {
            if (v < 1)
            {
                v = 1;
            }
        }
        return chunkDims;
    }();

    H5Pset_chunk(creationPropertyList, numDims, chunkDims.data());
    setNumericFillValue<ValueType>(creationPropertyList, dataType);

    /* It would be nice to have an option not to write full incomplete edge chunks,
     * but the closest option is:
     * H5Pset_chunk_opts(createPropertyList, H5D_CHUNK_DONT_FILTER_PARTIAL_CHUNKS);
     * But that only avoids compressing/decompressing the edge chunks.
     * Keep an eye open for alternatives.
     * Refs #5286: It is possible that it would be time-efficient
     * to avoid compressing edge chunks when writing checkpoints. Pros and cons for slightly
     * larger files vs slightly faster checkpoint writing must be evaluated.
     * Currently it seems like incomplete edge chunks are compressed even with this option.
     */
    H5Pset_chunk_opts(creationPropertyList, H5D_CHUNK_DONT_FILTER_PARTIAL_CHUNKS);

    /* Set a reasonable cache (in bytes) based on chunk sizes. The cache is not stored in file,
     * so must be set when opening a data set / Magnus
     * Refs #5286
     */
    size_t cacheSize = H5Tget_size(dataType);
    for (const auto& chunkSize : chunkDims)
    {
        cacheSize *= chunkSize;
    }

    const auto [accessPropertyList, accessPropertyListGuard] =
            makeH5mdPropertyListGuard(H5Pcreate(H5P_DATASET_ACCESS));
    throwUponH5mdError(
            H5Pset_chunk_cache(accessPropertyList, H5D_CHUNK_CACHE_NSLOTS_DEFAULT, cacheSize, H5D_CHUNK_CACHE_W0_DEFAULT)
                    < 0,
            "Cannot set chunk cache parameters.");

    const auto [dataSpace, dataSpaceGuard] =
            makeH5mdDataSpaceGuard(H5Screate_simple(numDims, dataSetDims.data(), dataSetMaxDims.data()));

    const hid_t dataSet = H5Dcreate(
            container, dataSetName.c_str(), dataType, dataSpace, H5P_DEFAULT, creationPropertyList, accessPropertyList);
    throwUponInvalidHid(dataSet, "Cannot create dataSet.");

    return dataSet;
}

} // namespace

namespace gmx
{

template<typename ValueType>
hid_t create1dFrameDataSet(const hid_t container, const std::string& dataSetName)
{
    constexpr int numDims = 1;

    const DataSetDims<numDims> dataSetDims    = { 0 };
    const DataSetDims<numDims> dataSetMaxDims = { H5S_UNLIMITED };
    const hid_t                dataType       = hdf5DataTypeFor<ValueType>();

    return createDataSet<ValueType, numDims>(container, dataSetName, dataType, dataSetDims, dataSetMaxDims);
}

hid_t openDataSet(const hid_t container, const std::string& dataSetName)
{
    const hid_t dataSet = H5Dopen(container, dataSetName.c_str(), H5P_DEFAULT);
    throwUponInvalidHid(dataSet, "Could not open data set.");

    return dataSet;
}

template<int numDims>
hsize_t getNumFrames(const hid_t dataSet)
{
    const DataSetDims<numDims> dataSetDims = getDataSetDims<numDims>(dataSet);

    return dataSetDims.at(0);
}

template<int numDims>
hid_t getFrameDataSpace(const hid_t dataSet, const hsize_t frameIndex)
{
    const DataSetDims<numDims> dataSetDims = getDataSetDims<numDims>(dataSet);
    const DataSetDims<numDims> chunkOffset = getFrameChunkOffset<numDims>(frameIndex);
    const DataSetDims<numDims> chunkDims   = getFrameChunkDims<numDims>(dataSetDims);

    const hid_t frameDataSpace = H5Dget_space(dataSet);
    throwUponH5mdError(
            H5Sselect_hyperslab(
                    frameDataSpace, H5S_SELECT_SET, chunkOffset.data(), nullptr, chunkDims.data(), nullptr)
                    < 0,
            "Could not select hyperslab for given frame index.");

    return frameDataSpace;
}

template<int numDims>
hid_t getFrameMemoryDataSpace(const hid_t dataSet)
{
    const DataSetDims<numDims> dataSetDims = getDataSetDims<numDims>(dataSet);
    const DataSetDims<numDims> chunkDims   = getFrameChunkDims<numDims>(dataSetDims);

    return H5Screate_simple(numDims, chunkDims.data(), nullptr);
}

template hid_t create1dFrameDataSet<int32_t>(const hid_t, const std::string&);

template hid_t create1dFrameDataSet<int64_t>(const hid_t, const std::string&);

template hid_t create1dFrameDataSet<float>(const hid_t, const std::string&);

template hid_t create1dFrameDataSet<double>(const hid_t, const std::string&);

template hsize_t getNumFrames<1>(const hid_t dataSet);

template hid_t getFrameDataSpace<1>(const hid_t dataSet, const hsize_t frameIndex);

template hid_t getFrameMemoryDataSpace<1>(const hid_t);

} // namespace gmx

CLANG_DIAGNOSTIC_RESET
