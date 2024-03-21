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

#include "gmxpre.h"

#include "h5md_datablock.h"

#include <string>

#include "gromacs/math/vectypes.h"
#include "gromacs/utility/exceptions.h"
#include "gromacs/utility/gmxassert.h"
#include "gromacs/utility/smalloc.h"

#include "h5md_util.h"

#define GMX_USE_HDF5 1 // FIXME: Temporary just for the editor

#if GMX_USE_HDF5
#    include <hdf5.h>
#endif

namespace gmx
{
namespace h5mdio
{

GmxH5mdTimeDataBlock::GmxH5mdTimeDataBlock(hid_t                container,
                                           const std::string    name,
                                           const std::string    unit,
                                           hsize_t              numFramesPerChunk,
                                           hsize_t              numEntries,
                                           hsize_t              numValuesPerEntry,
                                           hid_t                datatype,
                                           CompressionAlgorithm compression,
                                           double               compressionAbsoluteError)
{
    container_ = container;
    name_      = name;

    group_ = openOrCreateGroup(container_, name_.c_str());
    char tmpFullName[c_maxFullNameLength];
    H5Iget_name(group_, tmpFullName, c_maxFullNameLength);
    fullName_          = tmpFullName;
    readingFrameIndex_ = 0;
    writingFrameIndex_ = 0;

    static constexpr char c_valueName[] = "value";
    static constexpr char c_stepName[]  = "step";
    static constexpr char c_timeName[]  = "time";

    hsize_t chunkDims[3];
    /* With these default settings new data sets cannot be created. Just load existing from file (if any). */
    if (datatype == -1 && numEntries == 0)
    {
        /* The chunk cache size must be set when opening the data set. We must first open and check the chunk size. */
        hid_t tmpDataSet = H5Dopen(group_, c_valueName, H5P_DEFAULT);

        hid_t createPropertyList = H5Dget_create_plist(tmpDataSet);
        if (H5Pget_chunk(createPropertyList, DIM, chunkDims) < 0)
        {
            H5Eprint2(H5E_DEFAULT, nullptr);
            throw gmx::FileIOError("Error getting chunk size of data set.");
        }
        H5Dclose(tmpDataSet);

        size_t cacheSize = sizeof(real) * 2;
        for (int i = 0; i < DIM; i++)
        {
            cacheSize *= chunkDims[i];
        }
        hid_t accessPropertyList = H5Pcreate(H5P_DATASET_ACCESS);
        if (H5Pset_chunk_cache(accessPropertyList, H5D_CHUNK_CACHE_NSLOTS_DEFAULT, cacheSize, H5D_CHUNK_CACHE_W0_DEFAULT)
            < 0)
        {
            H5Eprint2(H5E_DEFAULT, nullptr);
            throw gmx::FileIOError("Error setting chunk size of data set.");
        }
        mainDataSet_ = H5Dopen(group_, c_valueName, accessPropertyList);
        stepDataSet_ = H5Dopen(group_, c_stepName, H5P_DEFAULT);
        timeDataSet_ = H5Dopen(group_, c_timeName, H5P_DEFAULT);

        updateUnitsFromFile();
    }
    else
    {
        mainUnit_ = unit;
        timeUnit_ = "ps";

        chunkDims[0] = numFramesPerChunk;
        chunkDims[1] = numEntries;
        chunkDims[2] = numValuesPerEntry;

        mainDataSet_ = openOrCreateDataSet<DIM>(
                group_, c_valueName, mainUnit_.c_str(), datatype, chunkDims, compression, compressionAbsoluteError);

        hsize_t chunkDimsTimeStep[1] = { numFramesPerChunk };

        stepDataSet_ = openOrCreateDataSet<1>(
                group_, c_stepName, nullptr, H5T_NATIVE_INT64, chunkDimsTimeStep, CompressionAlgorithm::None, 0);

#if GMX_DOUBLE
        const hid_t timeDatatype = H5Tcopy(H5T_NATIVE_DOUBLE);
#else
        const hid_t timeDatatype = H5Tcopy(H5T_NATIVE_FLOAT);
#endif
        timeDataSet_ = openOrCreateDataSet<1>(
                group_, c_timeName, timeUnit_.c_str(), timeDatatype, chunkDimsTimeStep, CompressionAlgorithm::None, 0);
    }
}

void GmxH5mdTimeDataBlock::closeAllDataSets()
{
    if (mainDataSet_ >= 0)
    {
        H5Dclose(mainDataSet_);
        mainDataSet_ = -1;
    }
    if (timeDataSet_ >= 0)
    {
        H5Dclose(timeDataSet_);
        timeDataSet_ = -1;
    }
    if (stepDataSet_ >= 0)
    {
        H5Dclose(stepDataSet_);
        stepDataSet_ = -1;
    }
}

bool GmxH5mdTimeDataBlock::operator==(const std::string fullNameComparison)
{
    if (fullName_ == fullNameComparison)
    {
        return true;
    }
    return false;
}

void GmxH5mdTimeDataBlock::writeFrame(const void* data, int64_t step, real time, int64_t frame)
{
    if (frame < 0)
    {
        frame = writingFrameIndex_;
    }
    GMX_ASSERT(frame >= 0, "Invalid frame when writing.");

    writeData<3, false>(mainDataSet_, data, frame);
    writeData<1, false>(stepDataSet_, &step, frame);
    writeData<1, false>(timeDataSet_, &time, frame);

    writingFrameIndex_ = frame + 1;
}

bool GmxH5mdTimeDataBlock::readFrame(real* data, int64_t frame)
{
    size_t totalNumElements;
    size_t dataTypeSize = getDataTypeSize(mainDataSet_);

    /* FIXME: Make the number of dimensions flexible */

#if GMX_DOUBLE
    if (dataTypeSize != 8)
    {
        void* buffer = nullptr;
        readData<3, false>(mainDataSet_, frame, dataTypeSize, &buffer, &totalNumElements);
        for (size_t i = 0; i < totalNumElements; i++)
        {
            data[i] = static_cast<float*>(buffer)[i];
        }
        H5free_memory(buffer);
    }
    else
    {
        readData<3, false>(
                mainDataSet_, frame, dataTypeSize, reinterpret_cast<void**>(&data), &totalNumElements);
    }
#else
    if (dataTypeSize != 4)
    {
        void* buffer = nullptr;
        readData<3, false>(mainDataSet_, frame, dataTypeSize, &buffer, &totalNumElements);
        for (size_t i = 0; i < totalNumElements; i++)
        {
            data[i] = static_cast<double*>(buffer)[i];
        }
        H5free_memory(buffer);
    }
    else
    {
        readData<3, false>(
                mainDataSet_, frame, dataTypeSize, reinterpret_cast<void**>(&data), &totalNumElements);
    }
#endif

    return true;
}

bool GmxH5mdTimeDataBlock::readNextFrame(real* data)
{
    if (readingFrameIndex_ >= writingFrameIndex_)
    {
        return false;
    }
    return readFrame(data, readingFrameIndex_++);
}

void GmxH5mdTimeDataBlock::updateUnitsFromFile()
{
    char* tmpStr = nullptr;
    if (getAttribute(mainDataSet_, "unit", &tmpStr))
    {
        mainUnit_ = tmpStr;
        H5free_memory(tmpStr);
        tmpStr = nullptr;
    }
    if (getAttribute(timeDataSet_, "unit", &tmpStr))
    {
        timeUnit_ = tmpStr;
        H5free_memory(tmpStr);
        tmpStr = nullptr;
    }
}

void GmxH5mdTimeDataBlock::updateNumWrittenFrames()
{
    GMX_ASSERT(stepDataSet_ >= 0,
               "There must be a data set with steps to determine the actual number "
               "of frames.");
    hid_t     dataSpace = H5Dget_space(stepDataSet_);
    const int numDims   = H5Sget_simple_extent_ndims(dataSpace);
    if (numDims != 1)
    {
        H5Eprint2(H5E_DEFAULT, nullptr);
        throw gmx::FileIOError("The step data set should be one-dimensional.");
    }
    hsize_t dimExtents;
    H5Sget_simple_extent_dims(dataSpace, &dimExtents, nullptr);
    hsize_t numValidFrames = dimExtents;

    hid_t            createPropertyList = H5Dget_create_plist(stepDataSet_);
    H5D_fill_value_t fillValueStatus;
    if (H5Pfill_value_defined(createPropertyList, &fillValueStatus) < 0
        || fillValueStatus == H5D_FILL_VALUE_UNDEFINED)
    {
        writingFrameIndex_ = numValidFrames;
        return;
    }
    hid_t   datatype = H5T_NATIVE_INT64;
    int64_t fillValue;
    H5Pget_fill_value(createPropertyList, datatype, &fillValue);

    int64_t stepData = fillValue;
    while (numValidFrames > 0 && stepData == fillValue)
    {
        numValidFrames--;
        hsize_t location = numValidFrames;
        H5Sselect_elements(dataSpace, H5S_SELECT_SET, 1, &location);
        const hsize_t memorySpaceSize = 1;
        hid_t         memSpace        = H5Screate_simple(1, &memorySpaceSize, nullptr);
        if (H5Dread(stepDataSet_, datatype, memSpace, dataSpace, H5P_DEFAULT, &stepData) < 0)
        {
            H5Eprint2(H5E_DEFAULT, nullptr);
            throw gmx::FileIOError(
                    "Error reading step data set when determining the number of frames.");
        }
    }
    writingFrameIndex_ = numValidFrames + 1;
}

size_t GmxH5mdTimeDataBlock::getNumParticles() const
{
    hid_t dataSpace = H5Dget_space(mainDataSet_);
    if (dataSpace < 0)
    {
        H5Eprint2(H5E_DEFAULT, nullptr);
        throw gmx::FileIOError(
                "The main data block of the time dependent data set cannot be found.");
    }
    const int dataSpaceNumDims = H5Sget_simple_extent_ndims(dataSpace);

    if (dataSpaceNumDims < 2)
    {
        throw gmx::FileIOError("No atoms in time dependent data set.");
    }
    hsize_t* dimExtents;
    snew(dimExtents, dataSpaceNumDims);
    H5Sget_simple_extent_dims(dataSpace, dimExtents, nullptr);

    return dimExtents[1];
}

int64_t GmxH5mdTimeDataBlock::getStepOfFrame(int64_t frame) const
{
    GMX_ASSERT(stepDataSet_ >= 0, "There must be a data set with steps to get the step of a frame.");

    size_t totalNumElements;
    size_t dataTypeSize = getDataTypeSize(stepDataSet_);

    if (dataTypeSize != 4 && dataTypeSize != 8)
    {
        throw gmx::FileIOError("Can only read 32- or 64-bit step data.");
    }

    if (dataTypeSize == 4)
    {
        void* buffer = nullptr;
        readData<1, false>(stepDataSet_, frame, dataTypeSize, &buffer, &totalNumElements);
        int* tmpIntData = static_cast<int*>(buffer);
        int  tmpValue   = *tmpIntData;
        H5free_memory(tmpIntData);
        tmpIntData = nullptr;
        return tmpValue;
    }
    else
    {
        int64_t  tmpValue;
        int64_t* tmpValuePtr = &tmpValue;
        readData<1, false>(
                stepDataSet_, frame, dataTypeSize, reinterpret_cast<void**>(&tmpValuePtr), &totalNumElements);
        return tmpValue;
    }
}

real GmxH5mdTimeDataBlock::getTimeOfFrame(int64_t frame) const
{
    GMX_ASSERT(timeDataSet_ >= 0, "There must be a data set with time to get the time of a frame.");

    size_t totalNumElements;
    size_t dataTypeSize = getDataTypeSize(timeDataSet_);

    if (dataTypeSize != 4 && dataTypeSize != 8)
    {
        throw gmx::FileIOError("Can only read float or double time data.");
    }

#if GMX_DOUBLE
    if (dataTypeSize != 8)
    {
        void* buffer = nullptr;
        readData<1, false>(timeDataSet_, frame, dataTypeSize, &buffer, &totalNumElements);
        float* tmpFloatData = static_cast<float*>(buffer);
        double tmpValue     = *tmpFloatData;
        H5free_memory(tmpFloatData);
        tmpFloatData = nullptr;
        return tmpValue;
    }
    else
    {
        double  tmpValue;
        double* tmpValuePtr = &tmpValue;
        readData<1, false>(
                timeDataSet_, frame, dataTypeSize, reinterpret_cast<void**>(&tmpValuePtr), &totalNumElements);
        return tmpValue;
    }
#else
    if (dataTypeSize != 4)
    {
        void* buffer = nullptr;
        readData<1, false>(timeDataSet_, frame, dataTypeSize, &buffer, &totalNumElements);
        double* tmpDoubleData = static_cast<double*>(buffer);
        float tmpValue = *tmpDoubleData;
        H5free_memory(tmpDoubleData);
        tmpDoubleData = nullptr;
        return tmpValue;
    }
    else
    {
        float tmpValue;
        float* tmpValuePtr = &tmpValue;
        readData<1, false>(
                timeDataSet_, frame, dataTypeSize, reinterpret_cast<void**>(&tmpValuePtr), &totalNumElements);
        return tmpValue;
    }
#endif
}

real GmxH5mdTimeDataBlock::getLossyCompressionError()
{
    return getDataSetSz3CompressionError(mainDataSet_);
}

extern template hid_t
openOrCreateDataSet<1>(hid_t, const char*, const char*, hid_t, const hsize_t*, CompressionAlgorithm, double);
extern template hid_t
openOrCreateDataSet<3>(hid_t, const char*, const char*, hid_t, const hsize_t*, CompressionAlgorithm, double);

extern template void writeData<1, false>(hid_t, const void*, hsize_t);
extern template void writeData<3, false>(hid_t, const void*, hsize_t);

extern template void readData<1, false>(hid_t, hsize_t, size_t, void**, size_t*);
extern template void readData<3, false>(hid_t, hsize_t, size_t, void**, size_t*);

} // namespace h5mdio
} // namespace gmx
