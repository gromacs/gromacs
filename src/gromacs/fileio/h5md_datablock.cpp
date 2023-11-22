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
/* This file was inspired by ch5md by Pierre de Buyl (BSD license). */

#include "gmxpre.h"

#include "h5md_datablock.h"

#include "config.h"

#include <string>

#include <sys/_types/_int64_t.h>

#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/gmxassert.h"
#include "gromacs/utility/smalloc.h"

#include "h5md_util.h"

#define GMX_USE_HDF5 1 // FIXME: Temporary just for the editor

#if GMX_USE_HDF5
#    include <hdf5.h>
#endif


GmxH5mdTimeDataBlock::GmxH5mdTimeDataBlock(hid_t                container,
                                           const std::string    name,
                                           const std::string    unit,
                                           int                  writingInterval,
                                           hsize_t              numFramesPerChunk,
                                           hsize_t              numEntries,
                                           hsize_t              numValuesPerEntry,
                                           hid_t                datatype,
                                           CompressionAlgorithm compression,
                                           double               compressionAbsoluteError)
{
    container_       = container;
    name_            = name;
    writingInterval_ = writingInterval;

    group_ = openOrCreateGroup(container_, name_.c_str());
    char tmpFullName[c_maxFullNameLength];
    H5Iget_name(group_, tmpFullName, c_maxFullNameLength);
    fullName_ = tmpFullName;

    static constexpr char c_valueName[] = "value";
    static constexpr char c_stepName[]  = "step";
    static constexpr char c_timeName[]  = "time";

    /* With these default settings new data sets cannot be created. Just load existing from file (if any). */
    if (datatype == -1 && numEntries == 0)
    {
        mainDataSet_ = H5Dopen(group_, c_valueName, H5P_DEFAULT);
        stepDataSet_ = H5Dopen(group_, c_stepName, H5P_DEFAULT);
        timeDataSet_ = H5Dopen(group_, c_timeName, H5P_DEFAULT);
    }
    else
    {
        hsize_t chunkDims[3] = { numFramesPerChunk, numEntries, numValuesPerEntry };

        mainDataSet_ = openOrCreateDataSet<3>(
                group_, c_valueName, unit.c_str(), datatype, chunkDims, compression, compressionAbsoluteError);

        hsize_t chunkDimsTimeStep[1] = { numFramesPerChunk };
        stepDataSet_                 = openOrCreateDataSet<1>(
                group_, c_stepName, nullptr, H5T_NATIVE_INT64, chunkDimsTimeStep, CompressionAlgorithm::None, 0);

#if GMX_DOUBLE
        const hid_t timeDatatype = H5Tcopy(H5T_NATIVE_DOUBLE);
#else
        const hid_t timeDatatype = H5Tcopy(H5T_NATIVE_FLOAT);
#endif
        static constexpr char c_timeUnit[] = "ps";
        timeDataSet_                       = openOrCreateDataSet<1>(
                group_, c_timeName, c_timeUnit, timeDatatype, chunkDimsTimeStep, CompressionAlgorithm::None, 0);
    }
    updateNumWrittenFrames();
}

void GmxH5mdTimeDataBlock::closeAllDataSets()
{
    if (mainDataSet_ >= 0)
    {
        H5Dclose(mainDataSet_);
    }
    if (timeDataSet_ >= 0)
    {
        H5Dclose(timeDataSet_);
    }
    if (stepDataSet_ >= 0)
    {
        H5Dclose(stepDataSet_);
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

void GmxH5mdTimeDataBlock::writeFrame(const void* data, int64_t step, real time)
{
    GMX_ASSERT(step >= 0, "Invalid step when writing frame.");

    /* If there is no specified writing interval for this data block, write after the previous output. */
    const int frameNumber = writingInterval_ > 0 ? step / writingInterval_ : numWrittenFrames_;

    writeData<3, false>(mainDataSet_, data, frameNumber);
    writeData<1, false>(stepDataSet_, &step, frameNumber);
    writeData<1, false>(timeDataSet_, &time, frameNumber);
    ++numWrittenFrames_;
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
        gmx_file("The step data set should be one-dimensional.");
    }
    hsize_t dimExtents;
    H5Sget_simple_extent_dims(dataSpace, &dimExtents, nullptr);
    hsize_t numValidFrames = dimExtents;

    hid_t            createPropertyList = H5Dget_create_plist(stepDataSet_);
    H5D_fill_value_t fillValueStatus;
    if (H5Pfill_value_defined(createPropertyList, &fillValueStatus) < 0
        || fillValueStatus == H5D_FILL_VALUE_UNDEFINED)
    {
        numWrittenFrames_ = numValidFrames;
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
            gmx_file("Error reading step data set when determining the number of frames.");
        }
    }
    numWrittenFrames_ = numValidFrames + 1;
}

real GmxH5mdTimeDataBlock::getTimeOfFrame(hsize_t frame) const
{
    GMX_ASSERT(timeDataSet_ >= 0, "There must be a data set with time to get the time of a frame.");

    void*  buffer;
    size_t dataSize;
    readData<1, false>(timeDataSet_, frame, &buffer, &dataSize);

    if (dataSize != 4 && dataSize != 8)
    {
        gmx_file("Can only read float or double time data.");
    }

    if (dataSize == 4)
    {
        float* tmpFloatData = static_cast<float*>(buffer);
        return *tmpFloatData;
    }
    double* tmpDoubleData = static_cast<double*>(buffer);
    return *tmpDoubleData;
}

extern template hid_t
openOrCreateDataSet<1>(hid_t, const char*, const char*, hid_t, const hsize_t*, CompressionAlgorithm, double);
extern template hid_t
openOrCreateDataSet<3>(hid_t, const char*, const char*, hid_t, const hsize_t*, CompressionAlgorithm, double);

extern template void writeData<1, false>(hid_t, const void*, hsize_t);
extern template void writeData<3, false>(hid_t, const void*, hsize_t);

extern template void readData<1, false>(hid_t, hsize_t, void**, size_t*);
