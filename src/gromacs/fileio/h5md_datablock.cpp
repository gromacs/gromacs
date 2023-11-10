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
                                           const char*          name,
                                           const char*          mainDataSetName,
                                           const char*          unit,
                                           int                  writingInterval,
                                           hsize_t              numFramesPerChunk,
                                           hsize_t              numEntries,
                                           hsize_t              numValuesPerEntry,
                                           hid_t                datatype,
                                           CompressionAlgorithm compression,
                                           double               compressionAbsoluteError)
{
    container_ = container;
    strncpy(name_, name, 128);
    writingInterval_ = writingInterval;

    group_ = openOrCreateGroup(container_, name_);
    H5Iget_name(group_, fullName_, c_maxFullNameLength);

    hsize_t chunkDims[3] = { numFramesPerChunk, numEntries, numValuesPerEntry };

    mainDataSet_ = openOrCreateDataSet<3>(
            group_, mainDataSetName, unit, datatype, chunkDims, compression, compressionAbsoluteError);

    if (writingInterval_ > 0)
    {
        hsize_t               chunkDimsTimeStep[1] = { numFramesPerChunk };
        static constexpr char c_stepName[]         = "step";
        stepDataSet_                               = openOrCreateDataSet<1>(
                group_, c_stepName, nullptr, H5T_NATIVE_INT64, chunkDimsTimeStep, CompressionAlgorithm::None, 0);

        static constexpr char c_timeName[] = "time";
        static constexpr char c_timeUnit[] = "ps";
#if GMX_DOUBLE
        const hid_t timeDatatype = H5Tcopy(H5T_NATIVE_DOUBLE);
#else
        const hid_t timeDatatype = H5Tcopy(H5T_NATIVE_FLOAT);
#endif
        timeDataSet_ = openOrCreateDataSet<1>(
                group_, c_timeName, c_timeUnit, timeDatatype, chunkDimsTimeStep, CompressionAlgorithm::None, 0);
    }
    else
    {
        timeDataSet_ = -1;
        stepDataSet_ = -1;
    }
}

GmxH5mdTimeDataBlock::GmxH5mdTimeDataBlock(const GmxH5mdTimeDataBlock& other) :
    container_(other.container_),
    group_(other.group_),
    mainDataSet_(other.mainDataSet_),
    timeDataSet_(other.timeDataSet_),
    stepDataSet_(other.stepDataSet_),
    writingInterval_(other.writingInterval_)
{
    strcpy(name_, other.name_);
    strcpy(fullName_, other.fullName_);
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

bool GmxH5mdTimeDataBlock::operator==(const char* fullNameComparison)
{
    if (strncmp(fullName_, fullNameComparison, c_maxFullNameLength) == 0)
    {
        return true;
    }
    return false;
}

void GmxH5mdTimeDataBlock::writeFrame(const void* data, int64_t step, real time)
{
    GMX_ASSERT(step == 0 || writingInterval_ > 0, "Invalid writing interval when writing frame.");

    int frameNumber = step > 0 ? step / writingInterval_ : 0;
    writeData<3, false>(mainDataSet_, data, frameNumber);
    writeData<1, false>(stepDataSet_, &step, frameNumber);
    writeData<1, false>(timeDataSet_, &time, frameNumber);
}

int64_t GmxH5mdTimeDataBlock::getNumberOfFrames()
{
    GMX_ASSERT(stepDataSet_ >= 0 || mainDataSet_ >= 0,
               "There must be a data set with steps or data values to determine the actual number "
               "of frames.");
    if (stepDataSet_ >= 0)
    {
        hid_t dataSpace = H5Dget_space(stepDataSet_);
        const int numDims = H5Sget_simple_extent_ndims(dataSpace);
        if (numDims != 1)
        {
            gmx_file("The step data set should be one-dimensional.");
        }
        hsize_t dimExtents;
        H5Sget_simple_extent_dims(dataSpace, &dimExtents, nullptr);
        int64_t numValidFrames = dimExtents;

        hid_t            createPropertyList = H5Dget_create_plist(stepDataSet_);
        H5D_fill_value_t fillValueStatus;
        if (H5Pfill_value_defined(createPropertyList, &fillValueStatus) < 0
            || fillValueStatus == H5D_FILL_VALUE_UNDEFINED)
        {
            return numValidFrames;
        }
        hid_t datatype = H5T_NATIVE_INT64;
        int64_t fillValue;
        H5Pget_fill_value(createPropertyList, datatype, &fillValue);

        int64_t stepData = fillValue;
        do {
            hsize_t location = numValidFrames - 1;
            H5Sselect_elements(dataSpace, H5S_SELECT_SET, 1, &location);
            if (H5Dread(stepDataSet_, datatype, dataSpace, H5S_ALL, H5P_DEFAULT, &stepData) < 0)
            {
                gmx_file("Error reading step data set when determining the number of frames.");
            }
            numValidFrames--;
        } while (numValidFrames > 0 && stepData == fillValue);

        return numValidFrames;
    }
    else
    {
        hid_t dataSpace = H5Dget_space(stepDataSet_);
        const int numDims = H5Sget_simple_extent_ndims(dataSpace);
        if (numDims != 3)
        {
            gmx_file("The time dependent data set should be three-dimensional.");
        }
        hsize_t dimExtents[3];
        H5Sget_simple_extent_dims(dataSpace, dimExtents, nullptr);
        int64_t numValidFrames = dimExtents[0];

        return numValidFrames;
    }
}

extern template hid_t
openOrCreateDataSet<1>(hid_t, const char*, const char*, hid_t, const hsize_t*, CompressionAlgorithm, double);
extern template hid_t
openOrCreateDataSet<3>(hid_t, const char*, const char*, hid_t, const hsize_t*, CompressionAlgorithm, double);

extern template void writeData<1, false>(hid_t, const void*, hsize_t);
extern template void writeData<3, false>(hid_t, const void*, hsize_t);
