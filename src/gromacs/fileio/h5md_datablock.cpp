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

#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/gmxassert.h"

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

    mainDataSet_ = openOrCreateDataSet(group_,
                                       mainDataSetName,
                                       unit,
                                       datatype,
                                       numFramesPerChunk,
                                       numEntries,
                                       numValuesPerEntry,
                                       compression,
                                       compressionAbsoluteError);

    if (writingInterval_ > 0)
    {
        static constexpr char c_stepName[] = "step";
        stepDataSet_                       = openOrCreateDataSet(
                group_, c_stepName, nullptr, H5T_NATIVE_INT64, numFramesPerChunk, 1, 1, CompressionAlgorithm::None, 0);

        static constexpr char c_timeName[] = "time";
        static constexpr char c_timeUnit[] = "ps";
#if GMX_DOUBLE
        const hid_t timeDatatype = H5Tcopy(H5T_NATIVE_DOUBLE);
#else
        const hid_t timeDatatype = H5Tcopy(H5T_NATIVE_FLOAT);
#endif
        timeDataSet_ = openOrCreateDataSet(
                group_, c_timeName, c_timeUnit, timeDatatype, numFramesPerChunk, 1, 1, CompressionAlgorithm::None, 0);
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
    writeData(mainDataSet_, data, frameNumber);
    writeData(stepDataSet_, &step, frameNumber);
    writeData(timeDataSet_, &time, frameNumber);
}
