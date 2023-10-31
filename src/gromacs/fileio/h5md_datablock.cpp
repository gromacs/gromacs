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

#include "h5md_util.h"

#define GMX_USE_HDF5 1 // FIXME: Temporary just for the editor

#if GMX_USE_HDF5
#include <hdf5.h>
#endif


// GmxH5mdDataBlock::GmxH5mdDataBlock()
// {
//     container_ = -1;
//     strcpy(name_, "");
//     numFramesPerChunk_ = 1;
//     numEntries_ = 0;
//     numValuesPerEntry_ = 1;
//     numWrittenFrames_ = 0;
//     compressionAlgorithm_ = CompressionAlgorithm::None;
//     compressionAbsoluteError_ = 0;
//     datatype_ = -1;
// }

GmxH5mdDataBlock::GmxH5mdDataBlock(hid_t container, const char* name, const char* unit, hsize_t numFramesPerChunk, hsize_t numEntries,
                                   hsize_t numValuesPerEntry, hid_t datatype, CompressionAlgorithm compression, double compressionError)
{
    container_ = container;
    strncpy(name_, name, 128);
    strncpy(unit_, unit, 64);
    numFramesPerChunk_ = numFramesPerChunk;
    numEntries_ = numEntries;
    numValuesPerEntry_ = numValuesPerEntry;
    numWrittenFrames_ = 0;
    datatype_ = datatype;
    compressionAlgorithm_ = compression;
    compressionAbsoluteError_ = compressionError;
}

void GmxH5mdDataBlock::initDataProperties(hsize_t numFramesPerChunk, hsize_t numEntries, hsize_t numValuesPerEntry)
{
    if (numWrittenFrames_ != 0)
    {
        gmx_file("Cannot change number of frames per chunk after writing.");
    }
    numFramesPerChunk_ = numFramesPerChunk;
    numEntries_ = numEntries;
    numValuesPerEntry_ = numValuesPerEntry;
}

void GmxH5mdDataBlock::writeFrame(const void* data, int64_t step, real time, const char* valueName)
{
#if GMX_DOUBLE
    const hid_t timeDatatype = H5Tcopy(H5T_NATIVE_DOUBLE);
#else
    const hid_t timeDatatype = H5Tcopy(H5T_NATIVE_FLOAT);
#endif
    char stepName[] = "step";
    char timeName[] = "time";
    char timeUnit[] = "ps";
    hid_t group = openOrCreateGroup(container_, name_);

    writeData(group, valueName, unit_, data, numFramesPerChunk_, numEntries_, numValuesPerEntry_, numWrittenFrames_, datatype_, compressionAlgorithm_, compressionAbsoluteError_);
    writeData(group, stepName, nullptr, &step, numFramesPerChunk_, 1, 1, numWrittenFrames_, H5T_NATIVE_INT64, CompressionAlgorithm::None, 0);
    writeData(group, timeName, timeUnit, &time, numFramesPerChunk_, 1, 1, numWrittenFrames_, timeDatatype, CompressionAlgorithm::None, 0);
    H5Gclose(group);
    numWrittenFrames_++;
}

