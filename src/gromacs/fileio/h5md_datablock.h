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

#ifndef GMX_FILEIO_H5MD_DATABLOCK_H
#define GMX_FILEIO_H5MD_DATABLOCK_H

#include <list>
#include <string>

#include <sys/_types/_int64_t.h>

#include "gromacs/math/vectypes.h"
#include "gromacs/utility/real.h"

#include "h5md_util.h"

typedef int64_t            hid_t;
typedef unsigned long long hsize_t;

constexpr int c_maxNameLength     = 128;
constexpr int c_maxFullNameLength = 256;

/*! \brief A class that handles H5MD data blocks with data can change during the MD trajectory. */
class GmxH5mdTimeDataBlock
{
private:
    hid_t       container_; //!< The HDF5 container of this data block.
    std::string name_;      //!< The name of the data block, e.g. "position".
    std::string fullName_;
    std::string mainUnit_;
    std::string timeUnit_;
    hid_t       group_;
    hid_t       mainDataSet_;
    hid_t       timeDataSet_;
    hid_t       stepDataSet_;
    int         writingInterval_; //!< The interval (in MD steps) between outputs.
    int         writingFrameIndex_;
    int         readingFrameIndex_;

public:
    GmxH5mdTimeDataBlock(hid_t                container         = -1,
                         const std::string    name              = "",
                         const std::string    unit              = "",
                         int                  writingInterval   = 0,
                         hsize_t              numFramesPerChunk = 1,
                         hsize_t              numEntries        = 0,
                         hsize_t              numValuesPerEntry = 1,
                         hid_t                datatype          = -1,
                         CompressionAlgorithm compression       = CompressionAlgorithm::None,
                         double               compressionError  = 0.001);

    void closeAllDataSets();

    bool operator==(const std::string fullSpecifier);

    /*! \brief Write a set of time independent data to the data block.
     *
     * \param[in] data The data that should be written.
     */
    void writeTimeIndependentData(const void* data);

    /*! \brief Write a frame of time dependent data to the data block.
     *
     * \param[in] data The data that should be written.
     * \param[in] step The MD simulation step of the data record.
     * \param[in] time The time stamp (in ps) of the data record.
     */
    void writeFrame(const void* data, int64_t step, real time);

    void writeFrame(const void* data, int64_t step, real time, int frame);

    bool readFrame(real* data, int frame);

    bool readNextFrame(real* data);

    void updateUnitsFromFile();

    void updateNumWrittenFrames();

    size_t getNumParticles() const;

    int64_t getStepOfFrame(hsize_t frame) const;

    real getTimeOfFrame(hsize_t frame) const;

    int         numberOfFrames() const { return writingFrameIndex_; }
    std::string name() const { return name_; }
    std::string fullName() const { return fullName_; }
    int         writingFrameIndex() const { return writingFrameIndex_; }
    int         readingFrameIndex() const { return readingFrameIndex_; }
    std::string mainUnit() const { return mainUnit_; }
    std::string timeUnit() const { return timeUnit_; }
};

#endif // GMX_FILEIO_H5MD_DATABLOCK_H
