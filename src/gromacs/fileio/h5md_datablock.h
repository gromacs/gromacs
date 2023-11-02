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

#include <string>

#include "h5md_util.h"
#include "gromacs/utility/real.h"

typedef int64_t hid_t;
typedef unsigned long long hsize_t;

/*! \brief A class that handles H5MD data blocks with data can change during the MD trajectory. */
class GmxH5mdDataBlock
{
private:
    hid_t    container_; //!< The HDF5 container of this data block.
    char     name_[128]; //!< The name of the data block, e.g. "position".
    char     datasetName_[64]; //!< The label of the dataset in the data block, e.g. "value".
    char     unit_[64]; //!< The unit of the data in the data block. The unit of the time records is automatically set to "ps".
    int      writingInterval_; //!< The interval (in MD steps) between outputs.
    hsize_t  numFramesPerChunk_; //!< Number of frames per HDF5 chunk, i.e. the blocks of data that are compressed, written and read together.
    hsize_t  numEntries_; //!< The number of entries per frame. For particle data this is the number of atoms. For box dimensions in GROMACS it is DIM.
    hsize_t  numValuesPerEntry_; //!< The dimensionality of each data value. For particle data it is often DIM or 1 (for e.g., masses).
    hid_t    datatype_; //!< The HDF5 datatype of this data block.
    CompressionAlgorithm compressionAlgorithm_; //!< The compression algorithm that should be used.
    double compressionAbsoluteError_; //!< The absolute error of lossy compression algorithms.
public:
    GmxH5mdDataBlock(hid_t container = -1, const char* name = "", const char* datasetName = "value", const char* unit = "", int writingInterval = 0, hsize_t numFramesPerChunk = 1, hsize_t numEntries = 0, hsize_t numValuesPerEntry = 1, hid_t datatype = -1,
                     CompressionAlgorithm compression = CompressionAlgorithm::None, double compressionError = 0.001);

    /*! \brief Write a frame of data to the data block.
     *
     * \param[in] data The data that should be written.
     * \param[in] step The MD simulation step of the data record.
     * \param[in] time The time stamp (in ps) of the data record.
     */
    void writeFrame(const void* data, int64_t step, real time);
};

#endif // GMX_FILEIO_H5MD_DATABLOCK_H
