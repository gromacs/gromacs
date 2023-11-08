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

#ifndef GMX_FILEIO_H5MD_UTIL_H
#define GMX_FILEIO_H5MD_UTIL_H

#include <string>

#include "gromacs/utility/real.h"


typedef int64_t            hid_t;
typedef unsigned long long hsize_t;
enum class PbcType : int;

/*! \brief An enumeration of compression options */
enum class CompressionAlgorithm
{
    None,                //!< No compression
    LosslessNoShuffle,   //!< Lossless, only gzip (deflate)
    LosslessWithShuffle, //!< Lossless, byte shuffle followed by gzip (deflate)
    LossySz3             //!< Lossy SZ3 compression.
};

/*! \brief Open an existing HDF5 group or create it if it did not exist already.
 *
 * \param[in] container  The container where the group is located, or should be created.
 * \param[in] name       The name of the group.
 * \returns the ID of the group.
 */
hid_t openOrCreateGroup(hid_t container, const char* name);

/*! \brief Registers the SZ3 filter by using the automatic registration mechanism by H5Pset_filter().
 * Must be done before appending (e.g. when restarting from acheckpoint) to a compressed dataset. */
void registerSz3FilterImplicitly();

/*! \brief Open an existing dataset (called name, in container). If it does not exist create a new dataset.
 *
 * \param[in] container The ID of the container of the data. This can be a group in the HDF5 or the HDF5 file itself.
 * \param[in] name The name of the data set.
 * \param[in] unit The unit of the data. See de Buyl et al., 2014 (https://www.sciencedirect.com/science/article/pii/S0010465514000447) for more information.
 * \param[in] datatype The HDF5 data type of the data.
 * \param[in] numFramesPerChunk The number of frames per chunk (compression unit) in the file.
 * \param[in] numEntries The number of particles, or similar.
 * \param[in] numValuesPerEntry The number of output values per entry (particle). This of often 1 or the number of dimensions, depending on the data.
 * \param[in] compression The compression algorithm to use.
 * \param[in] compressionError The required precision of lossy compression.
 * \returns The ID of the dataset.
 */

hid_t openOrCreateDataSet(hid_t                container,
                          const char*          name,
                          const char*          unit,
                          hid_t                dataType,
                          hsize_t              numFramesPerChunk,
                          hsize_t              numEntries,
                          hsize_t              numValuesPerEntry,
                          CompressionAlgorithm compression,
                          double               compressionError);

/*! \brief Writes an HDF5 data set, labelled by name, to the specified container.
 *
 * \param[in] dataSet The ID of the dataset to write to.
 * \param[in] data The data to write.
 * \param[in] positionToWrite The frame number to write.
 */
void writeData(hid_t dataSet, const void* data, hsize_t positionToWrite);

void setBoxGroupAttributes(hid_t boxGroup, PbcType pbcType);

template<typename T>
void setAttribute(hid_t container, const char* name, const T value, hid_t dataType);

void setAttribute(hid_t container, const char* name, const char* value);

template<hid_t numEntries, hid_t stringLength>
void setAttributeStringList(hid_t container, const char* name, const char value[numEntries][stringLength]);

#endif
