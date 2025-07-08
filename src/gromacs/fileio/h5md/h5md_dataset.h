
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

/*! \brief Declarations of H5md data set manipulation routines.
 *
 * \author Magnus Lundborg <lundborg.magnus@gmail.com>
 * \author Petter Johansson <pettjoha@kth.se>
 */

#ifndef GMX_FILEIO_H5MD_DATASET_H
#define GMX_FILEIO_H5MD_DATASET_H

#include <hdf5.h>

#include <string>
#include <vector>

namespace gmx
{

/*! \brief Dimensions of a HDF5 data set.
 */
using DataSetDims = std::vector<hsize_t>;

/*! \brief Open an existing data set (called \p dataSetName, in \p container).
 *
 * The returned handle must be closed with H5Dclose to avoid resource leaks.
 *
 * \param[in] container The ID of the container of the data.
 * \param[in] dataSetName The name of the data set.
 * \returns The ID of the data set.
 *
 * \throws gmx::FileIOError if data set cannot be opened.
 */
hid_t openDataSet(const hid_t container, const std::string& dataSetName);

/*! \brief Return the dimensions of a data set.
 *
 * \param[in] dataSet Handle to data set.
 * \returns Vector containing the dimensions.
 *
 * \throws gmx::FileIOError if the dimensions cannot be read.
 */
DataSetDims getDataSetDims(const hid_t dataSet);

/*! \brief Return the number of frames in a data set.
 *
 * Number of frames refers to the size of the major axis of the data set. For example,
 * a 1d data set with dimensions [30] obviously has 30 frames. A 3d data set with
 * dimensions [30, 150, 3] also has 30 frames.
 *
 * \param[in] dataSet Handle to data set.
 * \returns Number of frames.
 */
hsize_t getNumFrames(const hid_t dataSet);

/*! \brief Resize \p dataSet to store \p numFrames frames.
 *
 * This resizes the data set along the major axis (=0). It does not alter the size
 * along other axes.
 *
 * If the new number of frames is smaller than the current size, data will be discarded.
 *
 * \param[in] dataSet   Handle to data set.
 * \param[in] numFrames New number of frames for data set.
 */
void setNumFrames(const hid_t dataSet, hsize_t numFrames);

/*! \brief Return a data space for a given data set, with a selected hyperslab for the frame.
 *
 * The returned handle must be closed with H5Sclose to avoid resource leaks.
 *
 * When reading/writing data from/to a given frame we need to select the memory inside
 * of the data set where this data is placed. This is selected by creating a hyperslab
 * for the chunk of memory corresponding to this frame index, and creating a data space
 * which refers to this slab. This data space is then used by the reading/writing routines
 * of HDF5 to perform the action.
 *
 * \param[in] dataSet Handle to data set.
 * \param[in] frameIndex Index of frame to select hyperslab for.
 * \returns Handle to data space with the selected hyperslab.
 */
hid_t getFrameDataSpace(const hid_t dataSet, const hsize_t frameIndex);

/*! \brief Return a data space for storing a single frame from a data set.
 *
 * The returned handle must be closed with H5Sclose to avoid resource leaks.
 *
 * When reading/writing data from/to a given frame we need a data space of corresponding
 * size for handling. This function returns a space with the correct size for single
 * frame memory chunks, to be used by the reading/writing routines of HDF5 to perform
 * the action.
 *
 * \param[in] dataSet Handle to data set.
 * \returns Handle to data space with storage for a single frame.
 */
hid_t getFrameMemoryDataSpace(const hid_t dataSet);

} // namespace gmx

#endif // GMX_FILEIO_H5MD_DATASET_H
