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

#include <vector>

#include "gromacs/fileio/h5md/h5md_datasetbase.h"

namespace gmx
{

template<typename ValueType>
class BasicVector;

/*! \brief Dimensions of a HDF5 data set.
 */
using DataSetDims = std::vector<hsize_t>;

/*! \brief Return the number of frames in a data set.
 *
 * Number of frames refers to the size of the major axis of the data set. For example,
 * a 1d data set with dimensions [30] obviously has 30 frames. A 3d data set with
 * dimensions [30, 150, 3] also has 30 frames.
 *
 * \param[in] dataSet Handle to data set.
 * \returns Number of frames.
 */
template<typename ValueType>
hsize_t getNumFrames(const H5mdDataSetBase<ValueType>& dataSet);

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
template<typename ValueType>
void setNumFrames(const H5mdDataSetBase<ValueType>& dataSet, hsize_t numFrames);

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
template<typename ValueType>
hid_t getFrameDataSpace(const H5mdDataSetBase<ValueType>& dataSet, hsize_t frameIndex);

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
template<typename ValueType>
hid_t getFrameMemoryDataSpace(const H5mdDataSetBase<ValueType>& dataSet);

extern template hsize_t getNumFrames(const H5mdDataSetBase<int32_t>&);

extern template hsize_t getNumFrames(const H5mdDataSetBase<int64_t>&);

extern template hsize_t getNumFrames(const H5mdDataSetBase<float>&);

extern template hsize_t getNumFrames(const H5mdDataSetBase<double>&);

extern template hsize_t getNumFrames(const H5mdDataSetBase<gmx::BasicVector<float>>&);

extern template hsize_t getNumFrames(const H5mdDataSetBase<gmx::BasicVector<double>>&);

extern template void setNumFrames(const H5mdDataSetBase<int32_t>&, hsize_t);

extern template void setNumFrames(const H5mdDataSetBase<int64_t>&, hsize_t);

extern template void setNumFrames(const H5mdDataSetBase<float>&, hsize_t);

extern template void setNumFrames(const H5mdDataSetBase<double>&, hsize_t);

extern template void setNumFrames(const H5mdDataSetBase<gmx::BasicVector<float>>&, hsize_t);

extern template void setNumFrames(const H5mdDataSetBase<gmx::BasicVector<double>>&, hsize_t);

extern template hid_t getFrameDataSpace(const H5mdDataSetBase<int32_t>&, hsize_t);

extern template hid_t getFrameDataSpace(const H5mdDataSetBase<int64_t>&, hsize_t);

extern template hid_t getFrameDataSpace(const H5mdDataSetBase<float>&, hsize_t);

extern template hid_t getFrameDataSpace(const H5mdDataSetBase<double>&, hsize_t);

extern template hid_t getFrameDataSpace(const H5mdDataSetBase<BasicVector<float>>&, hsize_t);

extern template hid_t getFrameDataSpace(const H5mdDataSetBase<BasicVector<double>>&, hsize_t);

extern template hid_t getFrameMemoryDataSpace(const H5mdDataSetBase<int32_t>&);

extern template hid_t getFrameMemoryDataSpace(const H5mdDataSetBase<int64_t>&);

extern template hid_t getFrameMemoryDataSpace(const H5mdDataSetBase<float>&);

extern template hid_t getFrameMemoryDataSpace(const H5mdDataSetBase<double>&);

extern template hid_t getFrameMemoryDataSpace(const H5mdDataSetBase<BasicVector<float>>&);

extern template hid_t getFrameMemoryDataSpace(const H5mdDataSetBase<BasicVector<double>>&);

} // namespace gmx

#endif // GMX_FILEIO_H5MD_DATASET_H
