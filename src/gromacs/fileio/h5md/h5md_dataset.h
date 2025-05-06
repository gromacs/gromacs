
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

#include <array>
#include <string>

namespace gmx
{

/*! \brief Dimensions of a HDF5 data set.
 * \tparam numDims Number of dimensions of the data set.
 */
template<int numDims>
using DataSetDims = std::array<hsize_t, numDims>;

/*! \brief Create a 1d data set for per-frame data (called \p dataSetName, in \p container).
 *
 * The data set is empty and has no limit on the number of values it can contain. This informs
 * its purpose: to append values corresponding to different frames. Since it can grow it should
 * not be used to create a fixed-size data set.
 *
 * Use this to create data sets for e.g. storing time stamps, simulation steps or energies.
 *
 * Note that this should only be used for simple primitive \p ValueType types. More complicated
 * types such as strings, gmx::RVecs or similar should use specialized functions.
 *
 * The returned handle must be closed with H5Dclose to avoid resource leaks.
 *
 * \tparam ValueType The type of data to store.
 * \param[in] container The ID of the container of the data.
 * \param[in] dataSetName The name of the data set.
 * \returns The ID of the data set.
 *
 * \throws gmx::FileIOError if the data set cannot be created.
 */
template<typename ValueType>
hid_t create1dFrameDataSet(const hid_t container, const std::string& dataSetName);

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

/*! \brief Return the number of frames in a data set.
 *
 * Number of frames refers to the size of the major axis of the data set. For example,
 * a 1d data set with dimensions [30] obviously has 30 frames. A 3d data set with
 * dimensions [30, 150, 3] also has 30 frames.
 *
 * \tparam numDims Number of dimensions of data set.
 * \param[in] dataSet Handle to data set.
 * \returns Number of frames.
 */
template<int numDims>
hsize_t getNumFrames(const hid_t dataSet);

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
 * \tparam numDims Number of dimensions of data set.
 * \param[in] dataSet Handle to data set.
 * \param[in] frameIndex Index of frame to select hyperslab for.
 * \returns Handle to data space with the selected hyperslab.
 */
template<int numDims>
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
 * \tparam numDims Number of dimensions of data set.
 * \param[in] dataSet Handle to data set.
 * \returns Handle to data space with storage for a single frame.
 */
template<int numDims>
hid_t getFrameMemoryDataSpace(const hid_t dataSet);

extern template hid_t create1dFrameDataSet<int32_t>(const hid_t, const std::string&);

extern template hid_t create1dFrameDataSet<int64_t>(const hid_t, const std::string&);

extern template hid_t create1dFrameDataSet<float>(const hid_t, const std::string&);

extern template hid_t create1dFrameDataSet<double>(const hid_t, const std::string&);

extern template hsize_t getNumFrames<1>(const hid_t dataSet);

extern template hid_t getFrameDataSpace<1>(const hid_t dataSet, const hsize_t frameIndex);

extern template hid_t getFrameMemoryDataSpace<1>(const hid_t);

} // namespace gmx

#endif // GMX_FILEIO_H5MD_DATASET_H
