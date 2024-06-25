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

/*! \brief Declarations for organizing time-dependent H5MD data blocks.
 *
 * \author Magnus Lundborg <lundborg.magnus@gmail.com>
 */

#ifndef GMX_FILEIO_H5MD_DATABLOCK_H
#define GMX_FILEIO_H5MD_DATABLOCK_H

#include <list>
#include <string>

/* This is currently written for use from GROMACS. If this is used in other packades,
 * the datatype real cannot be used. */
#include "gromacs/utility/real.h"

#include "h5md_util.h"

namespace gmx
{
namespace h5mdio
{
typedef int64_t            hid_t;
typedef unsigned long long hsize_t;

/*! \brief A class that handles H5MD data blocks with data can change during the MD trajectory.
 * Data is stored in three data sets, grouped together: main (value), time and step.
 */
class GmxH5mdTimeDataBlock
{
private:
    hid_t container_;   //!< The HDF5 container of this HDF5 group, storing the data sets.
    hid_t group_;       //!< The HDF5 ID of the group storing the data sets.
    hid_t mainDataSet_; //!< The ID of the main data set (values).
    hid_t timeDataSet_; //!< The ID of the time data set.
    hid_t stepDataSet_; //!< The ID of the data set storing simulation step numbers.

    /*! The name of the data block, the HDF5 group containing the data sets, e.g. "position". */
    std::string name_;
    std::string fullName_; //!< The full HDF5 path of the group storing the data sets.
    std::string mainUnit_; //!< The physical unit of the main (value) data.
    std::string timeUnit_; //!< The unit of the time data.

    /*! The index of the next frame to write. 0 when no frames have been written. */
    int64_t writingFrameIndex_;

    /*! The index of the next frame to read, 0 or the frame after the previously read frame. */
    int64_t readingFrameIndex_;

public:
    /*! \brief Create a management entity for a time dependent set of data.
     * \param[in] container The ID of the container (HDF5 group or file) of the data.
     * \param[in] name The name of this set of time dependent data - the H5MD group.
     * \param[in] unit The unit of the time dependent values.
     * \param[in] numFramesPerChunk Number of frames per chunk of data, relevant for compressed data.
     * \param[in] numEntries Number of data entries per frame, e.g., the number of atoms.
     * \param[in] numValuesPerEntry Number of data values per entry, e.g. 3 for 3D data.
     * \param[in] compression The compression algorithm to use.
     * \param[in] comopressionError The absolute error for lossy compression algorithms.
     */
    GmxH5mdTimeDataBlock(hid_t                container         = -1,
                         const std::string    name              = "",
                         const std::string    unit              = "",
                         hsize_t              numFramesPerChunk = 1,
                         hsize_t              numEntries        = 0,
                         hsize_t              numValuesPerEntry = 1,
                         hid_t                datatype          = -1,
                         CompressionAlgorithm compression       = CompressionAlgorithm::None,
                         double               compressionError  = 0.001);

    /* Close the data sets: main (or value), step and time. */
    void closeAllDataSets();

    /* Overloaded operator to compare two objects. N.b., only compares the fullName_,
     * does not check if files are the same. */
    bool operator==(const GmxH5mdTimeDataBlock& comparison);

    /* Overloaded operator to compare fullName_ to the fullNameComparison string. */
    bool operator==(const std::string fullNameComparison);

    /*! \brief Write a frame of time dependent data to the data block.
     *
     * \param[in] data The data that should be written.
     * \param[in] step The MD simulation step of the data record.
     * \param[in] time The time stamp (in ps) of the data record.
     * \param[in] frame The frame number to write. If <= 0 the frame to write will be
     * the one after the previously written frame.
     */
    void writeFrame(const void* data, int64_t step, real time, int64_t frame = -1);

    /*! \brief Read a specific frame.
     *
     * \param[out] data The data that is read. Memory must be allocated beforehand.
     * \param[in] frame The frame number to read.
     */
    bool readFrame(real* data, int64_t frame);

    /*! \brief Read the next, or the first frame, frame.
     *
     * \param[out] data The data that is read. Memory must be allocated beforehand.
     */
    bool readNextFrame(real* data);

    /*! \brief Read the units properties from file and update mainUnit_ and timeUnit_ accordingly.
     * FIXME There are no unit conversions yet.
     */
    void updateUnitsFromFile();

    /*! \brief Find out how many frames are written, ignoring fill value frames at the end.
     * Update writingFrameIndex_ to keep track of what is the next frame to write.
     */
    void updateNumWrittenFrames();

    /*! \brief Return the number of particles in the data block. */
    size_t getNumParticles() const;

    /*! \brief Return the MD simulation step of a given frame. */
    int64_t getStepOfFrame(int64_t frame) const;

    /*! \brief Return the time of a given frame. */
    real getTimeOfFrame(int64_t frame) const;

    /*! \brief Return the MD simulation step of the next frame to read or -1 if there are no more frames to read.. */
    int64_t getStepOfNextReadingFrame() const;

    /*! \brief Return the time of the next frame to read or -1 if there are no more frames to read. */
    real getTimeOfNextReadingFrame() const;

    /*! \brief Returns the compression error of lossy SZ3 compression, or -1 if there is no lossy SZ3 compression. */
    double getLossyCompressionError() const;

    /*! \brief Get the number of frames in the data block. */
    int64_t numberOfFrames() const { return writingFrameIndex_; }

    /*! \brief Get the name of the data block. */
    std::string name() const { return name_; }

    /*! \brief Get the path to the data block. */
    std::string fullName() const { return fullName_; }
    /*! \brief Get the frame number of the next frame to write (in sequence). */
    int64_t writingFrameIndex() const { return writingFrameIndex_; }
    /*! \brief Get the frame number of the next frame to read (in sequence). */
    int64_t readingFrameIndex() const { return readingFrameIndex_; }
    /*! \brief Get the unit of the data stored in the data block. */
    std::string mainUnit() const { return mainUnit_; }
    /*! \brief Get the unit of the time stamps in the data block. */
    std::string timeUnit() const { return timeUnit_; }
};

} // namespace h5mdio
} // namespace gmx
#endif // GMX_FILEIO_H5MD_DATABLOCK_H
