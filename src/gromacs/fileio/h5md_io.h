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

#ifndef GMX_FILEIO_H5MD_IO_H
#define GMX_FILEIO_H5MD_IO_H

#include <string>
#include <vector>

#include <sys/_types/_int64_t.h>

#include "gromacs/math/vectypes.h"
#include "gromacs/utility/arrayref.h"
#include "gromacs/utility/real.h"

#include "h5md_datablock.h"

struct gmx_mtop_t;
typedef int64_t            hid_t;
typedef unsigned long long hsize_t;
enum class PbcType : int;
constexpr size_t c_atomNameLen = 17;


/*! \brief The container of the H5MD data. The class is designed to read/write data according to de Buyl et al., 2014
 * (https://www.sciencedirect.com/science/article/pii/S0010465514000447) and https://www.nongnu.org/h5md/h5md.html
 * The class contains a number of standard data blocks that are commonly used by GROMACS. */
class GmxH5mdIo
{
private:
    hid_t file_; //!< The HDF5 identifier of the file. This is the H5MD root.
    std::list<GmxH5mdTimeDataBlock> dataBlocks_; //!< A list of time dependent data blocks in the HDF5 file.

    std::string systemOutputName_; //!< The name of the selection of particles output. Defaults to "system"

    /*! \brief Sets the author (user) and creator (application name) properties in the h5md group (h5mdGroup_). */
    void setAuthorAndCreator();

public:
    /*! \brief Construct a GmxH5mdIo object and open a GmxHdf5 file.
     *
     * \param[in] fileName    Name of the file to open. The same as the file path.
     * \param[in] mode        The mode to open the file, described by a lower-case letter
     *                        'w' means writing (and reading), i.e. backup an existing file and
     * replace it. 'a' means appending (and reading), i.e., that existing files will be not be
     * overwritten, but extended. 'r' means only reading.
     * \throws FileIOError if fileName is specified and the file cannot be opened.
     */
    GmxH5mdIo(const std::string fileName = "", const char mode = '\0');

    ~GmxH5mdIo();

    /*! \brief Open an H5MD file.
     *
     * \param[in] fileName    Name of the file to open. The same as the file path.
     * \param[in] mode        The mode to open the file, described by a case-insensitive string of
     *                        letters, up to three characters long. Reading is always assumed.
     *                        'w' means writing, i.e. backup an existing file and replace it,
     *                        'a' means truncate, i.e., that existing files will be overwritten
     *                        'r' means only read.
     * \throws FileIOError    If the file cannot be opened.
     */
    void openFile(const std::string fileName, const char mode);

    /*! \brief Close the H5MD file.
     * \throws FileIOError    If the file cannot be closed.
     */
    void closeFile();

    /*! \brief Write all unwritten data to the file.
     * \throws FileIOError    If there were errors during flushing.
     */
    void flush();

    /*! \brief Is there an open file? */
    bool isFileOpen() const { return file_ > 0; }

    /*! \brief Create and initialize time dependent particles data block objects from the H5MD file.
     * \throws FileIOError    If the particles data could not be read from the file. Invalid H5MD file?
     */
    void initParticleDataBlocksFromFile();

    /*! \brief Set up data blocks related to particle data.
     *
     * This needs to be done before writing the particle data to the trajectory.
     *
     * \param[in] writeCoordinatesSteps The lossless coordinate output interval.
     * \param[in] writeVelocitiesSteps The lossless velocity output interval.
     * \param[in] writeForcesSteps The lossless force output interval.
     * \param[in] numParticles The number of particles/atoms in the system.
     * \param[in] pbcType The periodic boundary condition that is used.
     * \param[in] xCompressionError The accepted error of the lossy compression of coordinates.
     * \throws FileIOError If the data blocks could not be created.
     */
    void setUpParticlesDataBlocks(int64_t writeCoordinatesSteps,
                                  int64_t writeVelocitiesSteps,
                                  int64_t writeForcesSteps,
                                  int64_t numParticles,
                                  PbcType pbcType,
                                  double  xCompressionError);

    /*! \brief Write molecule system related data to the file.
     *
     * This is currently not updated during the trajectory. The data that is written are atom masses, atom charges and atom names.
     *
     * \param[in] topology The molecular topology describing the system.
     * \param[in] index    The selected atoms to include. If empty, use all atoms in the topology.
     * \param[in] index_group_name The name of the atom selection specified by index.
     * \throws FileIOError    If there is no file open.
     */
    void setupMolecularSystem(const gmx_mtop_t&        topology,
                              gmx::ArrayRef<const int> index            = {},
                              const std::string        index_group_name = "");

    std::vector<std::string> readAtomNames();

    /*! \brief Write a trajectory frame to the file. Only writes the data that is passed as input
     *
     * \param[in] step The simulation step.
     * \param[in] time The time stamp (in ps).
     * \param[in] lambda The lambda state. FIXME: Currently not written.
     * \param[in] box The box dimensions.
     * \param[in] numParticles The number of particles in the system.
     * \param[in] x The particle coordinates for lossless output.
     * \param[in] v The particle velocities for lossless output.
     * \param[in] f The particle forces for lossless output.
     * \param[in] xCompressionError The accepted error of the lossy compression of coordinates.
     * \throws FileIOError    If there is no file open or if errors occured during writing.
     */
    void writeFrame(int64_t     step,
                    real        time,
                    real        lambda,
                    const rvec* box,
                    int64_t     numParticles,
                    const rvec* x,
                    const rvec* v,
                    const rvec* f,
                    double      xCompressionError);

    /*! \brief Read the next frame of box, coordinates, velocities and forces. With next frame means
     * the lowest step/time reading from the previous read frame of that data type. If data is written
     * at different intervals the read data types will be different from one function call to the next.
     * \param[out] step    The step number of the read data blocks.
     * \param[out] time    The time of the read data blocks.
     * \param[out] box     Read box data. Memory must be allocated by the caller.
     * \param[out] x       Read coordinate data. Memory must be allocated by the caller.
     * \param[out] v       Read velocity data. Memory must be allocated by the caller.
     * \param[out] f       Read force data. Memory must be allocated by the caller.
     * \param[out] xCompressionError The error of lossy (SZ3) coordinate compression. -1 if no lossy SZ3 compression.
     * \param[out] readBox Whether box data was read or not, i.e. if there was box data matching step.
     * \param[out] readX   Whether coordinate data was read or not, i.e. if there was coordinate data matching step.
     * \param[out] readV   Whether velocity data was read or not, i.e. if there was velocity data matching step.
     * \param[out] readF   Whether force data was read or not, i.e. if there was force data matching step.
     * \returns Whether any frame was read or not.
     * \throws FileIOError If there was an error reading the next frame or if the data type of the read data was unknown.
     */
    bool readNextFrameOfStandardDataBlocks(int64_t* step,
                                           real*    time,
                                           rvec*    box,
                                           rvec*    x,
                                           rvec*    v,
                                           rvec*    f,
                                           real*    xCompressionError,
                                           bool*    readBox,
                                           bool*    readX,
                                           bool*    readV,
                                           bool*    readF);

    /*! \brief Get the number of frames of a particles data block
     * \param[in] dataBlockName The name of the data block.
     * \returns the number of frames of the data block or -1 if not found.
     */
    int64_t getNumberOfFrames(const std::string dataBlockName);

    /*! \brief Get the number of particles of a particles data block
     * \param[in] dataBlockName The name of the data block.
     * \returns the number of particles of the data block or -1 if not found.
     * \throws FileIOError    If there was an error reading the number of particles could not be
     * read, such as no atom data in particles data blocks.
     */
    int64_t getNumberOfParticles(const std::string dataBlockName);

    /*! \brief Get the first time stamp of a particles data block
     * \param[in] dataBlockName The name of the data block.
     * \returns the first time of the data block or -1 if not found.
     * \throws FileIOError    If there was an error determining the time of the first frame.
     */
    real getFirstTime(const std::string dataBlockName);

    /*! \brief Get the very first time stamp of all particles data blocks
     * \returns the first time of all data blocks or -1 if no data blocks were found.
     * \throws FileIOError    If there was an error determining the time of the first frame.
     */
    real getFirstTimeFromAllDataBlocks();

    /*! \brief Get the final time stamp of a particles data block
     * \param[in] dataBlockName The name of the data block.
     * \returns the final time of the data block or -1 if not found.
     * \throws FileIOError    If there was an error determining the time of the final frame.
     */
    real getFinalTime(const std::string dataBlockName);

    /*! \brief Get the very last time stamp of all particles data blocks
     * \returns the last time of all data blocks or -1 if no data blocks were found.
     * \throws FileIOError    If there was an error determining the time of the final frame.
     */
    real getFinalTimeFromAllDataBlocks();
};

#endif // GMX_FILEIO_H5MD_IO_H
