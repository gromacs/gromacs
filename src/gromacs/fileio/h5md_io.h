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

#include "gromacs/math/vectypes.h"
#include "gromacs/utility/arrayref.h"
#include "gromacs/utility/real.h"

#include "h5md_datablock.h"
#include "h5md_util.h"

enum class PbcType : int;
struct gmx_mtop_t;

namespace gmx
{
namespace h5mdio
{

typedef int64_t            hid_t;
typedef unsigned long long hsize_t;
constexpr size_t           c_atomNameLen = 17;


/*! \brief The container of the H5MD data. The class is designed to read/write data according to de Buyl et al., 2014
 * (https://www.sciencedirect.com/science/article/pii/S0010465514000447) and https://www.nongnu.org/h5md/h5md.html
 * The class contains a number of standard data blocks that are commonly used by GROMACS. */
class GmxH5mdIo
{
private:
    hid_t file_; //!< The HDF5 identifier of the file. This is the H5MD root.
    std::list<GmxH5mdTimeDataBlock> dataBlocks_; //!< A list of time dependent data blocks in the HDF5 file.

    /* FIXME: systemOutputName_ is a hacky solution. */
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

    /*! \brief Create and initialize time dependent data block objects from the H5MD file.
     * \param groupName       The name of the group with blocks to create and initialize.
     * \throws FileIOError    If there were errors reading data blocks from the file. Invalid file?
     * \returns The number of initialized data blocks.
     */
    int initGroupTimeDataBlocksFromFile(std::string groupName);

    /*! \brief Set the author name attribute in the H5MD file
     * \param[in] authorName The author name.
     * \throws FileIOError If the author name attribute could not be set.
     */
    void setAuthor(std::string authorName);

    /*! \brief Get the author name attribute from the H5MD file
     * \returns the author name.
     */
    std::string getAuthor();

    /*! \brief Set the name of the creating program as an attribute in the H5MD file.
     * \param[in] creatorName The creator name, i.e. the name of the program that created the file.
     * \throws FileIOError If the author name attribute could not be set.
     */
    void setCreatorProgramName(std::string creatorName);

    /*! \brief Get the name of the creating program attribute from the H5MD file.
     * \returns he creator name, i.e. the name of the program that created the file.
     */
    std::string getCreatorProgramName();

    /*! \brief Set the version of the creating program as an attribute in the H5MD file.
     * \param[in] version The version.
     * \throws FileIOError If the author name attribute could not be set.
     */
    void setCreatorProgramVersion(std::string version);

    /*! \brief Get the version of the creating program attribute from the H5MD file.
     * \returns the version.
     */
    std::string getCreatorProgramVersion();

    /* FIXME: This should be handled better. There should be ways to use multiple selections. */
    /*! Set the name of the selection used for writing output.
     * \param[in] systemOutputName The name to use for the selection.
     */
    void setSystemOutputName(std::string systemOutputName);

    /*! Get the name of the selection used for writing output.
     * \returns The name used for the output.
     */
    std::string getSystemOutputName();

    /*! \brief Set atom names, if not already set. Atom names cannot be modified after setting.
     *
     * \param[in] atomNames The names of the atoms.
     * \param[in] selectionName The name of the selection the atoms belong to.
     * \throws FileIOError If there was an error creating the data block or writing the data.
     */
    void setAtomNames(const std::vector<std::string>& atomNames,
                      std::string                     selectionName = "system");

    /*! \brief Set atom partial charges, if not already set. Partial charges cannot be modified after setting.
     *
     * \param[in] atomCharges The names of the atoms.
     * \param[in] selectionName The name of the selection the atoms belong to.
     * \throws FileIOError If there was an error creating the data block or writing the data.
     */
    void setAtomPartialCharges(const std::vector<real>& atomCharges,
                               std::string              selectionName = "system");

    /*! \brief Set atom masses, if not already set. Masses cannot be modified after setting.
     *
     * \param[in] atomMasses The names of the atoms.
     * \param[in] selectionName The name of the selection the atoms belong to.
     * \throws FileIOError If there was an error creating the data block or writing the data.
     */
    void setAtomMasses(const std::vector<real>& atomMasses, std::string selectionName = "system");

    std::vector<std::string> readAtomNames();

    /*! \brief Write a frame of data to the file.
     * \param[in] step The simulation step.
     * \param[in] time The time stamp (in ps).
     * \param[in] dataSetFullName Full path (including groups) to the data set.
     * \param[in] dataDimensionalityFirstDim The number of data entries, e.g., the number of particles.
     *                                       Must be > 0.
     * \param[in] dataDimensionalitySecondDim The number of dimesions per particles, e.g., 3 for coordinates.
     *                                        Must be > 0.
     * \param[in] data The data to write.
     * \param[in] unit The (SI) unit of the data.
     * \param[in] numberOfFramesPerChunk The number of frames per compression chunk. Cannot change.
     * \param[in] compressionAlgorithm The compression algorithm to use when writing the data.
     * \param[in] compressionError The accepted error of the lossy compression (if applicable).
     * \throws FileIOError    If there is no file open or if errors occured during writing.
     */
    void writeDataFrame(int64_t              step,
                        real                 time,
                        std::string          dataSetFullName,
                        int                  dataDimensionalityFirstDim,
                        int                  dataDimensionalitySecondDim,
                        const real*          data,
                        std::string          unit                   = "",
                        hsize_t              numberOfFramesPerChunk = 1,
                        CompressionAlgorithm compressionAlgorithm   = CompressionAlgorithm::None,
                        double               lossyCompressionError  = 0);

    /*! \brief Read the next frame of box, coordinates, velocities and forces. With next frame means
     * the lowest step/time reading from the previous read frame of that data type. If data is
     * written at different intervals the read data types will be different from one function call
     * to the next.
     * \param[out] step       The step number of the read data blocks.
     * \param[out] time       The time of the read data blocks.
     * \param[out] lambda     Read lambda data. Memory (1 value) must be allocated by the caller.
     * \param[out] box        Read box data. Memory must be allocated by the caller.
     * \param[out] x          Read coordinate data. Memory must be allocated by the caller.
     * \param[out] v          Read velocity data. Memory must be allocated by the caller.
     * \param[out] f          Read force data. Memory must be allocated by the caller.
     * \param[out] xCompressionError The error of lossy (SZ3) coordinate compression. -1 if no lossy SZ3.
     * \param[out] readLambda Whether lambda data was read or not, i.e. if there was lambda data at this step.
     * \param[out] readBox    Whether box data was read or not, i.e. if there was box data at this step.
     * \param[out] readX      Whether coordinate data was read or not, i.e. if there was coordinate at this step.
     * \param[out] readV      Whether velocity data was read or not, i.e. if there was velocity data at this step.
     * \param[out] readF      Whether force data was read or not, i.e. if there was force data at this step.
     * \returns Whether any frame was read or not.
     * \throws FileIOError If there was an error reading the next frame or if the data type of the data was unknown.
     */
    bool readNextFrameOfStandardDataBlocks(int64_t* step,
                                           real*    time,
                                           real*    lambda,
                                           rvec*    box,
                                           rvec*    x,
                                           rvec*    v,
                                           rvec*    f,
                                           real*    xCompressionError,
                                           bool*    readLambda,
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

} // namespace h5mdio

/*! \brief Set the author and creator fields based on user and program data available in gromacs.

 * \param[in] file     The H5MD file manager to use.
 */
void setH5mdAuthorAndCreator(h5mdio::GmxH5mdIo* file);

/*! \brief Write molecule system related data to the file.
 *
 * This is currently not updated during the trajectory. The data that is written are atom masses, atom charges and atom names.
 *
 * \param[in] file     The H5MD file manager to use.
 * \param[in] topology The molecular topology describing the system.
 * \param[in] index    The selected atoms to include. If empty, use all atoms in the topology.
 * \param[in] index_group_name The name of the atom selection specified by index.
 * \throws FileIOError    If there is no file open or if the data could not be written.
 */
void setupMolecularSystem(h5mdio::GmxH5mdIo*       file,
                          const gmx_mtop_t&        topology,
                          gmx::ArrayRef<const int> index            = {},
                          const std::string        index_group_name = "");

/*! \brief Write a trajectory frame to the file. Only writes the data that is passed as input
 *
 * \param[in] file The H5MD file manager to use.
 * \param[in] step The simulation step.
 * \param[in] time The time stamp (in ps).
 * \param[in] lambda The lambda state.
 * \param[in] box The box dimensions.
 * \param[in] numParticles The number of particles in the system.
 * \param[in] x The particle coordinates for lossless output.
 * \param[in] v The particle velocities for lossless output.
 * \param[in] f The particle forces for lossless output.
 * \param[in] xCompressionError The accepted error of the lossy compression of coordinates.
 * \throws FileIOError    If there is no file open or if errors occured during writing.
 */
void writeFrame(h5mdio::GmxH5mdIo* file,
                int64_t            step,
                real               time,
                real               lambda,
                const rvec*        box,
                int64_t            numParticles,
                const rvec*        x,
                const rvec*        v,
                const rvec*        f,
                double             xCompressionError);


} // namespace gmx
#endif // GMX_FILEIO_H5MD_IO_H
