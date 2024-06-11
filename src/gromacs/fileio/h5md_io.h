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
#include "gromacs/topology/atoms.h"
#include "gromacs/utility/arrayref.h"
#include "gromacs/utility/real.h"

#include "h5md_datablock.h"
#include "h5md_util.h"

enum class PbcType : int;
struct gmx_moltype_t;
struct gmx_mtop_t;
struct MoleculeBlockIndices;

namespace gmx
{
namespace h5mdio
{

typedef int64_t            hid_t;
typedef unsigned long long hsize_t;
constexpr size_t           c_atomResidueStringLen               = 20;
constexpr size_t           c_moleculeTypeStringLen              = 256;
constexpr size_t           c_provenanceRecordStringLen          = 512;
constexpr int              c_maxFullNameLength                  = 512;
constexpr int              c_h5mdMajorVersion                   = 1;
constexpr int              c_h5mdMinorVersion                   = 1;
constexpr int              c_gmxH5mdParametersGroupMajorVersion = 0;
constexpr int              c_gmxH5mdParametersGroupMinorVersion = 9;
constexpr int              c_gmxH5mdProvenanceGroupMajorVersion = 0;
constexpr int              c_gmxH5mdProvenanceGroupMinorVersion = 9;
static std::string         s_gromacsTopologyGroupName           = "/parameters/gromacs_topology";
static std::string         s_provenanceGroupName                = "/modules/provenance";


/*! \brief The container of the H5MD data. The class is designed to read/write data according to de Buyl et al., 2014
 * (https://www.sciencedirect.com/science/article/pii/S0010465514000447) and https://www.nongnu.org/h5md/h5md.html
 * The class contains a number of standard data blocks that are commonly used by GROMACS. */
class GmxH5mdIo
{
private:
    hid_t file_; //!< The HDF5 identifier of the file. This is the H5MD root.
    std::list<GmxH5mdTimeDataBlock> dataBlocks_; //!< A list of time dependent data blocks in the HDF5 file.
    char filemode_; //!< Whether the file is open for reading ('r'), writing ('w') or appending ('a')

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

    /*! \brief Return the version number of the H5MD root group. The version number is returned
     * as a string, composed as "<majorVersion>.<minorVersion>". If not version number could be
     * read, the returned string is empty. */
    std::string getH5mdRootVersionNumber();

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
     * \returns The creator name, i.e. the name of the program that created the file.
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

    /*! \brief Get the HDF5 ID of the group specified by \p fullName (by absolute path in the file).
     * \returns The ID of the group or < 0 if the group does not exist.
     */
    hid_t getGroupId(const std::string& fullName);

    /*! \brief Create an HDF5 group specified by \p fullName (by absolute path in the file).
     * \returns The ID of the group or < 0 if the group could not be created.
     */
    hid_t createGroup(const std::string& fullName);

    /*! \brief Create an HDF5 group in \p container specified by \p nameInContainer.
     * \returns The ID of the group or < 0 if the group could not be created.
     */
    hid_t createGroup(hid_t container, const std::string& nameInContainer);

    /*! \brief Set multiple values of a string data set. All entries in the propertyValues vector are set,
     * starting at the beginning of the dataset.
     * N.b., The maximum string length is set to c_atomResidueStringLen to avoid the overhead that comes
     * with flexible string lengths in HDF5.
     *
     * \param[in] containerName The name of the HDF5 group that contains the data set, e.g., "/particles/system".
     * \param[in] dataSetName The name of the property to set, e.g., "atomname".
     * \param[in] propertyValues A vector of strings.
     * \param[in] replaceExisting Whether to replace the property if it already exists.
     * \param[in] maxStringLength Is the length if fixed-length strings. If 0, the string will be
     * written as variable-length. Variable-length strings cannot be compressed, and are therefore
     * recommended against.
     * \throws FileIOError If there was an error creating the data block or writing the data.
     */
    void setStringDataSet(const std::string&              containerName,
                          const std::string&              dataSetName,
                          const std::vector<std::string>& propertyValues,
                          bool                            replaceExisting = false,
                          size_t                          maxStringLength = c_atomResidueStringLen);

    /*! \brief Set multiple values of a numeric data set. All entries in the propertyValues vector are set.
     *
     * \param[in] containerName The name of the HDF5 group that contains the data set, e.g., "/particles/system".
     * \param[in] dataSetName The name of the property to set, e.g., "position".
     * \param[in] propertyValues A vector of numeric values.
     * \param[in] replaceExisting Whether to replace the property if it already exists.
     * \throws FileIOError If there was an error creating the data block or writing the data.
     */
    template<typename T>
    void setNumericDataSet(const std::string&    containerName,
                           const std::string&    dataSetName,
                           const std::vector<T>& propertyValues,
                           const std::string&    unit,
                           bool                  replaceExisting = false);


    /*! \brief Read the contents of a string data set. The string can be either fixed-length or variable-length.
     *
     * \param[in] containerName The name of the HDF5 group that contains the data set, e.g., "/particles/system".
     * \param[in] dataSetName The name of the property to read, e.g., "atomname".
     * \returns A vector containing the strings. Empty if no strings could be read.
     * \throws FileIOError If there was an error reading the data.
     */
    std::vector<std::string> readStringDataSet(const std::string& containerName,
                                               const std::string& dataSetName);

    /*! \brief Read a numeric property.
     *
     * \param[in] containerName The name of the HDF5 group that contains the data set, e.g.,
     * "/particles/system".
     * \param[in] dataSetName The name of the property to read, e.g., "position".
     * \returns A vector containing the numeric values. Empty if no data could be read.
     * \throws FileIOError If there was an error reading the data.
     */
    template<typename T>
    std::vector<T> readNumericDataSet(const std::string& containerName, const std::string& dataSetName);

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

    /*! \brief Reads the next frame of the data block. If stepToRead >= 0 the step of the next frame
     * must match stepToRead.
     *
     * \param[in] dataBlockFullName The full path (in the H5MD file) to the data block group.
     * \param[in] data Then container in which the read data will be stored.
     * \param[in] stepToRead If >0 the step of the next frame must match this. Otherwise no frame
     *                       will be read.
     * \returns Whether a frame was read or not.
     */
    bool readNextFrameOfDataBlock(std::string dataBlockFullName, real* data, int64_t stepToRead = -1);

    /* \brief Returns the lossy compression error of the data block matching the name.
     *        Returns -1 if there is no lossy compression or if the data block was not found.
     */
    real getLossyCompressionErrorOfDataBlock(std::string dataBlockFullName);

    /*! \brief Get the number of frames of a particles data block
     * \param[in] dataBlockName The name of the data block.
     * \param[in] selectionName The name of the selection to use.
     * \returns the number of frames of the data block or -1 if not found.
     */
    int64_t getNumberOfFrames(const std::string dataBlockName,
                              const std::string selectionName = "system");

    /*! \brief Get the number of particles of a particles data block
     * \param[in] dataBlockName The name of the data block.
     * \param[in] selectionName The name of the selection to use.
     * \returns the number of particles of the data block or -1 if not found.
     * \throws FileIOError    If there was an error reading the number of particles could not be
     * read, such as no atom data in particles data blocks.
     */
    int64_t getNumberOfParticles(const std::string dataBlockName,
                                 const std::string selectionName = "system");

    /*! \brief Get the first time stamp of a particles data block
     * \param[in] dataBlockName The name of the data block.
     * \param[in] selectionName The name of the selection to use.
     * \returns the first time of the data block or -1 if not found.
     * \throws FileIOError    If there was an error determining the time of the first frame.
     */
    real getFirstTime(const std::string dataBlockName, const std::string selectionName = "system");

    /*! \brief Get the very first time stamp of all particles data blocks
     * \returns the first time of all data blocks or -1 if no data blocks were found.
     * \throws FileIOError    If there was an error determining the time of the first frame.
     */
    real getFirstTimeFromAllDataBlocks();

    /* \brief Looks through all data blocks to find which are the frames and steps of their
     * next frame to read.
     * \returns The lowest step and its corresponding time in a tuple.
     */
    std::tuple<int64_t, real> getNextStepAndTimeToRead();

    /*! \brief Get the final time stamp of a particles data block
     * \param[in] dataBlockName The name of the data block.
     * \param[in] selectionName The name of the selection to use.
     * \returns the final time of the data block or -1 if not found.
     * \throws FileIOError    If there was an error determining the time of the final frame.
     */
    real getFinalTime(const std::string dataBlockName, const std::string selectionName = "system");

    /*! \brief Get the very last time stamp of all particles data blocks
     * \returns the last time of all data blocks or -1 if no data blocks were found.
     * \throws FileIOError    If there was an error determining the time of the final frame.
     */
    real getFinalTimeFromAllDataBlocks();

    /*! \brief Add an entry to the provenance record.
     *
     * \param[in] commandLine The command line used when opening the file (for writing).
     * \param[in] version The version of the program used when opening the file.
     * \param[in] comment A comment to this provenance entry, e.g., "open file", "flush", "close file" etc.
     */
    void addToProvenanceRecord(const std::string& commandLine = "",
                               const std::string& version     = "",
                               const std::string& comment     = "");
};

} // namespace h5mdio

} // namespace gmx
#endif // GMX_FILEIO_H5MD_IO_H
