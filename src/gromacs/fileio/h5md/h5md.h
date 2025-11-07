/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2024- The GROMACS Authors
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

/*! \brief Declares the i/o interface to H5MD HDF5 files.
 *
 * \author Petter Johansson <pettjoha@kth.se>
 * \author Magnus Lundborg <lundborg.magnus@gmail.com>
 */

#ifndef GMX_FILEIO_H5MD_H
#define GMX_FILEIO_H5MD_H

#include "config.h" // To define GMX_USE_HDF5

#if GMX_USE_HDF5
#    include <hdf5.h>

#    include "h5md_particleblock.h"
#endif

#include <filesystem>
#include <optional>
#include <string>
#include <unordered_map>

#include "gromacs/utility/vectypes.h"

enum class PbcType : int;

struct gmx_mtop_t;
struct t_atoms;
struct t_inputrec;

namespace gmx
{

template<typename T>
class ArrayRef;
class MessageStringCollector;

#if !GMX_USE_HDF5
// Backup typedefs to be used when we are not compiling with the HDF5 library enabled.
// These are needed to ensure that the externally exposed interface from the module
// is available.
// TODO: These can be removed once the module implementation is fully separated
// from the exposed interface.
using hid_t = int64_t;
#endif

enum class H5mdFileMode : char
{
    Read  = 'r', //! Only read from the file.
    Write = 'w', //! Write to the file, replaces it if it exists. Also allows reading from the file.
    Append = 'a' //! Write to the file without replacing it if it exists. Also allows reading from the file.
};

/*! \brief Manager of an H5MD filehandle.
 * The class is designed to read/write data according to de Buyl et al., 2014
 * (https://doi.org/10.1016/j.cpc.2014.01.018) and https://www.nongnu.org/h5md/h5md.html
 */
class H5md
{
#if GMX_USE_HDF5
private:
    hid_t file_;            //!< The HDF5 identifier of the file. This is the H5MD root.
    H5mdFileMode filemode_; //!< Whether the file is open for reading ('r'), writing ('w') or appending ('a')

    /*! \brief Cursor for reading trajectory data from a particle block.
     *
     * Writing trajectory data always appends frames to the respective data sets. When reading
     * we want more up front control over the index to facilitate things like seeking. This
     * cursor will faciliate this behavior.
     *
     * TODO: This class is only sketched out for now. Expand functionality as needed
     * when implementing reading.
     */
    class TrajectoryReadCursor
    {
    public:
        //! \brief Constructor.
        TrajectoryReadCursor(H5mdParticleBlock&& block) : block_{ std::move(block) } {}
        //! \brief Return the particle block.
        H5mdParticleBlock& block() { return block_; }

    private:
        //! \brief Particle block.
        H5mdParticleBlock block_;
    };

    //! \brief List of particle blocks in /particles/, with the key being the block name
    std::unordered_map<std::string, TrajectoryReadCursor> particleBlocks_;
#endif

public:
    /*! \brief Open an H5MD file and manage its filehandle.
     *
     * \param[in] fileName    Name of the file to open. The same as the file path.
     * \param[in] mode        The mode to open the file.
     * \throws FileIOError if fileName is specified and the file cannot be opened.
     */
    H5md(const std::filesystem::path& fileName, H5mdFileMode mode);

    ~H5md();

    H5md(const H5md&)            = delete;
    H5md& operator=(const H5md&) = delete;
    H5md(H5md&&)                 = delete;
    H5md& operator=(H5md&&)      = delete;

    /*! \brief Return the HDF5 handle of the file.
     *
     * This handle is used as the root container of the file's data hierarchy, to which
     * the tree of groups and data sets can be added.
     *
     *  \returns The file handle.
     */
    hid_t fileid() const;

    /*! \brief Write all unwritten data to the file.
     * \param[in] throwExceptionUponError Whether to throw an exception if an error occurs.
     * Assumes a valid file_ identifier.
     * \throws FileIOError If there were errors during flushing (and throwExceptionUponError is true).
     */
    void flush(bool throwExceptionUponError = true);

    /*! \brief Set the author name attribute in the H5MD file
     * \param[in] authorName The author name.
     * \throws FileIOError If the author name attribute could not be set.
     */
    void setAuthor(const std::string& authorName);

    /*! \brief Get the author name attribute from the H5MD file
     * \returns the author name if the attribute was set.
     */
    std::optional<std::string> author();

    /*! \brief Set the name of the creating program as an attribute in the H5MD file.
     * \param[in] creatorName The creator name, i.e. the name of the program that created the file.
     * \throws FileIOError If the creator name attribute could not be set.
     */
    void setCreatorProgramName(const std::string& creatorName);

    /*! \brief Get the name of the creating program attribute from the H5MD file.
     * \returns The creator name, i.e. the name of the program that created the file,
     * if the attribute was set.
     */
    std::optional<std::string> creatorProgramName();

    /*! \brief Set the version of the creating program as an attribute in the H5MD file.
     * \param[in] version The version.
     * \throws FileIOError If the version attribute could not be set.
     */
    void setCreatorProgramVersion(const std::string& version);

    /*! \brief Get the version of the creating program attribute from the H5MD file.
     * \returns the version if the attribute was set.
     */
    std::optional<std::string> creatorProgramVersion();

    /*! \brief Set up the file from input data.
     *
     * \param[in] topology    Topology data of system.
     * \param[in] inputRecord Simulation parameters.
     */
    void setupFileFromInput(const gmx_mtop_t& topology, const t_inputrec& inputRecord);

    /*! \brief Set up from an existing file.
     *
     * Scans for trajectory data available in the /particles/system group
     * of the HDF5 file. For trajectories written by GROMACS this group
     * contains data for all atoms in the simulated system.
     *
     * \note Ignores trajectory data in other subgroups of /particles.
     */
    void setupFromExistingFile();

    /*! \brief Write input data as the next frame of the trajectory.
     *
     * \param[in] positions     Position data to write (or empty if not to write).
     * \param[in] velocities    Velocity data to write (or empty if not to write).
     * \param[in] forces        Force data to write (or empty if not to write).
     * \param[in] box           Simulation box for frame.
     * \param[in] step          Simulation step for frame.
     * \param[in] time          Simulation time for frame.
     *
     * \throws gmx::FileIOError if \p position, \p velocity or \p force data is given
     *     but the corresponding data set has not been created, or if the size of
     *     the data buffers do not match the number of atoms of the system.
     */
    void writeNextFrame(ArrayRef<const RVec> positions,
                        ArrayRef<const RVec> velocities,
                        ArrayRef<const RVec> forces,
                        const matrix         box,
                        int64_t              step,
                        double               time);

private:
    /*! \brief Set up particle blocks for trajectory writing according to the H5md specification.
     *
     * \param[in] topology         Molecular topology for the simulated system.
     * \param[in] selectionIndices Global indices for all atoms in the group.
     * \param[in] selectionName    Name of group.
     * \param[in] inputRecord      Simulation input record.
     */
    void setupParticleBlockForGroup(const gmx_mtop_t&   topology,
                                    ArrayRef<const int> selectionIndices,
                                    const std::string&  selectionName,
                                    const t_inputrec&   inputRecord);

    /*! \brief Set up particle block with name \p selectionName from the existing file structure.
     *
     * \param[in]  selectionName Name of group to set up particle block for
     */
    void setupParticleBlockForGroupFromExistingFile(const std::string& selectionName);

    /*! \brief Set up the H5MD metadata group according to the H5md specification.
     *
     * \note Also sets up the modules group with GROMACS and units modules. This
     * method must be called before adding any other H5MD module information.
     *
     * Creates the following groups and attributes (relative to the HDF5 root):
     *
     * /h5md                (H5md metadata group)
     * \++ version          (H5md specification version attribute)
     * \-- /author          (author group)
     *     \++ name         (current username, or a dummy name if we cannot get it)
     * \-- /creator         (creator group)
     *     \++ name         (GROMACS)
     *     \++ version      (GROMACS version)
     * \-- /modules         (modules group)
     *     \-- /gromacs     (GROMACS module group)
     *         \++ version  (GROMACS module version)
     *     \-- /units       (units module group)
     *         \++ version  (units module version)
     *         \++ system   (unit system string)
     *
     * These groups and attributes are required by the H5MD specification to exist
     * in the file. After setting these up we can proceed with the rest of the file.
     */
    void setupMetadataGroup();
};

} // namespace gmx
#endif // GMX_FILEIO_H5MD_H
