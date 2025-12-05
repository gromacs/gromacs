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
 * \author Yang Zhang <yang.zhang@scilifelab.se>
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
struct t_trxframe;

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
        TrajectoryReadCursor(H5mdParticleBlock&& block) :
            block_{ std::move(block) },
            nextPositionFrameToRead_{ 0 },
            nextVelocityFrameToRead_{ 0 },
            nextForceFrameToRead_{ 0 }
        {
        }

        //! \brief Return the particle block.
        H5mdParticleBlock& block() { return block_; }

        /*! \brief Check which data blocks are available for the next frame and return whether a next frame exists.
         *
         * \param[out] hasPosition Whether or not position data exist for the next frame.
         * \param[out] hasVelocity Whether or not velocity data exist for the next frame.
         * \param[out] hasForce    Whether or not force data exist for the next frame.
         * \param[out] hasBox      Whether or not simulation box data exist for the next frame.
         * \param[out] hasStep     Whether or not step data exist for the next frame.
         * \param[out] hasTime     Whether or not time data exist for the next frame.
         *
         * \returns True if a next frame exists to read, otherwise false.
         */

        bool nextFrameContents(bool* hasPosition,
                               bool* hasVelocity,
                               bool* hasForce,
                               bool* hasBox,
                               bool* hasStep,
                               bool* hasTime);

        /*! \brief Read the next frame into given data buffers.
         *
         * \note This attempts to read the next frame for all non-empty input data, without checking
         * that they are stored at the same simulation step and time. To read the data in sync,
         * first make a call to \c nextFrameContents, which sets flags corresponding to the the data
         * blocks found in the next frame. Then use those flags to select the output
         * blocks in a call to this function.
         *
         * \param[out] positions     Buffer to read position data into (or empty if not to read).
         * \param[out] velocities    Buffer to read velocity data into (or empty if not to read).
         * \param[out] forces        Buffer to read force data into (or empty if not to read).
         * \param[out] box           Buffer to read simulation box data into.
         * \param[out] step          Buffer to read step into.
         * \param[out] time          Buffer to read simulation time into.
         *
         * \returns True if a frame was read, otherwise false.
         *
         * \throws gmx::FileIOError if \p position, \p velocity or \p force data is given
         *     but the corresponding data set has not been created, or if the size of
         *     the data buffers do not match the number of atoms of the system.
         */
        bool readNextFrame(ArrayRef<RVec> position,
                           ArrayRef<RVec> velocity,
                           ArrayRef<RVec> force,
                           matrix         box,
                           int64_t*       step,
                           double*        time);

    private:
        //! \brief Particle block.
        H5mdParticleBlock block_;

        int64_t nextPositionFrameToRead_;
        int64_t nextVelocityFrameToRead_;
        int64_t nextForceFrameToRead_;
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

    /*! \brief Set up from an existing file when restarting a simulation from a checkpoint.
     *
     * Scans for trajectory data available in the /particles/system group
     * of the HDF5 file. For trajectories written by GROMACS this group
     * contains data for all atoms in the simulated system.
     *
     * \note Ignores trajectory data in other subgroups of /particles.
     *
     * \param[in] restartingFromStep Step which we are restarting the simulation from
     * \param[in] numParticles       Number of particles in system
     *
     * \throws gmx::FileIOError if the existing file contents are not
     * for a system of \p numParticles.
     */
    void setupFromExistingFileForAppending(int64_t restartingFromStep, int64_t numParticles);

    /*! \brief Read the next \p frame of the trajectory.
     *
     * \param[out] frame         Container to read data into.
     * \param[in]  selectionName Name of group to read frame for.
     *
     * \returns True if a frame was read, otherwise false.
     */
    bool readNextFrame(t_trxframe* frame, const std::string& selectionName = "system");

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
     * \param[in] selectionName        Name of group to set up particle block for
     * \param[in] restartingFromStep   (Optional) If restarting simulation, the step that we are
     *                                 restarting from.
     * \param[in] expectedNumParticles (Optional) Expected number of particles of block in file.
     *
     * \throws gmx::FileIOError if \p expectedNumParticles is given and the contents of the existing
     * particle block does not match it.
     */
    void setupParticleBlockForGroupFromExistingFile(const std::string&     selectionName,
                                                    std::optional<int64_t> restartingFromStep,
                                                    std::optional<int64_t> expectedNumParticles);

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

    /*! \brief Set up the module for GROMACS molecular topology
     *
     * This module store the GROMACS molecular topology for internal use.
     * The version of this module is determined by the h5md_topologyutils module.
     *
     * Create the following hierarchy (relative to the HDF5 root):
     *
     * /h5md/modules/gromacs_topology  (group for H5md GROMACS topology module)
     * \++ version                     (attribute for the version of GROMACS topology)
     * \++ molecule_names              (attribute for the names of molecules in the system)
     * \++ system_name                 (attribute for the name of the system)
     * \-- molecule1                   (group for molecule 1, same group name as molecule_names[0])
     *     \++ nr_particles            (attribute for the number of particles in molecule 1)
     *     \++ nr_residues             (attribute for the number of residues in molecule 1)
     *     \++ nr_blocks               (attribute for the number of blocks for molecule 1)
     *     \-- id                      (dataset for the atomic identifier in molecule 1)
     *     \-- mass                    (dataset for the atomic masses in molecule 1)
     *     \-- charge                  (dataset for the atomic charges in molecule 1)
     *     \-- species                 (dataset for the atomic species in molecule 1)
     *     \-- particle_name           (dataset for the atomic names (indices into particle name table) in molecule 1)
     *     \-- particle_name_table     (dataset for the atomic name lookup table in molecule 1)
     *     \-- residue_id              (dataset for the residue identifier in molecule 1)
     *     \-- sequence                (dataset for the residue sequence (indices into residue name table) in molecule 1)
     *     \-- residue_name            (dataset for the residue names (indices into residue name table) in molecule 1)
     *     \-- residue_name_table      (dataset for the residue name lookup table in molecule 1)
     *     \-- posres_xA               (dataset for the position restraint for topology A values in molecule 1)
     *     \-- posres_xB               (dataset for the position restraint for topology B values in molecule 1)
     *
     * \param[in] topology Molecular topology for the simulated system.
     */
    void setupGromacsTopology(const gmx_mtop_t& topology);

    /*! \brief Set up the H5MD metadata group for connectivity information.
     *
     * Create the following hierarchy (relative to the HDF5 root):
     *
     * /connectivity          (group for connectivity information)
     * \++ nr_bonds           (attribute for the number of bonds in the system)
     *     \-- bonds          (dataset for the bonds (shaped [numBond, 2]) in the system)
     *
     * \param[in] topology Molecular topology for the simulated system.
     */
    void setupBondConnectivity(const gmx_mtop_t& topology);
};

} // namespace gmx
#endif // GMX_FILEIO_H5MD_H
