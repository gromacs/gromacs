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

#ifndef GMX_FILEIO_GMXH5MD_IO_H
#define GMX_FILEIO_GMXH5MD_IO_H

#include <string>
#include <vector>

#include "gromacs/math/vectypes.h"
#include "gromacs/topology/atoms.h"
#include "gromacs/utility/arrayref.h"
#include "gromacs/utility/real.h"

#include "h5md_datablock.h"
#include "h5md_io.h"
#include "h5md_util.h"

enum class PbcType : int;
struct gmx_moltype_t;
struct gmx_mtop_t;
struct MoleculeBlockIndices;

namespace gmx
{
namespace h5mdio
{
class GmxH5mdIo;
}
/*! \brief Set the author and creator fields based on user and program data available in gromacs.

 * \param[in] file     The H5MD file manager to use.
 */
void setH5mdAuthorAndCreator(h5mdio::GmxH5mdIo* file);

/*! \brief Setup molecular system particle related data to the file. Write the data that is present from the beginning.
 *
 * This is currently not updated during the trajectory. The data that is written are atom masses, atom charges, atom
 * species and IDs.
 *
 * \param[in] file     The H5MD file manager to use.
 * \param[in] topology The molecular topology describing the system.
 * \param[in] index    The selected atoms to include. If empty, use all atoms in the topology.
 * \param[in] selectionName The name of the atom selection specified by index.
 * \throws FileIOError If there is no file open or if the data could not be written.
 */
void setupMolecularSystemParticleData(h5mdio::GmxH5mdIo*       file,
                                      const gmx_mtop_t&        topology,
                                      gmx::ArrayRef<const int> index         = {},
                                      const std::string        selectionName = "");

/*! \brief Get a specific GROMACS molecule block indices data structure (corresponding to a number
 * of molecule of a certain molecule type) from the H5MD GROMACS topology section in the file.
 *
 * \param[in] file               The H5MD file manager to use.
 * \param[in] molBlockIndex      The index of the molecule block.
 * \returns the MoleculeBlockIndices struct. If the index does not match a valid molecule block it
 * will be empty.
 * \throws FileIOError If there was en error reading the data.
 */
MoleculeBlockIndices getMoleculeBlockIndicesByIndex(h5mdio::GmxH5mdIo* file, size_t molBlockIndex);

/*! \brief Setup molecular system topology data.
 *
 * \param[in] file     The H5MD file manager to use.
 * \param[in] topology The molecular topology describing the system.
 * \param[in] index    The selected atoms to include. If empty, use all atoms in the topology.
 * \param[in] selectionName The name of the atom selection specified by index.
 * \param[in] abortIfPresent Do not set up the topology if it is already present in the file.
 * \param[in] writeVmdStructureData Whether to write information describing the system according
 * to VMD specifications in the H5MD group /parameters/vmd_structure/.
 * \throws FileIOError If there is no file open or if the data could not be written.
 */
void setupMolecularSystemTopology(h5mdio::GmxH5mdIo*       file,
                                  const gmx_mtop_t&        topology,
                                  gmx::ArrayRef<const int> index                 = {},
                                  bool                     abortIfPresent        = true,
                                  bool                     writeVmdStructureData = false);

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
 * \param[in] selectionName The name of the selection to use.
 * \throws FileIOError    If there is no file open or if errors occured during writing.
 */
void writeFrameToStandardDataBlocks(h5mdio::GmxH5mdIo* file,
                                    int64_t            step,
                                    real               time,
                                    real               lambda,
                                    const rvec*        box,
                                    int64_t            numParticles,
                                    const rvec*        x,
                                    const rvec*        v,
                                    const rvec*        f,
                                    double             xCompressionError,
                                    const std::string  selectionName = "system");

/*! \brief Read the next frame of box, coordinates, velocities and forces. With next frame means
 * the lowest step/time reading from the previous read frame of that data type. If data is
 * written at different intervals the read data types will be different from one function call
 * to the next.
 *
 * \param[in] file The H5MD file manager to use.
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
 * \param[in] selectionName The name of the selection to use.
 * \returns Whether any frame was read or not.
 * \throws FileIOError If there was an error reading the next frame or if the data type of the data was unknown.
 */
bool readNextFrameOfStandardDataBlocks(h5mdio::GmxH5mdIo* file,
                                       int64_t*           step,
                                       real*              time,
                                       real*              lambda,
                                       rvec*              box,
                                       rvec*              x,
                                       rvec*              v,
                                       rvec*              f,
                                       real*              xCompressionError,
                                       bool*              readLambda,
                                       bool*              readBox,
                                       bool*              readX,
                                       bool*              readV,
                                       bool*              readF,
                                       const std::string  selectionName = "system");

/*! \brief Copies the provenance group (in /modules) from \p srcFile to \p destFile, if it exists in \p srcFile.
 *
 * \param[in] srcFile The H5MD file manager of the source file.
 * \param[in] destFile The H5MD file manager of the destination file.
 * \returns Whether the group was copied.
 */
bool copyProvenanceRecords(h5mdio::GmxH5mdIo* srcFile, h5mdio::GmxH5mdIo* destFile);

} // namespace gmx
#endif // GMX_FILEIO_GMXH5MD_IO_H
