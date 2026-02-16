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

/*! \libinternal
 * \brief Declares wrappers to H5MD class.
 *
 * This avoids clashes between symbols in HDF5 headers
 * and GROMACS headers. */

#ifndef GMX_FILEIO_H5MD_H5MD_WRAPPER_H
#define GMX_FILEIO_H5MD_H5MD_WRAPPER_H

#include <filesystem>

namespace gmx
{

class H5md;

enum class H5mdFileMode : char;

/*! \brief Builder method to isolate creation from callers
 *
 * This should replaced by an impl class to hide the dependency
 * on hdf5 types. */
H5md* makeH5md(const std::filesystem::path& fileName, H5mdFileMode mode);

/*! \brief Set up the file from input data.
 *
 * \param[out] h5md       Handle to H5md object
 * \param[in] topology    Topology data of system.
 * \param[in] inputRecord Simulation parameters.
 */
void setupFileFromInput(H5md* h5md, const gmx_mtop_t& topology, const t_inputrec& inputRecord);

/*! \brief Set up from an existing \p h5md file for restarting with appending.
 *
 * Scans for trajectory data available in the /particles/system group
 * of the HDF5 file. For trajectories written by GROMACS this group
 * contains data for all atoms in the simulated system.
 *
 * \note Ignores trajectory data in other subgroups of /particles.
 *
 * \param[out] h5md              Handle to H5md object
 * \param[in] restartingFromStep Step which we are restarting the simulation from
 * \param[in]  numParticles      Number of particles in system
 *
 * \throws gmx::FileIOError if the existing file contents are not
 * for a system of \p numParticles.
 */
void setupFromExistingFileForAppending(H5md* h5md, int64_t restartingFromStep, int64_t numParticles);

/*! \brief Write input data as the next frame of the trajectory.
 *
 * \param[out] h5md         Handle to H5md object
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
void writeNextFrame(H5md*                h5md,
                    ArrayRef<const RVec> positions,
                    ArrayRef<const RVec> velocities,
                    ArrayRef<const RVec> forces,
                    const matrix         box,
                    int64_t              step,
                    double               time);

//! Deallocate \c h5md
void destroyH5md(H5md* h5md);

//! Flush \c h5md file.
void flushH5md(H5md* h5md);

} // namespace gmx

#endif
