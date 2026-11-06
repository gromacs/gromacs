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

/*! \brief Declarations of H5md particle data block routines.
 *
 * \author Petter Johansson <pettjoha@kth.se>
 */

#ifndef GMX_FILEIO_H5MD_PARTICLEBLOCK_H
#define GMX_FILEIO_H5MD_PARTICLEBLOCK_H

#include <hdf5.h>

#include <optional>

#include "gromacs/fileio/h5md/h5md_timedatablock.h"
#include "gromacs/utility/vectypes.h"

namespace gmx
{

/*! \brief Class which manages trajectory data inside of a H5md particle group.
 *
 * Per the H5md specification, trajectory data for subgroups is stored in
 * the tree hierarchy as:
 *
 * \-- /particles/group1
 *     \-- box
 *         \-- edges
 *             \-- value: real[numFrames][3][3]
 *     \-- position
 *         \-- value:  RVec[numFrames][numAtoms]
 *         \-- step:   int64_t[numFrames]
 *         \-- (time): double[numFrames] (optional)
 *     \-- velocity
 *         \-- value:  RVec[numFrames][numAtoms]
 *         \-- step:   int64_t[numFrames]
 *         \-- (time): double[numFrames] (optional)
 *     \-- force
 *         \-- value:  RVec[numFrames][numAtoms]
 *         \-- step:   int64_t[numFrames]
 *         \-- (time): double[numFrames] (optional)
 *
 * This class collects data sets for a single group into one object, providing
 * access to \c H5mdTimeDataBlocks<RVec> for position, velocity and force using getters,
 * and a \c H5mdFrameDataSet<real> for the simulation box.
 */
class H5mdParticleBlock
{
public:
    //! \brief Constructor for particle block containing \p numParticles particles.
    H5mdParticleBlock(int64_t numParticles) : numParticles_{ numParticles } {}

    /*! \brief Add and manage a data block for positions.
     *
     * \throws H5mdError if the block has already been set or if \p position
     * does not store \c numParticles_ particles.
     */
    void addPosition(H5mdTimeDataBlock<RVec>&& position);

    /*! \brief Add and manage a data block for velocities.
     *
     * \throws H5mdError if the block has already been set or if \p velocity
     * does not store \c numParticles_ particles.
     */
    void addVelocity(H5mdTimeDataBlock<RVec>&& velocity);

    /*! \brief Add and manage a data block for forces.
     *
     * \throws H5mdError if the block has already been set or if \p force
     * does not store \c numParticles_ particles.
     */
    void addForce(H5mdTimeDataBlock<RVec>&& force);

    /*! \brief Add and manage a data block for the simulation box.
     *
     * \throws H5mdError if the block has already been set or if its
     * data set frame dimensions are not { DIM, DIM }.
     */
    void addBox(H5mdFrameDataSet<real>&& box);

    //! \brief Return whether the block stores a per-position-frame box value, else false.
    bool hasBox() const { return box_.has_value(); }

    //! \brief Return whether the block stores a position block.
    bool hasPosition() const { return position_.has_value(); }

    //! \brief Return whether the block stores a velocity block.
    bool hasVelocity() const { return velocity_.has_value(); }

    //! \brief Return whether the block stores a force block.
    bool hasForce() const { return force_.has_value(); }

    //! \brief Return the number of particles for the managed blocks (0 if none are managed).
    int64_t numParticles() const { return numParticles_; }

    //! \brief Get the simulation box data set.
    std::optional<H5mdFrameDataSet<real>>& box() { return box_; }

    //! \brief Get the position data set.
    std::optional<H5mdTimeDataBlock<RVec>>& position() { return position_; }

    //! \brief Get the velocity data set.
    std::optional<H5mdTimeDataBlock<RVec>>& velocity() { return velocity_; }

    //! \brief Get the force data set.
    std::optional<H5mdTimeDataBlock<RVec>>& force() { return force_; }

    //! \brief Trim managed data sets to a maximum of \p maxStep - 1 steps.
    void trimDataSetsFromStep(int64_t maxStep);

private:
    //! \brief Positions of particles in block (if written to disk)
    std::optional<H5mdTimeDataBlock<RVec>> position_;

    //! \brief Velocities of particles in block (if written to disk)
    std::optional<H5mdTimeDataBlock<RVec>> velocity_;

    //! \brief Forces of particles in block (if written to disk)
    std::optional<H5mdTimeDataBlock<RVec>> force_;

    //! \brief Simulation box (if non-constant over the simulation)
    std::optional<H5mdFrameDataSet<real>> box_;

    //! \brief Number of particles for managed data sets.
    const int64_t numParticles_;
};

} // namespace gmx

#endif // GMX_FILEIO_H5MD_PARTICLEBLOCK_H
