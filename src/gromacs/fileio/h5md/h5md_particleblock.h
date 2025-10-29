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

class H5mdParticleBlockBuilder;

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
 *
 * Use the \c H5mdParticleBlockBuilder builder class to construct this object
 * by adding selected data sets at will.
 */
class H5mdParticleBlock
{
public:
    // Only make constructors available to our builder class.
    friend class H5mdParticleBlockBuilder;

    //! \brief Return whether the block stores a per-position-frame box value, else false.
    bool hasBox() const { return box_.has_value(); }

    //! \brief Return whether the block stores a position block.
    bool hasPosition() const { return position_.has_value(); }

    //! \brief Return whether the block stores a velocity block.
    bool hasVelocity() const { return velocity_.has_value(); }

    //! \brief Return whether the block stores a force block.
    bool hasForce() const { return force_.has_value(); }

    //! \brief Return whether time is stored in any of the managed data block sets.
    bool hasTime() const { return hasTime_; }

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

private:
    //! \brief Constructor.
    H5mdParticleBlock(std::optional<H5mdTimeDataBlock<RVec>>&& position,
                      std::optional<H5mdTimeDataBlock<RVec>>&& velocity,
                      std::optional<H5mdTimeDataBlock<RVec>>&& force,
                      std::optional<H5mdFrameDataSet<real>>&&  box,
                      int64_t                                  numParticles);

    //! \brief Positions of particles in block (if written to disk)
    std::optional<H5mdTimeDataBlock<RVec>> position_;

    //! \brief Velocities of particles in block (if written to disk)
    std::optional<H5mdTimeDataBlock<RVec>> velocity_;

    //! \brief Forces of particles in block (if written to disk)
    std::optional<H5mdTimeDataBlock<RVec>> force_;

    //! \brief Simulation box (if non-constant over the simulation)
    std::optional<H5mdFrameDataSet<real>> box_;

    //! \brief Whether or not any managed data set stores times.
    bool hasTime_;

    //! \brief Number of particles for managed data sets.
    int64_t numParticles_;
};

/*! \brief Builder class for \c H5mdParticleBlock objects.
 *
 * Asserts that all added data sets are consistent with our expectation
 * for trajectory data. Trajectory data sets must have a 1-dimensional
 * frame dimension (i.e. each frame has shape real[numParticles][3])
 * and created for the same number of particles. Adding a \c H5mdTimeDataBlock
 * which has inconsistent frame dimensions and number of particles results
 * in a throw.
 *
 * If no data sets are added, the constructed \c H5mdParticleBlock will
 * have numParticles() = 0. In this case no data sets are managed.
 */
class H5mdParticleBlockBuilder
{
public:
    //! Construct a new builder.
    H5mdParticleBlockBuilder() = default;

    /*! \brief Assign a per-frame box data set.
     *
     * \throws gmx::FileIOError if the data set frame dimensions are not { 3, 3 }.
     */
    H5mdParticleBlockBuilder& setBox(H5mdFrameDataSet<real>&& box);

    /*! \brief Assign a position data set.
     *
     * \throws gmx::FileIOError if the value frame dimensions are not 1-dimensional
     * or do not exactly match the frame dimensions of other already added data sets.
     */
    H5mdParticleBlockBuilder& setPosition(H5mdTimeDataBlock<RVec>&& position);

    /*! \brief Assign a velocity data set.
     *
     * \throws gmx::FileIOError if the value frame dimensions are not 1-dimensional
     * or do not exactly match the frame dimensions of other already added data sets.
     */
    H5mdParticleBlockBuilder& setVelocity(H5mdTimeDataBlock<RVec>&& velocity);

    /*! \brief Assign a force data set.
     *
     * \throws gmx::FileIOError if the value frame dimensions are not 1-dimensional
     * or do not exactly match the frame dimensions of other already added data sets.
     */
    H5mdParticleBlockBuilder& setForce(H5mdTimeDataBlock<RVec>&& force);

    //! \brief Build and return the \c H5mdParticleBlock.
    H5mdParticleBlock build();

private:
    /*! \brief Assert that an \c H5mdTimeDataBlock set \p dataBlock matches others.
     *
     * \note On first invocation this updates numParticles_ to the number of particles
     * in the added \p dataBlock. On subsequent invocations the number of particles
     * in \p dataBlock is compared against this value. Any inconsistency results
     * in a throw.
     *
     * \throws gmx::FileIOError if an inconsistency is found.
     */
    void ensureParticleCountConsistency(const H5mdTimeDataBlock<RVec>& dataBlock);

    //! \brief Positions of particles in block (if written to disk)
    std::optional<H5mdTimeDataBlock<RVec>> position_;

    //! \brief Velocities of particles in block (if written to disk)
    std::optional<H5mdTimeDataBlock<RVec>> velocity_;

    //! \brief Forces of particles in block (if written to disk)
    std::optional<H5mdTimeDataBlock<RVec>> force_;

    //! \brief Simulation box (if non-constant over the simulation)
    std::optional<H5mdFrameDataSet<real>> box_;

    //! \brief Number of particles for managed data sets.
    std::optional<int64_t> numParticles_ = std::nullopt;

    GMX_DISALLOW_COPY_MOVE_AND_ASSIGN(H5mdParticleBlockBuilder);
};

} // namespace gmx

#endif // GMX_FILEIO_H5MD_PARTICLEBLOCK_H
