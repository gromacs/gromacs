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

/*! \brief Definitions of H5md particle data block routines.
 *
 * \author Petter Johansson <pettjoha@kth.se>
 */

#include "h5md_particleblock.h"

#include "gromacs/fileio/h5md/exceptions.h"
#include "gromacs/utility/stringutil.h"

namespace gmx
{

void H5mdParticleBlock::addPosition(H5mdTimeDataBlock<RVec>&& position)
{
    GMX_H5MD_THROW_UPON_ERROR(position_.has_value(), "Cannot add position data block: already set");
    GMX_H5MD_THROW_UPON_ERROR(
            position.frameDims().size() != 1,
            formatString("Position data block must have a 1-dimensional frame (had %ld dimensions)",
                         gmx::ssize(position.frameDims())));

    const int64_t numParticlesInBlock = static_cast<int64_t>(position.frameDims()[0]);
    GMX_H5MD_THROW_UPON_ERROR(
            numParticles_ != numParticlesInBlock,
            formatString("Inconsistent number of particles: particle block is for %" PRId64
                         " particles, but trying to add position data block for %" PRId64,
                         numParticles_,
                         numParticlesInBlock));

    position_.emplace(std::move(position));
}

void H5mdParticleBlock::addVelocity(H5mdTimeDataBlock<RVec>&& velocity)
{
    GMX_H5MD_THROW_UPON_ERROR(velocity_.has_value(), "Cannot add velocity data block: already set");
    GMX_H5MD_THROW_UPON_ERROR(
            velocity.frameDims().size() != 1,
            formatString("Velocity data block must have a 1-dimensional frame (had %ld dimensions)",
                         gmx::ssize(velocity.frameDims())));

    const int64_t numParticlesInBlock = static_cast<int64_t>(velocity.frameDims()[0]);
    GMX_H5MD_THROW_UPON_ERROR(
            numParticles_ != numParticlesInBlock,
            formatString("Inconsistent number of particles: particle block is for %" PRId64
                         " particles, but trying to add velocity data block for %" PRId64,
                         numParticles_,
                         numParticlesInBlock));

    velocity_.emplace(std::move(velocity));
}

void H5mdParticleBlock::addForce(H5mdTimeDataBlock<RVec>&& force)
{
    GMX_H5MD_THROW_UPON_ERROR(force_.has_value(), "Cannot add force data block: already set");
    GMX_H5MD_THROW_UPON_ERROR(
            force.frameDims().size() != 1,
            formatString("Force data block must have a 1-dimensional frame (had %ld dimensions)",
                         gmx::ssize(force.frameDims())));

    const int64_t numParticlesInBlock = static_cast<int64_t>(force.frameDims()[0]);
    GMX_H5MD_THROW_UPON_ERROR(
            numParticles_ != numParticlesInBlock,
            formatString("Inconsistent number of particles: particle block is for %" PRId64
                         " particles, but trying to add force data block for %" PRId64,
                         numParticles_,
                         numParticlesInBlock));

    force_.emplace(std::move(force));
}

void H5mdParticleBlock::addBox(H5mdFrameDataSet<real>&& box)
{
    GMX_H5MD_THROW_UPON_ERROR(box_.has_value(), "Cannot add box data block: already set");
    GMX_H5MD_THROW_UPON_ERROR(
            box.frameDims().size() != 2 || box.frameDims()[0] != DIM || box.frameDims()[1] != DIM,
            "Simulation box does not have required dimensions { 3, 3 }");
    box_.emplace(std::move(box));
}

void H5mdParticleBlock::trimDataSetsFromStep(const int64_t maxStep)
{
    if (hasPosition())
    {
        position_->trimFramesFromMaxStep(maxStep);
    }
    if (hasVelocity())
    {
        velocity_->trimFramesFromMaxStep(maxStep);
    }
    if (hasForce())
    {
        force_->trimFramesFromMaxStep(maxStep);
    }
    if (hasBox())
    {
        box_->shrinkToNumFrames(position_->numFrames());
    }
}

} // namespace gmx
