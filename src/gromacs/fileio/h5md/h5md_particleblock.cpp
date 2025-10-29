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

#include "gromacs/utility/stringutil.h"

namespace gmx
{

H5mdParticleBlock::H5mdParticleBlock(std::optional<H5mdTimeDataBlock<RVec>>&& position,
                                     std::optional<H5mdTimeDataBlock<RVec>>&& velocity,
                                     std::optional<H5mdTimeDataBlock<RVec>>&& force,
                                     std::optional<H5mdFrameDataSet<real>>&&  box,
                                     int64_t                                  numParticles) :
    position_{ std::move(position) },
    velocity_{ std::move(velocity) },
    force_{ std::move(force) },
    box_{ std::move(box) },
    hasTime_{ (position_.has_value() && position_->hasTime())
              || (velocity_.has_value() && velocity_->hasTime())
              || (force_.has_value() && force_->hasTime()) },
    numParticles_{ numParticles }
{
}

H5mdParticleBlock H5mdParticleBlockBuilder::build()
{
    return H5mdParticleBlock(std::move(position_),
                             std::move(velocity_),
                             std::move(force_),
                             std::move(box_),
                             numParticles_.value_or(0));
}

void H5mdParticleBlockBuilder::ensureParticleCountConsistency(const H5mdTimeDataBlock<RVec>& dataBlock)
{
    throwUponH5mdError(dataBlock.frameDims().size() != 1,
                       "Trajectory data block does not have a 1d frame dimension");

    const int64_t numParticlesInBlock = static_cast<int64_t>(dataBlock.frameDims()[0]);
    if (!numParticles_.has_value())
    {
        numParticles_ = numParticlesInBlock;
    }
    else
    {
        throwUponH5mdError(numParticles_.value() != numParticlesInBlock,
                           formatString("Inconsistent number of particles: "
                                        "Trying to add data set for %lld particles, but already "
                                        "existing data sets are for %lld particles",
                                        static_cast<long long>(numParticlesInBlock),
                                        static_cast<long long>(numParticles_.value())));
    }
}

H5mdParticleBlockBuilder& H5mdParticleBlockBuilder::setBox(H5mdFrameDataSet<real>&& box)
{
    throwUponH5mdError(box.frameDims().size() != 2 || box.frameDims()[0] != DIM || box.frameDims()[1] != DIM,
                       "Simulation box does not have required dimensions { 3, 3 }");
    box_ = std::move(box);
    return *this;
}

H5mdParticleBlockBuilder& H5mdParticleBlockBuilder::setPosition(H5mdTimeDataBlock<RVec>&& position)
{
    ensureParticleCountConsistency(position);
    position_ = std::move(position);
    return *this;
}

H5mdParticleBlockBuilder& H5mdParticleBlockBuilder::setVelocity(H5mdTimeDataBlock<RVec>&& velocity)
{
    ensureParticleCountConsistency(velocity);
    velocity_ = std::move(velocity);
    return *this;
}

H5mdParticleBlockBuilder& H5mdParticleBlockBuilder::setForce(H5mdTimeDataBlock<RVec>&& force)
{
    ensureParticleCountConsistency(force);
    force_ = std::move(force);
    return *this;
}

} // namespace gmx
