/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2021, by the GROMACS development team, led by
 * Mark Abraham, David van der Spoel, Berk Hess, and Erik Lindahl,
 * and including many others, as listed in the AUTHORS file in the
 * top-level source directory and at http://www.gromacs.org.
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
 * http://www.gnu.org/licenses, or write to the Free Software Foundation,
 * Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA.
 *
 * If you want to redistribute modifications to GROMACS, please
 * consider that scientific software is very special. Version
 * control is crucial - bugs must be traceable. We will be happy to
 * consider code for inclusion in the official distribution, but
 * derived work must not be called official GROMACS. Details are found
 * in the README & COPYING files - if they are missing, get the
 * official version at http://www.gromacs.org.
 *
 * To help us fund GROMACS development, we humbly ask that you cite
 * the research papers on the package. Check out http://www.gromacs.org.
 */
/*! \inpublicapi \file
 * \brief
 * Implements nblib simulation box
 *
 * \author Victor Holanda <victor.holanda@cscs.ch>
 * \author Joe Jordan <ejjordan@kth.se>
 * \author Prashanth Kanduri <kanduri@cscs.ch>
 * \author Sebastian Keller <keller@cscs.ch>
 */
#ifndef NBLIB_SYSTEMDESCRIPTION_H
#define NBLIB_SYSTEMDESCRIPTION_H

#include "gromacs/ewald/ewald_utils.h"
#include "gromacs/gmxlib/nrnb.h"
#include "gromacs/gpu_utils/device_stream_manager.h"
#include "gromacs/mdlib/rf_util.h"
#include "gromacs/mdtypes/enerdata.h"
#include "gromacs/mdtypes/forcerec.h"
#include "gromacs/mdtypes/interaction_const.h"
#include "gromacs/mdtypes/simulation_workload.h"
#include "gromacs/nbnxm/gpu_data_mgmt.h"
#include "gromacs/nbnxm/nbnxm_gpu.h"
#include "gromacs/nbnxm/atomdata.h"
#include "gromacs/nbnxm/nbnxm.h"
#include "gromacs/nbnxm/pairlistset.h"
#include "gromacs/nbnxm/pairlistsets.h"
#include "gromacs/nbnxm/pairsearch.h"
#include "gromacs/timing/wallcycle.h"
#include "gromacs/utility/listoflists.h"
#include "gromacs/utility/range.h"

namespace nblib
{

//! \brief client-side provided system description data
struct SystemDescription
{
    SystemDescription() = default;

    SystemDescription(gmx::ArrayRef<int>     particleTypeIdOfAllParticles,
                      gmx::ArrayRef<real>    nonBondedParams,
                      gmx::ArrayRef<real>    charges,
                      gmx::ArrayRef<int64_t> particleInteractionFlags)
    {
        std::array inputSizes{ particleTypeIdOfAllParticles.size(),
                               charges.size(),
                               particleInteractionFlags.size() };
        if (static_cast<unsigned long>(std::count(begin(inputSizes), end(inputSizes), inputSizes[0]))
            != inputSizes.size())
        {
            throw InputException("input array size inconsistent");
        }

        int numParticleTypes = int(std::round(std::sqrt(nonBondedParams.size() / 2)));
        if (2 * numParticleTypes * numParticleTypes != int(nonBondedParams.size()))
        {
            throw InputException("Wrong size of nonBondedParams");
        }

        numParticles_     = particleTypeIdOfAllParticles.size();
        numParticleTypes_ = numParticleTypes;

        particleTypeIdOfAllParticles_ = std::vector<int>(particleTypeIdOfAllParticles.begin(),
                                                         particleTypeIdOfAllParticles.end());

        nonBondedParams_ = std::vector<real>(nonBondedParams.begin(), nonBondedParams.end());
        charges_         = std::vector<real>(charges.begin(), charges.end());
        particleInfo_ =
                std::vector<int64_t>(particleInteractionFlags.begin(), particleInteractionFlags.end());
    }

    //! number of particles
    size_t numParticles_{ 0 };

    //! number of particle types
    size_t numParticleTypes_{ 0 };

    //! particle type id of all particles
    std::vector<int> particleTypeIdOfAllParticles_;

    //! Storage for parameters for short range interactions.
    std::vector<real> nonBondedParams_;

    //! electrostatic charges
    std::vector<real> charges_;

    //! flag for each particle to set LJ and Q interactions
    std::vector<int64_t> particleInfo_;

    //! Legacy matrix for box
    Box box_{ 0 };
};


} // namespace nblib
#endif // NBLIB_SYSTEMDESCRIPTION_H
