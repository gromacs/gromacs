/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2004, The GROMACS development team.
 * Copyright (c) 2013,2014,2015,2016,2017,2018, by the GROMACS development team, led by
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

/*! \libinternal
 * \ingroup group_mdrun
 * \brief
 * Implements a part of a hybrid Monte Carlo integrator for md-vv.
 *
 * \author Sebastian Wingbermuehle
 */

/*! \libinternal \file
 *
 * \brief
 * Declares the class MetropolisStepMehlig.
 * TODO: Extend with other Metropolis criteria (e.g. for simulated tempering
 *       and expanded ensemble)
 *
 * This class implements the Metropolis criterion suggested by Mehlig et al. (1992)
 * for pure MD simulations. It calculates the value of the criterion, draws a random
 * number, and provides a boolean on whether to accept or reject the current
 * configuration.
 *
 * \author Sebastian Wingbermuehle
 * \inlibraryapi
 * \ingroup module_hybridMCMD
 */

#include "gromacs/hybridMCMD/acceptorrewind.h"
#include "gromacs/hybridMCMD/hybridmcmdvelocities.h"
#include "gromacs/random/threefry.h"
#include "gromacs/random/uniformrealdistribution.h"
#include "gromacs/utility/basedefinitions.h"
#include "gromacs/utility/real.h"

struct gmx_enerdata_t;

namespace gmx
{

class MetropolisStepMehlig : public IMetropolisStep, public IMetropolisStepVelocities
{
    public:
        //! \brief Constructor
        MetropolisStepMehlig(const real temperatureEnsemble, const real temperatureVelocities, const int seed, const int initStep, const bool startingFromCheckpoint);

        void setInitialKineticEnergy(const double initialKineticEnergy) override;
        bool accept(const int64_t step, const gmx_enerdata_t *enerd) override;

        //! \brief Destructor
        virtual ~MetropolisStepMehlig() override {}

    private:
        double calculateMetropolisCriterion(const double newPotentialEnergy, const double newKineticEnergy) const;
        double drawRandomNumber(const int64_t step);
        void updatePotentialEnergy(const double newPotentialEnergy);

        double     oldPotentialEnergy_;
        double     initialKineticEnergy_;
        const real temperatureEnsemble_;
        const real temperatureVelocities_;
        const int  initStep_;
        bool       startingFromCheckpoint_;
        ThreeFry2x64<0>                 rng_;
        UniformRealDistribution<real>   uniformRealDist_;
};

} // namespace
