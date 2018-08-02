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
 * \defgroup module_hybridMCMD Hybrid Monte Carlo integrator for md-vv
 * \ingroup group_mdrun
 * \brief
 * Implements the key part of a hybrid Monte Carlo integrator for md-vv.
 *
 * After drawing random initial velocities from a Boltzmann distribution,
 * a short MD simulation in the NVE ensemble is carried out. The resulting
 * configuration is accepted or rejected based on a Metropolis step.
 * In case of acceptance, the current simulation state (t_state) is saved;
 * in case of rejection, we rewind to the last back-up.
 *
 * \author Sebastian Wingbermuehle
 */

/*! \libinternal \file
 *
 * \brief
 * Declares the AcceptOrRewind class and an interface to a Metropolis step.
 *
 * This class provides the Metropolis step with the current energies and, in return,
 * gets a boolean on whether to accept or reject the current configuration. It
 * updates t_state accordingly (and restarts domain decomposition if it is used). 
 * TODO Extend this class to back-up the state of domain decomposition, dynamic
 *      loadbalancing, pme-tuning, and neighbour searching
 *
 * \author Sebastian Wingbermuehle
 * \inlibraryapi
 * \ingroup module_hybridMCMD
 */

#ifndef GMX_ACCEPT_OR_REWIND_H
#define GMX_ACCEPT_OR_REWIND_H

#include "gromacs/math/paddedvector.h"
#include "gromacs/mdtypes/state.h"
#include "gromacs/utility/arrayref.h"
#include "gromacs/utility/basedefinitions.h"
#include "gromacs/utility/real.h"

struct gmx_enerdata_t;
struct gmx_localtop_t;
struct gmx_mtop_t;
struct gmx_update_t;
struct gmx_vsite_t;
struct gmx_wallcycle;
struct t_commrec;
struct t_forcerec;
struct t_inputrec;
struct t_nrnb;

namespace gmx
{
class Constraints;
class MDAtoms;

class IMetropolisStep
{
    public:
        //! \brief Decide to accept or reject the new configuration
        virtual bool accept(const gmx_int64_t step, const gmx_enerdata_t *enerd) = 0;

        //! \brief Destructor
        virtual ~IMetropolisStep() {};
};

class AcceptOrRewind
{
    public:
        //! \brief Constructor
        AcceptOrRewind(IMetropolisStep *metropolisStep);
        bool updateTState(const gmx_int64_t     step,                                  // arguments for core functionality
                          t_state              *localState,
                          gmx_enerdata_t       *enerd,
                          FILE                 *fplog,                                 // arguments that are forwarded to domain decomposition functions
                          const t_commrec      *cr,
                          const gmx_mtop_t     *top_global,
                          const t_inputrec     *ir,
                          PaddedRVecVector     *f,
                          MDAtoms              *mdAtoms,
                          gmx_localtop_t       *top,
                          t_forcerec           *fr,
                          gmx_vsite_t          *vsite,
                          Constraints          *constr,
                          t_nrnb               *nrnb,
                          gmx_wallcycle        *wcycle,
                          gmx_update_t         *upd,
                          bool                 *shouldCheckNumberOfBondedInteractions); // bool required for signalling in md.cpp
        void printOutput(FILE *fplog, const int nstlog, const gmx_int64_t step);


    private:
        IMetropolisStep *metropolisStep_;
        t_state          localStateBackUp_;
        bool             accepted_;
        real             totalProposedConfigurations_;
        real             acceptedProposedConfigurations_;
};

} // namespace

#endif
