/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2018,2019, by the GROMACS development team, led by
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
 * Declares the AcceptOrRewind class.
 *
 * This class provides the Metropolis step with the current energies and, in return,
 * gets a boolean on whether to accept or reject the current configuration. It
 * updates t_state accordingly (and signals to restart domain decomposition if it is used).
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
#include "gromacs/utility/futil.h"
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
class IMetropolisStep;
class MDAtoms;

class AcceptOrRewind
{
    public:
        /*! \brief Constructor */
        AcceptOrRewind(IMetropolisStep &metropolisStep);

        /*! \brief To make sure globalState_ is initialised properly */
        void initialiseGlobalState(t_state *globalState, t_commrec *cr);

        /*! \brief Updates the class t_state according to the result of the Metropolis step.
         *         If the configuration is accepted, it makes a back-up; if the configuration is rejected,
         *         it rewinds to the last accepted configuration (local state without domain decomposition,
         *         global state with domain decomposition). It signals to start over new with domain
         *         decomposition in case of rejection in order to obtain a valid local state.
         */
        bool run(const int64_t         step,
                 gmx_enerdata_t       *enerd,
                 t_state              *localState,
                 t_state              *globalState,
                 const t_commrec      *cr);

        /*! \brief Prints the result of the Metropolis step(s) to the log file */
        void printOutput(FILE *fplog, const int nstlog, const int64_t step);

        /* gets and sets for acceptedProposedConfigurations_ and totalProposedConfigurations_*/
        int getAcceptedProposedConfigurations();
        int getTotalProposedConfigurations();
        void setAcceptedProposedConfigurations(int inputValue);
        void setTotalProposedConfigurations(int inputValue);

    private:
        IMetropolisStep &metropolisStep_;
        t_state          localStateBackUp_;
        t_state          globalStateBackUp_;
        bool             accepted_;
        int              totalProposedConfigurations_;
        int              acceptedProposedConfigurations_;
};

} // namespace

#endif
