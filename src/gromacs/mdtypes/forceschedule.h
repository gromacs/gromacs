/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2018, by the GROMACS development team, led by
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
/*! \internal \file
 *
 * \brief Describes the Schedule abstraction for Force Calculations
 * in each MD step
 *
 * \author Prashanth Kanduri <kanduri@cscs.ch>
 * \ingroup module_mdlib, gmxlib
 */

#ifndef GROMACS_FORCESCHEDULE_H
#define GROMACS_FORCESCHEDULE_H

#include "gromacs/mdlib/integrator.h"

enum class DdOpenBalanceRegionBeforeForceComputation;
enum class DdCloseBalanceRegionAfterForceComputation;

class t_state;
typedef struct gmx_edsam *gmx_edsam_t;
struct gmx_groups_t;

namespace gmx
{


struct MDLoopSharedPrimitives
{
    gmx_int64_t                               *step_;
    double                                    *t_;     //time

    gmx_bool                                  *bBornRadii_;

    DdOpenBalanceRegionBeforeForceComputation *ddOpenBalanceRegion_;
    DdCloseBalanceRegionAfterForceComputation *ddCloseBalanceRegion_;

    MDLoopSharedPrimitives(gmx_int64_t *step, double *t, gmx_bool *bBornRadii,
                           DdOpenBalanceRegionBeforeForceComputation *ddOpenBalanceRegion,
                           DdCloseBalanceRegionAfterForceComputation *ddCloseBalanceRegion) :

        step_(step), t_(t), bBornRadii_(bBornRadii), ddOpenBalanceRegion_(ddOpenBalanceRegion),
        ddCloseBalanceRegion_(ddCloseBalanceRegion)
    {
        // Initialize the struct with the values that can change during the simulation
    }
};

class AbstractNBSchedule
{
    protected:
        // Members
        FILE                           *log_;
        t_commrec                      *cr_;
        gmx_multisim_t                 *ms_;
        t_inputrec                     *inputrec_;

        struct t_nrnb                  *nrnb_;
        gmx_wallcycle_t                 wcycle_;
        gmx_localtop_t                 *top_;
        gmx_groups_t                   *groups_;
        t_state                        *state_;

        gmx::PaddedArrayRef<gmx::RVec> *force_;
        tensor                         *vir_force_;
        t_mdatoms                      *mdatoms_;
        gmx_enerdata_t                 *enerd_;
        t_fcdata                       *fcd_;

        t_graph                        *graph_;
        t_forcerec                     *fr_;
        gmx_vsite_t                    *vsite_;
        rvec                           *mu_tot_;

        struct gmx_edsam               *ed_;
        struct MDLoopSharedPrimitives  *sp_;

    public:

        AbstractNBSchedule();                                             // Null initialization

        ~AbstractNBSchedule();                                            // Default destructor

        AbstractNBSchedule(const AbstractNBSchedule &abstractNBSchedule); //copy constructor

        void init( FILE *fplog, t_commrec *cr, gmx_multisim_t *ms, t_inputrec *inputrec, t_nrnb *nrnb,
                   gmx_wallcycle_t wcycle, gmx_localtop_t *top, gmx_groups_t *groups,
                   t_state *state, gmx::PaddedArrayRef<gmx::RVec> *force, tensor *vir_force,
                   t_mdatoms *mdatoms, gmx_enerdata_t *enerd, t_fcdata *fcd, t_graph *graph, t_forcerec *fr,
                   gmx_vsite_t *vsite, rvec *mu_tot, gmx_edsam_t ed, MDLoopSharedPrimitives *sp);

        virtual void computeStep(int flags) = 0;
};



class DefaultNBSchedule : public AbstractNBSchedule
{

    public:

        void computeStep(int flags);

};

// TODO: Implement a schedule selector function. Also find a way to wrap the execution context
//       and problem details into nice, readable objects
// static AbstractNBSchedule selectSchedule(/* execution context and problem definition */);

}      // namespace gmx

#endif //GROMACS_FORCESCHEDULE_H
