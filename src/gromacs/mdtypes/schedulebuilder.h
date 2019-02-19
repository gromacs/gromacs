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

#ifndef GMX_MDTYPES_SCHEDULEBUILDER_H
#define GMX_MDTYPES_SCHEDULEBUILDER_H

#include "defaultschedule.h"

namespace gmx
{

/*! \libinternal \brief
 * Class to encapsulate the initialization of concrete schedules
 *
 * The class has the responsibility of initializing schedules either
 * by using a selector function to process the execution context or
 * manually specifying an identifier for the concrete schedules.
 *
 * \inlibraryapi
 * \ingroup module_mdtypes
 */

class ScheduleBuilder
{

    private:
        std::shared_ptr<IForceSchedule> root_;
        ScheduleType                    type_                  = None;
        bool                            experimentalSchedules_ = false;

        void setSchedule()
        {
            if (type_ == Verlet)
            {
                root_ = std::make_shared<DefaultVerletSchedule>();
            }
            else if (type_ == Group)
            {
                root_ = std::make_shared<DefaultGroupSchedule>();
            }
        }

    public:

        ScheduleBuilder() = default;

        /*! Invokes the use of experimental schedules under development.
         *  Uses the default VERLET scheme when FALSE.
         */

        explicit ScheduleBuilder(bool experimentalSchedules)
        {
            experimentalSchedules_ = experimentalSchedules;
        }

        //! Selects a schedule automatically based on inputs and execution context
        void selectSchedule(t_inputrec *inputrec
                            /* execution context and problem definition */)
        {
            switch (inputrec->cutoff_scheme)
            {
                case ecutsVERLET:
                    type_ = Verlet;
                    break;
                case ecutsGROUP:
                    type_ = Group;
                    break;
                default:
                    gmx_incons("Invalid cut-off scheme passed!");
            }

            setSchedule();

        }

        //! Allows to change the schedule type directly
        void setScheduleType(ScheduleType type)
        {
            type_ = type;
            setSchedule();
        }

        //! Wrapper to initialize a schedule with the current context
        void initSchedule(FILE *fplog, t_commrec *cr, const gmx_multisim_t *ms, t_inputrec *inputrec, t_nrnb *nrnb,
                          gmx_wallcycle *wcycle, gmx_localtop_t *top, gmx_groups_t *groups,
                          t_state *state, PaddedVector<gmx::RVec>* force, tensor *vir_force,
                          t_mdatoms *mdatoms, gmx_enerdata_t *enerd, t_fcdata *fcd, t_graph *graph, t_forcerec *fr,
                          gmx_vsite_t *vsite, rvec *mu_tot, gmx_edsam *ed, MDLoopSharedPrimitives *sp, Awh *awh,
                          gmx_enfrot *enforcedRotation, PpForceWorkload *ppForceWorkload)
        {
            if (!root_)
            {
                selectSchedule(inputrec);
            }
            root_->init(fplog, cr, ms, inputrec, nrnb, wcycle, top, groups, state, force,
                        vir_force, mdatoms, enerd, fcd, graph, fr, vsite, mu_tot, ed, sp, awh,
                        enforcedRotation, ppForceWorkload);
        }

        //! Returns a pointer to a schedule
        std::shared_ptr<IForceSchedule> getSchedule()
        {
            if (root_)
            {
                return root_;
            }
            else
            {
                gmx_incons("Schedule uninitialized!");
            }
        }

};

}      // namespace gmx

#endif //GMX_MDTYPES_SCHEDULEBUILDER_H
