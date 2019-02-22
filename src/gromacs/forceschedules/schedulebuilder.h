/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2019, by the GROMACS development team, led by
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

/*! \libinternal \file
 *
 * \brief Builds objects that schedules the operations necessary
 * for force calculations in particular contexts.
 *
 * \author Prashanth Kanduri <kanduri@cscs.ch>
 *
 * \ingroup module_forceschedules
 * \inlibraryapi
 */

#ifndef GMX_FORCESCHEDULES_SCHEDULEBUILDER_H
#define GMX_FORCESCHEDULES_SCHEDULEBUILDER_H

#include <cstdio>

#include <memory>

#include "gromacs/math/paddedvector.h"
#include "gromacs/math/vectypes.h"

class DDBalanceRegionHandler;
class t_state;
struct t_commrec;
struct t_fcdata;
struct t_forcerec;
struct t_graph;
struct t_inputrec;
struct t_mdatoms;
struct t_nrnb;
struct gmx_enerdata_t;
struct gmx_edsam;
struct gmx_enfrot;
struct gmx_groups_t;
struct gmx_localtop_t;
struct gmx_multisim_t;
struct gmx_wallcycle;
struct gmx_vsite_t;

namespace gmx
{

class Awh;
class IForceSchedule;
class PpForceWorkload;
struct MDLoopSharedPrimitives;

/*! \libinternal \brief
 * Class to encapsulate the initialization of concrete schedules
 *
 * The class has the responsibility of initializing schedules either
 * by using a selector function to process the execution context or
 * manually specifying an identifier for the concrete schedules.
 *
 * \todo Later this builder will get setter functions to become aware
 * of hardware, parallelism and method considerations, and then
 * generateSchedule would contain the logic that reacts to those by
 * building an object that implements IForceSchedule that e.g.
 * do_md can use transparently.
 */
class ScheduleBuilder
{
    public:
        /*! \brief Creates, initializes and returns a pointer to a schedule
         *
         * \todo Consider taking a templated parameter pack (like
         * makeConstraints()) so that the arguments for the
         * constructor of IForceSchedule do not have to be duplicated
         * here, but rather passed through with std::forward. */
        std::unique_ptr<IForceSchedule> generateSchedule(FILE                         *fplog,
                                                         t_commrec                    *cr,
                                                         const gmx_multisim_t         *ms,
                                                         t_inputrec                   *inputrec,
                                                         t_nrnb                       *nrnb,
                                                         gmx_wallcycle                *wcycle,
                                                         gmx_localtop_t               *top,
                                                         gmx_groups_t                 *groups,
                                                         t_state                      *state,
                                                         PaddedVector<gmx::RVec>     * force,
                                                         tensor                       *vir_force,
                                                         t_mdatoms                    *mdatoms,
                                                         gmx_enerdata_t               *enerd,
                                                         t_fcdata                     *fcd,
                                                         t_graph                      *graph,
                                                         t_forcerec                   *fr,
                                                         gmx_vsite_t                  *vsite,
                                                         rvec                         *mu_tot,
                                                         gmx_edsam                    *ed,
                                                         MDLoopSharedPrimitives       *sp,
                                                         Awh                          *awh,
                                                         gmx_enfrot                   *enforcedRotation,
                                                         PpForceWorkload              *ppForceWorkload,
                                                         const DDBalanceRegionHandler *ddBalanceRegionHandler);
};

}      // namespace gmx

#endif //GMX_FORCESCHEDULES_SCHEDULEBUILDER_H
