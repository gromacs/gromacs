/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2004, The GROMACS development team.
 * Copyright (c) 2013,2014,2015,2016,2017 by the GROMACS development team.
 * Copyright (c) 2018,2019,2020, by the GROMACS development team, led by
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
#ifndef GMX_MDLIB_FORCE_H
#define GMX_MDLIB_FORCE_H

#include "gromacs/math/arrayrefwithpadding.h"
#include "gromacs/math/vectypes.h"
#include "gromacs/utility/arrayref.h"

class DDBalanceRegionHandler;
struct gmx_edsam;
struct gmx_enerdata_t;
struct gmx_enfrot;
struct SimulationGroups;
struct gmx_localtop_t;
struct gmx_multisim_t;
struct gmx_vsite_t;
struct gmx_wallcycle;
class history_t;
class InteractionDefinitions;
struct pull_t;
struct t_commrec;
struct t_fcdata;
struct t_forcerec;
struct t_inputrec;
struct t_lambda;
struct t_mdatoms;
struct t_nrnb;

namespace gmx
{
class Awh;
class ForceOutputs;
class ForceWithVirial;
class ImdSession;
class MdrunScheduleWorkload;
class MDLogger;
class StepWorkload;
} // namespace gmx

void do_force(FILE*                               log,
              const t_commrec*                    cr,
              const gmx_multisim_t*               ms,
              const t_inputrec*                   inputrec,
              gmx::Awh*                           awh,
              gmx_enfrot*                         enforcedRotation,
              gmx::ImdSession*                    imdSession,
              pull_t*                             pull_work,
              int64_t                             step,
              t_nrnb*                             nrnb,
              gmx_wallcycle*                      wcycle,
              const gmx_localtop_t*               top,
              const matrix                        box,
              gmx::ArrayRefWithPadding<gmx::RVec> coordinates,
              history_t*                          hist,
              gmx::ArrayRefWithPadding<gmx::RVec> force,
              tensor                              vir_force,
              const t_mdatoms*                    mdatoms,
              gmx_enerdata_t*                     enerd,
              t_fcdata*                           fcd,
              gmx::ArrayRef<real>                 lambda,
              t_forcerec*                         fr,
              gmx::MdrunScheduleWorkload*         runScheduleWork,
              const gmx_vsite_t*                  vsite,
              rvec                                mu_tot,
              double                              t,
              gmx_edsam*                          ed,
              int                                 legacyFlags,
              const DDBalanceRegionHandler&       ddBalanceRegionHandler);

/* Communicate coordinates (if parallel).
 * Do neighbor searching (if necessary).
 * Calculate forces.
 * Communicate forces (if parallel).
 * Spread forces for vsites (if present).
 *
 * f is always required.
 */


/* Compute listed forces, Ewald, PME corrections add when (when used).
 *
 * xWholeMolecules only needs to contain whole molecules when orientation
 * restraints need to be computed and can be empty otherwise.
 */
void do_force_lowlevel(t_forcerec*                               fr,
                       const t_inputrec*                         ir,
                       const InteractionDefinitions&             idef,
                       const t_commrec*                          cr,
                       const gmx_multisim_t*                     ms,
                       t_nrnb*                                   nrnb,
                       gmx_wallcycle*                            wcycle,
                       const t_mdatoms*                          md,
                       gmx::ArrayRefWithPadding<const gmx::RVec> coordinates,
                       gmx::ArrayRef<const gmx::RVec>            xWholeMolecules,
                       history_t*                                hist,
                       gmx::ForceOutputs*                        forceOutputs,
                       gmx_enerdata_t*                           enerd,
                       t_fcdata*                                 fcd,
                       const matrix                              box,
                       const real*                               lambda,
                       const rvec*                               mu_tot,
                       const gmx::StepWorkload&                  stepWork,
                       const DDBalanceRegionHandler&             ddBalanceRegionHandler);
/* Call all the force routines */

#endif
