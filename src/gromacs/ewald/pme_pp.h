/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2020, by the GROMACS development team, led by
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
 * \brief This file contains function declarations necessary for
 * mananging the PP side of PME-only ranks.
 *
 * \author Berk Hess <hess@kth.se>
 * \author Mark Abraham <mark.j.abraham@gmail.com>
 * \ingroup module_ewald
 */

#ifndef GMX_EWALD_PME_PP_H
#define GMX_EWALD_PME_PP_H

#include "gromacs/math/vectypes.h"
#include "gromacs/utility/basedefinitions.h"

struct gmx_wallcycle;
struct interaction_const_t;
struct t_commrec;
struct t_forcerec;

class GpuEventSynchronizer;

namespace gmx
{
class ForceWithVirial;
class PmePpCommGpu;
} // namespace gmx

/*! \brief Send the charges and maxshift to out PME-only node. */
void gmx_pme_send_parameters(const t_commrec*           cr,
                             const interaction_const_t* ic,
                             gmx_bool                   bFreeEnergy_q,
                             gmx_bool                   bFreeEnergy_lj,
                             real*                      chargeA,
                             real*                      chargeB,
                             real*                      sqrt_c6A,
                             real*                      sqrt_c6B,
                             real*                      sigmaA,
                             real*                      sigmaB,
                             int                        maxshift_x,
                             int                        maxshift_y);

/*! \brief Send the coordinates to our PME-only node and request a PME calculation */
void gmx_pme_send_coordinates(t_forcerec*           fr,
                              const t_commrec*      cr,
                              const matrix          box,
                              const rvec*           x,
                              real                  lambda_q,
                              real                  lambda_lj,
                              gmx_bool              bEnerVir,
                              int64_t               step,
                              bool                  useGpuPmePpComms,
                              bool                  reinitGpuPmePpComms,
                              bool                  sendCoordinatesFromGpu,
                              GpuEventSynchronizer* coordinatesReadyOnDeviceEvent,
                              gmx_wallcycle*        wcycle);

/*! \brief Tell our PME-only node to finish */
void gmx_pme_send_finish(const t_commrec* cr);

/*! \brief Tell our PME-only node to reset all cycle and flop counters */
void gmx_pme_send_resetcounters(const t_commrec* cr, int64_t step);

/*! \brief PP nodes receive the long range forces from the PME nodes */
void gmx_pme_receive_f(gmx::PmePpCommGpu*    pmePpCommGpu,
                       const t_commrec*      cr,
                       gmx::ForceWithVirial* forceWithVirial,
                       real*                 energy_q,
                       real*                 energy_lj,
                       real*                 dvdlambda_q,
                       real*                 dvdlambda_lj,
                       bool                  useGpuPmePpComms,
                       bool                  receivePmeForceToGpu,
                       float*                pme_cycles);

/*! \brief Tell our PME-only node to switch to a new grid size */
void gmx_pme_send_switchgrid(const t_commrec* cr, ivec grid_size, real ewaldcoeff_q, real ewaldcoeff_lj);

#endif
