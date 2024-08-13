/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2012- The GROMACS Authors
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
/*! \libinternal \file
 *
 * \brief This file contains function declarations necessary for
 * managing automatic load balance of PME calculations (Coulomb and
 * LJ).
 *
 * \author Berk Hess <hess@kth.se>
 * \inlibraryapi
 * \ingroup module_ewald
 */

#ifndef GMX_EWALD_PME_LOAD_BALANCING_H
#define GMX_EWALD_PME_LOAD_BALANCING_H

#include "gromacs/math/vectypes.h"
#include "gromacs/timing/wallcycle.h"


struct t_commrec;
struct t_forcerec;
struct t_inputrec;
struct interaction_const_t;
struct gmx_pme_t;
class t_state;

namespace gmx
{
struct nonbonded_verlet_t;
class MDLogger;
template<typename T>
class ArrayRef;
} // namespace gmx

/*! \brief Object to manage PME load balancing */
struct pme_load_balancing_t;

/*! \brief Return whether PME load balancing is active */
bool pme_loadbal_is_active(const pme_load_balancing_t* pme_lb);

/*! \brief Initialize the PP-PME load balacing data and infrastructure
 *
 * Initialize the PP-PME load balacing data and infrastructure.
 * The actual load balancing might start right away, later or never.
 * The PME grid in pmedata is reused for smaller grids to lower the memory
 * usage.
 */
void pme_loadbal_init(pme_load_balancing_t**         pme_lb_p,
                      t_commrec*                     cr,
                      const gmx::MDLogger&           mdlog,
                      const t_inputrec&              ir,
                      const matrix                   box,
                      const interaction_const_t&     ic,
                      const gmx::nonbonded_verlet_t& nbv,
                      gmx_pme_t*                     pmedata,
                      gmx_bool                       bUseGPU);

/*! \brief Process cycles and PME load balance when necessary
 *
 * Process the cycles measured over the last nstlist steps and then
 * either continue balancing or check if we need to trigger balancing.
 * Should be called after the WallCycleCounter::Step cycle counter has been stopped.
 * Returns if the load balancing is printing to fp_err.
 */
void pme_loadbal_do(pme_load_balancing_t*          pme_lb,
                    struct t_commrec*              cr,
                    FILE*                          fp_err,
                    FILE*                          fp_log,
                    const gmx::MDLogger&           mdlog,
                    const t_inputrec&              ir,
                    t_forcerec*                    fr,
                    const matrix                   box,
                    gmx::ArrayRef<const gmx::RVec> x,
                    gmx_wallcycle*                 wcycle,
                    int64_t                        step,
                    int64_t                        step_rel,
                    gmx_bool*                      bPrinting,
                    bool                           useGpuPmePpCommunication);

/*! \brief Finish the PME load balancing and print the settings when fplog!=NULL */
void pme_loadbal_done(pme_load_balancing_t* pme_lb, FILE* fplog, const gmx::MDLogger& mdlog, gmx_bool bNonBondedOnGPU);

#endif
