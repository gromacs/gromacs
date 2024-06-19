/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 1991- The GROMACS Authors
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
#ifndef GMX_MDLIB_MD_SUPPORT_H
#define GMX_MDLIB_MD_SUPPORT_H

#include <cstdint>

#include "gromacs/math/vectypes.h"
#include "gromacs/mdlib/vcm.h"
#include "gromacs/mdtypes/md_enums.h"
#include "gromacs/timing/wallcycle.h"
#include "gromacs/utility/basedefinitions.h"

class gmx_ekindata_t;
struct gmx_enerdata_t;
struct gmx_global_stat;
struct t_extmass;
struct t_forcerec;
struct t_grpopts;
struct t_inputrec;
struct t_nrnb;
class t_state;
struct t_trxframe;
struct t_commrec;
struct t_mdatoms;

namespace gmx
{
template<typename T>
class ArrayRef;
class MDLogger;
class ObservablesReducer;
class SimulationSignaller;
} // namespace gmx

/* Define a number of flags to better control the information
 * passed to compute_globals in md.c and global_stat.
 */

/* we are computing the kinetic energy from average velocities */
#define CGLO_EKINAVEVEL (1u << 2u)
/* we are removing the center of mass momenta */
#define CGLO_STOPCM (1u << 3u)
/* bGStat is defined in do_md */
#define CGLO_GSTAT (1u << 4u)
/* Sum the energy terms in global computation */
#define CGLO_ENERGY (1u << 6u)
/* Sum the kinetic energy terms in global computation */
#define CGLO_TEMPERATURE (1u << 7u)
/* Sum the kinetic energy terms in global computation */
#define CGLO_PRESSURE (1u << 8u)
/* Sum the constraint term in global computation */
#define CGLO_CONSTRAINT (1u << 9u)
/* Reading ekin from the trajectory */
#define CGLO_READEKIN (1u << 10u)
/* we need to reset the ekin rescaling factor here */
#define CGLO_SCALEEKIN (1u << 11u)

/*! \brief Return the number of steps that will take place between
 * intra-simulation communications, given the constraints of the
 * inputrec. */
int computeGlobalCommunicationPeriod(const t_inputrec* ir);

/*! \brief Return the number of steps that will take place between
 * intra-simulation communications, given the constraints of the
 * inputrec, and write information to log.
 * Calls computeGlobalCommunicationPeriod(ir) internally. */
int computeGlobalCommunicationPeriod(const gmx::MDLogger& mdlog, const t_inputrec* ir, const t_commrec* cr);

void rerun_parallel_comm(t_commrec* cr, t_trxframe* fr, gmx_bool* bLastStep);

//! \brief Allocate and initialize node-local state entries
void set_state_entries(t_state* state, const t_inputrec* ir, bool useModularSimulator);

/* Compute global variables during integration
 *
 * Coordinates x are needed for kinetic energy calculation with cosine accelation
 * and for COM removal with rotational and acceleration correction modes.
 * Velocities v are needed for kinetic energy calculation and for COM removal.
 */
void compute_globals(gmx_global_stat*               gstat,
                     t_commrec*                     cr,
                     const t_inputrec*              ir,
                     t_forcerec*                    fr,
                     gmx_ekindata_t*                ekind,
                     gmx::ArrayRef<const gmx::RVec> x,
                     gmx::ArrayRef<const gmx::RVec> v,
                     const matrix                   box,
                     const t_mdatoms*               mdatoms,
                     t_nrnb*                        nrnb,
                     t_vcm*                         vcm,
                     gmx_wallcycle*                 wcycle,
                     gmx_enerdata_t*                enerd,
                     tensor                         force_vir,
                     tensor                         shake_vir,
                     tensor                         total_vir,
                     tensor                         pres,
                     gmx::SimulationSignaller*      signalCoordinator,
                     const matrix                   lastbox,
                     gmx_bool*                      bSumEkinhOld,
                     int                            flags,
                     int64_t                        step,
                     gmx::ObservablesReducer*       observablesReducer);

#endif
