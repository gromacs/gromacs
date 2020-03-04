/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2004, The GROMACS development team.
 * Copyright (c) 2013,2014,2015,2016,2017 The GROMACS development team.
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
#ifndef GMX_MDLIB_MD_SUPPORT_H
#define GMX_MDLIB_MD_SUPPORT_H

#include "gromacs/mdlib/vcm.h"
#include "gromacs/timing/wallcycle.h"

struct gmx_ekindata_t;
struct gmx_enerdata_t;
struct gmx_global_stat;
struct gmx_multisim_t;
struct gmx_signalling_t;
struct t_extmass;
struct t_forcerec;
struct t_grpopts;
struct t_inputrec;
struct t_lambda;
struct t_nrnb;
class t_state;
struct t_trxframe;

namespace gmx
{
template<typename T>
class ArrayRef;
class Constraints;
class MDLogger;
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
/* After a new DD partitioning, we need to set a flag to schedule
 * global reduction of the total number of bonded interactions that
 * will be computed, to check none are missing. */
#define CGLO_CHECK_NUMBER_OF_BONDED_INTERACTIONS (1u << 12u)


/*! \brief Return the number of steps that will take place between
 * intra-simulation communications, given the constraints of the
 * inputrec. */
int computeGlobalCommunicationPeriod(const gmx::MDLogger& mdlog, t_inputrec* ir, const t_commrec* cr);

/*! \brief Return true if the \p value is equal across the set of multi-simulations
 *
 * \todo This duplicates some of check_multi_int. Consolidate. */
bool multisim_int_all_are_equal(const gmx_multisim_t* ms, int64_t value);

void rerun_parallel_comm(t_commrec* cr, t_trxframe* fr, gmx_bool* bLastStep);

//! \brief Allocate and initialize node-local state entries
void set_state_entries(t_state* state, const t_inputrec* ir, bool useModularSimulator);

/* Set the lambda values in the global state from a frame read with rerun */
void setCurrentLambdasRerun(int64_t           step,
                            const t_lambda*   fepvals,
                            const t_trxframe* rerun_fr,
                            const double*     lam0,
                            t_state*          globalState);

/* Set the lambda values at each step of mdrun when they change */
void setCurrentLambdasLocal(int64_t             step,
                            const t_lambda*     fepvals,
                            const double*       lam0,
                            gmx::ArrayRef<real> lambda,
                            int                 currentFEPState);

int multisim_min(const gmx_multisim_t* ms, int nmin, int n);
/* Set an appropriate value for n across the whole multi-simulation */


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
                     real                           vdwLambda,
                     const t_mdatoms*               mdatoms,
                     t_nrnb*                        nrnb,
                     t_vcm*                         vcm,
                     gmx_wallcycle_t                wcycle,
                     gmx_enerdata_t*                enerd,
                     tensor                         force_vir,
                     tensor                         shake_vir,
                     tensor                         total_vir,
                     tensor                         pres,
                     gmx::Constraints*              constr,
                     gmx::SimulationSignaller*      signalCoordinator,
                     const matrix                   lastbox,
                     int*                           totalNumberOfBondedInteractions,
                     gmx_bool*                      bSumEkinhOld,
                     int                            flags);

#endif
