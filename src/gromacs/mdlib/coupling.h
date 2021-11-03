/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2004, The GROMACS development team.
 * Copyright (c) 2013,2014,2015,2016,2017 by the GROMACS development team.
 * Copyright (c) 2018,2019,2020,2021, by the GROMACS development team, led by
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
#ifndef GMX_MDLIB_COUPLING_H
#define GMX_MDLIB_COUPLING_H

#include <cstdio>

#include <array>
#include <vector>

#include "gromacs/math/vectypes.h"
#include "gromacs/mdtypes/md_enums.h"
#include "gromacs/utility/enumerationhelpers.h"
#include "gromacs/utility/real.h"

class gmx_ekindata_t;
struct gmx_enerdata_t;
struct t_commrec;
struct t_extmass;
struct t_grpopts;
struct t_inputrec;
struct t_nrnb;
class t_state;

enum class PbcType;

namespace gmx
{
class BoxDeformation;
class Constraints;
class Update;
template<typename>
class ArrayRef;
}; // namespace gmx

/* Update the size of per-atom arrays (e.g. after DD re-partitioning,
   which might increase the number of home atoms). */

void update_tcouple(int64_t                             step,
                    const t_inputrec*                   inputrec,
                    t_state*                            state,
                    gmx_ekindata_t*                     ekind,
                    const t_extmass*                    MassQ,
                    int                                 homenr,
                    gmx::ArrayRef<const unsigned short> cTC);

/* Update Parrinello-Rahman, to be called before the coordinate update */
void update_pcouple_before_coordinates(FILE*             fplog,
                                       int64_t           step,
                                       const t_inputrec* inputrec,
                                       t_state*          state,
                                       matrix            parrinellorahmanMu,
                                       matrix            M,
                                       bool              bInitStep);

/* Update the box, to be called after the coordinate update.
 * For Berendsen P-coupling, also calculates the scaling factor
 * and scales the coordinates.
 * When the deform option is used, scales coordinates and box here.
 */
void update_pcouple_after_coordinates(FILE*                               fplog,
                                      int64_t                             step,
                                      const t_inputrec*                   inputrec,
                                      int                                 homenr,
                                      gmx::ArrayRef<const unsigned short> cFREEZE,
                                      const matrix                        pressure,
                                      const matrix                        forceVirial,
                                      const matrix                        constraintVirial,
                                      matrix                              pressureCouplingMu,
                                      t_state*                            state,
                                      t_nrnb*                             nrnb,
                                      gmx::BoxDeformation*                boxDeformation,
                                      bool                                scaleCoordinates);

/* Return TRUE if OK, FALSE in case of Shake Error */

extern bool update_randomize_velocities(const t_inputrec*                   ir,
                                        int64_t                             step,
                                        const t_commrec*                    cr,
                                        int                                 homenr,
                                        gmx::ArrayRef<const unsigned short> cTC,
                                        gmx::ArrayRef<const real>           invMass,
                                        gmx::ArrayRef<gmx::RVec>            v,
                                        const gmx::Update*                  upd,
                                        const gmx::Constraints*             constr);

void berendsen_tcoupl(const t_inputrec*    ir,
                      gmx_ekindata_t*      ekind,
                      real                 dt,
                      std::vector<double>& therm_integral); //NOLINT(google-runtime-references)

void andersen_tcoupl(const t_inputrec*                   ir,
                     int64_t                             step,
                     const t_commrec*                    cr,
                     int                                 homenr,
                     gmx::ArrayRef<const unsigned short> cTC,
                     gmx::ArrayRef<const real>           invMass,
                     gmx::ArrayRef<gmx::RVec>            v,
                     real                                rate,
                     const std::vector<bool>&            randomize,
                     gmx::ArrayRef<const real>           boltzfac);

void nosehoover_tcoupl(const t_grpopts*      opts,
                       const gmx_ekindata_t* ekind,
                       real                  dt,
                       gmx::ArrayRef<double> xi,
                       gmx::ArrayRef<double> vxi,
                       const t_extmass*      MassQ);

void trotter_update(const t_inputrec*                   ir,
                    int64_t                             step,
                    gmx_ekindata_t*                     ekind,
                    const gmx_enerdata_t*               enerd,
                    t_state*                            state,
                    const tensor                        vir,
                    int                                 homenr,
                    gmx::ArrayRef<const unsigned short> cTC,
                    gmx::ArrayRef<const real>           invMass,
                    const t_extmass*                    MassQ,
                    gmx::ArrayRef<std::vector<int>>     trotter_seqlist,
                    TrotterSequence                     trotter_seqno);

gmx::EnumerationArray<TrotterSequence, std::vector<int>>
init_npt_vars(const t_inputrec* ir, t_state* state, t_extmass* Mass, bool bTrotter);

real NPT_energy(const t_inputrec* ir, const t_state* state, const t_extmass* MassQ);
/* computes all the pressure/tempertature control energy terms to get a conserved energy */

void vrescale_tcoupl(const t_inputrec*     ir,
                     int64_t               step,
                     gmx_ekindata_t*       ekind,
                     real                  dt,
                     gmx::ArrayRef<double> therm_integral);
/* Compute temperature scaling. For V-rescale it is done in update. */

void rescale_velocities(const gmx_ekindata_t*               ekind,
                        gmx::ArrayRef<const unsigned short> cTC,
                        int                                 start,
                        int                                 end,
                        gmx::ArrayRef<gmx::RVec>            v);
/* Rescale the velocities with the scaling factor in ekind */

/*!
 * \brief Compute the new annealing temperature for a temperature group
 *
 * \param inputrec          The input record
 * \param temperatureGroup  The temperature group
 * \param time              The current time
 * \return  The new reference temperature for the group
 */
real computeAnnealingTargetTemperature(const t_inputrec& inputrec, int temperatureGroup, real time);

//! Check whether we do simulated annealing.
bool doSimulatedAnnealing(const t_inputrec* ir);

//! Initialize simulated annealing.
bool initSimulatedAnnealing(t_inputrec* ir, gmx::Update* upd);

// TODO: This is the only function in update.h altering the inputrec
void update_annealing_target_temp(t_inputrec* ir, real t, gmx::Update* upd);
/* Set reference temp for simulated annealing at time t*/

real calc_temp(real ekin, real nrdf);
/* Calculate the temperature */

real calc_pres(PbcType pbcType, int nwall, const matrix box, const tensor ekin, const tensor vir, tensor pres);
/* Calculate the pressure tensor, returns the scalar pressure.
 * The unit of pressure is bar.
 */

void parrinellorahman_pcoupl(FILE*             fplog,
                             int64_t           step,
                             const t_inputrec* ir,
                             real              dt,
                             const tensor      pres,
                             const tensor      box,
                             tensor            box_rel,
                             tensor            boxv,
                             tensor            M,
                             matrix            mu,
                             bool              bFirstStep);

/*! \brief Calculate the pressure coupling scaling matrix
 *
 * Used by Berendsen and C-Rescale pressure coupling, this function
 * computes the current value of the scaling matrix. The template
 * parameter determines the pressure coupling algorithm.
 */
template<PressureCoupling pressureCouplingType>
void pressureCouplingCalculateScalingMatrix(FILE*             fplog,
                                            int64_t           step,
                                            const t_inputrec* ir,
                                            real              dt,
                                            const tensor      pres,
                                            const matrix      box,
                                            const matrix      force_vir,
                                            const matrix      constraint_vir,
                                            matrix            mu,
                                            double*           baros_integral);

/*! \brief Scale the box and coordinates
 *
 * Used by Berendsen and C-Rescale pressure coupling, this function scales
 * the box, the positions, and the velocities (C-Rescale only) according to
 * the scaling matrix mu. The template parameter determines the pressure
 * coupling algorithm.
 */
template<PressureCoupling pressureCouplingType>
void pressureCouplingScaleBoxAndCoordinates(const t_inputrec*                   ir,
                                            const matrix                        mu,
                                            matrix                              box,
                                            matrix                              box_rel,
                                            int                                 start,
                                            int                                 nr_atoms,
                                            gmx::ArrayRef<gmx::RVec>            x,
                                            gmx::ArrayRef<gmx::RVec>            v,
                                            gmx::ArrayRef<const unsigned short> cFREEZE,
                                            t_nrnb*                             nrnb,
                                            bool                                scaleCoordinates);

void pleaseCiteCouplingAlgorithms(FILE* fplog, const t_inputrec& ir);

/*! \brief Generate a new kinetic energy for the v-rescale thermostat
 *
 * Generates a new value for the kinetic energy, according to
 * Bussi et al JCP (2007), Eq. (A7)
 *
 * This is used by update_tcoupl(), and by the VRescaleThermostat of the modular
 * simulator.
 * \todo Move this to the VRescaleThermostat once the modular simulator becomes
 *       the default code path.
 *
 * \param[in] kk     present value of the kinetic energy of the atoms to be thermalized (in
 *                   arbitrary units)
 * \param[in] sigma  target average value of the kinetic energy (ndeg k_b T/2)  (in
 *                   the same units as kk)
 * \param[in] ndeg   number of degrees of freedom of the atoms to be thermalized
 * \param[in] taut   relaxation time of the thermostat, in units of 'how often this
 *                   routine is called'
 * \param[in] step   the time step this routine is called on
 * \param[in] seed   the random number generator seed
 * \return  the new kinetic energy
 */
real vrescale_resamplekin(real kk, real sigma, real ndeg, real taut, int64_t step, int64_t seed);

#endif // GMX_MDLIB_COUPLING_H
