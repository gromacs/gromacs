/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2004, The GROMACS development team.
 * Copyright (c) 2013,2014,2015,2016,2017, by the GROMACS development team, led by
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
#ifndef GMX_MDLIB_UPDATE_H
#define GMX_MDLIB_UPDATE_H

#include "gromacs/math/paddedvector.h"
#include "gromacs/math/vectypes.h"
#include "gromacs/timing/wallcycle.h"
#include "gromacs/utility/arrayref.h"
#include "gromacs/utility/basedefinitions.h"
#include "gromacs/utility/real.h"

class ekinstate_t;
struct gmx_constr;
struct gmx_ekindata_t;
struct gmx_enerdata_t;
struct t_extmass;
struct t_fcdata;
struct t_graph;
struct t_grpopts;
struct t_idef;
struct t_inputrec;
struct t_mdatoms;
struct t_nrnb;
class t_state;

/* Abstract type for update */
struct gmx_update_t;

/* Initialize the stochastic dynamics struct */
gmx_update_t *init_update(const t_inputrec *ir);

/* Update pre-computed constants that depend on the reference
 * temperature for coupling.
 *
 * This could change e.g. in simulated annealing. */
void update_temperature_constants(gmx_update_t *upd, const t_inputrec *ir);

/* Update the size of per-atom arrays (e.g. after DD re-partitioning,
   which might increase the number of home atoms). */
void update_realloc(gmx_update_t *upd, int natoms);

/* Store the box at step step
 * as a reference state for simulations with box deformation.
 */
void set_deform_reference_box(gmx_update_t *upd,
                              gmx_int64_t step, matrix box);

void update_tcouple(gmx_int64_t       step,
                    t_inputrec       *inputrec,
                    t_state          *state,
                    gmx_ekindata_t   *ekind,
                    t_extmass        *MassQ,
                    t_mdatoms        *md
                    );

/* Update Parrinello-Rahman, to be called before the coordinate update */
void update_pcouple_before_coordinates(FILE             *fplog,
                                       gmx_int64_t       step,
                                       const t_inputrec *inputrec,
                                       t_state          *state,
                                       matrix            parrinellorahmanMu,
                                       matrix            M,
                                       gmx_bool          bInitStep);

/* Update the box, to be called after the coordinate update.
 * For Berendsen P-coupling, also calculates the scaling factor
 * and scales the coordinates.
 * When the deform option is used, scales coordinates and box here.
 */
void update_pcouple_after_coordinates(FILE             *fplog,
                                      gmx_int64_t       step,
                                      const t_inputrec *inputrec,
                                      const t_mdatoms  *md,
                                      const matrix      pressure,
                                      const matrix      forceVirial,
                                      const matrix      constraintVirial,
                                      const matrix      parrinellorahmanMu,
                                      t_state          *state,
                                      t_nrnb           *nrnb,
                                      gmx_update_t     *upd);

void update_coords(FILE                          *fplog,
                   gmx_int64_t                    step,
                   t_inputrec                    *inputrec, /* input record and box stuff	*/
                   t_mdatoms                     *md,
                   t_state                       *state,
                   gmx::PaddedArrayRef<gmx::RVec> f, /* forces on home particles */
                   t_fcdata                      *fcd,
                   gmx_ekindata_t                *ekind,
                   matrix                         M,
                   gmx_update_t                  *upd,
                   int                            bUpdatePart,
                   t_commrec                     *cr, /* these shouldn't be here -- need to think about it */
                   gmx_constr                    *constr);

/* Return TRUE if OK, FALSE in case of Shake Error */

extern gmx_bool update_randomize_velocities(t_inputrec *ir, gmx_int64_t step, const t_commrec *cr, t_mdatoms *md, t_state *state, gmx_update_t *upd, gmx_constr *constr);

void update_constraints(FILE                    *fplog,
                        gmx_int64_t              step,
                        real                    *dvdlambda, /* FEP stuff */
                        t_inputrec              *inputrec,  /* input record and box stuff	*/
                        t_mdatoms               *md,
                        t_state                 *state,
                        gmx_bool                 bMolPBC,
                        t_graph                 *graph,
                        gmx::ArrayRef<gmx::RVec> force, /* forces on home particles */
                        t_idef                  *idef,
                        tensor                   vir_part,
                        t_commrec               *cr,
                        t_nrnb                  *nrnb,
                        gmx_wallcycle_t          wcycle,
                        gmx_update_t            *upd,
                        gmx_constr              *constr,
                        gmx_bool                 bFirstHalf,
                        gmx_bool                 bCalcVir);

/* Return TRUE if OK, FALSE in case of Shake Error */

void calc_ke_part(t_state *state, t_grpopts *opts, t_mdatoms *md,
                  gmx_ekindata_t *ekind, t_nrnb *nrnb, gmx_bool bEkinAveVel);
/*
 * Compute the partial kinetic energy for home particles;
 * will be accumulated in the calling routine.
 * The tensor is
 *
 * Ekin = SUM(i) 0.5 m[i] v[i] (x) v[i]
 *
 *     use v[i] = v[i] - u[i] when calculating temperature
 *
 * u must be accumulated already.
 *
 * Now also computes the contribution of the kinetic energy to the
 * free energy
 *
 */


void
init_ekinstate(ekinstate_t *ekinstate, const t_inputrec *ir);

void
update_ekinstate(ekinstate_t *ekinstate, gmx_ekindata_t *ekind);

/*! \brief Restores data from \p ekinstate to \p ekind, then broadcasts it
   to the rest of the simulation */
void
restore_ekinstate_from_state(t_commrec *cr,
                             gmx_ekindata_t *ekind, const ekinstate_t *ekinstate);

void berendsen_tcoupl(const t_inputrec *ir, const gmx_ekindata_t *ekind, real dt,
                      std::vector<double> &therm_integral);

void andersen_tcoupl(t_inputrec *ir, gmx_int64_t step,
                     const t_commrec *cr, const t_mdatoms *md, t_state *state, real rate, const gmx_bool *randomize, const real *boltzfac);

void nosehoover_tcoupl(t_grpopts *opts, gmx_ekindata_t *ekind, real dt,
                       double xi[], double vxi[], t_extmass *MassQ);

void trotter_update(t_inputrec *ir, gmx_int64_t step, gmx_ekindata_t *ekind,
                    gmx_enerdata_t *enerd, t_state *state, tensor vir, t_mdatoms *md,
                    t_extmass *MassQ, int **trotter_seqlist, int trotter_seqno);

int **init_npt_vars(t_inputrec *ir, t_state *state, t_extmass *Mass, gmx_bool bTrotter);

real NPT_energy(const t_inputrec *ir, const t_state *state, const t_extmass *MassQ);
/* computes all the pressure/tempertature control energy terms to get a conserved energy */

void NBaroT_trotter(t_grpopts *opts, real dt,
                    double xi[], double vxi[], real *veta, t_extmass *MassQ);

void vrescale_tcoupl(t_inputrec *ir, gmx_int64_t step,
                     gmx_ekindata_t *ekind, real dt,
                     double therm_integral[]);
/* Compute temperature scaling. For V-rescale it is done in update. */

void rescale_velocities(gmx_ekindata_t *ekind, t_mdatoms *mdatoms,
                        int start, int end, rvec v[]);
/* Rescale the velocities with the scaling factor in ekind */

void update_annealing_target_temp(t_inputrec *ir, real t, gmx_update_t *upd);
/* Set reference temp for simulated annealing at time t*/

real calc_temp(real ekin, real nrdf);
/* Calculate the temperature */

real calc_pres(int ePBC, int nwall, matrix box, tensor ekin, tensor vir,
               tensor pres);
/* Calculate the pressure tensor, returns the scalar pressure.
 * The unit of pressure is bar.
 */

void parrinellorahman_pcoupl(FILE *fplog, gmx_int64_t step,
                             const t_inputrec *ir, real dt, const tensor pres,
                             tensor box, tensor box_rel, tensor boxv,
                             tensor M, matrix mu,
                             gmx_bool bFirstStep);

void berendsen_pcoupl(FILE *fplog, gmx_int64_t step,
                      const t_inputrec *ir, real dt,
                      const tensor pres, const matrix box,
                      const matrix force_vir, const matrix constraint_vir,
                      matrix mu, double *baros_integral);

void berendsen_pscale(const t_inputrec *ir, const matrix mu,
                      matrix box, matrix box_rel,
                      int start, int nr_atoms,
                      rvec x[], const unsigned short cFREEZE[],
                      t_nrnb *nrnb);

void correct_ekin(FILE *log, int start, int end, rvec v[],
                  rvec vcm, real mass[], real tmass, tensor ekin);
/* Correct ekin for vcm */

#endif
