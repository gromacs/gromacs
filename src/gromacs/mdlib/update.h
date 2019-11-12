/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2004, The GROMACS development team.
 * Copyright (c) 2013,2014,2015,2016,2017,2018,2019, by the GROMACS development team, led by
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
#include "gromacs/mdtypes/md_enums.h"
#include "gromacs/timing/wallcycle.h"
#include "gromacs/utility/arrayref.h"
#include "gromacs/utility/basedefinitions.h"
#include "gromacs/utility/classhelpers.h"
#include "gromacs/utility/real.h"

class ekinstate_t;
struct gmx_ekindata_t;
struct gmx_enerdata_t;
struct t_extmass;
struct t_fcdata;
struct t_graph;
struct t_grpopts;
struct t_inputrec;
struct t_mdatoms;
struct t_nrnb;
class t_state;

/* Abstract type for update */
struct gmx_stochd_t;

namespace gmx
{
class BoxDeformation;
class Constraints;


/*! \libinternal
 * \brief Contains data for update phase */
class Update
{
public:
    //! Constructor
    Update(const t_inputrec* ir, BoxDeformation* boxDeformation);
    ~Update();
    // TODO Get rid of getters when more free functions are incorporated as member methods
    //! Returns handle to stochd_t struct
    gmx_stochd_t* sd() const;
    //! Returns pointer to PaddedVector xp
    PaddedVector<gmx::RVec>* xp();
    //! Returns handle to box deformation class
    BoxDeformation* deform() const;
    //! Resizes xp
    void setNumAtoms(int nAtoms);

private:
    //! Implementation type.
    class Impl;
    //! Implementation object.
    PrivateImplPointer<Impl> impl_;
};

}; // namespace gmx

/* Update pre-computed constants that depend on the reference
 * temperature for coupling.
 *
 * This could change e.g. in simulated annealing. */
void update_temperature_constants(gmx_stochd_t* sd, const t_inputrec* ir);

/* Update the size of per-atom arrays (e.g. after DD re-partitioning,
   which might increase the number of home atoms). */

void update_tcouple(int64_t           step,
                    const t_inputrec* inputrec,
                    t_state*          state,
                    gmx_ekindata_t*   ekind,
                    const t_extmass*  MassQ,
                    const t_mdatoms*  md);

/* Update Parrinello-Rahman, to be called before the coordinate update */
void update_pcouple_before_coordinates(FILE*             fplog,
                                       int64_t           step,
                                       const t_inputrec* inputrec,
                                       t_state*          state,
                                       matrix            parrinellorahmanMu,
                                       matrix            M,
                                       gmx_bool          bInitStep);

/* Update the box, to be called after the coordinate update.
 * For Berendsen P-coupling, also calculates the scaling factor
 * and scales the coordinates.
 * When the deform option is used, scales coordinates and box here.
 */
void update_pcouple_after_coordinates(FILE*             fplog,
                                      int64_t           step,
                                      const t_inputrec* inputrec,
                                      const t_mdatoms*  md,
                                      const matrix      pressure,
                                      const matrix      forceVirial,
                                      const matrix      constraintVirial,
                                      matrix            pressureCouplingMu,
                                      t_state*          state,
                                      t_nrnb*           nrnb,
                                      gmx::Update*      upd,
                                      bool              scaleCoordinates);

void update_coords(int64_t           step,
                   const t_inputrec* inputrec, /* input record and box stuff	*/
                   const t_mdatoms*  md,
                   t_state*          state,
                   gmx::ArrayRefWithPadding<const gmx::RVec> f, /* forces on home particles */
                   const t_fcdata*                           fcd,
                   const gmx_ekindata_t*                     ekind,
                   const matrix                              M,
                   gmx::Update*                              upd,
                   int                                       bUpdatePart,
                   const t_commrec* cr, /* these shouldn't be here -- need to think about it */
                   const gmx::Constraints* constr);

/* Return TRUE if OK, FALSE in case of Shake Error */

extern gmx_bool update_randomize_velocities(const t_inputrec*        ir,
                                            int64_t                  step,
                                            const t_commrec*         cr,
                                            const t_mdatoms*         md,
                                            gmx::ArrayRef<gmx::RVec> v,
                                            const gmx::Update*       upd,
                                            const gmx::Constraints*  constr);

void constrain_velocities(int64_t step,
                          real* dvdlambda, /* the contribution to be added to the bonded interactions */
                          t_state*          state,
                          tensor            vir_part,
                          gmx::Constraints* constr,
                          gmx_bool          bCalcVir,
                          bool              do_log,
                          bool              do_ene);

void constrain_coordinates(int64_t step,
                           real* dvdlambda, /* the contribution to be added to the bonded interactions */
                           t_state*          state,
                           tensor            vir_part,
                           gmx::Update*      upd,
                           gmx::Constraints* constr,
                           gmx_bool          bCalcVir,
                           bool              do_log,
                           bool              do_ene);

void update_sd_second_half(int64_t step,
                           real* dvdlambda, /* the contribution to be added to the bonded interactions */
                           const t_inputrec* inputrec, /* input record and box stuff */
                           const t_mdatoms*  md,
                           t_state*          state,
                           const t_commrec*  cr,
                           t_nrnb*           nrnb,
                           gmx_wallcycle_t   wcycle,
                           gmx::Update*      upd,
                           gmx::Constraints* constr,
                           bool              do_log,
                           bool              do_ene);

void finish_update(const t_inputrec*       inputrec,
                   const t_mdatoms*        md,
                   t_state*                state,
                   const t_graph*          graph,
                   t_nrnb*                 nrnb,
                   gmx_wallcycle_t         wcycle,
                   gmx::Update*            upd,
                   const gmx::Constraints* constr);

/* Return TRUE if OK, FALSE in case of Shake Error */

void calc_ke_part(const rvec*      x,
                  const rvec*      v,
                  const matrix     box,
                  const t_grpopts* opts,
                  const t_mdatoms* md,
                  gmx_ekindata_t*  ekind,
                  t_nrnb*          nrnb,
                  gmx_bool         bEkinAveVel);
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


void init_ekinstate(ekinstate_t* ekinstate, const t_inputrec* ir);

void update_ekinstate(ekinstate_t* ekinstate, const gmx_ekindata_t* ekind);

/*! \brief Restores data from \p ekinstate to \p ekind, then broadcasts it
   to the rest of the simulation */
void restore_ekinstate_from_state(const t_commrec* cr, gmx_ekindata_t* ekind, const ekinstate_t* ekinstate);

void berendsen_tcoupl(const t_inputrec*    ir,
                      gmx_ekindata_t*      ekind,
                      real                 dt,
                      std::vector<double>& therm_integral); //NOLINT(google-runtime-references)

void andersen_tcoupl(const t_inputrec*         ir,
                     int64_t                   step,
                     const t_commrec*          cr,
                     const t_mdatoms*          md,
                     gmx::ArrayRef<gmx::RVec>  v,
                     real                      rate,
                     const std::vector<bool>&  randomize,
                     gmx::ArrayRef<const real> boltzfac);

void nosehoover_tcoupl(const t_grpopts*      opts,
                       const gmx_ekindata_t* ekind,
                       real                  dt,
                       double                xi[],
                       double                vxi[],
                       const t_extmass*      MassQ);

void trotter_update(const t_inputrec*               ir,
                    int64_t                         step,
                    gmx_ekindata_t*                 ekind,
                    const gmx_enerdata_t*           enerd,
                    t_state*                        state,
                    const tensor                    vir,
                    const t_mdatoms*                md,
                    const t_extmass*                MassQ,
                    gmx::ArrayRef<std::vector<int>> trotter_seqlist,
                    int                             trotter_seqno);

std::array<std::vector<int>, ettTSEQMAX>
init_npt_vars(const t_inputrec* ir, t_state* state, t_extmass* Mass, gmx_bool bTrotter);

real NPT_energy(const t_inputrec* ir, const t_state* state, const t_extmass* MassQ);
/* computes all the pressure/tempertature control energy terms to get a conserved energy */

void vrescale_tcoupl(const t_inputrec* ir, int64_t step, gmx_ekindata_t* ekind, real dt, double therm_integral[]);
/* Compute temperature scaling. For V-rescale it is done in update. */

void rescale_velocities(const gmx_ekindata_t* ekind, const t_mdatoms* mdatoms, int start, int end, rvec v[]);
/* Rescale the velocities with the scaling factor in ekind */

//! Check whether we do simulated annealing.
bool doSimulatedAnnealing(const t_inputrec* ir);

//! Initialize simulated annealing.
bool initSimulatedAnnealing(t_inputrec* ir, gmx::Update* upd);

// TODO: This is the only function in update.h altering the inputrec
void update_annealing_target_temp(t_inputrec* ir, real t, gmx::Update* upd);
/* Set reference temp for simulated annealing at time t*/

real calc_temp(real ekin, real nrdf);
/* Calculate the temperature */

real calc_pres(int ePBC, int nwall, const matrix box, const tensor ekin, const tensor vir, tensor pres);
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
                             gmx_bool          bFirstStep);

void berendsen_pcoupl(FILE*             fplog,
                      int64_t           step,
                      const t_inputrec* ir,
                      real              dt,
                      const tensor      pres,
                      const matrix      box,
                      const matrix      force_vir,
                      const matrix      constraint_vir,
                      matrix            mu,
                      double*           baros_integral);

void berendsen_pscale(const t_inputrec*    ir,
                      const matrix         mu,
                      matrix               box,
                      matrix               box_rel,
                      int                  start,
                      int                  nr_atoms,
                      rvec                 x[],
                      const unsigned short cFREEZE[],
                      t_nrnb*              nrnb,
                      bool                 scaleCoordinates);

void pleaseCiteCouplingAlgorithms(FILE* fplog, const t_inputrec& ir);

/*! \brief Computes the atom range for a thread to operate on, ensuring SIMD aligned ranges
 *
 * \param[in]  numThreads   The number of threads to divide atoms over
 * \param[in]  threadIndex  The thread to get the range for
 * \param[in]  numAtoms     The total number of atoms (on this rank)
 * \param[out] startAtom    The start of the atom range
 * \param[out] endAtom      The end of the atom range, note that this is in general not a multiple of the SIMD width
 */
void getThreadAtomRange(int numThreads, int threadIndex, int numAtoms, int* startAtom, int* endAtom);

/*! \brief Generate a new kinetic energy for the v-rescale thermostat
 *
 * Generates a new value for the kinetic energy, according to
 * Bussi et al JCP (2007), Eq. (A7)
 *
 * This is used by update_tcoupl(), and by the VRescaleThermostat of the modular
 * simulator.
 * TODO: Move this to the VRescaleThermostat once the modular simulator becomes
 *       the default code path.
 *
 * @param kk     present value of the kinetic energy of the atoms to be thermalized (in arbitrary units)
 * @param sigma  target average value of the kinetic energy (ndeg k_b T/2)  (in the same units as kk)
 * @param ndeg   number of degrees of freedom of the atoms to be thermalized
 * @param taut   relaxation time of the thermostat, in units of 'how often this routine is called'
 * @param step   the time step this routine is called on
 * @param seed   the random number generator seed
 * @return  the new kinetic energy
 */
real vrescale_resamplekin(real kk, real sigma, real ndeg, real taut, int64_t step, int64_t seed);

#endif
