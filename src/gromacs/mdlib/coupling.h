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
#ifndef GMX_MDLIB_COUPLING_H
#define GMX_MDLIB_COUPLING_H

#include <cstdio>

#include <array>
#include <vector>

#include "gromacs/math/matrix.h"
#include "gromacs/math/vectypes.h"
#include "gromacs/mdtypes/md_enums.h"
#include "gromacs/utility/enumerationhelpers.h"
#include "gromacs/utility/real.h"

class gmx_ekindata_t;
struct t_commrec;
struct t_extmass;
struct t_grpopts;
struct t_inputrec;
struct t_nrnb;
class t_state;

enum class PbcType;
struct PressureCouplingOptions;

namespace gmx
{
class BoxDeformation;
class MDLogger;
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

/* Update Parrinello-Rahman, to be called before the coordinate update
 * Returns the box-scaling matrix mu and the coordinate-scaling matrix M */
void update_pcouple_before_coordinates(const gmx::MDLogger&           mdlog,
                                       int64_t                        step,
                                       const PressureCouplingOptions& pressureCoupling,
                                       const tensor                   deform,
                                       real                           delta_t,
                                       t_state*                       state,
                                       gmx::Matrix3x3*                parrinellorahmanMu,
                                       gmx::Matrix3x3*                M);

/*! \brief Implement aspects of pressure coupling.
 *
 * Always called after the coordinate update due to forces. This
 * ensures that when the coordinates should be updated by pressure
 * coupling, the box and the shift vectors remained consistent at the
 * time of the coordinate update.
 *
 * For Berendsen and c-rescale, on pressure-coupling steps:
 * - recalculate the scaling matrix mu,
 * - use mu to scale coordinates iff \p scaleCoordinates,
 *   doing positions always and velocites iff c-rescale
 * - use mu to update box size
 *
 * For Parrinello-Rahman, on the step one before pressure-coupling
 * steps:
 * - update the box according to the new velocity computed in
 *   parrinellorahman_pcoupl
 * - use mu to scale positions iff \p scaleCoordinates
 *
 * For MTTK, on pressure-coupling steps:
 * - update the box according to state->veta
 * - recalculate the box velocity according to state->veta and the new box
 * - the positions are NOT scaled, as they are already scaled during the integration
 *
 * When box deformation is active, always applies it.
 *
 * /param[out] pressureCouplingMu  Used later by the GPU update to scale the
 *                                 coordinates to suit the new box size
 */
void update_pcouple_after_coordinates(FILE*                               fplog,
                                      int64_t                             step,
                                      const PressureCouplingOptions&      pressureCoupling,
                                      int64_t                             ld_seed,
                                      real                                ensembleTemperature,
                                      const ivec*                         nFreeze,
                                      const tensor                        deform,
                                      real                                delta_t,
                                      int                                 homenr,
                                      gmx::ArrayRef<const unsigned short> cFREEZE,
                                      const matrix                        pressure,
                                      const matrix                        forceVirial,
                                      const matrix                        constraintVirial,
                                      gmx::Matrix3x3*                     pressureCouplingMu,
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

void trotter_update(const t_inputrec*                   ir,
                    int64_t                             step,
                    gmx_ekindata_t*                     ekind,
                    t_state*                            state,
                    const tensor                        vir,
                    int                                 homenr,
                    gmx::ArrayRef<const unsigned short> cTC,
                    gmx::ArrayRef<const real>           invMass,
                    const t_extmass*                    MassQ,
                    gmx::ArrayRef<std::vector<int>>     trotter_seqlist,
                    TrotterSequence                     trotter_seqno);

void init_npt_masses(const t_inputrec& ir, const gmx_ekindata_t& ekind, t_state* state, t_extmass* MassQ, bool bInit);

gmx::EnumerationArray<TrotterSequence, std::vector<int>> init_npt_vars(const t_inputrec*     ir,
                                                                       const gmx_ekindata_t& ekind,
                                                                       t_state*              state,
                                                                       t_extmass*            Mass,
                                                                       bool bTrotter);

real NPT_energy(const PressureCouplingOptions& pressureCoupling,
                TemperatureCoupling            etc,
                gmx::ArrayRef<const real>      degreesOfFreedom,
                const gmx_ekindata_t&          ekind,
                bool                           isTrotterWithConstantTemperature,
                const t_state*                 state,
                const t_extmass*               MassQ);
/* computes all the pressure/tempertature control energy terms to get a conserved energy */

void vrescale_tcoupl(const t_inputrec*     ir,
                     int64_t               step,
                     gmx_ekindata_t*       ekind,
                     real                  delta_t,
                     gmx::ArrayRef<double> therm_integral);
/* Compute temperature scaling. For V-rescale it is done in update. */

void rescale_velocities(const gmx_ekindata_t*               ekind,
                        gmx::ArrayRef<const unsigned short> cTC,
                        int                                 start,
                        int                                 end,
                        gmx::ArrayRef<gmx::RVec>            v);
/* Rescale the velocities with the scaling factor in ekind */

//! Initialize simulated annealing.
bool initSimulatedAnnealing(const t_inputrec& ir, gmx_ekindata_t* ekind, gmx::Update* upd);

//! Set reference temperature for simulated annealing
void update_annealing_target_temp(const t_inputrec& ir, real t, gmx_ekindata_t* ekind, gmx::Update* upd);

real calc_temp(real ekin, real nrdf);
/* Calculate the temperature */

real calc_pres(PbcType pbcType, int nwall, const matrix box, const tensor ekin, const tensor vir, tensor pres);
/* Calculate the pressure tensor, returns the scalar pressure.
 * The unit of pressure is bar.
 */

/*! \brief Initialize the Parrinello-Rahman mu and M tensors
 *
 * mu describes the relative change in box vectors introduced by the
 * coupling and applied at each step. This ensures there is no large
 * jump in box-vector values.
 *
 * M is used in the update code to couple the change in particle
 * velocities to that of the box vectors introduced by mu.
 *
 * \param[in]  pressureCouplingOptions  The pressure-coupling options
 * \param[in]  deform                   The matrix describing any box deformation
 * \param[in]  couplingTimePeriod       The period between changes to the pressure coupling
 * \param[in]  box                      The simulation box vectors
 * \param[in]  box_rel                  The relative box vectors
 * \param[in]  boxv                     The simulation box vector velocity for Parrinello-Rahman coupling
 * \param[out] M                        The scaling factor for particle velocities
 * \param[out] mu                       The scaling factor for box vectors
 */
void init_parrinellorahman(const PressureCouplingOptions& pressureCouplingOptions,
                           const tensor                   deform,
                           real                           couplingTimePeriod,
                           const tensor                   box,
                           tensor                         box_rel,
                           tensor                         boxv,
                           gmx::Matrix3x3*                M,
                           gmx::Matrix3x3*                mu);

/*! \brief Calculate the change in box vectors due to Parrinello-Rahman pressure coupling
 *
 * \param[in]  mdlog                    Log file for warning when pressure change is large
 * \param[in]  step                     The MD step count for the above warning
 * \param[in]  pressureCouplingOptions  The pressure-coupling options
 * \param[in]  deform                   The matrix describing any box deformation
 * \param[in]  couplingTimePeriod       The period between changes to the pressure coupling
 * \param[in]  pres                     The current pressure
 * \param[in]  box                      The simulation box vectors
 * \param[in]  box_rel                  The relative box vectors
 * \param[in]  boxv                     The simulation box vector velocity for Parrinello-Rahman coupling
 * \param[out] M                        The scaling factor for particle velocities
 * \param[out] mu                       The scaling factor for box vectors
 *
 * This doesn't do any coordinate updating. It just integrates the box
 * vector equations from the calculated acceleration due to pressure
 * difference. We also compute the tensor M which is used in update to
 * couple the particle coordinates to the box vectors.
 *
 * We also do NOT update the box vectors themselves here, since we
 * need them for shifting later. It is instead done last in the update
 * routines.
 *
 * In Nose and Klein (Mol.Phys 50 (1983) no 5., p 1055) this is given
 * as
 *
 *            -1    .           .     -1
 * M_nk = (h')   * (h' * h + h' h) * h
 *
 * with the dots denoting time derivatives, apostrophes denoting
 * transposition, and h is the transformation from the scaled frame to
 * the real frame, i.e. the TRANSPOSE of the box. This also goes for
 * the pressure and M tensors - they are transposed relative to
 * ours. Our equation thus becomes:
 *
 *                  -1       .    .           -1
 * M_gmx = M_nk' = b  * (b * b' + b * b') * b'
 *
 * where b is the gromacs box matrix. Our box accelerations are given
 * by
 *
 *   ..                                    ..
 *   b = vol/W inv(box') * (P-ref_P)     (=h')
 *
 */
void parrinellorahman_pcoupl(const gmx::MDLogger&           mdlog,
                             int64_t                        step,
                             const PressureCouplingOptions& pressureCouplingOptions,
                             const tensor                   deform,
                             real                           couplingTimePeriod,
                             const tensor                   pres,
                             const tensor                   box,
                             tensor                         box_rel,
                             tensor                         boxv,
                             gmx::Matrix3x3*                M,
                             gmx::Matrix3x3*                mu);

/*! \brief Calculate the pressure coupling scaling matrix
 *
 * Used by Berendsen and C-Rescale pressure coupling, this function
 * computes the current value of the scaling matrix. The template
 * parameter determines the pressure coupling algorithm.
 */
template<PressureCoupling pressureCouplingType>
void pressureCouplingCalculateScalingMatrix(FILE*                          fplog,
                                            int64_t                        step,
                                            const PressureCouplingOptions& pressureCoupling,
                                            int64_t                        ld_seed,
                                            real                           ensembleTemperature,
                                            real                           delta_t,
                                            const tensor                   pres,
                                            const matrix                   box,
                                            const matrix                   force_vir,
                                            const matrix                   constraint_vir,
                                            gmx::Matrix3x3*                mu,
                                            double*                        baros_integral);

/*! \brief Scale the box and coordinates
 *
 * Used by Berendsen and C-Rescale pressure coupling, this function scales
 * the box, the positions, and the velocities (C-Rescale only) according to
 * the scaling matrix mu. The template parameter determines the pressure
 * coupling algorithm.
 */
template<PressureCoupling pressureCouplingType>
void pressureCouplingScaleBoxAndCoordinates(const PressureCouplingOptions&      pressureCoupling,
                                            const tensor                        deform,
                                            const ivec*                         nFreeze,
                                            const gmx::Matrix3x3&               mu,
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
