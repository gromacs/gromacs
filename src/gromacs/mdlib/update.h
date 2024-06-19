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
#ifndef GMX_MDLIB_UPDATE_H
#define GMX_MDLIB_UPDATE_H

#include <cstdint>

#include <memory>
#include <vector>

#include "gromacs/math/matrix.h"
#include "gromacs/math/paddedvector.h"
#include "gromacs/math/vectypes.h"
#include "gromacs/mdtypes/md_enums.h"
#include "gromacs/timing/wallcycle.h"
#include "gromacs/utility/arrayref.h"
#include "gromacs/utility/real.h"

class ekinstate_t;
class gmx_ekindata_t;
struct gmx_enerdata_t;
enum class PbcType;
struct t_fcdata;
struct t_graph;
struct t_grpopts;
struct t_inputrec;
struct t_nrnb;
class t_state;
enum class ParticleType;
struct t_commrec;

namespace gmx
{
class BoxDeformation;
class Constraints;
template<typename T>
class ArrayRefWithPadding;


/*! \libinternal
 * \brief Contains data for update phase */
class Update
{
public:
    /*! \brief Constructor
     *
     * \param[in] inputRecord     Input record, used to construct SD object.
     * \param[in] ekind           Kinetic energy data
     * \param[in] boxDeformation  Periodic box deformation object.
     */
    Update(const t_inputrec& inputRecord, const gmx_ekindata_t& ekind, BoxDeformation* boxDeformation);
    //! Destructor
    ~Update();
    /*! \brief Get the pointer to updated coordinates
     *
     * Update saves the updated coordinates into separate buffer, so that constraints will have
     * access to both updated and not update coordinates. For that, update owns a separate buffer.
     * See finish_update(...) for details.
     *
     * \returns The pointer to the intermediate coordinates buffer.
     */
    PaddedVector<gmx::RVec>* xp();
    /*!\brief Getter to local copy of box deformation class.
     *
     * \returns handle to box deformation class
     */
    BoxDeformation* deform() const;
    /*! \brief Sets data that changes only at domain decomposition time.
     *
     * \param[in] numAtoms  Updated number of atoms.
     * \param[in] cFREEZE   Group index for freezing
     * \param[in] cTC       Group index for center of mass motion removal
     * \param[in] cAcceleration  Group index for constant acceleration groups
     */
    void updateAfterPartition(int                                 numAtoms,
                              gmx::ArrayRef<const unsigned short> cFREEZE,
                              gmx::ArrayRef<const unsigned short> cTC,
                              gmx::ArrayRef<const unsigned short> cAcceleration);

    /*! \brief Perform numerical integration step.
     *
     * Selects the appropriate integrator, based on the input record and performs a numerical integration step.
     *
     * \param[in]  inputRecord               Input record.
     * \param[in]  step                      Current timestep.
     * \param[in]  homenr                    The number of atoms on this processor.
     * \param[in]  havePartiallyFrozenAtoms  Whether atoms are frozen along 1 or 2 (not 3) dimensions?
     * \param[in]  ptype                     The list of particle types.
     * \param[in]  invMass                   Inverse atomic mass per atom, 0 for vsites and shells.
     * \param[in]  invMassPerDim             Inverse atomic mass per atom and dimension, 0 for vsites, shells and frozen dimensions
     * \param[in]  state                     System state object.
     * \param[in]  f                         Buffer with atomic forces for home particles.
     * \param[in]  fcdata                    Force calculation data to update distance and orientation restraints.
     * \param[in]  ekind                     Kinetic energy data (for temperature coupling, energy groups, etc.).
     * \param[in]  parrinelloRahmanM         Parrinello-Rahman velocity scaling matrix.
     * \param[in]  updatePart                What should be updated, coordinates or velocities. This enum only used in VV integrator.
     * \param[in]  cr                        Comunication record  (Old comment: these shouldn't be here -- need to think about it).
     * \param[in]  haveConstraints           If the system has constraints.
     */
    void update_coords(const t_inputrec&                                inputRecord,
                       int64_t                                          step,
                       int                                              homenr,
                       bool                                             havePartiallyFrozenAtoms,
                       gmx::ArrayRef<const ParticleType>                ptype,
                       gmx::ArrayRef<const real>                        invMass,
                       gmx::ArrayRef<const gmx::RVec>                   invMassPerDim,
                       t_state*                                         state,
                       const gmx::ArrayRefWithPadding<const gmx::RVec>& f,
                       t_fcdata*                                        fcdata,
                       const gmx_ekindata_t*                            ekind,
                       const Matrix3x3&                                 parrinelloRahmanM,
                       int                                              updatePart,
                       const t_commrec*                                 cr,
                       bool                                             haveConstraints);

    /*! \brief Finalize the coordinate update.
     *
     * Copy the updated coordinates to the main coordinates buffer for the atoms that are not frozen.
     *
     * \param[in]  inputRecord      Input record.
     * \param[in]  havePartiallyFrozenAtoms  Whether atoms are frozen along 1 or 2 (not 3) dimensions?
     * \param[in]  homenr                    The number of atoms on this processor.
     * \param[in]  state            System state object.
     * \param[in]  wcycle           Wall-clock cycle counter.
     * \param[in]  haveConstraints  If the system has constraints.
     */
    void finish_update(const t_inputrec& inputRecord,
                       bool              havePartiallyFrozenAtoms,
                       int               homenr,
                       t_state*          state,
                       gmx_wallcycle*    wcycle,
                       bool              haveConstraints);

    /*! \brief Secong part of the SD integrator.
     *
     * The first part of integration is performed in the update_coords(...) method.
     *
     * \param[in]  inputRecord  Input record.
     * \param[in]  step         Current timestep.
     * \param[in]  dvdlambda    Free energy derivative. Contribution to be added to
     *                          the bonded interactions.
     * \param[in]  homenr       The number of atoms on this processor.
     * \param[in]  ptype        The list of particle types.
     * \param[in]  invMass      Inverse atomic mass per atom, 0 for vsites and shells.
     * \param[in]  state        System state object.
     * \param[in]  cr           Comunication record.
     * \param[in]  nrnb         Cycle counters.
     * \param[in]  wcycle       Wall-clock cycle counter.
     * \param[in]  constr       Constraints object. The constraints are applied
     *                          on coordinates after update.
     * \param[in]  do_log       If this is logging step.
     * \param[in]  do_ene       If this is an energy evaluation step.
     */
    void update_sd_second_half(const t_inputrec&                 inputRecord,
                               int64_t                           step,
                               real*                             dvdlambda,
                               int                               homenr,
                               gmx::ArrayRef<const ParticleType> ptype,
                               gmx::ArrayRef<const real>         invMass,
                               t_state*                          state,
                               const t_commrec*                  cr,
                               t_nrnb*                           nrnb,
                               gmx_wallcycle*                    wcycle,
                               gmx::Constraints*                 constr,
                               bool                              do_log,
                               bool                              do_ene);

    /*! \brief Performs a leap-frog update without updating \p state so the constrain virial
     * can be computed.
     */
    void update_for_constraint_virial(const t_inputrec&              inputRecord,
                                      int                            homenr,
                                      bool                           havePartiallyFrozenAtoms,
                                      gmx::ArrayRef<const real>      invmass,
                                      gmx::ArrayRef<const gmx::RVec> invMassPerDim,
                                      const t_state&                 state,
                                      const gmx::ArrayRefWithPadding<const gmx::RVec>& f,
                                      const gmx_ekindata_t&                            ekind);

    /*! \brief Update pre-computed constants that depend on the reference temperature for coupling.
     *
     * This could change e.g. in simulated annealing.
     *
     * \param[in]  inputRecord  The input record
     * \param[in]  ekind        Kinetic energy data
     */
    void update_temperature_constants(const t_inputrec& inputRecord, const gmx_ekindata_t& ekind);

    /*!\brief Getter for the list of the randomize groups.
     *
     *  Needed for Andersen temperature control.
     *
     * \returns Reference to the groups from the SD data object.
     */
    const std::vector<bool>& getAndersenRandomizeGroup() const;
    /*!\brief Getter for the list of the Boltzmann factors.
     *
     *  Needed for Andersen temperature control.
     *
     * \returns Reference to the Boltzmann factors from the SD data object.
     */
    const std::vector<real>& getBoltzmanFactor() const;

private:
    //! Implementation type.
    class Impl;
    //! Implementation object.
    std::unique_ptr<Impl> impl_;
};

}; // namespace gmx

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

/*! \brief Updates \p ekinstate.
 *
 * Should be called on all ranks as a global reduction might be required.
 * This copies ekind->ekinh and ekind->dekindl, which are assumed to be computed
 * from the velocities at step at time t-dt/2.
 *
 * \param[out] ekinstate  The kinetic energy state to update
 * \param[in]  ekind      The kinetic energy data to store
 * \param[in]  sumEkin    Whether kinetic energy terms still need to be summed over all ranks
 * \param[in]  cr         Communication record, needed when sumEkin==true
 */
void update_ekinstate(ekinstate_t* ekinstate, const gmx_ekindata_t* ekind, bool sumEkin, const t_commrec* cr);

/*! \brief Restores data from \p ekinstate to \p ekind, then broadcasts it
   to the rest of the simulation */
void restore_ekinstate_from_state(const t_commrec* cr, gmx_ekindata_t* ekind, const ekinstate_t* ekinstate);

/*! \brief Computes the atom range for a thread to operate on, ensuring SIMD aligned ranges
 *
 * \param[in]  numThreads   The number of threads to divide atoms over
 * \param[in]  threadIndex  The thread to get the range for
 * \param[in]  numAtoms     The total number of atoms (on this rank)
 * \param[out] startAtom    The start of the atom range
 * \param[out] endAtom      The end of the atom range, note that this is in general not a multiple of the SIMD width
 */
void getThreadAtomRange(int numThreads, int threadIndex, int numAtoms, int* startAtom, int* endAtom);

#endif
