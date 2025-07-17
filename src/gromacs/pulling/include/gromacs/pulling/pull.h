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

/*! \libinternal \file
 *
 *
 * \brief
 * This file contains datatypes and function declarations necessary
   for mdrun to interface with the pull code.
 *
 * \author Berk Hess
 *
 * \inlibraryapi
 */

#ifndef GMX_PULLING_PULL_H
#define GMX_PULLING_PULL_H

#include <cstdio>

#include <optional>

#include "gromacs/mdtypes/pull_params.h"
#include "gromacs/utility/arrayref.h"
#include "gromacs/utility/basedefinitions.h"
#include "gromacs/utility/real.h"
#include "gromacs/utility/vectypes.h"

struct gmx_domdec_t;
struct gmx_mtop_t;
struct gmx_output_env_t;
struct pull_coord_work_t;
struct pull_params_t;
struct pull_t;
struct t_filenm;
struct t_inputrec;
struct t_pbc;
class t_state;
enum class PbcType;

namespace gmx
{
template<typename>
class ArrayRef;
class ForceWithVirial;
class MpiComm;
class LocalAtomSetManager;
} // namespace gmx

/*! \brief Returns the units of the pull coordinate.
 *
 * \param[in] pcrd The pull coordinate to query the units for.
 * \returns a string with the units of the coordinate.
 */
const char* pull_coordinate_units(const t_pull_coord& pcrd);

/*! \brief Returns the conversion factor from the pull coord init/rate unit to internal value unit.
 *
 * \param[in] pcrd The pull coordinate to get the conversion factor for.
 * \returns the conversion factor.
 */
double pull_conversion_factor_userinput2internal(const t_pull_coord& pcrd);

/*! \brief Returns the conversion factor from the pull coord internal value unit to the init/rate unit.
 *
 * \param[in] pcrd The pull coordinate to get the conversion factor for.
 * \returns the conversion factor.
 */
double pull_conversion_factor_internal2userinput(const t_pull_coord& pcrd);

/*! \brief Get the value for pull coord coord_ind.
 *
 * \param[in,out] pull        The pull struct.
 * \param[in]     coordIndex  Index of the pull coordinate in the list of coordinates
 * \param[in]     pbc         Information structure about periodicity.
 * \param[in]     t           The time
 * \returns the value of the pull coordinate.
 */
double get_pull_coord_value(pull_t* pull, int coordIndex, const t_pbc& pbc, double t);

/*! \brief Get the value for pull coord coord_ind.
 *
 * Should only be called when time is not allowed to be used as a transformation
 * coordinate variable. This condition is (release) asserted upon.
 *
 * \param[in,out] pull        The pull struct.
 * \param[in]     coordIndex  Index of the pull coordinate in the list of coordinates
 * \param[in]     pbc         Information structure about periodicity.
 * \returns the value of the pull coordinate.
 */
double get_pull_coord_value(pull_t* pull, int coordIndex, const t_pbc& pbc);

/*! \brief Registers the provider of an external potential for a coordinate.
 *
 * This function is only used for checking the consistency of the pull setup.
 * For each pull coordinate of type external-potential, selected by the user
 * in the mdp file, there has to be a module that provides this potential.
 * The module registers itself as the provider by calling this function.
 * The passed \p provider string has to match the string that the user
 * passed with the potential-provider pull coordinate mdp option.
 * This function should be called after init_pull has been called and before
 * pull_potential is called for the first time.
 * This function does many consistency checks and when it returns and the
 * first call to do_potential passes, the pull setup is guaranteed to be
 * correct (unless the module doesn't call apply_external_pull_coord_force
 * every step or calls it with incorrect forces). This registering function
 * will exit with a (release) assertion failure when used incorrely or
 * with a fatal error when the user (mdp) input in inconsistent.
 *
 * Thread-safe for simultaneous registration from multiple threads.
 *
 * \param[in,out] pull         The pull struct.
 * \param[in]     coord_index  The pull coordinate index to register the external potential for.
 * \param[in]     provider     Provider string, should match the potential-provider pull coordinate mdp option.
 */
void register_external_pull_potential(struct pull_t* pull, int coord_index, const char* provider);


/*! \brief Apply forces of an external potential to a pull coordinate.
 *
 * This function applies the external scalar force \p coord_force to
 * the pull coordinate. The corresponding potential energy
 * value should be added to the pull or the module's potential energy term
 * separately by the module itself.
 * This function should be called after pull_potential() has been called and,
 * before calling pull_apply_forces().
 *
 * \param[in,out] pull             The pull struct.
 * \param[in]     coord_index      The pull coordinate index to set the force for.
 * \param[in]     coord_force      The scalar force for the pull coordinate.
 */
void apply_external_pull_coord_force(pull_t* pull, int coord_index, double coord_force);

/*! \brief Set the all the pull forces to zero.
 *
 * \param pull              The pull group.
 */
void clear_pull_forces(pull_t* pull);


/*! \brief Computes the COM pull forces, returns the potential
 *
 * The function computes the COMs of the pull groups and the potentials and forces
 * acting on the pull groups, except for external potential coordinates, which forces
 * are set by calls to \p apply_external_pull_coord_force() after calling this function.
 * To finalize the pull application, a call to \p pull_apply_forces() is required to
 * distribute the forces on the COMs to the atoms.
 *
 * Note: performance global MPI communication, potentially on a subset of the MPI ranks.
 *
 * \param[in,out] pull     The pull struct.
 * \param[in]     masses   Atoms masses.
 * \param[in]     pbc      Information struct about periodicity.
 * \param[in]     mpiComm  Communication object for my group.
 * \param[in]     t        Time.
 * \param[in]     lambda   The value of lambda in FEP calculations.
 * \param[in]     x        Positions.
 * \param[out] dvdlambda   Pull contribution to dV/d(lambda).
 *
 * \returns The pull potential energy.
 */
real pull_potential(pull_t*                        pull,
                    gmx::ArrayRef<const real>      masses,
                    const t_pbc&                   pbc,
                    const gmx::MpiComm&            mpiComm,
                    double                         t,
                    real                           lambda,
                    gmx::ArrayRef<const gmx::RVec> x,
                    real*                          dvdlambda);

/*! \brief Applies the computed COM pull forces to the atoms and accumulates the virial
 *
 * When \p force!=nullptr, distributes the pull force on the COM of each normal pull
 * group to the atoms in the group (using mass weighting).
 *
 * Also performs the recursion for transformation pull coordinates, when present,
 * distributing the force on transformation coordinates to the COM of groups involved.
 *
 * This function should be called after calling \p pull_potential() and also after
 * other modules, e.g. AWH, have called \p apply_external_pull_coord_force().
 *
 * Note: this function is fully local and does not perform MPI communication.
 *
 * \param[in,out] pull     The pull struct.
 * \param[in]     masses   Atoms masses.
 * \param[in]     mpiComm  Communication object for my group.
 * \param[in,out] force    Forces and virial.
 */
void pull_apply_forces(struct pull_t*            pull,
                       gmx::ArrayRef<const real> masses,
                       const gmx::MpiComm&       mpiComm,
                       gmx::ForceWithVirial*     force);

/*! \brief Constrain the coordinates xp in the directions in x
 * and also constrain v when v != NULL.
 *
 * \param[in,out] pull   The pull data.
 * \param[in]     masses Atoms masses.
 * \param[in]     pbc    Information struct about periodicity.
 * \param[in]     mpiComm  Communication object for my group.
 * \param[in]     dt     The time step length.
 * \param[in]     t      The time.
 * \param[in]     x      Positions.
 * \param[in,out] xp     Updated x, can be NULL.
 * \param[in,out] v      Velocities, which may get a pull correction.
 * \param[in,out] vir    The virial, which, if != NULL, gets a pull correction.
 */
void pull_constraint(struct pull_t*            pull,
                     gmx::ArrayRef<const real> masses,
                     const t_pbc&              pbc,
                     const gmx::MpiComm&       mpiComm,
                     double                    dt,
                     double                    t,
                     gmx::ArrayRef<gmx::RVec>  x,
                     gmx::ArrayRef<gmx::RVec>  xp,
                     gmx::ArrayRef<gmx::RVec>  v,
                     tensor                    vir);


/*! \brief Make a selection of the home atoms for all pull groups.
 * Should be called at every domain decomposition.
 *
 * \param mpiComm  Communication object for my group.
 * \param dd       Domain decomposition struct, is nullptr when DD is not in use.
 * \param pull     The pull group.
 */
void dd_make_local_pull_groups(const gmx::MpiComm& mpiComm, gmx_domdec_t* dd, pull_t* pull);


/*! \brief Allocate, initialize and return a pull work struct.
 *
 * \param fplog       General output file, normally md.log.
 * \param pull_params The pull input parameters containing all pull settings.
 * \param ir          The inputrec.
 * \param mtop        The topology of the whole system.
 * \param mpiComm     Communication object for my group.
 * \param dd          Domain decomposition object, is nullptr when DD is not used.
 * \param atomSets    The manager that handles the pull atom sets
 * \param lambda      FEP lambda.
 */
struct pull_t* init_pull(FILE*                     fplog,
                         const pull_params_t*      pull_params,
                         const t_inputrec*         ir,
                         const gmx_mtop_t&         mtop,
                         const gmx::MpiComm&       mpiComm,
                         gmx_domdec_t*             dd,
                         gmx::LocalAtomSetManager* atomSets,
                         real                      lambda);


/*! \brief Close the pull output files and delete pull.
 *
 * \param pull       The pull data structure.
 */
void finish_pull(struct pull_t* pull);


/*! \brief Calculates centers of mass all pull groups.
 *
 * \param[in] mpiComm  Communication object for my group.
 * \param[in] pull     The pull data structure.
 * \param[in] masses   Atoms masses.
 * \param[in] pbc      Information struct about periodicity.
 * \param[in] t        Time, only used for cylinder ref.
 * \param[in] x        The local positions.
 * \param[in,out] xp   Updated x, can be NULL.
 *
 */
void pull_calc_coms(const gmx::MpiComm&            mpiComm,
                    pull_t*                        pull,
                    gmx::ArrayRef<const real>      masses,
                    const t_pbc&                   pbc,
                    double                         t,
                    gmx::ArrayRef<const gmx::RVec> x,
                    gmx::ArrayRef<gmx::RVec>       xp);

/*! \brief Margin for checking pull group PBC distances compared to half the box size */
static constexpr real c_pullGroupPbcMargin = 0.9;
/*! \brief Threshold (as a factor of half the box size) for accepting pull groups without explicitly set refatom */
static constexpr real c_pullGroupSmallGroupThreshold = 0.5;

/*! \brief Checks whether all groups that use a reference atom are within PBC restrictions
 *
 * Groups that use a reference atom for determining PBC should have all their
 * atoms within half the box size from the PBC atom. The box size is used
 * per dimension for rectangular boxes, but can be a combination of
 * dimensions for triclinic boxes, depending on which dimensions are
 * involved in the pull coordinates a group is involved in. A margin is specified
 * to ensure that atoms are not too close to the maximum distance.
 *
 * Should be called without MPI parallelization and after pull_calc_coms()
 * has been called at least once.
 *
 * \param[in] pull       The pull data structure
 * \param[in] x          The coordinates
 * \param[in] pbc        Information struct about periodicity
 * \param[in] pbcMargin  The minimum margin (as a fraction) to half the box size
 * \returns -1 when all groups obey PBC or the first group index that fails PBC
 */
int pullCheckPbcWithinGroups(const pull_t& pull, gmx::ArrayRef<const gmx::RVec> x, const t_pbc& pbc, real pbcMargin);

/*! \brief Checks whether a specific group that uses a reference atom is within PBC restrictions
 *
 * Groups that use a reference atom for determining PBC should have all their
 * atoms within half the box size from the PBC atom. The box size is used
 * per dimension for rectangular boxes, but can be a combination of
 * dimensions for triclinic boxes, depending on which dimensions are
 * involved in the pull coordinates a group is involved in. A margin is specified
 * to ensure that atoms are not too close to the maximum distance. Only one group is
 * checked.
 *
 * Should be called without MPI parallelization and after pull_calc_coms()
 * has been called at least once.
 *
 * \param[in] pull       The pull data structure
 * \param[in] x          The coordinates
 * \param[in] pbc        Information struct about periodicity
 * \param[in] groupNr    The index of the group (in pull.group[]) to check
 * \param[in] pbcMargin  The minimum margin (as a fraction) to half the box size
 * \returns true if the group obeys PBC otherwise false
 */
bool pullCheckPbcWithinGroup(const pull_t&                  pull,
                             gmx::ArrayRef<const gmx::RVec> x,
                             const t_pbc&                   pbc,
                             int                            groupNr,
                             real                           pbcMargin);

/*! \brief Returns if we have pull coordinates with potential pulling.
 *
 * \param[in] pull     The pull data structure.
 */
bool pull_have_potential(const pull_t& pull);


/*! \brief Returns if we have pull coordinates with constraint pulling.
 *
 * \param[in] pull     The pull data structure.
 */
bool pull_have_constraint(const pull_t& pull);

/*! \brief Returns if inputrec has pull coordinates with constraint pulling.
 *
 * \param[in] pullParameters  Pulling input parameters from input record.
 */
bool pull_have_constraint(const pull_params_t& pullParameters);

/*! \brief Returns the maxing distance for pulling
 *
 * For distance geometries, only dimensions with pcrd->params[dim]=1
 * are included in the distance calculation.
 * For directional geometries, only dimensions with pcrd->vec[dim]!=0
 * are included in the distance calculation.
 *
 * \param[in] pcrd Pulling data structure
 * \param[in] pbc  Information on periodic boundary conditions
 * \returns The maximume distance
 */
real max_pull_distance2(const pull_coord_work_t& pcrd, const t_pbc& pbc);

/*! \brief Sets the previous step COM in pull to the current COM, and optionally
 *         updates it in the provided ArrayRef
 *
 * \param[in] pull  The COM pull force calculation data structure
 * \param[in] comPreviousStep  The COM of the previous step of each pull group
 */
void updatePrevStepPullCom(pull_t* pull, std::optional<gmx::ArrayRef<double>> comPreviousStep);

/*! \brief Returns a copy of the previous step pull COM as flat vector
 *
 * Used for modular simulator checkpointing. Allows to keep the
 * implementation details of pull_t hidden from its users.
 *
 * \param[in] pull  The COM pull force calculation data structure
 * \return A copy of the previous step COM
 */
std::vector<double> prevStepPullCom(const pull_t* pull);

/*! \brief Set the previous step pull COM from a flat vector
 *
 * Used to restore modular simulator checkpoints. Allows to keep the
 * implementation details of pull_t hidden from its users.
 *
 * \param[in] pull  The COM pull force calculation data structure
 * \param[in] prevStepPullCom  The previous step COM to set
 */
void setPrevStepPullCom(pull_t* pull, gmx::ArrayRef<const double> prevStepPullCom);

/*! \brief Allocates, initializes and communicates the previous step pull COM (if that option is set to true).
 *
 * If ir->pull->bSetPbcRefToPrevStepCOM is not true nothing is done.
 *
 * \param[in] ir                     The input options/settings of the simulation.
 * \param[in] pull_work              The COM pull force calculation data structure
 * \param[in] masses                 Atoms masses.
 * \param[in] state                  The local (to this rank) state.
 * \param[in] state_global           The global state.
 * \param[in] mpiComm                Communication object for my group.
 * \param[in] startingFromCheckpoint Is the simulation starting from a checkpoint?
 */
void preparePrevStepPullCom(const t_inputrec*         ir,
                            pull_t*                   pull_work,
                            gmx::ArrayRef<const real> masses,
                            t_state*                  state,
                            const t_state*            state_global,
                            const gmx::MpiComm&       mpiComm,
                            bool                      startingFromCheckpoint);

/*! \brief Initializes the COM of the previous step (set to initial COM)
 *
 * \param[in] mpiComm  Communication object for my group.
 * \param[in] pull     The pull data structure.
 * \param[in] masses   Atoms masses.
 * \param[in] pbc      Information struct about periodicity.
 * \param[in] x        The local positions.
 */
void initPullComFromPrevStep(const gmx::MpiComm&            mpiComm,
                             pull_t*                        pull,
                             gmx::ArrayRef<const real>      masses,
                             const t_pbc&                   pbc,
                             gmx::ArrayRef<const gmx::RVec> x);

/*! \brief Initializes the previous step pull COM for new simulations (no reading from checkpoint).
 *
 * \param[in] mpiComm          Communication object for my group.
 * \param[in] pull_work        The COM pull force calculation data structure.
 * \param[in] masses           Atoms masses.
 * \param[in] x                The local positions.
 * \param[in] box              The current box matrix.
 * \param[in] pbcType          The type of periodic boundary conditions.
 * \param[in] comPreviousStep  The COM of the previous step of each pull group.
 */
void preparePrevStepPullComNewSimulation(const gmx::MpiComm&                    mpiComm,
                                         pull_t*                                pull_work,
                                         gmx::ArrayRef<const real>              masses,
                                         gmx::ArrayRef<const gmx::RVec>         x,
                                         const matrix                           box,
                                         PbcType                                pbcType,
                                         std::optional<gmx::ArrayRef<double>>&& comPreviousStep);

#endif
