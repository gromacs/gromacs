/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2004, The GROMACS development team.
 * Copyright (c) 2013,2014,2015,2016,2017,2018, by the GROMACS development team, led by
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

#include "gromacs/math/vectypes.h"
#include "gromacs/mdlib/mdrun.h"
#include "gromacs/mdtypes/pull-params.h"
#include "gromacs/utility/arrayref.h"
#include "gromacs/utility/basedefinitions.h"
#include "gromacs/utility/real.h"

struct ContinuationOptions;
struct gmx_mtop_t;
struct gmx_output_env_t;
struct pull_coord_work_t;
struct pull_params_t;
struct t_commrec;
struct t_filenm;
struct t_inputrec;
struct t_mdatoms;
struct t_pbc;
class t_state;

namespace gmx
{
class ForceWithVirial;
class LocalAtomSetManager;
}

/*! \brief Returns if the pull coordinate is an angle
 *
 * \param[in] pcrd The pull coordinate to query the type for.
 * \returns a boolean telling if the coordinate is of angle type.
 */
bool pull_coordinate_is_angletype(const t_pull_coord *pcrd);

/*! \brief Returns the units of the pull coordinate.
 *
 * \param[in] pcrd The pull coordinate to query the units for.
 * \returns a string with the units of the coordinate.
 */
const char *pull_coordinate_units(const t_pull_coord *pcrd);

/*! \brief Returns the conversion factor from the pull coord init/rate unit to internal value unit.
 *
 * \param[in] pcrd The pull coordinate to get the conversion factor for.
 * \returns the conversion factor.
 */
double pull_conversion_factor_userinput2internal(const t_pull_coord *pcrd);

/*! \brief Returns the conversion factor from the pull coord internal value unit to the init/rate unit.
 *
 * \param[in] pcrd The pull coordinate to get the conversion factor for.
 * \returns the conversion factor.
 */
double pull_conversion_factor_internal2userinput(const t_pull_coord *pcrd);

/*! \brief Get the value for pull coord coord_ind.
 *
 * \param[in,out] pull      The pull struct.
 * \param[in]     coord_ind Number of the pull coordinate.
 * \param[in]     pbc       Information structure about periodicity.
 * \returns the value of the pull coordinate.
 */
double get_pull_coord_value(struct pull_t      *pull,
                            int                 coord_ind,
                            const struct t_pbc *pbc);

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
void register_external_pull_potential(struct pull_t *pull,
                                      int            coord_index,
                                      const char    *provider);


/*! \brief Apply forces of an external potential to a pull coordinate.
 *
 * This function applies the external scalar force \p coord_force to
 * the pull coordinate, distributing it over the atoms in the groups
 * involved in the pull coordinate. The corresponding potential energy
 * value should be added to the pull or the module's potential energy term
 * separately by the module itself.
 * This function should be called after pull_potential has been called and,
 * obviously, before the coordinates are updated uses the forces.
 *
 * \param[in,out] pull             The pull struct.
 * \param[in]     coord_index      The pull coordinate index to set the force for.
 * \param[in]     coord_force      The scalar force for the pull coordinate.
 * \param[in]     mdatoms          Atom properties, only masses are used.
 * \param[in,out] forceWithVirial  Force and virial buffers.
 */
void apply_external_pull_coord_force(struct pull_t        *pull,
                                     int                   coord_index,
                                     double                coord_force,
                                     const t_mdatoms      *mdatoms,
                                     gmx::ForceWithVirial *forceWithVirial);


/*! \brief Set the all the pull forces to zero.
 *
 * \param pull              The pull group.
 */
void clear_pull_forces(struct pull_t *pull);


/*! \brief Determine the COM pull forces and add them to f, return the potential
 *
 * \param[in,out] pull   The pull struct.
 * \param[in]     md     All atoms.
 * \param[in]     pbc    Information struct about periodicity.
 * \param[in]     cr     Struct for communication info.
 * \param[in]     t      Time.
 * \param[in]     lambda The value of lambda in FEP calculations.
 * \param[in]     x      Positions.
 * \param[in,out] force  Forces and virial.
 * \param[out] dvdlambda Pull contribution to dV/d(lambda).
 *
 * \returns The pull potential energy.
 */
real pull_potential(struct pull_t *pull, const t_mdatoms *md, struct t_pbc *pbc,
                    const t_commrec *cr, double t, real lambda,
                    const rvec *x, gmx::ForceWithVirial *force, real *dvdlambda);


/*! \brief Constrain the coordinates xp in the directions in x
 * and also constrain v when v != NULL.
 *
 * \param[in,out] pull   The pull data.
 * \param[in]     md     All atoms.
 * \param[in]     pbc    Information struct about periodicity.
 * \param[in]     cr     Struct for communication info.
 * \param[in]     dt     The time step length.
 * \param[in]     t      The time.
 * \param[in]     x      Positions.
 * \param[in,out] xp     Updated x, can be NULL.
 * \param[in,out] v      Velocities, which may get a pull correction.
 * \param[in,out] vir    The virial, which, if != NULL, gets a pull correction.
 */
void pull_constraint(struct pull_t *pull, const t_mdatoms *md, struct t_pbc *pbc,
                     const t_commrec *cr, double dt, double t,
                     rvec *x, rvec *xp, rvec *v, tensor vir);


/*! \brief Make a selection of the home atoms for all pull groups.
 * Should be called at every domain decomposition.
 *
 * \param cr             Structure for communication info.
 * \param pull           The pull group.
 */
void dd_make_local_pull_groups(const t_commrec *cr, struct pull_t *pull);


/*! \brief Allocate, initialize and return a pull work struct.
 *
 * \param fplog       General output file, normally md.log.
 * \param pull_params The pull input parameters containing all pull settings.
 * \param ir          The inputrec.
 * \param mtop        The topology of the whole system.
 * \param cr          Struct for communication info.
 * \param atomSets    The manager that handles the pull atom sets
 * \param lambda      FEP lambda.
 */
struct pull_t *init_pull(FILE                      *fplog,
                         const pull_params_t       *pull_params,
                         const t_inputrec          *ir,
                         const gmx_mtop_t          *mtop,
                         const t_commrec           *cr,
                         gmx::LocalAtomSetManager  *atomSets,
                         real                       lambda);


/*! \brief Close the pull output files and delete pull.
 *
 * \param pull       The pull data structure.
 */
void finish_pull(struct pull_t *pull);


/*! \brief Calculates centers of mass all pull groups.
 *
 * \param[in] cr       Struct for communication info.
 * \param[in] pull     The pull data structure.
 * \param[in] md       All atoms.
 * \param[in] pbc      Information struct about periodicity.
 * \param[in] t        Time, only used for cylinder ref.
 * \param[in] x        The local positions.
 * \param[in,out] xp   Updated x, can be NULL.
 *
 */
void pull_calc_coms(const t_commrec *cr,
                    pull_t          *pull,
                    const t_mdatoms *md,
                    t_pbc           *pbc,
                    double           t,
                    const rvec       x[],
                    rvec            *xp);

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
int pullCheckPbcWithinGroups(const pull_t &pull,
                             const rvec   *x,
                             const t_pbc  &pbc,
                             real          pbcMargin);

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
bool pullCheckPbcWithinGroup(const pull_t                  &pull,
                             gmx::ArrayRef<const gmx::RVec> x,
                             const t_pbc                   &pbc,
                             int                            groupNr,
                             real                           pbcMargin);

/*! \brief Returns if we have pull coordinates with potential pulling.
 *
 * \param[in] pull     The pull data structure.
 */
gmx_bool pull_have_potential(const struct pull_t *pull);


/*! \brief Returns if we have pull coordinates with constraint pulling.
 *
 * \param[in] pull     The pull data structure.
 */
gmx_bool pull_have_constraint(const struct pull_t *pull);

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
real max_pull_distance2(const pull_coord_work_t *pcrd,
                        const t_pbc             *pbc);

/*! \brief Copies the COM from the previous step of all pull groups to the checkpoint state container
 *
 * \param[in]   pull  The COM pull force calculation data structure
 * \param[in]   state The global state container
 */
void setStatePrevStepPullCom(const struct pull_t *pull, t_state *state);

/*! \brief Copies the pull group COM of the previous step from the checkpoint state to the pull state
 *
 * \param[in]   pull  The COM pull force calculation data structure
 * \param[in]   state The global state container
 */
void setPrevStepPullComFromState(struct pull_t *pull, const t_state *state);

/*! \brief Sets the previous step COM to the current COM
 *
 * \param[in]   pull The COM pull force calculation data structure
 */
void updatePrevStepCom(struct pull_t *pull);

/*! \brief Resizes the vector, in the state container, containing the COMs from the previous step
 *
 * \param[in]   state The global state container
 * \param[in]   pull  The COM pull force calculation data structure
 */
void allocStatePrevStepPullCom(t_state *state, pull_t *pull);

/*! \brief Initializes the COM of the previous step (set to initial COM)
 *
 * \param[in] cr       Struct for communication info.
 * \param[in] pull     The pull data structure.
 * \param[in] md       All atoms.
 * \param[in] pbc      Information struct about periodicity.
 * \param[in] x        The local positions.
 */
void initPullComFromPrevStep(const t_commrec *cr,
                             pull_t          *pull,
                             const t_mdatoms *md,
                             t_pbc           *pbc,
                             const rvec       x[]);

#endif
