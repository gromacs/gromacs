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
#include "gromacs/mdtypes/pull-params.h"
#include "gromacs/pulling/pull_internal.h"
#include "gromacs/utility/basedefinitions.h"
#include "gromacs/utility/real.h"

struct ContinuationOptions;
struct gmx_mtop_t;
struct gmx_output_env_t;
struct pull_params_t;
struct t_commrec;
struct t_filenm;
struct t_inputrec;
struct t_mdatoms;
struct t_pbc;

namespace gmx
{
class ForceWithVirial;
}

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
 * \param[out]    value     The value of the pull coordinate.
 */
void get_pull_coord_value(struct pull_t      *pull,
                          int                 coord_ind,
                          const struct t_pbc *pbc,
                          double             *value);


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
 * \param[in,out] pull           The pull struct.
 * \param[in]     coord_index    The pull coordinate index to set the force for.
 * \param[in]     coord_force    The scalar force for the pull coordinate.
 * \param[in]     mdatoms        Atom properties, only masses are used.
 * \param[in,out] force          The force buffer.
 * \param[in,out] virial         The virial, can be NULL.
 */
void apply_external_pull_coord_force(struct pull_t   *pull,
                                     int              coord_index,
                                     double           coord_force,
                                     const t_mdatoms *mdatoms,
                                     rvec            *force,
                                     tensor           virial);


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
real pull_potential(struct pull_t *pull, t_mdatoms *md, struct t_pbc *pbc,
                    t_commrec *cr, double t, real lambda,
                    rvec *x, gmx::ForceWithVirial *force, real *dvdlambda);


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
void pull_constraint(struct pull_t *pull, t_mdatoms *md, struct t_pbc *pbc,
                     t_commrec *cr, double dt, double t,
                     rvec *x, rvec *xp, rvec *v, tensor vir);


/*! \brief Make a selection of the home atoms for all pull groups.
 * Should be called at every domain decomposition.
 *
 * \param cr             Structure for communication info.
 * \param pull           The pull group.
 * \param md             All atoms.
 */
void dd_make_local_pull_groups(t_commrec *cr,
                               struct pull_t *pull, t_mdatoms *md);


/*! \brief Allocate, initialize and return a pull work struct.
 *
 * \param fplog       General output file, normally md.log.
 * \param pull_params The pull input parameters containing all pull settings.
 * \param ir          The inputrec.
 * \param nfile       Number of files.
 * \param fnm         Standard filename struct.
 * \param mtop        The topology of the whole system.
 * \param cr          Struct for communication info.
 * \param oenv        Output options.
 * \param lambda      FEP lambda.
 * \param bOutFile    Open output files?
 * \param continuationOptions  Options for continuing from checkpoint file
 */
struct pull_t *init_pull(FILE                      *fplog,
                         const pull_params_t       *pull_params,
                         const t_inputrec          *ir,
                         int                        nfile,
                         const t_filenm             fnm[],
                         const gmx_mtop_t          *mtop,
                         t_commrec                * cr,
                         const gmx_output_env_t    *oenv,
                         real                       lambda,
                         gmx_bool                   bOutFile,
                         const ContinuationOptions &continuationOptions);


/*! \brief Close the pull output files.
 *
 * \param pull       The pull group.
 */
void finish_pull(struct pull_t *pull);


/*! \brief Print the pull output (x and/or f)
 *
 * \param pull     The pull data structure.
 * \param step     Time step number.
 * \param time     Time.
 */
void pull_print_output(struct pull_t *pull, gmx_int64_t step, double time);


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
void pull_calc_coms(t_commrec        *cr,
                    struct pull_t    *pull,
                    t_mdatoms        *md,
                    struct t_pbc     *pbc,
                    double            t,
                    rvec              x[],
                    rvec             *xp);


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

#endif
