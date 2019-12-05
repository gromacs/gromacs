/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2005 - 2014, The GROMACS development team.
 * Copyright (c) 2015,2016,2017,2018,2019, by the GROMACS development team, led by
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
/*! \defgroup module_domdec Spatial domain decomposition (for parallelization over MPI)
 * \ingroup group_mdrun
 *
 * \brief Manages the decomposition of the simulation volume over MPI
 * ranks to try to distribute work evenly with minimal communication
 * overheads.
 *
 * \todo Get domdec stuff out of mdtypes/commrec.h
 *
 * \author Berk Hess <hess@kth.se>
 *
 */

/*! \libinternal \file
 *
 * \brief This file declares functions for mdrun to call to manage the
 * details of its domain decomposition.
 *
 * \author Berk Hess <hess@kth.se>
 * \inlibraryapi
 * \ingroup module_domdec
 */

#ifndef GMX_DOMDEC_DOMDEC_H
#define GMX_DOMDEC_DOMDEC_H

#include <vector>

#include "gromacs/math/vectypes.h"
#include "gromacs/utility/arrayref.h"
#include "gromacs/utility/basedefinitions.h"
#include "gromacs/utility/real.h"

struct cginfo_mb_t;
struct gmx_domdec_t;
struct gmx_ddbox_t;
struct gmx_domdec_zones_t;
struct gmx_localtop_t;
struct gmx_mtop_t;
struct gmx_vsite_t;
struct t_block;
struct t_blocka;
struct t_commrec;
struct t_forcerec;
struct t_inputrec;
struct t_mdatoms;
struct t_nrnb;
struct gmx_wallcycle;
class t_state;

namespace gmx
{
class ForceWithShiftForces;
class MDLogger;
class RangePartitioning;
} // namespace gmx

/*! \brief Returns the global topology atom number belonging to local atom index i.
 *
 * This function is intended for writing ASCII output
 * and returns atom numbers starting at 1.
 * When dd=NULL returns i+1.
 */
int ddglatnr(const gmx_domdec_t* dd, int i);

/*! \brief Returns a list of update group partitioning for each molecule type or empty when update groups are not used */
gmx::ArrayRef<const gmx::RangePartitioning> getUpdateGroupingPerMoleculetype(const gmx_domdec_t& dd);

/*! \brief Store the global cg indices of the home cgs in state,
 *
 * This means it can be reset, even after a new DD partitioning.
 */
void dd_store_state(struct gmx_domdec_t* dd, t_state* state);

/*! \brief Returns a pointer to the gmx_domdec_zones_t struct */
struct gmx_domdec_zones_t* domdec_zones(struct gmx_domdec_t* dd);

/*! \brief Returns the range for atoms in zones*/
int dd_numAtomsZones(const gmx_domdec_t& dd);

/*! \brief Returns the number of home atoms */
int dd_numHomeAtoms(const gmx_domdec_t& dd);

/*! \brief Returns the atom range in the local state for atoms that need to be present in mdatoms */
int dd_natoms_mdatoms(const gmx_domdec_t* dd);

/*! \brief Returns the atom range in the local state for atoms involved in virtual sites */
int dd_natoms_vsite(const gmx_domdec_t* dd);

/*! \brief Sets the atom range for atom in the local state for atoms received in constraints communication */
void dd_get_constraint_range(const gmx_domdec_t* dd, int* at_start, int* at_end);

/*! \libinternal \brief Struct for passing around the number of PME domains */
struct NumPmeDomains
{
    int x; //!< The number of PME domains along dimension x
    int y; //!< The number of PME domains along dimension y
};

/*! \brief Returns the number of PME domains, can be called with dd=NULL */
NumPmeDomains getNumPmeDomains(const gmx_domdec_t* dd);

/*! \brief Returns the set of DD ranks that communicate with pme node cr->nodeid */
std::vector<int> get_pme_ddranks(const t_commrec* cr, int pmenodeid);

/*! \brief Returns the maximum shift for coordinate communication in PME, dim x */
int dd_pme_maxshift_x(const gmx_domdec_t* dd);

/*! \brief Returns the maximum shift for coordinate communication in PME, dim y */
int dd_pme_maxshift_y(const gmx_domdec_t* dd);

/*! \brief Return whether constraints, not including settles, cross domain boundaries */
bool ddHaveSplitConstraints(const gmx_domdec_t& dd);

/*! \brief Return whether update groups are used */
bool ddUsesUpdateGroups(const gmx_domdec_t& dd);

/*! \brief Return whether the DD has a single dimension with a single pulse
 *
 * The GPU halo exchange code requires a 1D single-pulse DD, and its
 * setup code can use the returned value to understand what it should
 * do. */
bool is1DAnd1PulseDD(const gmx_domdec_t& dd);

/*! \brief Initialize data structures for bonded interactions */
void dd_init_bondeds(FILE*              fplog,
                     gmx_domdec_t*      dd,
                     const gmx_mtop_t*  mtop,
                     const gmx_vsite_t* vsite,
                     const t_inputrec*  ir,
                     gmx_bool           bBCheck,
                     cginfo_mb_t*       cginfo_mb);

/*! \brief Returns whether molecules are always whole, i.e. not broken by PBC */
bool dd_moleculesAreAlwaysWhole(const gmx_domdec_t& dd);

/*! \brief Returns if we need to do pbc for calculating bonded interactions */
gmx_bool dd_bonded_molpbc(const gmx_domdec_t* dd, int ePBC);

/*! \brief Change the DD non-bonded communication cut-off.
 *
 * This could fail when trying to increase the cut-off,
 * then FALSE will be returned and the cut-off is not modified.
 *
 * \param[in] cr               Communication recrod
 * \param[in] box              Box matrix, used for computing the dimensions of the system
 * \param[in] x                Position vector, used for computing the dimensions of the system
 * \param[in] cutoffRequested  The requested atom to atom cut-off distance, usually the pair-list cutoff distance
 */
gmx_bool change_dd_cutoff(t_commrec* cr, const matrix box, gmx::ArrayRef<const gmx::RVec> x, real cutoffRequested);

/*! \brief Set up communication for averaging GPU wait times over domains
 *
 * When domains (PP MPI ranks) share a GPU, the individual GPU wait times
 * are meaningless, as it depends on the order in which tasks on the same
 * GPU finish. Therefore there wait times need to be averaged over the ranks
 * sharing the same GPU. This function sets up the communication for that.
 */
void dd_setup_dlb_resource_sharing(const t_commrec* cr, int gpu_id);

/*! \brief Cycle counter indices used internally in the domain decomposition */
enum
{
    ddCyclStep,
    ddCyclPPduringPME,
    ddCyclF,
    ddCyclWaitGPU,
    ddCyclPME,
    ddCyclNr
};

/*! \brief Add the wallcycle count to the DD counter */
void dd_cycles_add(const gmx_domdec_t* dd, float cycles, int ddCycl);

/*! \brief Communicate the coordinates to the neighboring cells and do pbc. */
void dd_move_x(struct gmx_domdec_t* dd, const matrix box, gmx::ArrayRef<gmx::RVec> x, gmx_wallcycle* wcycle);

/*! \brief Sum the forces over the neighboring cells.
 *
 * When fshift!=NULL the shift forces are updated to obtain
 * the correct virial from the single sum including f.
 */
void dd_move_f(struct gmx_domdec_t* dd, gmx::ForceWithShiftForces* forceWithShiftForces, gmx_wallcycle* wcycle);

/*! \brief Communicate a real for each atom to the neighboring cells. */
void dd_atom_spread_real(struct gmx_domdec_t* dd, real v[]);

/*! \brief Sum the contributions to a real for each atom over the neighboring cells. */
void dd_atom_sum_real(struct gmx_domdec_t* dd, real v[]);

/*! \brief Reset all the statistics and counters for total run counting */
void reset_dd_statistics_counters(struct gmx_domdec_t* dd);

/* In domdec_con.c */

/*! \brief Communicates the virtual site forces, reduces the shift forces when \p fshift != NULL */
void dd_move_f_vsites(struct gmx_domdec_t* dd, rvec* f, rvec* fshift);

/*! \brief Clears the forces for virtual sites */
void dd_clear_f_vsites(struct gmx_domdec_t* dd, rvec* f);

/*! \brief Move x0 and also x1 if x1!=NULL. bX1IsCoord tells if to do PBC on x1 */
void dd_move_x_constraints(struct gmx_domdec_t* dd, const matrix box, rvec* x0, rvec* x1, gmx_bool bX1IsCoord);

/*! \brief Communicates the coordinates involved in virtual sites */
void dd_move_x_vsites(struct gmx_domdec_t* dd, const matrix box, rvec* x);

/*! \brief Returns the local atom count array for all constraints
 *
 * The local atom count for a constraint, possible values 2/1/0, is needed
 * to avoid not/double-counting contributions linked to the Lagrange
 * multiplier, such as the virial and free-energy derivatives.
 *
 * \note When \p dd = nullptr, an empty reference is returned.
 */
gmx::ArrayRef<const int> dd_constraints_nlocalatoms(const gmx_domdec_t* dd);

/* In domdec_top.c */

/*! \brief Print error output when interactions are missing */
[[noreturn]] void dd_print_missing_interactions(const gmx::MDLogger&  mdlog,
                                                t_commrec*            cr,
                                                int                   local_count,
                                                const gmx_mtop_t*     top_global,
                                                const gmx_localtop_t* top_local,
                                                const rvec*           x,
                                                const matrix          box);

/*! \brief Generate and store the reverse topology */
void dd_make_reverse_top(FILE*              fplog,
                         gmx_domdec_t*      dd,
                         const gmx_mtop_t*  mtop,
                         const gmx_vsite_t* vsite,
                         const t_inputrec*  ir,
                         gmx_bool           bBCheck);

/*! \brief Generate the local topology and virtual site data */
void dd_make_local_top(struct gmx_domdec_t*       dd,
                       struct gmx_domdec_zones_t* zones,
                       int                        npbcdim,
                       matrix                     box,
                       rvec                       cellsize_min,
                       const ivec                 npulse,
                       t_forcerec*                fr,
                       rvec*                      cgcm_or_x,
                       const gmx_mtop_t&          top,
                       gmx_localtop_t*            ltop);

/*! \brief Sort ltop->ilist when we are doing free energy. */
void dd_sort_local_top(gmx_domdec_t* dd, const t_mdatoms* mdatoms, gmx_localtop_t* ltop);

/*! \brief Initialize local topology
 *
 * \param[in] top_global Reference to global topology.
 * \param[in,out] top Pointer to new local topology
 */
void dd_init_local_top(const gmx_mtop_t& top_global, gmx_localtop_t* top);

/*! \brief Construct local state */
void dd_init_local_state(struct gmx_domdec_t* dd, const t_state* state_global, t_state* local_state);

/*! \brief Generate a list of links between atoms that are linked by bonded interactions
 *
 * Also stores whether atoms are linked in \p cginfo_mb.
 */
t_blocka* makeBondedLinks(const gmx_mtop_t* mtop, cginfo_mb_t* cginfo_mb);

/*! \brief Calculate the maximum distance involved in 2-body and multi-body bonded interactions */
void dd_bonded_cg_distance(const gmx::MDLogger& mdlog,
                           const gmx_mtop_t*    mtop,
                           const t_inputrec*    ir,
                           const rvec*          x,
                           const matrix         box,
                           gmx_bool             bBCheck,
                           real*                r_2b,
                           real*                r_mb);

#endif
