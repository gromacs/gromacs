/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2005,2006,2007,2008,2009,2010,2012,2013,2014,2015,2016,2017, by the GROMACS development team, led by
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

#include <stdio.h>

#include "gromacs/gmxlib/nrnb.h"
#include "gromacs/math/vectypes.h"
#include "gromacs/mdlib/vsite.h"
#include "gromacs/mdtypes/forcerec.h"
#include "gromacs/mdtypes/mdatom.h"
#include "gromacs/timing/wallcycle.h"
#include "gromacs/topology/block.h"
#include "gromacs/topology/idef.h"
#include "gromacs/topology/topology.h"
#include "gromacs/utility/basedefinitions.h"
#include "gromacs/utility/real.h"

struct gmx_domdec_t;
struct gmx_ddbox_t;
struct gmx_domdec_zones_t;
struct MdrunOptions;
struct t_commrec;
struct t_inputrec;
class t_state;

/*! \brief Returns the global topology atom number belonging to local atom index i.
 *
 * This function is intended for writing ASCII output
 * and returns atom numbers starting at 1.
 * When dd=NULL returns i+1.
 */
int ddglatnr(const gmx_domdec_t *dd, int i);

/*! \brief Return a block struct for the charge groups of the whole system */
t_block *dd_charge_groups_global(struct gmx_domdec_t *dd);

/*! \brief Store the global cg indices of the home cgs in state,
 *
 * This means it can be reset, even after a new DD partitioning.
 */
void dd_store_state(struct gmx_domdec_t *dd, t_state *state);

/*! \brief Returns a pointer to the gmx_domdec_zones_t struct */
struct gmx_domdec_zones_t *domdec_zones(struct gmx_domdec_t *dd);

/*! \brief Sets the j-charge-group range for i-charge-group \p icg */
void dd_get_ns_ranges(const gmx_domdec_t *dd, int icg,
                      int *jcg0, int *jcg1, ivec shift0, ivec shift1);

/*! \brief Returns the atom range in the local state for atoms that need to be present in mdatoms */
int dd_natoms_mdatoms(const gmx_domdec_t *dd);

/*! \brief Returns the atom range in the local state for atoms involved in virtual sites */
int dd_natoms_vsite(const gmx_domdec_t *dd);

/*! \brief Sets the atom range for atom in the local state for atoms received in constraints communication */
void dd_get_constraint_range(const gmx_domdec_t *dd,
                             int *at_start, int *at_end);

/*! \brief Get the number of PME nodes along x and y, can be called with dd=NULL */
void get_pme_nnodes(const struct gmx_domdec_t *dd,
                    int *npmenodes_x, int *npmenodes_y);

/*! \brief Returns the set of DD nodes that communicate with pme node cr->nodeid */
void get_pme_ddnodes(struct t_commrec *cr, int pmenodeid,
                     int *nmy_ddnodes, int **my_ddnodes, int *node_peer);

/*! \brief Returns the maximum shift for coordinate communication in PME, dim x */
int dd_pme_maxshift_x(const gmx_domdec_t *dd);

/*! \brief Returns the maximum shift for coordinate communication in PME, dim y */
int dd_pme_maxshift_y(const gmx_domdec_t *dd);

/*! \brief The options for the domain decomposition MPI task ordering. */
enum class DdRankOrder
{
    select,     //!< First value (needed to cope with command-line parsing)
    interleave, //!< Interleave the PP and PME ranks
    pp_pme,     //!< First all PP ranks, all PME rank at the end
    cartesian,  //!< Use Cartesian communicators for PP, PME and PP-PME
    nr          //!< The number of options
};

/*! \brief The options for the dynamic load balancing. */
enum class DlbOption
{
    select,           //!< First value (needed to cope with command-line parsing)
    turnOnWhenUseful, //!< Turn on DLB when we think it would improve performance
    no,               //!< Never turn on DLB
    yes,              //!< Turn on DLB from the start and keep it on
    nr                //!< The number of options
};

/*! \libinternal \brief Structure containing all (command line) options for the domain decomposition */
struct DomdecOptions
{
    /*! \brief Constructor */
    DomdecOptions();

    //! If true, check that all bonded interactions have been assigned to exactly one domain/rank.
    gmx_bool          checkBondedInteractions;
    //! If true, don't communicate all atoms between the non-bonded cut-off and the larger bonded cut-off, but only those that have non-local bonded interactions. This significantly reduces the communication volume.
    gmx_bool          useBondedCommunication;
    //! The domain decomposition grid cell count, 0 means let domdec choose based on the number of ranks.
    ivec              numCells;
    //! The number of separate PME ranks requested, -1 = auto.
    int               numPmeRanks;
    //! Ordering of the PP and PME ranks, values from enum above.
    DdRankOrder       rankOrder;
    //! The minimum communication range, used for extended the communication range for bonded interactions (nm).
    real              minimumCommunicationRange;
    //! Communication range for atom involved in constraints (P-LINCS) (nm).
    real              constraintCommunicationRange;
    //! Dynamic load balancing option, values from enum above.
    DlbOption         dlbOption;
    /*! \brief Fraction in (0,1) by whose reciprocal the initial
     * DD cell size will be increased in order to provide a margin
     * in which dynamic load balancing can act, while preserving
     * the minimum cell size. */
    real              dlbScaling;
    //! String containing a vector of the relative sizes in the x direction of the corresponding DD cells.
    const char       *cellSizeX;
    //! String containing a vector of the relative sizes in the y direction of the corresponding DD cells.
    const char       *cellSizeY;
    //! String containing a vector of the relative sizes in the z direction of the corresponding DD cells.
    const char       *cellSizeZ;
};

/*! \brief Initialized the domain decomposition, chooses the DD grid and PME ranks, return the DD struct */
gmx_domdec_t *init_domain_decomposition(FILE                *fplog,
                                        t_commrec           *cr,
                                        const DomdecOptions &options,
                                        const MdrunOptions  &mdrunOptions,
                                        const gmx_mtop_t    *mtop,
                                        const t_inputrec    *ir,
                                        const matrix         box,
                                        const rvec          *xGlobal,
                                        gmx_ddbox_t         *ddbox,
                                        int                 *npme_x,
                                        int                 *npme_y);

/*! \brief Initialize data structures for bonded interactions */
void dd_init_bondeds(FILE              *fplog,
                     gmx_domdec_t      *dd,
                     const gmx_mtop_t  *mtop,
                     const gmx_vsite_t *vsite,
                     const t_inputrec  *ir,
                     gmx_bool           bBCheck,
                     cginfo_mb_t       *cginfo_mb);

/*! \brief Returns if we need to do pbc for calculating bonded interactions */
gmx_bool dd_bonded_molpbc(const gmx_domdec_t *dd, int ePBC);

/*! \brief Change the DD non-bonded communication cut-off.
 *
 * This could fail when trying to increase the cut-off,
 * then FALSE will be returned and the cut-off is not modified.
 */
gmx_bool change_dd_cutoff(struct t_commrec *cr,
                          t_state *state, const t_inputrec *ir,
                          real cutoff_req );

/*! \brief Limit DLB to preserve the option of returning to the current cut-off.
 *
 * Domain boundary changes due to the DD dynamic load balancing can limit
 * the cut-off distance that can be set in change_dd_cutoff. This function
 * sets/changes the DLB limit such that using the passed (pair-list) cut-off
 * should still be possible after subsequently setting a shorter cut-off
 * with change_dd_cutoff.
 */
void set_dd_dlb_max_cutoff(struct t_commrec *cr, real cutoff);

/*! \brief Return if we are currently using dynamic load balancing */
gmx_bool dd_dlb_is_on(const struct gmx_domdec_t *dd);

/*! \brief Return if the DLB lock is set */
gmx_bool dd_dlb_is_locked(const struct gmx_domdec_t *dd);

/*! \brief Set a lock such that with DLB=auto DLB cannot get turned on */
void dd_dlb_lock(struct gmx_domdec_t *dd);

/*! \brief Clear a lock such that with DLB=auto DLB may get turned on later */
void dd_dlb_unlock(struct gmx_domdec_t *dd);

/*! \brief Set up communication for averaging GPU wait times over domains
 *
 * When domains (PP MPI ranks) share a GPU, the individual GPU wait times
 * are meaningless, as it depends on the order in which tasks on the same
 * GPU finish. Therefore there wait times need to be averaged over the ranks
 * sharing the same GPU. This function sets up the communication for that.
 */
void dd_setup_dlb_resource_sharing(t_commrec           *cr,
                                   int                  gpu_id);

/*! \brief Collects local rvec arrays \p lv to \p v on the master rank */
void dd_collect_vec(struct gmx_domdec_t    *dd,
                    const t_state          *state_local,
                    const PaddedRVecVector *lv,
                    rvec                   *v);

/*! \brief Collects local rvec arrays \p lv to \p v on the master rank */
void dd_collect_vec(struct gmx_domdec_t    *dd,
                    const t_state          *state_local,
                    const PaddedRVecVector *lv,
                    PaddedRVecVector       *v);

/*! \brief Collects the local state \p state_local to \p state on the master rank */
void dd_collect_state(struct gmx_domdec_t *dd,
                      const t_state *state_local, t_state *state);

/*! \brief Cycle counter indices used internally in the domain decomposition */
enum {
    ddCyclStep, ddCyclPPduringPME, ddCyclF, ddCyclWaitGPU, ddCyclPME, ddCyclNr
};

/*! \brief Add the wallcycle count to the DD counter */
void dd_cycles_add(const gmx_domdec_t *dd, float cycles, int ddCycl);

/*! \brief Start the force flop count */
void dd_force_flop_start(struct gmx_domdec_t *dd, t_nrnb *nrnb);

/*! \brief Stop the force flop count */
void dd_force_flop_stop(struct gmx_domdec_t *dd, t_nrnb *nrnb);

/*! \brief Return the PME/PP force load ratio, or -1 if nothing was measured.
 *
 * Should only be called on the DD master node.
 */
float dd_pme_f_ratio(struct gmx_domdec_t *dd);

/*! \brief Communicate the coordinates to the neighboring cells and do pbc. */
void dd_move_x(struct gmx_domdec_t *dd, matrix box, rvec x[]);

/*! \brief Sum the forces over the neighboring cells.
 *
 * When fshift!=NULL the shift forces are updated to obtain
 * the correct virial from the single sum including f.
 */
void dd_move_f(struct gmx_domdec_t *dd, rvec f[], rvec *fshift);

/*! \brief Communicate a real for each atom to the neighboring cells. */
void dd_atom_spread_real(struct gmx_domdec_t *dd, real v[]);

/*! \brief Sum the contributions to a real for each atom over the neighboring cells. */
void dd_atom_sum_real(struct gmx_domdec_t *dd, real v[]);

/*! \brief Partition the system over the nodes.
 *
 * step is only used for printing error messages.
 * If bMasterState==TRUE then state_global from the master node is used,
 * else state_local is redistributed between the nodes.
 * When f!=NULL, *f will be reallocated to the size of state_local.
 */
void dd_partition_system(FILE                *fplog,
                         gmx_int64_t          step,
                         t_commrec           *cr,
                         gmx_bool             bMasterState,
                         int                  nstglobalcomm,
                         t_state             *state_global,
                         const gmx_mtop_t    *top_global,
                         const t_inputrec    *ir,
                         t_state             *state_local,
                         PaddedRVecVector    *f,
                         t_mdatoms           *mdatoms,
                         gmx_localtop_t      *top_local,
                         t_forcerec          *fr,
                         gmx_vsite_t         *vsite,
                         struct gmx_constr   *constr,
                         t_nrnb              *nrnb,
                         gmx_wallcycle_t      wcycle,
                         gmx_bool             bVerbose);

/*! \brief Reset all the statistics and counters for total run counting */
void reset_dd_statistics_counters(struct gmx_domdec_t *dd);

/*! \brief Print statistics for domain decomposition communication */
void print_dd_statistics(struct t_commrec *cr, const t_inputrec *ir, FILE *fplog);

/* In domdec_con.c */

/*! \brief Communicates the virtual site forces, reduces the shift forces when \p fshift != NULL */
void dd_move_f_vsites(struct gmx_domdec_t *dd, rvec *f, rvec *fshift);

/*! \brief Clears the forces for virtual sites */
void dd_clear_f_vsites(struct gmx_domdec_t *dd, rvec *f);

/*! \brief Move x0 and also x1 if x1!=NULL. bX1IsCoord tells if to do PBC on x1 */
void dd_move_x_constraints(struct gmx_domdec_t *dd, matrix box,
                           rvec *x0, rvec *x1, gmx_bool bX1IsCoord);

/*! \brief Communicates the coordinates involved in virtual sites */
void dd_move_x_vsites(struct gmx_domdec_t *dd, matrix box, rvec *x);

/*! \brief Returns the local atom count array for all constraints
 *
 * The local atom count for a constraint, possible values 2/1/0, is needed
 * to avoid not/double-counting contributions linked to the Lagrange
 * multiplier, such as the virial and free-energy derivatives.
 */
int *dd_constraints_nlocalatoms(struct gmx_domdec_t *dd);

/* In domdec_top.c */

/*! \brief Print error output when interactions are missing */
void dd_print_missing_interactions(FILE *fplog, struct t_commrec *cr,
                                   int local_count,
                                   const gmx_mtop_t *top_global,
                                   const gmx_localtop_t *top_local,
                                   t_state *state_local);

/*! \brief Generate and store the reverse topology */
void dd_make_reverse_top(FILE *fplog,
                         gmx_domdec_t *dd, const gmx_mtop_t *mtop,
                         const gmx_vsite_t *vsite,
                         const t_inputrec *ir, gmx_bool bBCheck);

/*! \brief Store the local charge group index in \p lcgs */
void dd_make_local_cgs(struct gmx_domdec_t *dd, t_block *lcgs);

/*! \brief Generate the local topology and virtual site data */
void dd_make_local_top(struct gmx_domdec_t *dd, struct gmx_domdec_zones_t *zones,
                       int npbcdim, matrix box,
                       rvec cellsize_min, ivec npulse,
                       t_forcerec *fr,
                       rvec *cgcm_or_x,
                       gmx_vsite_t *vsite,
                       const gmx_mtop_t *top, gmx_localtop_t *ltop);

/*! \brief Sort ltop->ilist when we are doing free energy. */
void dd_sort_local_top(gmx_domdec_t *dd, const t_mdatoms *mdatoms,
                       gmx_localtop_t *ltop);

/*! \brief Construct local topology */
gmx_localtop_t *dd_init_local_top(const gmx_mtop_t *top_global);

/*! \brief Construct local state */
void dd_init_local_state(struct gmx_domdec_t *dd,
                         t_state *state_global, t_state *local_state);

/*! \brief Generate a list of links between charge groups that are linked by bonded interactions */
t_blocka *make_charge_group_links(const gmx_mtop_t *mtop, gmx_domdec_t *dd,
                                  cginfo_mb_t *cginfo_mb);

/*! \brief Calculate the maximum distance involved in 2-body and multi-body bonded interactions */
void dd_bonded_cg_distance(FILE *fplog, const gmx_mtop_t *mtop,
                           const t_inputrec *ir,
                           const rvec *x, const matrix box,
                           gmx_bool bBCheck,
                           real *r_2b, real *r_mb);

/*! \brief Dump a pdb file with the current DD home + communicated atoms.
 *
 * When natoms=-1, dump all known atoms.
 */
void write_dd_pdb(const char *fn, gmx_int64_t step, const char *title,
                  const gmx_mtop_t *mtop,
                  t_commrec *cr,
                  int natoms, rvec x[], matrix box);


/* In domdec_setup.c */

/*! \brief Returns the volume fraction of the system that is communicated */
real comm_box_frac(const ivec dd_nc, real cutoff, const gmx_ddbox_t *ddbox);

/*! \brief Determines the optimal DD cell setup dd->nc and possibly npmenodes
 * for the system.
 *
 * On the master node returns the actual cellsize limit used.
 */
real dd_choose_grid(FILE *fplog,
                    t_commrec *cr, gmx_domdec_t *dd,
                    const t_inputrec *ir,
                    const gmx_mtop_t *mtop,
                    const matrix box, const gmx_ddbox_t *ddbox,
                    int nPmeRanks,
                    gmx_bool bDynLoadBal, real dlb_scale,
                    real cellsize_limit, real cutoff_dd,
                    gmx_bool bInterCGBondeds);


/* In domdec_box.c */

/*! \brief Set the box and PBC data in \p ddbox */
void set_ddbox(gmx_domdec_t *dd, gmx_bool bMasterState, t_commrec *cr_sum,
               const t_inputrec *ir, const matrix box,
               gmx_bool bCalcUnboundedSize, const t_block *cgs, const rvec *x,
               gmx_ddbox_t *ddbox);

/*! \brief Set the box and PBC data in \p ddbox */
void set_ddbox_cr(t_commrec *cr, const ivec *dd_nc,
                  const t_inputrec *ir, const matrix box,
                  const t_block *cgs, const rvec *x,
                  gmx_ddbox_t *ddbox);

#endif
