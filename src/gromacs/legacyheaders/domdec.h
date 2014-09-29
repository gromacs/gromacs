/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2005,2006,2007,2008,2009,2010,2012,2013,2014, by the GROMACS development team, led by
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

#ifndef _domdec_h
#define _domdec_h

#include "typedefs.h"
#include "vsite.h"
#include "genborn.h"

#ifdef __cplusplus
extern "C" {
#endif

int ddglatnr(gmx_domdec_t *dd, int i);
/* Returns the global topology atom number belonging to local atom index i.
 * This function is intended for writing ascii output
 * and returns atom numbers starting at 1.
 * When dd=NULL returns i+1.
 */

t_block *dd_charge_groups_global(gmx_domdec_t *dd);
/* Return a block struct for the charge groups of the whole system */

gmx_bool dd_filled_nsgrid_home(gmx_domdec_t *dd);
/* Is the ns grid already filled with the home particles? */

void dd_store_state(gmx_domdec_t *dd, t_state *state);
/* Store the global cg indices of the home cgs in state,
 * so it can be reset, even after a new DD partitioning.
 */

gmx_domdec_zones_t *domdec_zones(gmx_domdec_t *dd);

void dd_get_ns_ranges(gmx_domdec_t *dd, int icg,
                      int *jcg0, int *jcg1, ivec shift0, ivec shift1);

int dd_natoms_vsite(gmx_domdec_t *dd);

void dd_get_constraint_range(gmx_domdec_t *dd,
                             int *at_start, int *at_end);

real dd_cutoff_mbody(gmx_domdec_t *dd);

real dd_cutoff_twobody(gmx_domdec_t *dd);

void get_pme_nnodes(const gmx_domdec_t *dd,
                    int *npmenodes_x, int *npmenodes_y);
/* Get the number of PME nodes along x and y, can be called with dd=NULL */

gmx_bool gmx_pmeonlynode(t_commrec *cr, int nodeid);
/* Return if nodeid in cr->mpi_comm_mysim is a PME-only node */

void get_pme_ddnodes(t_commrec *cr, int pmenodeid,
                     int *nmy_ddnodes, int **my_ddnodes, int *node_peer);
/* Returns the set of DD nodes that communicate with pme node cr->nodeid */

int dd_pme_maxshift_x(gmx_domdec_t *dd);
/* Returns the maximum shift for coordinate communication in PME, dim x */

int dd_pme_maxshift_y(gmx_domdec_t *dd);
/* Returns the maximum shift for coordinate communication in PME, dim y */

void make_dd_communicators(FILE *fplog, t_commrec *cr, int dd_node_order);

gmx_domdec_t *
init_domain_decomposition(FILE *fplog,
                          t_commrec *cr,
                          unsigned long Flags,
                          ivec nc,
                          real comm_distance_min, real rconstr,
                          const char *dlb_opt, real dlb_scale,
                          const char *sizex, const char *sizey, const char *sizez,
                          gmx_mtop_t *mtop, t_inputrec *ir,
                          matrix box, rvec *x,
                          gmx_ddbox_t *ddbox,
                          int *npme_x, int *npme_y);

void dd_init_bondeds(FILE *fplog,
                     gmx_domdec_t *dd, gmx_mtop_t *mtop,
                     gmx_vsite_t *vsite,
                     t_inputrec *ir, gmx_bool bBCheck, cginfo_mb_t *cginfo_mb);
/* Initialize data structures for bonded interactions */

gmx_bool dd_bonded_molpbc(gmx_domdec_t *dd, int ePBC);
/* Returns if we need to do pbc for calculating bonded interactions */

void set_dd_parameters(FILE *fplog, gmx_domdec_t *dd, real dlb_scale,
                       t_inputrec *ir,
                       gmx_ddbox_t *ddbox);
/* Set DD grid dimensions and limits,
 * should be called after calling dd_init_bondeds.
 */

gmx_bool change_dd_cutoff(t_commrec *cr, t_state *state, t_inputrec *ir,
                          real cutoff_req );
/* Change the DD non-bonded communication cut-off.
 * This could fail when trying to increase the cut-off,
 * then FALSE will be returned and the cut-off is not modified.
 */

void change_dd_dlb_cutoff_limit(t_commrec *cr);
/* Domain boundary changes due to the DD dynamic load balancing can limit
 * the cut-off distance that can be set in change_dd_cutoff. This function
 * limits the DLB such that using the currently set cut-off should still be
 * possible after subsequently setting a shorter cut-off with change_dd_cutoff.
 */

gmx_bool dd_dlb_is_locked(const gmx_domdec_t *dd);
/* Return if the DLB lock is set */

void dd_dlb_set_lock(gmx_domdec_t *dd, gmx_bool bValue);
/* Set a lock such that with DLB=auto DLB can (not) get turned on */

void dd_setup_dlb_resource_sharing(t_commrec           *cr,
                                   const gmx_hw_info_t *hwinfo,
                                   const gmx_hw_opt_t  *hw_opt);
/* When domains (PP MPI ranks) share a GPU, the individual GPU wait times
 * are meaningless, as it depends on the order in which tasks on the same
 * GPU finish. Therefore there wait times need to be averaged over the ranks
 * sharing the same GPU. This function sets up the communication for that.
 */

void setup_dd_grid(FILE *fplog, gmx_domdec_t *dd);

void dd_collect_vec(gmx_domdec_t *dd,
                    t_state *state_local, rvec *lv, rvec *v);

void dd_collect_state(gmx_domdec_t *dd,
                      t_state *state_local, t_state *state);

enum {
    ddCyclStep, ddCyclPPduringPME, ddCyclF, ddCyclWaitGPU, ddCyclPME, ddCyclNr
};

void dd_cycles_add(gmx_domdec_t *dd, float cycles, int ddCycl);
/* Add the wallcycle count to the DD counter */

void dd_force_flop_start(gmx_domdec_t *dd, t_nrnb *nrnb);
/* Start the force flop count */

void dd_force_flop_stop(gmx_domdec_t *dd, t_nrnb *nrnb);
/* Stop the force flop count */

float dd_pme_f_ratio(gmx_domdec_t *dd);
/* Return the PME/PP force load ratio, or -1 if nothing was measured.
 * Should only be called on the DD master node.
 */

void dd_move_x(gmx_domdec_t *dd, matrix box, rvec x[]);
/* Communicate the coordinates to the neighboring cells and do pbc. */

void dd_move_f(gmx_domdec_t *dd, rvec f[], rvec *fshift);
/* Sum the forces over the neighboring cells.
 * When fshift!=NULL the shift forces are updated to obtain
 * the correct virial from the single sum including f.
 */

void dd_atom_spread_real(gmx_domdec_t *dd, real v[]);
/* Communicate a real for each atom to the neighboring cells. */

void dd_atom_sum_real(gmx_domdec_t *dd, real v[]);
/* Sum the contributions to a real for each atom over the neighboring cells. */

void dd_partition_system(FILE                *fplog,
                         gmx_int64_t          step,
                         t_commrec           *cr,
                         gmx_bool             bMasterState,
                         int                  nstglobalcomm,
                         t_state             *state_global,
                         gmx_mtop_t          *top_global,
                         t_inputrec          *ir,
                         t_state             *state_local,
                         rvec               **f,
                         t_mdatoms           *mdatoms,
                         gmx_localtop_t      *top_local,
                         t_forcerec          *fr,
                         gmx_vsite_t         *vsite,
                         gmx_shellfc_t        shellfc,
                         gmx_constr_t         constr,
                         t_nrnb              *nrnb,
                         gmx_wallcycle_t      wcycle,
                         gmx_bool             bVerbose);
/* Partition the system over the nodes.
 * step is only used for printing error messages.
 * If bMasterState==TRUE then state_global from the master node is used,
 * else state_local is redistributed between the nodes.
 * When f!=NULL, *f will be reallocated to the size of state_local.
 */

void reset_dd_statistics_counters(gmx_domdec_t *dd);
/* Reset all the statistics and counters for total run counting */

void print_dd_statistics(t_commrec *cr, t_inputrec *ir, FILE *fplog);

/* In domdec_con.c */

void dd_move_f_vsites(gmx_domdec_t *dd, rvec *f, rvec *fshift);

void dd_clear_f_vsites(gmx_domdec_t *dd, rvec *f);

void dd_move_x_constraints(gmx_domdec_t *dd, matrix box,
                           rvec *x0, rvec *x1, gmx_bool bX1IsCoord);
/* Move x0 and also x1 if x1!=NULL. bX1IsCoord tells if to do PBC on x1 */

void dd_move_x_vsites(gmx_domdec_t *dd, matrix box, rvec *x);

int *dd_constraints_nlocalatoms(gmx_domdec_t *dd);

void dd_clear_local_constraint_indices(gmx_domdec_t *dd);

void dd_clear_local_vsite_indices(gmx_domdec_t *dd);

int dd_make_local_vsites(gmx_domdec_t *dd, int at_start, t_ilist *lil);

int dd_make_local_constraints(gmx_domdec_t *dd, int at_start,
                              const gmx_mtop_t *mtop,
                              const int *cginfo,
                              gmx_constr_t constr, int nrec,
                              t_ilist *il_local);

void init_domdec_constraints(gmx_domdec_t *dd,
                             gmx_mtop_t   *mtop);

void init_domdec_vsites(gmx_domdec_t *dd, int n_intercg_vsite);


/* In domdec_top.c */

void dd_print_missing_interactions(FILE *fplog, t_commrec *cr,
                                   int local_count,  gmx_mtop_t *top_global, t_state *state_local);

void dd_make_reverse_top(FILE *fplog,
                         gmx_domdec_t *dd, gmx_mtop_t *mtop,
                         gmx_vsite_t *vsite,
                         t_inputrec *ir, gmx_bool bBCheck);

void dd_make_local_cgs(gmx_domdec_t *dd, t_block *lcgs);

void dd_make_local_top(gmx_domdec_t *dd, gmx_domdec_zones_t *zones,
                       int npbcdim, matrix box,
                       rvec cellsize_min, ivec npulse,
                       t_forcerec *fr,
                       rvec *cgcm_or_x,
                       gmx_vsite_t *vsite,
                       gmx_mtop_t *top, gmx_localtop_t *ltop);

void dd_sort_local_top(gmx_domdec_t *dd, t_mdatoms *mdatoms,
                       gmx_localtop_t *ltop);
/* Sort ltop->ilist when we are doing free energy. */

gmx_localtop_t *dd_init_local_top(gmx_mtop_t *top_global);

void dd_init_local_state(gmx_domdec_t *dd,
                         t_state *state_global, t_state *local_state);

t_blocka *make_charge_group_links(gmx_mtop_t *mtop, gmx_domdec_t *dd,
                                  cginfo_mb_t *cginfo_mb);

void dd_bonded_cg_distance(FILE *fplog, gmx_mtop_t *mtop,
                           t_inputrec *ir, rvec *x, matrix box,
                           gmx_bool bBCheck,
                           real *r_2b, real *r_mb);

void write_dd_pdb(const char *fn, gmx_int64_t step, const char *title,
                  gmx_mtop_t *mtop,
                  t_commrec *cr,
                  int natoms, rvec x[], matrix box);
/* Dump a pdb file with the current DD home + communicated atoms.
 * When natoms=-1, dump all known atoms.
 */


/* In domdec_setup.c */

real comm_box_frac(ivec dd_nc, real cutoff, gmx_ddbox_t *ddbox);
/* Returns the volume fraction of the system that is communicated */

real dd_choose_grid(FILE *fplog,
                    t_commrec *cr, gmx_domdec_t *dd, t_inputrec *ir,
                    gmx_mtop_t *mtop, matrix box, gmx_ddbox_t *ddbox,
                    gmx_bool bDynLoadBal, real dlb_scale,
                    real cellsize_limit, real cutoff_dd,
                    gmx_bool bInterCGBondeds);
/* Determines the optimal DD cell setup dd->nc and possibly npmenodes
 * for the system.
 * On the master node returns the actual cellsize limit used.
 */


/* In domdec_box.c */

void set_ddbox(gmx_domdec_t *dd, gmx_bool bMasterState, t_commrec *cr_sum,
               t_inputrec *ir, matrix box,
               gmx_bool bCalcUnboundedSize, t_block *cgs, rvec *x,
               gmx_ddbox_t *ddbox);

void set_ddbox_cr(t_commrec *cr, ivec *dd_nc,
                  t_inputrec *ir, matrix box, t_block *cgs, rvec *x,
                  gmx_ddbox_t *ddbox);

#ifdef __cplusplus
}
#endif

#endif  /* _domdec_h */
