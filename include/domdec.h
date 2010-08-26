/* -*- mode: c; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; c-file-style: "stroustrup"; -*-
 *
 * 
 * This file is part of Gromacs        Copyright (c) 1991-2008
 * David van der Spoel, Erik Lindahl, Berk Hess, University of Groningen.
 *
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 *
 * To help us fund GROMACS development, we humbly ask that you cite
 * the research papers on the package. Check out http://www.gromacs.org
 * 
 * And Hey:
 * Gnomes, ROck Monsters And Chili Sauce
 */

#ifndef _domdec_h
#define _domdec_h

#include "typedefs.h"
#include "vsite.h"
#include "genborn.h"

#ifdef GMX_LIB_MPI
#include <mpi.h>
#endif
#ifdef GMX_THREADS
#include "tmpi.h"
#endif

#ifdef __cplusplus
extern "C" {
#endif

int ddglatnr(gmx_domdec_t *dd,int i);
/* Returns the global topology atom number belonging to local atom index i.
 * This function is intended for writing ascii output
 * and returns atom numbers starting at 1.
 * When dd=NULL returns i+1.
 */

t_block *dd_charge_groups_global(gmx_domdec_t *dd);
/* Return a block struct for the charge groups of the whole system */

gmx_bool dd_filled_nsgrid_home(gmx_domdec_t *dd);
/* Is the ns grid already filled with the home particles? */

void dd_store_state(gmx_domdec_t *dd,t_state *state);
/* Store the global cg indices of the home cgs in state,
 * so it can be reset, even after a new DD partitioning.
 */

gmx_domdec_zones_t *domdec_zones(gmx_domdec_t *dd);

void dd_get_ns_ranges(gmx_domdec_t *dd,int icg,
                             int *jcg0,int *jcg1,ivec shift0,ivec shift1);

int dd_natoms_vsite(gmx_domdec_t *dd);

void dd_get_constraint_range(gmx_domdec_t *dd,
				    int *at_start,int *at_end);

real dd_cutoff_mbody(gmx_domdec_t *dd);

real dd_cutoff_twobody(gmx_domdec_t *dd);

gmx_bool gmx_pmeonlynode(t_commrec *cr,int nodeid);
/* Return if nodeid in cr->mpi_comm_mysim is a PME-only node */

void get_pme_ddnodes(t_commrec *cr,int pmenodeid,
			    int *nmy_ddnodes,int **my_ddnodes,int *node_peer);
/* Returns the set of DD nodes that communicate with pme node cr->nodeid */

int dd_pme_maxshift_x(gmx_domdec_t *dd);
/* Returns the maximum shift for coordinate communication in PME, dim x */

int dd_pme_maxshift_y(gmx_domdec_t *dd);
/* Returns the maximum shift for coordinate communication in PME, dim y */

void make_dd_communicators(FILE *fplog,t_commrec *cr,int dd_node_order);

gmx_domdec_t *
init_domain_decomposition(FILE *fplog,
                          t_commrec *cr,
                          unsigned long Flags,
                          ivec nc,
                          real comm_distance_min,real rconstr,
                          const char *dlb_opt,real dlb_scale,
                          const char *sizex,const char *sizey,const char *sizez,
                          gmx_mtop_t *mtop,t_inputrec *ir,
                          matrix box,rvec *x,
                          gmx_ddbox_t *ddbox,
                          int *npme_x, int *npme_y);

void dd_init_bondeds(FILE *fplog,
                            gmx_domdec_t *dd,gmx_mtop_t *mtop,
                            gmx_vsite_t *vsite,gmx_constr_t constr,
                            t_inputrec *ir,gmx_bool bBCheck,cginfo_mb_t *cginfo_mb);
/* Initialize data structures for bonded interactions */

void set_dd_parameters(FILE *fplog,gmx_domdec_t *dd,real dlb_scale,
                              t_inputrec *ir,t_forcerec *fr,
                              gmx_ddbox_t *ddbox);
/* Set DD grid dimensions and limits,
 * should be called after calling dd_init_bondeds.
 */

void setup_dd_grid(FILE *fplog,gmx_domdec_t *dd);

void dd_collect_vec(gmx_domdec_t *dd,
                           t_state *state_local,rvec *lv,rvec *v);

void dd_collect_state(gmx_domdec_t *dd,
                             t_state *state_local,t_state *state);

enum { ddCyclStep, ddCyclPPduringPME, ddCyclF, ddCyclPME, ddCyclNr };

void dd_cycles_add(gmx_domdec_t *dd,float cycles,int ddCycl);
/* Add the wallcycle count to the DD counter */

void dd_force_flop_start(gmx_domdec_t *dd,t_nrnb *nrnb);
/* Start the force flop count */

void dd_force_flop_stop(gmx_domdec_t *dd,t_nrnb *nrnb);
/* Stop the force flop count */

void dd_move_x(gmx_domdec_t *dd,matrix box,rvec x[]);
/* Communicate the coordinates to the neighboring cells and do pbc. */

void dd_move_f(gmx_domdec_t *dd,rvec f[],rvec *fshift);
/* Sum the forces over the neighboring cells.
 * When fshift!=NULL the shift forces are updated to obtain
 * the correct virial from the single sum including f.
 */

void dd_atom_spread_real(gmx_domdec_t *dd,real v[]);
/* Communicate a real for each atom to the neighboring cells. */

void dd_atom_sum_real(gmx_domdec_t *dd,real v[]);
/* Sum the contributions to a real for each atom over the neighboring cells. */

void dd_partition_system(FILE            *fplog,
                                gmx_large_int_t      step,
                                t_commrec       *cr,
                                gmx_bool            bMasterState,
                                int             nstglobalcomm,
                                t_state         *state_global,
                                gmx_mtop_t      *top_global,
                                t_inputrec      *ir,
                                t_state         *state_local,
                                rvec            **f,
                                t_mdatoms       *mdatoms,
                                gmx_localtop_t  *top_local,
                                t_forcerec      *fr,
                                gmx_vsite_t     *vsite,
                                gmx_shellfc_t   shellfc,
                                gmx_constr_t    constr,
                                t_nrnb          *nrnb,
                                gmx_wallcycle_t wcycle,
                                gmx_bool            bVerbose);
/* Partition the system over the nodes.
 * step is only used for printing error messages.
 * If bMasterState==TRUE then state_global from the master node is used,
 * else state_local is redistributed between the nodes.
 * When f!=NULL, *f will be reallocated to the size of state_local.
 */

void reset_dd_statistics_counters(gmx_domdec_t *dd);
/* Reset all the statistics and counters for total run counting */

void print_dd_statistics(t_commrec *cr,t_inputrec *ir,FILE *fplog);

/* In domdec_con.c */

void dd_move_f_vsites(gmx_domdec_t *dd,rvec *f,rvec *fshift);

void dd_clear_f_vsites(gmx_domdec_t *dd,rvec *f);

void dd_move_x_constraints(gmx_domdec_t *dd,matrix box,
				  rvec *x0,rvec *x1);
/* Move x0 and also x1 if x1!=NULL */

void dd_move_x_vsites(gmx_domdec_t *dd,matrix box,rvec *x);

int *dd_constraints_nlocalatoms(gmx_domdec_t *dd);

void dd_clear_local_constraint_indices(gmx_domdec_t *dd);

void dd_clear_local_vsite_indices(gmx_domdec_t *dd);

int dd_make_local_vsites(gmx_domdec_t *dd,int at_start,t_ilist *lil);

int dd_make_local_constraints(gmx_domdec_t *dd,int at_start,
                                     gmx_mtop_t *mtop,
                                     gmx_constr_t constr,int nrec,
                                     t_ilist *il_local);

void init_domdec_constraints(gmx_domdec_t *dd,
                                    int natoms,gmx_mtop_t *mtop,
                                    gmx_constr_t constr);

void init_domdec_vsites(gmx_domdec_t *dd,int natoms);


/* In domdec_top.c */

void dd_print_missing_interactions(FILE *fplog,t_commrec *cr,
                                          int local_count,  gmx_mtop_t *top_global, t_state *state_local);

void dd_make_reverse_top(FILE *fplog,
                                gmx_domdec_t *dd,gmx_mtop_t *mtop,
                                gmx_vsite_t *vsite,gmx_constr_t constr,
                                t_inputrec *ir,gmx_bool bBCheck);

void dd_make_local_cgs(gmx_domdec_t *dd,t_block *lcgs);

void dd_make_local_top(FILE *fplog,
                              gmx_domdec_t *dd,gmx_domdec_zones_t *zones,
                              int npbcdim,matrix box,
                              rvec cellsize_min,ivec npulse,
                              t_forcerec *fr,gmx_vsite_t *vsite,
                              gmx_mtop_t *top,gmx_localtop_t *ltop);

void dd_sort_local_top(gmx_domdec_t *dd,t_mdatoms *mdatoms,
                              gmx_localtop_t *ltop);
/* Sort ltop->ilist when we are doing free energy. */

gmx_localtop_t *dd_init_local_top(gmx_mtop_t *top_global);

void dd_init_local_state(gmx_domdec_t *dd,
                                t_state *state_global,t_state *local_state);

t_blocka *make_charge_group_links(gmx_mtop_t *mtop,gmx_domdec_t *dd,
                                         cginfo_mb_t *cginfo_mb);

void dd_bonded_cg_distance(FILE *fplog,
                                  gmx_domdec_t *dd,gmx_mtop_t *mtop,
                                  t_inputrec *ir,rvec *x,matrix box,
                                  gmx_bool bBCheck,
                                  real *r_2b,real *r_mb);

void write_dd_pdb(const char *fn,gmx_large_int_t step,const char *title,
                         gmx_mtop_t *mtop,
                         t_commrec *cr,
                         int natoms,rvec x[],matrix box);
/* Dump a pdb file with the current DD home + communicated atoms.
 * When natoms=-1, dump all known atoms.
 */


/* In domdec_setup.c */

real comm_box_frac(ivec dd_nc,real cutoff,gmx_ddbox_t *ddbox);
/* Returns the volume fraction of the system that is communicated */

real dd_choose_grid(FILE *fplog,
                           t_commrec *cr,gmx_domdec_t *dd,t_inputrec *ir,
                           gmx_mtop_t *mtop,matrix box,gmx_ddbox_t *ddbox,
                           gmx_bool bDynLoadBal,real dlb_scale,
                           real cellsize_limit,real cutoff_dd,
                           gmx_bool bInterCGBondeds,gmx_bool bInterCGMultiBody);
/* Determines the optimal DD cell setup dd->nc and possibly npmenodes
 * for the system.
 * On the master node returns the actual cellsize limit used.
 */


/* In domdec_box.c */

void set_ddbox(gmx_domdec_t *dd,gmx_bool bMasterState,t_commrec *cr_sum,
                      t_inputrec *ir,matrix box,
                      gmx_bool bCalcUnboundedSize,t_block *cgs,rvec *x,
                      gmx_ddbox_t *ddbox);

void set_ddbox_cr(t_commrec *cr,ivec *dd_nc,
                         t_inputrec *ir,matrix box,t_block *cgs,rvec *x,
                         gmx_ddbox_t *ddbox);

#ifdef __cplusplus
}
#endif

#endif	/* _domdec_h */
