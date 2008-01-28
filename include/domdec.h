#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#ifndef _domdec_h
#define _domdec_h

#include "typedefs.h"
#include "vsite.h"

#ifdef GMX_MPI
#include <mpi.h>
#endif

extern int glatnr(gmx_domdec_t *dd,int i);
/* Returns the global topology atom number belonging to local atom index i.
 * This function is intended for writing ascii output
 * and returns atom numbers starting at 1.
 * When dd=NULL returns i+1.
 */

extern bool dd_filled_nsgrid_home(gmx_domdec_t *dd);
/* Is the ns grid already filled with the home particles? */

extern void dd_store_state(gmx_domdec_t *dd,t_state *state);
/* Store the global cg indices of the home cgs in state,
 * so it can be reset, even after a new DD partitioning.
 */

extern void dd_get_ns_ranges(gmx_domdec_t *dd,int icg,
			     int *jcg0,int *jcg1,ivec shift0,ivec shift1);

extern int dd_natoms_vsite(gmx_domdec_t *dd);

extern void dd_get_constraint_range(gmx_domdec_t *dd,
				    int *at_start,int *at_end);

extern void dd_make_reverse_top(FILE *fplog,
				gmx_domdec_t *dd,t_topology *top,
				gmx_vsite_t *vsite,gmx_constr_t constr,
				bool bDynamics,int eeltype);

extern int gmx_ddcoord2pmeslab(t_commrec *cr,int x,int y,int z);
/* Returns the pme slab for DD cell x,y,z */

extern bool gmx_pmeonlynode(t_commrec *cr,int nodeid);
/* Return if nodeid in cr->mpi_comm_mysim is a PME-only node */

extern void get_pme_ddnodes(t_commrec *cr,int pmenodeid,
			    int *nmy_ddnodes,int **my_ddnodes,int *node_peer);
/* Returns the set of DD nodes that communicate with pme node cr->nodeid */

extern int dd_pme_maxshift(gmx_domdec_t *dd);
/* Returns the maximum shift for coordinate communication in PME */

extern void make_dd_communicators(FILE *fplog,t_commrec *cr,int dd_node_order);

extern gmx_domdec_t *init_domain_decomposition(FILE *fplog,
					       t_commrec *cr,ivec nc,
					       real comm_distance_min,
					       real rconstr,
					       bool bDynLoadBal,
					       char *sizex,
					       char *sizey,
					       char *sizez,
					       t_topology *top,matrix box,
					       t_inputrec *ir);

extern void set_dd_parameters(FILE *fplog,gmx_domdec_t *dd,
			      t_topology *top,t_inputrec *ir,t_forcerec *fr,
			      matrix box);

extern void setup_dd_grid(FILE *fplog,gmx_domdec_t *dd);

extern void dd_collect_vec(gmx_domdec_t *dd,t_block *cgs,rvec *lv,rvec *v);

extern void dd_collect_state(gmx_domdec_t *dd,t_block *cgs,
			     t_state *state_local,t_state *state);

enum { ddCyclStep, ddCyclPPduringPME, ddCyclF, ddCyclPME, ddCyclNr };

extern void dd_cycles_add(gmx_domdec_t *dd,float cycles,int ddCycl);
/* Add the wallcycle count to the DD counter */

extern void dd_force_flop_start(gmx_domdec_t *dd,t_nrnb *nrnb);
/* Start the force flop count */

extern void dd_force_flop_stop(gmx_domdec_t *dd,t_nrnb *nrnb);
/* Stop the force flop count */

enum {
  ddForward,ddBackward
};

extern void dd_sendrecv_int(const gmx_domdec_t *dd,
			    int ddim,int direction,
			    int *buf_s,int n_s,
			    int *buf_r,int n_r);

extern void dd_sendrecv_rvec(const gmx_domdec_t *dd,
			     int ddim,int direction,
			     rvec *buf_s,int n_s,
			     rvec *buf_r,int n_r);
/* Move data (int/rvec) one cell in decomposition dimension ddim
 * over the DD grid.
 * For direction see the enum above.
 */

extern void dd_sendrecv2_rvec(const gmx_domdec_t *dd,
			      int ddim,
			      rvec *buf_s_fw,int n_s_fw,
			      rvec *buf_r_fw,int n_r_fw,
			      rvec *buf_s_bw,int n_s_bw,
			      rvec *buf_r_bw,int n_r_bw);
/* Move data (rvec) forward and backward 
 * one cell in decomposition dimension ddim over the DD grid.
 * If possible, this is done simultaneously.
 */

extern void dd_move_x(gmx_domdec_t *dd,matrix box,rvec x[],rvec buf[]);
/* buf should should have size natoms (of the whole system)
 * although in most cases far less will be used.
 */

extern void dd_move_f(gmx_domdec_t *dd,rvec f[],rvec buf[],rvec *fshift);
/* buf should should have size natoms (of the whole system)
 * although in most cases far less will be used.
 * When fshift!=NULL the shift forces are updated to obtain
 * the correct virial from the single sum including f.
 */

extern void dd_bcast(gmx_domdec_t *dd,int nbytes,void *data);
/* Broadcasts nbytes bytes from the DD master nodes to all DD nodes */

extern void dd_partition_system(FILE            *fplog,
				int             step,
				t_commrec       *cr,
                                bool            bMasterState,
				t_state         *state_global,
				t_topology      *top_global,
				t_inputrec      *ir,
				t_state         *state_local,
				rvec            **f,
				rvec            **buf,
				t_mdatoms       *mdatoms,
				t_topology      *top_local,
				t_forcerec      *fr,
				gmx_vsite_t     *vsite,
				gmx_shellfc_t   shellfc,
				gmx_constr_t    constr,
				t_nrnb          *nrnb,
				gmx_wallcycle_t wcycle,
				bool            bVerbose);
/* Partition the system over the nodes.
 * step is only used for printing error messages.
 * If bMasterState==TRUE then state_global from the master node is used,
 * else state_local is redistributed between the nodes.
 */

void print_dd_statistics(t_commrec *cr,t_inputrec *ir,FILE *fplog);

/* In domdec_con.c */

extern void dd_move_f_vsites(gmx_domdec_t *dd,rvec *f,rvec *fshift);

extern void dd_clear_f_vsites(gmx_domdec_t *dd,rvec *f);

extern void dd_move_x_constraints(gmx_domdec_t *dd,matrix box,
				  rvec *x0,rvec *x1);
/* Move x0 and also x1 if x1!=NULL */

extern void dd_move_x_vsites(gmx_domdec_t *dd,matrix box,rvec *x);

extern void dd_clear_local_constraint_indices(gmx_domdec_t *dd);

extern void dd_clear_local_vsite_indices(gmx_domdec_t *dd);

extern int dd_make_local_vsites(gmx_domdec_t *dd,int at_start,t_ilist *lil);

extern int dd_make_local_constraints(gmx_domdec_t *dd,int at_start,t_iatom *ia,
				     gmx_constr_t constr,int nrec);

extern void init_domdec_constraints(gmx_domdec_t *dd,
				    int natoms,t_idef *idef,
				    gmx_constr_t constr);

extern void init_domdec_vsites(gmx_domdec_t *dd,int natoms);


/* In domdec_top.c */

extern void dd_print_missing_interactions(FILE *fplog,t_commrec *cr,
					  int local_count);

extern void dd_make_local_cgs(gmx_domdec_t *dd,t_block *lcgs);

extern void dd_make_local_top(FILE *fplog,gmx_domdec_t *dd,
			      matrix box,real rc,ivec npulse,
			      t_forcerec *fr,gmx_vsite_t *vsite,
			      t_topology *top,t_topology *ltop);

extern t_topology *dd_init_local_top(t_topology *top_global);

extern t_state *dd_init_local_state(gmx_domdec_t *dd,t_state *state_global);

#endif	/* _domdec_h */


