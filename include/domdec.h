#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#ifndef _domdec_h
#define _domdec_h

#include "typedefs.h"

#ifdef GMX_MPI
#include <mpi.h>
#endif

extern int glatnr(gmx_domdec_t *dd,int i);
/* Returns the global topology atom number belonging to local atom index i.
 * This function is intended for writing ascii output
 * and returns atom numbers starting at 1.
 * When dd=NULL returns i+1.
 */

extern int dd_nicg(gmx_domdec_t *dd);

extern int dd_ncg_tot(gmx_domdec_t *dd);

extern void dd_get_ns_ranges(gmx_domdec_t *dd,int icg,
			     int *jcg0,int *jcg1,ivec shift0,ivec shift1);

extern void dd_make_reverse_top(gmx_domdec_t *dd,
				int natoms,t_idef *idef,bool bDynamics);

extern bool dd_node2pme_cart_coords(t_commrec *cr,int nodeid,int *coords);
/* Returns TRUE if nodeid is a PP node.
 * Set coords to the Cartesian coordinates of the PME-only node
 * that the PP node nodeid communicates with.
 */

extern int dd_node2pmenode(t_commrec *cr,int nodeid);
/* Returns the nodeid in cr->mpi_comm_mysim of the PME-only node
 * that nodeid communicates with.
 * Returns -1 if nodeid is a PME-only node.
 */

extern void make_dd_communicators(FILE *fplog,t_commrec *cr);

extern gmx_domdec_t *init_domain_decomposition(FILE *fplog,
					       t_commrec *cr,ivec nc);

extern void setup_dd_grid(FILE *fplog,matrix box,gmx_domdec_t *dd);

extern void dd_collect_vec(gmx_domdec_t *dd,t_block *cgs,rvec *lv,rvec *v);

extern void dd_collect_state(gmx_domdec_t *dd,t_block *cgs,
			     t_state *state_local,t_state *state);

enum {
  ddForward,ddBackward
};

extern void dd_sendrecv_int(const gmx_domdec_t *dd,
			    int dim,int direction,
			    int *buf_s,int n_s,
			    int *buf_r,int n_r);

extern void dd_sendrecv_rvec(const gmx_domdec_t *dd,
			     int dim,int direction,
			     rvec *buf_s,int n_s,
			     rvec *buf_r,int n_r);
/* Move data (int/rvec) one cell in dimension dim over the DD grid.
 * For direction see the enum above.
 */

extern void dd_move_x(gmx_domdec_t *dd,rvec x[],rvec buf[]);
/* buf should should have size natoms (of the whole system)
 * although in most cases far less will be used.
 */

extern void dd_move_f(gmx_domdec_t *dd,rvec f[],rvec buf[]);
/* buf should should have size natoms (of the whole system)
 * although in most cases far less will be used.
 */

extern void dd_partition_system(FILE         *fplog,
				gmx_domdec_t *dd,
				bool         bMasterState,
				t_state      *state_global,
				t_topology   *top_global,
				t_inputrec   *ir,
				t_state      *state_local,
				rvec         *buf,
				t_mdatoms    *mdatoms,
				t_topology   *top_local,
				t_nsborder   *nsb,
				t_forcerec   *fr,
				t_nrnb       *nrnb);
/* Partition the system over the nodes.
 * If bMasterState==TRUE then state_global from the master node is used,
 * else state_local is redistributed between the nodes.
 */

/* In domdec_con.c */

extern void dd_move_x_constraints(gmx_domdec_t *dd,rvec *x);

extern void clear_local_constraint_indices(gmx_domdec_t *dd);

extern void make_local_constraints(gmx_domdec_t *dd,t_iatom *ia,int nrec,
				   gmx_domdec_constraints_t *dc);

extern gmx_domdec_constraints_t *init_domdec_constraints(int natoms,
							 t_idef *idef,
							 bool bDynamics);

#endif	/* _domdec_h */
