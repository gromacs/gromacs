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

extern void gmx_ddindex2xyz(ivec nc,int ind,ivec xyz);
/* Returns the DD cell coordinates xyz for DD index ind,
 * nc should contain the number of cells.
 */

extern void dd_get_ns_ranges(gmx_domdec_t *dd,int icg,
			     int *jcg0,int *jcg1,ivec shift0,ivec shift1);

extern void dd_make_reverse_top(FILE *fplog,
				gmx_domdec_t *dd,t_topology *top,
				bool bDynamics,int eeltype);

extern int gmx_ddindex2pmeslab(t_commrec *cr,int ddindex);
/* Returns the pme slab for DD index ddindex */

extern int gmx_ddindex2nodeid(t_commrec *cr,int ddindex);
/* Returns the nodeid in cr->mpi_comm_mysim for ddindex */

extern bool gmx_pmeonlynode(t_commrec *cr,int nodeid);
/* Return if nodeid in cr->mpi_comm_mysim is a PME-only node */

extern void make_dd_communicators(FILE *fplog,t_commrec *cr,bool bCartesian);

extern gmx_domdec_t *init_domain_decomposition(FILE *fplog,
					       t_commrec *cr,ivec nc,
					       real comm_distance_min,
					       bool bDynLoadBal,
					       char *loadx,
					       char *loady,
					       char *loadz);

extern void set_dd_parameters(FILE *fplog,gmx_domdec_t *dd,
			      t_topology *top,t_inputrec *ir,t_forcerec *fr);

extern void setup_dd_grid(FILE *fplog,gmx_domdec_t *dd);

extern void dd_collect_vec(gmx_domdec_t *dd,t_block *cgs,rvec *lv,rvec *v);

extern void dd_collect_state(gmx_domdec_t *dd,t_block *cgs,
			     t_state *state_local,t_state *state);

enum { ddCyclMoveX, ddCyclF, ddCyclMoveF, ddCyclPME, ddCyclNr };

extern void dd_cycles_add(gmx_domdec_t *dd,float cycles,int ddCycl);
/* Add the wallcycle count to the DD counter */

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

extern void dd_partition_system(FILE         *fplog,
				int          step,
				t_commrec    *cr,
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
				t_nrnb       *nrnb,
				gmx_wallcycle_t wcycle,
				bool         bVerbose);
/* Partition the system over the nodes.
 * step is only used for printing error messages.
 * If bMasterState==TRUE then state_global from the master node is used,
 * else state_local is redistributed between the nodes.
 */


/* In domdec_con.c */

extern void dd_move_f_vsites(gmx_domdec_t *dd,rvec *f,rvec *fshift);

extern void dd_move_x_constraints(gmx_domdec_t *dd,matrix box,rvec *x);

extern void dd_move_x_vsites(gmx_domdec_t *dd,matrix box,rvec *x);

extern void clear_local_constraint_indices(gmx_domdec_t *dd);

extern void clear_local_vsite_indices(gmx_domdec_t *dd);

extern void make_local_vsites(gmx_domdec_t *dd,t_ilist *lil);

extern void make_local_constraints(gmx_domdec_t *dd,t_iatom *ia,int nrec);

extern void init_domdec_constraints(gmx_domdec_t *dd,
				    int natoms,t_idef *idef,t_block *cgs,
				    bool bDynamics);

extern void init_domdec_vsites(gmx_domdec_t *dd,int natoms);


/* In domdec_top.c */

extern void dd_print_missing_interactions(FILE *fplog,t_commrec *cr,
					  int local_count);

extern void make_local_cgs(gmx_domdec_t *dd,t_block *lcgs);

extern void make_local_top(FILE *fplog,gmx_domdec_t *dd,
			   t_forcerec *fr,t_topology *top,t_topology *ltop);

#endif	/* _domdec_h */
