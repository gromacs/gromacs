#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#ifndef _domdec_h
#define _domdec_h

#include "typedefs.h"

#ifdef GMX_MPI
#include <mpi.h>
#endif

extern int dd_nicg(gmx_domdec_t *dd);

extern int dd_ncg_tot(gmx_domdec_t *dd);

extern void dd_get_ns_ranges(gmx_domdec_t *dd,int icg,
			     int *jcg0,int *jcg1,ivec shift0,ivec shift1);

extern gmx_domdec_t *init_domain_decomposition(FILE *fplog,t_commrec *cr,
					       ivec nc,int ncg,int natoms,
					       matrix box);

extern void setup_dd_grid(FILE *fplog,matrix box,gmx_domdec_t *dd);

extern void get_cg_distribution(FILE *fplog,gmx_domdec_t *dd,
				matrix box,t_block *cgs,rvec pos[]);

extern void setup_dd_communication(FILE *fplog,gmx_domdec_t *dd,t_block *cgs,
				   rvec cg_cm[],matrix box,real r_comm);

extern void dd_collect_vec(gmx_domdec_t *dd,t_block *cgs,rvec *lv,rvec *v);

extern void dd_collect_state(gmx_domdec_t *dd,t_block *cgs,
			     t_state *state_local,t_state *state);

extern void dd_distribute_state(gmx_domdec_t *dd,t_block *cgs,
				t_state *state,t_state *state_local);

extern void dd_redistribute_cg(FILE *fplog,
			       gmx_domdec_t *dd,t_block *cgs,
			       t_state *state,rvec cg_cm[]);

extern void dd_calc_cgcm_home(gmx_domdec_t *dd,t_block *cgs,
			      rvec x[],rvec cg_cm[]);

extern void dd_move_x(gmx_domdec_t *dd,t_block *cgs,rvec x[],rvec buf[]);

extern void dd_move_f(gmx_domdec_t *dd,t_block *cgs,rvec f[],rvec buf[]);

extern void make_local_top(FILE *fplog,gmx_domdec_t *dd,
			   t_topology *top,t_topology *ltop,
			   t_nsborder *nsb);

#endif	/* _domdec_h */
