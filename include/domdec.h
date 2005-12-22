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

extern int dd_cg_j0(gmx_domdec_t *dd,int icg);

extern int dd_cg_j1(gmx_domdec_t *dd,int icg);

extern gmx_domdec_t *init_domain_decomposition(FILE *fplog,t_commrec *cr,
					       ivec nc,int ncg,int natoms,
					       matrix box);

extern void setup_dd_grid(FILE *fplog,matrix box,gmx_domdec_t *dd);

extern void get_cg_distribution(FILE *fplog,gmx_domdec_t *dd,
				matrix box,t_block *cgs,rvec pos[]);

extern void dd_collect_vec(gmx_domdec_t *dd,t_block *cgs,rvec *lv,rvec *v);

extern void dd_collect_state(gmx_domdec_t *dd,t_block *cgs,
			     t_state *state_local,t_state *state);

extern void dd_distribute_state(gmx_domdec_t *dd,t_block *cgs,
				t_state *state,t_state *state_local);

extern void dd_move_x(gmx_domdec_t *dd,rvec x[]);

extern void dd_move_f(gmx_domdec_t *dd,rvec f[],rvec buf[]);

extern void make_local_top(FILE *fplog,gmx_domdec_t *dd,
			   t_topology *top,t_topology *ltop,
			   t_nsborder *nsb);

#endif	/* _domdec_h */
