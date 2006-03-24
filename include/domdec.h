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
					       matrix box,
					       t_idef *idef);

extern void setup_dd_grid(FILE *fplog,matrix box,gmx_domdec_t *dd);

extern void dd_collect_vec(gmx_domdec_t *dd,t_block *cgs,rvec *lv,rvec *v);

extern void dd_collect_state(gmx_domdec_t *dd,t_block *cgs,
			     t_state *state_local,t_state *state);

extern void dd_move_x(gmx_domdec_t *dd,rvec x[],rvec buf[]);

extern void dd_start_move_x(gmx_domdec_t *dd,rvec x[],rvec buf[]);

extern void dd_finish_move_x(gmx_domdec_t *dd);
/* buf should should have size natoms (of the whole system)
 * although in most cases far less will be used.
 *
 * dd_move_x waits with returning until the communication has finished.
 *
 * With dd_start_move_x, a call to dd_finish_move_x is required
 * to make sure the communication has finished.
 * Before that only the home coordinates in x can be used
 * and buf can not be used.
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

#endif	/* _domdec_h */
