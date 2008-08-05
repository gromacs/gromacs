#ifdef HAVE_CONFIG_H
#include<config.h>
#endif

#include "typedefs.h"
 
/* Initialization function, also predicts the initial shell postions.
 * If x!=NULL, the shells are predict for the global coordinates x.
 */
extern gmx_shellfc_t init_shell_flexcon(FILE *log,
					gmx_mtop_t *mtop,int nflexcon,
					rvec *x);

/* Get the local shell with domain decomposition */
extern void make_local_shells(t_commrec *cr,t_mdatoms *md,
			      gmx_shellfc_t shfc);

/* Optimize shell positions */
extern int relax_shell_flexcon(FILE *log,t_commrec *cr,bool bVerbose,
			       int mdstep,t_inputrec *inputrec,
			       bool bDoNS,bool bStopCM,
			       t_topology *top,gmx_constr_t constr,
			       real ener[],t_fcdata *fcd,
			       t_state *state,rvec f[],
			       rvec buf[],tensor force_vir,
			       t_mdatoms *md,
			       t_nrnb *nrnb,gmx_wallcycle_t wcycle,
			       t_graph *graph,
			       gmx_groups_t *groups,t_groups *grps,
			       gmx_shellfc_t shfc,
			       t_forcerec *fr,
			       real t,rvec mu_tot,
			       int natoms,bool *bConverged,
			       gmx_vsite_t *vsite,
			       FILE *fp_field);
