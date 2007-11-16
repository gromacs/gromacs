#ifdef HAVE_CONFIG_H
#include<config.h>
#endif

#include "typedefs.h"
 
/* Initialization function, also predicts the initial shell postions.
 * With domain decomposition the prediction is for the global x.
 */
extern gmx_shellfc_t init_shell_flexcon(FILE *log,t_commrec *cr,
					t_topology *top,int nflexcon,
					bool bContinuation,rvec *x);

/* Get the local shell with domain decomposition */
extern void make_local_shells(gmx_domdec_t *dd,t_mdatoms *md,
			      gmx_shellfc_t shfc);

/* Optimize shell positions */
extern int relax_shell_flexcon(FILE *log,t_commrec *cr,bool bVerbose,
			       int mdstep,t_inputrec *inputrec,
			       bool bDoNS,bool bStopCM,
			       t_topology *top,gmx_constr_t constr,
			       real ener[],t_fcdata *fcd,
			       t_state *state,rvec f[],
			       rvec buf[],t_mdatoms *md,
			       t_nrnb *nrnb,gmx_wallcycle_t wcycle,
			       t_graph *graph,t_groups *grps,
			       gmx_shellfc_t shfc,
			       t_forcerec *fr,
			       real t,rvec mu_tot,
			       int natoms,bool *bConverged,
			       gmx_vsite_t *vsite,
			       FILE *fp_field);
