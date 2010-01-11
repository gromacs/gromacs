#ifndef _qhop_h
#define _qhop_h

#ifdef HAVE_IDENT
#ident	"@(#) qhop.h 1 28/2/01"
#endif /* HAVE_IDENT */
#include "typedefs.h"
#include "pbc.h"
#include "network.h"
#include "tgroup.h"
#include "gmx_qhop_db.h"
#include "resall.h"

#define DOO   0.35
#define BOND  0.15
#define HBOND  0.2
#define HOPANG 120

t_qhoprec *mk_qhoprec(void);
extern int init_qhop(t_commrec *cr, gmx_mtop_t *mtop, t_inputrec *ir, t_forcerec *fr,
	      rvec *x,matrix box, t_mdatoms *md, qhop_db_t *db);

extern void do_qhop(FILE *fplog, 
	     t_commrec *cr,
	     t_inputrec *ir, 
	     t_nrnb *nrnb,
	     gmx_wallcycle_t wcycle, 
	     gmx_localtop_t *top,
	     gmx_mtop_t *mtop, 
	     gmx_groups_t *groups,
	     t_state *state,
	     t_mdatoms *md, 
	     t_fcdata *fcd,
	     t_graph *graph, 
	     t_forcerec *fr,
	     gmx_vsite_t *vsite,
	     rvec mu_tot,
	     /*gmx_genborn_t *born,*/ 
	     bool bBornRadii,
	     real T,
	     int step,
	     tensor force_vir,
	     qhop_db_t db);

#endif	/* _qhop_h */
