#ifndef _TITRATION_H
#define _TITRATION_H

#ifdef HAVE_IDENT
#ident	"@(#) qhop.h 1 28/2/01"
#endif /* HAVE_IDENT */
#include <stdio.h>
#include "typedefs.h"
#include "pbc.h"
#include "network.h"
#include "tgroup.h"
#include "titrationrec.h"
#include "resall.h"

/** \brief Initializes the titration structure inside forcerec.
 * 
 * What it does:
 *   - Initializes titration structure inside forcerec
 *   - Picks out the atoms that are titratable
 *   - Figures out which residue subtypes to use at step 0,
 *     based on the H existence map.
 *   - Creates a bqhopdonor (andbqhopacceptor) array in
 *     the t_mdatoms to be used by nbsearch,
 *   - Completes the qhop residues array
 *   - Reads in the hopping and force field parameters.
 * \return nothing.
 */
extern void init_titration(FILE *fplog,const char *ff,
			   t_commrec *cr, gmx_mtop_t *mtop, 
			   t_inputrec *ir, t_forcerec *fr,
			   matrix box, t_mdatoms *md);

/** \brief Identifies potential hops calculates probabilities.
 *
 * do_qhop() identifies acceptors in hopping range of the donors,
 * calculates relevant energies and subsequently probabilities for
 * every hop. The list of hops is scrambled according to the qhopmode
 * and the hops are tested against random numbers. The actual hopping
 * will be done later.
 * \return the QMMM energy add into the system.
 */
extern real do_titration(FILE *fplog, 
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
			 gmx_constr_t constr,
			 gmx_vsite_t *vsite,
			 rvec mu_tot,
			 gmx_bool bBornRadii,
			 real Temperature,
			 gmx_large_int_t step,
			 gmx_ekindata_t *ekindata,
			 tensor force_vir,
			 rvec *f_old);

extern void fold_inactive_protons(titration_t T, rvec x[], rvec v[]);

/* Sets the interactions according to the hydrogen exstence map.
 * This requires a finalized t_mdatoms. */
extern void finalize_titration(titration_t T, gmx_localtop_t *top, t_mdatoms *md,t_commrec *cr);


#endif	/* _qhop_h */
