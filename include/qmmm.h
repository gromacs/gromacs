#ifndef _QMMM_h
#define _QMMM_h


#ifdef HAVE_IDENT
#ident	"@(#) QMMM.h 1 28/2/01"
#endif /* HAVE_IDENT */
#include "typedefs.h"
#include "pbc.h"
#include "network.h"
#include "tgroup.h"

void 
atomic_number(int nr, char ***atomtype, int *nucnum);

extern t_QMMMrec *mk_QMMMrec(void);
/* allocates memory for QMMMrec */

extern void init_QMMMrec(t_commrec *cr,
			 matrix box,
			 t_topology *top,
			 t_inputrec *ir,
			 t_forcerec *fr);

/* init_QMMMrec initializes the QMMM record. From
 * topology->atoms.atomname and topology->atoms.atomtype the atom
 * names and types are read; from inputrec->QMcharge
 * resp. inputrec->QMmult the nelecs and multiplicity are determined
 * and md->cQMMM gives numbers of the MM and QM atoms 
 */

extern void update_QMMMrec(t_commrec *cr,
			   t_forcerec *fr,
			   rvec x[],
			   t_mdatoms *md,
			   matrix box,
			   t_topology *top);

/* update_QMMMrec fills the MM stuff in QMMMrec. The MM atoms are
 * taken froom the neighbourlists of the QM atoms. In a QMMM run this
 * routine should be called at every step, since it updates the MM
 * elements of the t_QMMMrec struct.  
 */

extern real calculate_QMMM(t_commrec *cr,
			   rvec x[], rvec f[],
			   t_forcerec *fr,
			   t_mdatoms *md);

/* QMMM computes the QM forces. This routine makes either function
 * calls to gmx QM routines (derived from MOPAC7 (semi-emp.) and MPQC
 * (ab initio)) or generates input files for an external QM package
 * (listed in QMMMrec.QMpackage). The binary of the QM package is
 * called by system().
 */

#endif	/* _QMMM_h */

