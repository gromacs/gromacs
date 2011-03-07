#ifndef _qhoprec_h
#define _qhoprec_h


#include "types/qhoprec.h"
#include "types/topology.h"
#include "types/mdatom.h"
extern t_qhoprec *mk_qhoprec(void);

/* Sets the interactions according to the hydrogen exstence map.
 * This requires a finalized t_mdatoms. */
extern void finalize_qhoprec(t_qhoprec *qr, gmx_localtop_t *top, t_mdatoms *md,t_commrec *cr);

#endif	/* _qhoprec_h */
