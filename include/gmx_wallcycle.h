#ifndef _gmx_wallcycle_h
#define _gmx_wallcycle_h

#include <stdio.h>
#include "typedefs.h"

enum { ewcRUN, ewcDOMDEC, ewcNS, ewcFORCE, ewcPMEMESH, ewcPMEMESH_SEP, ewcUPDATE, ewcNR };

extern gmx_wallcycle_t init_wallcycle(void);
/* Returns the wall cycle structure.
 * Returns NULL when cycle counting is not supported.
 */

extern void wallcycle_start(gmx_wallcycle_t wc, int ewc);
/* Set the start cycle count for ewc */

extern void wallcycle_stop(gmx_wallcycle_t wc, int ewc);
/* Stop the cycle count for ewc */

extern void wallcycle_sum(t_commrec *cr, gmx_wallcycle_t wc,double cycles[]);
/* Sum the cycles over the nodes in cr->mpi_comm_mysim */

extern void wallcycle_print(FILE *fplog, int nnodes, int npme, double realtime,
			    gmx_wallcycle_t wc, double cycles[]);
/* Print the cycle and time accounting */

#endif /* _gmx_wallcycle_h */
