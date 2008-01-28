#ifndef _gmx_wallcycle_h
#define _gmx_wallcycle_h

#include <stdio.h>
#include "typedefs.h"

enum { ewcRUN, ewcSTEP, ewcPPDURINGPME, ewcDOMDEC, ewcVSITECONSTR, ewcPP_PMESENDX, ewcMOVEX, ewcNS, ewcFORCE, ewcMOVEF, ewcPMEMESH, ewcPMEMESH_SEP, ewcPMEWAITCOMM, ewcPP_PMEWAITRECVF, ewcVSITESPREAD, ewcTRAJ, ewcUPDATE, ewcCONSTR, ewcMoveE, ewcTEST, ewcNR };

extern bool wallcycle_have_counter(void);
/* Returns if cycle counting is supported */

extern gmx_wallcycle_t wallcycle_init(FILE *fplog,t_commrec *cr);
/* Returns the wall cycle structure.
 * Returns NULL when cycle counting is not supported.
 */

extern void wallcycle_start(gmx_wallcycle_t wc, int ewc);
/* Set the start cycle count for ewc */

extern double wallcycle_stop(gmx_wallcycle_t wc, int ewc);
/* Stop the cycle count for ewc, returns the last cycle count */

extern void wallcycle_sum(t_commrec *cr, gmx_wallcycle_t wc,double cycles[]);
/* Sum the cycles over the nodes in cr->mpi_comm_mysim */

extern void wallcycle_print(FILE *fplog, int nnodes, int npme, double realtime,
			    gmx_wallcycle_t wc, double cycles[]);
/* Print the cycle and time accounting */

#endif /* _gmx_wallcycle_h */
