/*
 * 
 *                This source code is part of
 * 
 *                 G   R   O   M   A   C   S
 * 
 *          GROningen MAchine for Chemical Simulations
 * 
 * Written by David van der Spoel, Erik Lindahl, Berk Hess, and others.
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2008, The GROMACS development team,
 * check out http://www.gromacs.org for more information.
 
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 * 
 * If you want to redistribute modifications, please consider that
 * scientific software is very special. Version control is crucial -
 * bugs must be traceable. We will be happy to consider code for
 * inclusion in the official distribution, but derived work must not
 * be called official GROMACS. Details are found in the README & COPYING
 * files - if they are missing, get the official version at www.gromacs.org.
 * 
 * To help us fund GROMACS development, we humbly ask that you cite
 * the papers on the package - you can find them in the top README file.
 * 
 * For more info, check our website at http://www.gromacs.org
 * 
 * And Hey:
 * Gallium Rubidium Oxygen Manganese Argon Carbon Silicon
 */


#ifndef _gmx_wallcycle_h
#define _gmx_wallcycle_h

#include <stdio.h>
#include "typedefs.h"

#ifdef __cplusplus
extern "C" {
#endif

enum { ewcRUN, ewcSTEP, ewcPPDURINGPME, ewcDOMDEC, ewcDDCOMMLOAD,
       ewcDDCOMMBOUND, ewcVSITECONSTR, ewcPP_PMESENDX, ewcNS, ewcLAUNCH_GPU_NB,
       ewcMOVEX, ewcGB, ewcFORCE, ewcMOVEF, ewcPMEMESH,
       ewcPME_REDISTXF, ewcPME_SPREADGATHER, ewcPME_FFT, ewcPME_FFTCOMM, ewcPME_SOLVE,
       ewcPMEWAITCOMM, ewcPP_PMEWAITRECVF, ewcWAIT_GPU_NB_NL, ewcWAIT_GPU_NB_L, ewcNB_XF_BUF_OPS,
       ewcVSITESPREAD, ewcTRAJ, ewcUPDATE, ewcCONSTR, ewcMoveE, ewcROT, ewcROTadd,
       ewcTEST, ewcNR };

enum { ewcsDD_REDIST, ewcsDD_GRID, ewcsDD_SETUPCOMM,
       ewcsDD_MAKETOP, ewcsDD_MAKECONSTR, ewcsDD_TOPOTHER,
       ewcsNBS_GRID_LOCAL, ewcsNBS_GRID_NONLOCAL,
       ewcsNBS_SEARCH_LOCAL, ewcsNBS_SEARCH_NONLOCAL,
       ewcsBONDED, ewcsNONBONDED, ewcsEWALD_CORRECTION,
       ewcsNB_X_BUF_OPS, ewcsNB_F_BUF_OPS,
       ewcsNR };

gmx_bool wallcycle_have_counter(void);
/* Returns if cycle counting is supported */

gmx_wallcycle_t wallcycle_init(FILE *fplog, int resetstep, t_commrec *cr, 
                               int nthreads_pp, int nthreads_pme);
/* Returns the wall cycle structure.
 * Returns NULL when cycle counting is not supported.
 */

void wallcycle_start(gmx_wallcycle_t wc, int ewc);
/* Starts the cycle counter (and increases the call count) */

void wallcycle_start_nocount(gmx_wallcycle_t wc, int ewc);
/* Starts the cycle counter without increasing the call count */

double wallcycle_stop(gmx_wallcycle_t wc, int ewc);
/* Stop the cycle count for ewc, returns the last cycle count */

void wallcycle_reset_all(gmx_wallcycle_t wc);
/* Resets all cycle counters to zero */

void wallcycle_sum(t_commrec *cr, gmx_wallcycle_t wc);
/* Sum the cycles over the nodes in cr->mpi_comm_mysim */

void wallcycle_print(FILE *fplog, int nnodes, int npme, double realtime,
			    gmx_wallcycle_t wc, wallclock_gpu_t *gpu_t);
/* Print the cycle and time accounting */

gmx_large_int_t wcycle_get_reset_counters(gmx_wallcycle_t wc);
/* Return reset_counters from wc struct */

void wcycle_set_reset_counters(gmx_wallcycle_t wc, gmx_large_int_t reset_counters);
/* Set reset_counters */

/* Uncomment the next line to get extra cycle counters of sub parts */
/* #define GMX_CYCLE_SUB */

#ifdef GMX_CYCLE_SUB

void wallcycle_sub_start(gmx_wallcycle_t wc, int ewcs);
/* Set the start sub cycle count for ewcs */

void wallcycle_sub_stop(gmx_wallcycle_t wc, int ewcs);
/* Stop the sub cycle count for ewcs */

#else
/* Define the counter call to nothing to avoid any effect on performance */
#define wallcycle_sub_start(wc, ewcs)
#define wallcycle_sub_stop(wc, ewcs)
#endif

#ifdef __cplusplus
}
#endif

#endif /* _gmx_wallcycle_h */
