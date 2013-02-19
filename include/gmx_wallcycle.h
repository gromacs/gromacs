/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2008, The GROMACS development team,
 * check out http://www.gromacs.org for more information.
 * Copyright (c) 2012,2013, by the GROMACS development team, led by
 * David van der Spoel, Berk Hess, Erik Lindahl, and including many
 * others, as listed in the AUTHORS file in the top-level source
 * directory and at http://www.gromacs.org.
 *
 * GROMACS is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public License
 * as published by the Free Software Foundation; either version 2.1
 * of the License, or (at your option) any later version.
 *
 * GROMACS is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with GROMACS; if not, see
 * http://www.gnu.org/licenses, or write to the Free Software Foundation,
 * Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA.
 *
 * If you want to redistribute modifications to GROMACS, please
 * consider that scientific software is very special. Version
 * control is crucial - bugs must be traceable. We will be happy to
 * consider code for inclusion in the official distribution, but
 * derived work must not be called official GROMACS. Details are found
 * in the README & COPYING files - if they are missing, get the
 * official version at http://www.gromacs.org.
 *
 * To help us fund GROMACS development, we humbly ask that you cite
 * the research papers on the package. Check out http://www.gromacs.org.
 */


#ifndef _gmx_wallcycle_h
#define _gmx_wallcycle_h

#include <stdio.h>
#include "visibility.h"
#include "typedefs.h"
#include "types/commrec.h"

#ifdef __cplusplus
extern "C" {
#endif

enum {
    ewcRUN, ewcSTEP, ewcPPDURINGPME, ewcDOMDEC, ewcDDCOMMLOAD,
    ewcDDCOMMBOUND, ewcVSITECONSTR, ewcPP_PMESENDX, ewcNS, ewcLAUNCH_GPU_NB,
    ewcMOVEX, ewcGB, ewcFORCE, ewcMOVEF, ewcPMEMESH,
    ewcPME_REDISTXF, ewcPME_SPREADGATHER, ewcPME_FFT, ewcPME_FFTCOMM, ewcPME_SOLVE,
    ewcPMEWAITCOMM, ewcPP_PMEWAITRECVF, ewcWAIT_GPU_NB_NL, ewcWAIT_GPU_NB_L, ewcNB_XF_BUF_OPS,
    ewcVSITESPREAD, ewcTRAJ, ewcUPDATE, ewcCONSTR, ewcMoveE, ewcROT, ewcROTadd,
    ewcTEST, ewcNR
};

enum {
    ewcsDD_REDIST, ewcsDD_GRID, ewcsDD_SETUPCOMM,
    ewcsDD_MAKETOP, ewcsDD_MAKECONSTR, ewcsDD_TOPOTHER,
    ewcsNBS_GRID_LOCAL, ewcsNBS_GRID_NONLOCAL,
    ewcsNBS_SEARCH_LOCAL, ewcsNBS_SEARCH_NONLOCAL,
    ewcsBONDED, ewcsNONBONDED, ewcsEWALD_CORRECTION,
    ewcsNB_X_BUF_OPS, ewcsNB_F_BUF_OPS,
    ewcsNR
};

gmx_bool wallcycle_have_counter(void);
/* Returns if cycle counting is supported */

GMX_LIBMD_EXPORT
gmx_wallcycle_t wallcycle_init(FILE *fplog, int resetstep, t_commrec *cr,
                               int nthreads_pp, int nthreads_pme);
/* Returns the wall cycle structure.
 * Returns NULL when cycle counting is not supported.
 */

GMX_LIBMD_EXPORT
void wallcycle_start(gmx_wallcycle_t wc, int ewc);
/* Starts the cycle counter (and increases the call count) */

void wallcycle_start_nocount(gmx_wallcycle_t wc, int ewc);
/* Starts the cycle counter without increasing the call count */

GMX_LIBMD_EXPORT
double wallcycle_stop(gmx_wallcycle_t wc, int ewc);
/* Stop the cycle count for ewc, returns the last cycle count */

GMX_LIBMD_EXPORT
void wallcycle_reset_all(gmx_wallcycle_t wc);
/* Resets all cycle counters to zero */

void wallcycle_sum(t_commrec *cr, gmx_wallcycle_t wc);
/* Sum the cycles over the nodes in cr->mpi_comm_mysim */

void wallcycle_print(FILE *fplog, int nnodes, int npme, double realtime,
                     gmx_wallcycle_t wc, wallclock_gpu_t *gpu_t);
/* Print the cycle and time accounting */

GMX_LIBMD_EXPORT
gmx_large_int_t wcycle_get_reset_counters(gmx_wallcycle_t wc);
/* Return reset_counters from wc struct */

GMX_LIBMD_EXPORT
void wcycle_set_reset_counters(gmx_wallcycle_t wc, gmx_large_int_t reset_counters);
/* Set reset_counters */

/* Uncomment the next line to get extra cycle counters of sub parts */
/* #define GMX_CYCLE_SUBCOUNTERS */

#ifdef GMX_CYCLE_SUBCOUNTERS

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
