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
* Copyright (c) 2001-2004, The GROMACS development team,
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
* Gromacs Runs On Most of All Computer Systems
*/

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#ifdef __cplusplus
extern "C" {
#endif

extern void init_swapcoords(
            FILE             *fplog, /* general output file md.log */
            gmx_bool         bVerbose,
            t_inputrec       *ir,
            const char       *fn,    /* output file name for swap data */
            gmx_mtop_t       *mtop,
            rvec             x[],
            matrix           box,
            swapstate_t      *swapstate,
            t_commrec        *cr,
            const output_env_t oenv,
            unsigned long    Flags);

extern void dd_make_local_swap_groups(gmx_domdec_t *dd, t_swapcoords *si_pub);


extern gmx_bool do_swapcoords(
        t_commrec        *cr,
        gmx_large_int_t  step,
        real             t,
        t_inputrec       *ir,
        rvec             x[],            /* positions of home particles */
        matrix           box,
        gmx_mtop_t       *mtop,
        gmx_bool         bVerbose,
        gmx_bool         bRerun);

extern void swapcoords_get_charges(
        t_commrec        *cr,
        FILE             *fplog,
        gmx_large_int_t  step,
        t_inputrec       *ir,
        gmx_mtop_t       *mtop);

#ifdef __cplusplus
}
#endif
