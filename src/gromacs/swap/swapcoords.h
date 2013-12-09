/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2013, The GROMACS development team.
 * Copyright (c) 2013, by the GROMACS development team, led by
 * Mark Abraham, David van der Spoel, Berk Hess, and Erik Lindahl,
 * and including many others, as listed in the AUTHORS file in the
 * top-level source directory and at http://www.gromacs.org.
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

#ifndef _swapcoords_h
#define _swapcoords_h

#include "typedefs.h"
#include "types/commrec.h"

extern void init_swapcoords(
        FILE             *fplog,    /* general output file md.log */
        gmx_bool         bVerbose,
        t_inputrec       *ir,
        const char       *fn,       /* output file name for swap data */
        gmx_mtop_t       *mtop,
        rvec             x[],       /* the initial positions */
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

#endif
