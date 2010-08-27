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
 * Copyright (c) 2001-2012, The GROMACS development team,
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

#ifndef _gmx_membed_h
#define _gmx_membed_h

#include "typedefs.h"

#ifdef __cplusplus
extern "C" {
#endif

/* information about scaling center */
typedef struct {
    rvec    xmin;         /* smallest coordinates of all embedded molecules */
    rvec    xmax;         /* largest coordinates of all embedded molecules */
    rvec    *geom_cent;   /* scaling center of each independent molecule to embed */
    int     pieces;       /* number of molecules to embed independently */
    int     *nidx;        /* n atoms for every independent embedded molecule (index in subindex) */
    atom_id **subindex;   /* atomids for independent molecule *
                           * atoms of piece i run from subindex[i][0] to subindex[i][nidx[i]] */
} pos_ins_t;

/* variables needed in do_md */
typedef struct {
    int   it_xy;          /* number of iterations (steps) used to grow something in the xy-plane */
    int   it_z;           /* same, but for z */
    real  xy_step;        /* stepsize used during resize in xy-plane */
    real  z_step;         /* same, but in z */
    rvec  fac;            /* initial scaling of the molecule to grow into the membrane */
    rvec  *r_ins;         /* final coordinates of the molecule to grow  */
    pos_ins_t *pos_ins;   /* scaling center for each piece to embed */
} gmx_membed_t;

/* initialisation of membed code */
void init_membed(FILE *fplog, gmx_membed_t *membed, int nfile, const t_filenm fnm[], 
                        gmx_mtop_t *mtop, t_inputrec *inputrec, t_state *state, t_commrec *cr,
                        real *cpt);

/* rescaling the coordinates voor de membed code */
void rescale_membed(int step_rel, gmx_membed_t *membed, rvec *x);

#ifdef __cplusplus
}
#endif

#endif
