/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2004, The GROMACS development team,
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

#ifndef _pull_h
#define _pull_h
#include "visibility.h"
#include "vec.h"
#include "typedefs.h"

#ifdef __cplusplus
extern "C" {
#endif


/* This file contains datatypes and function declarations necessary
   for mdrun to interface with the pull code */

/* Get the distance to the reference and deviation for pull group g */
GMX_LIBMD_EXPORT
void get_pullgrp_distance(t_pull *pull, t_pbc *pbc, int g, double t,
                          dvec dr, dvec dev);

/* Set the all the pull forces to zero */
void clear_pull_forces(t_pull *pull);

/* Determine the COM pull forces and add them to f, return the potential */
real pull_potential(int ePull, t_pull *pull, t_mdatoms *md, t_pbc *pbc,
                    t_commrec *cr, double t, real lambda,
                    rvec *x, rvec *f, tensor vir, real *dvdlambda);

/* Constrain the coordinates xp in the directions in x
 * and also constrain v when v!=NULL.
 */
void pull_constraint(t_pull *pull, t_mdatoms *md, t_pbc *pbc,
                     t_commrec *cr, double dt, double t,
                     rvec *x, rvec *xp, rvec *v, tensor vir);

/* Make a selection of the home atoms for all pull groups.
 * Should be called at every domain decomposition.
 */
GMX_LIBMD_EXPORT
void dd_make_local_pull_groups(gmx_domdec_t *dd,
                               t_pull *pull, t_mdatoms *md);

/* get memory and initialize the fields of pull that still need it, and
   do runtype specific initialization */
GMX_LIBMD_EXPORT
void init_pull(FILE              *fplog,
               t_inputrec        *ir,       /* the inputrec */
               int                nfile,
               const t_filenm     fnm[],    /* standard filename struct */
               gmx_mtop_t        *mtop,     /* the topology of the whole system */
               t_commrec        * cr,       /* struct for communication info */
               const output_env_t oenv,     /* output options */
               real               lambda,   /* FEP lambda */
               gmx_bool           bOutFile, /* open output files */
               unsigned long      Flags);

/* Close the pull output files */
GMX_LIBMD_EXPORT
void finish_pull(FILE *fplog, t_pull *pull);

/* Print the pull output (x and/or f) */
GMX_LIBMD_EXPORT
void pull_print_output(t_pull *pull, gmx_large_int_t step, double time);

/* In pullutil.c */

/* Calculates centers of mass all pull groups */
GMX_LIBMD_EXPORT
void pull_calc_coms(t_commrec *cr,
                    t_pull    *pull, /* the pull group */
                    t_mdatoms *md,   /* all atoms */
                    t_pbc     *pbc,
                    double     t,    /* only used for cylinder ref. */
                    rvec       x[],  /* local coordinates */
                    rvec      *xp    /* updated x, can be NULL */
                    );

#ifdef __cplusplus
}
#endif

#endif
