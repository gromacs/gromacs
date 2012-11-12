/*
 *
 *                This source code is part of
 *
 *                 G   R   O   M   A   C   S
 *
 *          GROningen MAchine for Chemical Simulations
 *
 *                        VERSION 3.2.0
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

#ifndef _nsfactor_h
#define _nsfactor_h

#include "index.h"
#include "types/simple.h"
#include "gmxcomplex.h"
#include "oenv.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef struct gmx_neutron_atomic_structurefactors_t {
    int     nratoms;
    int    *p;       /* proton number */
    int    *n;       /* neuton number */
    double *slength; /* scattering length in fm */
    char  **atomnm;  /* atom symbol */
} gmx_neutron_atomic_structurefactors_t;

typedef struct gmx_sans_t {
    t_topology *top;     /* topology */
    double     *slength; /* scattering length for this topology */
} gmx_sans_t;

typedef struct gmx_radial_distribution_histogram_t {
    int     grn;      /* number of bins */
    double  binwidth; /* bin size */
    double *r;        /* Distances */
    double *gr;       /* Probability */
} gmx_radial_distribution_histogram_t;

typedef struct gmx_static_structurefactor_t {
    int     qn;    /* number of items */
    double *s;     /* scattering */
    double *q;     /* q vectors */
    double  qstep; /* q increment */
} gmx_static_structurefactor_t;

void check_binwidth(real binwidth);

void check_mcover(real mcover);

void normalize_probability(int n, double *a);

gmx_neutron_atomic_structurefactors_t *gmx_neutronstructurefactors_init(const char *datfn);

gmx_sans_t *gmx_sans_init(t_topology *top, gmx_neutron_atomic_structurefactors_t *gnsf);

gmx_radial_distribution_histogram_t *calc_radial_distribution_histogram  (gmx_sans_t  *gsans,
                                                                          rvec        *x,
                                                                          matrix       box,
                                                                          atom_id     *index,
                                                                          int          isize,
                                                                          double       binwidth,
                                                                          gmx_bool     bMC,
                                                                          gmx_bool     bNORM,
                                                                          real         mcover,
                                                                          unsigned int seed);

gmx_static_structurefactor_t *convert_histogram_to_intensity_curve (gmx_radial_distribution_histogram_t *pr, double start_q, double end_q, double q_step);


#ifdef __cplusplus
}
#endif
#endif
