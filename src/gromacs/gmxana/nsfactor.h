/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2012,2013,2014, by the GROMACS development team, led by
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

#ifndef _nsfactor_h
#define _nsfactor_h

#include "index.h"
#include "types/simple.h"
#include "oenv.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef struct gmx_neutron_atomic_structurefactors_t {
    int       nratoms;
    int      *p;       /* proton number */
    int      *n;       /* neuton number */
    double   *slength; /* scattering length in fm */
    char    **atomnm;  /* atom symbol */
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
    int      qn;    /* number of items */
    double  *s;     /* scattering */
    double  *q;     /* q vectors */
    double   qstep; /* q increment */
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
