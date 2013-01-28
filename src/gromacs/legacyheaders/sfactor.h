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

#ifndef _sfactor_h
#define _sfactor_h

 
#include "index.h"
#include "types/simple.h"
#include "gmxcomplex.h"
#include "oenv.h"



#ifdef __cplusplus
extern "C" {
#endif

/*
 * Definitions of data structures
 */

typedef struct gmx_structurefactors_t {
    int nratoms;
    int *p;        /* proton number */
    int *n;        /* neutron number */
    /* Parameters for the Cromer Mann fit */
    real **a;    /* parameter a */
    real **b;    /* parameter b */
    real *c;       /* parameter c */
    /* Parameters for Neutrons */
    real  *coh_b; /* scattering length in fm */
    /* generic stuff */
    char **atomnm;  /* atomname */

} gmx_structurefactors_t;

typedef struct reduced_atom_t {
    rvec x;
    int  t;
} reduced_atom_t ;

typedef struct structure_factor_t
{
  int     n_angles;
  int     n_groups;
  double  lambda;
  double  energy;
  double  momentum;
  double  ref_k;
  double  **F;
  int     nSteps;
  int     total_n_atoms;
} structure_factor_t;

typedef struct gmx_sans_t {
    t_topology *top; /* topology */
    double *slength; /* scattering length for this topology */
} gmx_sans_t;

typedef struct gmx_radial_distribution_histogram_t {
    int     grn; /* number of bins */
    double binwidth; /* bin size */
    double *r; /* Distances */
    double *gr; /* Probability */
} gmx_radial_distribution_histogram_t;

typedef struct gmx_sans_structurefactor_t {
    int     qn; /* number of items */
    double  *s; /* scattering */
    double  *q; /* q vectors */
    double  qstep; /* q increment */
} gmx_sans_structurefactor_t;

typedef struct gmx_nse_structurefcator_t {
    int     tn; /*  number of frames */
    double   q; /*  q for this s(q=const(t)) */
    double  *s; /*  scattering */
} gmx_nse_structurefactor_t;

typedef struct gmx_nse_t {
    gmx_sans_t                            *sans; /*  scattering params */
    gmx_radial_distribution_histogram_t   **gr;  /*  array of gr */
    gmx_sans_structurefactor_t           **sq; /*  array of sq */
    gmx_nse_structurefactor_t          **sqt; /*  s(q(t))  array */
    rvec       **x; /* md coordinates */
    real        *t; /*  md time */
    real        *dt; /* time for correlation */
    matrix      *box; /* box for current coordinates */
    int        nrframes;
    int        sqtn;
} gmx_nse_t;

/*
 * Common functions for SANS and SAXS
 */

gmx_structurefactors_t *gmx_structurefactors_init(const char *datfn);

void done_gmx_structurefactors(gmx_structurefactors_t *gsf);

int *create_indexed_atom_type (reduced_atom_t * atm, int size);

int return_atom_type (const char *name,gmx_structurefactors_t *gsf);

void rearrange_atoms (reduced_atom_t * positions, t_trxframe *fr, atom_id * index,
              int isize, t_topology * top, gmx_bool flag,gmx_structurefactors_t *gsf);

void check_binwidth(real binwidth);

void check_mcover(real mcover);

void normalize_probability(int n, double *a);

/*
 * SAXS
 */

void compute_structure_factor (structure_factor_t * sft, matrix box,
			       reduced_atom_t * red, int isize, real start_q,
			       real end_q, int group,real **sf_table);

int gmx_structurefactors_get_sf(gmx_structurefactors_t *gsf, int elem, real a[4], real b[4], real *c);

real **gmx_structurefactors_table(gmx_structurefactors_t *gsf,real momentum, real ref_k,
        real lambda, int n_angles);

void save_data (structure_factor_t * sft, const char *file, int ngrps,
                real start_q, real end_q, const output_env_t oenv);

double CMSF (gmx_structurefactors_t *gsf,int type,int nh,double lambda, double sin_theta);

t_complex *** rc_tensor_allocation(int x, int y, int z);

real **compute_scattering_factor_table (gmx_structurefactors_t *gsf,structure_factor_t * sft,int *nsftable);

/*
 * SANS
 */

gmx_sans_t *gmx_sans_init(t_topology *top, gmx_structurefactors_t *sf);

void done_sans(gmx_sans_t *gsans);

gmx_radial_distribution_histogram_t *calc_radial_distribution_histogram  (gmx_sans_t *gsans,
                            rvec *x, rvec *xf,
                            matrix box, matrix boxf,
                            atom_id *index,
                            int isize,
                            double binwidth,
                            gmx_bool bMC, gmx_bool bNSE,
                            real mcover,
                            unsigned int seed);

gmx_sans_structurefactor_t *convert_histogram_to_intensity_curve (gmx_radial_distribution_histogram_t *pr, double start_q, double end_q, double q_step);

#ifdef __cplusplus
}
#endif
#endif
