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


typedef struct gmx_structurefactors gmx_structurefactors_t;

typedef struct structure_factor structure_factor_t;

typedef struct reduced_atom reduced_atom_t;

int * create_indexed_atom_type (reduced_atom_t * atm, int size);

void compute_structure_factor (structure_factor_t * sft, matrix box,
			       reduced_atom_t * red, int isize, real start_q,
			       real end_q, int group,real **sf_table);

gmx_structurefactors_t *gmx_structurefactors_init(const char *datfn);

void gmx_structurefactors_done(gmx_structurefactors_t *gsf);

int gmx_structurefactors_get_sf(gmx_structurefactors_t *gsf, int elem, real a[4], real b[4], real *c);

real **gmx_structurefactors_table(gmx_structurefactors_t *gsf,real momentum, real ref_k,
        real lambda, int n_angles);

void save_data (structure_factor_t * sft, const char *file, int ngrps,
                real start_q, real end_q, const output_env_t oenv);

double CMSF (gmx_structurefactors_t *gsf,int type,int nh,double lambda, double sin_theta);

int return_atom_type (const char *name,gmx_structurefactors_t *gsf);

void rearrange_atoms (reduced_atom_t * positions, t_trxframe *fr, atom_id * index,
		      int isize, t_topology * top, gmx_bool flag,gmx_structurefactors_t *gsf);

int do_scattering_intensity (const char* fnTPS, const char* fnNDX,
                             const char* fnXVG, const char *fnTRX,
                             const char* fnDAT,
                             real start_q,real end_q,
                             real energy,int ng,const output_env_t oenv);

t_complex *** rc_tensor_allocation(int x, int y, int z);

real **compute_scattering_factor_table (gmx_structurefactors_t *gsf,structure_factor_t * sft,int *nsftable);

#ifdef __cplusplus
}
#endif
#endif
