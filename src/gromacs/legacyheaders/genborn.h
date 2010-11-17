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


#ifndef _genborn_h
#define _genborn_h

#include "typedefs.h"
#include "grompp.h"

#ifdef __cplusplus
extern "C" {
#endif

/* Still parameters - make sure to edit in genborn_sse.c too if you change these! */
#define STILL_P1  0.073*0.1              /* length        */
#define STILL_P2  0.921*0.1*CAL2JOULE    /* energy*length */
#define STILL_P3  6.211*0.1*CAL2JOULE    /* energy*length */
#define STILL_P4  15.236*0.1*CAL2JOULE
#define STILL_P5  1.254 

#define STILL_P5INV (1.0/STILL_P5)
#define STILL_PIP5  (M_PI*STILL_P5)


/* Initialise GB stuff */
int init_gb(gmx_genborn_t **p_born,
	    const t_commrec *cr, t_forcerec *fr, const t_inputrec *ir,
	    const gmx_mtop_t *mtop, real rgbradii, int gb_algorithm);


/* Born radii calculations, both with and without SSE acceleration */
int calc_gb_rad(t_commrec *cr, t_forcerec *fr, t_inputrec *ir,gmx_localtop_t *top,
				const t_atomtypes *atype, rvec x[], t_nblist *nl, gmx_genborn_t *born,t_mdatoms *md, t_nrnb     *nrnb);



/* Bonded GB interactions */								
real gb_bonds_tab(rvec x[], rvec f[], rvec fshift[], real *charge, real *p_gbtabscale,
				  real *invsqrta, real *dvda, real *GBtab, t_idef *idef, real epsilon_r,
				  real gb_epsilon_solvent, real facel, const t_pbc *pbc, const t_graph *graph);



void gb_pd_send(t_commrec *cr, real *send_data, int nr);


/* Functions for setting up the F_GB12,13,14 lists in grompp */
int 
init_gb_plist(t_params *p_list);

int 
convert_gb_params(gmx_ffparams_t *ffparams, t_functype ftype, t_params *gb_plist, t_ilist *il);



/* Functions for calculating adjustments due to ie chain rule terms */
void
calc_gb_forces(t_commrec *cr, t_mdatoms *md, gmx_genborn_t *born, gmx_localtop_t *top, const t_atomtypes *atype,
			   rvec x[], rvec f[], t_forcerec *fr,t_idef *idef,int gb_algorithm, int sa_algorithm, t_nrnb *nrnb, gmx_bool bRad,
			   const t_pbc *pbc, const t_graph *graph, gmx_enerdata_t *enerd);


int
make_gb_nblist(t_commrec *cr, int gb_algorithm, real gbcut,
	       rvec x[], matrix box,
	       t_forcerec *fr, t_idef *idef, t_graph *graph, gmx_genborn_t *born);

void 
make_local_gb(const t_commrec *cr, gmx_genborn_t *born, int gb_algorithm);

#ifdef __cplusplus
}
#endif

#endif /* _genborn_h */

