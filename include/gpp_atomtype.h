/*
 * $Id$
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

#ifndef _gpp_atomtype_h
#define _gpp_atomtype_h

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <stdio.h>
#include "typedefs.h"
#include "macros.h"

typedef struct gpp_atomtype *t_atomtype;

extern int get_atomtype_type(char *str,t_atomtype at);
/* Return atomtype corresponding to case-insensitive str
   or NOTSET if not found */

extern int get_atomtype_ntypes(t_atomtype at);
/* Return number of atomtypes */

extern char *get_atomtype_name(int nt, t_atomtype at);
/* Return name corresponding to atomtype nt, or NULL if not found */

extern real get_atomtype_massA(int nt,t_atomtype at);
extern real get_atomtype_massB(int nt,t_atomtype at);
extern real get_atomtype_qA(int nt,t_atomtype at);
extern real get_atomtype_qB(int nt,t_atomtype at);
extern real get_atomtype_radius(int nt,t_atomtype at);
extern real get_atomtype_vol(int nt,t_atomtype at);
extern real get_atomtype_surftens(int nt,t_atomtype at);
extern real get_atomtype_gb_radius(int nt,t_atomtype at);
extern real get_atomtype_S_hct(int nt,t_atomtype at);
extern int get_atomtype_ptype(int nt,t_atomtype at);
extern int get_atomtype_batype(int nt,t_atomtype at);
extern int get_atomtype_atomnumber(int nt,t_atomtype at);

/* Return the above variable for atomtype nt, or NOTSET if not found */

extern real get_atomtype_nbparam(int nt,int param,t_atomtype at);
/* Similar to the previous but returns the paramth parameter or NOTSET */

extern t_atomtype init_atomtype(void);
/* Return a new atomtype structure */

extern void done_atomtype(t_atomtype *at);
/* Free the memory in the structure */

extern int set_atomtype(int nt,t_atomtype at,t_symtab *tab,
			t_atom *a,char *name,t_param *nb,
			int bondatomtype,
			real radius,real vol,real surftens,int atomnumber,
			real gb_radius, real S_hct);
/* Set the values of an existing atom type nt. Returns nt on success or
   NOTSET on error. */	

int
set_atomtype_gbparam(t_atomtype at, int i,
		     real radius,real vol,real surftens,
		     real gb_radius, real S_hct);

extern int add_atomtype(t_atomtype at,t_symtab *tab,
			t_atom *a,char *name,t_param *nb,
			int bondatomtype,
			real radius,real vol,real surftens,real atomnumber,
			real gb_radius, real S_hct);
/* Add a complete new atom type to an existing atomtype structure. Returns
   the number of the atom type. */

extern void print_at (FILE * out, t_atomtype at);
/* Print an atomtype record to a text file */

extern void renum_atype(t_params plist[],gmx_mtop_t *mtop,
			int *wall_atomtype,
			t_atomtype at,bool bVerbose);
			
extern void copy_atomtype_atomtypes(t_atomtype atype,t_atomtypes *atypes);
/* Copy from one structure to another */
#endif	/* _gpp_atomtype_h */








