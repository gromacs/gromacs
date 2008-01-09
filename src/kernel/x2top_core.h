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
 * Gallium Rubidium Oxygen Manganese Argon Carbon Silicon
 */

#ifndef _x2top_h
#define _x2top_h
	
#include <stdio.h>
#include "grompp.h"
	
extern void calc_angles_dihs(t_params *ang,t_params *dih,rvec x[],bool bPBC,
			     matrix box);
			     
extern void set_force_const(t_params plist[],real kb,real kt,real kp,
			    bool bRound,bool bParam);
			    
extern int *set_cgnr(t_atoms *atoms,bool bUsePDBcharge,real *qtot,real *mtot);

extern real calc_dip(t_atoms *atoms,rvec x[]);

extern void delete_shell_interactions(t_params plist[F_NRE],t_atoms *atoms,
				      t_atomtype *atype,t_nextnb *nnb,
				      t_excls excls[]);
				      
extern void dump_hybridization(FILE *fp,t_atoms *atoms,int nbonds[]);

extern void reset_q(t_atoms *atoms);

extern void print_rtp(char *filenm,char *title,t_atoms *atoms,
		      t_params plist[],int cgnr[],int nbts,int bts[]);
		      
extern void mk_bonds(int nnm,t_nm2type nmt[],
		     t_atoms *atoms,rvec x[],t_params *bond,int nbond[],char *ff,
		     bool bPBC,matrix box,void *atomprop,real tol);
		     
extern t_atomtype *set_atom_type(t_symtab *tab,t_atoms *atoms,t_params *bonds,
				 int *nbonds,int nnm,t_nm2type nm2t[]);
		     
#endif
