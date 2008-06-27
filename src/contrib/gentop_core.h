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

#ifndef _gentop_core_h
#define _gentop_core_h
	
#include <stdio.h>
#include "typedefs.h"
#include "pdbio.h"
#include "gentop_nm2type.h"
#include "gpp_nextnb.h"
#include "gpp_atomtype.h"

extern void calc_angles_dihs(t_params *ang,t_params *dih,rvec x[],bool bPBC,
			     matrix box);
			     
extern void set_force_const(t_params plist[],real kb,real kt,real kp,
			    bool bRound,bool bParam);
			    
extern real calc_dip(t_atoms *atoms,rvec x[]);

extern void delete_shell_interactions(t_params plist[F_NRE],t_atoms *atoms,
				      t_atomtype atype,t_nextnb *nnb,
				      t_excls excls[]);
				      
extern void dump_hybridization(FILE *fp,t_atoms *atoms,int nbonds[]);

extern void reset_q(t_atoms *atoms);

extern void print_rtp(char *filenm,char *title,t_atoms *atoms,
		      t_params plist[],int cgnr[],int nbts,int bts[]);
		      
extern void mk_bonds(gentop_nm2t nmt,t_atoms *atoms,rvec x[],
		     gmx_conect gc,t_params *bond,int nbond[],char *ff,
		     bool bPBC,matrix box,gmx_atomprop_t aps,real tol);
		     
extern t_atomtype set_atom_type(t_symtab *tab,t_atoms *atoms,t_params *bonds,
				int *nbonds,gentop_nm2t nm2t,
				gmx_atomprop_t aps);
		     
extern void add_shells(gentop_nm2t nm2t,t_atoms **atoms,
		       t_atomtype atype,t_params plist[],
		       rvec **x,t_symtab *symtab,t_excls **excls);
		       
extern void symmetrize_charges(t_atoms *atoms,t_atomtype atype,
			       t_params *bonds,gmx_atomprop_t aps);

enum { ecgGroup, ecgAtom, ecgNeutral, ecgNR };

extern int *generate_charge_groups(int cgtp,t_atoms *atoms,t_atomtype atype,
				   t_params *bonds,t_params *pols,
				   bool bUsePDBcharge,real *qtot,real *mutot);

extern void sort_on_charge_groups(int *cgnr,t_atoms *atoms,t_params plist[],
				  t_atomtype type,rvec x[]);

#endif
