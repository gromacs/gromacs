/* -*- mode: c; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; c-file-style: "stroustrup"; -*-
 * $Id: gentop_core.h,v 1.8 2009/02/02 21:11:11 spoel Exp $
 * 
 *                This source code is part of
 * 
 *                 G   R   O   M   A   C   S
 * 
 *          GROningen MAchine for Chemical Simulations
 * 
 *                        VERSION 4.0.99
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
 * Groningen Machine for Chemical Simulation
 */

#ifndef _gentop_core_h
#define _gentop_core_h

	
#include <stdio.h>
#include "typedefs.h"
#include "pdbio.h"
#include "gentop_nm2type.h"
#include "gpp_nextnb.h"
#include "gpp_atomtype.h"
#include "poldata.h"

extern void calc_angles_dihs(t_params *ang,t_params *dih,rvec x[],gmx_bool bPBC,
                             matrix box);
			     
extern real calc_dip(t_atoms *atoms,rvec x[]);

extern void dump_hybridization(FILE *fp,t_atoms *atoms,int nbonds[]);

extern void reset_q(t_atoms *atoms);

extern void print_rtp(char *filenm,char *title,t_atoms *atoms,
                      t_params plist[],int cgnr[],int nbts,int bts[]);
		      
extern void mk_bonds(gmx_poldata_t pd,t_atoms *atoms,rvec x[],
                     gmx_conect gc,t_params plist[],int nbond[],
                     gmx_bool bH14,gmx_bool bAllDihedrals,gmx_bool bRemoveDoubleDihedrals,
                     int nexcl,t_excls **excls,
                     gmx_bool bPBC,matrix box,gmx_atomprop_t aps,real tol);
		     
extern gpp_atomtype_t set_atom_type(FILE *fp,char *molname,
                                    t_symtab *tab,t_atoms *atoms,t_params *bonds,
                                    int nbonds[],char **smnames,
                                    gmx_poldata_t pd,gmx_atomprop_t aps,
                                    rvec x[],t_pbc *pbc,real th_toler,
                                    real ph_toler,gentop_vsite_t gvt);
		     
extern void add_shells(gmx_poldata_t pd,int maxatom,t_atoms *atoms,
                       gpp_atomtype_t atype,t_params plist[],
                       rvec *x,t_symtab *symtab,t_excls **excls,
                       char **smnames);
		       
extern int *symmetrize_charges(gmx_bool bQsym,
                               t_atoms *atoms,t_params *bonds,gmx_poldata_t pd,
                               gmx_atomprop_t aps,char *symm_string);

enum { ecgAtom, ecgGroup, ecgNeutral, ecgNR };

extern int *generate_charge_groups(int cgtp,t_atoms *atoms,
                                   t_params *bonds,t_params *pols,
                                   gmx_bool bUsePDBcharge,
                                   real *qtot,real *mtot);

extern void sort_on_charge_groups(int *cgnr,t_atoms *atoms,t_params plist[],
                                  rvec x[],t_excls excls[],
                                  char *smnames[],const char *ndxout,
                                  int nmol);

#endif
