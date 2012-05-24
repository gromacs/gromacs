/*
 * $Id: molprop_util.h,v 1.18 2009/05/27 13:44:55 spoel Exp $
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

#ifndef _molprop_util_h
#define _molprop_util_h

#include "grompp.h"
#include "atomprop.h"
#include "molprop.h"
#include "poldata.h"
#include "molselect.h"

enum { iqmExp, iqmBoth, iqmQM, iqmNR };

/* Returns 1 when OK, 0 otherwise. Level of theory and/or conformation
   and/or type (elec, ESP, RESP) may be set. */
extern int get_val(gmx_molprop_t mp,int expref,char *type,int emp,
		   double *value,double *error,double vec[3],
		   tensor quadrupole);
		   
extern int mp_get_prop(gmx_molprop_t mp,int eGP,int iQM,char *lot,
		       char *conf,char *type,double *value);

/* mylot returns the used level of theory. mylot and ref may be NULL */
extern int mp_get_prop_ref(gmx_molprop_t mp,int emp,int iQM,char *lot,
			   char *conf,char *type,double *value,double *err,
			   char **ref,char **mylot,
			   double vec[3],tensor quadrupole);

extern int gmx_molprop_get_calc_lot(gmx_molprop_t mpt,char *lot);

extern void generate_formula(int nmol,gmx_molprop_t mp[],gmx_atomprop_t ap);

/* Returns TRUE if the present molecule has the indicated composition */
extern gmx_bool gmx_molprop_support_composition(gmx_molprop_t mp,char *composition);

extern void generate_composition(int nmol,gmx_molprop_t mp[],gmx_poldata_t pd,
                                 gmx_atomprop_t ap,gmx_bool bForceGenComp,
                                 double th_toler,double phi_toler);

extern gmx_molprop_t atoms_2_molprop(char *molname,t_atoms*atoms,char **smnames,
				     gmx_atomprop_t ap,gmx_poldata_t pd,gmx_bool bForce,
				     double th_toler,double ph_toler);

/* Return number of atoms, 0 means failure */
extern int molprop_2_atoms(gmx_molprop_t mp,gmx_atomprop_t ap,
			   t_symtab *tab,const char *lot,
			   t_atoms *atoms,const char *q_algorithm,
			   rvec **x);
				     				     
extern gmx_molprop_t *merge_xml(int nfile,char **infiles,char *outf,
				char *sorted,char *doubles,int *nmolprop,
				gmx_atomprop_t ap,gmx_poldata_t pd,
				gmx_bool bForceMerge,gmx_bool bForceGenComp,
				double th_toler,double ph_toler);

enum { empSORT_Molname, empSORT_Formula, empSORT_Composition, empSORT_Selection, empSORT_NR };
	       
extern void gmx_molprop_sort(int np,gmx_molprop_t mp[],int alg,
			     gmx_atomprop_t ap,gmx_molselect_t gms);
/* gms may be NULL if the alg is not empSORT_Selection */

typedef struct {
  int  n;
  int  *count;
  char **method,**basis,**type,**lot;
  int  nconf;
  char **conf;
} t_qmcount;

/* Check the available molprops to see what kind of calculations are stored in there */
extern t_qmcount *find_calculations(int np,gmx_molprop_t mp[],int emp,char *fc_str);
				
/* Routines for sending and receiving a molprop */
extern void gmx_molprop_send(t_commrec *cr,int dest,gmx_molprop_t mp);

extern gmx_molprop_t gmx_molprop_receive(t_commrec *cr,int from);

#endif
