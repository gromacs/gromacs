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

#ifndef _molprop_util_hpp
#define _molprop_util_hpp

#include "grompp.h"
#include "atomprop.h"
#include "molprop.hpp"
#include "poldata.h"
#include "molselect.h"

enum { iqmExp, iqmBoth, iqmQM, iqmNR };

/* Returns 1 when OK, 0 otherwise. Level of theory and/or conformation
   and/or type (elec, ESP, RESP) may be set. */
//extern int get_val(alexandria::ExperimentIterator mp,char *type,int emp,
//                   double *value,double *error,double vec[3],
//                   tensor quadrupole);
		   
extern int mp_get_prop(alexandria::MolProp mp,int eGP,int iQM,char *lot,
                       char *conf,char *type,double *value);

/* mylot returns the used level of theory. mylot and ref may be NULL */
extern int mp_get_prop_ref(alexandria::MolProp mp,int emp,int iQM,char *lot,
                           const char *conf,const char *type,
                           double *value,double *err,
                           char **ref,char **mylot,
                           double vec[3],tensor quadrupole);

extern int gmx_molprop_get_calc_lot(alexandria::MolProp mp,char *lot);

extern void generate_formula(std::vector<alexandria::MolProp> mp,gmx_atomprop_t ap);

/* Returns TRUE if the present molecule has the indicated composition */
extern gmx_bool gmx_molprop_support_composition(alexandria::MolProp mp,char *composition);

void generate_composition(std::vector<alexandria::MolProp> mp,
                          gmx_poldata_t pd,
                          gmx_atomprop_t ap,gmx_bool bForceGenComp,
                          double th_toler,double phi_toler);

extern alexandria::MolProp atoms_2_molprop(char *molname,t_atoms*atoms,char **smnames,
                                           gmx_atomprop_t ap,gmx_poldata_t pd,gmx_bool bForce,
                                           double th_toler,double ph_toler);

/* Return number of atoms, 0 means failure */
extern int molprop_2_topology(alexandria::MolProp mp,gmx_atomprop_t ap,
                              gmx_poldata_t pd,
                              t_symtab *tab,const char *lot,
                              t_topology *top,const char *q_algorithm,
                              rvec **x,t_params plist[F_NRE],
                              int nexcl,t_excls **excls);

extern void merge_doubles(std::vector<alexandria::MolProp> mp,
                          char *doubles,gmx_bool bForceMerge);
				     				     
extern std::vector<alexandria::MolProp> merge_xml(int nfile,char **infiles,char *outf,
                                                  char *sorted,char *doubles,
                                                  gmx_atomprop_t ap,gmx_poldata_t pd,
                                                  gmx_bool bForceMerge,gmx_bool bForceGenComp,
                                                  double th_toler,double ph_toler);

enum { empSORT_Molname, empSORT_Formula, empSORT_Composition, empSORT_Selection, empSORT_NR };
	       
extern void gmx_molprop_sort(std::vector<alexandria::MolProp> mp,
                             int alg,gmx_atomprop_t ap,gmx_molselect_t gms);
/* gms may be NULL if the alg is not empSORT_Selection */

typedef struct {
    int  n;
    int  *count;
    char **method,**basis,**type,**lot;
    int  nconf;
    char **conf;
} t_qmcount;

/* Check the available molprops to see what kind of calculations are stored in there */
extern t_qmcount *find_calculations(std::vector<alexandria::MolProp> mp,
                                    int emp,const char *fc_str);
				
extern void MolPropSort(std::vector<alexandria::MolProp> mp,
                        int alg,gmx_atomprop_t apt,
                        gmx_molselect_t gms);

#endif
