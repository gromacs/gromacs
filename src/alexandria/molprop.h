/*
 * $Id: molprop.h,v 1.16 2009/05/29 15:01:18 spoel Exp $
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

#ifndef _molprop_h
#define _molprop_h

typedef struct gmx_molprop *gmx_molprop_t;

enum { eMOLPROP_Exp, eMOLPROP_Calc, eMOLPROP_Any, eMOLPROP_NR };

enum { empPOTENTIAL, empDIPOLE, empQUADRUPOLE, empPOLARIZABILITY, 
       empENERGY, empNR };

#ifdef __cplusplus
extern "C" {
#endif

extern const char *emp_name[empNR];

extern gmx_molprop_t gmx_molprop_init();

extern void gmx_molprop_delete(gmx_molprop_t mpt);

extern void gmx_molprop_check_consistency(gmx_molprop_t mpt);

extern void gmx_molprop_set_mass(gmx_molprop_t mpt,double mass);

extern double gmx_molprop_get_mass(gmx_molprop_t mpt);

extern void gmx_molprop_set_charge(gmx_molprop_t mpt,int charge);

extern int gmx_molprop_get_charge(gmx_molprop_t mpt);

extern void gmx_molprop_set_multiplicity(gmx_molprop_t mpt,int multiplicity);

extern int gmx_molprop_get_multiplicity(gmx_molprop_t mpt);

extern void gmx_molprop_set_formula(gmx_molprop_t mpt,char *formula);

extern char *gmx_molprop_get_formula(gmx_molprop_t mpt);

extern void gmx_molprop_set_molname(gmx_molprop_t mpt,char *molname);

extern char *gmx_molprop_get_molname(gmx_molprop_t mpt);

/* Return IUPAC name or, if not found, the molname */
extern void gmx_molprop_set_iupac(gmx_molprop_t mpt,char *iupac);

extern char *gmx_molprop_get_iupac(gmx_molprop_t mpt);

extern void gmx_molprop_set_cas(gmx_molprop_t mpt,char *cas);

extern char *gmx_molprop_get_cas(gmx_molprop_t mpt);

extern void gmx_molprop_set_cid(gmx_molprop_t mpt,char *cid);

extern char *gmx_molprop_get_cid(gmx_molprop_t mpt);

extern void gmx_molprop_set_inchi(gmx_molprop_t mpt,char *inchi);

extern char *gmx_molprop_get_inchi(gmx_molprop_t mpt);

extern void gmx_molprop_delete_composition(gmx_molprop_t mpt,const char *compname);

extern void gmx_molprop_add_composition(gmx_molprop_t mpt,const char *compname);

extern void gmx_molprop_reset_composition(gmx_molprop_t mpt);

extern char *gmx_molprop_get_composition(gmx_molprop_t mpt);

extern int gmx_molprop_get_natom(gmx_molprop_t mpt);
  
extern void gmx_molprop_add_composition_atom(gmx_molprop_t mpt,const char *compname,
					     char *atomname,int natom);
					     
extern void gmx_molprop_del_composition_atom(gmx_molprop_t mpt,char *compname,
					     char *atomname);

extern void gmx_molprop_replace_composition_atom(gmx_molprop_t mpt,char *compname,char *oldatom,char *newatom);

extern int gmx_molprop_get_composition_atom(gmx_molprop_t mpt,char *compname,
					    char **atomname,int *natom);

extern int gmx_molprop_count_composition_atoms(gmx_molprop_t mpt,
					       char *compname,char *atom);
				
/* Returns an internal reference that can be used when adding data values */
extern void gmx_molprop_add_experiment(gmx_molprop_t mpt,const char *reference,
				       const char *conformation,int *expref);

extern int gmx_molprop_get_experiment(gmx_molprop_t mpt,char **reference,
				      char **conformation,int *expref);

extern void gmx_molprop_reset_experiment(gmx_molprop_t mpt);

/* Polarizabilities. Ref is either exp or calc ref. from add_exper,
   add_calc, get_exper or get_calc. */	    
extern void gmx_molprop_add_polar(gmx_molprop_t mpt,int ref,
				  const char *type,const char *unit,
				  double xx,double yy,double zz,
				  double aver,double error);
				       
extern int gmx_molprop_get_polar(gmx_molprop_t mpt,int ref,
				 char **type,char **unit,
				 double *xx,double *yy,double *zz,
				 double *aver,double *error);

/* Energies. Ref is either exp or calc ref. from add_exper,
   add_calc, get_exper or get_calc. */	    
extern void gmx_molprop_add_energy(gmx_molprop_t mpt,int ref,
				   char *type,char *unit,
				   double value,double error);
				       
extern int gmx_molprop_get_energy(gmx_molprop_t mpt,int ref,
				  char **type,char **unit,
				  double *value,double *error);

/* Potential. Ref is a calc ref. from add_calc, or get_calc. */	    
extern void gmx_molprop_add_potential(gmx_molprop_t mpt,int ref,
				      char *xyz_unit,char *V_unit,
				      int espid,
				      double x,double y,double z,
				      double V);
				       
extern int gmx_molprop_get_potential(gmx_molprop_t mpt,int ref,
				     char **xyz_unit,char **V_unit,
				     int *espid,
				     double *x,double *y,double *z,
				     double *V);

/* Dipoles. Ref is either exp or calc ref. from add_exper,
   add_calc, get_exper or get_calc. */	    
extern void gmx_molprop_add_dipole(gmx_molprop_t mpt,int ref,
				   char *type,char *unit,
				   double x,double y,double z,
				   double aver,double error);
				       
extern int gmx_molprop_get_dipole(gmx_molprop_t mpt,int ref,
				  char **type,char **unit,
				  double *x,double *y,double *z,
				  double *aver,double *error);

/* Quadrupoles. Ref is either exp or calc ref. from add_exper,
   add_calc, get_exper or get_calc. */	    
extern void gmx_molprop_add_quadrupole(gmx_molprop_t mpt,int ref,
				       char *type,char *unit,
				       double xx,double yy,double zz,
				       double xy,double xz,double yz);
				       
extern int gmx_molprop_get_quadrupole(gmx_molprop_t mpt,int ref,
				      char **type,char **unit,
				      double *xx,double *yy,double *zz,
				      double *xy,double *xz,double *yz);

extern void gmx_molprop_add_category(gmx_molprop_t mpt,char *category);

/* Returns one category at a time. If NULL, you got them all previously. */
extern char *gmx_molprop_get_category(gmx_molprop_t mpt);

extern void gmx_molprop_reset_category(gmx_molprop_t mpt);

extern int gmx_molprop_search_category(gmx_molprop_t mpt,char *catname);

extern void gmx_molprop_reset_calculation(gmx_molprop_t mpt);

extern void gmx_molprop_reset(gmx_molprop_t mpt);

/* Returns calcref that can be used to add properties later on */
extern void gmx_molprop_add_calculation(gmx_molprop_t mpt,
					const char *program,char *method,
					char *basisset,char *reference,
					char *conformation,int *calcref);

extern int gmx_molprop_get_calculation(gmx_molprop_t mpt,char **program,char **method,
				       char **basisset,char **reference,
				       char **conformation,int *calcref);

extern int gmx_molprop_get_calc_lot(gmx_molprop_t mpt,char *lot);

/* Returns atomref that can be used to set coordinates and charges */
extern void gmx_molprop_calc_add_atom(gmx_molprop_t mpt,int calcref,
				      char *atomname,int atomid,int *atomref);

extern int gmx_molprop_calc_get_atom(gmx_molprop_t mpt,int calcref,
				     char **atomname,int *atomid,int *atomref);
				     
extern void gmx_molprop_calc_set_atomcoords(gmx_molprop_t mpt,int calcref,int atomref,
					    char *unit,double x,double y,double z);
					   
extern int gmx_molprop_calc_get_atomcoords(gmx_molprop_t mpt,int calcref,int atomref,
					   char **unit,double *x,double *y,double *z);
					   
extern void gmx_molprop_calc_set_atomcharge(gmx_molprop_t mpt,int calcref,int atomref,
					   char *type,char *unit,double q);
					   
extern int gmx_molprop_calc_get_atomcharge(gmx_molprop_t mpt,int calcref,int atomref,
					   char **type,char **unit,double *q);

/* Generic routines */					   
extern gmx_molprop_t gmx_molprop_copy(gmx_molprop_t mpt);

extern void gmx_molprop_merge(gmx_molprop_t dst,gmx_molprop_t src);

#ifdef __cplusplus
}
#endif

#endif
