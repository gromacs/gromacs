/* -*- mode: c; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; c-file-style: "stroustrup"; -*-
 * $Id: poldata.h,v 1.13 2009/02/03 16:08:37 spoel Exp $
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

#ifndef _poldata_h
#define _poldata_h

/* This source code file is part of the Alexandria project */
	
enum { eqgNone, eqgYang, eqgBultinck, eqgRappe,
       eqgAXp, eqgAXs, eqgAXg, eqgESP, eqgRESP, eqgNR };

#ifdef __cplusplus
extern "C" {
#endif

typedef struct gmx_poldata *gmx_poldata_t;

extern gmx_poldata_t gmx_poldata_init();

extern int gmx_poldata_get_natypes(gmx_poldata_t pd);

extern void gmx_poldata_add_ffatype(gmx_poldata_t pd,char *elem,char *desc,
                                    char *gt_name,char *gt_type,
                                    char *miller_equiv,
                                    char *charge,char *geometry,int numbonds,
                                    char *neighbors,
                                    double polarizability,double sig_pol,
                                    char *vdwparams);

extern void gmx_poldata_set_ffatype(gmx_poldata_t pd,char *gt_type,
                                    double polarizability,double sig_pol);
		
extern void gmx_poldata_set_ffatype_ff(gmx_poldata_t pd,char *forcefield);

extern void gmx_poldata_set_ffatype_unit(gmx_poldata_t pd,char *polar_unit);

extern void gmx_poldata_set_length_unit(gmx_poldata_t pd,char *length_unit);

extern char *gmx_poldata_get_ffatype_unit(gmx_poldata_t pd);

extern char *gmx_poldata_get_ffatype_ff(gmx_poldata_t pd);

extern char *gmx_poldata_get_length_unit(gmx_poldata_t pd);

extern char *gmx_poldata_get_gt_atom(gmx_poldata_t pd,char *elem,
                                    int nbond,char *neighbors[],
                                    const char *geometry);
				    
extern char *gmx_poldata_get_geometry(gmx_poldata_t pd,char *gt_atom);

extern char *gmx_poldata_get_gt_type(gmx_poldata_t pd,char *gt_atom);

extern char *gmx_poldata_get_desc(gmx_poldata_t pd,char *gt_atom);

/* Get the charge from the gentop.dat file */
extern char *gmx_poldata_get_charge(gmx_poldata_t pd,char *gt_atom);

/* Returns name or NULL if last or not found */
extern char *gmx_poldata_get_ffatype(gmx_poldata_t pd,char *name,
                                     char **elem,char **desc,
                                     char **gt_type,char **miller_equiv,
                                     char **charge,char **geometry,
                                     int *numbonds,char **neighbors,
                                     double *polarizability,double *sig_pol,
                                     char **vdwparams);

/* Return 1 if OK, 0 if not found */				 
extern int gmx_poldata_gt_type_polarizability(gmx_poldata_t pd,char *gt_type,
                                              double *polarizability,double *sig_pol);
 
extern void gmx_poldata_add_miller(gmx_poldata_t pd,char *name,
                                   int atomnumber,
                                   double tau_ahc,double alpha_ahp,
                                   char *alexandria_equiv);

extern void gmx_poldata_set_miller_units(gmx_poldata_t pd,char *tau_unit,
                                         char *ahp_unit);

extern void gmx_poldata_get_miller_units(gmx_poldata_t pd,char **tau_unit,
                                         char **ahp_unit);

/* Returns name or NULL if last or not found */
extern char *gmx_poldata_get_miller(gmx_poldata_t pd,char *name,
                                    int *atomnumber,
                                    double *tau_ahc,double *alpha_ahp,
                                    char **alexandria_equiv);

extern char *gmx_poldata_get_miller_equiv(gmx_poldata_t pd,char *gt_name);

extern void gmx_poldata_add_bosque(gmx_poldata_t pd,char *elem,
                                   double polarizability);
				  
extern void gmx_poldata_set_bosque_unit(gmx_poldata_t pd,char *polar_unit);

extern char *gmx_poldata_get_bosque_unit(gmx_poldata_t pd);

/* Returns elem or NULL if last or not found */
extern char *gmx_poldata_get_bosque(gmx_poldata_t pd,char *name,char *elem,
                                    double *polarizability);
				  
extern void gmx_poldata_add_gt_bond(gmx_poldata_t pd,char *atom1,char *atom2,
                                    double length,double sigma,double bondorder,
                                    char *params);
				   
extern int gmx_poldata_get_gt_bond(gmx_poldata_t pd,char **atom1,char **atom2,
                                   double *length,double *sigma,double *bondorder,
                                   char **params);

extern int gmx_poldata_search_gt_bond(gmx_poldata_t pd,char *atom1,char *atom2,
                                      double *length,double *sigma,char **params);

/* Returns the bondorder */
extern double gmx_poldata_get_bondorder(gmx_poldata_t pd,char *atom1,char *atom2,
                                        double distance,double toler);
				     
extern void gmx_poldata_add_gt_angle(gmx_poldata_t pd,char *atom1,char *atom2,
                                     char *atom3,double angle,double sigma,char *params);
				   
extern int gmx_poldata_get_gt_angle(gmx_poldata_t pd,char **atom1,char **atom2,
                                    char **atom3,double *angle,double *sigma,
                                    char **params);

extern int gmx_poldata_search_gt_angle(gmx_poldata_t pd,char *atom1,char *atom2,
                                       char *atom3,double *angle,double *sigma,
                                       char **params);

extern void gmx_poldata_set_angle_unit(gmx_poldata_t pd,char *angle_unit);

extern char *gmx_poldata_get_angle_unit(gmx_poldata_t pd);

extern void gmx_poldata_add_gt_dihedral(gmx_poldata_t pd,char *atom1,char *atom2,
                                        char *atom3,char *atom4,
                                        double dihedral,double sigma,char *params);
				   
extern int gmx_poldata_get_gt_dihedral(gmx_poldata_t pd,char **atom1,char **atom2,
                                       char **atom3,char **atom4,
                                       double *dihedral,double *sigma,char **params);

extern int gmx_poldata_search_gt_dihedral(gmx_poldata_t pd,char *atom1,char *atom2,
                                          char *atom3,char *atom4,
                                          double *dihedral,double *sigma,char **params);

extern void gmx_poldata_set_dihedral_unit(gmx_poldata_t pd,char *dihedral_unit);

extern char *gmx_poldata_get_dihedral_unit(gmx_poldata_t pd);

extern void gmx_poldata_add_symcharges(gmx_poldata_t pd,char *central,
                                       char *attached,int numattach);

extern int gmx_poldata_get_symcharges(gmx_poldata_t pd,char **central,
                                      char **attached,int *numattach);
				      
extern int gmx_poldata_search_symcharges(gmx_poldata_t pd,char *central,
                                         char *attached,int numattach);

extern int name2eemtype(const char *name);

extern const char *get_eemtype_name(int eem);

extern char *gmx_poldata_get_eemref(gmx_poldata_t pd,int eqg_model);

extern int gmx_poldata_get_numprops(gmx_poldata_t pd,int eqg_model);

extern int gmx_poldata_have_pol_support(gmx_poldata_t pd,char *gt_type);

extern int gmx_poldata_have_eem_support(gmx_poldata_t pd,int eqg_model,char *name,
                                        gmx_bool bAllowZeroParameters);

extern double gmx_poldata_get_j00(gmx_poldata_t pd,int eqg_model,char *name);

extern int gmx_poldata_get_nzeta(gmx_poldata_t pd,int eqg_model,char *name);

extern double gmx_poldata_get_zeta(gmx_poldata_t pd,int eqg_model,char *name,int zz);

extern char *gmx_poldata_get_qstr(gmx_poldata_t pd,int eqg_model,char *name);

extern char *gmx_poldata_get_rowstr(gmx_poldata_t pd,int eqg_model,char *name);

extern double gmx_poldata_get_q(gmx_poldata_t pd,int eqg_model,char *name,int zz);

extern int gmx_poldata_get_row(gmx_poldata_t pd,int eqg_model,char *name,int zz);

extern double gmx_poldata_get_chi0(gmx_poldata_t pd,int eqg_model,char *name);

extern char *gmx_poldata_get_opts(gmx_poldata_t pd,int eqg_model,char *name);

extern void gmx_poldata_set_eemprops(gmx_poldata_t pd,
                                     int eqg_model,char *name,
                                     double J0,double chi0,
                                     char *zeta,char *q,char *row);

extern int gmx_poldata_get_eemprops(gmx_poldata_t pd,
                                    int *eqg_model,char **name,
                                    double *J0,double *chi0,
                                    char **zeta,char **q,char **row);
				    
extern void gmx_poldata_set_epref(gmx_poldata_t pd,int eqg_model,char *epref);
				    
extern char *gmx_poldata_get_epref(gmx_poldata_t pd,int eqg_model);

extern int gmx_poldata_list_epref(gmx_poldata_t pd,int *eqg_model,char **epref);

extern void gmx_poldata_comm_eemprops(gmx_poldata_t pd,t_commrec *cr);

#ifdef __cplusplus
}
#endif

#endif
