/*
 * $Id: molprop_tables.h,v 1.10 2009/05/29 15:01:18 spoel Exp $
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

#ifndef _molprop_tables_h
#define _molprop_tables_h

#include "grompp.h"
#include "atomprop.h"
#include "molprop.h"
#include "poldata.h"
#include "molprop_util.h"
#include "molselect.h"

extern char *ftoa(double f);

extern char *itoa(int f);

/* Outlier indicates a level (in units of sigma, one standard
   deviation). Calculations that deviate more than this level from the
   experiment are not taken into account when computing
   statistics. Moreover, the outliers are printed to the standard
   error. If outlier is 0, no action is take. */
extern void gmx_molprop_stats_table(FILE *fp,int eprop,
				    int nmol,gmx_molprop_t mp[],int ntot,
				    t_qmcount *qmc,int iQM,char *lot,
				    double outlier,gmx_molselect_t gms,int ims);

extern void gmx_molprop_composition_table(FILE *fp,int nmol,gmx_molprop_t mp[],
					  gmx_molselect_t gms,int ims);

extern void gmx_molprop_category_table(FILE *fp,int np,gmx_molprop_t mp[],
                                       gmx_molselect_t gms,int ims);

extern void gmx_molprop_prop_table(FILE *fp,int eprop,real rtoler,real atoler,
				   int np,gmx_molprop_t mp[],int bDS,t_qmcount *qmc,
				   gmx_bool bPrintAll,gmx_molselect_t gms,int ims);
				   
/* Calling this with bPolar TRUE will print an atomtype table with polarizability information.
 * With bPolar FALSE it will print the same table with EEM parameters.
 */
extern void gmx_molprop_atomtype_table(FILE *fp,gmx_bool bPolar,
				       int npd,gmx_poldata_t pd[],
				       gmx_poldata_t pd_aver, /* Output! */
				       int nmol,gmx_molprop_t mp[],
				       int iQM,char *lot,output_env_t oenv,
				       const char *histo);

#endif
