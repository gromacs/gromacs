/* 
 * $Id: gauss_io.h,v 1.8 2009/02/02 21:11:11 spoel Exp $
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

#ifndef _gauss_io_h
#define _gauss_io_h

#include "molprop.hpp"

typedef struct gau_atomprop *gau_atomprop_t;

/* read composite method atom data */
gau_atomprop_t read_gauss_data(void);

void done_gauss_data(gau_atomprop_t gaps);

int gau_atomprop_get_value(gau_atomprop_t gaps,const char *element,
                           const char *method,
                           char *desc,double temp,double *value);

void ReadGauss(const char *g98,
               alexandria::MolProp& mpt,
               gmx_bool bBabel,
               gmx_atomprop_t aps,gmx_poldata_t pd,
               char *molnm,char *iupac,char *conformation,
               char *basisset,gau_atomprop_t gaps,
               real th_toler,real ph_toler,
               int maxpot,gmx_bool bVerbose);

void translate_atomtypes(t_atoms *atoms,t_symtab *tab,const char *forcefield);

#endif
