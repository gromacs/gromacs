/*
 * $Id$
 * 
 *                This source code is part of
 * 
 *                 G   R   O   M   A   C   S
 * 
 *          GROningen MAchine for Chemical Simulations
 * 
 *                        VERSION 3.0
 * 
 * Copyright (c) 1991-2001
 * BIOSON Research Institute, Dept. of Biophysical Chemistry
 * University of Groningen, The Netherlands
 * 
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
 * Do check out http://www.gromacs.org , or mail us at gromacs@gromacs.org .
 * 
 * And Hey:
 * Giving Russians Opium May Alter Current Situation
 */

#ifndef _ewald_util_h
#define _ewald_util_h

static char *SRCID_ewald_util_h = "$Id$";
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <math.h>
#include "typedefs.h"
#include "gmxcomplex.h"


extern real ewald_LRcorrection(FILE *fp,t_nsborder *nsb,
			       t_commrec *cr,t_forcerec *fr,
			       real charge[],t_block *excl,rvec x[],
			       matrix box,rvec mu_tot, real qsum,
			       real surface_eps,matrix lrvir);
/* Calculate the Long range correction to ewald, due to 
 * 1-4 interactions, surface dipole term and charge terms
 */

extern real calc_ewaldcoeff(real rc,real dtol);
/* Determines the Ewald parameter, both for Ewald and PME */
#endif
