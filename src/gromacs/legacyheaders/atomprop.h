/*
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
 * Gromacs Runs On Most of All Computer Systems
 */

#ifndef _atomprop_h
#define _atomprop_h

#ifdef __cplusplus
extern "C" {
#endif

#include "index.h"

/* Abstract type for the atom property database */
typedef struct gmx_atomprop *gmx_atomprop_t;

enum { epropMass, epropVDW, epropDGsol, epropElectroneg, epropElement, 
       epropNR };

gmx_atomprop_t gmx_atomprop_init(void);
/* Initializes and returns the atom properties struct */

void gmx_atomprop_destroy(gmx_atomprop_t aps);
/* Get rid of memory after use */

char *gmx_atomprop_element(gmx_atomprop_t aps,int atomnumber);

int gmx_atomprop_atomnumber(gmx_atomprop_t aps,const char *element);

gmx_bool gmx_atomprop_query(gmx_atomprop_t aps,
                        int eprop,const char *resnm,const char *atomnm,
                        real *value);
/* Extract a value from the database. Returns TRUE on succes,
 * FALSE otherwise. In the latter case, value is a deafult value.
 * The first time this function is called for this property
 * the database will be read.
 */

#ifdef __cplusplus
}
#endif


#endif
