/* -*- mode: c; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; c-file-style: "stroustrup"; -*-
 *
 * 
 *                This source code is part of
 * 
 *                 G   R   O   M   A   C   S
 * 
 *          GROningen MAchine for Chemical Simulations
 * 
 * Written by David van der Spoel, Erik Lindahl, Berk Hess, and others.
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2012, The GROMACS development team,
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

#ifndef _tables_h
#define _tables_h


#ifdef __cplusplus
extern "C" {
#endif


void table_spline3_fill_ewald_lr(real *tabf,real *tabv,
                                 int ntab,int tableformat,
                                 real dr,real beta);
/* Fill table tabf of size ntab with spacing dr with the ewald long-range
 * (mesh) force and with tableformatF and tabv!=NULL, fill tabv energy.
 * With tableformatFDV0 the size of the tabf array should be ntab*4, tabv=NULL.
 * This function interpolates the Ewald mesh potential contribution
 * with coefficient beta using a quadratic spline.
 * The force can then be interpolated linearly.
 */

real ewald_spline3_table_scale(real ewaldcoeff,real rc);
/* Return the scaling for the Ewald quadratic spline tables. */


#ifdef __cplusplus
}
#endif

#endif	/* _tables\_h */
