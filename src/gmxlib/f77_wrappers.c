/*
 * $Id$
 * 
 *                This source code is part of
 * 
 *                 G   R   O   M   A   C   S
 * 
 *          GROningen MAchine for Chemical Simulations
 * 
 *                        VERSION 3.1
 * Copyright (c) 1991-2001, University of Groningen, The Netherlands
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
 * Great Red Owns Many ACres of Sand 
 */
static char *SRCID_f77_wrappers_c = "$Id$";
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <math.h>
#include <typedefs.h>
#include <callf77.h>

/* These are wrapper routines for f77 to use the C
 * math library routines.
 */

#ifdef USE_FORTRAN
real F77_FUNC(cerfc,CERFC)(real *x)
{
	return erfc(*x);
}

real F77_FUNC(cpow,CPOW)(real *x,real *y)
{
	return pow(*x,*y);
}
#endif
