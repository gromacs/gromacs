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
 * Good ROcking Metal Altar for Chronical Sinners
 */
static char *SRCID_bfunc_h = "$Id$";
/*
 *	bfunc.h
 *
 *	Bcopy/Memcpy patch.
 *
$Log$
Revision 1.6  2001/05/14 17:58:06  lindahl

Tagged files with gromacs 3.0 header

Revision 1.5  1999/11/03 12:45:47  hess
copyrgted

Revision 1.4  1998/12/10 07:43:44  spoel
Trying to get everything in synch again, Makefiles remain problematic.
For instance the shared libraries do not work any longer...

 * Revision 1.3  1997/12/23  11:52:07  anton
 * Edited by Copyright -> 2.0
 *
 * Revision 1.2  1997/11/27  16:29:42  anton
 * Edited by copyrgt -> v1.6; fixed loads of inconsistent formatting of .h files
 *
 * Revision 1.1.1.1  1997/11/03  16:08:02  spoel
 * Generated_by_makecvs
 *
 * Revision 1.1  1993/08/30  23:26:46  manchek
 * Initial revision
 *
 */


#if defined(SYSVBFUNC)
#include <memory.h>
#define BZERO(d,n)      memset(d,0,n)
#define BCMP(s,d,n)     memcmp(d,s,n)
#define BCOPY(s,d,n)    memcpy(d,s,n)

#else
#define BZERO(d,n)      bzero(d,n)
#define BCMP(s,d,n)     bcmp(s,d,n)
#define BCOPY(s,d,n)    bcopy(s,d,n)

#endif

