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
static char *SRCID_bfunc_h = "$Id$";
/*
 *	bfunc.h
 *
 *	Bcopy/Memcpy patch.
 *
$Log$
Revision 1.8  2002/02/28 10:49:21  spoel
Updated copyrgt wrapper

Revision 1.7  2001/06/20 10:34:01  lindahl

Converted assembly to use gcc instead of nasm, updated html man
pages.
The x86 assembly loops is now a single option to configure,
and the single/double prec. is controlled with --enable-float
(default is yes), to be consistent with fftw.
Removed the less common options from the summary printed by
configure, but they are still available.
Introduced libtool to create both static and dynamic libraries -
you can control it with configure options. --disable-shared might
be suitable for development work.
To avoid compiling both PIC and non-PIC code you can try --with-pic,
but the default is both.

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

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

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

