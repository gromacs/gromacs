/*
 *       $Id$
 *
 *       This source code is part of
 *
 *        G   R   O   M   A   C   S
 *
 * GROningen MAchine for Chemical Simulations
 *
 *            VERSION 2.0
 * 
 * Copyright (c) 1991-1997
 * BIOSON Research Institute, Dept. of Biophysical Chemistry
 * University of Groningen, The Netherlands
 * 
 * Please refer to:
 * GROMACS: A message-passing parallel molecular dynamics implementation
 * H.J.C. Berendsen, D. van der Spoel and R. van Drunen
 * Comp. Phys. Comm. 91, 43-56 (1995)
 *
 * Also check out our WWW page:
 * http://rugmd0.chem.rug.nl/~gmx
 * or e-mail to:
 * gromacs@chem.rug.nl
 *
 * And Hey:
 * GRowing Old MAkes el Chrono Sweat
 */
static char *SRCID_bfunc_h = "$Id$";

/*
 *	bfunc.h
 *
 *	Bcopy/Memcpy patch.
 *
$Log$
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

