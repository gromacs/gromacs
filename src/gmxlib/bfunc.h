/*
 *       @(#) copyrgt.c 1.12 9/30/97
 *
 *       This source code is part of
 *
 *        G   R   O   M   A   C   S
 *
 * GROningen MAchine for Chemical Simulations
 *
 *            VERSION 2.0b
 * 
 * Copyright (c) 1990-1997,
 * BIOSON Research Institute, Dept. of Biophysical Chemistry,
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
 * GROwing Monsters And Cloning Shrimps
 */

/*
 *	bfunc.h
 *
 *	Bcopy/Memcpy patch.
 *
$Log$
Revision 1.1  1997/11/03 16:08:02  spoel
Initial revision

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

