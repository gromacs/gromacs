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
 * Great Red Oystrich Makes All Chemists Sane
 */

#ifndef	_maths_h
#define	_maths_h

#ifdef HAVE_IDENT
#ident	"@(#) maths.h 1.11 11/24/92"
#endif /* HAVE_IDENT */

#include <math.h>
#include "typedefs.h"

#ifndef M_PI
#define	M_PI		3.14159265358979323846
#endif

#ifndef M_PI_2
#define	M_PI_2		1.57079632679489661923
#endif

#ifndef M_2PI
#define	M_2PI		6.28318530718
#endif

#ifndef M_SQRT2
#define M_SQRT2 sqrt(2.0)
#endif

#ifdef CPLUSPLUS
extern "C" {
#endif

extern	int		gmx_nint(real a);
extern	real		sqr(real x);
extern  real            sign(real x,real y);

#ifdef CPLUSPLUS
}
#endif

#endif	/* _maths_h */
