/*
 * $Id$
 * 
 *       This source code is part of
 * 
 *        G   R   O   M   A   C   S
 * 
 * GROningen MAchine for Chemical Simulations
 * 
 *               VERSION 2.0
 * 
 * Copyright (c) 1991-1999
 * BIOSON Research Institute, Dept. of Biophysical Chemistry
 * University of Groningen, The Netherlands
 * 
 * Please refer to:
 * GROMACS: A message-passing parallel molecular dynamics implementation
 * H.J.C. Berendsen, D. van der Spoel and R. van Drunen
 * Comp. Phys. Comm. 91, 43-56 (1995)
 * 
 * Also check out our WWW page:
 * http://md.chem.rug.nl/~gmx
 * or e-mail to:
 * gromacs@chem.rug.nl
 * 
 * And Hey:
 * Good ROcking Metal Altar for Chronical Sinners
 */

#ifndef _macros_h
#define _macros_h

static char *SRCID_macros_h = "$Id$";

#ifdef HAVE_IDENT
#ident	"@(#) macros.h 1.8 11/23/92"
#endif /* HAVE_IDENT */

#ifdef CPLUSPLUS
extern "C" { 
#endif


#include "typedefs.h" /* for real definition only */

/* 
 * With the macros below you don't
 * have to use an index if you don't wan't to. You can eg. use
 * angle.C0[23] instead if angle.c[0][23].
 * In a similar fashion, you can use angle.AI[3] instead of
 * angle.a[0][3]
 */
#define AI 	a[0]
#define AJ 	a[1]
#define AK 	a[2]
#define AL 	a[3]
#define AM      a[4]
#define C0 	c[0]
#define C1 	c[1]
#define C2 	c[2]
#define C3 	c[3]
#define C4 	c[4]
#define C5 	c[5]

#define min(a,b) (((a) < (b)) ? (a) : (b) )
#define max(a,b) (((a) > (b)) ? (a) : (b) )
#define even(a) ( ( (a+1) / 2) == (a / 2) )

/* This macro calculates the size of a array */
#define asize(a) (sizeof(a)/sizeof((a)[0]))

extern real ZERO;
extern real THIRD;
extern real HALF;
extern real ONE;
extern real TWO;
extern real THREE;
extern real SIX;
extern real TEN;
extern real TWELVE;

#ifdef CPLUSPLUS
}
#endif


#endif	/* _macros_h */

