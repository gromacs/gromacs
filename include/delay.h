/*
 *       $Id$
 *
 *       This source code is part of
 *
 *        G   R   O   M   A   C   S
 *
 * GROningen MAchine for Chemical Simulations
 *
 *            VERSION 1.6
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
 * Great Red Owns Many ACres of Sand 
 */

#ifndef _delay_h
#define _delay_h

static char *SRCID_delay_h = "$Id$";

#ifdef HAVE_IDENT
#ident	"@(#) delay.h 1.5 11/23/92"
#endif /* HAVE_IDENT */

extern void delay(int ms);
     /* 
      * Delay for ms milliseconds within 1 percent accurate.
      */

extern void delay01(int ms01);
     /* 
      * Delay for ms 0.01 milliseconds within 1 percent accurate.
      */

#endif	/* _delay_h */
