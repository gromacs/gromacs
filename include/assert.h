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
 * Gnomes, ROck Monsters And Chili Sauce
 */

#ifndef _assert_h
#define _assert_h

static char *SRCID_assert_h = "$Id$";

#ifdef HAVE_IDENT
#ident	"@(#) assert.h 1.12 11/23/92"
#endif /* HAVE_IDENT */

#include <ctype.h>
#include "sysstuff.h"

#ifdef NDEBUG
#define assert(EX)	((void)0)
#else

#define assert(EX)  if ((EX) != 0); \
  else  \
    do { \
      fprintf(stderr,"Assertion failed for \"%s\" in file %s," \
      	      " line %d\n",#EX, __FILE__, __LINE__); \
      fprintf(stderr,"dump core ? (y/n):"); fflush(stderr); \
      if (toupper(getc(stdin))=='Y') \
         (void) abort(); else exit(-1); \
    } while (0)
#endif	/* NDEBUG */
#endif	/* _assert_h */
