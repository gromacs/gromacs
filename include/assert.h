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
 * GROningen Mixture of Alchemy and Childrens' Stories
 */

#ifndef _assert_h
#define _assert_h

static char *SRCID_assert_h = "$Id$";

#ifdef HAVE_IDENT
#ident	"@(#) assert.h 1.12 11/23/92"
#endif /* HAVE_IDENT */

#include <ctype.h>
#include "sysstuff.h"

#ifdef assert
#undef assert
#endif

#define assert(EXPRESSION)  \
  if (!(EXPRESSION)) { \
    fprintf(stderr,"Assertion failed for \"%s\" in file %s, " \
	    "line %d\ndump core ? (y/n):",#EXPRESSION, __FILE__, __LINE__); \
    fflush(stderr); \
    if (toupper(getc(stdin))=='Y') \
      (void) abort(); else exit(-1); \
  } else 

#endif	/* _assert_h */
