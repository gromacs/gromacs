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
 * Gromacs Runs On Most of All Computer Systems
 */

#ifndef _rwtop_h
#define _rwtop_h

static char *SRCID_rwtop_h = "$Id$";

#ifdef HAVE_IDENT
#ident	"@(#) rwtop.h 1.5 12/16/92"
#endif /* HAVE_IDENT */

#include <stdio.h>
#include "typedefs.h"

/*
 * This module handles topolgy manipulation, including (binary) writing
 * to and reading from a file, freeing the allocated space and textual 
 * representation of the complete topology structure.
 */

extern long wr_top(FILE *fp,t_topology *top);
/*
 * Writes the topology to the file, specified by fp. The function 
 * returns the number of bytes written. The topology is not modified!
 */

extern long rd_top(FILE *fp,t_topology *top);
/*
 * Reads the topology from the file, specified by fp. This will 
 * include allocating the needed space. The function returns the 
 * number of bytes read.
 */

extern void rm_top(t_topology *top);
/*
 * Frees the space allocated by the topology. This is only
 * guaranteed to work when the same allocation strategy is used as
 * in rd_top().
 */

extern void pr_energies(FILE *fp,int indent,char *title,t_energy *e,int n);
/*
 * This routine prints out a (human) readable representation of
 * an array of energy structs to the file fp. Ident specifies the
 * number of spaces the text should be indented. Title is used to
 * print a header text.
 */

extern void pr_inputrec(FILE *fp,int indent,char *title,t_inputrec *ir);
/*
 * This routine prints out a (human) readable representation of
 * an input record to the file fp. Ident specifies the number of spaces
 * the text should be indented. Title is used to print a header text.
 */
 
#endif	/* _rwtop_h */
