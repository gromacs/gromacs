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

#ifndef _rwtop_h
#define _rwtop_h

static char *SRCID_rwtop_h = "$Id$";
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

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
