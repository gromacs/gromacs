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
 * Getting the Right Output Means no Artefacts in Calculating Stuff
 */

#ifndef _viewit_h
#define _viewit_h

static char *SRCID_viewit_h = "$Id$";
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#ifdef HAVE_IDENT
#ident "$Id$"
#endif /* HAVE_IDENT */

#ifdef CPLUSPLUS
extern "C" {
#endif

#include "typedefs.h"

extern void do_view(char *fn, char *opts);
/* forks off appropriate command to view file.
 * currently eps, xpm, xvg and pdb are supported 
 * defaults are provided, can be overriden with environment vars 
 */

extern void view_all(int nf, t_filenm fnm[]);
/* calls do_view for all viewable output files in fnm[] */
 
#ifdef CPLUSPLUS
}
#endif

#endif	/* _maths_h */
