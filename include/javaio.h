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

#ifndef _javaio_h
#define _javaio_h

static char *SRCID_javaio_h = "$Id$";

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "typedefs.h"

void write_java(FILE *out,
                char *program,
                int nldesc,char *desc[],
                int nfile,t_filenm fnm[],
                int npargs,t_pargs pa[],
                int nbug,char *bugs[]);

#endif  /* _javaio_h */
 
