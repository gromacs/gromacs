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
 * GROningen MAchine for Chemical Simulation
 */

#ifndef _comlib_h
#define _comlib_h

static char *SRCID_comlib_h = "$Id$";

#ifdef HAVE_IDENT
#ident	"@(#) comlib.h 1.3 11/23/92"
#endif /* HAVE_IDENT */

extern void put_serverbyte(unsigned char data);
extern unsigned char get_serverbyte();
extern void get_serverdata(void *data,int size);
extern void put_serverdata(void *data,int size);

#endif	/* _comlib_h */
