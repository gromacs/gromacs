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
 * Giant Rising Ordinary Mutants for A Clerical Setup
 */

#ifndef	_runner_h
#define	_runner_h

#ifdef HAVE_IDENT
#ident	"@(#) runner.h 1.7 02 Nov 1995"
#endif /* HAVE_IDENT */
#include "typedefs.h"
#include "network.h"
#include "filenm.h"
	
extern void mdrunner(t_commrec *cr,int nfile,t_filenm fnm[],bool bVerbose,
		     bool bCompact,int nDlb);

extern void nmrunner(t_commrec *cr,int nfile,t_filenm fnm[],bool bVerbose,
		     bool bCompact,int nDlb);
#endif	/* _runner_h */
