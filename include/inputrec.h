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
 * Great Red Owns Many ACres of Sand 
 */

#ifndef	_inputrec_h
#define	_inputrec_h

#ifdef HAVE_IDENT
#ident	"@(#) inputrec.h 1.6 11/23/92"
#endif /* HAVE_IDENT */

#include "typedefs.h"
#include "string2.h"

enum { eshNONE, eshHBONDS, eshALLBONDS, eshHANGLES, eshALLANGLES, eshNR };

typedef struct {
  int warnings;
  int nshake;
  int nprocs;
  int splitalg;
  char *title;
  char *cpp;
  char *include;
  char *define;
  bool bGen;
  real tempi;
  real epsr;
  int  seed;
} t_gromppopts;

extern void read_ir(char *mdparin,t_inputrec *ir,t_gromppopts *opts);

extern void write_ir(char *outfile,t_inputrec *ir,t_gromppopts *opts);

extern void init_ir(t_inputrec *ir,t_gromppopts *opts);

extern void check_ir(t_inputrec *ir,t_gromppopts *opts);

#endif	/* _inputrec_h */
