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
 * Good gRace! Old Maple Actually Chews Slate
 */

#ifndef _readir_h
#define _readir_h

static char *SRCID_readir_h = "$Id$";

#ifdef HAVE_IDENT
#ident	"@(#) readir.h 1.20 9/30/97"
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
  char *watertype;
  bool bGenVel;
  bool bGenPairs;
  real tempi;
  int  seed;
} t_gromppopts;


extern void init_ir(t_inputrec *ir, t_gromppopts *opts);
/* Initiate stuff */

extern void check_ir(t_inputrec *ir, t_gromppopts *opts);
/* Validate inputrec data */

extern void get_ir(char *mdparin,char *mdparout,
		   t_inputrec *ir,t_gromppopts *opts);
/* Read the input file, and retrieve data for inputrec.
 * More data are read, but the are only evaluated when the next
 * function is called. Also prints the input file back to mdparout
 */
 
extern void do_index(char *ndx,
		     t_symtab   *symtab,
		     t_atoms    *atoms,bool bVerbose,
		     t_inputrec *ir,t_idef *idef,int *forward);
/* Read the index file and assign grp numbers to atoms */

#endif	/* _readir_h */
