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
 * GROningen MAchine for Chemical Simulation
 */

#ifndef	_enerio_h
#define	_enerio_h

#ifdef HAVE_IDENT
#ident	"@(#) enerio.h 1.13 2/2/97"
#endif /* HAVE_IDENT */

#ifdef CPLUSPLUS
external "C" {
#endif

  /**************************************************************
   *
   * The routines in the corresponding c-file enxio.c
   * are based on the lower level routines in gmxfio.c
   * The integer file pointer returned from open_enx
   * can also be used with the routines in gmxfio.h
   *
   **************************************************************/

#include "sysstuff.h"
#include "typedefs.h"
#include "xdrf.h"
  
  /* New energy reading and writing interface */
  extern int open_enx(char *fn,char *mode);
  
  extern void close_enx(int fp_ene);
  
  extern void do_enxnms(int fp_ene,int *nre,char ***nms);
  
  extern bool do_enx(int fp_ene,real *t,int *step,int *nre,
		     t_energy ener[],t_drblock *drblock);
  
  


 
#ifdef CPLUSPLUS
}
#endif

#endif	/* _enerio_h */
