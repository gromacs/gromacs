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

#include "sysstuff.h"
#include "typedefs.h"
#include "xdrf.h"

typedef struct {
  real t;		/* Time frame				*/
  int  step;		/* MD step				*/
  int  nre;		/* Number of energies			*/
  int  ndisre;		/* Number of disre blocks		*/
  int  nuser;		/* User definable number		*/
  int  e_size;		/* Size (in bytes) of energies		*/
  int  d_size;		/* Size (in bytes) of disre blocks	*/
  int  u_size;		/* Size (in bytes) of user blocks	*/
} t_eheader;

extern void wr_ener_nms(FILE *out,int nre,char *nms[]);
/* Write energy names to file */

extern void wr_ener(FILE *out,real t,int step,int nre,t_energy ener[],
		    t_drblock *drblock);
/* Write energy terms and distance restraint data to file. */

extern void rd_ener_nms(FILE *in,int *nre,char ***nm);
/* Read energy names and number. */

extern bool rd_ener(FILE *in,real *t,int *step,
		    t_energy ener[],t_drblock *drblock);
/* Read energies and distance restraint data. 
 * If memory for distance restraints has not been allocated
 * (drblock->ndr == 0) the routine will do so, and set
 * drblock->ndr.
 * Return FALSE on eof, TRUE on succes. 
 */

/************************************************
 *
 *   X D R format energies
 *
 ************************************************/

extern void edr_nms(XDR *xdr,int *nre,char ***nms);
/* Read/Write energy names to file */

extern bool edr_io(XDR *xdr,real *t,int *step,int *nre,t_energy ener[],
		   t_drblock *drblock);
/* Read/Write energies to portable file format, return TRUE if succesfull
 * FALSE otherwise.
 */

#ifdef CPLUSPLUS
}
#endif

#endif	/* _enerio_h */
