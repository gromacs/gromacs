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
#ifndef _shake_h
#define _shake_h

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "force.h"
#include "typedefs.h"

extern int bshakef(FILE *log,		/* Log file			*/
		   int natoms,		/* Total number of atoms	*/
		   real invmass[],	/* Atomic masses		*/
		   int nblocks,		/* The number of shake blocks	*/
		   int sblock[],        /* The shake blocks             */
		   t_idef *idef,	/* The interaction def		*/
		   t_inputrec *ir,	/* Input record		        */
		   matrix box,		/* The box			*/
		   rvec x_s[],		/* Coords before update		*/
		   rvec xp[],		/* Output coords		*/
		   t_nrnb *nrnb);        /* Performance measure          */
/* Shake all the atoms blockwise. It is assumed that all the constraints
 * in the idef->shakes field are sorted, to ascending block nr. The
 * sblock array points into the idef->shakes.iatoms field, with block 0 
 * starting
 * at sblock[0] and running to ( < ) sblock[1], block n running from 
 * sblock[n] to sblock[n+1]. Array sblock should be large enough.
 * Return 0 when OK, -1 when shake-error
 */
extern void csettle(FILE *log,
		    int nshake,		/* Number of water molecules 	*/
		    int owptr[],	/* pointer to Oxygen in b4 & after */
		    real b4[],		/* Old coordinates		*/
		    real after[],	/* New coords, to be settled	*/
		    real dOH,		/* Constraint length Ox-Hyd	*/
		    real dHH, 		/* Constraint length Hyd-Hyd	*/
		    real mO,  		/* Mass of Oxygen		*/
		    real mH  		/* Mass of Hydrogen		*/
		    );

extern void cshake(atom_id iatom[],int ncon,int *nnit,int maxnit,
		   real dist2[],real xp[],real rij[],real m2[],
		   real invmass[],real tt[],int *nerror);
/* Regular iterative shake */

/* Fortran versions of shake and settle */
DECLAREF77(fsettle) (int *nshake,int owptr[],
		     real b4[],real after[],
		     real *dOH,real *dHH,real *mO,real *mH);
		     
DECLAREF77(fshake)  (atom_id iatom[],int *ncon,int *nit,int *maxnit,
		     real dist2[],real xp[],real rij[],real m2[],
		     real invmass[],real tt[],int *error);

#endif
