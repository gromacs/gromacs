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
 * Gyas ROwers Mature At Cryogenic Speed
 */

#ifndef	_mvdata_h
#define	_mvdata_h

#ifdef HAVE_IDENT
#ident	"@(#) mvdata.h 1.6 11/23/92"
#endif /* HAVE_IDENT */

#include "typedefs.h"
#include "nsb.h"

extern void ld_data(int left,int right,t_parm *parm,t_nsborder *nsb,
		    t_topology *top,rvec **x,rvec **v);

extern void mv_data(int left,int right,t_parm *parm,t_nsborder *nsb,
		    t_topology *top,rvec x[],rvec v[]);

extern void move_cgcm(FILE *log,t_commrec *cr,rvec cg_cm[],int nload[]);
		     
extern void move_rvecs(FILE *log,bool bForward,bool bSum,
		       int left,int right,rvec vecs[],rvec buf[],
		       int shift,t_nsborder *nsb,t_nrnb *nrnb);

extern void move_x(FILE *log,
		   int left,int right,rvec x[],t_nsborder *nsb,t_nrnb *nrnb);
		    
extern void move_f(FILE *log,
		   int left,int right,rvec f[],rvec fadd[],
		   t_nsborder *nsb,t_nrnb *nrnb);
		    
#endif	/* _mvdata_h */
