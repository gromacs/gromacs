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

#ifndef _mvdata_h
#define _mvdata_h

static char *SRCID_mvdata_h = "$Id$";

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

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
