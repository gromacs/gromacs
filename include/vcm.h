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
 * Green Red Orange Magenta Azure Cyan Skyblue
 */

#ifndef _vcm_h
#define _vcm_h

static char *SRCID_vcm_h = "$Id$";

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#ifdef HAVE_IDENT
#ident	"@(#) vcm.h 1.9 9/29/97"
#endif /* HAVE_IDENT */
#include "sysstuff.h"
#include "typedefs.h"

typedef struct {
  int  nr;
  rvec *group_mvcm;
  real *group_mass;
  char **group_name;         /* These two are copies to pointers in */
  unsigned short *group_id;  /* other structures.                   */
} t_vcm;

t_vcm *init_vcm(FILE *fp,t_topology *top,t_mdatoms *md,
		int start,int homenr,int nstcomm);

extern void calc_vcm(FILE *log,int homenr,int start,
		     real mass[],rvec v[],rvec vcm);

extern void do_stopcm(FILE *log,int homenr,int start,
		      rvec v[],rvec mvcm,real tm,real invmass[]);

extern void check_cm(FILE *log,rvec mvcm,real tm);

/* remove global rotation of system by fitting to structure of nstcomm
   steps ago */
extern void do_stoprot(FILE *log, int natoms, rvec box, rvec x[], 
		       real mass[]);

/* Do a per group center of mass things */
extern void calc_vcm_grp(FILE *log,int homenr,int start,real mass[],rvec v[],
			 t_vcm *vcm);

extern void do_stopcm_grp(FILE *log,int homenr,int start,rvec v[],
			  t_vcm *vcm,real invmass[]);

extern void check_cm_grp(FILE *log,t_vcm *vcm);
			 
#endif /* _vcm_h */
