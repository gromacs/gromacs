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

#ifndef _do_gct_h
#define _do_gct_h

static char *SRCID_do_gct_h = "$Id$";

#include <stdio.h>
#include "typedefs.h"
#include "filenm.h"
#include "network.h"

enum { eoPres, eoEpot, eoVir, eoPolarizability, eoDipole, 
       eoMemory, eoInter, eoUseVirial, eoNR };
extern char *eoNames[eoNR];

typedef struct {
  int  at_i,at_j;   	/* Atom type # for i and j                   	*/
  int  eObs;		/* Observable to couple to              	*/
  bool bPrint;		/* Does this struct have to be printed		*/
  real c6,c12;		/* Actual value of params			*/
  real xi_6,xi_12;	/* Constants for coupling C6 and C12 		*/
} t_coupl_LJ;

typedef struct {
  int  at_i,at_j;   	/* Atom type # for i and j                   	*/
  int  eObs;		/* Observable to couple to              	*/
  bool bPrint;		/* Does this struct have to be printed		*/
  real a,b,c;		/* Actual value of params			*/
  real xi_a,xi_b,xi_c;	/* Constants for coupling A, B and C 		*/
} t_coupl_BU;

typedef struct {
  int  at_i;		/* Atom type					*/
  int  eObs;		/* Observable to couple to              	*/
  bool bPrint;		/* Does this struct have to be printed		*/
  real Q;		/* Actual value of charge			*/
  real xi_Q;		/* Constant for coupling Q			*/
} t_coupl_Q;

typedef struct {
  int       type;	/* Type number in the iparams struct	*/
  int       eObs;       /* Observable to couple to              */
  t_iparams xi;	        /* Parameters that need to be changed	*/
  t_iparams iprint;
} t_coupl_iparams;

typedef struct {
  real       pres0,pres;
  real       vir0,vir;
  real       epot0,epot;
  int        nLJ,nBU,nQ,nIP;
  t_coupl_LJ *tcLJ;
  t_coupl_BU *tcBU;
  t_coupl_Q  *tcQ;
  t_coupl_iparams *tIP;
  real       polarizability;
  real       dipole;
  int        nmemory;
  bool       bInter;
  bool       bVirial;
} t_coupl_rec;

extern void write_gct(char *fn,t_coupl_rec *tcr,t_idef *idef);

extern void read_gct(char *fn,t_coupl_rec *tcr);

extern void comm_tcr(FILE *log,t_commrec *cr,t_coupl_rec **tcr);

extern void copy_ff(t_coupl_rec *tcr,t_forcerec *fr,t_mdatoms *md,
		    t_idef *idef);

extern t_coupl_rec *init_coupling(FILE *log,int nfile,t_filenm fnm[],
				  t_commrec *cr,t_forcerec *fr,t_mdatoms *md,
				  t_idef *idef);

extern void do_coupling(FILE *log,int nfile,t_filenm fnm[],
			t_coupl_rec *tcr,real t,int step,real ener[],
			t_forcerec *fr,t_inputrec *ir,bool bMaster,
			t_mdatoms *md,t_idef *idef,real mu_aver,int nmols,
			t_commrec *cr,matrix box,tensor virial);
		     
#endif
