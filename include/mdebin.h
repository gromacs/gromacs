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

#ifndef _mdebin_h
#define _mdebin_h

static char *SRCID_mdebin_h = "$Id$";

#ifdef HAVE_IDENT
#ident	"@(#) mdebin.h 1.12 2/2/97"
#endif /* HAVE_IDENT */
#include "typedefs.h"
#include "sysstuff.h"
#include "ebin.h"
#include "enxio.h"

typedef struct {
  t_ebin *ebin;
  int    ie,ib,isvir,ifvir,ipres,ivir,isurft,itc,iu,imu,ivcos,ivisc;
  int    nE,nEg,nEc,nTC,nU;
  int    *igrp;
} t_mdebin;

extern t_mdebin *init_mdebin(int fp_ene,t_groups *grps,t_atoms *atoms,
			     t_idef *idef,bool bLR,bool BLJLR,bool bBHAM,
			     bool b14,bool bFEP,bool bPcoupl,bool bDispCorr,
			     t_commrec *cr);
/* Initiate MD energy bin and write header to energy file. */

extern void upd_mdebin(t_mdebin *md,FILE *fp_dgdl,
		       real tmass,int step,real time,
		       real ener[],
		       matrix box,
		       tensor svir,
		       tensor fvir,
		       tensor vir,
		       tensor pres,
		       t_groups *grps,
		       rvec mu_tot);
     
extern void print_ebin(int fp_ene,bool bEne,bool bDR,
		       FILE *log,int steps,real time,real lamb,
		       real SAfactor,int mode,bool bCompact,
		       t_mdebin *md,t_groups *grps,t_atoms *atoms);

#endif	/* _mdebin_h */

