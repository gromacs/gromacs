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
 * Gromacs Runs On Most of All Computer Systems
 */

#ifndef	_mdebin_h
#define	_mdebin_h

#ifdef HAVE_IDENT
#ident	"@(#) mdebin.h 1.12 2/2/97"
#endif /* HAVE_IDENT */
#include "typedefs.h"
#include "sysstuff.h"
#include "ebin.h"

typedef struct {
  t_ebin *ebin;
  int    ie,ib,isvir,ifvir,ipres,ivir,itc,iu;
  int    nE,nEg,nEc,nTC,nU;
  int    *igrp;
} t_mdebin;

extern t_mdebin *init_mdebin(FILE *ene,t_groups *grps,t_atoms *atoms,
			     bool bLR,bool bBHAM,bool b14);
/* Initiate MD energy bin and write header to energy file. */

extern void upd_mdebin(t_mdebin *md,real tmass,int step,
		       real ener[],
		       matrix box,
		       tensor svir,
		       tensor fvir,
		       tensor vir,
		       tensor pres,
		       t_groups *grps);
     
extern void print_ebin(FILE *ene,FILE *log,int steps,real time,real lamb,
		       real SAfactor,int mode,bool bCompact,
		       t_mdebin *md,t_groups *grps,t_atoms *atoms);

#endif	/* _mdebin_h */

