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
 * Great Red Oystrich Makes All Chemists Sane
 */

#ifndef _hxprops_h
#define _hxprops_h

static char *SRCID_hxprops_h = "$Id$";

#include <stdio.h>
#include "typedefs.h"

#define PHI_AHX (-55.0)
#define PSI_AHX (-45.0)
/* Canonical values of the helix phi/psi angles */


typedef struct {
  real phi,psi,pprms2;
  real jcaha;
  real d3,d4,d5,rmsa;
  bool bHelix;
  int  nhx;
  int  nrms,resno;
  int  Cprev,N,H,CA,C,O,Nnext;
  char label[32];
} t_bb;

enum { 
  efhRAD,  efhTWIST, efhRISE, efhLEN,  
  efhDIP,  efhRMS,   efhRMSA, efhCD222,
  efhPPRMS,efhCPHI,  efhPHI,  efhPSI,  
  efhHB3,  efhHB4,   efhHB5,  efhJCA,
  efhAHX,  efhNR 
};

extern real ahx_len(int gnx,atom_id index[],rvec x[],matrix box);
/* Assume we have a list of Calpha atoms only! */

extern real ellipticity(int nres,t_bb bb[]);

extern real radius(FILE *fp,int nca,atom_id ca_index[],rvec x[]);
/* Assume we have calphas */

extern real twist(FILE *fp,int nca,atom_id caindex[],rvec x[]);
/* Calculate the twist of the helix */

extern real pprms(FILE *fp,int nbb,t_bb bb[]);
/* Calculate the average RMS from canonical phi/psi values
 * and the distance per residue
 */
 
extern real ca_phi(int gnx,atom_id index[],rvec x[],matrix box);
/* Assume we have a list of Calpha atoms only! */

extern real dip(int nbb,atom_id bbind[],rvec x[],t_atom atom[]);

extern real rise(int gnx,atom_id index[],rvec x[]);
/* Assume we have a list of Calpha atoms only! */

extern void av_hblen(FILE *fp3,FILE *fp3a,
		     FILE *fp4,FILE *fp4a,
		     FILE *fp5,FILE *fp5a,
		     real t,int nres,t_bb bb[]);
		     
extern void av_phipsi(FILE *fphi,FILE *fpsi,FILE *fphi2,FILE *fpsi2,
		      real t,int nres,t_bb bb[]);

extern t_bb *mkbbind(char *fn,int *nres,int *nbb,int res0,
		     int *nall,atom_id **index,
		     char ***atomname,t_atom atom[],
		     char ***resname);
		     
extern void do_start_end(int nres,t_bb bb[],rvec x[],int *nbb,
			 atom_id bbindex[],int *nca,atom_id caindex[],
			 bool bRange,int rStart,int rEnd);
		     
extern void calc_hxprops(int nres,t_bb bb[],rvec x[],matrix box);

extern void pr_bb(FILE *fp,int nres,t_bb bb[]);

#endif
