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
static char *SRCID_anaf_c = "$Id$";

#include <stdio.h>
#include <string.h>
#include <math.h>
#include "main.h"
#include "macros.h"
#include "futil.h"
#include "statutil.h"
#include "copyrite.h"
#include "sysstuff.h"
#include "txtdump.h"
#include "fatal.h"
#include "xtcio.h"
#include "enxio.h"
#include "assert.h"
#include "smalloc.h"
#include "gmxfio.h"
#include "tpxio.h"
#include "trnio.h"
#include "txtdump.h"
#include "vec.h"

static char *nm[5]  = { "OW", "HW1", "HW2", "DW", "SW" };
  
static void list_trn(char *fn)
{
  static real mass[5] = { 15.9994, 1.008, 1.008, 0.0, 0.0 };
  int         i,j,m,fpread,fpwrite,nframe;
  rvec        *x,*v,*f,fmol[2],xcm[2],torque[j],dx;
  real        mmm,len;
  matrix      box;
  t_trnheader trn;
  bool        bOK;

  printf("Going to open %s\n",fn);
  fpread  = open_trn(fn,"r"); 
  fpwrite = open_tpx(NULL,"w");
  fio_setdebug(fpwrite,TRUE);
  
  mmm=mass[0]+2*mass[1];
  for(i=0; (i<5); i++) 
    mass[i] /= mmm;
  
  nframe = 0;
  while (fread_trnheader(fpread,&trn,&bOK)) {
    snew(x,trn.natoms);
    snew(v,trn.natoms);
    snew(f,trn.natoms);
    if (fread_htrn(fpread,&trn,
		   trn.box_size ? box : NULL,
		   trn.x_size   ? x : NULL,
		   trn.v_size   ? v : NULL,
		   trn.f_size   ? f : NULL)) {
		   
      if (trn.x_size && trn.f_size) {
	printf("There are %d atoms\n",trn.natoms);
	for(j=0; (j<2); j++) {
	  clear_rvec(xcm[j]);
	  clear_rvec(fmol[j]);
	  clear_rvec(torque[j]);
	  for(i=5*j; (i<5*j+5); i++) {
	    rvec_inc(fmol[j],f[i]);
	    for(m=0; (m<DIM); m++)
	      xcm[j][m] += mass[i%5]*x[i][m];
	  }
	  for(i=5*j; (i<5*j+5); i++) {
	    rvec_dec(x[i],xcm[j]);
	    oprod(x[i],f[i],dx);
	    rvec_inc(torque[j],dx);
	    rvec_inc(x[i],xcm[j]);
	  }
	}
	pr_rvecs(stdout,0,"FMOL  ",fmol,2);
	pr_rvecs(stdout,0,"TORQUE",torque,2);
	printf("Distance matrix Water1-Water2\n%5s","");
	for(j=0; (j<5); j++) 
	  printf("  %10s",nm[j]);
	printf("\n");
	for(j=0; (j<5); j++) {
	  printf("%5s",nm[j]);
	  for(i=5; (i<10); i++) {
	    rvec_sub(x[i],x[j],dx);
	    len = sqrt(iprod(dx,dx));
	    printf("  %10.7f",len);
	  }
	  printf("\n");
	}
      }
    }
    sfree(x);
    sfree(v);
    sfree(f);
    nframe++;
  }
  if (!bOK)
    fprintf(stderr,"\nWARNING: Incomplete frame header: nr %d, t=%g\n",
	    nframe,trn.t);
  close_tpx(fpwrite);
  close_trn(fpread);
}

int main(int argc,char *argv[])
{
  static char *desc[] = {
    "gmxdump reads a run input file ([TT].tpa[tt]/[TT].tpr[tt]/[TT].tpb[tt]),",
    "a trajectory ([TT].trj[tt]/[TT].trr[tt]/[TT].xtc[tt]) or an energy",
    "file ([TT].ene[tt]/[TT].edr[tt]) and prints that to standard",
    "output in a readable format. This program is essential for",
    "checking your run input file in case of problems.[PAR]"
  };
  t_filenm fnm[] = {
    { efTRN, "-f", NULL, ffOPTRD }
  };
#define NFILE asize(fnm)
  char *fn;
  
  /* Command line options */
  
  CopyRight(stdout,argv[0]);
  parse_common_args(&argc,argv,0,FALSE,NFILE,fnm,0,NULL,
		    asize(desc),desc,0,NULL);
  
  if (ftp2bSet(efTRN,NFILE,fnm)) {
    fn = ftp2fn(efTRN,NFILE,fnm);
    printf("Going to open %s\n",fn);
    list_trn(fn);
  }
  
  thanx(stdout);

  return 0;
}
