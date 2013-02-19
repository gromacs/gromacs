/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2006, The GROMACS development team,
 * check out http://www.gromacs.org for more information.
 * Copyright (c) 2012,2013, by the GROMACS development team, led by
 * David van der Spoel, Berk Hess, Erik Lindahl, and including many
 * others, as listed in the AUTHORS file in the top-level source
 * directory and at http://www.gromacs.org.
 *
 * GROMACS is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public License
 * as published by the Free Software Foundation; either version 2.1
 * of the License, or (at your option) any later version.
 *
 * GROMACS is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with GROMACS; if not, see
 * http://www.gnu.org/licenses, or write to the Free Software Foundation,
 * Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA.
 *
 * If you want to redistribute modifications to GROMACS, please
 * consider that scientific software is very special. Version
 * control is crucial - bugs must be traceable. We will be happy to
 * consider code for inclusion in the official distribution, but
 * derived work must not be called official GROMACS. Details are found
 * in the README & COPYING files - if they are missing, get the
 * official version at http://www.gromacs.org.
 *
 * To help us fund GROMACS development, we humbly ask that you cite
 * the research papers on the package. Check out http://www.gromacs.org.
 */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

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
#include "gmx_fatal.h"
#include "xtcio.h"
#include "enxio.h"
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
  int         i,j=0,m,fpread,fpwrite,nframe;
  rvec        *x,*v,*f,fmol[2],xcm[2],torque[j],dx;
  real        mmm,len;
  matrix      box;
  t_trnheader trn;
  gmx_bool        bOK;

  printf("Going to open %s\n",fn);
  fpread  = open_trn(fn,"r"); 
  fpwrite = open_tpx(NULL,"w");
  gmx_fio_setdebug(fpwrite,TRUE);
  
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
	    cprod(x[i],f[i],dx);
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
    "[TT]gmxdump[tt] reads a run input file ([TT].tpa[tt]/[TT].tpr[tt]/[TT].tpb[tt]),",
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
  parse_common_args(&argc,argv,0,NFILE,fnm,0,NULL,
		    asize(desc),desc,0,NULL);
  
  if (ftp2bSet(efTRN,NFILE,fnm)) {
    fn = ftp2fn(efTRN,NFILE,fnm);
    printf("Going to open %s\n",fn);
    list_trn(fn);
  }
  
  thanx(stderr);

  return 0;
}
