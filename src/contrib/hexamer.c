/*
 * 
 *                This source code is part of
 * 
 *                 G   R   O   M   A   C   S
 * 
 *          GROningen MAchine for Chemical Simulations
 * 
 *                        VERSION 3.3.99_development_20071104
 * Written by David van der Spoel, Erik Lindahl, Berk Hess, and others.
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2006, The GROMACS development team,
 * check out http://www.gromacs.org for more information.

 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 * 
 * If you want to redistribute modifications, please consider that
 * scientific software is very special. Version control is crucial -
 * bugs must be traceable. We will be happy to consider code for
 * inclusion in the official distribution, but derived work must not
 * be called official GROMACS. Details are found in the README & COPYING
 * files - if they are missing, get the official version at www.gromacs.org.
 * 
 * To help us fund GROMACS development, we humbly ask that you cite
 * the papers on the package - you can find them in the top README file.
 * 
 * For more info, check our website at http://www.gromacs.org
 * 
 * And Hey:
 * Groningen Machine for Chemical Simulation
 */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <math.h>
#include <string.h>

#include "gromacs/fileio/pdbio.h"
#include "gromacs/fileio/confio.h"
#include "gromacs/utility/smalloc.h"
#include "gromacs/topology/symtab.h"
#include "macros.h"
#include "copyrite.h"
#include "gromacs/commandline/pargs.h"
#include "gromacs/math/vec.h"
#include "typedefs.h"
#include "gbutil.h"
#include "gromacs/math/units.h"

void copy_atom(t_symtab *tab,t_atoms *a1,int i1,t_atoms *a2,int i2,
	       rvec xin[],rvec xout[],rvec vin[],rvec vout[])
{
  a2->atom[i2]     = a1->atom[i1];
  a2->atomname[i2] = put_symtab(tab,*a1->atomname[i1]);
  a2->resname[a2->atom[i2].resnr] =
    put_symtab(tab,*a1->resname[a1->atom[i1].resnr]);
  copy_rvec(xin[i1],xout[i2]);
  copy_rvec(vin[i1],vout[i2]);
}

static void rotate_x(int natom,rvec xin[],real angle,rvec xout[],
		     gmx_bool bZ,gmx_bool bUpsideDown,real dz)
{
  int i;
  matrix mat;
  
  angle *= DEG2RAD;
  clear_mat(mat);
  if (bZ) {
    mat[XX][XX] = cos(angle);
    mat[XX][YY] = sin(angle);
    mat[YY][XX] = -sin(angle);
    mat[YY][YY] = cos(angle);
    mat[ZZ][ZZ] = 1;
  }
  else {
    mat[XX][XX] = 1;
    mat[YY][YY] = cos(angle);
    mat[YY][ZZ] = sin(angle);
    mat[ZZ][YY] = -sin(angle);
    mat[ZZ][ZZ] = cos(angle);
  }
    
  for(i=0; (i<natom); i++) {
    mvmul(mat,xin[i],xout[i]);
    if (bUpsideDown)
      xout[i][ZZ] *= -1;
    xout[i][ZZ] += dz;
  }
}

static void prep_x(int natom,rvec x[],real rDist,real rAngleZ,real rAngleX)
{
  int  i;
  rvec xcm;
  rvec *xx;
  
  /* Center on Z-axis */
  clear_rvec(xcm);
  for(i=0; (i<natom); i++) {
    xcm[XX] += x[i][XX];
    xcm[YY] += x[i][YY];
    xcm[ZZ] += x[i][ZZ];
  }
  xcm[XX] /= natom;
  xcm[YY] /= natom;
  xcm[ZZ] /= natom;
  for(i=0; (i<natom); i++) {
    x[i][XX] -= xcm[XX];
    x[i][YY] -= xcm[YY];
    x[i][ZZ] -= xcm[ZZ];
  }
  if (rAngleZ != 0) {
    snew(xx,natom);
    rotate_x(natom,x,rAngleZ,xx,TRUE,FALSE,0);
    for(i=0; (i<natom); i++) 
      copy_rvec(xx[i],x[i]);
    sfree(xx);
  }
  if (rAngleX != 0) {
    snew(xx,natom);
    rotate_x(natom,x,rAngleX,xx,FALSE,FALSE,0);
    for(i=0; (i<natom); i++) 
      copy_rvec(xx[i],x[i]);
    sfree(xx);
  }
  if (rDist > 0) {
    for(i=0; (i<natom); i++) 
      x[i][XX] += rDist;
  }
}

int main(int argc, char *argv[])
{
  t_symtab tab;
  static char *desc[] = {
    "[TT]hexamer[tt] takes a single input coordinate file and makes five symmetry",
    "related copies."
  };
#define NPA asize(pa)
  t_filenm fnm[] = {
    { efSTX, "-f", NULL, ffREAD },
    { efPDB, "-o", NULL, ffWRITE }
  };
#define NFILE asize(fnm)
  gmx_bool bCenter    = FALSE;
  gmx_bool bTrimer    = FALSE;
  gmx_bool bAlternate = FALSE;
  real rDist = 0,rAngleZ = 0,rAngleX = 0, alterz = 0;
  t_pargs pa[] = {
    { "-center",   FALSE, etBOOL,  {&bCenter}, 
      "Center molecule on Z-axis first" },
    { "-trimer",   FALSE, etBOOL,  {&bTrimer},
      "Make trimer rather than hexamer" },
    { "-alternate",FALSE, etBOOL,  {&bAlternate},
      "Turn every other molecule upside down" },
    { "-alterz",   FALSE, etREAL,  {&alterz},
      "Add this amount to Z-coordinate in every other molecule" },
    { "-radius",   FALSE, etREAL,  {&rDist},
      "Distance of protein axis from Z-axis (implies [TT]-center[tt])" },
    { "-anglez",   FALSE, etREAL,  {&rAngleZ},
      "Initial angle of rotation around Z-axis of protein" },
    { "-anglex",   FALSE, etREAL,  {&rAngleX},
      "Initial angle of rotation around X-axis of protein" }
  };
#define NPA asize(pa)
  FILE    *fp;
  int     i,iout,now,natom;
  rvec    *xin,*vin,*xout;
  matrix  box;
  t_atoms atoms,aout;
  char    *infile,*outfile,title[256],buf[32];
  
  CopyRight(stderr,argv[0]);
  parse_common_args(&argc,argv,0,NFILE,fnm,NPA,pa,
		    asize(desc),desc,0,NULL);
  bCenter = bCenter || (rDist > 0) || bAlternate;
  
  infile  = ftp2fn(efSTX,NFILE,fnm);
  outfile = ftp2fn(efPDB,NFILE,fnm);
  
  get_stx_coordnum(infile,&natom);
  init_t_atoms(&atoms,natom,TRUE);
  snew(xin,natom);
  snew(xout,natom);
  snew(vin,natom);
  read_stx_conf(infile,title,&atoms,xin,vin,box);
  printf("Read %d atoms\n",atoms.nr); 
  
  if (bCenter) 
    prep_x(atoms.nr,xin,rDist,rAngleZ,rAngleX);
  
  fp = gmx_ffopen(outfile,"w");
  for(i=0; (i<(bTrimer ? 3 : 6)); i++) {
    rotate_x(atoms.nr,xin,i*(bTrimer ? 120.0 : 60.0),xout,TRUE,
	     bAlternate && ((i % 2) != 0),alterz*(((i % 2) == 0) ? 0 : 1));
    sprintf(buf,"Rotated %d degrees",i*(bTrimer ? 120 : 60));
    write_pdbfile(fp,buf,&atoms,xout,box,'A'+i,1+i);
  }
  gmx_ffclose(fp);
  
  gmx_thanx(stderr);
  
  return 0;
}
