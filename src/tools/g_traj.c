/*
 * $Id$
 * 
 *                This source code is part of
 * 
 *                 G   R   O   M   A   C   S
 * 
 *          GROningen MAchine for Chemical Simulations
 * 
 *                        VERSION 3.0
 * 
 * Copyright (c) 1991-2001
 * BIOSON Research Institute, Dept. of Biophysical Chemistry
 * University of Groningen, The Netherlands
 * 
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
 * Do check out http://www.gromacs.org , or mail us at gromacs@gromacs.org .
 * 
 * And Hey:
 * GROup of MAchos and Cynical Suckers
 */
static char *SRCID_g_traj_c = "$Id$";
#include <math.h>
#include <string.h>
#include "statutil.h"
#include "sysstuff.h"
#include "typedefs.h"
#include "smalloc.h"
#include "macros.h"
#include "vec.h"
#include "pbc.h"
#include "copyrite.h"
#include "futil.h"
#include "statutil.h"
#include "rdgroup.h"
#include "mshift.h"
#include "xvgr.h"
#include "tpxio.h"
#include "rmpbc.h"
#include "physics.h"

static void low_print_data(FILE *fp,real time,rvec x[],int n,atom_id *index,
		    bool bDim[])
{
  int i,d;
  
  fprintf(fp," %g",time);
  for(i=0; i<n; i++) {
    for(d=0; d<DIM; d++)
      if (bDim[d]) {
	if (index)
	  fprintf(fp,"\t%g",x[index[i]][d]);
	else
	  fprintf(fp,"\t%g",x[i][d]);
      }
    if (bDim[DIM]) {
      if (index)
	  fprintf(fp,"\t%g",norm(x[index[i]]));
	else
	  fprintf(fp,"\t%g",norm(x[i]));
    }
  }
  fprintf(fp,"\n");
}

static void average_data(rvec x[],rvec xav[],real *mass,
			 int ngrps,int isize[],atom_id **index)
{
  int  g,i;
  real m,mtot;
  rvec tmp;

  for(g=0; g<ngrps; g++) {
    clear_rvec(xav[g]);
    mtot = 0;
    for(i=0; i<isize[g]; i++)
      if (mass) {
	m = mass[index[g][i]];
	svmul(m,x[index[g][i]],tmp);
	rvec_inc(xav[g],tmp);
	mtot += m;
      } else
	rvec_inc(xav[g],x[index[g][i]]);
    if (mass)
      svmul(1/mtot,xav[g],xav[g]);
  }
}
  
static void print_data(FILE *fp,real time,rvec x[],real *mass,bool bCom,
		       int ngrps,int isize[],atom_id **index,bool bDim[])
{
  static rvec *xav=NULL;
  
  if (bCom) {
    if (xav==NULL)
      snew(xav,ngrps);
    average_data(x,xav,mass,ngrps,isize,index);
    low_print_data(fp,time,xav,ngrps,NULL,bDim);
  } else
    low_print_data(fp,time,x,isize[0],index[0],bDim);
}

static void make_legend(FILE *fp,int ngrps,int isize,atom_id index[],
			char **name,bool bCom,bool bMol,bool bDim[])
{
  char **leg;
  char *dimtxt[] = { " X", " Y", " Z", "" };
  int n,i,j,d;
  
  if (bCom)
    n = ngrps;
  else
    n = isize;

  snew(leg,4*n);
  j=0;
  for(i=0; i<n; i++)
    for(d=0; d<=DIM; d++)
      if (bDim[d]) {
	snew(leg[j],STRLEN);
	if (bMol) 
	  sprintf(leg[j],"mol %d%s",index[i]+1,dimtxt[d]);
	else if (bCom)
	  sprintf(leg[j],"%s%s",name[i],dimtxt[d]);
	else
	  sprintf(leg[j],"atom %d%s",index[i]+1,dimtxt[d]);
	j++;
      }
  xvgr_legend(fp,j,leg);
  for(i=0; i<j; i++)
    sfree(leg[i]);
  sfree(leg);
}

static real ekrot(rvec x[],rvec v[],real mass[],int isize,atom_id index[])
{
  matrix TCM,L;
  real   tm,m0,lxx,lxy,lxz,lyy,lyz,lzz,ekrot;
  rvec   dx,a0,ocm;
  rvec   xcm,vcm,acm;
  int    i,j,m,n;

  clear_rvec(xcm);
  clear_rvec(vcm);
  clear_rvec(acm);
  tm=0.0;
  for(i=0; i<isize; i++) {
    j = index[i];
    m0=mass[j];
    tm+=m0;
    oprod(x[j],v[j],a0);
    for(m=0; (m<DIM); m++) {
      xcm[m]+=m0*x[j][m]; /* c.o.m. position */
      vcm[m]+=m0*v[j][m]; /* c.o.m. velocity */
      acm[m]+=m0*a0[m];   /* rotational velocity around c.o.m. */
    }
  }
  oprod(xcm,vcm,a0);
  for(m=0; (m<DIM); m++) {
    xcm[m]/=tm;
    vcm[m]/=tm;
    acm[m]-=a0[m]/tm;
  }

  lxx=lxy=lxz=lyy=lyz=lzz=0.0;
  for(i=0; i<isize; i++) {
    j = index[i];
    m0=mass[j];
    for(m=0; (m<DIM); m++)
      dx[m]=x[j][m]-xcm[m];
    lxx+=dx[XX]*dx[XX]*m0;
    lxy+=dx[XX]*dx[YY]*m0;
    lxz+=dx[XX]*dx[ZZ]*m0;
    lyy+=dx[YY]*dx[YY]*m0;
    lyz+=dx[YY]*dx[ZZ]*m0;
    lzz+=dx[ZZ]*dx[ZZ]*m0;
  }
  clear_mat(L);
  
  L[XX][XX]=lyy+lzz;
  L[YY][XX]=-lxy;
  L[ZZ][XX]=-lxz;
  L[XX][YY]=-lxy;
  L[YY][YY]=lxx+lzz;
  L[ZZ][YY]=-lyz;
  L[XX][ZZ]=-lxz;
  L[YY][ZZ]=-lyz;
  L[ZZ][ZZ]=lxx+lyy;
  
  m_inv(L,TCM);

  /* Compute omega (hoeksnelheid) */
  clear_rvec(ocm);
  ekrot=0;
  for(m=0; (m<DIM); m++) {
    for(n=0; (n<DIM); n++)
      ocm[m]+=TCM[m][n]*acm[n];
    ekrot+=ocm[m]*acm[m];
  }
  return ekrot;
}

static real temp(rvec v[],real mass[],int isize,atom_id index[])
{
  matrix TCM,L;
  real   m,ekin=0;
  int    i,j;

  for(i=0; i<isize; i++) {
    j = index[i];
    ekin += mass[j]*norm2(v[j]);
  }

  return ekin/(2*isize*BOLTZ);
}

static void remove_jump(matrix box,int natoms,rvec xp[],rvec x[])
{
  rvec hbox;
  int d,i,m;

  for(d=0; d<DIM; d++)
    hbox[d] = 0.5*box[d][d];
  for(i=0; i<natoms; i++)
    for(m=DIM-1; m>=0; m--) {
      while (x[i][m]-xp[i][m] <= -hbox[m])
	for(d=0; d<=m; d++)
	  x[i][d] += box[m][d];
      while (x[i][m]-xp[i][m] > hbox[m])
	for(d=0; d<=m; d++)
	  x[i][d] -= box[m][d];
    }
}

int main(int argc,char *argv[])
{
  static char *desc[] = {
    "g_traj plots coordinates, velocities, forces and/or the box.",
    "With [TT]-com[tt] the coordinates, velocities and forces are",
    "calculated for the center of mass of each group.",
    "When [TT]-mol[tt] is set, the numbers in the index file are",
    "interpreted as molecule numbers and the same procedure as with",
    "[TT]-com[tt] is used for each molecule.[PAR]",
    "Option [TT]-ot[tt] plots the temperature of each group.",
    "This implies [TT]-com[tt].[PAR]",
    "Option [TT]-ekr[tt] plots the rotational kinetic energy of each group.",
    "This implies [TT]-com[tt]."
  };
  static bool bMol=FALSE,bCom=FALSE,bNoJump=FALSE;
  static bool bX=TRUE,bY=TRUE,bZ=TRUE,bNorm=FALSE;
  t_pargs pa[] = {
    { "-com", FALSE, etBOOL, {&bCom},
      "Plot data for the com of each group" },
    { "-mol", FALSE, etBOOL, {&bMol},
      "Index contains molecule numbers iso atom numbers" },
    { "-nojump", FALSE, etBOOL, {&bNoJump},
      "Remove jumps of atoms across the box" },
    { "-x", FALSE, etBOOL, {&bX},
      "Plot X-component" },
    { "-y", FALSE, etBOOL, {&bY},
      "Plot Y-component" },
    { "-z", FALSE, etBOOL, {&bZ},
      "Plot Z-component" },
    { "-len", FALSE, etBOOL, {&bNorm},
      "Plot vector length" }
    
  };
  FILE       *outx=NULL,*outv=NULL,*outf=NULL,*outb=NULL,*outt=NULL,*outekr=NULL;
  t_topology top;
  real       *mass;
  char       title[STRLEN],*indexfn;
  t_trxframe fr;
  int        flags;
  rvec       *xtop,*xp=NULL;
  matrix     topbox;
  int        status;
  int        i,j,n;
  int        ngrps;
  char       **grpname;
  int        *isize0,*isize;
  atom_id    **index0,**index;
  atom_id    *a,*atndx;
  t_block    *mols;
  bool       bTop,bOX,bOV,bOF,bOB,bOT,bEKR,bDim[4],bDum[4];
  char       *box_leg[6] = { "XX", "YY", "ZZ", "YX", "ZX", "ZY" };

  t_filenm fnm[] = {
    { efTRX, "-f", NULL, ffREAD },
    { efTPS, NULL, NULL, ffREAD },
    { efNDX, NULL, NULL, ffOPTRD },
    { efXVG, "-ox", "coord.xvg", ffOPTWR },
    { efXVG, "-ov", "veloc.xvg", ffOPTWR },
    { efXVG, "-of", "force.xvg", ffOPTWR },
    { efXVG, "-ob", "box.xvg",   ffOPTWR },
    { efXVG, "-ot", "temp.xvg",  ffOPTWR },
    { efXVG, "-ekr","ekrot.xvg", ffOPTWR }
  };
#define NFILE asize(fnm)

  CopyRight(stderr,argv[0]);
  parse_common_args(&argc,argv,PCA_CAN_TIME | PCA_CAN_VIEW,TRUE,
		    NFILE,fnm,asize(pa),pa,asize(desc),desc,0,NULL);

  if (bMol)
    fprintf(stderr,"Interpreting indexfile entries as molecules.\n"
	    "Using center of mass.\n");
  
  bOX  = opt2bSet("-ox",NFILE,fnm);
  bOV  = opt2bSet("-ov",NFILE,fnm);
  bOF  = opt2bSet("-of",NFILE,fnm);
  bOB  = opt2bSet("-ob",NFILE,fnm);
  bOT  = opt2bSet("-ot",NFILE,fnm);
  bEKR = opt2bSet("-ekr",NFILE,fnm);
  if (bMol || bOT || bEKR)
    bCom = TRUE;

  bDim[XX] = bX;
  bDim[YY] = bY;
  bDim[ZZ] = bZ;
  bDim[DIM] = bNorm;

  bTop = read_tps_conf(ftp2fn(efTPS,NFILE,fnm),title,&top,&xtop,NULL,topbox,
		       bCom && (bOX || bOV || bOT || bEKR));
  sfree(xtop);
  if (bMol && !bTop)
    fatal_error(0,"Need a run input file for option -mol");

  if (bMol)
    indexfn = ftp2fn(efNDX,NFILE,fnm);
  else
    indexfn = ftp2fn_null(efNDX,NFILE,fnm);

  if (bCom && !bMol) {
    fprintf(stderr,"How many groups do you want to analyze? ");
    scanf("%d",&ngrps);
  } else
    ngrps = 1;
  snew(grpname,ngrps);
  snew(isize0,ngrps);
  snew(index0,ngrps);
  get_index(&(top.atoms),indexfn,ngrps,isize0,index0,grpname);
  
  if (bMol) {
    mols=&(top.blocks[ebMOLS]);
    a = mols->a;
    atndx = mols->index;
    ngrps = isize0[0];
    snew(isize,ngrps);
    snew(index,ngrps);
    for (i=0; i<ngrps; i++) {
      isize[i] = atndx[index0[0][i]+1] - atndx[index0[0][i]];
      snew(index[i],isize[i]);
      for(j=0; j<isize[i]; j++)
	index[i][j] = a[atndx[index0[0][i]]+j];
    }
  } else {
    isize = isize0;
    index = index0;
  }
  if (bCom) {
    snew(mass,top.atoms.nr);
    for(i=0; i<top.atoms.nr; i++)
      mass[i] = top.atoms.atom[i].m;
  } else
    mass = NULL;

  flags = 0;
  if (bOX) {
    flags = flags | TRX_READ_X;
    outx = xvgropen(opt2fn("-ox",NFILE,fnm),
		    bCom ? "Center of mass" : "Coordinate",
		    "Time (ps)","Coordinate (nm)");
    make_legend(outx,ngrps,isize0[0],index0[0],grpname,bCom,bMol,bDim);
  }
  if (bOV) {
    flags = flags | TRX_READ_V;
    outv = xvgropen(opt2fn("-ov",NFILE,fnm),
		    bCom ? "Center of mass velocity" : "Velocity",
		    "Time (ps)","Velocity (nm/ps)");
   make_legend(outv,ngrps,isize0[0],index0[0],grpname,bCom,bMol,bDim); 
  }
  if (bOF) {
    flags = flags | TRX_READ_F;
    outf = xvgropen(opt2fn("-of",NFILE,fnm),"Force",
		    "Time (ps)","Force (kJ mol\\S-1\\N nm\\S-1\\N)");
    make_legend(outf,ngrps,isize0[0],index0[0],grpname,bCom,bMol,bDim);
  }
  if (bOB) {
    outb = xvgropen(opt2fn("-ob",NFILE,fnm),"Box vector elements",
		    "Time (ps)","(nm)");
   
    xvgr_legend(outb,6,box_leg);
  }
  if (bOT) {
    bDum[XX] = FALSE;
    bDum[YY] = FALSE;
    bDum[ZZ] = FALSE;
    bDum[DIM] = TRUE;
    flags = flags | TRX_READ_V;
    outt = xvgropen(opt2fn("-ot",NFILE,fnm),"Temperature","Time (ps)","(K)");
    make_legend(outt,ngrps,isize[0],index[0],grpname,bCom,bMol,bDum);
  }
  if (bEKR) {
    bDum[XX] = FALSE;
    bDum[YY] = FALSE;
    bDum[ZZ] = FALSE;
    bDum[DIM] = TRUE;
    flags = flags | TRX_READ_X | TRX_READ_V;
    outekr = xvgropen(opt2fn("-ekr",NFILE,fnm),"Center of mass rotation",
		      "Time (ps)","Energy (kJ mol\\S-1\\N)");
    make_legend(outekr,ngrps,isize[0],index[0],grpname,bCom,bMol,bDum);
  }
  if (flags == 0 && !bOB) {
    fprintf(stderr,"Please select one or more output file options\n");
    exit(0);
  }

  read_first_frame(&status,ftp2fn(efTRX,NFILE,fnm),&fr,flags);

  do {
    if (fr.bX && bNoJump && fr.bBox) {
      if (xp)
	remove_jump(fr.box,fr.natoms,xp,fr.x);
      else 
	snew(xp,fr.natoms);
      for(i=0; i<fr.natoms; i++)
	copy_rvec(fr.x[i],xp[i]);
    }
    
    if (fr.bX && bCom)
      rm_pbc(&(top.idef),fr.natoms,fr.box,fr.x,fr.x);

    if (bOX && fr.bX)
      print_data(outx,fr.time,fr.x,mass,bCom,ngrps,isize,index,bDim);
    if (bOV && fr.bV)
      print_data(outv,fr.time,fr.v,mass,bCom,ngrps,isize,index,bDim);
    if (bOF && fr.bF)
      print_data(outf,fr.time,fr.f,NULL,bCom,ngrps,isize,index,bDim);
    if (bOB && fr.bBox)
      fprintf(outb,"\t%g\t%g\t%g\t%g\t%g\t%g\t%g\n",fr.time,
	      fr.box[XX][XX],fr.box[YY][YY],fr.box[ZZ][ZZ],
	      fr.box[YY][XX],fr.box[ZZ][XX],fr.box[ZZ][YY]);
    if (bOT && fr.bV) {
      fprintf(outt," %g",fr.time);
      for(i=0; i<ngrps; i++)
	fprintf(outt,"\t%g",temp(fr.v,mass,isize[i],index[i]));
      fprintf(outt,"\n");
    }
    if (bEKR && fr.bX && fr.bV) {
      fprintf(outekr," %g",fr.time);
      for(i=0; i<ngrps; i++)
	fprintf(outekr,"\t%g",ekrot(fr.x,fr.v,mass,isize[i],index[i]));
      fprintf(outekr,"\n");
    }
  } while(read_next_frame(status,&fr));
  

  /* clean up a bit */
  close_trj(status);
  if (bOX) {
    fclose(outx);
    do_view(opt2fn("-ox",NFILE,fnm), NULL);
  }
  if (bOV) {
    fclose(outv);
    do_view(opt2fn("-ov",NFILE,fnm), NULL);
  }
  if (bOF) {
    fclose(outf);
    do_view(opt2fn("-of",NFILE,fnm), NULL);
  }
  if (bOB) {
    fclose(outb);
    do_view(opt2fn("-ob",NFILE,fnm), NULL);
  }
  if (bOT) {
    fclose(outt);
    do_view(opt2fn("-ot",NFILE,fnm), NULL);
  }
  if (bEKR) {
    fclose(outekr);
    do_view(opt2fn("-ekr",NFILE,fnm), NULL);
  } 

  thanx(stderr);
  
  return 0;
}

