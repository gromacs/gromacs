/*
 * $Id$
 * 
 *                This source code is part of
 * 
 *                 G   R   O   M   A   C   S
 * 
 *          GROningen MAchine for Chemical Simulations
 * 
 *                        VERSION 3.2.0
 * Written by David van der Spoel, Erik Lindahl, Berk Hess, and others.
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2004, The GROMACS development team,
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
 * Green Red Orange Magenta Azure Cyan Skyblue
 */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

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
#include "index.h"
#include "mshift.h"
#include "xvgr.h"
#include "tpxio.h"
#include "rmpbc.h"
#include "physics.h"
#include "nrjac.h"
#include "confio.h"

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
  int  g,i,ind,d;
  real m,mtot;
  rvec tmp;
  double sum[DIM];

  for(g=0; g<ngrps; g++) {
    for(d=0; d<DIM; d++)
      sum[d] = 0;
    clear_rvec(xav[g]);
    mtot = 0;
    for(i=0; i<isize[g]; i++) {
      ind = index[g][i];
      if (mass) {
	m = mass[ind];
	svmul(m,x[ind],tmp);
	for(d=0; d<DIM; d++)
	  sum[d] += tmp[d];
	mtot += m;
      } else
	for(d=0; d<DIM; d++)
	  sum[d] += x[ind][d];
    }
    if (mass) {
      for(d=0; d<DIM; d++)
	xav[g][d] = sum[d]/mtot;
    } else {
      for(d=0; d<DIM; d++)
	xav[g][d] = sum[d]/isize[g];
    }
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
  static real **TCM=NULL,**L;
  real   tm,m0,lxx,lxy,lxz,lyy,lyz,lzz,ekrot;
  rvec   dx,a0,ocm;
  rvec   xcm,vcm,acm;
  int    i,j,m,n;

  if (TCM == NULL) {
    snew(TCM,DIM);
    for(i=0; i<DIM; i++)
      snew(TCM[i],DIM);
    snew(L,DIM);
    for(i=0; i<DIM; i++)
      snew(L[i],DIM);
  }

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
  
  L[XX][XX]=lyy+lzz;
  L[YY][XX]=-lxy;
  L[ZZ][XX]=-lxz;
  L[XX][YY]=-lxy;
  L[YY][YY]=lxx+lzz;
  L[ZZ][YY]=-lyz;
  L[XX][ZZ]=-lxz;
  L[YY][ZZ]=-lyz;
  L[ZZ][ZZ]=lxx+lyy;

  m_inv_gen(L,DIM,TCM);

  /* Compute omega (hoeksnelheid) */
  clear_rvec(ocm);
  ekrot=0;
  for(m=0; (m<DIM); m++) {
    for(n=0; (n<DIM); n++)
      ocm[m]+=TCM[m][n]*acm[n];
    ekrot+=0.5*ocm[m]*acm[m];
  }
  return ekrot;
}

static real ektrans(rvec v[],real mass[],int isize,atom_id index[])
{
  rvec   mvcom;
  real   mtot=0;
  int    i,j,d;
  
  clear_rvec(mvcom);
  for(i=0; i<isize; i++) {
    j = index[i];
    for(d=0; d<DIM; d++)
      mvcom[d] += mass[j]*v[j][d];
    mtot += mass[j];
  }

  return norm2(mvcom)/(mtot*2);
}

static real temp(rvec v[],real mass[],int isize,atom_id index[])
{
  real ekin2=0;
  int  i,j;

  for(i=0; i<isize; i++) {
    j = index[i];
    ekin2 += mass[j]*norm2(v[j]);
  }

  return ekin2/(3*isize*BOLTZ);
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

static void write_pdb_bfac(char *fname,char *title,t_atoms *atoms,matrix box,
			   int isize,atom_id *index,int nfr,rvec *x,rvec *sum,
			   bool bDim[],real scale_factor)
{
  FILE    *fp;
  real    max,len2,scale;
  atom_id maxi; 
  int     i,m,onedim;
  bool    bOne;

  if (nfr == 0) {
    fprintf(stderr,"No frames found for %s, will not write %s\n",title,fname);
  } else {
    fprintf(stderr,"Used %d frames for %s\n",nfr,title);
    onedim = -1;
    if (!bDim[DIM]) {
      m = 0;
      for(i=0; i<DIM; i++)
	if (bDim[i]) {
	  onedim = i;
	  m++;
	}
      if (m != 1)
	onedim = -1;
    }
    scale = 1.0/nfr;
    for(i=0; i<isize; i++)
      svmul(scale,x[index[i]],x[index[i]]);
    max  = 0;
    maxi = 0;
    for(i=0; i<isize; i++) {
      len2 = 0;
      for(m=0; m<DIM; m++)
	if (bDim[m] || bDim[DIM])
	  len2 += sqr(sum[index[i]][m]);
      if (len2 > max) {
	max  = len2;
	maxi = index[i];
      }
    }
    if (scale_factor != 0) {
      scale = scale_factor;
    } else {
      if (max == 0)
	scale = 1;
      else
	scale = 10.0/sqrt(max);
    }
    
    fprintf(stdout,"Maximum %s is %g on atom %d %s, res. %s %d\n",
	    title,sqrt(max)/nfr,maxi+1,*(atoms->atomname[maxi]),
	    *(atoms->resname[atoms->atom[maxi].resnr]),
	    atoms->atom[maxi].resnr+1);
    
    if (atoms->pdbinfo == NULL)
      snew(atoms->pdbinfo,atoms->nr);
    if (onedim == -1) {
      for(i=0; i<isize; i++) {
	len2 = 0;
	for(m=0; m<DIM; m++) 
	  if (bDim[m] || bDim[DIM])
	    len2 += sqr(sum[index[i]][m]);
	atoms->pdbinfo[index[i]].bfac = sqrt(len2)*scale;
      }
    } else {
      for(i=0; i<isize; i++)
	atoms->pdbinfo[index[i]].bfac = sum[index[i]][onedim]*scale;
    }
    write_sto_conf_indexed(fname,title,atoms,x,NULL,box,isize,index);
  }
}

int gmx_traj(int argc,char *argv[])
{
  static char *desc[] = {
    "g_traj plots coordinates, velocities, forces and/or the box.",
    "With [TT]-com[tt] the coordinates, velocities and forces are",
    "calculated for the center of mass of each group.",
    "When [TT]-mol[tt] is set, the numbers in the index file are",
    "interpreted as molecule numbers and the same procedure as with",
    "[TT]-com[tt] is used for each molecule.[PAR]",
    "Option [TT]-ot[tt] plots the temperature of each group,",
    "provided velocities are present in the trajectory file.",
    "No corrections are made for constrained degrees of freedom!",
    "This implies [TT]-com[tt].[PAR]",
    "Options [TT]-ekt[tt] and [TT]-ekr[tt] plot the translational and",
    "rotational kinetic energy of each group,", 
    "provided velocities are present in the trajectory file.",
    "This implies [TT]-com[tt].[PAR]",
    "Options [TT]-cv[tt] and [TT]-cf[tt] write the average velocities",
    "and average forces as temperature factors to a pdb file with",
    "the average coordinates. The temperature factors are scaled such",
    "that the maximum is 10. The scaling can be changed with the option",
    "[TT]-scale[tt]. To get the velocities or forces of one",
    "frame set both [TT]-b[tt] and [TT]-e[tt] to the time of",
    "desired frame. When averaging over frames you might need to use",
    "the [TT]-nojump[tt] option to obtain the correct average coordinates."
  };
  static bool bMol=FALSE,bCom=FALSE,bNoJump=FALSE;
  static bool bX=TRUE,bY=TRUE,bZ=TRUE,bNorm=FALSE;
  static real scale=0;
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
      "Plot vector length" },
    { "-scale", FALSE, etREAL, {&scale},
      "Scale factor for pdb output, 0 is autoscale" }
    
  };
  FILE       *outx=NULL,*outv=NULL,*outf=NULL,*outb=NULL,*outt=NULL;
  FILE       *outekt=NULL,*outekr=NULL;
  t_topology top;
  real       *mass,time;
  char       title[STRLEN],*indexfn;
  t_trxframe fr;
  int        flags;
  rvec       *xtop,*xp=NULL;
  rvec       *sumxv=NULL,*sumv=NULL,*sumxf=NULL,*sumf=NULL;
  matrix     topbox;
  int        status;
  int        i,j,n;
  int        ngrps,nr_vfr,nr_ffr;
  char       **grpname;
  int        *isize0,*isize;
  atom_id    **index0,**index;
  atom_id    *a,*atndx;
  t_block    *mols;
  bool       bTop,bOX,bOV,bOF,bOB,bOT,bEKT,bEKR,bCV,bCF,bDim[4],bDum[4];
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
    { efXVG, "-ekt","ektrans.xvg", ffOPTWR },
    { efXVG, "-ekr","ekrot.xvg", ffOPTWR },
    { efPDB, "-cv", "veloc.pdb", ffOPTWR },
    { efPDB, "-cf", "force.pdb", ffOPTWR }
  };
#define NFILE asize(fnm)

  CopyRight(stderr,argv[0]);
  parse_common_args(&argc,argv,
		    PCA_CAN_TIME | PCA_TIME_UNIT | PCA_CAN_VIEW | PCA_BE_NICE,
		    NFILE,fnm,asize(pa),pa,asize(desc),desc,0,NULL);

  if (bMol)
    fprintf(stderr,"Interpreting indexfile entries as molecules.\n"
	    "Using center of mass.\n");
  
  bOX  = opt2bSet("-ox",NFILE,fnm);
  bOV  = opt2bSet("-ov",NFILE,fnm);
  bOF  = opt2bSet("-of",NFILE,fnm);
  bOB  = opt2bSet("-ob",NFILE,fnm);
  bOT  = opt2bSet("-ot",NFILE,fnm);
  bEKT = opt2bSet("-ekt",NFILE,fnm);
  bEKR = opt2bSet("-ekr",NFILE,fnm);
  bCV  = opt2bSet("-cv",NFILE,fnm);
  bCF  = opt2bSet("-cf",NFILE,fnm);
  if (bMol || bOT || bEKT || bEKR)
    bCom = TRUE;

  bDim[XX] = bX;
  bDim[YY] = bY;
  bDim[ZZ] = bZ;
  bDim[DIM] = bNorm;

  bTop = read_tps_conf(ftp2fn(efTPS,NFILE,fnm),title,&top,&xtop,NULL,topbox,
		       bCom && (bOX || bOV || bOT || bEKT || bEKR));
  sfree(xtop);
  if ((bMol || bCV || bCF) && !bTop)
    gmx_fatal(FARGS,"Need a run input file for option -mol, -cv or -cf");
  

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
		    xvgr_tlabel(),"Coordinate (nm)");
    make_legend(outx,ngrps,isize0[0],index0[0],grpname,bCom,bMol,bDim);
  }
  if (bOV) {
    flags = flags | TRX_READ_V;
    outv = xvgropen(opt2fn("-ov",NFILE,fnm),
		    bCom ? "Center of mass velocity" : "Velocity",
		    xvgr_tlabel(),"Velocity (nm/ps)");
   make_legend(outv,ngrps,isize0[0],index0[0],grpname,bCom,bMol,bDim); 
  }
  if (bOF) {
    flags = flags | TRX_READ_F;
    outf = xvgropen(opt2fn("-of",NFILE,fnm),"Force",
		    xvgr_tlabel(),"Force (kJ mol\\S-1\\N nm\\S-1\\N)");
    make_legend(outf,ngrps,isize0[0],index0[0],grpname,bCom,bMol,bDim);
  }
  if (bOB) {
    outb = xvgropen(opt2fn("-ob",NFILE,fnm),"Box vector elements",
		    xvgr_tlabel(),"(nm)");
   
    xvgr_legend(outb,6,box_leg);
  }
  if (bOT) {
    bDum[XX] = FALSE;
    bDum[YY] = FALSE;
    bDum[ZZ] = FALSE;
    bDum[DIM] = TRUE;
    flags = flags | TRX_READ_V;
    outt = xvgropen(opt2fn("-ot",NFILE,fnm),"Temperature",xvgr_tlabel(),"(K)");
    make_legend(outt,ngrps,isize[0],index[0],grpname,bCom,bMol,bDum);
  }
  if (bEKT) {
    bDum[XX] = FALSE;
    bDum[YY] = FALSE;
    bDum[ZZ] = FALSE;
    bDum[DIM] = TRUE;
    flags = flags | TRX_READ_V;
    outekt = xvgropen(opt2fn("-ekt",NFILE,fnm),"Center of mass translation",
		      xvgr_tlabel(),"Energy (kJ mol\\S-1\\N)");
    make_legend(outekt,ngrps,isize[0],index[0],grpname,bCom,bMol,bDum);
  }
  if (bEKR) {
    bDum[XX] = FALSE;
    bDum[YY] = FALSE;
    bDum[ZZ] = FALSE;
    bDum[DIM] = TRUE;
    flags = flags | TRX_READ_X | TRX_READ_V;
    outekr = xvgropen(opt2fn("-ekr",NFILE,fnm),"Center of mass rotation",
		      xvgr_tlabel(),"Energy (kJ mol\\S-1\\N)");
    make_legend(outekr,ngrps,isize[0],index[0],grpname,bCom,bMol,bDum);
  }
  if (bCV)
    flags = flags | TRX_READ_X | TRX_READ_V;
  if (bCF)
    flags = flags | TRX_READ_X | TRX_READ_F;
  if (flags == 0 && !bOB) {
    fprintf(stderr,"Please select one or more output file options\n");
    exit(0);
  }

  read_first_frame(&status,ftp2fn(efTRX,NFILE,fnm),&fr,flags);

  if (bCV) {
    snew(sumxv,fr.natoms);
    snew(sumv,fr.natoms);
  }
  if (bCF) {
    snew(sumxf,fr.natoms);
    snew(sumf,fr.natoms);
  }
  nr_vfr = 0;
  nr_ffr = 0;
  
  do {
    time = convert_time(fr.time);

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
      print_data(outx,time,fr.x,mass,bCom,ngrps,isize,index,bDim);
    if (bOV && fr.bV)
      print_data(outv,time,fr.v,mass,bCom,ngrps,isize,index,bDim);
    if (bOF && fr.bF)
      print_data(outf,time,fr.f,NULL,bCom,ngrps,isize,index,bDim);
    if (bOB && fr.bBox)
      fprintf(outb,"\t%g\t%g\t%g\t%g\t%g\t%g\t%g\n",fr.time,
	      fr.box[XX][XX],fr.box[YY][YY],fr.box[ZZ][ZZ],
	      fr.box[YY][XX],fr.box[ZZ][XX],fr.box[ZZ][YY]);
    if (bOT && fr.bV) {
      fprintf(outt," %g",time);
      for(i=0; i<ngrps; i++)
	fprintf(outt,"\t%g",temp(fr.v,mass,isize[i],index[i]));
      fprintf(outt,"\n");
    }
    if (bEKT && fr.bV) {
      fprintf(outekt," %g",time);
      for(i=0; i<ngrps; i++)
	fprintf(outekt,"\t%g",ektrans(fr.v,mass,isize[i],index[i]));
      fprintf(outekt,"\n");
    }
    if (bEKR && fr.bX && fr.bV) {
      fprintf(outekr," %g",time);
      for(i=0; i<ngrps; i++)
	fprintf(outekr,"\t%g",ekrot(fr.x,fr.v,mass,isize[i],index[i]));
      fprintf(outekr,"\n");
    }
    if (bCV && fr.bX && fr.bV) {
      for(i=0; i<fr.natoms; i++)
	rvec_inc(sumxv[i],fr.x[i]);
      for(i=0; i<fr.natoms; i++)
	rvec_inc(sumv[i],fr.v[i]);
      nr_vfr++;
    }
    if (bCF && fr.bX && fr.bF) {
      for(i=0; i<fr.natoms; i++)
	rvec_inc(sumxf[i],fr.x[i]);
      for(i=0; i<fr.natoms; i++)
	rvec_inc(sumf[i],fr.f[i]);
      nr_ffr++;
    }

  } while(read_next_frame(status,&fr));
  

  /* clean up a bit */
  close_trj(status);
  
  if (bOX) fclose(outx);
  if (bOV) fclose(outv);
  if (bOF) fclose(outf);
  if (bOB) fclose(outb);
  if (bOT) fclose(outt);
  if (bEKT) fclose(outekt);
  if (bEKR) fclose(outekr);

  if ((bCV || bCF) && (nr_vfr>1 || nr_ffr>1) && !bNoJump)
    fprintf(stderr,"WARNING: More than one frame was used for option -cv or -cf\n"
	    "If atoms jump across the box you should use the -nojump option\n");

  if (bCV)
    write_pdb_bfac(opt2fn("-cv",NFILE,fnm),"average velocity",&(top.atoms),
		   topbox,isize[0],index[0],nr_vfr,sumxv,sumv,bDim,scale);
  if (bCF)
    write_pdb_bfac(opt2fn("-cf",NFILE,fnm),"average force",&(top.atoms),
		   topbox,isize[0],index[0],nr_ffr,sumxf,sumf,bDim,scale);

  /* view it */
  view_all(NFILE, fnm);
  
  thanx(stderr);
  
  return 0;
}

