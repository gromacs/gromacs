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
 * GROningen Mixture of Alchemy and Childrens' Stories
 */
static char *SRCID_g_com_c = "$Id$";

#include "sysstuff.h"
#include "smalloc.h"
#include "macros.h"
#include "statutil.h"
#include "random.h"
#include "names.h"
#include "vec.h"
#include "futil.h"
#include "copyrite.h"
#include "xvgr.h"
#include "string2.h"
#include "rdgroup.h"
#include "tpxio.h"

real calc_ekrot(int natoms,real mass[],rvec x[],rvec v[])
{
  matrix TCM,L;
  real   tm,m0,lxx,lxy,lxz,lyy,lyz,lzz,ekrot;
  rvec   dx,a0,ocm;
  rvec   xcm,vcm,acm;
  int    i,m,n;

  clear_rvec(xcm);
  clear_rvec(vcm);
  clear_rvec(acm);
  tm=0.0;
  for(i=0; (i<natoms); i++) {
    m0=mass[i];
    tm+=m0;
    oprod(x[i],v[i],a0);
    for(m=0; (m<DIM); m++) {
      xcm[m]+=m0*x[i][m]; /* c.o.m. position */
      vcm[m]+=m0*v[i][m]; /* c.o.m. velocity */
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
  for(i=0; (i<natoms); i++) {
    m0=mass[i];
    for(m=0; (m<DIM); m++)
      dx[m]=x[i][m]-xcm[m];
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

real calc_cm_group(real mass[],rvec x[],rvec cm,
		   int isize,atom_id *index)
{
  real tm,m0;
  int  i,m;
  atom_id aid;
  
  clear_rvec(cm);

  tm=0.0;
  for(i=0; (i<isize); i++) {
    aid=index[i];
    m0=mass[aid];
    tm+=m0;
    for(m=0; (m<DIM); m++) {
      cm[m]+=m0*x[aid][m];
    }
  }

  for(m=0; (m<DIM); m++) {
    cm[m]/=tm;
  }

  return tm;
}
	
int main(int argc,char *argv[])
{
  static char *desc[] = {
    "g_com computes the translational and rotational motion ",
    "of a group of atoms (i.e. a protein) as a function of time."
  };
  t_filenm fnm[] = {
    { efTRX,  "-f",  NULL, ffREAD },
    { efTPS,  NULL,  NULL, ffREAD },
    { efNDX,  NULL,  NULL, ffOPTRD },
    { efXVG, "-ox", "xcm", ffWRITE },
    { efXVG, "-oe", "ekrot",ffOPTWR }
  };
#define NFILE asize(fnm)

  static char  *axisX[]={ "Xx", "Xy", "Xz", "Xtot" };
  
  /* index stuff */
  int      ngrps;       /* the number of groups */
  int      *isize;      /* the size of each group */
  char     **grpnames;  /* the name of each group */
  atom_id  **index;     /* the index array of each group */
  t_topology top;
  int      g;           /* group counter */
  char     format[STRLEN],filename[STRLEN],title[STRLEN];
  FILE     **outX,*outek;
  int      status,ftpout;
  int      i,j,idum,step,natoms;
  real     t,rdum;
  rvec     *x,*v;
  real     *mass;
  rvec     xcm,acm;
  matrix   L;
  int      d,m,n;
  matrix   box;
  atom_id  *sysindex;
  bool     bHaveV,bReadV;

  CopyRight(stderr,argv[0]);
  parse_common_args(&argc,argv,PCA_CAN_TIME,TRUE,NFILE,fnm,
		    0,NULL,asize(desc),desc,0,NULL);
  ftpout=fn2ftp(ftp2fn(efTRX,NFILE,fnm));
  bHaveV=((ftpout==efTRJ) || (ftpout==efTRR));
  bReadV=opt2bSet("-oe",NFILE,fnm);
  if ( bReadV && !bHaveV ) {
    fprintf(stderr,"No velocities in input file, "
	    "will not calculate rotational energy\n");
    bReadV=FALSE;
  }
  
  /* open input files, read topology and index */
  read_tps_conf(ftp2fn(efTPS,NFILE,fnm),title,&top,&x,NULL,box,TRUE);
  sfree(x);
  
  fprintf(stderr,"How many groups do you want to calc com of ? ");
  scanf("%d",&ngrps);
  fprintf(stderr,"OK. I will calc com of %d groups\n",ngrps);
  
  snew(grpnames,ngrps);
  snew(index,ngrps);
  snew(isize,ngrps);
  
  get_index(&(top.atoms),ftp2fn_null(efNDX,NFILE,fnm),
	    ngrps,isize,index,grpnames);
  
  if ( bReadV )
    natoms=read_first_x_or_v(&status,ftp2fn(efTRX,NFILE,fnm),&t,&x,&v,box);
  else
    natoms=read_first_x(&status,ftp2fn(efTRX,NFILE,fnm),&t,&x,box);
  if ( natoms > top.atoms.nr )
    fatal_error(0,"Topology (%d atoms) does not match trajectory (%d atoms)",
		top.atoms.nr,natoms);
  
  snew(mass,natoms);
  for(i=0; (i<natoms); i++)
    mass[i]=top.atoms.atom[i].m;
  
  /* open output files */
  snew(outX,ngrps);
  if (ngrps==1) {
    outX[0]=xvgropen(opt2fn("-ox",NFILE,fnm),"COM : ","Time(ps)","x (nm)");
    xvgr_legend(outX[0],asize(axisX),axisX);
  } else {
    strcpy(format,opt2fn("-ox",NFILE,fnm));
    format[strlen(format)-4]='\0';
    strcat(format,"_%s.xvg");
    for(g=0;(g<ngrps);g++) {
      /* coordinates */
      sprintf(filename,format,grpnames[g]);
      outX[g+1]=xvgropen(filename,"COM : ","Time(ps)","x (nm)");
      xvgr_legend(outX[g],asize(axisX),axisX);
    }
  }
  if (bReadV)
    outek=xvgropen(opt2fn("-oe",NFILE,fnm),"EK Rot","Time (ps)","E (kJ/mole)");

  do {
    if (bReadV) {
      bHaveV=FALSE;
      for (i=0; ((i<natoms) && !bHaveV); i++)
	for (d=0; ((d<DIM) && !bHaveV); d++)
	  bHaveV=(v[i][d]!=0);
      if (bHaveV)
	fprintf(outek,"%10g  %10g\n",t,calc_ekrot(natoms,mass,x,v));
    }
 
    /* IF COORDINATES ARE PRESENT */
    for(g=0;(g<ngrps);g++) {
      calc_cm_group(mass,x,xcm,isize[g],index[g]);
      for (j = 0; j < DIM; j++) {
	while (xcm[j] < 0) xcm[j]+= box[j][j];
	while (xcm[j] > box[j][j]) xcm[j]-=box[j][j];
      }
      fprintf(outX[g],"%10g  %10g  %10g  %10g  %10g\n",
	      t,xcm[XX],xcm[YY],xcm[ZZ],norm(xcm));
    }
    if (bReadV)
      for(i=0; (i<natoms); i++)
	for(m=0; (m<DIM); m++)
	  v[i][m]=-x[i][m];
  } while ((bReadV && read_next_x_or_v(status,&t,natoms,x,v,box)) ||
	   (!bReadV && read_next_x(status,&t,natoms,x,box)));
  sfree(x);
  if (bReadV)
    sfree(v);
  sfree(mass);
  
  close_trj(status);
  if (bReadV)
    fclose(outek);
  for(g=0;(g<ngrps);g++) {
    fclose(outX[g]);
  }
  
  thanx(stdout);
  
  return 0;
}

