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
 * GROningen MAchine for Chemical Simulation
 */
#include "sysstuff.h"
#include "smalloc.h"
#include "statusio.h"
#include "statutil.h"
#include "random.h"
#include "names.h"
#include "vec.h"
#include "futil.h"
#include "copyrite.h"
#include "xvgr.h"
#include "string2.h"
#include "rdgroup.h"

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
	
void do_com(char *trj,char *topf,char *outfX,char *outekrot,
	    int ngrps,int *isize,atom_id **index,char **grpnames)
{
  static char  *axisX[]={ "Xx", "Xy", "Xz", "Xtot" };
  FILE         **outX,*outek;
  t_statheader sh;
  t_topology   *top;
  int          status;
  int          i,j,idum,step,natoms;
  real         t,rdum;
  rvec         *x,*v;
  real         *mass;
  rvec         xcm,acm;
  matrix       L;
  int          g,d,m,n;
  matrix       box;
  atom_id *sysindex;
  
  /* malloc the output files */
  snew(outX,1+ngrps);

  top=read_top(topf);
  natoms=top->atoms.nr;
  snew(mass,natoms);
  snew(x,natoms);
  snew(v,natoms);
  snew(sysindex,natoms);
  for(i=0; (i<top->atoms.nr); i++) {
    mass[i]=top->atoms.atom[i].m;
    sysindex[i]=i;
  }

  outX[0]=xvgropen(outfX,"COM : ","Time(ps)","x (nm)");
  xvgr_legend(outX[0],asize(axisX),axisX);
  for(g=0;(g<ngrps);g++) {
    char filename[STRLEN];

    /* coordinates */
    sprintf(filename,"xcm_%s.xvg",grpnames[g]);
    outX[g+1]=xvgropen(filename,"COM : ","Time(ps)","x (nm)");    
    xvgr_legend(outX[g+1],asize(axisX),axisX);
  }

  outek=xvgropen(outekrot,"EK Rot","Time (ps)","E (kJ/mole)");
  if ((natoms=read_first_x(&status,trj,&t,&x,box))!=top->atoms.nr) 
    fatal_error(0,"Topology (%d atoms) does not match trajectory (%d atoms)",
		top->atoms.nr,natoms);
  
  j=0;
  do {
    if (j>0) {
      for(i=0; (i<natoms); i++)
	for(m=0; (m<DIM); m++) {
	  v[i][m]=(x[i][m]+v[i][m])*2.0;
	}
      fprintf(outek,"%10g  %10g\n",t,calc_ekrot(natoms,mass,x,v));
    }
    fprintf(stderr,"\rFrame: %d",j++);
 
    /* IF COORDINATES ARE PRESENT */
    calc_cm_group(mass,x,xcm,natoms,sysindex);
    fprintf(outX[0],"%10g  %10g  %10g  %10g  %10g\n",
	    t,xcm[XX],xcm[YY],xcm[ZZ],norm(xcm));
    fflush(outX[0]);
    for(g=0;(g<ngrps);g++) {
      calc_cm_group(mass,x,xcm,isize[g],index[g]);
      fprintf(outX[g+1],"%10g  %10g  %10g  %10g  %10g\n",
	      t,xcm[XX],xcm[YY],xcm[ZZ],norm(xcm));
      fflush(outX[g+1]);
    }
    for(i=0; (i<natoms); i++)
      for(m=0; (m<DIM); m++)
	v[i][m]=-x[i][m];
  }  while (read_next_x(status,&t,natoms,x,box));
  close_trj(status);
  fclose(outek);
  
  sfree(x);
  sfree(mass);
  for(g=0;(g<=ngrps);g++) {
    fclose(outX[g]);
  }
}

int main(int argc,char *argv[])
{
  static char *desc[] = {
    "g_com computes the translational and rotational motion of the center",
    "of mass of a molecule as a function of time."
  };
  t_filenm fnm[] = {
    { efTRX,  "-f",  NULL, ffREAD },
    { efTPB,  NULL,  NULL, ffREAD },
    { efNDX,  NULL,  NULL, ffOPTRD },
    { efXVG, "-ox", "xcm", ffWRITE },
    { efXVG, "-oe", "ekrot",ffWRITE }
  };
#define NFILE asize(fnm)

  /* index stuff */
  int      ngrps;       /* the number of groups */
  int      *isize;      /* the size of each group */
  char     **grpnames;  /* the name of each group */
  atom_id  **index;     /* the index array of each group */
  t_block   *block;     /* the total index file */

  int g; /* group counter */

  CopyRight(stderr,argv[0]);
  parse_common_args(&argc,argv,PCA_CAN_TIME,TRUE,NFILE,fnm,
		    0,NULL,asize(desc),desc,0,NULL);

  /* check if ndx file is given */
  if (ftp2bSet(efNDX,NFILE,fnm)) {
    block = (t_block *)init_index(ftp2fn(efNDX,NFILE,fnm),&grpnames);
    ngrps = block->nr;          
    snew(isize,ngrps);
    snew(index,ngrps);
    for(g=0;(g<ngrps);g++) {
      int i;
      isize[g]=block->index[g+1]-block->index[g];
      index[g]=&(block->a[block->index[g]]);
    }
  } else {
    /* no groups means all atoms for vcm and xcm */
    ngrps=0;
  }

  do_com(ftp2fn(efTRX,NFILE,fnm),
	 ftp2fn(efTPB,NFILE,fnm),
	 opt2fn("-ox",NFILE,fnm),
	 opt2fn("-oe",NFILE,fnm),
	 ngrps,isize,index,grpnames);

  thanx(stdout);
  
  return 0;
}
