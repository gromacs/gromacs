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
 * Grunge ROck MAChoS
 */
#include "sysstuff.h"
#include "smalloc.h"
#include "typedefs.h"
#include "copyrite.h"
#include "statutil.h"
#include "futil.h"
#include "confio.h"
#include "pbc.h"
#include "rdgroup.h"
#include "macros.h"
#include "rmpbc.h"
#include <vec.h>
#include <vdw.h>
#include <string2.h>

real distance2(rvec v1,rvec v2,matrix box)
{
  rvec dx;
  pbc_dx(box,v1,v2,dx);
  return ( iprod(dx,dx) );  
}


static real min_distance(char resname[],char atomname[],
			 t_vdw **vdw,int *max_i,real r_distance)
{
  int i,match=0;
  real distance=0;
  for (i=0; ((i<*max_i) && (match < 1)); i++ )
    if ( strcmp(atomname,(*vdw)[i].atomname) == 0){
      distance = (*vdw)[i].distance;
      match=1;
    } 
  if (match==0) {
    (*max_i)++;
    srenew(*vdw,*max_i);
    (*vdw)[*max_i-1].distance=r_distance;
    strcpy((*vdw)[*max_i-1].atomname,atomname);
    distance = (*vdw)[*max_i-1].distance;
    fprintf(stderr,"distance of %s %s is set to %f\n",
	    resname,(*vdw)[*max_i-1].atomname,
	    (*vdw)[*max_i-1].distance);
  }
  return distance;
} /*min_distance()*/


static void mk_vdw(t_atoms *a,
		   real rvdw[],t_vdw **vdw,int *max_vdw,
		   real r_distance)
{
  int i;
  
  /*initialise van der waals arrays of configuration */
  fprintf(stderr,"Initialising van der waals distances...\n");
  for(i=0; (i < a->nr); i++)
    rvdw[i]=min_distance(*(a->resname[a->atom[i].resnr]),
			 *(a->atomname[i]),
			 vdw,max_vdw,r_distance);
}



int main(int argc,char *argv[])
{
  static char *desc[] = {
    "g_probesurface is a secret."
  };
  static char *opts[] = {
    "-T Time",
    "-R Probe Radius",
    "-S Grid spacing"
  };
  static char *odesc[] = {
    "How",
    "Are",
    "You?"
  };
  t_manual man = { asize(desc),desc,asize(opts),opts,odesc,0,NULL};

  /* variables checked */
  int natoms;            /* number of atoms */
  int d;                 /* dimension counter */
  rvec xmin,xmax;        /* minimum and maximum of the surface */
  real radius=0.7;       /* probe radius */
  real radius2=radius*radius; /* squared probe radius */
  real tradius=radius + 0.25;  /* needed for quick scan */
  real tradius2=tradius*tradius; /* squared probe tradius for quick scan */
  rvec pstart,pstop;     /* start coordinates of the probe */
  real stepsize=0.1;     /* stepsize (nm) of the probe sphere */
  int  nx,ny,nz;         /* probe counters */
  int  nsteps[DIM];      /* number of steps for probe */
  rvec probe;            /* prbe coordinates */
  real **probemin;       /* lower side of the surface */
  real **probemax;       /* upper side of the surface */
  FILE *fp;              /* file pointer to output files */
  char filename[255];    /* output filename */
  
  int countmin;
  int countmax;

  /* van der waals data */
  int    max_vdw;
  t_vdw  *vdw;
  real   *rvdw;
  real vdwdist;

  /* variables */

  FILE         *in;
  bool         bStop;
  t_topology   *top;
  t_statheader sh;
  rvec         *x,*v,*xx,*vv;
  matrix       box;
  int          i,j,n;
  real         t,t0=0,r;
  int          time=0;
  int          isize[1];
  atom_id      *index[1];
  char         *grpnames[1];
  t_atoms      atoms,*useatoms;
  t_filenm fnm[] = {
    { efTRJ,  "-f",       NULL, ffREAD},
    { efTPX,  NULL,       NULL, ffREAD},
    { efNDX,  NULL,       NULL, ffOPTRD},
    { efVDW, "-wi",       NULL, ffOPTRD},
    { efVDW, "-wo", "outradii", ffWRITE}
  };
#define NFILE asize(fnm)

  CopyRight(stderr,argv[0]);
  parse_common_args(&argc,argv,0,NFILE,fnm,FALSE,&man);
  
  t0=-1;
  for(i=1; (i<argc); i++) {
    if (argv[i][0] != '-')
      usage(argv[0],argv[i]);
    switch(argv[i][1]) {
    case 'S':
      stepsize=dscan(argc,argv,&i);
      break;
    case 'R':
      radius=dscan(argc,argv,&i);
      tradius=radius + 0.25;
      tradius2=tradius*tradius;
      radius2 = radius * radius;
      break;
    case 'T':
      t0=dscan(argc,argv,&i);
      time = (int)t0;
      break;
    default:
      usage(argv[0],argv[i]);
    }
  } 
  
  /*read van der waals distances for all the existing atoms*/
  max_vdw=read_vdw(opt2fn("-wi",NFILE,fnm),&vdw);


  top=read_top(ftp2fn(efTPX,NFILE,fnm));
    
  snew(x,top->atoms.nr);
  snew(v,top->atoms.nr);
  snew(rvdw,top->atoms.nr);

  in=ftp2FILE(efTRJ,NFILE,fnm,"r");
  
  bStop=FALSE;
  j=0;
  do {
    fprintf(stderr,"\rFrame: %d",j++);
    bStop=eof(in);
    if (!bStop) {
      rd_header(in,&sh);
      xx=NULL;
      vv=NULL;
      if (sh.x_size)
	xx=x;
      else
	xx=NULL;
      if (sh.v_size)
	vv=v;
      else
	vv=NULL;
      
      rd_hstatus(in,&sh,&n,&t,
		 &r,NULL,box,NULL,NULL,&n,
		 xx,vv,NULL,&n,NULL,NULL);
      if (t0!=-1)
	bStop=(t>=t0);
      time = t;
    }
  } while (!bStop);
  fclose(in);
  fprintf(stderr,"\n");

    
  if (ftp2bSet(efNDX,NFILE,fnm)) {
    snew(atoms.atom,top->atoms.nr);
    snew(atoms.atomname,top->atoms.nr);
    
    /* Copy a pointer... */
    atoms.resname=top->atoms.resname;
    
    rd_index(ftp2fn(efNDX,NFILE,fnm),1,isize,index,grpnames);
    atoms.nr=isize[0];
    atoms.nres=top->atoms.nres;
    for(i=0;(i<isize[0]);i++) {
      atoms.atomname[i]=top->atoms.atomname[index[0][i]];
      atoms.atom[i].resnr=top->atoms.atom[index[0][i]].resnr;
      x[i][XX]=x[index[0][i]][XX];
      x[i][YY]=x[index[0][i]][YY];
      x[i][ZZ]=x[index[0][i]][ZZ];
    }
    useatoms=&atoms;
  }
  else
    useatoms=&(top->atoms);
    


  natoms = atoms.nr;
  /* find the initialisation values for x */
  for(d=0;(d<DIM);d++) {
    xmin[d]=x[0][d];
    xmax[d]=x[0][d];
  }
  for(n=1;(n<natoms);n++) 
    for(d=0;(d<DIM);d++) {
      if ( x[n][d] < xmin[d] )
	xmin[d]=x[n][d];
      if( x[n][d] > xmax[d] )
	xmax[d]=x[n][d];
    }
  
  for(d=0;(d<DIM);d++) {
    pstart[d]=xmin[d]-radius*1.01;
    pstop[d]=xmax[d]+radius*1.01;
    fprintf(stderr,"DIM: %5d PSTART: %8.3f PSTOP: %8.3f\n",d,pstart[d],pstop[d]); 
  }
  
  /* determione the number of steps */
  for(d=0;(d<DIM);d++) 
    nsteps[d]=(int)box[d][d]/stepsize;

  fprintf(stderr,"NSTEPS[XX] = %5d\n",nsteps[XX]);
  fprintf(stderr,"NSTEPS[YY] = %5d\n",nsteps[YY]);
  fprintf(stderr,"NSTEPS[ZZ] = %5d\n",nsteps[ZZ]);
  
  /* allocate memory for the output arrays */
  snew(probemin,nsteps[YY]);
  snew(probemax,nsteps[YY]);
  for(ny=0;(ny<nsteps[YY]);ny++) {
    snew(probemin[ny],nsteps[ZZ]);
    snew(probemax[ny],nsteps[ZZ]);
  }
  
  /* load van der waals data */
  mk_vdw(&(top->atoms),rvdw,&vdw,&max_vdw,0.1);
  

  /* init pbc */
  init_pbc(box,FALSE);

  /* let a probe initialise on several yz coordinates*/
  for(ny=0;(ny<nsteps[YY]);ny++) {
    fprintf(stderr,"\r%5d    ",ny);
    for(nz=0;(nz<nsteps[ZZ]);nz++) {
      probe[YY] = box[YY][YY] * (real)ny / (real)nsteps[YY];
      probe[ZZ] = box[ZZ][ZZ] * (real)nz / (real)nsteps[ZZ];

      probemin[ny][nz]=pstop[XX];
      probemax[ny][nz]=pstart[XX];
      
      /* let probe walk through the surface */
      for(probe[XX]=pstart[XX];(probe[XX]<pstop[XX]);probe[XX]+=stepsize) {
	static real mindist;

	/* find the nearest atom */
	mindist=distance2(probe,x[0],box);
	for(n=1;(n<natoms);n++) {
	  /* calculate the vander waals distance */
	  vdwdist = rvdw[n];

	  if (fabs(x[n][YY]-probe[YY])<tradius)
	    if (fabs(x[n][ZZ]-probe[ZZ])<tradius)
	      if (fabs(x[n][XX]-probe[XX])<tradius) {
		static real tdist;
		tdist = distance2(probe,x[n],box) - vdwdist;
		if ( tdist < mindist )
		  mindist = tdist;
	      }
	}
	
	if ( mindist < radius2 ) {
	  if ( probe[XX] > probemax[ny][nz] )
	    probemax[ny][nz]=probe[XX];
	  if ( probe[XX] < probemin[ny][nz] )
	    probemin[ny][nz]=probe[XX];
	}

	
	
      }
    }
  }


  /* calculate the area which is still in the minimum or maximum */
  countmin = 0;
  countmax = 0;
  for(ny=0;(ny<nsteps[YY]);ny++) {
     for(nz=0;(nz<nsteps[ZZ]);nz++) {
       if ( probemin[ny][nz] == pstop[XX] )
	 countmin++;
       if ( probemax[ny][nz] == pstart[XX] )
	 countmax++;
     }
  }
  fprintf(stderr,"countmin = %5d %8.3f\n",countmin , countmin * stepsize * stepsize );
  fprintf(stderr,"countmax = %5d %8.3f\n",countmax , countmax * stepsize * stepsize );
  

  /* print maxima */
  fprintf(stderr,"writing maxima\n");
  sprintf(filename,"probemax_%d.gnp",time);
  fp = ffopen(filename,"w");
  for(ny=0;(ny<nsteps[YY]);ny++) {
    for(nz=0;(nz<nsteps[ZZ]);nz++) 
      fprintf(fp,"%8.3f\n",probemax[ny][nz]);
    fprintf(fp,"\n");
  }
  fflush(fp);
  fclose(fp);



  /* print minima */
  fprintf(stderr,"writing minima\n");
  sprintf(filename,"probemin_%d.gnp",time);
  fp = ffopen(filename,"w");
    for(ny=0;(ny<nsteps[YY]);ny++) {
    for(nz=0;(nz<nsteps[ZZ]);nz++) 
      fprintf(fp,"%8.3f\n",probemin[ny][nz]);
    fprintf(fp,"\n");
  }
  fflush(fp);
  fclose(fp);

  /* write van der waals output */
  write_vdw(opt2fn("-wo",NFILE,fnm),vdw,max_vdw);

  
}







