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
#include "smalloc.h"
#include "math.h"
#include "typedefs.h"
#include "xvgr.h"
#include "copyrite.h"
#include "statutil.h"
#include "statusio.h"
#include "string2.h"
#include "vec.h"
#include "rdgroup.h"
#include "pbc.h"
#include "fatal.h"
#include "futil.h"
#include "do_fit.h"
#include "rmpbc.h"

typedef struct {
  real    dist;
  atom_id ai;
  atom_id aj;
} tclose;

rvec *x;

void check_close(atom_id ai,atom_id aj,matrix box,tclose *close,int nclose)
{
  int this;
  rvec dx;
  real dist;
  bool bCont;
  int n;

  pbc_dx(box,x[ai],x[aj],dx);
  dist = iprod(dx,dx);
  this = nclose - 1;
  
  if ( dist < close[this].dist ) {
    close[this].dist=dist;
    close[this].ai=ai;
    close[this].aj=aj;
    bCont = TRUE;
  }
 
  while ( bCont ) {
    if ( this == 0 ) 
      return; 
    if ( close[this].dist < close[this - 1].dist ) {
      real tdist;
      atom_id ta;
      tdist = close[this].dist;
      close[this].dist = close[this - 1].dist;
      close[this - 1].dist = tdist;
      ta = close[this].ai;
      close[this].ai = close[this - 1].ai;
      close[this - 1].ai = ta;
      ta = close[this].aj;
      close[this].aj = close[this - 1].aj;
      close[this - 1].aj = ta;
      this--;
    } else 
      bCont = FALSE;
  }
}


int main (int argc,char *argv[])
{
  static char *desc[] = {
    "g_close computes the "
  };
  static char *opts[] = {
    "-N"
  };
  static char *odesc[] = {
    "set the number of nearest pairs"
  };
  t_manual man = { asize(desc),desc,asize(opts),opts,odesc,0,NULL};
#define MAXFRAME 10000
  int          step,nre,natom,natoms,i,j,m,teller=0;
  real         t,lambda;
  bool         bTruncOct;
  
  t_statheader header;
  t_inputrec   ir;
  t_topology   top;
  bool         bPBC=TRUE;

  matrix       box;
  rvec         *v;
  int          status;
  char         buf[256];
  
  int          nrgroups;
  char         **groupname;
  int          *isize;
  atom_id      **index;
  int          nclose=1;
  tclose       *close;
  FILE         *fp;
  
  t_filenm fnm[] = {
    { efTPB, NULL, NULL, ffREAD },
    { efTRX, "-f", NULL, ffREAD },
    { efNDX, NULL, NULL, ffREAD },
    { efOUT, "-o", NULL, ffWRITE },
  };
#define NFILE asize(fnm)

  CopyRight(stderr,argv[0]);
  parse_common_args(&argc,argv,PCA_CAN_TIME,NFILE,fnm,TRUE,&man);
  
  for(i=1; (i<argc); i++) {
    if (strcmp(argv[i],"-N") == 0)
      nclose=iscan(argc,argv,&i);
    else
      usage(argv[0],argv[i]);
  }
  snew(close,nclose);

  read_status_header(ftp2fn(efTPB,NFILE,fnm),&header);

  snew(x,header.natoms);
  snew(v,header.natoms);

  read_status(ftp2fn(efTPB,NFILE,fnm),
              &step,&t,&lambda,&ir,
              box,NULL,NULL,
              &natom,
              x,NULL,NULL,&nre,NULL,
              &top);

  /*set box type*/
  if (ir.eBox==ebtTRUNCOCT)
    bTruncOct=TRUE;
  else
    bTruncOct=FALSE;
  init_pbc(box,bTruncOct);


  nrgroups = 2;
  snew(groupname,nrgroups);
  snew(index,nrgroups);
  snew(isize,nrgroups);
 
  rd_index(ftp2fn(efNDX,NFILE,fnm),nrgroups,isize,index,groupname);

  rm_pbc(&(top.idef),top.atoms.nr,box,x,x);


  fp = fopen(ftp2fn(efOUT,NFILE,fnm),"w");

  /* do a first step */
  if ((natoms=read_first_x(&status,ftp2fn(efTRX,NFILE,fnm),&t,&x,box)) 
      != top.atoms.nr) 
    fatal_error(0,"Topology (%d atoms) does not match trajectory (%d atoms)",
		top.atoms.nr,natoms);
  do {
    int i,j,n;
    atom_id ai,aj;
    
    for(n=0;(n<nclose);n++) {
      close[n].ai=0;
      close[n].aj=0;
      close[n].dist=box[0][0]+box[1][1]+box[2][2];
    }
    
    rm_pbc(&(top.idef),top.atoms.nr,box,x,x);
    
    
    /* do everything here */
    for(i=0;(i<isize[0]);i++) {
      ai = index[0][i];
      for(j=0;(j<isize[1]);j++) {
	aj = index[1][j];
	check_close(ai,aj,box,close,nclose);
      }
    }
   
    fprintf(fp,"%8.3f ",t);
    for(n=0;(n<nclose);n++) {
      fprintf(fp,"%8.3f %5d %5d",sqrt(close[n].dist),close[n].ai + 1,close[n].aj + 1);
    }
    fprintf(fp,"\n");
    
    /*print time of frame*/
    if ((teller % 10) == 0)
      fprintf(stderr,"\r %5.2f",t);
      
    teller++;
  } while (read_next_x(status,&t,natom,x,box));
  close_trj(status);
    
  fclose(fp);
  
  thanx(stdout);
  
  return 0;
}
