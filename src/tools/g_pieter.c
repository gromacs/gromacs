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
#include <typedefs.h>
#include "smalloc.h"
#include "math.h"
#include "xvgr.h"
#include "copyrite.h"
#include "statutil.h"
#include "string2.h"
#include "vec.h"
#include "rdgroup.h"
#include "pbc.h"
#include "fatal.h"
#include "futil.h"
#include "gstat.h"

void main(int argc,char *argv[])
{
  static char *desc[] = {
    "g_rms computes the root mean square deviation (RMSD) of a structure",
    "from a trajectory x(t) with respect to a reference structure from a",
    "topology x(ref) by LSQ fitting the structures on top of each other.",
    "The reference structure is taken from the binary topology file.[PAR]",
    "Option -a produces time averaged RMSD per group (e.g. residues in",
    "a protein)"
  };
  static char *opts[] = {
    "-k",
    "-nf",
    "-nfc",
    "-C"
  };
  static char *odesc[] = {
    "no PBC check",
    "don't fit to reference structure, just center on the origin",
    "don't fit to reference structure and don't center trajectory on origin",
    "calc average RMS for each element of normal RMS group (max 1 group)"
  };
  t_manual man = { asize(desc),desc,asize(opts),opts,odesc,0,NULL};

  t_topology *top;
  real t;
  rvec *x;
  matrix box;
  int status;
  int natoms;

  int i,j,teller;
  atom_id aid;
  int score;

  bool **bOcc;   /* 2-D boolean array of the */
  real gsize = 0.1;    /* the size of the grid */
  int  nx,ny,nz;    /* number of cells in two dimensions */

  int     ngrps;     /* the number of index groups */
  atom_id **index;   /* the index for the atom numbers */
  int     *isize;    /* the size of each group */
  char    **grpname; /* the name of each group */
  int step,s;
  int surf;
  bool *bRes;
  int nres;
  int *nindex;
  
  t_filenm fnm[] = {
    { efTPX, NULL, NULL, ffREAD },
    { efTRX, "-f", NULL, ffREAD },
    { efNDX, NULL, NULL, ffREAD },
  };
#define NFILE asize(fnm)


  CopyRight(stdout,argv[0]);

  parse_common_args(&argc,argv,PCA_CAN_TIME,NFILE,fnm,TRUE,&man);
  
  for(i=1; (i<argc); i++) {
    if (strcmp(argv[i],"-g") == 0)
      gsize=dscan(argc,argv,&i);
    else
      usage(argv[0],argv[i]);
  }
  

  top=read_top(ftp2fn(efTPX,NFILE,fnm));
  
  natoms=read_first_x(&status,ftp2fn(efTRX,NFILE,fnm),&t,&x,box);
  
  
  /* read index files */
  ngrps = 2;
  snew(grpname,ngrps);
  snew(index,ngrps);
  snew(isize,ngrps);
  rd_index(ftp2fn(efNDX,NFILE,fnm),ngrps,isize,index,grpname);
  

  /* find out how many residues are in group 1 */
  snew(bRes,top->atoms.nres);
  for(i=0;(i<top->atoms.nres);i++) 
    bRes[i]=FALSE;

  for(i=0;(i<isize[1]);i++)
    bRes[top->atoms.atom[index[1][i]].resnr]=TRUE;

  nres = 0;
  for(i=0;(i<top->atoms.nres);i++)
    if ( bRes[i] )
      nres++;
  fprintf(stderr,"nres = %5d\n",nres);
  
  teller = 0;
  
  /* initialise grid */
  nz = (int)(box[ZZ][ZZ] / gsize) + 100;
  ny = (int)(box[YY][YY] / gsize) + 100;
  
  snew(bOcc,nz);
  for(i=0;(i<nz);i++) {
    snew(bOcc[i],ny);
    for(j=0;(j<ny);j++) {
      bOcc[i][j] = FALSE;
    }
  }

  do {
    if ((teller % 10) == 0)
      fprintf(stderr,"\rFrame %d",teller);


    /* initialise grid with group one */
    for(i=0;(i<isize[0]);i++) {
      aid = index[0][i];
      if ( x[aid][ZZ] < 0 )
	x[aid][ZZ] += box[ZZ][ZZ];
      if ( x[aid][YY] < 0 )
	x[aid][YY] += box[YY][YY];
      bOcc[(int)(x[aid][ZZ]/gsize)][(int)(x[aid][YY]/gsize)] = TRUE;
    }
   
    /* fill holes */
    for(s=0;(s<5);s++) 
      for(i=1;(i<nz-1);i++)
	for(j=1;(j<nz-1);j++) 
	  if ( bOcc[i][j] == FALSE ) {
	    int yes = 0;
	    if ( bOcc[i-1][j-1] )
	      yes++;
	    if ( bOcc[i-1][j] )
	      yes++;
	    if ( bOcc[i-1][j+1] )
	      yes++;
	    if ( bOcc[i][j-1] )
	      yes++;
	    if ( bOcc[i][j+1] )
	      yes++;
	    if ( bOcc[i+1][j-1] )
	      yes++;
	    if ( bOcc[i+1][j] )
	      yes++;
	    if ( bOcc[i+1][j+1] )
	      yes++;
	    if ( yes > 4 )
	      bOcc[i][j] = TRUE;
	  }
    

    printf("%8.3f ",t);

    score = 0;

    surf = 0;
    for(i=0;(i<nz);i++)
      for(j=0;(j<ny);j++)
	if ( bOcc[i][j] )
	  surf++;
    printf("%8.3f ", surf * gsize * gsize);
    
    for(i=0;(i<top->atoms.nres);i++) 
      bRes[i]=FALSE;
  
    for(i=0;(i<isize[1]);i++) {
      aid = index[1][i];
      if ( x[aid][ZZ] < 0 ) 
	x[aid][ZZ] += box[ZZ][ZZ];
      if ( x[aid][YY] < 0 ) 
	x[aid][YY] += box[YY][YY];
      if ( bOcc[(int)(x[aid][ZZ]/gsize)][(int)(x[aid][YY]/gsize)] == TRUE ) {
	score++;
	bRes[top->atoms.atom[aid].resnr]=TRUE;
      }
    }
    printf("%5d ",score); 
    
    score = 0;
    for(i=0;(i<top->atoms.nres);i++) 
      if ( bRes[i] )
	score++;
    printf("%5d ",score); 
    
    score = 0;
    
    
    printf("\n");
    fflush(stdout);
    
    teller++;
  } while (read_next_x(status,&t,natoms,x,box));
  fprintf(stderr,"\n");
  close_trj(status);

}
