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
#include "gstat.h"

/* calculates center of mass in xcm and returns total mass */

int main(int argc,char *argv[])
{
  static char *desc[] = {
    "g_coord computes the coordinates of atoms or the center of mass of",
    "molecules as function of time. It takes an index file with molecule",
    "number and generates three output files with resp. the x-, y- and z-",
    "coordinate for each entry in the index file"
  };
  static bool bMolecular = FALSE;  /* default is to use atoms    */
  t_pargs pa[] = {
    { "-m", FALSE, etBOOL, &bMolecular,
      "If given, index is interpreted as containing molecule, not atom numbers"}
  };
  FILE    *xout,*yout,*zout; /* xvgr files with coordinates*/
  t_topology *top;           /* topology                   */
  rvec       *x,*x_s;        /* coord, with and without pbc*/
  rvec       xcm;            /* center of mass of molecule */
  matrix     box;            /* box matrix (3x3)           */
  real       t,tm;           /* time, total mass molecule  */
  int        natoms;         /* number of atoms in system  */
  int        status;
  int        i,j,k=0;        /* loopcounters               */
  char       *grpname;       /* name of the group          */
  int        gnx;            /* number of elements in group*/
  int        moleculesize;   /* size of molecule in numbers*/
  atom_id    *index;         /* molecule numbers in index  */
  atom_id    *a, *atndx;     /* indices from block.        */
  t_block *mols;             /* all molecules in system    */
  t_filenm fnm[] = {
    { efTRX, "-f", NULL, ffREAD },
    { efTPB, NULL, NULL, ffREAD },
    { efXVG, "-ox", "x-coor.xvg", ffWRITE },
    { efXVG, "-oy", "y-coor.xvg", ffWRITE },
    { efXVG, "-oz", "z-coor.xvg", ffWRITE },
    { efNDX, NULL, NULL, ffREAD }
  };
#define NFILE asize(fnm)

  CopyRight(stderr,argv[0]);
  parse_common_args(&argc,argv,PCA_CAN_TIME | PCA_CAN_VIEW,TRUE,
		    NFILE,fnm,asize(pa),pa,asize(desc),desc,0,NULL);

  if (bMolecular)
    fprintf(stderr,"Interpreting indexfile entries as molecules."
	    "Using center of mass.\n");
  
  rd_index(ftp2fn(efNDX,NFILE,fnm),1,&gnx,&index,&grpname);
  top=read_top(ftp2fn(efTPB,NFILE,fnm));
  
  mols=&(top->blocks[ebMOLS]);
  a = mols->a;
  atndx = mols->index;

  natoms=read_first_x(&status,ftp2fn(efTRX,NFILE,fnm),&t,&x,box);
  snew(x_s,natoms);

  if (bMolecular)
  {
    xout=xvgropen(opt2fn("-ox",NFILE,fnm),
		  "X-coordinates of center of mass","Time (ps)","Coord (nm)");
    yout=xvgropen(opt2fn("-oy",NFILE,fnm),
		  "Y-coordinates of center of mass","Time (ps)","Coord (nm)");
    zout=xvgropen(opt2fn("-oz",NFILE,fnm),
		  "Z-coordinates of center of mass","Time (ps)","Coord (nm)");
  } else {
    xout=xvgropen(opt2fn("-ox",NFILE,fnm),
		  "X-coordinates of atoms","Time (ps)","Coord (nm)");
    yout=xvgropen(opt2fn("-oy",NFILE,fnm),
		  "Y-coordinates of atoms","Time (ps)","Coord (nm)");
    zout=xvgropen(opt2fn("-oz",NFILE,fnm),
		  "Z-coordinates of atoms","Time (ps)","Coord (nm)");
  };

  do {
    fprintf(xout,"%10g  ", t);
    fprintf(yout,"%10g  ", t);
    fprintf(zout,"%10g  ", t);

    if (bMolecular)
    {
      rm_pbc(&(top->idef),top->atoms.nr,box,x,x_s); /* remove pbc. */
    
      /* Loop over all molecules. Calculate the center of mass for each
	 molecule. To do so, give an index with all atoms in the molecule
	 to calc_xcm. Index number for atom j of molecule i is 
	 a[atndx[index[i]]+j]. Don't you just love Gromacs? See block.h
	 for the definition of t_block and how to use it. 
	 */
      for (i = 0; i < gnx; i++)
      {
	moleculesize = atndx[index[i]+1] - atndx[index[i]];
	tm=calc_xcm(x_s, moleculesize, &a[atndx[index[i]]], 
		  top->atoms.atom, xcm, FALSE);
	
	/* We used non-pbc coordinates. Now put cm back in the box */
	for (j = 0; j < 3; j++)
	{
	  if (xcm[j] < 0) xcm[j]+= box[j][j];
	  if (xcm[j] > box[j][j]) xcm[j]-=box[j][j];
	}
	fprintf(xout,"%10g\t", xcm[0]);
	fprintf(yout,"%10g\t", xcm[1]);
	fprintf(zout,"%10g\t", xcm[2]);
      }
      /* End loop over all molecules */
      fprintf(xout, "\n");
      fprintf(yout, "\n");
      fprintf(zout, "\n");

    } else {
      /* We are not interested in molecules, just in the coordinates of atoms */
      for(i=0; i<gnx; i++)
      {
	fprintf(xout,"%10g\t", x[index[i]][0]);
	fprintf(yout,"%10g\t", x[index[i]][1]);
	fprintf(zout,"%10g\t", x[index[i]][2]);
      }
      /* End loop over all atoms */
      fprintf(xout, "\n");
      fprintf(yout, "\n");
      fprintf(zout, "\n");
    } 
     
    if ((k % 10) == 0)
      fprintf(stderr,"\rt = %6.1f",t);
    k++;
  } while(read_next_x(status,&t,natoms,x,box));
  

  /* clean up a bit */
  fprintf(stderr,"\n");
  close_trj(status);
  fclose(xout);
  fclose(yout);
  fclose(zout);

  xvgr_file(opt2fn("-ox",NFILE,fnm), NULL);
  xvgr_file(opt2fn("-oy",NFILE,fnm), NULL);
  xvgr_file(opt2fn("-oz",NFILE,fnm), NULL);
  
  thanx(stdout);
  
  return 0;
}

