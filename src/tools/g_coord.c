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
static char *SRCID_g_coord_c = "$Id$";

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
    "molecules as function of time. It takes an index file with atom, ",
    "or optionally molecule, numbers and generates an output file with the ",
    "coordinates for each entry in the index file"
  };
  static bool bMolecular = FALSE;  /* default is to use atoms    */
  t_pargs pa[] = {
    { "-m", FALSE, etBOOL, &bMolecular,
      "index contains molecule numbers i.s.o. atom"}
  };
  FILE       *out;           /* xvgr file with coordinates */
  t_topology *top;           /* topology                   */
  rvec       *x,*x_s;        /* coord, with and without pbc*/
  rvec       xcm;            /* center of mass of molecule */
  matrix     box;            /* box matrix (3x3)           */
  real       t;              /* time, total mass molecule  */
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
    { efTPX, NULL, NULL, ffREAD },
    { efXVG, "-o", "coord.xvg", ffWRITE },
    { efNDX, NULL, NULL, ffREAD }
  };
#define NFILE asize(fnm)

  CopyRight(stderr,argv[0]);
  parse_common_args(&argc,argv,PCA_CAN_TIME | PCA_CAN_VIEW,TRUE,
		    NFILE,fnm,asize(pa),pa,asize(desc),desc,0,NULL);

  if (bMolecular)
    fprintf(stderr,"Interpreting indexfile entries as molecules."
	    "Using center of mass.\n");
  
  top=read_top(ftp2fn(efTPX,NFILE,fnm));
  get_index(&(top->atoms),ftp2fn_null(efNDX,NFILE,fnm),
	    1,&gnx,&index,&grpname);
  
  mols=&(top->blocks[ebMOLS]);
  a = mols->a;
  atndx = mols->index;

  natoms=read_first_x(&status,ftp2fn(efTRX,NFILE,fnm),&t,&x,box);
  snew(x_s,natoms);

  out=xvgropen(opt2fn("-o",NFILE,fnm),
	       bMolecular ? "Center of mass" : "Coordinate",
	       "Time (ps)","Coord (nm)");

  do {
    fprintf(out,"%10g", t);

    if (bMolecular) {
      rm_pbc(&(top->idef),top->atoms.nr,box,x,x_s); /* remove pbc. */
    
      /* Loop over all molecules. Calculate the center of mass for each
	 molecule. To do so, give an index with all atoms in the molecule
	 to calc_xcm. Index number for atom j of molecule i is 
	 a[atndx[index[i]]+j]. Don't you just love Gromacs? See block.h
	 for the definition of t_block and how to use it. 
	 */
      for (i = 0; i < gnx; i++) {
	moleculesize = atndx[index[i]+1] - atndx[index[i]];
	calc_xcm(x_s, moleculesize, &a[atndx[index[i]]], 
		 top->atoms.atom, xcm, FALSE);
	
	/* We used non-pbc coordinates. Now put cm back in the box */
	for (j = 0; j < DIM; j++) {
	  if (xcm[j] < 0) xcm[j]+= box[j][j];
	  if (xcm[j] > box[j][j]) xcm[j]-=box[j][j];
	}
	fprintf(out,"\t%10g\t%10g\t%10g", xcm[XX], xcm[YY], xcm[ZZ]);
      }
      /* End loop over all molecules */
      fprintf(out, "\n");

    } else {
      /* We are not interested in molecules, just in the coordinates of atoms */
      for(i=0; i<gnx; i++) {
	fprintf(out,"\t%10g\t%10g\t%10g", 
		x[index[i]][0], x[index[i]][1], x[index[i]][2]);
      }
      /* End loop over all atoms */
      fprintf(out, "\n");
    } 
     
    k++;
  } while(read_next_x(status,&t,natoms,x,box));
  

  /* clean up a bit */
  fprintf(stderr,"\n");
  close_trj(status);
  fclose(out);

  xvgr_file(opt2fn("-o",NFILE,fnm), NULL);
  
  thanx(stdout);
  
  return 0;
}

