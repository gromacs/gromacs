/*
 *       $Id$
 *
 *       This source code is part of
 *
 *        G   R   O   M   A   C   S
 *
 * GROningen MAchine for Chemical Simulations
 *
 *            VERSION 1.6
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
 * Great Red Owns Many ACres of Sand 
 */
static char *SRCID_gro2pdb_c = "$Id$";

#include <math.h>
#include "smalloc.h"
#include "typedefs.h"
#include "macros.h"
#include "string2.h"
#include "confio.h"
#include "vec.h"
#include "statutil.h"
#include "copyrite.h"
#include "pdbio.h"

int main(int argc, char *argv[])
{
  static char *desc[] = {
    "gro2pdb converts a gromos coordinate file to a pdb (protein databank)",
    "type file, suitable for displaying using e.g. rasmol, reads B-factor",
    "data from file with with following format: first line states number of",
    "residues or atoms in the file, next lines state either residue or atom",
    "number followed by B-factor. Obiously, any type of numeric data can",
    "be displayed in stead of B-factors."
  };

  int n,n_bfac;
  
  char **bfac_lines;
  int *bfac_nr;
  double *bfac_val,bfac_max,bfac_min;
  bool doBfact;
  
  FILE    *out;
  char    title[STRLEN];
  int     natoms,i;
  rvec    *x,*v;
  char    ch;
  matrix  box;
  t_atoms atoms;
  t_pdbatom *ptr;
  t_filenm fnm[] = {
    { efGRO, "-f", NULL, ffREAD },
    { efPDB, "-o", NULL, ffWRITE },
    { efDAT, "-bf", "bfact", ffOPTRD }
  };
#define NFILE asize(fnm)
  static char *chbuf=NULL;
  static bool peratom=FALSE;
  t_pargs pa[] = {
    { "-chain", FALSE, etSTR, &chbuf, "Chain name" },
    { "-atom", FALSE, etBOOL, &peratom, "B-factors in input file per atom (default per residue)" }
  };
  
  CopyRight(stderr,argv[0]);
  parse_common_args(&argc,argv,0,FALSE,NFILE,fnm,asize(pa),pa,
		    asize(desc),desc,0,NULL);
  if (chbuf)
    ch = chbuf[0];
  else
    ch=' ';

  doBfact = opt2bSet("-bf",NFILE,fnm);
  
  atoms.nr=0;
  atoms.nres=0;
  get_coordnum(fnm[0].fn,&natoms);
  snew(atoms.atomname,natoms);
  snew(atoms.resname,natoms);
  snew(atoms.atom,natoms);
  snew(x,natoms);
  snew(v,natoms);
  read_whole_conf(ftp2fn(efGRO,NFILE,fnm),title,&atoms,x,v,box);
  ptr=atoms2pdba(&atoms,x);  
  pdba_trimnames(natoms,ptr);

  if (doBfact) {  
    n_bfac = get_lines(ftp2fn(efDAT,NFILE,fnm), &bfac_lines);
    snew(bfac_val, n_bfac);
    snew(bfac_nr, n_bfac);
    fprintf(stderr, "Reading %d B-factors from %s\n",n_bfac,ftp2fn(efDAT,NFILE,fnm));
    bfac_max=-1e10;
    bfac_min=1e10;
    for(i=0; (i<n_bfac); i++) {
      /*fprintf(stderr, "Line %d: %s\n",i,bfac_lines[i]);*/
      sscanf(bfac_lines[i],"%d %lf",&bfac_nr[i],&bfac_val[i]);
      if (bfac_val[i] > bfac_max) 
	bfac_max = bfac_val[i];
      if (bfac_val[i] < bfac_min) 
	bfac_min = bfac_val[i];
    }
    while ( (bfac_max > 99.99) || (bfac_min < -99.99) ) {
      fprintf(stderr,"Range of values for B-factors to large (min %g, max %g) "
	      "will scale down a factor 10\n",bfac_min,bfac_max);
      for(i=0; (i<n_bfac); i++)
	bfac_val[i] /= 10;
      bfac_max /= 10;
      bfac_min /= 10;
    }
    
    for(i=0; (i<natoms); i++)
      ptr[i].bfac=0;
    
    if (peratom) {
      if ( n_bfac != natoms )
	warning("Number of B-factors in inputfile does not match number of atoms.");
      for(i=0; (i<n_bfac); i++)
	ptr[bfac_nr[i]].bfac=bfac_val[i];
    } else {
      if ( n_bfac != atoms.nres ) 
	warning("Number of B-factors in inputfile does not match number of residues.");
      for(i=0; (i<n_bfac); i++)
	for(n=0; (n<natoms); n++)
	  if ( bfac_nr[i] == (ptr[n].resnr+1) )
	    ptr[n].bfac=bfac_val[i];
    }
  }
  for(i=0; (i<natoms); i++)
    ptr[i].chain=ch;
  
  out=ftp2FILE(efPDB,NFILE,fnm,"w");
  print_pdbatoms(out,natoms,ptr,box);
  fclose(out);
  
  thanx(stdout);

  return 0;
}

