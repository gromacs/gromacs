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
 * Gyas ROwers Mature At Cryogenic Speed
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
#include "gstat.h"

int main(int argc,char *argv[])
{
  static char *desc[] = {
    "getconf extracts a coordinate frame from a trajectory file",
    "at a user specified time in the simulation.",
    "The periodic boundary conditions may be removed",
    "which means that all molecules are continuous in space",
    "which may be important for displaying programs.",
    "As an extra option coordinates may be in PDB format rather",
    "than gromos which is the default.",
    "An indexfile may be passed to select a subset of the particles",
    "only (eg. a protein from a protein/water system)."
  };
  static char *bugs[] = {
    "The program will crash when the infile is corrupt and no -t option",
    "has been given."
  };
  
  FILE         *in;
  bool         bStop,bPDB;
  t_topology   *top;
  t_statheader sh;
  rvec         *x,*v,*xx,*vv;
  matrix       box;
  int          i,j,n;
  real         t,r;
  int          isize[1];
  atom_id      *index[1];
  char         *grpnames[1];
  t_atoms      atoms,*useatoms;
  
  static real t0=-1;
  static bool remove_pbc=FALSE;
  t_pargs pa[] = {
    { "-t", FALSE, etREAL, &t0, 
      "Time frame of which you want the configuration" },
    { "-x", FALSE, etBOOL, &remove_pbc,
      "Remove periodic boundary conditions" }
  };
  t_filenm fnm[] = {
    { efTRJ, "-f", NULL, ffREAD },
    { efTPX, NULL, NULL, ffREAD },
    { efGRO, "-o", NULL, ffOPTWR },
    { efPDB, "-q", NULL, ffOPTWR },
    { efNDX, NULL, NULL, ffOPTRD }
  };
#define NFILE asize(fnm)

  CopyRight(stderr,argv[0]);
  parse_common_args(&argc,argv,0,FALSE,NFILE,fnm,asize(pa),pa,
		    asize(desc),desc,asize(bugs),bugs);
  
  top=read_top(ftp2fn(efTPX,NFILE,fnm));
    
  snew(x,top->atoms.nr);
  snew(v,top->atoms.nr);
  
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
    }
  } while (!bStop);
  fclose(in);
  fprintf(stderr,"\n");

  if (remove_pbc == TRUE)
    rm_pbc(&(top->idef),top->atoms.nr,box,x,x);
    
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
      
      v[i][XX]=v[index[0][i]][XX];
      v[i][YY]=v[index[0][i]][YY];
      v[i][ZZ]=v[index[0][i]][ZZ];
      
    }
    useatoms=&atoms;
  }
  else
    useatoms=&(top->atoms);
    
  bPDB=ftp2bSet(efPDB,NFILE,fnm);
  if (bPDB) {
    printf("Dumping pdb...\n");
    write_pdb_conf(ftp2fn(efPDB,NFILE,fnm),useatoms,x,box,FALSE);
  }
  if (!bPDB || ftp2bSet(efGRO,NFILE,fnm)) {
    printf("Dumping Gromos...\n");
    write_conf(ftp2fn(efGRO,NFILE,fnm),*(top->name),useatoms,x,v,box);
  }
  thanx(stdout);
  
  return 0;
}
