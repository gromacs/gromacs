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
 * Green Red Orange Magenta Azure Cyan Skyblue
 */
static char *SRCID_gen_xh_c = "$Id$";

#include "macros.h"
#include "smalloc.h"
#include "vec.h"
#include "confio.h"
#include "statutil.h"
#include "copyrite.h"

#define AA 0.081649
#define BB 0.0
#define CC 0.0577350

void gen_one(rvec x[3])
{
  const  rvec   XH1[6] = {
    { AA,     BB,     CC },
    { AA,     BB,     CC },
    { AA,     BB,     CC },
    { -AA,    BB,     CC },
    { -AA,    BB,     CC },
    { BB,     AA,    -CC }
  };
  const  rvec   XH2[6] = {
    { -AA,   BB,   CC },
    { BB,    AA,  -CC },
    { BB,   -AA,  -CC },
    { BB,    AA,  -CC },
    { BB,   -AA,  -CC },
    { BB,   -AA,  -CC }
  };
  static int l=0;
  int        m;
  
  /* This was stolen from Gromos */
  
  for(m=0; (m<DIM); m++) {
    x[1][m]=x[0][m]+XH1[l][m];
    x[2][m]=x[0][m]+XH2[l][m];
  }
      
  l=(l+1) % 6;
}

int main(int argc,char *argv[])
{
  static char *desc[] = {
    "gen_xh generates hydrogen atoms on crystal waters. Normally",
    "in pdb files no hydrogen atoms are included. This program can",
    "add these hydrogens. It is assumed that *ONLY* water oxygen",
    "atoms are in the input file."
  };
  t_manual man = {asize(desc),desc,0,NULL,NULL,0,NULL};

  t_filenm fnm[] = {
    { efGRO, "-f", "confin",  ffREAD },
    { efGRO, "-o", "confout", ffWRITE }
  };
#define NFILE asize(fnm)

  FILE    *out;
  rvec    *xin,xout[3],*v;
  matrix  box;
  int     natoms,i,m,n;
  char    title[256];
  char    *nms[3] = { "OW1", "HW2", "HW3" };
  
  CopyRight(stdout,argv[0]);
  parse_common_args(&argc,argv,0,NFILE,fnm,FALSE,&man);
  
  get_coordnum(opt2fn("-f",NFILE,fnm),&natoms);
  snew(xin,natoms);
  snew(v,natoms);
  
  read_conf(opt2fn("-f",NFILE,fnm),title,&natoms,xin,v,box);
  
  out=ffopen(opt2fn("-o",NFILE,fnm),"w");
  fprintf(out,"%s\n",title);
  fprintf(out,"%5d\n",3*natoms);  
  
  for(i=0; (i<natoms); i++) {
    fprintf(stderr,"\rAtom %d",i);
    copy_rvec(xin[i],xout[0]);
    gen_one(xout);
    for(m=0; (m<DIM); m++) {
      fprintf(out,"%5d%5s%5s%5d",i+1,"SOL",nms[m],3*i+m);
      for(n=0; (n<DIM); n++)
	fprintf(out,"%8.3f",xout[m][n]);
      for(n=0; (n<DIM); n++)
	fprintf(out,"%8.4f",v[i][n]);
      fprintf(out,"\n");
    }
  }
  fprintf(stderr,"\n");
  for(n=0; (n<DIM); n++)
    fprintf(out,"%10.5f",box[n][n]);
  fprintf(out,"\n");
  fclose(out);
  
  thanx(stdout);
  
  return 0;
}
