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
static char *SRCID_ipgro_c = "$Id$";

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

int main(int argc,char *argv[])
{
  static char *desc[] = {
    "ipgro interpolates two gromos coordinate files and writes a trajectory."
  };
  
  FILE         *out;
  char         title1[256],title2[256];
  rvec         *x1,*x2,*v1,*v2,*xm;
  matrix       box1,box2;
  int          i,j,np,n1,n2,nat,m;
  real         frac,frac1;
  t_filenm fnm[] = {
    { efTRJ, "-o", NULL,  ffWRITE },
    { efGRO, "-c1", NULL, ffREAD },
    { efGRO, "-c2", NULL, ffREAD },
  };
#define NFILE asize(fnm)

  CopyRight(stderr,argv[0]);
  parse_common_args(&argc,argv,0,TRUE,
		    NFILE,fnm,0,NULL,asize(desc),desc,0,NULL);
  if (argc < 2)
    fatal_error(0,"Usage: %s npoints",argv[0]);
  np=atoi(argv[1]);
  if (np < 2)
    fatal_error(0,"np < 2");
    
  get_coordnum(opt2fn("-c1",NFILE,fnm),&n1);
  get_coordnum(opt2fn("-c2",NFILE,fnm),&n2);
  if (n1 != n2)
    fatal_error(0,"n1 != n2");
  nat=n1;
  
  snew(x1,nat);
  snew(v1,nat);
  snew(x2,nat);
  snew(v2,nat);
  snew(xm,nat);
  read_conf(opt2fn("-c1",NFILE,fnm),title1,&n1,x1,v1,box1);
  fprintf(stderr,"%s\n",title1);
  read_conf(opt2fn("-c2",NFILE,fnm),title2,&n2,x2,v2,box2);
  fprintf(stderr,"%s\n",title2);
  
  out=ftp2FILE(efTRJ,NFILE,fnm,"w");
  frac=1.0/(np - 1);
  frac1=1.0-frac;
  for(i=0; (i<np); i++) {
    for(j=0; (j<nat); j++)
      for(m=0; (m<DIM); m++) {
	xm[j][m]=frac1*x1[j][m]+frac*x2[j][m];
      }
    wr_status(out,i,i,0,NULL,box1,NULL,NULL,nat,xm,NULL,NULL,0,NULL,NULL);
  }
  fclose(out);
  
  thanx(stdout);
  
  return 0;
}

