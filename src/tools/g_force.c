/*
 * $Id$
 * 
 *       This source code is part of
 * 
 *        G   R   O   M   A   C   S
 * 
 * GROningen MAchine for Chemical Simulations
 * 
 *               VERSION 2.0
 * 
 * Copyright (c) 1991-1999
 * BIOSON Research Institute, Dept. of Biophysical Chemistry
 * University of Groningen, The Netherlands
 * 
 * Please refer to:
 * GROMACS: A message-passing parallel molecular dynamics implementation
 * H.J.C. Berendsen, D. van der Spoel and R. van Drunen
 * Comp. Phys. Comm. 91, 43-56 (1995)
 * 
 * Also check out our WWW page:
 * http://md.chem.rug.nl/~gmx
 * or e-mail to:
 * gromacs@chem.rug.nl
 * 
 * And Hey:
 * GRowing Old MAkes el Chrono Sweat
 */
static char *SRCID_g_com_c = "$Id$";

#include "sysstuff.h"
#include "smalloc.h"
#include "macros.h"
#include "statutil.h"
#include "random.h"
#include "names.h"
#include "vec.h"
#include "futil.h"
#include "copyrite.h"
#include "xvgr.h"
#include "string2.h"
#include "rdgroup.h"
#include "tpxio.h"

void calc_ftot(rvec f[],rvec ftot,int isize,atom_id *index)
{
  int  i,m;
  atom_id aid;
  
  clear_rvec(ftot);
  for(i=0; (i<isize); i++) {
    aid=index[i];
    for(m=0; (m<DIM); m++) {
      ftot[m]+=f[aid][m];
    }
  }
}
	
int main(int argc,char *argv[])
{
  static char *desc[] = {
    "g_force computes the net force on a group of atoms by summing up forces",
    "from the trajectory file over a user specified index."
  };
  t_filenm fnm[] = {
    { efTRX,  "-f",  NULL, ffREAD },
    { efNDX,  NULL,  NULL, ffREAD },
    { efXVG, "-o",  "f",   ffWRITE },
  };
#define NFILE asize(fnm)

  static char  *axisX[]={ "Fx", "Fy", "Fz", "Ftot" };
 
  /* index stuff */
  int      ngrps;       /* the number of groups */
  int      *isize;      /* the size of each group */
  char     **grpnames;  /* the name of each group */
  atom_id  **index;     /* the index array of each group */
  t_trxframe fr;
  int      flags;
  int      g;           /* group counter */
  char     format[STRLEN],filename[STRLEN],title[STRLEN];
  FILE     **outX,*outek=NULL;
  int      status,ftpout;
  int      i,j,idum,step,natoms;
  rvec     ftot;
  int      d,m,n;
  matrix   box;

  CopyRight(stderr,argv[0]);
  parse_common_args(&argc,argv,PCA_CAN_TIME,TRUE,NFILE,fnm,
		    0,NULL,asize(desc),desc,0,NULL);

  flags = TRX_NEED_F | TRX_DONT_SKIP;
  read_first_frame(&status,ftp2fn(efTRX,NFILE,fnm),&fr,flags);
  natoms = fr.natoms;
  
  fprintf(stderr,"How many groups do you want to calc com of ? ");
  scanf("%d",&ngrps);
  fprintf(stderr,"OK. I will calc net force on %d groups\n",ngrps);
  
  snew(grpnames,ngrps);
  snew(index,ngrps);
  snew(isize,ngrps);
  
  rd_index(ftp2fn_null(efNDX,NFILE,fnm),ngrps,isize,index,grpnames);
  
  /* open output files */
  snew(outX,ngrps);
  if (ngrps==1) {
    outX[0]=xvgropen(opt2fn("-o",NFILE,fnm),"Force","Time(ps)","F (kJ/mol nm)");
    xvgr_legend(outX[0],asize(axisX),axisX);
  } 
  else {
    strcpy(format,opt2fn("-o",NFILE,fnm));
    format[strlen(format)-4]='\0';
    strcat(format,"_%s.xvg");
    for(g=0;(g<ngrps);g++) {
      /* coordinates */
      sprintf(filename,format,grpnames[g]);
      outX[g]=xvgropen(filename,"Force","Time(ps)","F (kJ/mol nm)");
      xvgr_legend(outX[g],asize(axisX),axisX);
    }
  }
  
  do {
    for(g=0;(g<ngrps);g++) {
      calc_ftot(fr.f,ftot,isize[g],index[g]);
      fprintf(outX[g],"%10g  %10g  %10g  %10g  %10g\n",
	      fr.time,ftot[XX],ftot[YY],ftot[ZZ],norm(ftot));
    }
  } while (read_next_frame(status,&fr));

  sfree(fr.f);
  
  close_trj(status);
  for(g=0;(g<ngrps);g++) {
    fclose(outX[g]);
  }
  
  thanx(stderr);
  
  return 0;
}

