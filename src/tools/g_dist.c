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
static char *SRCID_g_dist_c = "$Id$";

#include <typedefs.h>
#include "smalloc.h"
#include "macros.h"
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
#include "pbc.h"

int main(int argc,char *argv[])
{
  static char *desc[] = {
    "g_dist can calculate the distance between the centers of mass of two",
    "groups of atoms as a function of time. The total distance and its",
    "x, y and z components are plotted.[PAR]",
    "Or when [TT]-dist[tt] is set, print all the atoms in group 2 that are",
    "closer than a certain distance to the center of mass of group 1."
  };
  
  t_topology *top=NULL;
  real t,cut2,dist2;
  rvec *x=NULL,*v=NULL,dx;
  matrix box;
  int status;
  int natoms;

  int g,d,i,j,res,teller=0;
  atom_id aid;

  int     ngrps;     /* the number of index groups */
  atom_id **index,max;   /* the index for the atom numbers */
  int     *isize;    /* the size of each group */
  char    **grpname; /* the name of each group */
  rvec    *com;
  real    *mass;
  FILE    *fp=NULL;
  bool    bCutoff;

  char    *leg[4] = { "|d|","d\\sx\\N","d\\sy\\N","d\\sz\\N" };

  static real cut=0;

  static t_pargs pa[] = {
    { "-dist",      FALSE, etREAL, {&cut},
      "Print all atoms in group 2 closer than dist to the center of mass of group 1" },
  };
#define NPA asize(pa)

  t_filenm fnm[] = {
    { efTRX, "-f", NULL, ffREAD },
    { efTPX, NULL, NULL, ffREAD },
    { efNDX, NULL, NULL, ffOPTRD },
    { efXVG, NULL, "dist", ffOPTWR },
  };
#define NFILE asize(fnm)


  CopyRight(stdout,argv[0]);

  parse_common_args(&argc,argv,PCA_CAN_TIME,TRUE,
		    NFILE,fnm,NPA,pa,asize(desc),desc,0,NULL);
  
  bCutoff=opt2parg_bSet("-dist",NPA,pa);
  cut2=cut*cut;
  
  top=read_top(ftp2fn(efTPX,NFILE,fnm));
  
  /* read index files */
  ngrps = 2;
  snew(com,ngrps);
  snew(grpname,ngrps);
  snew(index,ngrps);
  snew(isize,ngrps);
  get_index(&top->atoms,ftp2fn(efNDX,NFILE,fnm),ngrps,isize,index,grpname);
  
  /* calculate mass */
  max=0;
  snew(mass,ngrps);
  for(g=0;(g<ngrps);g++) {
    mass[g]=0;
    for(i=0;(i<isize[g]);i++) {
      if (index[g][i]>max)
	max=index[g][i];
      if (index[g][i] >= top->atoms.nr)
	fatal_error(0,"Atom number %d, item %d of group %d, is larger than number of atoms in the topolgy (%d)\n",index[g][i]+1,i+1,g+1,top->atoms.nr+1);
      mass[g]+=top->atoms.atom[index[g][i]].m;
    }
  }

  natoms=read_first_x(&status,ftp2fn(efTRX,NFILE,fnm),&t,&x,box);

  if (max>=natoms)
    fatal_error(0,"Atom number %d in an index group is larger than number of atoms in the trajectory (%d)\n",(int)max+1,natoms);

  if (!bCutoff) {
    /* open output file */
    fp = xvgropen(ftp2fn(efXVG,NFILE,fnm),
		  "Distance","Time (ps)","Distance (nm)");
    xvgr_legend(fp,4,leg);
  } else
    ngrps=1;
  
  do {
    /* initialisation for correct distance calculations */
    init_pbc(box,FALSE);
    /* make molecules whole again */
    rm_pbc(&top->idef,natoms,box,x,x);

    /* calculate center of masses */
    for(g=0;(g<ngrps);g++) {
      for(d=0;(d<DIM);d++) {
	com[g][d]=0;
	for(i=0;(i<isize[g]);i++) {
	  com[g][d] += x[index[g][i]][d] * top->atoms.atom[index[g][i]].m;
	}
	com[g][d] /= mass[g];
      }
    }
    
    if (!bCutoff) {
      /* write to output */
      fprintf(fp,"%8.3f ",t);
      for(g=0;(g<ngrps/2);g++) {
	pbc_dx(com[2*g],com[2*g+1],dx);
	fprintf(fp,"%10.5f %10.5f %10.5f %10.5f",
		norm(dx),dx[XX],dx[YY],dx[ZZ]);
      }
      fprintf(fp,"\n");
    } else {
      for(i=0;(i<isize[1]);i++) { 
	j=index[1][i];
	pbc_dx(x[j],com[0],dx);
	dist2 = norm2(dx);
	if (dist2<cut2) {
	  res=top->atoms.atom[j].resnr;
	  fprintf(stdout,"\rt: %g  %d %s %d %s  %g (nm)\n",
		  t,res+1,*top->atoms.resname[res],
		  j+1,*top->atoms.atomname[j],sqrt(dist2));     
	} 
      }
    }
    
    teller++;
  } while (read_next_x(status,&t,natoms,x,box));

  if (!bCutoff)
    fclose(fp);

  close_trj(status);
  
  thanx(stderr);
  return 0;
}
