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
static char *SRCID_g_gyrate_c = "$Id$";

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
#include "princ.h"
#include "rmpbc.h"
#include "txtdump.h"
#include "tpxio.h"

typedef struct {
  atom_id i;
  real    d;
} t_order;

t_order *order;

static int ocomp(const void *a,const void *b)
{
  t_order *oa,*ob;
  
  oa = (t_order *)a;
  ob = (t_order *)b;
  
  if (oa->d < ob->d)
    return -1;
  else
    return 1;  
}

int main(int argc,char *argv[])
{
  static char *desc[] = {
    "trjorder orders molecules according to the smallest distance",
    "to atoms in a reference group. It will ask for a group of reference",
    "atoms and a group of molecules. For each frame of the trajectory",
    "the selected molecules will be reordered according to the shortest",
    "distance between atom number [TT]-da[tt] in the molecule and all the",
    "atoms in the reference group. All atoms in the trajectory are written",
    "to the output trajectory.[PAR]",
    "trjorder can be useful for e.g. analyzing the n waters closest to a",
    "protein.",
    "In that case the reference group would be the protein and the group",
    "of molecules would consist of all the water atoms. When an index group",
    "of the first n waters is made, the ordered trajectory can be used",
    "with any Gromacs program to analyze the n closest waters."
  };
  static int na=3,ref_a=1;
  t_pargs pa[] = {
    { "-na", FALSE, etINT, {&na},
      "Number of atoms in a molecule" },
    { "-da", FALSE, etINT, {&ref_a},
      "Atom used for the distance calculation" }
  };
  int        status,out;
  t_topology top;
  rvec       *x,dx;
  matrix     box;
  real       t;
  int        natoms,nwat;
  char       **grpname,title[256];
  int        i,j,*isize;
  atom_id    sa,*swi,**index;
#define NLEG asize(leg) 
  t_filenm fnm[] = { 
    { efTRX, "-f", NULL, ffREAD  }, 
    { efTPS, NULL, NULL, ffREAD  }, 
    { efNDX, NULL, NULL, ffOPTRD },
    { efTRX, "-o", "ordered", ffWRITE } 
  }; 
#define NFILE asize(fnm) 

  CopyRight(stderr,argv[0]); 
  parse_common_args(&argc,argv,PCA_CAN_TIME,TRUE,
		    NFILE,fnm,asize(pa),pa,asize(desc),desc,0,NULL); 

  read_tps_conf(ftp2fn(efTPS,NFILE,fnm),title,&top,&x,NULL,box,TRUE);
  sfree(x);

  /* get index groups */
  printf("Select a group of reference atoms and a group of molecules to be ordered:\n"); 
  snew(grpname,2);
  snew(index,2);
  snew(isize,2);
  get_index(&top.atoms,ftp2fn_null(efNDX,NFILE,fnm),2,isize,index,grpname);

  natoms=read_first_x(&status,ftp2fn(efTRX,NFILE,fnm),&t,&x,box); 
  if (natoms > top.atoms.nr)
    fatal_error(0,"Number of atoms in the run input file is larger than in the trjactory");
  for(i=0; i<2; i++)
    for(j=0; j<isize[i]; j++)
      if (index[i][j] > natoms)
	fatal_error(0,"An atom number in group %s is larger than the number of atoms in the trajectory");
  
  if (isize[1] % na)
    fatal_error(0,"Number of atoms in the molecule group (%d) is not a multiple of na (%d)",isize[1],na);
  nwat = isize[1]/na;
  if (ref_a > na)
    fatal_error(0,"The reference atom can not be larger than the number of atoms in a molecule");
  ref_a--;
  snew(order,nwat);
  snew(swi,natoms);
  for(i=0; i<natoms; i++)
    swi[i] = i;
  
  out=open_trx(opt2fn("-o",NFILE,fnm),"w");
  do {
    rm_pbc(&top.idef,natoms,box,x,x);
    init_pbc(box,FALSE);
    
    for(i=0; i<nwat; i++) {
      sa = index[1][na*i];
      pbc_dx(x[index[0][0]],x[sa+ref_a],dx);
      order[i].i = sa;
      order[i].d = norm2(dx); 
    }
    for(j=1; j<isize[0]; j++)
      for(i=0; i<nwat; i++) {
	sa = index[1][na*i];
	pbc_dx(x[index[0][j]],x[sa+ref_a],dx);
	if (norm2(dx) < order[i].d)
	  order[i].d = norm2(dx);
      }

    qsort(order,nwat,sizeof(*order),ocomp);
    for(i=0; i<nwat; i++)
      for(j=0; j<na; j++)
	swi[index[1][na*i]+j] = order[i].i+j;

    write_trx(out,natoms,swi,&top.atoms,0,t,box,x,NULL);
    
  } while(read_next_x(status,&t,natoms,x,box));
  close_trj(status);
  close_trx(out);

  thanx(stderr);
  
  return 0;
}
