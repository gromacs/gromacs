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
static char *SRCID_g_sas_c = "$Id$";

#include <math.h>
#include <stdlib.h>
#include "sysstuff.h"
#include "string.h"
#include "typedefs.h"
#include "smalloc.h"
#include "macros.h"
#include "vec.h"
#include "xvgr.h"
#include "pbc.h"
#include "copyrite.h"
#include "futil.h"
#include "statutil.h"
#include "rdgroup.h"
#include "nsc.h"
#include "vdw.h"
#include "pdbio.h"
#include "rmpbc.h"

void connelly_plot(char *fn,int ndots,real dots[],rvec x[],t_atoms *atoms,
		   matrix box)
{
  static char *atomnm="DOT";
  static char *resnm ="DOT";
  FILE *fp;
  t_pdbatom *pdba;
  int  i,i0,k;
  
  pdba = atoms2pdba(atoms,x);
  i0 = atoms->nr;
  srenew(pdba,atoms->nr+ndots);
  for(i=k=0; (i<ndots); i++) {
    strcpy(pdba[i0+i].resnm,resnm);
    strcpy(pdba[i0+i].atomnm,atomnm);
    sprintf(pdba[i0+i].pdbresnr,"1");
    pdba[i0+i].pdbtp = epdbATOM;
    pdba[i0+i].chain = ' ';
    pdba[i0+i].atomnr= i0+i;
    pdba[i0+i].resnr = 1;
    pdba[i0+i].x[XX] = dots[k++];
    pdba[i0+i].x[YY] = dots[k++];
    pdba[i0+i].x[ZZ] = dots[k++];
    pdba[i0+i].bfac  = 0.0;
    pdba[i0+i].dummy = 0.0;
  }
  fp=ffopen(fn,"w");
  print_pdbatoms(fp,i0+ndots,pdba,box);
  fclose(fp);
  
  sfree(pdba);
}

real calc_radius(char *atom)
{
  real r;
  
  switch (atom[0]) {
  case 'C':
    r = 0.16;
    break;
  case 'O':
    r = 0.13;
    break;
  case 'N':
    r = 0.14;
    break;
  case 'S':
    r = 0.2;
    break;
  case 'H':
    r = 0.1;
    break;
  default:
    r = 1e-3;
  }
  return r;
}

void sas_plot(int nfile,t_filenm fnm[],real solsize,real defvdw,int ndots,
	      real qcut)
{
  FILE         *fp;
  real         t;
  int          nvdw,status;
  int          i,j,k,natoms,flag,nsurfacedots;
  rvec         *x;
  matrix       box;
  t_vdw        *vdw;
  t_topology   *top;
  bool         *bPhobic;
  bool         bConnelly;
  real         *radius,*area,*surfacedots;
  real         totarea,totvolume,harea;
    
  if ((natoms=read_first_x(&status,ftp2fn(efTRX,nfile,fnm),
			   &t,&x,box))==0) {
    fprintf(stderr,"Could not read coordinates from statusfile\n");
    exit(1);
  }
  top=read_top(ftp2fn(efTPX,nfile,fnm));
  
  nvdw = read_vdw(ftp2fn(efVDW,nfile,fnm),&vdw);
  fprintf(stderr,"There are %d VanderWaals radii\n",nvdw);
  
  /* Now comput atomic readii including solvent probe size */
  snew(radius,natoms);
  snew(bPhobic,natoms);
  for(i=0; (i<natoms); i++) {
    radius[i]  = calc_radius(*(top->atoms.atomname[i])) + solsize;
    /*get_vdw(nvdw,vdw,*(top->atoms.atomname[i])) + solsize;*/
    bPhobic[i] = fabs(top->atoms.atom[i].q) <= qcut;
    /*(*(top->atoms.atomname[i]))[0] == 'C';*/
  }
  fp=xvgropen(ftp2fn(efXVG,nfile,fnm),"Solvent Accessible Surface","Time (ps)",
	      "Area (nm\\S2\\N)");
  
  j=0;
  do {
    if ((j++ % 10) == 0)
      fprintf(stderr,"\rframe: %5d",j-1);
      
    rm_pbc(&top->idef,natoms,box,x,x);

    bConnelly = ((j == 1) && (opt2bSet("-q",nfile,fnm)));
    if (bConnelly)
      flag = FLAG_ATOM_AREA | FLAG_DOTS;
    else
      flag = FLAG_ATOM_AREA;
    
    if (NSC(x[0],radius,natoms,ndots,flag,&totarea,
	    &area,&totvolume,&surfacedots,&nsurfacedots))
      fatal_error(0,"Something wrong in NSC");
      
    if (bConnelly)
      connelly_plot(ftp2fn(efPDB,nfile,fnm),
		    nsurfacedots,surfacedots,x,&(top->atoms),box);
      
    harea = 0;
    for(i=0; (i<natoms); i++) {
      if (bPhobic[i])
	harea += area[i];
    }
    fprintf(fp,"%10g  %10g  %10g\n",t,harea,totarea);
    
    sfree(area);
    sfree(surfacedots);
      
  } while (read_next_x(status,&t,natoms,x,box));
  
  fprintf(stderr,"\n");
  close_trj(status);

  sfree(x);
}

int main(int argc,char *argv[])
{
  static char *desc[] = {
    "g_sas computes hydrophobic and total solvent accessible surface area."
  };

  static real solsize = 0.14;
  static real defvdw  = 0.15;
  static int  ndots   = 24;
  static real qcut    = 0.2;
  t_pargs pa[] = {
    { "-solsize", FALSE, etREAL, &solsize,
	"Radius of the solvent probe (nm)" },
    { "-defvdw",  FALSE, etREAL, &defvdw,
	"Default Van der Waals radius for unknown atoms" },
    { "-ndots",   FALSE, etINT,  &ndots,
	"Number of dots per sphere, more dots means more accuracy" },
    { "-qmax",    FALSE, etREAL, &qcut,
	"The maximum charge (e, absolute value) of a hydrophobic atom" }
  };
  t_filenm  fnm[] = {
    { efTRX, "-f",   NULL,       ffREAD },
    { efTPX, "-s",   NULL,       ffREAD },
    { efVDW, "-vdw", NULL,       ffREAD },
    { efXVG, "-o",   "area",     ffWRITE },
    { efPDB, "-q",   "connelly", ffOPTWR }
  };
#define NFILE asize(fnm)

  CopyRight(stderr,argv[0]);
  parse_common_args(&argc,argv,PCA_CAN_VIEW | PCA_CAN_TIME,TRUE,
		    NFILE,fnm,asize(pa),pa,asize(desc),desc,0,NULL);
  if (solsize <= 0) {
    solsize=1e-3;
    fprintf(stderr,"Solsize too small, setting it to %g\n",solsize);
  }
  if (ndots < 20) {
    ndots = 20;
    fprintf(stderr,"Ndots too small, setting it to %d\n",ndots);
  }
  
  sas_plot(NFILE,fnm,solsize,defvdw,ndots,qcut);
  
  xvgr_file(opt2fn("-o",NFILE,fnm),"-nxy");
  
  thanx(stdout);
  
  return 0;
}

