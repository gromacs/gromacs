/*
 * $Id$
 * 
 *                This source code is part of
 * 
 *                 G   R   O   M   A   C   S
 * 
 *          GROningen MAchine for Chemical Simulations
 * 
 *                        VERSION 3.0
 * 
 * Copyright (c) 1991-2001
 * BIOSON Research Institute, Dept. of Biophysical Chemistry
 * University of Groningen, The Netherlands
 * 
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 * 
 * If you want to redistribute modifications, please consider that
 * scientific software is very special. Version control is crucial -
 * bugs must be traceable. We will be happy to consider code for
 * inclusion in the official distribution, but derived work must not
 * be called official GROMACS. Details are found in the README & COPYING
 * files - if they are missing, get the official version at www.gromacs.org.
 * 
 * To help us fund GROMACS development, we humbly ask that you cite
 * the papers on the package - you can find them in the top README file.
 * 
 * Do check out http://www.gromacs.org , or mail us at gromacs@gromacs.org .
 * 
 * And Hey:
 * GROup of MAchos and Cynical Suckers
 */
static char *SRCID_g_bundle_c = "$Id$";
#include <math.h>
#include <string.h>
#include "statutil.h"
#include "sysstuff.h"
#include "typedefs.h"
#include "smalloc.h"
#include "macros.h"
#include "vec.h"
#include "copyrite.h"
#include "futil.h"
#include "statutil.h"
#include "rdgroup.h"
#include "xvgr.h"
#include "rmpbc.h"
#include "tpxio.h"
#include "physics.h"

typedef struct{
  int n;
  rvec *end[2];
  rvec *mid;
  rvec *dir;
  real *len;
} t_bundle;

static void rotate_ends(t_bundle *bun,rvec axis,int c0,int c1)
{
  int  end,i;
  rvec ax,tmp;

  unitv(axis,ax);
  for(end=0; end<2; end++)
    for(i=0; i<bun->n; i++) {
      copy_rvec(bun->end[end][i],tmp);
      bun->end[end][i][c0] = ax[c1]*tmp[c0] - ax[c0]*tmp[c1];
      bun->end[end][i][c1] = ax[c0]*tmp[c0] + ax[c1]*tmp[c1];
    }
  copy_rvec(axis,tmp);
  axis[c0] = ax[c1]*tmp[c0] - ax[c0]*tmp[c1];
  axis[c1] = ax[c0]*tmp[c0] + ax[c1]*tmp[c1];
}

static void calc_axes(rvec x[],t_atom atom[],int gnx[],atom_id *index[],
		      bool bRot,t_bundle *bun)
{
  int  end,i,div,d;
  real *mtot,m;
  rvec axis[2],cent;
  
  snew(mtot,bun->n);

  for(end=0; end<2; end++) {
    for(i=0; i<bun->n; i++) {
      clear_rvec(bun->end[end][i]);
      mtot[i] = 0;
    }
    div = gnx[end]/bun->n;
    for(i=0; i<gnx[end]; i++) {
      m = atom[index[end][i]].m;
      for(d=0; d<DIM; d++)
	bun->end[end][i/div][d] += m*x[index[end][i]][d];
      mtot[i/div] += m;
    }
    clear_rvec(axis[end]);
    for(i=0; i<bun->n; i++) {
      svmul(1.0/mtot[i],bun->end[end][i],bun->end[end][i]);
      rvec_inc(axis[end],bun->end[end][i]);
    }
    svmul(1.0/bun->n,axis[end],axis[end]);
  }
  sfree(mtot);

  rvec_add(axis[0],axis[1],cent);
  svmul(0.5,cent,cent);
  /* center the bundle on the origin */
  for(end=0; end<2; end++) {
    rvec_dec(axis[end],cent);
    for(i=0; i<bun->n; i++)
      rvec_dec(bun->end[end][i],cent);
  }
  if (bRot) {
    /* rotate the axis parallel to the z-axis */
    rotate_ends(bun,axis[0],YY,ZZ);
    rotate_ends(bun,axis[0],XX,ZZ);
  }
  for(i=0; i<bun->n; i++) {
    rvec_add(bun->end[0][i],bun->end[1][i],bun->mid[i]);
    svmul(0.5,bun->mid[i],bun->mid[i]);
    rvec_sub(bun->end[0][i],bun->end[1][i],bun->dir[i]);
    bun->len[i] = norm(bun->dir[i]);
    unitv(bun->dir[i],bun->dir[i]);
  }
}

static void dump_axes(int fp,t_trxframe *fr,t_atoms *outat,t_bundle *bun)
{
  t_trxframe  frout;
  static rvec *xout=NULL;
  int         i,end;
  
  if (xout==NULL)
    snew(xout,outat->nr);

  for(i=0; i<bun->n; i++) {
    copy_rvec(bun->end[0][i],xout[3*i]);
    copy_rvec(bun->mid[i],   xout[3*i+1]);
    copy_rvec(bun->end[1][i],xout[3*i+2]);
  }
  frout = *fr;
  frout.bV     = FALSE;
  frout.bF     = FALSE;
  frout.bBox   = FALSE;
  frout.bAtoms = TRUE;
  frout.natoms = outat->nr;
  frout.atoms  = outat;
  frout.x      = xout;
  write_trxframe(fp,&frout);
}

int main(int argc,char *argv[])
{
  static char *desc[] = {
    "g_bundle analyzes bundles of axes. The axes can be for instance",
    "helix axes. The program reads two index groups and divides both",
    "of them in [TT]-na[tt] parts. The groups define the tops and",
    "bottoms of the axes. Several quantities are written to file:",
    "the axis length, the distance from the axis centers to the average",
    "center of all axes, the total tilt, the radial tilt and the lateral",
    "tilt with respect to the average axis.",
    "[PAR]",
    "With option [TT]-oa[tt] the top, mid and bottom points of each axis",
    "are written to a pdb file each frame. The residue numbers correspond",
    "to the axis numbers. When viewing this file with [TT]rasmol[tt] use the",
    "command line option [TT]-nmrpdb[tt], and type [TT]set axis true[tt] to",
    "display the reference axis."
  };
  static int  n=0;
  static bool bZ=FALSE;
  t_pargs pa[] = {
    { "-na", FALSE, etINT, {&n},
	"Number of axes" },
    { "-z", FALSE, etBOOL, {&bZ},
	"Use the Z-axis as reference iso the average axis" }
  };
  FILE       *out,*flen,*fdist,*ftilt,*ftiltr,*ftiltl;
  int        status,fpdb;
  t_topology top;
  rvec       *xtop;
  matrix     box;
  t_trxframe fr;
  t_atoms    outatoms;
  real       t,comp;
  int        natoms;
  char       *grpname[2],title[256],*anm="CA",*rnm="GLY";
  int        i,j,gnx[2];
  atom_id    *index[2];
  t_bundle   bun;
#define NLEG asize(leg) 
  t_filenm fnm[] = { 
    { efTRX, "-f", NULL, ffREAD }, 
    { efTPS, NULL, NULL, ffREAD }, 
    { efNDX, NULL, NULL, ffOPTRD },
    { efXVG, "-ol", "bun_len", ffWRITE },
    { efXVG, "-od", "bun_dist", ffWRITE },
    { efXVG, "-ot", "bun_tilt", ffWRITE },
    { efXVG, "-otr", "bun_tiltr", ffWRITE },
    { efXVG, "-otl", "bun_tiltl", ffWRITE },
    { efPDB, "-oa", "axes", ffOPTWR }
  }; 
#define NFILE asize(fnm) 

  CopyRight(stderr,argv[0]); 
  parse_common_args(&argc,argv,PCA_CAN_TIME | PCA_BE_NICE,
		    NFILE,fnm,asize(pa),pa,asize(desc),desc,0,NULL); 

  read_tps_conf(ftp2fn(efTPS,NFILE,fnm),title,&top,&xtop,NULL,box,TRUE);

  fprintf(stderr,"Select a group of top and a group of bottom atoms\n");
  get_index(&top.atoms,ftp2fn_null(efNDX,NFILE,fnm),2,gnx,index,grpname);

  if (n<=0 || gnx[0] % n || gnx[1] % n)
    fatal_error(0,
		"The size of one of your index groups is not a multiple of n");
  bun.n = n;
  snew(bun.end[0],n);
  snew(bun.end[1],n);
  snew(bun.mid,n);
  snew(bun.dir,n);
  snew(bun.len,n);

  flen  = xvgropen(opt2fn("-ol",NFILE,fnm),"Axis lengths",
		   "Time (ps)","(nm)");
  fdist = xvgropen(opt2fn("-od",NFILE,fnm),"Distance of axis centers",
		   "Time (ps)","(nm)");
  ftilt = xvgropen(opt2fn("-ot",NFILE,fnm),"Axis tilts",
		   "Time (ps)","(degrees)");
  ftiltr = xvgropen(opt2fn("-otr",NFILE,fnm),"Radial axis tilts",
		   "Time (ps)","(degrees)");
  ftiltl = xvgropen(opt2fn("-otl",NFILE,fnm),"Lateral axis tilts",
		   "Time (ps)","(degrees)");

  if (opt2bSet("-oa",NFILE,fnm)) {
    init_t_atoms(&outatoms,3*n,FALSE);
    outatoms.nr = 3*n;
    for(i=0; i<3*n; i++) {
      outatoms.atomname[i] = &anm;
      outatoms.atom[i].resnr = i/3;
      outatoms.resname[i/3] = &rnm;
    }
    fpdb = open_trx(opt2fn("-oa",NFILE,fnm),"w");
  } else
    fpdb = -1;
  
  read_first_frame(&status,ftp2fn(efTRX,NFILE,fnm),&fr,TRX_NEED_X); 
  
  do {
    rm_pbc(&top.idef,fr.natoms,fr.box,fr.x,fr.x);
    calc_axes(fr.x,top.atoms.atom,gnx,index,!bZ,&bun);
    fprintf(flen," %10g",fr.time);
    fprintf(fdist," %10g",fr.time);
    fprintf(ftilt," %10g",fr.time);
    fprintf(ftiltr," %10g",fr.time);
    fprintf(ftiltl," %10g",fr.time);

    for(i=0; i<bun.n; i++) {
      fprintf(flen," %6g",bun.len[i]);
      fprintf(fdist," %6g",norm(bun.mid[i]));
      fprintf(ftilt," %6g",RAD2DEG*acos(bun.dir[i][ZZ]));
      comp = bun.mid[i][XX]*bun.dir[i][XX]+bun.mid[i][YY]*bun.dir[i][YY];
      fprintf(ftiltr," %6g",RAD2DEG*
	      asin(comp/sqrt(sqr(comp)+sqr(bun.dir[i][ZZ]))));
      comp = bun.mid[i][YY]*bun.dir[i][XX]-bun.mid[i][XX]*bun.dir[i][YY];
      fprintf(ftiltl," %6g",RAD2DEG*
	      asin(comp/sqrt(sqr(comp)+sqr(bun.dir[i][ZZ]))));
    }
    fprintf(flen,"\n");
    fprintf(fdist,"\n");
    fprintf(ftilt,"\n");
    fprintf(ftiltr,"\n");
    fprintf(ftiltl,"\n");
    if (fpdb >= 0)
      dump_axes(fpdb,&fr,&outatoms,&bun);
  } while(read_next_frame(status,&fr));

  close_trx(status);
  
  if (fpdb >= 0)
    close_trx(fpdb);
  fclose(fdist);
  fclose(flen);
  fclose(ftilt);
  fclose(ftiltr);
  fclose(ftiltl);
  
  do_view(ftp2fn(efXVG,NFILE,fnm),NULL);
  
  thanx(stderr);
  
  return 0;
}
