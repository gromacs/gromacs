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
 * S  C  A  M  O  R  G
 */
static char *SRCID_calcvir_c = "$Id$";

#include "sysstuff.h"
#include "force.h"
#include "assert.h"
#include "vec.h"
#include "mshift.h"
#include "macros.h"
	
static void dprod1(rvec vir,real x,rvec f)
{
  vir[XX]+=x*f[XX];
  vir[YY]+=x*f[YY];
  vir[ZZ]+=x*f[ZZ];
}

static void upd_vir(rvec vir,real dvx,real dvy,real dvz)
{
  vir[XX]-=0.5*dvx;
  vir[YY]-=0.5*dvy;
  vir[ZZ]-=0.5*dvz;
}

void calc_vir(FILE *log,int nxf,rvec x[],rvec f[],tensor vir,
	      t_commrec *cr)
{
  int    i;
  real   dvxx=0,dvxy=0,dvxz=0,dvyx=0,dvyy=0,dvyz=0,dvzx=0,dvzy=0,dvzz=0;
    
  for(i=0; (i<nxf); i++) {
    dvxx+=x[i][XX]*f[i][XX];
    dvxy+=x[i][XX]*f[i][YY];
    dvxz+=x[i][XX]*f[i][ZZ];
    
    dvyx+=x[i][YY]*f[i][XX];
    dvyy+=x[i][YY]*f[i][YY];
    dvyz+=x[i][YY]*f[i][ZZ];
    
    dvzx+=x[i][ZZ]*f[i][XX];
    dvzy+=x[i][ZZ]*f[i][YY];
    dvzz+=x[i][ZZ]*f[i][ZZ];
  }
  
  upd_vir(vir[XX],dvxx,dvxy,dvxz);
  upd_vir(vir[YY],dvyx,dvyy,dvyz);
  upd_vir(vir[ZZ],dvzx,dvzy,dvzz);
}


static void lo_fcv(int i0,int i1,int g0,
		   real x[],real f[],tensor vir,
		   t_ishift is[],real shift_vec[])
{
  int      i,i3,gg,t,t3;
  real     xx,yy,zz;
  real     dvxx=0,dvxy=0,dvxz=0,dvyx=0,dvyy=0,dvyz=0,dvzx=0,dvzy=0,dvzz=0;

  for(i=i0,gg=g0; (i<i1); i++,gg++) {
    i3=DIM*i;
    t=is[gg];
    t3=DIM*t;
    
    xx=x[i3+XX]-shift_vec[t3+XX];
    dvxx+=xx*f[i3+XX];
    dvxy+=xx*f[i3+YY];
    dvxz+=xx*f[i3+ZZ];
    
    yy=x[i3+YY]-shift_vec[t3+YY];
    dvyx+=yy*f[i3+XX];
    dvyy+=yy*f[i3+YY];
    dvyz+=yy*f[i3+ZZ];
    
    zz=x[i3+ZZ]-shift_vec[t3+ZZ];
    dvzx+=zz*f[i3+XX];
    dvzy+=zz*f[i3+YY];
    dvzz+=zz*f[i3+ZZ];
  }
  
  upd_vir(vir[XX],dvxx,dvxy,dvxz);
  upd_vir(vir[YY],dvyx,dvyy,dvyz);
  upd_vir(vir[ZZ],dvzx,dvzy,dvzz);
}

static void lo_fcv2(int i0,int i1,
		    rvec x[],rvec f[],tensor vir,
		    t_ishift is[],rvec shift_vec[])
{
  int      i,gg,t;
  real     xx,yy,zz;
  real     dvxx=0,dvxy=0,dvxz=0,dvyx=0,dvyy=0,dvyz=0,dvzx=0,dvzy=0,dvzz=0;

  for(i=i0,gg=0; (i<i1); i++,gg++) {
    t=is[gg];
    
    xx=x[i][XX]-shift_vec[t][XX];
    dvxx+=xx*f[i][XX];
    dvxy+=xx*f[i][YY];
    dvxz+=xx*f[i][ZZ];
    
    yy=x[i][YY]-shift_vec[t][YY];
    dvyx+=yy*f[i][XX];
    dvyy+=yy*f[i][YY];
    dvyz+=yy*f[i][ZZ];
    
    zz=x[i][ZZ]-shift_vec[t][ZZ];
    dvzx+=zz*f[i][XX];
    dvzy+=zz*f[i][YY];
    dvzz+=zz*f[i][ZZ];
  }
  
  upd_vir(vir[XX],dvxx,dvxy,dvxz);
  upd_vir(vir[YY],dvyx,dvyy,dvyz);
  upd_vir(vir[ZZ],dvzx,dvzy,dvzz);
}

void f_calc_vir(FILE *log,int i0,int i1,rvec x[],rvec f[],tensor vir,
		t_commrec *cr,t_graph *g,rvec shift_vec[])
{
  int start,end;
  
  if (g->nnodes > 0) {
    /* Calculate virial for bonded forces only when they belong to
     * this processor.
     */
    start = max(i0,g->start);
    end   = min(i1,g->end+1);
#ifdef SAFE
    lo_fcv2(start,end,x,f,vir,g->ishift,shift_vec);
#else
    lo_fcv(start,end,0,x[0],f[0],vir,g->ishift,shift_vec[0]);
#endif
    
    /* If not all atoms are bonded, calculate their virial contribution 
     * anyway, without shifting back their coordinates
     */
    if (start > i0) 
      calc_vir(log,start-i0,&(x[i0]),&(f[i0]),vir,cr);
    if (end < i1)
      calc_vir(log,i1-end,&(x[end]),&(f[end]),vir,cr);
  }
  else
    calc_vir(log,i1-i0,&(x[i0]),&(f[i0]),vir,cr);
}
