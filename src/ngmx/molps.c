/*
 * 
 *                This source code is part of
 * 
 *                 G   R   O   M   A   C   S
 * 
 *          GROningen MAchine for Chemical Simulations
 * 
 *                        VERSION 3.2.0
 * Written by David van der Spoel, Erik Lindahl, Berk Hess, and others.
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2004, The GROMACS development team,
 * check out http://www.gromacs.org for more information.

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
 * For more info, check our website at http://www.gromacs.org
 * 
 * And Hey:
 * Gyas ROwers Mature At Cryogenic Speed
 */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <math.h>
#include "sysstuff.h"
#include "string.h"
#include "smalloc.h"
#include "macros.h"
#include "xutil.h"
#include "3dview.h"
#include "gmx_fatal.h"
#include "buttons.h"
#include "manager.h"
#include "nmol.h"
#include "writeps.h"
#include "nleg.h"

#define MSIZE 4

static void ps_draw_atom(t_psdata ps,atom_id ai,iv2 vec2[],char **atomnm[])
{
  int xi,yi;
  
  xi=vec2[ai][XX];
  yi=vec2[ai][YY];
  ps_rgb(ps,Type2RGB(*atomnm[ai]));
  ps_line(ps,xi-MSIZE,yi,xi+MSIZE+1,yi);
  ps_line(ps,xi,yi-MSIZE,xi,yi+MSIZE+1);
}

/* Global variables */
static rvec gl_fbox,gl_hbox,gl_mhbox;

static void init_pbc(matrix box)
{
  int i;

  for(i=0; (i<DIM); i++) {
    gl_fbox[i]  =  box[i][i];
    gl_hbox[i]  =  gl_fbox[i]*0.5;
    gl_mhbox[i] = -gl_hbox[i];
  }
}

static gmx_bool local_pbc_dx(rvec x1, rvec x2)
{
  int  i;
  real dx;
  
  for(i=0; (i<DIM); i++) {
    dx=x1[i]-x2[i];
    if (dx > gl_hbox[i])
      return FALSE;
    else if (dx <= gl_mhbox[i])
      return FALSE;
  }
  return TRUE;
}

static void ps_draw_bond(t_psdata ps,
			 atom_id ai,atom_id aj,iv2 vec2[],
			 rvec x[],char **atomnm[],int size[],gmx_bool bBalls)
{
  char    *ic,*jc;
  int     xi,yi,xj,yj;
  int     xm,ym;

  if (bBalls) {
    ps_draw_atom(ps,ai,vec2,atomnm);
    ps_draw_atom(ps,aj,vec2,atomnm);
  }
  else {
    if (local_pbc_dx(x[ai],x[aj])) {
      ic=*atomnm[ai];
      jc=*atomnm[aj];
      xi=vec2[ai][XX];
      yi=vec2[ai][YY];
      xj=vec2[aj][XX];
      yj=vec2[aj][YY];
      
      if (ic != jc) {
	xm=(xi+xj) >> 1;
	ym=(yi+yj) >> 1;
      
	ps_rgb(ps,Type2RGB(ic));
	ps_line(ps,xi,yi,xm,ym);
	ps_rgb(ps,Type2RGB(jc));
	ps_line(ps,xm,ym,xj,yj);
      }
      else {
	ps_rgb(ps,Type2RGB(ic));
	ps_line(ps,xi,yi,xj,yj);
      }
    }
  }
}

void ps_draw_objects(t_psdata ps,int nobj,t_object objs[],iv2 vec2[],rvec x[],
		     char **atomnm[],int size[],gmx_bool bShowHydro,int bond_type,
		     gmx_bool bPlus)
{
  gmx_bool     bBalls;
  int      i;
  t_object *obj;

  bBalls=FALSE;
  for(i=0; (i<nobj); i++) {
    obj=&(objs[i]);
    switch (obj->eO) {
    case eOSingle:
      ps_draw_atom(ps,obj->ai,vec2,atomnm);
      break;
    case eOBond:
      ps_draw_bond(ps,obj->ai,obj->aj,vec2,x,atomnm,size,bBalls);
      break;
    case eOHBond:
      if (bShowHydro)
	ps_draw_bond(ps,obj->ai,obj->aj,vec2,x,atomnm,size,bBalls);
      break;
    default:
      break;
    }
  }
}

static void v4_to_iv2(vec4 x4,iv2 v2,int x0,int y0,real sx,real sy)
{
  real inv_z;

  inv_z=1.0/x4[ZZ];
  v2[XX]=x0+sx*x4[XX]*inv_z;
  v2[YY]=y0-sy*x4[YY]*inv_z;
}

static void draw_box(t_psdata ps,t_3dview *view,matrix box,
		     int x0,int y0,real sx,real sy)
{
  int  ivec[8][4] =  { 
    { 0,0,0,1 }, { 1,0,0,1 }, { 1,1,0,1 }, { 0,1,0,1 },
    { 0,0,1,1 }, { 1,0,1,1 }, { 1,1,1,1 }, { 0,1,1,1 }
  };
  int  bonds[12][2] = {
    { 0,1 }, { 1,2 }, { 2,3 }, { 3,0 }, 
    { 4,5 }, { 5,6 }, { 6,7 }, { 7,4 },
    { 0,4 }, { 1,5 }, { 2,6 }, { 3,7 }
  };
  int  i,j;
  rvec corner[8];
  vec4 x4;
  iv2  vec2[12];

  for (i=0; (i<8); i++) {
    for (j=0; (j<DIM); j++)
      corner[i][j] = ivec[i][j]*box[j][j];
    m4_op(view->proj,corner[i],x4);
    v4_to_iv2(x4,vec2[i],x0,y0,sx,sy);
  }
  ps_color(ps,0,0,0.5);
  for (i=0; (i<12); i++)
    ps_line(ps,
	    vec2[bonds[i][0]][XX],vec2[bonds[i][0]][YY],
	    vec2[bonds[i][1]][XX],vec2[bonds[i][1]][YY]);
}

void ps_draw_mol(t_psdata ps,t_manager *man)
{
  static char tstr[2][20];
  static int  ntime=0;
  t_windata *win;
  t_3dview  *view;
  t_molwin  *mw;
  int       i,x0,y0,nvis;
  iv2       *vec2;
  real      sx,sy;
  vec4      x4;

  if (!man->status)
    return;

  view=man->view;
  mw=man->molw;

  win=&(mw->wd);

  vec2=man->ix;
  x0=win->width/2;
  y0=win->height/2;
  sx=win->width/2*view->sc_x;
  sy=win->height/2*view->sc_y;

  init_pbc(man->box);

  /* create_visibility(man); */

  for(i=0; (i<man->natom); i++) {
    if (man->bVis[i]) {
      m4_op(view->proj,man->x[i],x4);
      man->zz[i]=x4[ZZ];
      v4_to_iv2(x4,vec2[i],x0,y0,sx,sy);
    }
  }
  set_sizes(man,sx,sy);
  
  z_fill (man,man->zz);
  
  /* Start drawing 
  XClearWindow(x11->disp,win->self); */

  /* Draw Time 
  sprintf(tstr[ntime],"Time: %.3 ps",man->time);
  if (strcmp(tstr[ntime],tstr[1-ntime]) != 0) {
    set_vbtime(x11,man->vbox,tstr[ntime]);
    ntime=1-ntime;
  }*/

  if (mw->boxtype != esbNone)
    draw_box(ps,view,man->box,x0,y0,sx,sy);

  /* Should sort on Z-Coordinates here! */
  nvis=filter_vis(man);
  if (nvis && man->bSort)
    qsort(man->obj,nvis,sizeof(man->obj[0]),compare_obj);
  
  /* Draw the objects */
  ps_draw_objects(ps,
		  nvis,man->obj,man->ix,man->x,man->top.atoms.atomname,
		  man->size,
		  mw->bShowHydrogen,mw->bond_type,man->bPlus);

  /* Draw the labels */
  ps_color(ps,0,0,0);
  for(i=0; (i<man->natom); i++) 
    if (man->bLabel[i] && man->bVis[i]) 
      ps_text(ps,vec2[i][XX]+2,vec2[i][YY]-2,man->szLab[i]);
}

