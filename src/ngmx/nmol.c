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
 * Great Red Oystrich Makes All Chemists Sane
 */
static char *SRCID_nmol_c = "$Id$";

#include <math.h>
#include "sysstuff.h"
#include "string.h"
#include "smalloc.h"
#include "macros.h"
#include "xutil.h"
#include "3dview.h"
#include "fatal.h"
#include "buttons.h"
#include "manager.h"
#include "nmol.h"

#define MSIZE 4

static bool MWCallBack(t_x11 *x11,XEvent *event, Window w, void *data)
{
  t_molwin *mw;
  Window   To;
  XEvent   letter;

  mw=(t_molwin *)data;
  To=mw->wd.Parent;
  letter.type=ClientMessage;
  letter.xclient.display=x11->disp;
  letter.xclient.window=To;
  letter.xclient.message_type=0;
  letter.xclient.format=32;
  switch(event->type) {
  case Expose:
    /* Do Not draw anything, but signal parent instead, he will
     * coordinate drawing.
     */
    letter.xclient.data.l[0]=IDDRAWMOL;
    letter.xclient.data.l[1]=Button1;
    XSendEvent(x11->disp,To,True,0,&letter);
    break;
  case ButtonPress:
#ifdef DEBUG
    printf("Molwindow: Buttonpress\n");
#endif
    letter.xclient.data.l[0]=IDLABEL;
    letter.xclient.data.l[1]=(long)event->xbutton.button;
    letter.xclient.data.l[2]=event->xbutton.x;
    letter.xclient.data.l[3]=event->xbutton.y;
    XSendEvent(x11->disp,To,True,0,&letter);
    break;
  case ConfigureNotify:
    mw->wd.width=event->xconfigure.width;
    mw->wd.height=event->xconfigure.height;
    break;
  default:
    break;
  }
  return FALSE;
}

void set_def (t_molwin *mw)
{
  mw->bShowHydrogen=TRUE;
  mw->bond_type=eBFat;
  mw->bBoxSelect=TRUE;
}

t_molwin *init_mw(t_x11 *x11,Window Parent,
		  int x,int y,int width,int height,
		  ulong fg,ulong bg)
{
  t_molwin *mw;

  snew(mw,1);
  set_def(mw);
  
  InitWin(&mw->wd,x,y,width,height,1,"Mol Window");

  mw->wd.Parent=Parent;
  mw->wd.self=XCreateSimpleWindow(x11->disp,Parent,x,y,width,height,1,fg,bg);
  x11->RegisterCallback(x11,mw->wd.self,Parent,MWCallBack,mw);
  x11->SetInputMask(x11,mw->wd.self,
		    ExposureMask | StructureNotifyMask |
		    ButtonPressMask);
  return mw;
}

void map_mw(t_x11 *x11,t_molwin *mw)
{
  XMapWindow(x11->disp,mw->wd.self);
}

bool toggle_hydrogen(t_x11 *x11,t_molwin *mw)
{
  mw->bShowHydrogen=!mw->bShowHydrogen;
  ExposeWin(x11->disp,mw->wd.self);

  return mw->bShowHydrogen;
}

void set_bond_type(t_x11 *x11,t_molwin *mw,int bt)
{ 
  if (bt != mw->bond_type) {
    mw->bond_type=bt;
    ExposeWin(x11->disp,mw->wd.self);
  }
}

bool toggle_box (t_x11 *x11,t_molwin *mw)
{ 
  mw->bBoxSelect =! mw->bBoxSelect;
  ExposeWin(x11->disp,mw->wd.self);

  return mw->bBoxSelect;
}

void done_mw(t_x11 *x11,t_molwin *mw)
{
  x11->UnRegisterCallback(x11,mw->wd.self);
  sfree(mw);
}

static void draw_atom(Display *disp,Window w,GC gc,
		      atom_id ai,iv2 vec2[],ulong col[],int size[],
		      bool bBall,bool bPlus)
{
  int xi,yi;
  
  xi=vec2[ai][XX];
  yi=vec2[ai][YY];
  XSetForeground(disp,gc,col[ai]);
  if (bBall) {
    XFillCircle(disp,w,gc,xi,yi,size[ai]-1);
    XSetForeground(disp,gc,BLACK);
    XDrawCircle(disp,w,gc,xi,yi,size[ai]);
    /*    XSetForeground(disp,gc,WHITE);
	  XFillCircle(disp,w,gc,xi+4,yi-4,4); */
  }
  else if (bPlus) {
    XDrawLine(disp,w,gc,xi-MSIZE,yi,xi+MSIZE+1,yi);
    XDrawLine(disp,w,gc,xi,yi-MSIZE,xi,yi+MSIZE+1);
  }
  else 
    XDrawLine(disp,w,gc,xi-1,yi,xi+1,yi);
  
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

static bool local_pbc_dx(rvec x1, rvec x2)
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

static void draw_bond(Display *disp,Window w,GC gc,
		      atom_id ai,atom_id aj,iv2 vec2[],
		      rvec x[],ulong col[],int size[],bool bBalls)
{
  ulong   ic,jc;
  int     xi,yi,xj,yj;
  int     xm,ym;

  if (bBalls) {
    draw_atom(disp,w,gc,ai,vec2,col,size,TRUE,FALSE);
    draw_atom(disp,w,gc,aj,vec2,col,size,TRUE,FALSE);
  }
  else {
    if (local_pbc_dx(x[ai],x[aj])) {
      ic=col[ai];
      jc=col[aj];
      xi=vec2[ai][XX];
      yi=vec2[ai][YY];
      xj=vec2[aj][XX];
      yj=vec2[aj][YY];
      
      if (ic != jc) {
	xm=(xi+xj) >> 1;
	ym=(yi+yj) >> 1;
      
	XSetForeground(disp,gc,ic);
	XDrawLine(disp,w,gc,xi,yi,xm,ym);
	XSetForeground(disp,gc,jc);
	XDrawLine(disp,w,gc,xm,ym,xj,yj);
      }
      else {
	XSetForeground(disp,gc,ic);
	XDrawLine(disp,w,gc,xi,yi,xj,yj);
      }
    }
  }
}

int compare_obj(const void *a,const void *b)
{ 
  t_object *oa,*ob;
  real     z;
  
  oa=(t_object *)a;
  ob=(t_object *)b;

  z=oa->z-ob->z;

  if (z < 0)
    return 1;
  else if (z > 0)
    return -1;
  else 
    return 0;
}

void create_visibility(t_manager *man)
{ 
  t_object *obj;
  int       i;

  for(i=0,obj=man->obj; (i<man->nobj); i++,obj++)
    if (obj->eV != eVHidden) {
      man->bVis[obj->ai]=TRUE;
      switch (obj->eO) {
      case eOBond:
      case eOHBond:
	man->bVis[obj->aj]=TRUE;
	break;
      default:
	break;
      }
    }
}

void z_fill(t_manager *man, real *zz)
{ 
  t_object *obj;
  int       i;
  
  for(i=0,obj=man->obj; (i<man->nobj); i++,obj++)
    switch (obj->eO) {
    case eOSingle:
      obj->z = zz[obj->ai];
      break;
    case eOBond:
    case eOHBond:
      obj->z = (zz[obj->ai] + zz[obj->aj]) * 0.5;
      break;
    default:
      break;
    }
}

int filter_vis(t_manager *man)
{
  int      i,nobj,nvis,nhide;
  atom_id  ai;
  bool     bAdd,*bVis;
  t_object *obj;
  t_object *newobj;

  nobj=man->nobj;
  snew(newobj,nobj);
  obj=man->obj;
  bVis=man->bVis;
  nvis=0;
  nhide=nobj-1;
  for(i=0; (i<nobj); i++,obj++) {
    ai=obj->ai;
    bAdd=bVis[ai];
    if (bAdd)
      if (obj->eO != eOSingle)
	bAdd=bVis[obj->aj];
    if (bAdd)
      newobj[nvis++]=*obj;
    else
      newobj[nhide--]=*obj;
  }
  sfree(man->obj);
  man->obj=newobj;
  
  return nvis;
}

void draw_objects(Display *disp,Window w,GC gc,int nobj,
		  t_object objs[],iv2 vec2[],rvec x[],
		  ulong col[],int size[],bool bShowHydro,int bond_type,
		  bool bPlus)
{
  bool     bBalls;
  int      i;
  t_object *obj;

  bBalls=FALSE;
  switch (bond_type) {
  case eBThin:
    XSetLineAttributes(disp,gc,1,LineSolid,CapNotLast,JoinRound);
    break;
  case eBFat:
    XSetLineAttributes(disp,gc,3,LineSolid,CapNotLast,JoinRound);
    break;
  case eBVeryFat:
    XSetLineAttributes(disp,gc,5,LineSolid,CapNotLast,JoinRound);
    break;
  case eBSpheres:
    bBalls=TRUE;
    bPlus=FALSE;
    break;
  default:
    fatal_error(0,"Invalid bond_type selected: %d\n",bond_type);
    break;
  }
  for(i=0; (i<nobj); i++) {
    obj=&(objs[i]);
    switch (obj->eO) {
    case eOSingle:
      draw_atom(disp,w,gc,obj->ai,vec2,col,size,bBalls,bPlus);
      break;
    case eOBond:
      draw_bond(disp,w,gc,obj->ai,obj->aj,vec2,x,col,size,bBalls);
      break;
    case eOHBond:
      if (bShowHydro)
	draw_bond(disp,w,gc,obj->ai,obj->aj,vec2,x,col,size,bBalls);
      break;
    default:
      break;
    }
  }
  XSetLineAttributes(disp,gc,1,LineSolid,CapNotLast,JoinRound);
}

static void v4_to_iv2(vec4 x4,iv2 v2,int x0,int y0,real sx,real sy)
{
  real inv_z;

  inv_z=1.0/x4[ZZ];
  v2[XX]=x0+sx*x4[XX]*inv_z;
  v2[YY]=y0-sy*x4[YY]*inv_z;
}

static void draw_box(t_x11 *x11,Window w,t_3dview *view,matrix box,
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
  XSetForeground(x11->disp,x11->gc,YELLOW);
  for (i=0; (i<12); i++)
    XDrawLine(x11->disp,w,x11->gc,
	      vec2[bonds[i][0]][XX],vec2[bonds[i][0]][YY],
	      vec2[bonds[i][1]][XX],vec2[bonds[i][1]][YY]);
}

void set_sizes(t_manager *man,real sx,real sy)
{
  int  i;

  for(i=0; (i<man->natom); i++)
    if (man->bVis[i])
      man->size[i]=180*man->vdw[i];
}

void draw_mol(t_x11 *x11,t_manager *man)
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

  if (man->status == -1)
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
  
  /* Start drawing */
  XClearWindow(x11->disp,win->self);

  /* Draw Time */
  sprintf(tstr[ntime],"Time: %.3f ps",man->time);
  if (strcmp(tstr[ntime],tstr[1-ntime]) != 0) {
    set_vbtime(x11,man->vbox,tstr[ntime]);
    ntime=1-ntime;
  }

  if (mw->bBoxSelect)
    draw_box(x11,win->self,view,man->box,x0,y0,sx,sy);

  /* Should sort on Z-Coordinates here! */
  nvis=filter_vis(man);
  if (nvis && man->bSort)
    qsort(man->obj,nvis,sizeof(man->obj[0]),compare_obj);
  
  /* Draw the objects */
  draw_objects(x11->disp,win->self,x11->gc,
	       nvis,man->obj,man->ix,man->x,man->col,man->size,
	       mw->bShowHydrogen,mw->bond_type,man->bPlus);

  /* Draw the labels */
  XSetForeground(x11->disp,x11->gc,WHITE);
  for(i=0; (i<man->natom); i++) {
    if (man->bLabel[i] && man->bVis[i]) {
      XDrawString(x11->disp,win->self,x11->gc,vec2[i][XX]+2,vec2[i][YY]-2,
		  man->szLab[i],strlen(man->szLab[i]));
    }
  }

  XSetForeground(x11->disp,x11->gc,x11->fg);
}

