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
static char *SRCID_xutil_c = "$Id$";

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <Xstuff.h>
#include <xutil.h>
#include "assert.h"
#include "smalloc.h"
#include "typedefs.h"
#include "string2.h"

int CheckWin(Window win,char *file, int line)
{
  typedef struct {
    int  n;
    char *s;
  } t_winerr;
  t_winerr winerr[] = {
    { BadAlloc,  "Bad Alloc" },
    { BadColor,  "Bad Color" },
    { BadCursor, "Bad Cursor"},
    { BadMatch,  "Bad Match" },
    { BadPixmap, "Bad Pixmap"},
    { BadValue,  "Bad Value" },
    { BadWindow, "Bad Window"}
  };  
#define NERR (sizeof(winerr)/sizeof(winerr[0]))  
  int i;

  for(i=0; (i<NERR); i++)
    if (win==winerr[i].n) {
      fprintf(stderr,"%s",winerr[i].s);
      break;
    }
  if (i==NERR) return 1;

  fprintf(stderr," in file %s, line %d\n",file,line);
  return 0;
}

void LightBorder(Display *disp, Window win, unsigned long color)
{
  XSetWindowAttributes attributes;

  attributes.border_pixel = color;
  XChangeWindowAttributes(disp,win,CWBorderPixel,&attributes);
}

void SpecialTextInRect(t_x11 *x11,XFontStruct *font,Drawable win,
		       char *s,int x,int y,int width,int height,
		       eXPos eX,eYPos eY)
{
  int         fw,fh,x0,y0;
  XFontStruct *f;

  if (font) {
    XSetFont(x11->disp,x11->gc,font->fid);
    f=font;
  }
  else
    f=x11->font;

  fw=XTextWidth(f,s,strlen(s));
  fh=XTextHeight(f);
  switch (eX) {
  case eXLeft:
    x0=x;
    break;
  case eXRight:
    x0=x+width-fw;
    break;
  case eXCenter:
  default:
    x0=x+(width-fw)/2;
    break;
  }
  switch (eY) {
  case eYTop:
    y0=y+f->ascent;
    break;
  case eYBottom:
    y0=y+height-f->descent;
    break;
  case eYCenter:
  default:
    y0=y+(height-fh)/2+f->ascent;
    break;
  }
  XDrawString(x11->disp,win,x11->gc,x0,y0,s,strlen(s));
  if (font)
    XSetFont(x11->disp,x11->gc,x11->font->fid);
}

void TextInRect(t_x11 *x11,Drawable win,
		char *s,int x,int y,int width,int height,
		eXPos eX,eYPos eY)
{
  SpecialTextInRect(x11,NULL,win,s,x,y,width,height,eX,eY);
}

void TextInWin(t_x11 *x11, t_windata *win,
	       char *s, eXPos eX, eYPos eY)
{
  TextInRect(x11,win->self,s,0,0,win->width,win->height,eX,eY);
}

void InitWin(t_windata *win, int x0,int y0, int w, int h, int bw, char *text)
{
  win->self=0;
  win->color=0;
  win->x=x0;
  win->y=y0;
  win->width=w;
  win->height=h;
  win->bwidth=bw;
  win->bFocus=FALSE;
  win->cursor=0;
  if (text)
    win->text=strdup(text);
  else
    win->text=NULL;
#ifdef DEBUG
  printf("%s: %d x %d at %d, %d\n",text,w,h,x0,y0);
#endif  
}

void FreeWin(Display *disp, t_windata *win)
{
  if (win->text)
    sfree(win->text);
  if (win->cursor)
    XFreeCursor(disp,win->cursor);
}

void ExposeWin(Display *disp,Window win)
{
  XEvent event;

  event.type = Expose;
  event.xexpose.send_event=True;
  event.xexpose.window = win;
  event.xexpose.x=0;
  event.xexpose.y=0;
  event.xexpose.width=1000;
  event.xexpose.height=1000;
  event.xexpose.count = 0;
  XSendEvent(disp,win,False,ExposureMask,&event);
}

void XDrawRoundRect(Display *disp, Window win, GC gc, 
		    int x, int y, int w, int h)
{
#define RAD (OFFS_X/2)
#define SetPoint(pn,x0,y0) pn.x=x0; pn.y=y0

  if ((w<10) || (h<10))
    XDrawRectangle(disp,win,gc,x,y,w,h);
  else {
    XPoint p[9];

    SetPoint(p[0],x+RAD,y);
    SetPoint(p[1],w-2*RAD,0);
    SetPoint(p[2],RAD,RAD);
    SetPoint(p[3],0,h-2*RAD);
    SetPoint(p[4],-RAD,RAD);
    SetPoint(p[5],2*RAD-w,0);
    SetPoint(p[6],-RAD,-RAD);
    SetPoint(p[7],0,2*RAD-h);
    SetPoint(p[8],RAD,-RAD);
    XDrawLines(disp,win,gc,p,9,CoordModePrevious);
  }
}

void RoundRectWin(Display *disp, GC gc, t_windata *win, 
		  int offsx, int offsy,unsigned long color)
{
  XSetLineAttributes(disp,gc,1,LineOnOffDash,CapButt,JoinRound);
  XSetForeground(disp,gc,color);
  XDrawRoundRect(disp,win->self,gc,offsx,offsy,
		 win->width-2*offsx-1,win->height-2*offsy-1);
  XSetLineAttributes(disp,gc,1,LineSolid,CapButt,JoinRound);
}

void RectWin(Display *disp, GC gc, t_windata *win, unsigned long color)
{
  int bw=1; /*2*w.bwidth;*/

  XSetForeground(disp,gc,color);
  XDrawRoundRect(disp,win->self,gc,0,0,win->width-bw,win->height-bw);
}

typedef struct t_mpos {
  int    x,y;
  struct t_mpos *prev;
} t_mpos;

static t_mpos *mpos=NULL;

void PushMouse(Display *disp, Window dest, int x, int y)
{
  Window root,child;
  int    root_x,root_y;
  int    win_x,win_y;
  unsigned int   keybut;
  t_mpos *newpos;
  
  snew(newpos,1);
  XQueryPointer(disp,DefaultRootWindow(disp),&root,&child,&root_x,&root_y,
		&win_x,&win_y,&keybut);
  newpos->x=root_x;
  newpos->y=root_y;
  newpos->prev=mpos;
  mpos=newpos;
  XWarpPointer(disp,None,dest,0,0,0,0,x,y);
#ifdef DEBUG
  fprintf(stderr,"Pushmouse %d, %d\n",x,y);
#endif
}

void PopMouse(Display *disp)
{
  t_mpos *old;
  
  old=mpos;
  if (!old) return;

  XWarpPointer(disp,None,DefaultRootWindow(disp),0,0,0,0,old->x,old->y);
#ifdef DEBUG
  fprintf(stderr,"Popmouse %d, %d\n",old->x,old->y);
#endif
  mpos=old->prev;
  sfree(old);
}

bool HelpPressed(XEvent *event)
{
#define BUFSIZE 24
  char           buf[BUFSIZE+1];
  XComposeStatus compose;
  KeySym         keysym;

  (void)XLookupString(&(event->xkey),buf,BUFSIZE,&keysym,&compose);

  return (keysym == XK_F1);
}

bool GrabOK(FILE *out, int err)
{
  switch (err) {
  case GrabSuccess:
    return TRUE;
  case GrabNotViewable:
    fprintf(out,"GrabNotViewable\n"); break;
  case AlreadyGrabbed:
    fprintf(out,"AlreadyGrabbed\n"); break;
  case GrabFrozen:
    fprintf(out,"GrabFrozen\n"); break;
  case GrabInvalidTime:
    fprintf(out,"GrabInvalidTime\n"); break;
  default:
    break;
  }
  return FALSE;
}

