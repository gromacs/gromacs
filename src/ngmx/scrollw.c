/*
 *       @(#) copyrgt.c 1.12 9/30/97
 *
 *       This source code is part of
 *
 *        G   R   O   M   A   C   S
 *
 * GROningen MAchine for Chemical Simulations
 *
 *            VERSION 2.0b
 * 
 * Copyright (c) 1990-1997,
 * BIOSON Research Institute, Dept. of Biophysical Chemistry,
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
 * Gromacs Runs On Most of All Computer Systems
 */
#include <sysstuff.h>
#include <Xstuff.h>
#include <xutil.h>
#include <smalloc.h>
#include <macros.h>
#include <futil.h>
#include <string2.h>

#define YSPACE 2

typedef struct {
  t_windata   wd;		/* Window structure			*/
  int         nlines,top;	/* Number of lines, current top line	*/
  char        **lines;		/* The strings				*/
  int         wheight,wwidth;	/* The size of the window in chars	*/
  XFontStruct *font;		/* Font					*/
  ulong       fg,bg;		/* Colours				*/
} t_scrollw;

static void calc_scrollw(t_scrollw *sw,int w,int h)
{
  sw->wd.width=w;
  sw->wd.height=h;
  sw->wheight=h/(YSPACE+XTextHeight(sw->font));
  sw->wwidth=w/XTextWidth(sw->font,"W",1);
}

static bool SWCallback(t_x11 *x11,XEvent *event,Window w,void *data)
{
  t_scrollw *sw;
  int       i,y,nl,barw,btop,bheight;
  real      h,frac;

  sw=(t_scrollw *)data;

  /* Calc some bar data */
  barw=20;
  h=XTextHeight(sw->font)+YSPACE;
  frac=min(1.0,((real)sw->wheight)/((real)sw->nlines));
  btop=((((real)sw->top)/((real)sw->nlines)))*(sw->wd.height);
  bheight=frac*sw->wd.height;
  bheight-=bheight/h;

  switch(event->type) {
  case Expose:
    nl=min(sw->nlines,sw->top+sw->wheight);
    y=0;

    XClearWindow(x11->disp,w);
#ifdef DEBUG
    printf("btop: %d, bheight: %d, frac: %e, h: %e\n",btop,bheight,frac,h);
#endif
    /* Draw the bar */
    XSetForeground(x11->disp,x11->gc,LIGHTGREY);
    XFillRectangle(x11->disp,w,x11->gc,2,btop+2,barw-4,bheight-4);
    XDrawLine(x11->disp,w,x11->gc,barw,0,barw,sw->wd.height);

    /* Draw the text */
    XSetForeground(x11->disp,x11->gc,sw->fg);
    for(i=sw->top; (i<nl); i++) {
      SpecialTextInRect(x11,sw->font,w,
			sw->lines[i],barw+2,y,sw->wd.width-barw-4,(int)h,
			eXLeft,eYCenter);
      y+=h;
    }
    XSetForeground(x11->disp,x11->gc,x11->fg);
    break;
  case ConfigureNotify:
    calc_scrollw(sw,event->xconfigure.width,event->xconfigure.height);
    break;
  case ButtonPress:
    if (event->xbutton.x < barw) {
      int y=event->xbutton.y;

      if (sw->nlines > sw->wheight) {
	if (y<btop) 
	  sw->top=max(0,sw->top-1);
	else if (y>btop+bheight) 
	  sw->top=min(sw->nlines-sw->wheight,sw->top+1);
	else
	  break;
	ExposeWin(x11->disp,sw->wd.self);
      }
    }
    break;
  default:
    break;
  }

  return FALSE;
}

t_scrollw *init_scrollw(t_x11 *x11,Window parent,int x,int y,int w,int h,
			ulong fg,ulong bg)
{
  t_scrollw *sw;

  snew(sw,1);

  InitWin(&sw->wd,x,y,w,h,1,"Scroll Window");
  sw->fg=fg;
  sw->bg=bg;
  sw->font=x11->font;
  sw->wd.self=XCreateSimpleWindow(x11->disp,parent,x,y,w,h,
				  sw->wd.bwidth,fg,bg);
  x11->RegisterCallback(x11,sw->wd.self,parent,SWCallback,sw);
  x11->SetInputMask(x11,sw->wd.self,ExposureMask | ButtonPressMask |
		    StructureNotifyMask);
  calc_scrollw(sw,w,h);

  return sw;
}

void show_scrollw(t_x11 *x11,t_scrollw *sw)
{
  XMapWindow(x11->disp,sw->wd.self);
}

char *tab2spc(char *buf)
{
  char *buf2;
  char *new;
  int  i,j;

  snew(buf2,8*strlen(buf)+1);
  for(i=j=0; (buf[i]!='\0'); i++)
    if (buf[i]=='\t') 
      do {
	buf2[j++]=' ';
      } while ((j % 8)!=0);
    else
      buf2[j++]=buf[i];
  buf2[j]='\0';
  new=strdup(buf2);
  sfree(buf2);
  return new;
}

void read_lines(FILE *in,t_scrollw *sw)
{
  char buf[1024];

  while (fgets2(buf,1023,in)) {
    sw->nlines++;
    srenew(sw->lines,sw->nlines);
    sw->lines[sw->nlines-1]=tab2spc(buf);
  }
}

void main(int argc, char *argv[])
{
  t_x11     *x11;
  t_scrollw *sw;

  if ((x11=GetX11(&argc,argv))==NULL) {
    fprintf(stderr,"No X!\n");
    exit(1);
  }
  sw=init_scrollw(x11,x11->root,0,0,600,200,WHITE,BLACK);
  read_lines(stdin,sw);
  show_scrollw(x11,sw);
  x11->MainLoop(x11);

  x11->CleanUp(x11);
}
