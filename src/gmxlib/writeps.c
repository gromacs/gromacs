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
#include <stdio.h>
#include "futil.h"
#include "fatal.h"
#include "copyrite.h"
#include "writeps.h"

#define MAXRGB 200
static int   nrgb=0;
static t_rgb rgb[MAXRGB];

FILE *ps_open(char *fn,real x1,real y1,real x2,real y2)
{
  FILE *ps;
  
  ps=ffopen(fn,"w");
  fprintf(ps,"%%!PS-Adobe-2.0 EPSF-1.2\n");
  fprintf(ps,"%%%%Creator: GROMACS\n");
  fprintf(ps,"%%%%Title: %s\n",fn);
  fprintf(ps,"%%%%BoundingBox: %g %g %g %g\n",x1,y1,x2,y2);
  fprintf(ps,"%%%%EndComments\n");
  fprintf(ps,"/m {moveto} bind def\n");
  fprintf(ps,"/l {lineto} bind def\n");
  fprintf(ps,"/rm {rmoveto} bind def\n");
  fprintf(ps,"/r  {rlineto} bind def\n");
  fprintf(ps,"/f {fill} bind def\n");
  fprintf(ps,"/s {stroke} bind def\n");

  if (nrgb > 0) {
    fprintf(stderr,"Warning: resetting color table in %s when opening %s\n",
	    __FILE__,fn);
    nrgb=0;
  }
      
  return ps;
}

static void ps_defcolor(FILE *ps,real r,real g,real b,char *cname)
{
  fprintf(ps,"/%s {%g %g %g setrgbcolor} bind def\n",cname,r,g,b);
}

static void ps_selcolor(FILE *ps,char *cname)
{
  fprintf(ps,"%s\n",cname);
}

static char *i2a(int i)
{
  static char buf[12];
  
  sprintf(buf,"C%d",i);
  
  return buf;
}

static int search_col(FILE *ps,real r,real g,real b)
{
  int  i;
  
  for(i=0; (i<nrgb); i++) {
    if ((rgb[i].r == r) && (rgb[i].g == g) && (rgb[i].b == b))
      return i;
  }
  if (nrgb < MAXRGB) {
    ps_defcolor(ps,r,g,b,i2a(nrgb));
    rgb[i].r=r;
    rgb[i].g=g;
    rgb[i].b=b;
    nrgb++;
    
    return nrgb-1;
  }
  else {
    fclose(ps);
    fatal_error(0,"Colormap full! Increase MAXRGB");
  }
    
  return -1;
}

void ps_color(FILE *ps,real r,real g,real b)
{
  ps_selcolor(ps,i2a(search_col(ps,r,g,b)));
}

void ps_rgb(FILE *ps,t_rgb *rgb)
{
  ps_color(ps,rgb->r,rgb->g,rgb->b);
}

void ps_lineto(FILE *ps,real x,real y)
{
  fprintf(ps,"%g %g l\n",x,y);
}

void ps_linerel(FILE *ps,real dx,real dy)
{
  fprintf(ps,"%g %g r\n",dx,dy);
}

void ps_moveto(FILE *ps,real x,real y)
{
  fprintf(ps,"%g %g m\n",x,y);
}

void ps_moverel(FILE *ps,real dx,real dy)
{
  fprintf(ps,"%g %g rm\n",dx,dy);
}

void ps_line(FILE *ps,real x1,real y1,real x2,real y2)
{
  ps_moveto(ps,x1,y1);
  ps_lineto(ps,x2,y2);
  fprintf(ps,"s\n");
}

static void do_box(FILE *ps,real x1,real y1,real x2,real y2)
{
  ps_moveto(ps,x1,y1);
  ps_linerel(ps,0,y2-y1);
  ps_linerel(ps,x2-x1,0);
  ps_linerel(ps,0,y1-y2);
  ps_linerel(ps,x1-x2,0);
}

void ps_box(FILE *ps,real x1,real y1,real x2,real y2)
{
  do_box(ps,x1,y1,x2,y2);
  fprintf(ps,"s\n");
}

void ps_fillbox(FILE *ps,real x1,real y1,real x2,real y2)
{
  do_box(ps,x1,y1,x2,y2);
  fprintf(ps,"f\n");
}

void ps_arc(FILE *ps,real x1,real y1,real rad,real a0,real a1)
{
  fprintf(ps,"%g %g %g %g %g arc s\n",x1,y1,rad,a0,a1);
}

void ps_fillarc(FILE *ps,real x1,real y1,real rad,real a0,real a1)
{
  fprintf(ps,"%g %g %g %g %g arc f\n",x1,y1,rad,a0,a1);
}

void ps_arcslice(FILE *ps,real xc,real yc,
		 real rad1,real rad2,real a0,real a1)
{
  fprintf(ps,"newpath %g %g %g %g %g arc %g %g %g %g %g arcn closepath s\n",
	  xc,yc,rad1,a0,a1,xc,yc,rad2,a1,a0);
}
  
void ps_fillarcslice(FILE *ps,real xc,real yc,
		     real rad1,real rad2,real a0,real a1)
{
  fprintf(ps,"newpath %g %g %g %g %g arc %g %g %g %g %g arcn closepath f\n",
	  xc,yc,rad1,a0,a1,xc,yc,rad2,a1,a0);
}
  
void ps_circle(FILE *ps,real x1,real y1,real rad)
{
  ps_arc(ps,x1,y1,rad,0,360);
}

char *fontnm[efontNR] = { 
  "Times-Roman","Times-Italic",     "Times-Bold",    "Times-BoldItalic",
  "Helvetica",  "Helvetica-Oblique","Helvetica-Bold","Helvetica-BoldOblique",
  "Courier",    "Courier-Oblique",  "Courier-Bold",  "Courier-BoldOblique"
};

void ps_font(FILE *ps,int font,real size)
{
  
  if ((font < 0) || (font > efontNR)) {
    fprintf(stderr,"Invalid Font: %d, using %s\n",font,fontnm[0]);
    font=0;
  }
  fprintf(ps,"/%s findfont\n",fontnm[font]);
  fprintf(ps,"%g scalefont setfont\n",size);
}

void ps_strfont(FILE *ps,char *font,real size)
{
  fprintf(ps,"/%s findfont\n",font);
  fprintf(ps,"%g scalefont setfont\n",size);
}

void ps_text(FILE *ps,real x1,real y1,char *str)
{
  ps_moveto(ps,x1,y1);
  fprintf(ps,"(%s) show\n",str);
}

void ps_rotate(FILE *ps,bool bPlus)
{
  if (bPlus) 
    fprintf(ps,"612.5 0 translate 90 rotate\n");
  else
    fprintf(ps,"-90 rotate -612.5 0 translate\n");
}

void ps_ctext(FILE *ps,real x1,real y1,char *str,int expos)
{
  if (expos == eXLeft) {
    ps_text(ps,x1,y1,str);
    return;
  }
  ps_moveto(ps,x1,y1);
  fprintf(ps,"(%s) stringwidth\n",str);
  switch (expos) {
  case eXLeft:
    fprintf(ps,"exch 0 exch pop exch\n");
    break;
  case eXCenter:
    fprintf(ps,"exch 2 div neg exch\n");
    break;
  case eXRight:
    fprintf(ps,"exch neg exch\n");
    break;
  default:
    fatal_error(0,"invalid position index (expos=%d)",expos);
  }
  fprintf(ps,"rmoveto (%s) show\n",str);
}

void ps_translate(FILE *ps,real x,real y)
{
  fprintf(ps,"%g %g translate\n",x,y);
}

static int ostack=0;

void ps_setorigin(FILE *ps)
{
  fprintf(ps,"currentpoint dup 3 -1 roll dup 4 1 roll exch translate\n");
  ostack++;
}

void ps_unsetorigin(FILE *ps)
{
  if (ostack <= 0)
    fatal_error(0,"No origin on stack!\n");
  fprintf(ps,"neg exch neg exch translate\n");
  ostack--;
}

void ps_close(FILE *ps)
{
  fprintf(ps,"%%showpage\n");
  fprintf(ps,"%%%%EOF\n");
  fclose(ps);
}

void viewps(char *fn)
{
  char buf[256];
  
  sprintf(buf,"ghostview %s &",fn);
  system(buf);
}

void ps_comment(FILE *ps,char *s)
{
  fprintf(ps,"%%%% %s\n",s);
}
