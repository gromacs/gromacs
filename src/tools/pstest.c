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
 * GROtesk MACabre and Sinister
 */
#include "string2.h"
#include "simple.h"
#include "macros.h"
#include "writeps.h"
#include "futil.h"
#include "fatal.h"
#include "smalloc.h"

typedef struct {
  real r,g,b;
} t_rgb;

static  char  dssp_val[]=" EBSTHIG";
static  char *dssp_str[]={ 
  "None", "B-Sheet", "Bend", "S", "Turn", "Helix", "5-Helix", "3-Helix"
};
#define NVAL asize(dssp_val)-1


static t_rgb pl_int(char SSTP)
{
  t_rgb cols[] = {
    1,    1,    1,
    1,    0,    0,
    0.8,  0.8,  0,
    0,    0.5,  0.5,
    0,    1,    1,
    0,    0,    1,
    0.5,  0,    0.5,
    0.5,  0.5,  0.5
  };
  int Col;
  
  for(Col=0; (Col < NVAL); Col++)
    if (dssp_val[Col] == SSTP)
      return cols[Col];
      
  return cols[0];
}

void ps_dssp_leg(FILE *ps,int x0,int y0)
{
#define  FONTSIZE 15
#define  DDD      2
#define  BOXHH  (FONTSIZE+2*DDD)
  int   i,yy;
  t_rgb col;
  
  for(i=0; (i<NVAL); i++) {
    col=pl_int(dssp_val[i]);
    ps_color(ps,col.r,col.g,col.b);
    yy=y0+i*BOXHH;
    ps_fillbox(ps,x0+DDD,yy+DDD,x0+DDD+FONTSIZE,yy+BOXHH-DDD);
    ps_color(ps,0,0,0);
    ps_box(ps,x0+DDD,yy+DDD,x0+DDD+FONTSIZE,yy+BOXHH-DDD);
    ps_ctext(ps,x0+BOXHH+2*DDD,yy+DDD+FONTSIZE/3,dssp_str[i],eXLeft);
  }
}

void ps_dssp(char *outf,int nres,int nframes,real t[],char *ss[],
	     bool bLabels,int ymod,int r0,int xmod,bool bLeg)
{
#define  BOXS     4
#define  HBOXS    (BOXS/2)
#define  TICK     8
  static char *res="Residue";
  FILE   *out;
  int    W,H;
  int    font;
  int    tw,th;
  int    x,y,w,h,x0,y0,xx,yy;
  t_rgb  col;
  char   buf[20];
  
  /* Set up size of box for nice colors */
  w=nres*BOXS;
  h=nframes*BOXS;
  tw=12;
  th=18;
  
  /* Set up bounding box */
  W=w+100;
  H=h+100;
  
  /* Layout is 
   *
   *        Residue
   *          nr
   *    0
   * T  1
   * i  2
   * m  3
   * e  4
   *
   */
  /* Border is */
  x0=W-w-1;
  y0=20;
  out=ps_open(outf,0,0,W+6*BOXS,h+6*BOXS);
  
  font=efontHELV;
  ps_font(out,font,FONTSIZE);
  ps_color(out,0,0,0);
  ps_box(out,x0-1,y0-1,x0+w+1,y0+h+1);
  if (bLabels) {
    yy=y0+h+1;
    ps_ctext(out,x0+w/2,yy+TICK+th,res,eXCenter);
    for(x=0; (x<nres); x++) {
      if ((x==0) || (((x+r0) % xmod)==0)) {
	xx=x0+HBOXS+(x*BOXS);
	ps_line(out,xx,yy,xx,yy+TICK);
	sprintf(buf,"%d",x+r0);
	ps_ctext(out,xx,yy+TICK,buf,eXCenter);
      }
    }
  }
  for(y=0; (y<nframes); y++) {
    yy=y0+h-(y+1)*BOXS;
    for(x=0; (x<nres); x++) {
      col=pl_int(ss[y][x]);
      ps_color(out,col.r,col.g,col.b);
      ps_fillbox(out,x0+x*BOXS,yy,x0+(x+1)*BOXS,yy+BOXS);
    }
    if (bLabels) {
      if ((y % ymod) == 0) {
	ps_color(out,0,0,0);
	xx=x0-TICK;
	ps_line(out,xx,yy+HBOXS,x0-1,yy+HBOXS);
	sprintf(buf,"%4.0f",t[y]);
	ps_ctext(out,xx-HBOXS,yy+HBOXS-FONTSIZE/3,buf,eXRight);
      }
    }
  }
  if (bLeg)
    ps_dssp_leg(out,x0+w+FONTSIZE,y0);
  ps_rotate(out,TRUE);
  ps_ctext(out,y0+h/2,612.5-2*FONTSIZE,"Time (ps)",eXCenter);
  ps_close(out);
}

void dssp_test()
{
  char *ss[] = {
    "HHHHHHHHHH",
    "HHHGGHHHIH",
    "HHHEEEETTT",
    "EEE HH  EE",
    "HHHGGHHHIH",
    "GIHSBT   H",
    "HHHEEEETTT",
    "HHHSSHHHIH",
    "HHHGGGETTT",
    "HHHHHHHHHH",
    "HHHGGHHHIH",
    "HHHEEEETTT",
    "EEE HH  EE",
    "HHHGGHHHIH",
    "GIHSBT   H",
    "HHHEEEETTT",
    "HHHSSHHHIH",
    "HHHGGGETTT",
    "HHHHHHHHHH",
    "HHHGGHHHIH",
    "HHHEEEETTT",
    "EEE HH  EE",
    "HHHGGHHHIH",
    "GIHSBT   H",
    "HHHEEEETTT",
    "HHHSSHHHIH"
  };
#define NF   asize(ss)
#define NRES 10/*asize(ss[0])*/
  real time[NF];
  int  i;
  
  for(i=0; (i<NF); i++)
    time[i]=10*i;
  ps_dssp("dssp.ps",NRES,NF,time,ss,TRUE,5,21,10,FALSE);
}

void real_dssp_test(char *fn)
{
  FILE   *in;
  char   buf[2048],sss[2048];
  real   *time=NULL;
  char   **ss=NULL;
  double tt;
  int    nres=0,nframes=0;
    
  in=ffopen(fn,"r");
  while (fgets2(buf,2048-1,in) != NULL) {
    sscanf(buf,"%lf %s",&tt,sss);
    if (nres==0)
      nres=strlen(sss)-2;
    else if (nres != strlen(sss)-2)
      fatal_error(0,"Your input file sucks");
    srenew(time,nframes+1);
    srenew(ss,nframes+1);
    time[nframes]=tt;
    sss[nres+2]='\0';
    ss[nframes]=strdup(sss+1);
    nframes++;
  }
  fclose(in);
  
  ps_dssp("dssp.ps",nres,nframes,time,ss,TRUE,5,21,10,TRUE);
}

void font_test()
{
  FILE *ps;
  int i;
  
  ps=ps_open("fonts.ps",0,0,600,500);
  for(i=0; (i<efontNR); i++) {
    ps_font(ps,i,20);
    ps_text(ps,10+i,40*(i+1),fontnm[i]);
  }

  ps_close(ps);
}

void align_test()
{
  FILE *ps;
#define L 100
#define R 500
#define C ((L+R)/2)
  char *x[]={"Center", "Left", "Right" };
  char *y[]={"Center", "Top",  "Bottom"  };
  int  ii[]={C,L,R};
  char buf[128];
  int  i,j;
  
  ps=ps_open("align.ps",L-100,L-100,R+100,R+100);
  ps_font(ps,0,20);
  ps_color(ps,0,0,0);
  ps_box(ps,L,L,R,R);
  ps_line(ps,L,C,R,C);
  ps_line(ps,C,L,C,R);
  for(i=0; (i<3); i++) {
    for(j=0; (j<3); j++) {
      ps_rotate(ps,TRUE);
      sprintf(buf,"%s-%s",x[i],y[j]);
      ps_ctext(ps,ii[i],ii[j],buf,i);
    }
  }
  ps_close(ps);
}

void main()
{
  /*align_test();*/
  real_dssp_test("ss.xvg");
}
