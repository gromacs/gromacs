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
static char *SRCID_xpm2ps_c = "$Id$";

#include <math.h>
#include "string2.h"
#include "copyrite.h"
#include "typedefs.h"
#include "macros.h"
#include "statutil.h"
#include "writeps.h"
#include "futil.h"
#include "fatal.h"
#include "smalloc.h"
#include "string2.h"
#include "matio.h"

#define FUDGE 1.2
#define DDD   2

typedef struct {
  real major;
  real minor;
  real offset;
  bool first;
  int  lineatzero;
  real majorticklen;
  real minorticklen;
  char label[STRLEN];
  real fontsize;
  char font[STRLEN];
  real tickfontsize;
  char tickfont[STRLEN];
} t_axisdef;

typedef struct {
  int  bw;
  real linewidth;
  real xoffs,yoffs;
  bool bTitle;
  real titfontsize;
  char titfont[STRLEN];
  bool legend;
  real legfontsize;
  char legfont[STRLEN];
  char leglabel[STRLEN];
  char leg2label[STRLEN];
  real xboxsize;
  real yboxsize;
  real boxspacing;
  real boxlinewidth;
  real ticklinewidth;
  real zerolinewidth;
  t_axisdef X,Y;
} t_psrec;

void get_params(char *mpin,char *mpout,t_psrec *psr)
{
  static char *bools[BOOL_NR+1]  = { "no", "yes", NULL };
  /* this must correspond to t_rgb *linecolors[] below */
  static char *colors[] = { "none", "black", "white", NULL };
  t_inpfile *inp;
  char      *tmp;
  int       ninp;
  
  inp=read_inpfile(mpin,&ninp);
  ETYPE("black&white",    psr->bw,             bools);
  RTYPE("linewidth",      psr->linewidth,      2.0);
  STYPE("titlefont",      psr->titfont,        "Helvetica");
  RTYPE("titlefontsize",  psr->titfontsize,    20.0);
  ETYPE("legend",         psr->legend,         bools);
  STYPE("legendfont",     psr->legfont,        psr->titfont);
  STYPE("legendlabel",    psr->leglabel,       "");
  STYPE("legend2label",   psr->leg2label,      psr->leglabel);
  RTYPE("legendfontsize", psr->legfontsize,    14.0);
  RTYPE("xbox",           psr->xboxsize,       2.0);
  RTYPE("ybox",           psr->yboxsize,       psr->xboxsize);
  RTYPE("matrixspacing",  psr->boxspacing,     20.0);
  RTYPE("xoffset",        psr->xoffs,          0.0);
  RTYPE("yoffset",        psr->yoffs,          psr->xoffs);
  RTYPE("boxlinewidth",   psr->boxlinewidth,   psr->linewidth);
  RTYPE("ticklinewidth",  psr->ticklinewidth,  psr->linewidth);
  RTYPE("zerolinewidth",  psr->zerolinewidth,  psr->ticklinewidth);
  ETYPE("x-lineat0value", psr->X.lineatzero,   colors);
  RTYPE("x-major",        psr->X.major,        20.0);
  RTYPE("x-minor",        psr->X.minor,        5.0);
  RTYPE("x-firstmajor",   psr->X.offset,       0.0);
  ETYPE("x-majorat0",     psr->X.first,        bools);
  RTYPE("x-majorticklen", psr->X.majorticklen, 8.0);
  RTYPE("x-minorticklen", psr->X.minorticklen, 4.0);
  STYPE("x-label",        psr->X.label,        "");
  RTYPE("x-fontsize",     psr->X.fontsize,     16.0);
  STYPE("x-font",         psr->X.font,         psr->titfont);
  RTYPE("x-tickfontsize", psr->X.tickfontsize, 10.0);
  STYPE("x-tickfont",     psr->X.tickfont,     psr->X.font);
  ETYPE("y-lineat0value", psr->Y.lineatzero,   colors);
  RTYPE("y-major",        psr->Y.major,        psr->X.major);
  RTYPE("y-minor",        psr->Y.minor,        psr->X.minor);
  RTYPE("y-firstmajor",   psr->Y.offset,       psr->X.offset);
  ETYPE("y-majorat0",     psr->Y.first,        bools);
  RTYPE("y-majorticklen", psr->Y.majorticklen, psr->X.majorticklen);
  RTYPE("y-minorticklen", psr->Y.minorticklen, psr->X.minorticklen);
  STYPE("y-label",        psr->Y.label,        psr->X.label);
  RTYPE("y-fontsize",     psr->Y.fontsize,     psr->X.fontsize);
  STYPE("y-font",         psr->Y.font,         psr->X.font);
  RTYPE("y-tickfontsize", psr->Y.tickfontsize, psr->X.tickfontsize);
  STYPE("y-tickfont",     psr->Y.tickfont,     psr->Y.font);
  if (mpout)
    write_inpfile(mpout,ninp,inp);
}

t_rgb black={ 0.0, 0.0, 0.0 };
t_rgb white={ 1.0, 1.0, 1.0 };
#define BLACK (&black)
/* this must correspond to *colors[] in get_params */
t_rgb *linecolors[] = { NULL, &black, &white, NULL };

bool diff_maps(int nmap1,t_mapping *map1,int nmap2,t_mapping *map2)
{
  int i;
  bool bDiff,bColDiff=FALSE;

  if (nmap1 != nmap2) 
      bDiff=TRUE;
    else {
      bDiff=FALSE;
      for(i=0; i<nmap1; i++) {
	if (!matelmt_cmp(map1[i].code, map2[i].code)) bDiff=TRUE;
	if (strcmp(map1[i].desc,map2[i].desc) != 0) bDiff=TRUE;
	if ((map1[i].rgb.r!=map2[i].rgb.r) ||
	    (map1[i].rgb.g!=map2[i].rgb.g) ||
	    (map1[i].rgb.b!=map2[i].rgb.b))
	  bColDiff=TRUE;
      }
      if (!bDiff && bColDiff)
	fprintf(stderr,"Warning: two colormaps differ only in RGB value, using one colormap.\n");
    }
  
  return bDiff;
}
  
void leg_discrete(FILE *ps,real x0,real y0,char *label,
		  real fontsize,char *font,int nmap,t_mapping map[])
{
  int   i;
  real  yhh;
  real  boxhh;
  
  boxhh=fontsize+DDD;
  /* LANDSCAPE */
  ps_rgb(ps,BLACK);
  ps_strfont(ps,font,fontsize);
  yhh=y0+fontsize+3*DDD;
  if ((int)strlen(label) > 0)
    ps_ctext(ps,x0,yhh,label,eXLeft);
  ps_moveto(ps,x0,y0);
  for(i=0; (i<nmap); i++) {
    ps_setorigin(ps);
    ps_rgb(ps,&(map[i].rgb));
    ps_fillbox(ps,DDD,DDD,DDD+fontsize,boxhh-DDD);
    ps_rgb(ps,BLACK);
    ps_box(ps,DDD,DDD,DDD+fontsize,boxhh-DDD);
    ps_ctext(ps,boxhh+2*DDD,fontsize/3,map[i].desc,eXLeft);
    ps_unsetorigin(ps);
    ps_moverel(ps,DDD,-fontsize/3);
  }
}

void leg_continuous(FILE *ps,real x0,real x,real y0,char *label,
		    real fontsize,char *font,
		    int nmap,t_mapping map[])
{
  int   i;
  real  xx0;
  real  yhh,boxxh,boxyh;
  
  boxyh=fontsize;
  if (x<8*fontsize)
    x=8*fontsize;
  boxxh=(real)x/(real)nmap;
  if (boxxh>fontsize)
    boxxh=fontsize;
  
  /* LANDSCAPE */
  xx0=x0-(nmap*boxxh)/2.0;
  
  for(i=0; (i<nmap); i++) {
    ps_rgb(ps,&(map[i].rgb));
    ps_fillbox(ps,xx0+i*boxxh,y0,xx0+(i+1)*boxxh,y0+boxyh);
  }
  ps_strfont(ps,font,fontsize);
  ps_rgb(ps,BLACK);
  ps_box(ps,xx0,y0,xx0+nmap*boxxh,y0+boxyh);
  
  yhh=y0+boxyh+3*DDD;
  ps_ctext(ps,xx0+boxxh/2,yhh,map[0].desc,eXCenter);
  if ((int)strlen(label) > 0)
    ps_ctext(ps,x0,yhh,label,eXCenter);
  ps_ctext(ps,xx0+(nmap*boxxh)-boxxh/2,yhh,map[nmap-1].desc,eXCenter);
}

void leg_bicontinuous(FILE *ps,real x0,real x,real y0,char *label1,char *label2,
		      real fontsize,char *font,
		      int nmap1,t_mapping map1[],int nmap2,t_mapping map2[])
{
  real xx1,xx2,x1,x2;
  
  x1=x/(nmap1+nmap2)*nmap1;/* width of legend 1 */
  x2=x/(nmap1+nmap2)*nmap2;/* width of legend 2 */
  xx1=x0-(x2/2.0)-fontsize;/* center of legend 1 */
  xx2=x0+(x1/2.0)+fontsize;/* center of legend 2 */
  x1-=fontsize/2;/* adjust width */
  x2-=fontsize/2;/* adjust width */
  leg_continuous(ps,xx1,x1,y0,label1,fontsize,font,nmap1,map1);
  leg_continuous(ps,xx2,x2,y0,label2,fontsize,font,nmap2,map2);
}

static real box_height(t_matrix *mat,t_psrec *psr)
{
  return mat->ny*psr->yboxsize; 
}

static real box_dh(t_psrec *psr)
{
  return psr->boxspacing;
}

static real box_dh_top(t_psrec *psr)
{
  real dh;

  if (psr->bTitle)
    dh=2*psr->titfontsize;
  else
    dh=0;

  return  dh;
}

static bool box_do_all_x_maj_ticks(t_psrec *psr)
{
  return (psr->boxspacing>(1.5*psr->X.majorticklen));
}

static bool box_do_all_x_min_ticks(t_psrec *psr)
{
  return (psr->boxspacing>(1.5*psr->X.minorticklen));
}

static void draw_boxes(FILE *out,real x0,real y0,real w,
		       int nmat,t_matrix mat[],t_psrec *psr)
{
  char   buf[12];
  real   xxx;
  char   **xtick,**ytick;
  real   xx,yy,dy,xx00,yy00;
  int    i,j,x,y,strlength;
  
  /* Only necessary when there will be no y-labels */ 
  strlength = 0;
  
  /* Draw the box */
  ps_rgb(out,BLACK);
  ps_linewidth(out,psr->boxlinewidth);
  yy00=y0;
  for(i=0; (i<nmat); i++) {
    dy=box_height(&(mat[i]),psr);
    ps_box(out,x0-1,yy00-1,x0+w+1,yy00+dy+1);
    yy00+=dy+box_dh(psr)+box_dh_top(psr);
  }
  
  /* Draw the ticks on the axes */
  ps_linewidth(out,psr->ticklinewidth);
  xx00=x0-1;
  yy00=y0-1;
  for (i=0; (i<nmat); i++) {
    if (mat[i].axis_x==NULL) {
      snew(mat[i].axis_x,mat[i].nx);
      for(j=0; (j<mat[i].nx); j++)
	mat[i].axis_x[j]=j;
    }
    /* Make labels for x axis */
    snew(xtick,mat[i].nx);
    for(j=0; (j<mat[i].nx); j++) {
      sprintf(buf,"%g",mat[i].axis_x[j]);
      xtick[j]=strdup(buf);
    }
    ps_strfont(out,psr->X.tickfont,psr->X.tickfontsize);
    for(x=0; (x<mat[i].nx); x++) {
      xx=xx00+(x+0.7)*psr->xboxsize;
      if ( ( bRmod(mat[i].axis_x[x] - psr->X.offset, psr->X.major) || 
	     (psr->X.first && (x==0))) &&
	   ( (i == 0) || box_do_all_x_maj_ticks(psr) ) ) {
	/* Longer tick marks */
	ps_line (out,xx,yy00,xx,yy00-psr->X.majorticklen);
	/* Plot label on lowest graph only */
	if (i == 0)
	  ps_ctext(out,xx,
		   yy00-DDD-psr->X.majorticklen-psr->X.tickfontsize*0.8,
		   xtick[x],eXCenter);
      } else if ( bRmod(mat[i].axis_x[x] - psr->X.offset, psr->X.minor) &&
		( (i == 0) || box_do_all_x_min_ticks(psr) ) ){
	/* Shorter tick marks */
	ps_line(out,xx,yy00,xx,yy00-psr->X.minorticklen);
      } else if ( bRmod(mat[i].axis_x[x] - psr->X.offset, psr->X.major) ) {
	/* Even shorter marks, only each X.major */
	ps_line(out,xx,yy00,xx,yy00-(psr->boxspacing/2));
      }
    }
    ps_strfont(out,psr->Y.tickfont,psr->Y.tickfontsize);
    /* Make labels for Y axis */
    if (mat[i].axis_y==NULL) {
      snew(mat[i].axis_y,mat[i].ny);
      for(j=0; (j<mat[i].ny); j++)
	mat[i].axis_y[j]=j;
    }
    snew(ytick,mat[i].ny);
    for(j=0; (j<mat[i].ny); j++) {
      sprintf(buf,"%g",mat[i].axis_y[j]);
      ytick[j]=strdup(buf);
    }

    for(y=0; (y<mat[i].ny); y++) {
      yy=yy00+(y+0.7)*psr->yboxsize;
      if ( bRmod(mat[i].axis_y[y] - psr->Y.offset, psr->Y.major) || 
	   (psr->Y.first && (y==0))) {
	/* Major ticks */
	strlength=max(strlength,(int)strlen(ytick[y]));
	ps_line (out,xx00,yy,xx00-psr->Y.majorticklen,yy);
	ps_ctext(out,xx00-psr->Y.majorticklen-DDD,
		 yy-psr->Y.tickfontsize/3.0,ytick[y],eXRight);
      }
      else if ( bRmod(mat[i].axis_y[y] - psr->Y.offset, psr->Y.minor) ) {
	/* Minor ticks */
	ps_line(out,xx00,yy,xx00-psr->Y.minorticklen,yy);
      }
    }
    sfree(xtick);
    sfree(ytick);
    
    /* Label on Y-axis */
    ps_strfont(out,psr->Y.font,psr->Y.fontsize);
    ps_rotate(out,TRUE);
    xxx=x0-psr->X.majorticklen-psr->X.tickfontsize*strlength-DDD;
    ps_ctext(out,yy00+box_height(&mat[i],psr)/2.0,612.5-xxx,
	     mat[i].label_y,eXCenter);
    ps_rotate(out,FALSE);

    yy00+=box_height(&(mat[i]),psr)+box_dh(psr)+box_dh_top(psr);
  }
  /* Label on X-axis */  
  ps_strfont(out,psr->X.font,psr->X.fontsize);
  if (strlen(mat[0].label_x) == 0) 
    ps_ctext(out,x0+w/2,y0-DDD-psr->X.majorticklen-psr->X.tickfontsize*FUDGE-
	     psr->X.fontsize,psr->X.label,eXCenter);
  else
    ps_ctext(out,x0+w/2,y0-DDD-psr->X.majorticklen-psr->X.tickfontsize*FUDGE-
	     psr->X.fontsize,mat[0].label_x,eXCenter);
}

static void draw_zerolines(FILE *out,real x0,real y0,real w,
			   int nmat,t_matrix mat[],t_psrec *psr)
{
  real   xx,yy,dy,xx00,yy00;
  int    i,x,y;
  
  xx00=x0-1.5;
  yy00=y0-1.5;
  ps_linewidth(out,psr->zerolinewidth);
  for (i=0; (i<nmat); i++) {
    dy=box_height(&(mat[i]),psr);
    /* mat[i].axis_x and _y were already set by draw_boxes */
    if (psr->X.lineatzero) {
      ps_rgb(out,linecolors[psr->X.lineatzero]);
      for(x=0; (x<mat[i].nx); x++) {
	xx=xx00+(x+0.7)*psr->xboxsize;
      /* draw lines whenever tick label almost zero (e.g. next trajectory) */
	if ( x!=0 && x<mat[i].nx-1 &&
	     abs(mat[i].axis_x[x]) < 
	     0.1*abs(mat[i].axis_x[x+1]-mat[i].axis_x[x]) ) {
	  ps_line (out,xx,yy00,xx,yy00+dy+2);
	}
      }
    }
    if (psr->Y.lineatzero) {
      ps_rgb(out,linecolors[psr->Y.lineatzero]);
      for(y=0; (y<mat[i].ny); y++) {
	yy=yy00+(y+0.7)*psr->yboxsize;
	/* draw lines whenever tick label almost zero (e.g. next trajectory) */
	if ( y!=0 && y<mat[i].ny-1 && 
	     abs(mat[i].axis_y[y]) < 
	     0.1*abs(mat[i].axis_y[y+1]-mat[i].axis_y[y]) ) {
	  ps_line (out,xx00,yy,xx00+w+2,yy);
	}
      }
    }
    yy00+=box_height(&(mat[i]),psr)+box_dh(psr)+box_dh_top(psr);
  }
}

static void box_dim(int nmat,t_matrix mat[],t_matrix *mat2,t_psrec *psr,
		    char w_legend,bool bFrame,
		    real *w,real *h,real *dw,real *dh)
{
  int i,maxytick;
  real ww,hh,dww,dhh;
  
  hh=dww=dhh=0;
  maxytick=0;
  
  ww=0;
  for(i=0; (i<nmat); i++) {
    ww=max(ww,mat[i].nx*psr->xboxsize);
    hh+=box_height(&(mat[i]),psr);
    maxytick=max(maxytick,mat[i].nx);
  }
  if (bFrame) {
    if (mat[0].label_y[0])
      dww+=2.0*(psr->Y.fontsize+DDD);
    if (psr->Y.major > 0) 
      dww += psr->Y.majorticklen + DDD + 
	psr->Y.tickfontsize*(log(maxytick)/log(10.0));
    else if (psr->Y.minor > 0)
      dww+=psr->Y.minorticklen;
    
    if (mat[0].label_x[0])
      dhh+=psr->X.fontsize+2*DDD;
    if (((w_legend == 'b') && (mat[0].legend[0] || mat2[0].legend[0])) ||
	((w_legend == 'f') && mat[0].legend[0]) ||
    ((w_legend == 's') && mat2[0].legend[0]))
      dhh+=2*(psr->legfontsize*FUDGE+2*DDD);
    else 
      dhh+=psr->legfontsize*FUDGE+2*DDD;
    if (psr->X.major > 0)
    dhh+=psr->X.tickfontsize*FUDGE+2*DDD+psr->X.majorticklen;
    else if (psr->X.minor > 0)
      dhh+=psr->X.minorticklen;
    
    hh+=(nmat-1)*box_dh(psr)+nmat*box_dh_top(psr);
  }
  *w=ww;
  *h=hh;
  *dw=dww;
  *dh=dhh;
}

void xpm_mat(char *outf,
	     int nmat,t_matrix *mat,t_matrix *mat2,bool bDiag,bool bFirstDiag)
{
  FILE   *out;
  char   buf[100];
  int    i,j,k,x,y,col;
  static char mapper[] = "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789!@#$%^&*()-_=+{}|;:',<.>/?"; 
#define NMAP strlen(mapper)
  int       nmap,nmap1,nmap2;
  t_mapping *map;

   out=ffopen(outf,"w");
   
   for(i=0; i<nmat; i++) {
     if (!mat2 || !diff_maps(mat[i].nmap,mat[i].map,mat2[i].nmap,mat2[i].map))
       write_xpm_m(out,mat[0]);
     else {
       nmap1=mat[i].nmap;
       nmap2=mat2[i].nmap;
       nmap=nmap1+nmap2;
       snew(map,nmap);
       if (nmap > NMAP*NMAP) 
	 fatal_error(0,"Not enough symbols to merge the two colormaps\n");
       for(j=0; j<nmap1; j++) {
	 map[j].code.c1=mapper[j % NMAP];
	 if (nmap > NMAP)
	   map[j].code.c2=mapper[j/NMAP];
	 map[j].rgb.r=mat[i].map[j].rgb.r;
	 map[j].rgb.g=mat[i].map[j].rgb.g;
	 map[j].rgb.b=mat[i].map[j].rgb.b;
	 map[j].desc=mat[i].map[j].desc;
       }
       for(j=0; j<nmap2; j++) {
	 k=j+nmap1;
	 map[k].code.c1=mapper[k % NMAP];
	 if (nmap > NMAP)
	   map[k].code.c2=mapper[k/NMAP];
	 map[k].rgb.r=mat2[i].map[j].rgb.r;
	 map[k].rgb.g=mat2[i].map[j].rgb.g;
	 map[k].rgb.b=mat2[i].map[j].rgb.b;
	 map[k].desc=mat2[i].map[j].desc;
       }
       for(x=0; (x<mat[i].nx); x++) {
	 for(y=0; (y<mat[i].nx); y++) {
	   if ((x<y) || ((x==y) && bFirstDiag)) /* upper left  -> map1 */
	     col=mat[i].matrix[x][y];
	   else /* lower right -> map2 */
	     col=nmap1+mat[i].matrix[x][y];
	   if ((bDiag) || (x!=y))
	     mat[i].matrix[x][y]=col;
	   else
	     mat[i].matrix[x][y]=0;
	 }
       }
       mat[i].nmap=nmap;
       mat[i].map=map;
       if (mat2 && (strcmp(mat[i].title,mat2[i].title) != 0))
	 sprintf(mat[i].title,"%s / %s",mat[i].title,mat2[i].title);
       if (mat2 && (strcmp(mat[i].legend,mat2[i].legend) != 0))
	 sprintf(mat[i].legend,"%s / %s",mat[i].legend,mat2[i].legend); 
       write_xpm_m(out,mat[i]);
       }
   }
   fclose(out);
}

void ps_mat(char *outf,int nmat,t_matrix mat[],t_matrix mat2[],
	    bool bFrame,
	    bool bDiag,bool bFirstDiag,bool bTitle,char w_legend,
	    real boxx,real boxy,char *m2p,char *m2pout)
{
  char   buf[256],*legend;
  FILE   *out;
  t_psrec  psrec,*psr;
  int    W,H;
  int    i,x,y,col,leg=0;
  real   x0,y0,xx;
  real   w,h,dw,dh;
  int       nmap1=0,nmap2=0,leg_nmap;
  t_mapping *map1=NULL,*map2=NULL,*leg_map;
  bool   bMap1,bNextMap1,bDiscrete;
  
  get_params(libfn(m2p),m2pout,&psrec);

  psr=&psrec;

  if (boxx>0) {
    psr->xboxsize=boxx;
    psr->yboxsize=boxx;
  }
  if (boxy>0)
    psr->yboxsize=boxy;  

  nmap1=0;
  for (i=0; (i<nmat); i++)
    if (mat[i].nmap>nmap1) {
      nmap1=mat[i].nmap;
      map1=mat[i].map;
      leg=i+1;
    }
  if (leg!=1)
    printf("Selected legend of matrix # %d for display\n",leg);
  if (mat2) {
    nmap2=0;
    for (i=0; (i<nmat); i++)
      if (mat2[i].nmap>nmap2) {
	nmap2=mat2[i].nmap;
	map2=mat2[i].map;
	leg=i+1;
  }
    if (leg!=1)
      printf("Selected legend of matrix # %d for second display\n",leg);
  }
  if ( (mat[0].legend[0]==0) && psr->legend )
    strcpy(mat[0].legend, psr->leglabel);

  bTitle = bTitle && mat[nmat-1].title[0];
  psr->bTitle = bTitle;

  /* Set up size of box for nice colors */
  box_dim(nmat,mat,mat2,psr,w_legend,bFrame,&w,&h,&dw,&dh);
  
  /* Set up bounding box */
  W=w+dw;
  H=h+dh;
  
  /* Start box at */
  x0=dw;
  y0=dh;
  x = W+psr->xoffs;
  y = H+psr->yoffs;
  if (bFrame) {
    x += 5*DDD;
    y += 4*DDD;
  }
  out=ps_open(outf,0,0,x,y);
  ps_linewidth(out,psr->linewidth);
  ps_init_rgb_box(out,psr->xboxsize,psr->yboxsize);
  ps_init_rgb_nbox(out,psr->xboxsize,psr->yboxsize);
  ps_translate(out,psr->xoffs,psr->yoffs);

  if (bFrame) {
    ps_comment(out,"Here starts the BOX drawing");  
    draw_boxes(out,x0,y0,w,nmat,mat,psr);
  }

  for(i=0; (i<nmat); i++) {
    if (bTitle) {
      /* Print title, if any */
      ps_rgb(out,BLACK);
      ps_strfont(out,psr->titfont,psr->titfontsize); 
      if (!mat2 || (strcmp(mat[i].title,mat2[i].title) == 0))
	strcpy(buf,mat[i].title);
      else
	sprintf(buf,"%s / %s",mat[i].title,mat2[i].title);
      ps_ctext(out,x0+w/2,y0+box_height(&(mat[i]),psr)+psr->titfontsize,
	       buf,eXCenter);
    }
    sprintf(buf,"Here starts the filling of box #%d",i);
    ps_comment(out,buf);
    for(x=0; (x<mat[i].nx); x++) {
      int nexty;
      int nextcol;
      
      xx=x0+x*psr->xboxsize;
      ps_moveto(out,xx,y0);
      y=0;
      bMap1 = (!mat2 || (x<y || (x==y && bFirstDiag)));
      if ((bDiag) || (x!=y))
	col = mat[i].matrix[x][y];
      else
	col = -1;
      for(nexty=1; (nexty<=mat[i].ny); nexty++) {
	bNextMap1 = (!mat2 || (x<nexty || (x==nexty && bFirstDiag)));
	  /* TRUE:  upper left  -> map1 */
	  /* FALSE: lower right -> map2 */
	if ((nexty==mat[i].ny) || (!bDiag && (x==nexty)))
	  nextcol = -1;
	else
	  nextcol=mat[i].matrix[x][nexty];
	if ( (nexty==mat[i].ny) || (col!=nextcol) || (bMap1!=bNextMap1) ) {
	  if (col >= 0)
	    if (bMap1)
	      ps_rgb_nbox(out,&(mat[i].map[col].rgb),nexty-y);
	    else
	      ps_rgb_nbox(out,&(mat2[i].map[col].rgb),nexty-y);
	  else
	    ps_moverel(out,0,psr->yboxsize);
	  y=nexty;
	  bMap1=bNextMap1;
	  col=nextcol;
	  }
	}
    }
    y0+=box_height(&(mat[i]),psr)+box_dh(psr)+box_dh_top(psr);
  }
  
  if (psr->X.lineatzero || psr->Y.lineatzero) {
    /* reset y0 for first box */
    y0=dh;
    ps_comment(out,"Here starts the zero lines drawing");  
    draw_zerolines(out,x0,y0,w,nmat,mat,psr);
  }
  
  if (w_legend != 'n') {
    ps_comment(out,"Now it's legend time!");
    ps_linewidth(out,psr->linewidth);
    if ((mat2==NULL) || (w_legend != 's')) {
      bDiscrete = mat[0].bDiscrete;
      legend    = mat[0].legend;
      leg_nmap  = nmap1;
      leg_map   = map1;
    } else {
      bDiscrete = mat2[0].bDiscrete;
      legend    = mat2[0].legend;
      leg_nmap  = nmap2;
      leg_map   = map2;
    }
    if (bDiscrete)
      leg_discrete(out,psr->legfontsize,DDD,legend,
		   psr->legfontsize,psr->legfont,leg_nmap,leg_map);
    else {
      if (w_legend != 'b')
	leg_continuous(out,x0+w/2,w/2,DDD,legend,
		       psr->legfontsize,psr->legfont,leg_nmap,leg_map);
      else
	leg_bicontinuous(out,x0+w/2,w,DDD,mat[0].legend,mat2[0].legend,
			 psr->legfontsize,psr->legfont,nmap1,map1,nmap2,map2);
    }
    ps_comment(out,"Were there, dude");
  }
  
  ps_close(out);
}

void do_mat(int nmat,t_matrix *mat,t_matrix *mat2,
	    bool bFrame,
	    bool bDiag,bool bFirstDiag,bool bTitle,char w_legend,
	    real boxx,real boxy,
	    char *epsfile,char *xpmfile,char *m2p,char *m2pout)
{
  int      i,j,k,copy_start;

  if (mat2) {
    for (k=0; (k<nmat); k++) {
      if ((mat2[k].nx!=mat[k].nx) || (mat2[k].ny!=mat[k].ny)) 
	fatal_error(0,"WAKE UP!! Size of frame %d in 2nd matrix file (%dx%d) does not match size of 1st matrix (%dx%d) or the other way around.\n",
		    k,mat2[k].nx,mat2[k].ny,mat[k].nx,mat[k].ny);
      for (j=0; (j<mat[k].ny); j++) {
	if (bFirstDiag)
	  copy_start = j+1;
	else
	  copy_start = j;
	for (i=copy_start; (i<mat[k].nx); i++)
	  mat[k].matrix[i][j]=mat2[k].matrix[i][j];
      }
    }
  }
  for(i=0; (i<nmat); i++) 
    fprintf(stderr,"Matrix %d is %d x %d\n",i,mat[i].nx,mat[i].ny);

  
  if (epsfile!=NULL)
    ps_mat(epsfile,nmat,mat,mat2,bFrame,bDiag,bFirstDiag,
	   bTitle,w_legend,boxx,boxy,m2p,m2pout);
  if (xpmfile!=NULL)
    xpm_mat(xpmfile,nmat,mat,mat2,bDiag,bFirstDiag);
}

void rainbow_map(char *rainbow, int nmat, t_matrix mat[])
{
  t_mapping *map;
  int m,i;
  real c,r,g,b;

  for(m=0; m<nmat; m++) {
    map = mat[m].map;
    for(i=0; i<mat[m].nmap; i++) {
      c = (map[i].rgb.r + map[i].rgb.g + map[i].rgb.b)/3;
      if (c > 1)
	  c = 1;
      if (rainbow[0] == 'b')
	c = 1 - c;
      if (c <= 0.25) {
	r = 0;
	g = pow(4*c,0.666);
	b = 1;
      } else if (c <= 0.5) {
	r = 0;
	g = 1;
	b = pow(2-4*c,0.666);
      } else if (c <= 0.75) {
	r = pow(4*c-2,0.666);
	g = 1;
	b = 0;
      } else {
	r = 1;
	g = pow(4-4*c,0.666);
	b = 0;
      }
      map[i].rgb.r = r;
      map[i].rgb.g = g;
      map[i].rgb.b = b;
    }
  }
}

int main(int argc,char *argv[])
{
  static char *desc[] = {
    "xpm2ps makes a beautiful color plot of an XPixelMap file.",
    "Labels and axis can be displayed, when they are supplied",
    "in the correct matrix format.",
    "Matrix data may be generated by programs such as do_dssp, g_rms or",
    "g_mdmat.[PAR]",
    "Parameters are set in the [TT]m2p[tt] file optionally supplied with",
    "[TT]-di[tt]. Reasonable defaults are provided. Settings for the y-axis",
    "default to those for the x-axis. Font names have a defaulting hierarchy:",
    "titlefont -> legendfont; titlefont -> (xfont -> yfont -> ytickfont)",
    "-> xtickfont, e.g. setting titlefont sets all fonts, setting xfont",
    "sets yfont, ytickfont and xtickfont.[PAR]",
    "With [TT]-f2[tt] a 2nd matrix file can be supplied, both matrix",
    "files will be read simultaneously and the upper left half of the",
    "first one ([TT]-f[tt]) is plotted together with the lower right",
    "half of the second one ([TT]-f2[tt]). The diagonal will contain",
    "values from the matrix file selected with [TT]-diag[tt].",
    "Plotting of the diagonal values can be suppressed altogether by",
    "setting [TT]-diag[tt] to [TT]none[tt].[PAR]",
    "If the color coding and legend labels of both matrices are identical,",
    "only one legend will be displayed, else two separate legends are",
    "displayed.[PAR]",
    "[TT]-title[tt] can be set to [TT]none[tt] to suppress the title, or to",
    "[TT]ylabel[tt] to show the title in the Y-label position (alongside",
    "the Y-axis).[PAR]",
    "With the [TT]-rainbow[tt] option dull grey-scale matrices can be turned",
    "into attractive color pictures.[PAR]",
    "Merged or rainbowed matrices can be written to an XPixelMap file with",
    "the [TT]-xpm[tt] option."
  };

  char      *fn,*epsfile=NULL,*xpmfile=NULL;
  char      w_legend;
  int       i,nmat,nmat2;
  t_matrix *mat=NULL,*mat2=NULL;
  bool      bTitle,bDiag,bFirstDiag;
  static bool bFrame=TRUE;
  static real boxx=0,boxy=0;
  static char *title[]   = { NULL, "top", "ylabel", "none", NULL };
  static char *legend[]  = { NULL, "both", "first", "second", "none", NULL };
  static char *diag[]    = { NULL, "first", "second", "none", NULL };
  static char *rainbow[] = { NULL, "no", "blue", "red", NULL };
  t_pargs pa[] = {
    { "-frame",   FALSE, etBOOL, {&bFrame}, "Display frame, ticks, labels, title and legend" },
    { "-title",   FALSE, etENUM, {title},   "Show title at" },
    { "-legend",  FALSE, etENUM, {legend},  "Show legend" },
    { "-diag",    FALSE, etENUM, {diag},    "Diagonal" },
    { "-bx",      FALSE, etREAL, {&boxx},
      "Box x-size (also y-size when -by is not set)" },
    { "-by",      FALSE, etREAL, {&boxy},   "Box y-size" },
    { "-rainbow", FALSE, etENUM, {rainbow}, "Rainbow colors, convert white to" }
  };
  t_filenm  fnm[] = {
    { efXPM, "-f",  NULL,      ffREAD },
    { efXPM, "-f2", "root2",   ffOPTRD },
    { efM2P, "-di", NULL,      ffLIBRD },
    { efM2P, "-do", "out",     ffOPTWR },
    { efEPS, "-o",  NULL,      ffOPTWR },
    { efXPM, "-xpm",NULL,      ffOPTWR }
  };
#define NFILE asize(fnm)
  
  CopyRight(stdout,argv[0]);
  parse_common_args(&argc,argv,PCA_CAN_VIEW,FALSE,
		    NFILE,fnm,asize(pa),pa,
		    asize(desc),desc,0,NULL);

  if (!bFrame) {
    title[0]  = "none";
    legend[0] = "none";
  }

  if (ftp2bSet(efEPS,NFILE,fnm))
    epsfile=ftp2fn(efEPS,NFILE,fnm);
  if (opt2bSet("-xpm",NFILE,fnm))
    xpmfile=opt2fn("-xpm",NFILE,fnm);
  if ((epsfile==NULL) && (xpmfile==NULL))
    epsfile=ftp2fn(efEPS,NFILE,fnm);

  bDiag      = (diag[0][0] != 'n');
  bFirstDiag = (diag[0][0] != 's');
  
  fn=opt2fn("-f",NFILE,fnm);
  nmat=read_xpm_matrix(fn,&mat);
  fprintf(stderr,"There are %d matrices in %s\n",nmat,fn);
  if (opt2bSet("-f2",NFILE,fnm)) {
    fn=opt2fn("-f2",NFILE,fnm);
    nmat2=read_xpm_matrix(fn,&mat2);
    fprintf(stderr,"There are %d matrices in %s\n",nmat2,fn);
    if (nmat != nmat2) {
      fprintf(stderr,"Different number of matrices, using the smallest number.\n");
      nmat=nmat2=min(nmat,nmat2);
    }
  }
  else {
    nmat2=0;
  }
  bTitle = (title[0][0] == 't');
  if (title[0][0] == 'y') {
    bTitle=FALSE; /* don't print title in two places at once */
    for (i=0; (i<nmat); i++) {
      strcpy(mat[i].label_y, mat[i].title);
      strcpy(mat2[i].label_y, mat2[i].title);
    }
  }
  if (rainbow[0][0] != 'n') {
    rainbow_map(rainbow[0],nmat,mat);
    rainbow_map(rainbow[0],nmat,mat);
  }

  w_legend = legend[0][0];
  if ((mat2 == NULL) && (w_legend != 'n'))
    w_legend = 'f';

  do_mat(nmat,mat,mat2,bFrame,bDiag,bFirstDiag,bTitle,w_legend,
	 boxx,boxy,epsfile,xpmfile,
	 opt2fn_null("-di",NFILE,fnm),opt2fn_null("-do",NFILE,fnm));
  
  if (bDoView())
    viewps(ftp2fn(efEPS,NFILE,fnm));
    
  thanx(stderr);
  
  return 0;
}
