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
 * Gyas ROwers Mature At Cryogenic Speed
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
  int  major;
  int  minor;
  int  offset;
  bool first;
  real majorticklen;
  real minorticklen;
  char label[STRLEN];
  real fontsize;
  char font[STRLEN];
  real tickfontsize;
  char tickfont[STRLEN];
} t_axisdef;

typedef struct {
  int       bw;
  real      xoffs,yoffs;
  bool      bTitle;
  real      titfontsize;
  char      titfont[STRLEN];
  bool      legend;
  real      legfontsize;
  char      legfont[STRLEN];
  char      leglabel[STRLEN];
  char      leg2label[STRLEN];
  real      xboxsize;
  real      yboxsize;
  real      boxspacing;
  t_axisdef X,Y;
} t_psrec;

void get_params(char *mpin,char *mpout,t_psrec *psr)
{
  static char *bools[BOOL_NR+1]  = { "no", "yes", NULL };
  t_inpfile *inp;
  char      *tmp;
  int       ninp;
  
  inp=read_inpfile(mpin,&ninp);
  ETYPE("black&white",		psr->bw,          bools);
  STYPE("titlefont",		psr->titfont,     "Times-Roman");
  RTYPE("titlefontsize",	psr->titfontsize, 20.0);
  ETYPE("legend",		psr->legend,     bools);
  STYPE("legendfont",		psr->legfont,     "Times-Roman");
  STYPE("legendlabel",		psr->leglabel,	 "");
  STYPE("legend2label",		psr->leg2label,	 psr->leglabel);
  RTYPE("legendfontsize",	psr->legfontsize, 14.0);
  RTYPE("xbox",         	psr->xboxsize,    2.0);
  RTYPE("ybox",         	psr->yboxsize,    2.0);
  RTYPE("matrixspacing",        psr->boxspacing,  20.0);
  RTYPE("xoffset",              psr->xoffs,       0.0);
  RTYPE("yoffset",              psr->yoffs,       0.0);
  ITYPE("x-major",      	psr->X.major,     20);
  ITYPE("x-minor",      	psr->X.minor,     5);
  ITYPE("x-firstmajor",         psr->X.offset,    0);
  ETYPE("x-majorat0",           psr->X.first,     bools);
  RTYPE("x-majorticklen", 	psr->X.majorticklen,	8.0);
  RTYPE("x-minorticklen", 	psr->X.minorticklen,	4.0);
  STYPE("x-label",		psr->X.label,		"");
  RTYPE("x-fontsize",     	psr->X.fontsize, 	16.0);
  STYPE("x-font",      		psr->X.font,      	"Times-Roman");
  RTYPE("x-tickfontsize",     	psr->X.tickfontsize, 	10.0);
  STYPE("x-tickfont",      	psr->X.tickfont,		"Helvetica");
  ITYPE("y-major",      	psr->Y.major,      	20);
  ITYPE("y-minor",      	psr->Y.minor,      	5);
  ITYPE("y-firstmajor",         psr->Y.offset,    0);
  ETYPE("y-majorat0",           psr->Y.first,     bools);
  RTYPE("y-majorticklen", 	psr->Y.majorticklen,  	8.0);
  RTYPE("y-minorticklen", 	psr->Y.minorticklen,  	4.0);
  STYPE("y-label",		psr->Y.label,		"");
  RTYPE("y-fontsize",     	psr->Y.fontsize, 	16.0);
  STYPE("y-font",      		psr->Y.font,      	"Times-Roman");
  RTYPE("y-tickfontsize",     	psr->Y.tickfontsize, 	10.0);
  STYPE("y-tickfont",      	psr->Y.tickfont,		"Helvetica");
  if (mpout)
    write_inpfile(mpout,ninp,inp);
}

t_rgb black={ 0.0, 0.0, 0.0 };
#define BLACK (&black)

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

static void draw_boxes(FILE *out,real x0,real y0,real w,real h,
		       int nmat,t_matrix mat[],t_psrec *psr)
{
  char   buf[12];
  real   xxx;
  char   **xtick,**ytick;
  real   xx,yy,dy,xx00,yy00;
  int    i,j,x,y,strlength=0;
  
  /* Draw the box */
  ps_rgb(out,BLACK);
  yy00=y0;
  for(i=0; (i<nmat); i++) {
    dy=box_height(&(mat[i]),psr);
    ps_box(out,x0-1,yy00-1,x0+w+1,yy00+dy+1);
    yy00+=dy+box_dh(psr)+box_dh_top(psr);
  }
  
  /* Draw the ticks on the axes */
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
      if ( ( (x % psr->X.major == psr->X.offset) || 
	     (psr->X.first && (x==0))) &&
	   ( (i == 0) || box_do_all_x_maj_ticks(psr) ) ) {
	/* Longer tick marks */
	ps_line (out,xx,yy00,xx,yy00-psr->X.majorticklen);
	/* Plot label on lowest graph only */
	if (i == 0)
	  ps_ctext(out,xx,
		   yy00-DDD-psr->X.majorticklen-psr->X.tickfontsize*0.8,
		   xtick[x],eXCenter);
      } else if ( ( (x-psr->X.offset) % psr->X.minor == 0) &&
		( (i == 0) || box_do_all_x_min_ticks(psr) ) ){
	/* Shorter tick marks */
	ps_line(out,xx,yy00,xx,yy00-psr->X.minorticklen);
      } else if (x % psr->X.major == psr->X.offset) {
	/* Even shorter marks, only each X.offset */
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
    
    strlength=max(strlength,(int)strlen(ytick[mat[i].ny-1]));
    for(y=0; (y<mat[i].ny); y++) {
      yy=yy00+(y+0.7)*psr->yboxsize;
      if ( (y % psr->Y.major == psr->Y.offset) || (psr->Y.first && (y==0))) {
	/* Major ticks */
	ps_line (out,xx00,yy,xx00-psr->Y.majorticklen,yy);
	ps_ctext(out,xx00-psr->Y.majorticklen-DDD,
		 yy-psr->Y.tickfontsize/3.0,ytick[y],eXRight);
      }
      else if ((y-psr->Y.offset) % psr->Y.minor == 0) {
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

static void box_dim(int nmat,t_matrix mat[],t_psrec *psr,
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
  dww=0;
  if (mat[0].label_y[0])
    dww+=2.0*(psr->Y.fontsize+DDD);
  if (psr->Y.major > 0) 
    dww+=psr->Y.majorticklen+DDD+psr->Y.tickfontsize*(log(maxytick)/log(10.0));
  else if (psr->Y.minor > 0)
    dww+=psr->Y.minorticklen;
  
  dhh=0;
  if (mat[0].label_x[0])
    dhh+=psr->X.fontsize+2*DDD;
  if (mat[0].legend[0]==0)
    dhh+=psr->legfontsize*FUDGE+2*DDD;
  else
    dhh+=2*(psr->legfontsize*FUDGE+2*DDD);
  if (psr->X.major > 0)
    dhh+=psr->X.tickfontsize*FUDGE+2*DDD+psr->X.majorticklen;
  else if (psr->X.minor > 0)
    dhh+=psr->X.minorticklen;
    
  hh+=(nmat-1)*box_dh(psr)+nmat*box_dh_top(psr);
  
  *w=ww;
  *h=hh;
  *dw=dww;
  *dh=dhh;
}

void xpm_mat(char *outf,
	     int nmat,t_matrix *mat,t_matrix *mat2,bool bDiag)
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
	   if (x<=y) { /* upper left  -> map1 */
	     col=mat[i].matrix[x][y];
	   } else  {   /* lower right -> map2 */
	     col=nmap1+mat[i].matrix[x][y];
	   }
	   if ((bDiag) || (x!=y))
	     mat[i].matrix[x][y]=col;
	   else
	     mat[i].matrix[x][y]=0;
	 }
       }
       mat[i].nmap=nmap;
       mat[i].map=map;
       if (mat2 && (strcmp(mat[i].title,mat2[i].title) != 0))
	 sprintf(mat[i].title,"%s / %s\0",mat[i].title,mat2[i].title);
       if (mat2 && (strcmp(mat[i].legend,mat2[i].legend) != 0))
	 sprintf(mat[i].legend,"%s / %s\0",mat[i].legend,mat2[i].legend); 
       write_xpm_m(out,mat[i]);
       }
   }
   fclose(out);
}

void ps_mat(char *outf,int nmat,t_matrix mat[],t_matrix mat2[],
	    bool bDiag,bool bTitle,bool bLegend,bool bLegSet,
	    real boxx,real boxy,char *m2p,char *m2pout)
{
  char   buf[256];
  FILE   *out;
  t_psrec  psrec,*psr;
  int    W,H;
  int    i,x,y,col,leg;
  real   x0,y0,xx;
  real   w,h,dw,dh;
  int       nmap1,nmap2;
  t_mapping *map1,*map2;
  bool   bMap1,bNextMap1;
  
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
  box_dim(nmat,mat,psr,&w,&h,&dw,&dh);
  
  /* Set up bounding box */
  W=w+dw;
  H=h+dh;
  
  /* Start box at */
  x0=dw;
  y0=dh;
  /*
  if (bTitle) 
    out=ps_open(outf,0,0,W+psr->xoffs+5*DDD,H+psr->yoffs+4*DDD+
		2*psr->titfontsize);
  else
  */
  out=ps_open(outf,0,0,W+psr->xoffs+5*DDD,H+psr->yoffs+4*DDD);
  ps_init_rgb_box(out,psr->xboxsize,psr->yboxsize);
  ps_init_rgb_nbox(out,psr->xboxsize,psr->yboxsize);
  ps_translate(out,psr->xoffs,psr->yoffs);
    
  ps_comment(out,"Here starts the BOX drawing");  
  draw_boxes(out,x0,y0,w,h,nmat,mat,psr);
  /*
  if (bTitle) {
    ps_strfont(out,psr->titfont,psr->titfontsize); 
    if (!mat2 || (strcmp(mat[nmat-1].title,mat2[nmat-1].title) == 0))
      strcpy(buf,mat[nmat-1].title);
    else
      sprintf(buf,"%s / %s\0",mat[nmat-1].title,mat2[nmat-1].title);
    ps_ctext(out,x0+w/2,y0+h+2*DDD+psr->titfontsize,
	     buf,eXCenter);
  }
  */

  /* LANDSCAPE */
  for(i=0; (i<nmat); i++) {
    if (bTitle) {
      /* Print title, if any */
      ps_rgb(out,BLACK);
      ps_strfont(out,psr->titfont,psr->titfontsize); 
      if (!mat2 || (strcmp(mat[i].title,mat2[i].title) == 0))
	strcpy(buf,mat[i].title);
      else
	sprintf(buf,"%s / %s\0",mat[i].title,mat2[i].title);
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
      bMap1=(!mat2 || (x<=y));
      if ((bDiag) || (x!=y))
	col = mat[i].matrix[x][y];
      else
	col = -1;
      for(nexty=1; (nexty<=mat[i].ny); nexty++) {
	bNextMap1=(!mat2 || (x<=nexty));
	  /* TRUE:  upper left  -> map1 */
	  /* FALSE: lower right -> map2 */
	if ((bDiag) || (x!=nexty))
	  nextcol=mat[i].matrix[x][nexty];
	else
	  nextcol = -1;
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
  
  if ((!bLegSet && psr->legend) || (bLegSet && bLegend)) {
    ps_comment(out,"Now it's legend time!");
    if (mat[0].bDiscrete)
      leg_discrete(out,psr->legfontsize,DDD,mat[0].legend,
		   psr->legfontsize,psr->legfont,nmap1,map1);
    else {
      if ((mat2==NULL) || (!diff_maps(nmap1,map1,nmap2,map2)))
	leg_continuous(out,x0+w/2,w/2,DDD,mat[0].legend,
		       psr->legfontsize,psr->legfont,nmap1,map1);
      else
	leg_bicontinuous(out,x0+w/2,w,DDD,mat[0].legend,mat2[0].legend,
			 psr->legfontsize,psr->legfont,nmap1,map1,nmap2,map2);
    }
    ps_comment(out,"Were there, dude");
  }
  
  ps_close(out);
}

void do_mat(int nmat,t_matrix *mat,int nmat2,t_matrix *mat2,
	    bool bDiag,bool bTitle,bool bLegend,bool bLegSet,
	    real boxx,real boxy,
	    char *epsfile,char *xpmfile,char *m2p,char *m2pout)
{
  int      i,j,k;

  if (mat2) {
    for (k=0; (k<nmat); k++)
      for (j=0; (j<mat[k].ny); j++)
	for (i=j+1; (i<mat[k].nx); i++) {
	  if ((mat2[k].nx!=mat[k].nx) || (mat2[k].ny!=mat[k].ny)) 
	    fatal_error(0,"WAKE UP!! Size of frame %d in 2nd matrix file (%dx%d) does not match size of 1st matrix (%dx%d) or the other way around.\n",
			k,mat2[k].nx,mat2[k].ny,mat[k].nx,mat[k].ny);
	  mat[k].matrix[i][j]=mat2[k].matrix[i][j];
	}
  }
  for(i=0; (i<nmat); i++) 
    fprintf(stderr,"Matrix %d is %d x %d\n",i,mat[i].nx,mat[i].ny);

  
  if (epsfile!=NULL)
    ps_mat(epsfile,nmat,mat,mat2,bDiag,bTitle,bLegend,bLegSet,
	   boxx,boxy,m2p,m2pout);
  if (xpmfile!=NULL)
    xpm_mat(xpmfile,nmat,mat,mat2,bDiag);
}

int main(int argc,char *argv[])
{
  static char *desc[] = {
    "xpm2ps makes a beatiful color plot of an XPixelMap file.",
    "Labels and axis can be displayed, when there are supplied",
    "in the correct matrix format.",
    "Matrix data may be generated by programs such as do_dssp, sas2mat",
    "g_mdmat or g_rms.[PAR]",
    "Parameters are set in the [BB]m2p[bb] file optionally supplied with",
    "[BB]-di[bb]. Reasonable defaults are supplied in a library file.[PAR]",
    "With [BB]-f2[bb] a 2nd matrix file can be supplied, both matrix",
    "files will be read simultaneously and the upper left half of the",
    "first one ([BB]-f[bb]) is plotted together with the lower right",
    "half of the second one ([BB]-f2[bb]). The diagonal will contain",
    "values from the first matrix file ([BB]-f[bb]). Plotting of the",
    "diagonal values can be suppressed altogether using [BB]-nodiag[bb].[PAR]",
    "If the colour",
    "coding and legend labels of both matrices are identical, only",
    "one legend will be output, else two separate legends are output.[PAR]",
    "[BB]-ytitle[bb] shows the title string in the Y-label position",
    "(alongside the Y-axes) and suppresses the title in the normal",
    "position (above the matrix/matrices).[PAR]"
  };

  char      *fn,*epsfile=NULL,*xpmfile=NULL;
  int       i,nmat,nmat2;
  bool      bLegSet;
  t_matrix *mat=NULL,*mat2=NULL;
  static bool bTitle=TRUE,bDiag=TRUE,bLegend=TRUE,bTitLab=FALSE;
  static real boxx=0,boxy=0;
  t_pargs pa[] = {
    { "-title",  FALSE, etBOOL, &bTitle, "show title" },
    { "-ytitle",  FALSE, etBOOL, &bTitLab,"show title in Y-label position"},
    { "-legend", FALSE, etBOOL, &bLegend,"show legend" },
    { "-diag",   FALSE, etBOOL, &bDiag,  "plot diagonal" },
    { "-bx",     FALSE, etREAL, &boxx,
      "box x-size (also y-size when -by is not set)"},
    { "-by",     FALSE, etREAL, &boxy,   "box y-size"}
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

  if (ftp2bSet(efEPS,NFILE,fnm))
    epsfile=ftp2fn(efEPS,NFILE,fnm);
  if (opt2bSet("-xpm",NFILE,fnm))
    xpmfile=opt2fn("-xpm",NFILE,fnm);
  if ((epsfile==NULL) && (xpmfile==NULL))
    epsfile=ftp2fn(efEPS,NFILE,fnm);
  bLegSet=opt2parg_bSet("-legend",asize(pa),pa);
  
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
  if (bTitLab) {
    bTitle=FALSE; /* don't print title on two places at once */
    for (i=0; (i<nmat); i++)
      strcpy(mat[i].label_y, mat[i].title);
    for (i=0; (i<nmat2); i++)
      strcpy(mat2[i].label_y, mat2[i].title);
  }

  do_mat(nmat,mat,nmat2,mat2,bDiag,bTitle,bLegend,bLegSet,
	 boxx,boxy,epsfile,xpmfile,
	 opt2fn_null("-di",NFILE,fnm),opt2fn_null("-do",NFILE,fnm));
  
  if (bDoView())
    viewps(ftp2fn(efEPS,NFILE,fnm));
    
  thanx(stdout);
  
  return 0;
}
