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
#include "readcmap.h"
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

enum { elOff, elDiscrete, elContinuous, elNR };
static char *legs[] = { "Off", "Discrete", "Continuous", NULL };

typedef struct {
  int       bw;
  real      xoffs,yoffs;
  bool      bTitle;
  real      titfontsize;
  char      titfont[STRLEN];
  int       legtype;
  real      legfontsize;
  char      legfont[STRLEN];
  char      leglabel[STRLEN];
  char      leg2label[STRLEN];
  real      xboxsize;
  real      yboxsize;
  t_axisdef X,Y;
} t_psrec;

void get_params(char *mpin,t_psrec *psr)
{
  static char *bools[BOOL_NR+1]  = { "no", "yes", NULL };
  t_inpfile *inp;
  char      *tmp;
  int       ninp;
  
  inp=read_inpfile(mpin,&ninp);
  ETYPE("black&white",		psr->bw,          bools);
  STYPE("titlefont",		psr->titfont,     "Times-Roman");
  RTYPE("titlefontsize",	psr->titfontsize, 20.0);
  ETYPE("legend",		psr->legtype,     legs);
  STYPE("legendfont",		psr->legfont,     "Times-Roman");
  STYPE("legendlabel",		psr->leglabel,	 "");
  STYPE("legend2label",		psr->leg2label,	 psr->leglabel);
  RTYPE("legendfontsize",	psr->legfontsize, 14.0);
  RTYPE("xbox",         	psr->xboxsize,    2.0);
  RTYPE("ybox",         	psr->yboxsize,    2.0);
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
	if (map1[i].code != map2[i].code) bDiff=TRUE;
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
  
void leg_discrete(FILE *ps,int x0,int y0,char *label,
		  real fontsize,char *font,int nmap,t_mapping map[])
{
  int   i,yy,yhh;
  real  boxhh;
  
  boxhh=fontsize+DDD;
  /* LANDSCAPE */
  ps_rgb(ps,BLACK);
  ps_strfont(ps,font,fontsize);
  yhh=y0+fontsize+3*DDD;
  if ((int)strlen(label) > 0)
    ps_ctext(ps,x0,yhh,label,eXLeft);
  yy=y0;
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

void leg_continuous(FILE *ps,int x0,int x,int y0,char *label,
		    real fontsize,char *font,int nmap,t_mapping map[])
{
  int   i,xx0;
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

void leg_bicontinuous(FILE *ps,int x0,int x,int y0,char *label1,char *label2,
		      real fontsize,char *font,
		      int nmap1,t_mapping map1[],int nmap2,t_mapping map2[])
{
  int xx1,xx2,x1,x2;
  x1=x/(nmap1+nmap2)*nmap1;/* width of legend 1 */
  x2=x/(nmap1+nmap2)*nmap2;/* width of legend 2 */
  xx1=x0-(x2/2.0)-fontsize;/* center of legend 1 */
  xx2=x0+(x1/2.0)+fontsize;/* center of legend 2 */
  x1-=fontsize/2;/* adjust width */
  x2-=fontsize/2;/* adjust width */
  leg_continuous(ps,xx1,x1,y0,label1,fontsize,font,nmap1,map1);
  leg_continuous(ps,xx2,x2,y0,label2,fontsize,font,nmap2,map2);
}

static int box_height(t_matrix *mat,t_psrec *psr)
{
  return mat->ny*psr->yboxsize; 
}

static int box_dh(t_psrec *psr)
{
  return 3*psr->Y.majorticklen;
}

static int box_dh_top(t_psrec *psr)
{
  int dh;

  if (psr->bTitle)
    dh=2*psr->titfontsize;
  else
    dh=0;

  return  dh;
}

static void draw_boxes(FILE *out,real x0,real y0,real w,real h,
		       int nmat,t_matrix mat[],t_psrec *psr)
{
  char   buf[12];
  int    xxx;
  char   **xtick,**ytick;
  real   xx,yy,dy,xx00,yy00;
  int    i,j,x,itx,y,strlength=0;
  
  /* Draw the box */
  ps_rgb(out,BLACK);
  y=y0;
  for(i=0; (i<nmat); i++) {
    dy=box_height(&(mat[i]),psr);
    ps_box(out,x0-1,y-1,x0+w+1,y+dy+1);
    y+=dy+box_dh(psr)+box_dh_top(psr);
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
      /*
      itx=mat[i].axis_x[0];
      */
      if ((x % psr->X.major == psr->X.offset) || (psr->X.first && (x==0))) {
	/* Longer tick marks */
	ps_line (out,xx,yy00,xx,yy00-psr->X.majorticklen);
	/* Plot label on lowest graph only */
	if (i == 0)
	  ps_ctext(out,xx,
		   yy00-DDD-psr->X.majorticklen-psr->X.tickfontsize*0.8,
		   xtick[x],eXCenter);
      }
      else if ((x-psr->X.offset) % psr->X.minor == 0) {
	/* Shorter tick marks */
	ps_line(out,xx,yy00,xx,yy00-psr->X.minorticklen);
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
  int       nmap1,nmap2;
  t_mapping map[100];

   out=ffopen(outf,"w");
   
   for(i=0; i<nmat; i++) {
     if (!mat2 || !diff_maps(mat[i].nmap,mat[i].map,mat2[i].nmap,mat2[i].map))
       write_xpm_m(out,mat[0]);
     else {
       nmap1=mat[i].nmap;
       nmap2=mat2[i].nmap;
       if (nmap1+nmap2 > sizeof(mapper)) 
	 fatal_error(0,"Not enough symbols to merge the two colormaps\n");
       for(j=0; j<nmap1; j++) {
	 map[j].code=mapper[j];
	 map[j].rgb.r=mat[i].map[j].rgb.r;
	 map[j].rgb.g=mat[i].map[j].rgb.g;
	 map[j].rgb.b=mat[i].map[j].rgb.b;
	 map[j].desc=mat[i].map[j].desc;
       }
       for(j=0; j<nmap2; j++) {
	 k=j+nmap1;
	 map[k].code=mapper[k];
	 map[k].rgb.r=mat2[i].map[j].rgb.r;
	 map[k].rgb.g=mat2[i].map[j].rgb.g;
	 map[k].rgb.b=mat2[i].map[j].rgb.b;
	 map[k].desc=mat2[i].map[j].desc;
       }
       mat[i].nmap=nmap1+nmap2;
       mat[i].map=map;
       for(x=0; (x<mat[i].nx); x++) {
	 for(y=0; (y<mat[i].nx); y++) {
	   if (x<=y) { /* upper left  -> map1 */
	     col=searchcmap(mat[i].nmap,mat[i].map,mat[i].matrix[x][y]);
	   } else  {   /* lower right -> map2 */
	     col=nmap1+searchcmap(mat2[i].nmap,mat2[i].map,mat[i].matrix[x][y]);
	   }
	   if ((bDiag) || (x!=y))
	     mat[i].matrix[x][y]=mapper[col];
	   else
	     mat[i].matrix[x][y]=mapper[0];
	 }
       }
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
	    real boxx,real boxy,char *m2p)
{
  char   buf[256];
  FILE   *out;
  t_psrec  psrec,*psr;
  int    W,H;
  int    i,x,y,x0,y0,xx,col;
  real   w,h,dw,dh;
  int       nmap1,nmap2;
  t_mapping *map1,*map2;
  
  get_params(libfn(m2p),&psrec);

  psr=&psrec;

  if (boxx>0) {
    psr->xboxsize=boxx;
    psr->yboxsize=boxx;
  }
  if (boxy>0)
    psr->yboxsize=boxy;  

  nmap1=mat[0].nmap;
  map1=mat[0].map;
  if (mat2==NULL) {
    nmap2=nmap1;
    map2=map1;
  }
  else {
    nmap2=mat2[0].nmap;
    map2=mat2[0].map;
  }

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
      xx=x0+x*psr->xboxsize;
      for(y=0; (y<mat[i].ny); y++) {
	if (x<=y) { /* upper left  -> map1 */
	  col=searchcmap(nmap1,map1,mat[i].matrix[x][y]);
	  ps_rgb(out,&(map1[col].rgb));
	} else  {   /* lower right -> map2 */
	  col=searchcmap(nmap2,map2,mat[i].matrix[x][y]);
	  ps_rgb(out,&(map2[col].rgb));
	}
	if ((bDiag) || (x!=y)) { /* no diagonal */
	  ps_fillbox(out,xx,y0+(y*psr->yboxsize),xx+psr->xboxsize,
		     y0+(y+1)*psr->yboxsize);
	}
      }
    }
    y0+=box_height(&(mat[i]),psr)+box_dh(psr)+box_dh_top(psr);
  }
  
  if (!bLegSet && (psr->legtype != elOff) || (bLegSet && bLegend)) {
    ps_comment(out,"Now it's legend time!");
    if (mat[0].bDiscrete)
      leg_discrete(out,psr->legfontsize,DDD,mat[0].legend,
		   psr->legfontsize,psr->legfont,nmap1,map1);
    else {
      if (!diff_maps(nmap1,map1,nmap2,map2))
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
	    char *epsfile,char *xpmfile,char *m2p)
{
  int      i,j,k;

  if (mat2!=NULL) {
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
	   boxx,boxy,m2p);
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
    "With [BB]-f2[bb] a 2nd matrix file can be supplied, both matrix",
    "files will be read simultaneously and the upper left half of the",
    "first one ([BB]-f[bb]) is plotted together with the lower right",
    "half of the second one ([BB]-f2[bb]). The diagonal will contain",
    "values from the first matrix file ([BB]-f[bb]).[PAR]"
  };
  static char *bugs[] = {
    "Output files can be huge"
  };

  char      *fn,*epsfile=NULL,*xpmfile=NULL;
  int       nmat,nmat2;
  bool      bLegSet;
  t_matrix *mat=NULL,*mat2=NULL;
  static bool bTitle=TRUE,bDiag=TRUE,bLegend=TRUE;
  static real boxx=0,boxy=0;
  t_pargs pa[] = {
    { "-title", FALSE, etBOOL, &bTitle,"show title" },
    { "-legend", FALSE, etBOOL, &bLegend,"show legend" },
    { "-diag", FALSE, etBOOL, &bDiag,"plot diagonal" },
    { "-bx", FALSE, etREAL, &boxx,
      "box x-size (also y-size when -by is not set)"},
    { "-by", FALSE, etREAL, &boxy,"box y-size"}
  };
  t_filenm  fnm[] = {
    { efXPM, "-f",  NULL,      ffREAD },
    { efXPM, "-f2", "root2",   ffOPTRD },
    { efM2P, "-di", NULL,      ffLIBRD },
    { efEPS, "-o",  NULL,      ffOPTWR },
    { efXPM, "-xpm",NULL,      ffOPTWR }
  };
#define NFILE asize(fnm)
  
  CopyRight(stdout,argv[0]);
  parse_common_args(&argc,argv,PCA_CAN_VIEW,FALSE,
		    NFILE,fnm,asize(pa),pa,
		    asize(desc),desc,asize(bugs),bugs);

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

  do_mat(nmat,mat,nmat2,mat2,bDiag,bTitle,bLegend,bLegSet,
	 boxx,boxy,epsfile,xpmfile,
	 opt2fn_null("-di",NFILE,fnm));
  
  if (bDoView())
    viewps(ftp2fn(efEPS,NFILE,fnm));
    
  thanx(stdout);
  
  return 0;
}
