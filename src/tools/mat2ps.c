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
 * Giant Rising Ordinary Mutants for A Clerical Setup
 */
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
  real majorticklen;
  real minorticklen;
  char label[STRLEN];
  real fontsize;
  char font[STRLEN];
  real tickfontsize;
  char tickfont[STRLEN];
} t_axisdef;

enum { elOff, elDiscrete, elContinuous, elBicontinuous, elNR };

typedef struct {
  int       bw;
  real      xoffs,yoffs;
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

void get_params(char *mpin,char *mpout,
		t_psrec *psr)
{
  static char *bools[BOOL_NR+1]  = { "no", "yes", NULL };
  static char *legs[] = { "Off", "Discrete", "Continuous", "Bicontinuous", NULL };
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
  ITYPE("x-major",      	psr->X.major,     25);
  ITYPE("x-minor",      	psr->X.minor,     5);
  RTYPE("x-majorticklen", 	psr->X.majorticklen,	8.0);
  RTYPE("x-minorticklen", 	psr->X.minorticklen,	4.0);
  STYPE("x-label",		psr->X.label,		"");
  RTYPE("x-fontsize",     	psr->X.fontsize, 	16.0);
  STYPE("x-font",      		psr->X.font,      	"Times-Roman");
  RTYPE("x-tickfontsize",     	psr->X.tickfontsize, 	10.0);
  STYPE("x-tickfont",      	psr->X.tickfont,		"Helvetica");
  ITYPE("y-major",      	psr->Y.major,      	25);
  ITYPE("y-minor",      	psr->Y.minor,      	5);
  RTYPE("y-majorticklen", 	psr->Y.majorticklen,  	8.0);
  RTYPE("y-minorticklen", 	psr->Y.minorticklen,  	4.0);
  STYPE("y-label",		psr->Y.label,		"");
  RTYPE("y-fontsize",     	psr->Y.fontsize, 	16.0);
  STYPE("y-font",      		psr->Y.font,      	"Times-Roman");
  RTYPE("y-tickfontsize",     	psr->Y.tickfontsize, 	10.0);
  STYPE("y-tickfont",      	psr->Y.tickfont,		"Helvetica");

  write_inpfile(mpout,ninp,inp);
}

t_rgb black={ 0.0, 0.0, 0.0 };
#define BLACK (&black)
  
void leg_discrete(FILE *ps,int x0,int y0,char *label,
		  real fontsize,char *font,int nmap,t_mapping map[])
{
  int   i,yy;
  real  boxhh;
  
  boxhh=fontsize+DDD;
  ps_strfont(ps,font,fontsize);
  /* LANDSCAPE */
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
    y+=dy+box_dh(psr);
  }
  
  /* Draw the ticks on the axes */
  xx00=x0-1;
  yy00=y0-1;
  for (i=0; (i<nmat); i++) {
    /* Make labels for x axis */
    snew(xtick,mat[i].nx);
    for(j=0; (j<mat[i].nx); j++) {
      sprintf(buf,"%.0f",mat[i].axis_x[j]);
      xtick[j]=strdup(buf);
    }
    ps_strfont(out,psr->X.tickfont,psr->X.tickfontsize);
    for(x=0; (x<mat[i].nx); x++) {
      xx=xx00+(x+0.7)*psr->xboxsize;
      itx=mat[i].axis_x[0];
      if ((x == 0) || (((x+itx) % psr->X.major) == 0)) {
	/* Longer tick marks */
	ps_line (out,xx,yy00,xx,yy00-psr->X.majorticklen);
	/* Plot label on lowest graph only */
	if (i == 0)
	  ps_ctext(out,xx,
		   yy00-DDD-psr->X.majorticklen-psr->X.tickfontsize*0.8,
		   xtick[x],eXCenter);
      }
      else if (((x+itx) % psr->X.minor) == 0) {
	/* Shorter tick marks */
	ps_line(out,xx,yy00,xx,yy00-psr->X.minorticklen);
      }
    }
    ps_strfont(out,psr->Y.tickfont,psr->Y.tickfontsize);
    /* Make labels for Y axis */
    snew(ytick,mat[i].ny);
    for(j=0; (j<mat[i].ny); j++) {
      if (mat[0].axis_y==NULL) 
	sprintf(buf,"%d",j+mat[i].y0);
      else 
	sprintf(buf,"%.0f",mat[i].axis_y[j]);
      ytick[j]=strdup(buf);
    }
    
    strlength=max(strlength,(int)strlen(ytick[mat[i].ny-1]));
    for(y=0; (y<mat[i].ny); y++) {
      yy=yy00+(y+0.7)*psr->yboxsize;
      if ((y == 0) || (((y+mat[i].y0) % psr->Y.major) == 0)) {
	/* Major ticks */
	ps_line (out,xx00,yy,xx00-psr->Y.majorticklen,yy);
	ps_ctext(out,xx00-psr->Y.majorticklen-DDD,
		 yy-psr->Y.tickfontsize/3.0,ytick[y],eXRight);
      }
      else if (((y+mat[i].y0) % psr->Y.minor) == 0) {
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

    yy00+=box_height(&(mat[i]),psr)+box_dh(psr);
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
  if ((int)strlen(psr->Y.label) > 0)
    dww+=2*(psr->Y.fontsize+DDD);
  if (psr->Y.major > 0) 
    dww+=psr->Y.majorticklen+DDD+psr->Y.tickfontsize*(log(maxytick)/log(10.0));
  else if (psr->Y.minor > 0)
    dww+=psr->Y.minorticklen;
  
  dhh=0;
  if ((int)strlen(psr->X.label) > 0)
    dhh+=psr->X.fontsize+2*DDD;
  switch (psr->legtype) {
  case elDiscrete:
    dhh+=psr->legfontsize*FUDGE+2*DDD;
    break;
  case elContinuous:
  case elBicontinuous:
    dhh+=2*(psr->legfontsize*FUDGE+2*DDD);
    break;
  default:
    break;
  }
  if (psr->X.major > 0)
    dhh+=psr->X.tickfontsize*FUDGE+2*DDD+psr->X.majorticklen;
  else if (psr->X.minor > 0)
    dhh+=psr->X.minorticklen;
    
  hh+=(nmat-1)*box_dh(psr);
    
  *w=ww;
  *h=hh;
  *dw=dww;
  *dh=dhh;
}

void xpm_mat(char *outf,
	     int nmat,t_matrix mat[],
	     int nmap1,t_mapping map1[],
	     int nmap2,t_mapping map2[],bool bDiag)
{
  FILE   *out;
  int    i,j,x,y,col,offset;
  static char mapper[] = "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789!@#$%^&*()-_=+{}|;:',<.>/?"; 

  i=0;
  if (map1 == map2) offset=0; else offset=nmap1;
  out=ffopen(outf,"w");
  fprintf(out,"/* XPM */\n");
  fprintf(out,"static char * gv_xpm[] = {\n");
  fprintf(out,"\"%d %d   %d 1\",\n",mat[i].nx,mat[i].ny,offset+nmap2);
  for(j=0;j<nmap1;j++) {
    fprintf(out,"\"%c  c #%02X%02X%02X \",\n",mapper[j],
	    (int)(255*map1[j].rgb.r),
            (int)(255*map1[j].rgb.g),
            (int)(255*map1[j].rgb.b));
  }
  if (map1 != map2) {
    for(j=0;j<nmap2;j++) {
      fprintf(out,"\"%c  c #%02X%02X%02X \",\n",mapper[offset+j],
	    (int)(255*map2[j].rgb.r),
            (int)(255*map2[j].rgb.g),
            (int)(255*map2[j].rgb.b));
    }
  }
  for(y=mat[i].ny-1;(y>=0); y--) {
    fprintf(out,"\"");
    for(x=0; (x<mat[i].nx); x++) {
      if (x<=y) { /* upper left  -> map1 */
	col=searchcmap(nmap1,map1,mat[i].matrix[x][y]);
      } else  {   /* lower right -> map2 */
	col=offset+searchcmap(nmap2,map2,mat[i].matrix[x][y]);
      }
      if ((bDiag) || (x!=y))
	fprintf(out,"%c",mapper[col]);
      else
	fprintf(out,"%c",mapper[0]); 
    }
    if (y>0) 
      fprintf(out,"\",\n");
    else
      fprintf(out,"\"\n}\n");

  }
  fclose(out);
}

void ps_mat(char *outf,t_psrec *psr,
	     int nmat,t_matrix mat[],
	     int nmap1,t_mapping map1[],
	     int nmap2,t_mapping map2[],bool bDiag)
{
  char   buf[256];
  FILE   *out;
  int    W,H;
  int    i,x,y,x0,y0,xx,col;
  real   w,h,dw,dh;
  
  /* Set up size of box for nice colors */
  box_dim(nmat,mat,psr,&w,&h,&dw,&dh);
  
  /* Set up bounding box */
  W=w+dw;
  H=h+dh;
  
  /* Start box at */
  x0=dw;
  y0=dh;
  out=ps_open(outf,0,0,W+psr->xoffs+5*DDD,H+psr->yoffs+4*DDD);
  ps_translate(out,psr->xoffs,psr->yoffs);

  ps_comment(out,"Here starts the BOX drawing");  
  draw_boxes(out,x0,y0,w,h,nmat,mat,psr);

  /* LANDSCAPE */
  for(i=0; (i<nmat); i++) {
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
    y0+=box_height(&(mat[i]),psr)+box_dh(psr);
  }
  
  ps_comment(out,"Now it's legend time!");
  switch (psr->legtype) {
  case elDiscrete:
    leg_discrete(out,psr->legfontsize,DDD,psr->leglabel,
		 psr->legfontsize,psr->legfont,nmap1,map1);
    break;
  case elContinuous:
    leg_continuous(out,x0+w/2,w,DDD,psr->leglabel,
		   psr->legfontsize,psr->legfont,nmap1,map1);
    break;
  case elBicontinuous:
    leg_bicontinuous(out,x0+w/2,w,DDD,psr->leglabel,psr->leg2label,
		     psr->legfontsize,psr->legfont,nmap1,map1,nmap2,map2);
    break;
  default:
    fprintf(stderr,"No legend\n");
  }
  ps_comment(out,"Were there, dude");
  
  ps_close(out);
}

void do_mat(char *fn,char *fn2,char *epsfile,char *xpmfile,t_psrec *psr,
	    int n1,t_mapping map1[],
	    int n2,t_mapping map2[],bool bDiag)
{
  int      i,j,k,nmat,nmat2;
  t_matrix *mat,*mat2;

  mat=read_matrix(fn,&nmat);
  fprintf(stderr,"There are %d matrices in %s\n",nmat,fn);
  if (fn2!=NULL) {
    mat2=read_matrix(fn2,&nmat2);
    fprintf(stderr,"There are %d matrices in %s\n",nmat2,fn2);
    if (nmat2!=nmat) {
      if (nmat>nmat2) 
	nmat=nmat2;
      fprintf(stderr,"Will use only first %d frames.\n",nmat);
    }
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
    ps_mat(epsfile,psr,nmat,mat,n1,map1,n2,map2,bDiag);
  if (xpmfile!=NULL)
    xpm_mat(xpmfile,nmat,mat,n1,map1,n2,map2,bDiag);
  
}

int main(int argc,char *argv[])
{
  static char *desc[] = {
    "mat2ps makes a beatiful color plot of a matrix. The matrix has",
    "to be formatted in a special way, it must have been discretized",
    "into a given number of levels that corresponds to the number of",
    "levels in the colormap file.[PAR]",
    "Matrix data may be generated by programs such as do_dssp, sas2mat",
    "g_mdmat or g_rms.[PAR]",
    "With [BB]-f2[bb] a 2nd matrix file can be supplied, both matrix",
    "files will be read simultaneously and the upper left half of the",
    "first one ([BB]-f[bb]) is plotted together with the lower right",
    "half of the second one ([BB]-f2[bb]). The diagonal will contain",
    "values from the first matrix file ([BB]-f[bb]).[PAR]",
    "If legend type 'Bicontinuous' is specified in the [BB]m2p[bb] file",
    "two consecutive colour mapping portions should be present in the",
    "[BB]map[bb] file, one for each legend. This is especially usefull",
    "together with supplying two separate matrix files."
  };
  static char *bugs[] = {
    "Output files can be huge"
  };
  
  char      *epsfile=NULL,*xpmfile=NULL;
  int       n,n1,n2,i;
  t_mapping *map,*map1,*map2;
  t_psrec   psr;
  static bool bDiag=TRUE;
  t_pargs pa[] = {
    { "-diag", FALSE, etBOOL, &bDiag,
      "plot diagonal" }
  };
  t_filenm  fnm[] = {
    { efM2P, "-di", NULL,      ffOPTRD },
    { efM2P, "-do", "m2pout",  ffWRITE },
    { efMAP, "-ci", NULL,      ffLIBRD },
    { efMAP, "-ci2","ss2",     ffOPTRD },
    { efMAP, "-co", "cmapout", ffWRITE },
    { efMAT, "-f",  NULL,      ffREAD },
    { efMAT, "-f2", "ss2",     ffOPTRD },
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
  if (ftp2bSet(efXPM,NFILE,fnm))
    xpmfile=ftp2fn(efXPM,NFILE,fnm);
  if ((epsfile==NULL) && (xpmfile==NULL))
    epsfile=ftp2fn(efEPS,NFILE,fnm);

  if (epsfile != NULL) {
    if (!fexist(opt2fn("-di",NFILE,fnm)))
      fatal_error(0,"%s does not exist, need an m2p file for EPS\n",
		  opt2fn("-di",NFILE,fnm));
    get_params(opt2fn("-di",NFILE,fnm),opt2fn("-do",NFILE,fnm),&psr);
  }

  n1=readcmap(opt2fn("-ci",NFILE,fnm),&map1);

  /*
  n2=(n1-1)/2;
  for(i=0;i<n1;i++) {
    if (i<n2) {
      map1[i].rgb.r=0.0;
      map1[i].rgb.g=sqrt((float)i/n2);
      map1[i].rgb.b=sqrt(1.0-(float)i/n2);
    }
    else {
      map1[i].rgb.r=sqrt((float)(i-n2)/n2);
      map1[i].rgb.g=sqrt(1.0-(float)(i-n2)/n2);
      map1[i].rgb.b=0.0;
    }
  }
  */

  if (opt2bSet("-ci2",NFILE,fnm)) {
    if (opt2bSet("-di",NFILE,fnm))
      if (psr.legtype != elBicontinuous)
	fprintf(stderr,"WARNING!! Will use two colour maps (%s and %s) but no bicontinuous legend!\n",opt2fn("-ci",NFILE,fnm),opt2fn("-ci2",NFILE,fnm));
    if (!opt2bSet("-f2",NFILE,fnm))
      fprintf(stderr,"WARNING!! Two colour maps (%s and %s) specified but only one matrix!\n",opt2fn("-ci",NFILE,fnm),opt2fn("-ci2",NFILE,fnm));
    n2=readcmap(opt2fn("-ci2",NFILE,fnm),&map2);
  } else {
    if (opt2bSet("-di",NFILE,fnm))
      if (psr.legtype == elBicontinuous)
	fprintf(stderr,"WARNING!! Will make a bicontinuous legend with only one colour map (%s)!\n",opt2fn("-ci",NFILE,fnm));
    n2=n1;
    map2=map1;
    writecmap(opt2fn("-co",NFILE,fnm),n2,map2);
  }

  if (opt2bSet("-f2",NFILE,fnm)) {
    do_mat(opt2fn("-f",NFILE,fnm),opt2fn("-f2",NFILE,fnm),
	   epsfile,xpmfile,&psr,n1,map1,n2,map2,bDiag);
  } else {
    do_mat(opt2fn("-f",NFILE,fnm),NULL,
	   epsfile,xpmfile,&psr,n1,map1,n2,map2,bDiag);
  }
  
  if (bDoView())
    viewps(ftp2fn(efEPS,NFILE,fnm));
    
  thanx(stdout);
  
  return 0;
}
