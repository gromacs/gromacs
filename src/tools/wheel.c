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
static char *SRCID_wheel_c = "$Id$";

#include <math.h>
#include "sysstuff.h"
#include "physics.h"
#include "string2.h"
#include "typedefs.h"
#include "smalloc.h"
#include "macros.h"
#include "vec.h"
#include "xvgr.h"
#include "pbc.h"
#include "copyrite.h"
#include "futil.h"
#include "strdb.h"
#include "statutil.h"
#include "pbc.h"
#include "rdgroup.h"
#include "gstat.h"
#include "fatal.h"
#include "writeps.h"
#include "strdb.h"

bool *bPhobics(int nres,char *resnm[])
{
  int  i,nb;
  char **cb;
  bool *bb;
  
  nb=get_strings("phbres.dat",&cb);
  snew(bb,nres);
  
  for(i=0; (i<nres); i++) {
    if (search_str(nb,cb,resnm[i]) != -1)
      bb[i]=TRUE;
  }
  return bb;
}
 
void wheel(char *fn,int nres,char *resnm[],int r0,real rot0,char *title)
{
  const real fontsize  = 16;
  const real gray      = 0.9;
  const real fontasp   = 0.6;
  const real fontwidth = fontsize*fontasp;
  
  FILE *out;
  int  i,sl,slen;
  real ring,inner,outer;
  real xc,yc,box;
  bool *bPh;
  char **rnms;
  char sign;
  
  inner=75.0;
  slen=0;
  snew(rnms,nres);
  for(i=0; (i<nres); i++) {
    snew(rnms[i],256);
    sl=strlen(resnm[i]);
    sign=resnm[i][sl-1];
    if ((sign == '+') || (sign == '-'))
      resnm[i][sl-1] = '\0';
    sprintf(rnms[i],"%s-%d",resnm[i],i+r0);
    if ((sign == '+') || (sign == '-')) {
      sl=strlen(rnms[i]);
      rnms[i][sl]=sign;
      rnms[i][sl+1]='\0';
    }
    
    slen=max(slen,(int)strlen(rnms[i]));
  }
  ring=(2+slen)*fontwidth;
  outer=inner+ring;
  box=inner*1.5+(1+(nres / 18))*ring;
  
  bPh=bPhobics(nres,resnm);

  out=ps_open(fn,0,0,2.0*box,2.0*box);
  xc=box;
  yc=box;

  ps_font(out,efontHELV,1.5*fontsize);  
  ps_translate(out,xc,yc);
  if (title) 
    ps_ctext(out,0,-fontsize*1.5/2.0,title,eXCenter);
  ps_font(out,efontHELV,fontsize);  
  
  fprintf(out,"%g rotate\n",rot0);
  for(i=0; (i<nres); ) {
    if (bPh[i]) {
      ps_color(out,gray,gray,gray);
      ps_fillarcslice(out,0,0,inner,outer,-10,10);
      ps_color(out,0,0,0);
    }
    ps_arcslice(out,0,0,inner,outer,-10,10);
    
    ps_ctext(out,inner+fontwidth,-fontsize/2.0,rnms[i],eXLeft);
    fprintf(out,"-100 rotate\n");
    i++;
    
    if ((i % 18) == 0) {
      inner=outer;
      outer+=ring;
    }
  }
  ps_close(out);
}

void wheel2(char *fn,int nres,char *resnm[],int r0,real rot0,char *title)
{
  const real fontsize  = 14;
  const real gray      = 0.9;
  const real fontasp   = 0.45;
  const int  angle     = 9;
  const real fontwidth = fontsize*fontasp;
  
  FILE *out;
  int  i,slen;
  real ring,inner,outer;
  real xc,yc,box;
  
  inner=60.0;
  slen=0;
  for(i=0; (i<nres); i++) {
    slen=max(slen,(int)strlen(resnm[i]));
  }
  fprintf(stderr,"slen = %d\n",slen);
  ring=(slen)*fontwidth;
  outer=inner+ring;
  box=(1+(nres / (2*angle)))*outer;
  
  out=ps_open(fn,0,0,2.0*box,2.0*box);
  xc=box;
  yc=box;

  ps_font(out,efontHELV,1.5*fontsize);  
  ps_translate(out,xc,yc);
  ps_color(out,0,0,0);
  if (title) 
    ps_ctext(out,0,-fontsize*1.5/2.0,title,eXCenter);
  ps_font(out,efontHELV,fontsize);  
  
  fprintf(out,"%g rotate\n",rot0);
  for(i=0; (i<nres); ) {
    if ((i % 5) == 4) {
      ps_color(out,gray,gray,1.0);
      ps_fillarcslice(out,0,0,inner,outer,-angle,angle);
      ps_color(out,0,0,0);
    }
    ps_arcslice(out,0,0,inner,outer,-angle,angle);
    
    ps_ctext(out,inner+fontwidth,-fontsize/2.0,resnm[i],eXLeft);
    fprintf(out,"%d rotate\n",-2*angle);
    i++;
    
    if ((i % (2*angle)) == 0) {
      inner=outer;
      outer+=ring;
    }
  }
  ps_close(out);
}

int main(int argc,char *argv[])
{
  static char *desc[] = {
    "wheel plots a helical wheel representation of your sequence."
    "The input sequence is in the .dat file where the first line contains",
    "the number of residues and each consecutive line contains a residue"
    "name."
  };
  static real rot0=0;
  static bool bNum=TRUE;
  static char *title=NULL;
  static int  r0=1;
  t_pargs pa [] = {
    { "-r0",  FALSE, etINT, {&r0},
      "The first residue number in the sequence" },
    { "-rot0",FALSE, etREAL,{&rot0},
      "Rotate around an angle initially (90 degrees makes sense)" },
    { "-T",   FALSE, etSTR, {&title},
      "Plot a title in the center of the wheel (must be shorter than 10 characters, or it will overwrite the wheel)" },
    { "-nn",  FALSE, etBOOL,{&bNum},
      "Toggle numbers" }
  };
  t_filenm  fnm[] = {
    { efDAT, "-f", NULL,  ffREAD  },
    { efEPS, "-o", NULL,  ffWRITE }
  };
#define NFILE asize(fnm)
  
  int  i,nres;
  char **resnm;
  
  CopyRight(stderr,argv[0]);
  parse_common_args(&argc,argv,0,TRUE,NFILE,fnm,asize(pa),pa,
		    asize(desc),desc,0,NULL);
  
  for(i=1; (i<argc); i++) {
    if (strcmp(argv[i],"-r0") == 0) {
      r0=atoi(argv[++i]);
      fprintf(stderr,"First residue is %d\n",r0);
    }
    else if (strcmp(argv[i],"-rot0") == 0) {
      rot0=atof(argv[++i]);
      fprintf(stderr,"Initial rotation is %g\n",rot0);
    }
    else if (strcmp(argv[i],"-T") == 0) {
      title=strdup(argv[++i]);
      fprintf(stderr,"Title will be '%s'\n",title);
    }
    else if (strcmp(argv[i],"-nn") == 0) {
      bNum=FALSE;
      fprintf(stderr,"No residue numbers\n");
    }
    else
      fatal_error(0,"Incorrect usage of option %s",argv[i]);
  }
    
  nres=get_lines(ftp2fn(efDAT,NFILE,fnm),&resnm);
  if (bNum)
    wheel(ftp2fn(efEPS,NFILE,fnm),nres,resnm,r0,rot0,title);
  else
    wheel2(ftp2fn(efEPS,NFILE,fnm),nres,resnm,r0,rot0,title);
    
  thanx(stdout);
  
  return 0;
}
