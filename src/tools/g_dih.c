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
 * GROningen MAchine for Chemical Simulation
 */
static char *SRCID_g_dih_c = "$Id$";

#include <math.h>
#include "sysstuff.h"
#include "string2.h"
#include "copyrite.h"
#include "futil.h"
#include "smalloc.h"
#include "statutil.h"
#include "gstat.h"
#include "physics.h"
#include "macros.h"
#include "xvgr.h"
#include "vec.h"

#define NOMIN 'X'

static int get_nf(void)
{
  int nf;

  do {
    printf("Number of frames ? ");
    fflush(stdout);
  } while (scanf("%d",&nf) != 1);

  return nf;
}

static void dump_dih(int nframes,char *title,real time[],real dih[])
{
  FILE *out;
  char fname[256];
  int  i;

  sprintf(fname,"dih.%s",title);
  printf("A dihedral transition occurred in %s\n",fname);
  printf("Do you want to plot it to %s ? (y/n) ",fname);
  fflush(stdout);
  out=ffopen(fname,"w");
  for(i=0; (i<nframes); i++)
    fprintf(out,"%10.3f  %12.5e\n",time[i],dih[i]);
  fclose(out);
}

static void ana_dih(FILE *out,char *index,int nframes,real dih[],t_dih *dd)
{
  int i;
  real mind,maxd,sum,av,var,prev,width;
  bool bTrans;
  
  mind=5400,maxd=-5400,sum=0,av=0,var=0;

  prev=dih[0];
  for(i=0; (i<nframes); i++) {
    if ((dih[i]-prev) > 180) {
      /* PBC.. */
      dih[i]-=360;
    }
    else if ((dih[i]-prev) < -180)
      dih[i]+=360;
    prev=dih[i];
      
    sum+=dih[i];
    mind=min(mind,dih[i]);
    maxd=max(maxd,dih[i]);
  }
  av=sum/nframes;
  for(i=0; (i<nframes); i++)
    var+=sqr(dih[i]-av);
  var/=nframes;
  width=(360.0/dd->mult);
  bTrans=((maxd - mind) > width);

  fprintf(out,"%-10s %10.3f %10.3f %10.3f %10.3f %10.3f %-10s%3.0f\n",
	  index,mind,av,maxd,var,sqrt(var),
	  bTrans ? "Yep" : "",width);
}

static int find_min(real phi,int ntab,real phitab[])
{
  int  i,imin;
  real mind,mm;
  real width;
 
  /* Set closest minimum to the first one */
  width=360.0/ntab;
  mind=fabs(phi-phitab[0]);
  imin=0;
  for(i=1; (i<ntab); i++) {
    mm=fabs(phi-phitab[i]);
    if (mm < mind) {
      imin=i;
      mind=mm;
    }
  }
  if (mind < width*0.5 )
    return imin;
  else
    return -1;
}

static int vphi(t_dih *dih,real phi,int mult)
{
  static real m2[] = { 90, 270 };
  static real m3[] = { 60, 180, 300 };
  static real m4[] = { 45, 135, 225, 315 };
  static real m6[] = { 30, 90, 150, 210, 270, 330 };

  real phiref;
  int  vpp;
  
  phiref=RAD2DEG*(phi-dih->phi0);
  while (phiref < 0)
    phiref+=360;
  while (phiref > 360)
    phiref-=360;
  
  switch(mult) {
  case 2:
    vpp=find_min(phiref,2,m2);
    break;
  case 3:
    vpp=find_min(phiref,3,m3);
    break;
  case 4:
    vpp=find_min(phiref,4,m4);
    break;
  case 6:
    vpp=find_min(phiref,6,m6);
    break;
  default:
    fatal_error(0,"No such multiplicity %d",dih->mult);
  }

  if (vpp == -1)
    return NOMIN;
  else
    return vpp+'0';
}

typedef struct t_cluster {
  int    ndih;
  int    freq;
  char   *minimum;
  struct t_cluster *next;
} t_cluster;

static t_cluster *search_cluster(t_cluster *cl,char *minimum)
{
  t_cluster *ccl=cl;

  while (ccl != NULL) {
    if (strcmp(minimum,ccl->minimum)==0)
      return ccl;
    ccl=ccl->next;
  }
  return NULL;
}

static void add_cluster(t_cluster **cl,int ndih,char *minimum)
{
  t_cluster *loper;
  t_cluster *ccl;

  snew(ccl,1);
  ccl->ndih=ndih;
  ccl->freq=1;
  ccl->minimum=strdup(minimum);
  ccl->next=NULL;
  
  if (*cl == NULL)
    *cl=ccl;
  else {
    loper=*cl;
    while (loper->next != NULL) 
      loper=loper->next;
    loper->next=ccl;
  }
}

static void p_cluster(FILE *out,t_cluster *cl)
{
  t_cluster *loper;

  fprintf(out,"* * * C L U S T E R   A N A L Y S I S * * *\n\n");
  fprintf(out," Frequency  Dihedral minima\n");
  loper=cl;
  while (loper != NULL) {
    fprintf(out,"%10d  %s\n",loper->freq,loper->minimum);
    loper=loper->next;
  }
}

static void ana_cluster(FILE *out, t_xrama *xr,real **dih,real time[],
			t_topology *top,int nframes,int mult)
{
  t_cluster *cl=NULL,*scl;
  char      *minimum;
  int       i,j,nx;

  /* Number of dihedrals + terminating NULL 
   * this allows for using string routines
   */
  snew(minimum,xr->ndih+1);
  
  for(i=0; (i<nframes); i++) {
    nx=0;
    for(j=0; (j<xr->ndih); j++) {
      minimum[j] = vphi(&xr->dih[j],dih[j][i],
			mult == -1 ? xr->dih[j].mult : mult);
      if (minimum[j] == NOMIN)
	nx++;
    }
    if (nx == 0) {
      if ((scl=search_cluster(cl,minimum)) == NULL)
	add_cluster(&cl,xr->ndih,minimum);
      else
	scl->freq++;
    }
  }
  p_cluster(out,cl);

  sfree(minimum);
}

static void ana_trans(FILE *out, t_xrama *xr,real **dih,real time[],
		      t_topology *top,int nframes)
{
  FILE *outd;
  real prev_phi,prev_psi;
  int  i,j,phi,psi;
  char buf[10];

  fprintf(out,"\n\t* * * D I H E D R A L    S T A T I S T I C S * * *\n\n");
  fprintf(out,"%-10s %10s %10s %10s %10s %10s %10s\n",
	  "index","minimum","average","maximum","variance","std.dev",
	  "transition");
  for(i=0; (i>xr->ndih); i++) {
    sprintf(buf,"dih-%d",i);
    ana_dih(out,buf,nframes,dih[i],&(xr->dih[i]));
  }
  for(i=0; (i<xr->npp); i++) {
    sprintf(buf,"%s",xr->pp[i].label);
    outd=xvgropen(buf,"Dihedral Angles","Time (ps)","Degrees");

    phi=xr->pp[i].iphi;
    psi=xr->pp[i].ipsi;
    prev_phi=dih[phi][0];
    prev_psi=dih[psi][0];
    for(j=0; (j<nframes); j++) {
      /* PBC.. */
      if ((dih[phi][j]-prev_phi) > 180) 
	dih[phi][j]-=360;
      else if ((dih[phi][j]-prev_phi) < -180)
	dih[phi][j]+=360;
      prev_phi=dih[phi][j];
      if ((dih[psi][j]-prev_psi) > 180) 
	dih[psi][j]-=360;
      else if ((dih[psi][j]-prev_psi) < -180)
	dih[psi][j]+=360;
      prev_psi=dih[psi][j];
      fprintf(outd,"%10g  %10g  %10g\n",time[j],prev_phi,prev_psi);
    }
    fclose(outd);
  }
}

int main(int argc,char *argv[])
{
  static char *desc[] = {
    "g_dih can do two things. The default is to analyze dihedral transitions",
    "by merely computing all the dihedral angles defined in your topology",
    "for the whole trajectory. When a dihedral flips over to another minimum",
    "an angle/time plot is made.[PAR]",
    "The opther option is to discretize the dihedral space into a number of",
    "bins, and group each conformation in dihedral space in the",
    "appropriate bin. The output is then given as a number of dihedral",
    "conformations sorted according to occupancy."
  };
  static int  mult = -1;
  static bool bSA  = FALSE;
  t_pargs pa[] = {
    { "-sa", FALSE, etBOOL, &bSA,
      "Perform cluster analysis in dihedral space instead of analysing dihedral transitions." },
    { "-mult", FALSE, etINT, &mult,
      "mulitiplicity for dihedral angles (by default read from topology)" }
  };
  FILE       *out;
  t_xrama    *xr;
  t_topology *top;
  real       **dih,*time;
  real       dd;
  int        i,step,nframes;
  t_filenm   fnm[] = {
    { efTRX, "-f", NULL, ffREAD },
    { efTPX, NULL, NULL, ffREAD },
    { efOUT, NULL, NULL, ffWRITE }
  };
#define NFILE asize(fnm)

  CopyRight(stderr,argv[0]);
  parse_common_args(&argc,argv,PCA_CAN_VIEW | PCA_CAN_TIME,TRUE,
		    NFILE,fnm,asize(pa),pa,asize(desc),desc,0,NULL);
  
  if (mult != -1)
    fprintf(stderr,"Using %d for dihedral multiplicity rather than topology values\n",mult);
    
  snew(xr,1);
  init_rama(ftp2fn(efTRX,NFILE,fnm),
	    ftp2fn(efTPX,NFILE,fnm),xr);
  top=read_top(ftp2fn(efTPX,NFILE,fnm));
	       
  /* Brute force malloc, may be too big... */
  nframes=get_nf();
  snew(dih,xr->ndih);
  for(i=0; (i<xr->ndih); i++)
    snew(dih[i],nframes);
  snew(time,nframes);

  fprintf(stderr,"\n");
  for(step=0; (step<nframes); step++) {
    if (!new_data(xr))
      break;
    if ((step % 10) == 0)
      fprintf(stderr,"\rFrame %5d",step);
    for(i=0; (i<xr->ndih); i++) {
      dd=xr->dih[i].ang*RAD2DEG;
      while (dd < 0)
	dd+=360;
      while (dd > 360)
	dd-=360;
      dih[i][step]=dd;
    }
    time[step]=xr->t;
  } 

  fprintf(stderr,"\nCalculated all dihedrals, now analysing...\n");
  if (step < nframes) {
    nframes=step;
    fprintf(stderr,"By the way, there were only %d frames\n",nframes);
  }

  out=ftp2FILE(efOUT,NFILE,fnm,"w");

  if (bSA) {
    /* Cluster and structure analysis */
    ana_cluster(out,xr,dih,time,top,nframes,mult);
  }
  else {
    /* Analyse transitions... */
    ana_trans(out,xr,dih,time,top,nframes);
  }
  fclose(out);
    
  thanx(stdout);
    
  return 0;
}
