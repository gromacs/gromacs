#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "typedefs.h"
#include "smalloc.h"
#include "macros.h"
#include "copyrite.h"
#include "statutil.h"
#include "fatal.h"
#include "random.h"
#include "pdbio.h"
#include "futil.h"
#include "physics.h"
#include "xvgr.h"
#include "vec.h"
#include "names.h"
#include "ehdata.h"

t_histo *init_histo(int np,real minx,real maxx)
{
  t_histo *h;
  
  snew(h,1);
  snew(h->y,np+1);
  snew(h->nh,np+1);
  h->np   = np;
  if (maxx <= minx)
    fatal_error(0,"minx (%f) should be less than maxx (%f) in init_histo",minx,maxx);
  h->minx = minx;
  h->maxx = maxx;
  h->dx_1 = np/(maxx-minx);
  h->dx   = (maxx-minx)/np;
  
  return h;
}

void done_histo(t_histo *h)
{
  sfree(h->y);
  sfree(h->nh);
  h->np = 0;
}

void add_histo(t_histo *h,real x,real y)
{
  int n;
  
  n = (x-h->minx)*h->dx_1;
  if ((n < 0) || (n > h->np)) 
    fatal_error(0,"Invalid x (%f) in add_histo. SHould be in %f - %f",x,h->minx,h->maxx);
  h->y[n] += y;
  h->nh[n]++;
}

void dump_histo(t_histo *h,char *fn,char *title,char *xaxis,char *yaxis,
		int enorm,real norm_fac)
{
  FILE *fp;
  int  i,nn;
  
  for(nn=h->np; (nn > 0); nn--)
    if (h->nh[nn] != 0) 
      break;
  for(i=0; (i<nn); i++)
    if (h->nh[i] > 0)
      break;
  fp = xvgropen(fn,title,xaxis,yaxis);
  for(  ; (i<nn); i++) {
    switch (enorm) {
    case enormNO:
      fprintf(fp,"%12f  %12f  %d\n",h->minx+h->dx*i,h->y[i],h->nh[i]);
      break;
    case enormFAC:
      fprintf(fp,"%12f  %12f  %d\n",h->minx+h->dx*i,h->y[i]*norm_fac,h->nh[i]);
      break;
    case enormNP:
      if (h->nh[i] > 0)
	fprintf(fp,"%12f  %12f  %d\n",h->minx+h->dx*i,h->y[i]*norm_fac/h->nh[i],h->nh[i]);
      break;
    default:
      fatal_error(0,"Wrong value for enorm (%d)",enorm);
    }
  }
  fclose(fp);
}

/*******************************************************************
 *
 * Functions to analyse and monitor scattering
 *
 *******************************************************************/	

void add_scatter_event(t_ana_scat *scatter,rvec pos,bool bInel,
		       real t,real ekin)
{
  int np = scatter->np;
  
  if (np == scatter->maxp) {
    scatter->maxp += 32;
    srenew(scatter->time,scatter->maxp);
    srenew(scatter->ekin,scatter->maxp);
    srenew(scatter->bInel,scatter->maxp);
    srenew(scatter->pos,scatter->maxp);
  }
  scatter->time[np]  = t;
  scatter->bInel[np] = np;
  scatter->ekin[np]  = ekin;
  copy_rvec(pos,scatter->pos[np]);
  scatter->np++;
}

void reset_ana_scat(t_ana_scat *scatter)
{
  scatter->np = 0;
}

void done_scatter(t_ana_scat *scatter)
{
  scatter->np = 0;
  sfree(scatter->time);
  sfree(scatter->ekin);
  sfree(scatter->bInel);
  sfree(scatter->pos);
}

void analyse_scatter(t_ana_scat *scatter,t_histo *hmfp)
{
  int   i,n;
  rvec  dx;
  
  if (scatter->np > 1) {
    for(i=1; (i<scatter->np); i++) {
      rvec_sub(scatter->pos[i],scatter->pos[i-1],dx);
      add_histo(hmfp,scatter->ekin[i],norm(dx));
    }
  }
}

/*******************************************************************
 *
 * Functions to analyse structural changes
 *
 *******************************************************************/	

t_ana_struct *init_ana_struct(int nstep,int nsave,real timestep)
{
  t_ana_struct *anal;
  
  snew(anal,1);
  anal->nanal = (nstep / nsave)+1;
  anal->index = 0;
  anal->dt    = nsave*timestep;
  snew(anal->t,anal->nanal);
  snew(anal->maxdist,anal->nanal);
  snew(anal->averdist,anal->nanal);
  snew(anal->ad2,anal->nanal);
  snew(anal->nion,anal->nanal);
  
  return anal;
}

void done_ana_struct(t_ana_struct *anal)
{
  sfree(anal->t);
  sfree(anal->maxdist);
  sfree(anal->averdist);
  sfree(anal->ad2);
  sfree(anal->nion);
}

void reset_ana_struct(t_ana_struct *anal)
{
  int i;
  
  for(i=0; (i<anal->nanal); i++) {
    anal->t[i] = 0;
    anal->maxdist[i] = 0;
    anal->averdist[i] = 0;
    anal->ad2[i] = 0;
    anal->nion[i] = 0;
  }
  anal->index = 0;
}

void add_ana_struct(t_ana_struct *total,t_ana_struct *add)
{
  int i;
  
  if (total->index == 0)
    total->index = add->index;
  else if (total->index != add->index)
    fatal_error(0,"Analysis incompatible %s, %d",__FILE__,__LINE__);
  for(i=0; (i<total->index); i++) {
    if (total->t[i] == 0)
      total->t[i] = add->t[i];
    else if (total->t[i] != add->t[i])
      fatal_error(0,"Inconsistent times in analysis (%f-%f) %s, %d",
		  total->t[i],add->t[i],__FILE__,__LINE__);
    total->maxdist[i]  += add->maxdist[i];
    total->averdist[i] += add->averdist[i];
    total->ad2[i]      += add->ad2[i];
    total->nion[i]     += add->nion[i];
  }
}

void analyse_structure(t_ana_struct *anal,real t,rvec center,
		       rvec x[],int nparticle,real charge[])
{
  int  i,j,n=0;
  rvec dx;
  real dx2,dx1;
  
  j = anal->index;
  anal->t[j]       = t;
  anal->maxdist[j] = 0;
  for(i=0; (i<nparticle); i++) {
    if (charge[i] < 0) {
      rvec_sub(x[i],center,dx);
      dx2 = iprod(dx,dx);
      dx1 = sqrt(dx2);
      anal->ad2[j] += dx2;
      anal->averdist[j]  += dx1;
      if (dx1 > anal->maxdist[j])
	anal->maxdist[j] = dx1;
      n++;
    }
  }
  anal->nion[j] = n;
  anal->index++;
}

void dump_ana_struct(char *rmax,char *nion,char *gyr,
		     t_ana_struct *anal,int nsim)
{
  FILE *fp,*gp,*hp;
  int  i;
  real t;
  
  fp = xvgropen(rmax,"rmax","Time (fs)","r (nm)");
  gp = xvgropen(nion,"N ion","Time (fs)","N ions");
  hp = xvgropen(gyr,"Radius of gyration","Time (fs)","Rg (nm)");
  for(i=0; (i<anal->index); i++) {
    t = 1000*anal->t[i];
    fprintf(fp,"%12g  %12.3f\n",t,anal->maxdist[i]/nsim);
    fprintf(gp,"%12g  %12.3f\n",t,(1.0*anal->nion[i])/nsim);
    if (anal->nion[i] > 0)
      fprintf(hp,"%12g  %12.3f\n",t,sqrt(anal->ad2[i]/anal->nion[i]));
  }
  fclose(hp);
  fclose(gp);
  fclose(fp);
}
