#include <math.h>
#include "typedefs.h"
#include "smalloc.h"
#include "vec.h"
#include "nrjac.h"
#include "network.h"
#include "orires.h"

void init_orires(FILE *log,int nfa,t_iatom forceatoms[],t_iparams ip[],
		 t_inputrec *ir,t_commrec *mcr,t_fcdata *fcd)
{
  int i,j,ex,nr,*nr_ex;
  t_oriresdata *od;
  
  od = &(fcd->orires);
  od->fc  = ir->orires_fc;
  od->nex = 0;
  od->S   = NULL;

  if (ir->orires_tau > 0)
    od->edt = exp(-ir->delta_t/ir->orires_tau);
  else
    od->edt = 0;
  od->edt1 = 1 - od->edt;
  od->exp_min_t_tau = 1.0;
  od->nr = nfa/3;
  
  if (od->nr == 0)
    return;

  nr_ex = NULL;

  for(i=0; i<nfa; i+=3) {
    ex = ip[forceatoms[i]].orires.ex;
    if (ex >= od->nex) {
      srenew(nr_ex,ex+1);
      for(j=od->nex; j<ex+1; j++)
	nr_ex[j] = 0;
      od->nex = ex+1;
    }
    nr_ex[ex]++;
  }
  snew(od->S,od->nex);
  /* When not doing time averaging, the instaneous and time averaged data
   * are indentical and the pointers can point to the same memory.
   */
  snew(od->Dins,od->nr);
  if (ir->orires_tau == 0)
    od->Dtav = od->Dins;
  else
    snew(od->Dtav,od->nr);
  if (ir->orires_tau == 0) {
    snew(od->oins,od->nr);
    od->otav = od->oins;
  } else {
    snew(od->oins,2*od->nr);
    od->otav = &(od->oins[od->nr]);
  }
  /* When not ensemble averaging the local orientations are identical to
   * the ensemble averaged ones and the pointer can point to the same memory.
   */
  if (mcr)
    snew(od->oinsl,od->nr);
  else
    od->oinsl = od->oins;
  snew(od->tmp,od->nex);
  snew(od->TMP,od->nex);
  for(ex=0; ex<od->nex; ex++) {
    snew(od->TMP[ex],5);
    for(i=0; i<5; i++)
      snew(od->TMP[ex][i],5);
  }

  fprintf(log,"Found %d orientation experiments\n",od->nex);
  for(i=0; i<od->nex; i++)
    fprintf(log,"  experiment %d has %d restraints\n",i+1,nr_ex[i]);
  if (mcr)
    fprintf(log,"The orientation restraints are ensemble averaged over %d systems\n",mcr->nnodes);

  sfree(nr_ex);
}

void print_orires_log(FILE *log,t_fcdata *fcd)
{
  int           ex,i,j,nrot;
  bool          bZero;
  t_oriresdata  *od;
  static double **M=NULL,*eig,**v;
  
  od = &(fcd->orires);

  if (M == NULL) {
    snew(M,DIM);
    for(i=0; i<DIM; i++)
      snew(M[i],DIM);
    snew(eig,DIM);
    snew(v,DIM);
    for(i=0; i<DIM; i++)
      snew(v[i],DIM);
  }

  for(ex=0; ex<od->nex; ex++) {
    bZero = TRUE;
    for(i=0; i<DIM; i++)
      for(j=0; j<DIM; j++) {
	bZero = bZero && od->S[ex][i][j]==0;
	M[i][j] = od->S[ex][i][j];
      }
    if (!bZero) {
      jacobi(M,DIM,eig,v,&nrot);
      
      j=0;
      for(i=1; i<DIM; i++)
	if (sqr(eig[i]) > sqr(eig[j]))
	  j=i;
      
      fprintf(log,"  Orientation experiment %d:\n",ex+1);
      fprintf(log,"    order parameter: %g\n",eig[j]);
      for(i=0; i<DIM; i++)
	fprintf(log,"    eig: %6.3f   %6.3f %6.3f %6.3f\n",
		eig[i]/eig[j],v[XX][i],v[YY][i],v[ZZ][i]);
      fprintf(log,"\n");
    }
  }
}

real calc_orires_dev(t_commrec *mcr,
		     int nfa,t_iatom forceatoms[],t_iparams ip[],
		     rvec x[],t_forcerec *fr,t_fcdata *fcd)
{
  int          fa,d,i,j,type,ex;
  real         edt,edt1,pfac,r2,invr,corrfac,weight,wsv2,sw,dev;
  real         two_thn;
  tensor       *S;
  rvec5        *Dins,*Dtav,*rhs;
  real         ***T;
  rvec         r;
  t_oriresdata *od;
  bool         bTAV;

  od = &(fcd->orires);

  bTAV = (od->edt != 0);
  edt  = od->edt;
  edt1 = od->edt1;
  S    = od->S;
  Dins = od->Dins;
  Dtav = od->Dtav;
  T    = od->TMP;
  rhs  = od->tmp;
  
  od->exp_min_t_tau *= edt;

  d = 0;
  for(fa=0; fa<nfa; fa+=3) {
    type = forceatoms[fa];
    rvec_sub(x[forceatoms[fa+1]],x[forceatoms[fa+2]],r);
    r2   = norm2(r);
    invr = invsqrt(r2);
    /* Calculate the prefactor for the D tensor, this includes the factor 3! */
    pfac = ip[type].orires.c*invr*invr*3;
    for(i=0; i<ip[type].orires.pow; i++)
      pfac *= invr;
    Dins[d][0] = pfac*(2*r[0]*r[0] + r[1]*r[1] - r2);
    Dins[d][1] = pfac*(2*r[0]*r[1]);
    Dins[d][2] = pfac*(2*r[0]*r[2]);
    Dins[d][3] = pfac*(2*r[1]*r[1] + r[0]*r[0] - r2);
    Dins[d][4] = pfac*(2*r[1]*r[2]);
    if (bTAV)
      for(i=0; i<5; i++)
	Dtav[d][i] = edt*Dtav[d][i] + edt1*Dins[d][i];
    d++;
  }

  /* Correction factor to correct for the lack of history for short times */
  corrfac = 1.0/(1.0-od->exp_min_t_tau);

  /* Calculate the order tensor S for each experiment via optimization */
  for(ex=0; ex<od->nex; ex++)
    for(i=0; i<5; i++) {
      rhs[ex][i] = 0;
      for(j=0; j<=i; j++)
	T[ex][i][j] = 0;
    }
  d = 0;
  for(fa=0; fa<nfa; fa+=3) {
    type   = forceatoms[fa];
    ex     = ip[type].orires.ex;
    weight = ip[type].orires.kfac;
    /* Calculate the vector rhs and half the matrix T for the 5 equations */
    for(i=0; i<5; i++) {
      rhs[ex][i] += Dtav[d][i]*ip[type].orires.obs*weight;
      for(j=0; j<=i; j++)
	T[ex][i][j] += Dtav[d][i]*Dtav[d][j]*weight;
    }
    d++;
  }
  /* Now we have all the data we can calculate S */
  for(ex=0; ex<od->nex; ex++) {
    /* Correct corrfac and copy one half of T to the other half */
    for(i=0; i<5; i++) {
      rhs[ex][i]  *= corrfac;
      T[ex][i][i] *= sqr(corrfac);
      for(j=0; j<i; j++) {
	T[ex][i][j] *= sqr(corrfac);
	T[ex][j][i]  = T[ex][i][j];
      }
    }
    m_inv_gen(T[ex],5,T[ex]);
    /* Calculate the orientation tensor S for this experiment */
    S[ex][0][0] = 0;
    S[ex][0][1] = 0;
    S[ex][0][2] = 0;
    S[ex][1][1] = 0;
    S[ex][1][2] = 0;
    for(i=0; i<5; i++) {
      S[ex][0][0] += 1.5*T[ex][0][i]*rhs[ex][i];
      S[ex][0][1] += 1.5*T[ex][1][i]*rhs[ex][i];
      S[ex][0][2] += 1.5*T[ex][2][i]*rhs[ex][i];
      S[ex][1][1] += 1.5*T[ex][3][i]*rhs[ex][i];
      S[ex][1][2] += 1.5*T[ex][4][i]*rhs[ex][i];
    }
    S[ex][1][0] = S[ex][0][1];
    S[ex][2][0] = S[ex][0][2];
    S[ex][2][1] = S[ex][1][2];
    S[ex][2][2] = -S[ex][0][0] - S[ex][1][1];
  }

  if (mcr)
    two_thn = 2.0/(3.0*mcr->nnodes);
  else
    two_thn = 2.0/3.0;
  
  d = 0;
  for(fa=0; fa<nfa; fa+=3) {
    type = forceatoms[fa];
    ex = ip[type].orires.ex;

    od->otav[d] = two_thn*
      corrfac*(S[ex][0][0]*Dtav[d][0] + S[ex][0][1]*Dtav[d][1] +
	       S[ex][0][2]*Dtav[d][2] + S[ex][1][1]*Dtav[d][3] +
	       S[ex][1][2]*Dtav[d][4]);
    if (bTAV)
      od->oins[d] = two_thn*(S[ex][0][0]*Dins[d][0] + S[ex][0][1]*Dins[d][1] +
			     S[ex][0][2]*Dins[d][2] + S[ex][1][1]*Dins[d][3] +
			     S[ex][1][2]*Dins[d][4]);
    if (mcr)
      /* When ensemble averaging is used recalculate the local orientation
       * for output to the energy file.
       */
      od->oinsl[d] = mcr->nnodes*od->oins[d];

    d++;
  }

  if (mcr)
    /* Parallel call to ensemble average the orientations by summing them
     * over the processors.
     * To minimize communication when using time averaging the time averages
     * are stored in memory behind the instantaneous values:
     * od->otav = &(od->oins[od->nr])
     */
    gmx_sum(bTAV ? 2*od->nr : od->nr,od->oins,mcr);
  
  wsv2 = 0;
  sw   = 0;
  
  d = 0;
  for(fa=0; fa<nfa; fa+=3) {
    type = forceatoms[fa];
    
    dev = od->otav[d] - ip[type].orires.obs;
    
    wsv2 += ip[type].orires.kfac*sqr(dev);
    sw   += ip[type].orires.kfac;
    
    d++;
  }
  od->rmsdev = sqrt(wsv2/sw);
  
  return od->rmsdev;
  
  /* Approx. 120*nfa/3 flops */
}

real orires(int nfa,t_iatom forceatoms[],t_iparams ip[],
	    rvec x[],rvec f[],t_forcerec *fr,t_graph *g,
	    matrix box,real lambda,real *dvdlambda,
	    t_mdatoms *md,int ngrp,real egnb[],real egcoul[],
	    t_fcdata *fcd)
{
  atom_id      ai,aj;
  int          fa,d,i,type,ex,power,ki;
  ivec         dt;
  real         r2,invr,invr2,fc,smooth_fc,dev,devins,pfac;
  rvec         r,Sr,fij;
  real         vtot;
  t_oriresdata *od;
  bool         bTAV;

  vtot = 0;
  od = &(fcd->orires);

  if (od->fc != 0) {
    bTAV = (od->edt != 0);
    
    /* Smoothly switch on the restraining when time averaging is used */
    smooth_fc = od->fc*(1.0 - od->exp_min_t_tau);
    
    d = 0;
    for(fa=0; fa<nfa; fa+=3) {
      type  = forceatoms[fa];
      ai    = forceatoms[fa+1];
      aj    = forceatoms[fa+2];
      rvec_sub(x[ai],x[aj],r);
      r2    = norm2(r);
      invr  = invsqrt(r2);
      invr2 = invr*invr;
      ex    = ip[type].orires.ex;
      power = ip[type].orires.pow;
      fc    = smooth_fc*ip[type].orires.kfac;
      dev   = od->otav[d] - ip[type].orires.obs;
      
      /* NOTE: there is no real potential when time averaging is applied */
      vtot += 0.5*fc*sqr(dev);
      
      if (bTAV) {
	/* Calculate the force as the sqrt of tav times instantaneous */
	devins = od->oins[d] - ip[type].orires.obs;
	if (dev*devins <= 0)
	  dev = 0;
	else {
	  dev = sqrt(dev*devins);
	  if (dev < 0)
	    dev = -dev;
	}
      }
      
      pfac  = fc*ip[type].orires.c*invr2;
      for(i=0; i<power; i++)
	pfac *= invr;
      mvmul(od->S[ex],r,Sr);
      for(i=0; i<DIM; i++)
	fij[i] = -pfac*dev*(4*Sr[i] - 2*(2+power)*invr2*iprod(Sr,r)*r[i]);
      
      ivec_sub(SHIFT_IVEC(g,ai),SHIFT_IVEC(g,aj),dt);
      ki=IVEC2IS(dt);
      
      for(i=0; i<DIM; i++) {
	f[ai][i]               += fij[i];
	f[aj][i]               -= fij[i];
	fr->fshift[ki][i]      += fij[i];
	fr->fshift[CENTRAL][i] -= fij[i];
      }
      d++;
    }
  }
  
  return vtot;
  
  /* Approx. 80*nfa/3 flops */
}
