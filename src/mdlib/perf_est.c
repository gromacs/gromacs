#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "perf_est.h"
#include "physics.h"
#include "vec.h"
#include "mtop_util.h"

int n_bonded_dx(gmx_mtop_t *mtop,bool bExcl)
{
  int mb,nmol,ftype,ndxb,ndx_excl;
  int ndx;
  gmx_moltype_t *molt;

  /* Count the number of pbc_rvec_sub calls required for bonded interactions.
   * This number is also roughly proportional to the computational cost.
   */
  ndx = 0;
  ndx_excl = 0;
  for(mb=0; mb<mtop->nmolblock; mb++) {
    molt = &mtop->moltype[mtop->molblock[mb].type];
    nmol = mtop->molblock[mb].nmol;
    for(ftype=0; ftype<F_NRE; ftype++) {
      if (interaction_function[ftype].flags & IF_BOND) {
	switch (ftype) {
	F_POSRES:    ndxb = 1; break;
	F_CONNBONDS: ndxb = 0; break;
	default:     ndxb = NRAL(ftype) - 1; break;
	}
	ndx += nmol*ndxb*molt->ilist[ftype].nr/(1 + NRAL(ftype));
      }
    }
    if (bExcl) {
      ndx_excl += nmol*(molt->excls.nra - molt->atoms.nr)/2;
    } else {
      ndx_excl = 0;
    }
  }

  if (debug)
    fprintf(debug,"ndx bonded %d exclusions %d\n",ndx,ndx_excl);

  ndx += ndx_excl;

  return ndx;
}

float pme_load_estimate(gmx_mtop_t *mtop,t_inputrec *ir,matrix box)
{
  t_atom *atom;
  int  mb,nmol,atnr,cg,a,a0,ncqlj,ncq,nclj;
  bool bBHAM,bLJcut,bWater,bQ,bLJ;
  double nw,nqlj,nq,nlj,cost_bond,cost_pp,cost_spread,cost_fft;
  float fq,fqlj,flj,fljtab,fqljw,fqw,fqspread,ffft,fbond;
  float ratio;
  t_iparams *iparams;
  gmx_moltype_t *molt;

  bBHAM = (mtop->ffparams.functype[0] == F_BHAM);

  bLJcut = ((ir->vdwtype == evdwCUT) && !bBHAM);

  /* Computational cost relative to a tabulated q-q interaction.
   * This will be machine dependent.
   * The numbers here are accurate for Intel Core2 and AMD Athlon 64
   * in single precision. In double precision PME mesh is slightly cheaper,
   * although not so much that the numbers need to be adjusted.
   */
  fq    = 1.0;
  fqlj  = (bLJcut ? 1.5  : 2.0 );
  flj   = (bLJcut ? 0.5  : 1.5 );
  /* Cost of 1 water with one Q/LJ atom */
  fqljw = (bLJcut ? 1.75 : 2.25);
  /* Cost of 1 water with one Q atom or with 1/3 water (LJ negligible) */
  fqw   = 1.5;
  /* Cost of q spreading and force interpolation per charge */
  fqspread = 25.0;
  /* Cost of fft's + pme_solve, will be multiplied with N log(N) */
  ffft     =  0.4;
  /* Cost of a bonded interaction divided by the number of (pbc_)dx required */
  fbond = 5.0;

  iparams = mtop->ffparams.iparams;
  atnr = mtop->ffparams.atnr;
  nw   = 0;
  nqlj = 0;
  nq   = 0;
  nlj  = 0;
  for(mb=0; mb<mtop->nmolblock; mb++) {
    molt = &mtop->moltype[mtop->molblock[mb].type];
    atom = molt->atoms.atom;
    nmol = mtop->molblock[mb].nmol;
    a = 0;
    for(cg=0; cg<molt->cgs.nr; cg++) {
      bWater = !bBHAM;
      ncqlj = 0;
      ncq   = 0;
      nclj  = 0;
      a0    = a;
      while (a < molt->cgs.index[cg+1]) {
	bQ  = (atom[a].q != 0 || atom[a].qB != 0);
	bLJ = (iparams[(atnr+1)*atom[a].type].lj.c6  != 0 ||
	       iparams[(atnr+1)*atom[a].type].lj.c12 != 0);
	/* This if this atom fits into water optimization */
	if (!((a == a0   &&  bQ &&  bLJ) ||
	      (a == a0+1 &&  bQ && !bLJ) ||
	      (a == a0+2 &&  bQ && !bLJ && atom[a].q == atom[a-1].q) ||
	      (a == a0+3 && !bQ &&  bLJ)))
	  bWater = FALSE;
	if (bQ && bLJ) {
	  ncqlj++;
	} else {
	  if (bQ)
	    ncq++;
	  if (bLJ)
	    nclj++;
	}
	a++;
      }
      if (bWater) {
	nw   += nmol;
      } else {
	nqlj += nmol*ncqlj;
	nq   += nmol*ncq;
	nlj  += nmol*nclj;
      }
    }
  }
  if (debug)
    fprintf(debug,"nw %g nqlj %g nq %g nlj %g\n",nw,nqlj,nq,nlj);

  cost_bond = fbond*n_bonded_dx(mtop,TRUE);

  cost_pp = 0.5*(fqljw*nw*nqlj +
		 fqw  *nw*(3*nw + nq) +
		 fqlj *nqlj*nqlj +
		 fq   *nq*(3*nw + nq) +
		 flj  *nlj*(nw + nlj))
    *4/3*M_PI*ir->rlist*ir->rlist*ir->rlist/det(box);
  
  cost_spread = fqspread*(3*nw + nqlj + nq);
  cost_fft    = ffft*ir->nkx*ir->nky*ir->nkz*log(ir->nkx*ir->nky*ir->nkz);
  
  ratio =
    (cost_spread + cost_fft)/(cost_bond + cost_pp + cost_spread + cost_fft);

  if (debug) {
    fprintf(debug,
	    "cost_bond   %f\n"
	    "cost_pp     %f\n"
	    "cost_spread %f\n"
	    "cost_fft    %f\n",
	    cost_bond,cost_pp,cost_spread,cost_fft);

    fprintf(debug,"Estimate for relative PME load: %.3f\n",ratio);
  }

  return ratio;
}
