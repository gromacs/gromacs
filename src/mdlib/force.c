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
 * GRowing Old MAkes el Chrono Sweat
 */
static char *SRCID_force_c = "$Id$";

#include <math.h>
#include "sysstuff.h"
#include "typedefs.h"
#include "macros.h"
#include "smalloc.h"
#include "assert.h"
#include "physics.h"
#include "macros.h"
#include "vec.h"
#include "force.h"
#include "invblock.h"
#include "confio.h"
#include "nsb.h"
#include "names.h"
#include "network.h"
#include "wnblist.h"
#include "pbc.h"
#include "ns.h"
#include "nrnb.h"
#include "bondf.h"
#include "mshift.h"
#include "txtdump.h"
#include "ewald_util.h"
#include "shift_util.h"
#include "pppm.h"
#include "poisson.h"
#include "ewald.h"
#include "pme.h"
#include "copyrite.h"

t_forcerec *mk_forcerec(void)
{
  t_forcerec *fr;
  
  snew(fr,1);
  
  return fr;
}

#ifdef DEBUG
static void pr_nbfp(FILE *fp,real *nbfp,bool bBHAM,int atnr)
{
  int i,j;
  
  for(i=0; (i<atnr); i++) {
    for(j=0; (j<atnr); j++) {
      fprintf(fp,"%2d - %2d",i,j);
      if (bBHAM)
	fprintf(fp,"  a=%10g, b=%10g, c=%10g\n",BHAMA(nbfp,atnr,i,j),
		BHAMB(nbfp,atnr,i,j),BHAMC(nbfp,atnr,i,j));
      else
	fprintf(fp,"  c6=%10g, c12=%10g\n",C6(nbfp,atnr,i,j),
		C12(nbfp,atnr,i,j));
    }
  }
}
#endif

static real *mk_nbfp(t_idef *idef,bool bBHAM)
{
  real *nbfp;
  int  i,j,k,atnr;
  
  atnr=idef->atnr;
  if (bBHAM) {
    snew(nbfp,3*atnr*atnr);
    for(i=k=0; (i<atnr); i++) {
      for(j=0; (j<atnr); j++,k++) {
	BHAMA(nbfp,atnr,i,j) = idef->iparams[k].bham.a;
	BHAMB(nbfp,atnr,i,j) = idef->iparams[k].bham.b;
	BHAMC(nbfp,atnr,i,j) = idef->iparams[k].bham.c;
      }
    }
  }
  else {
    snew(nbfp,2*atnr*atnr);
    for(i=k=0; (i<atnr); i++) {
      for(j=0; (j<atnr); j++,k++) {
	C6(nbfp,atnr,i,j)   = idef->iparams[k].lj.c6;
	C12(nbfp,atnr,i,j)  = idef->iparams[k].lj.c12;
      }
    }
  }
  return nbfp;
}

static void calc_rffac(FILE *log,int eel,real eps,real Rc,real Temp,
		       real zsq,matrix box,
		       real *kappa,real *epsfac,real *krf,real *crf)
{
  /* Compute constants for Generalized reaction field */
  static bool bFirst=TRUE;
  real   k1,k2,I,vol,rmin;
  
  if ((eel == eelRF) || (eel == eelGRF)) {
    vol     = det(box);
    I       = zsq/vol;
    if (eel == eelGRF) {
      /* Consistency check */
      if (Temp <= 0.0)
	fatal_error(0,"Temperature is %f while using"
		    " Generalized Reaction Field\n",Temp);
      /* Ionic strength (only needed for eelGRF */
      *kappa  = sqrt(2*I/(EPSILON0*eps*BOLTZ*Temp));
    }
    else
      *kappa = 0;

    /* eps == 0 signals infinite dielectric */
    if (eps == 0) {
      *krf = 1/(2*Rc*Rc*Rc);
      *crf = 0;
    }
    else {
      k1      = (1+*kappa*Rc);
      k2      = eps*sqr((real)(*kappa*Rc));
      
      *krf    = (((eps-1)*k1+k2)/((2*eps+1)*k1+2*k2)/(Rc*Rc*Rc));
      *crf    = 1/Rc + *krf*Rc*Rc;
    }
    *epsfac = ONE_4PI_EPS0;
    rmin    = pow(*krf*2.0,-1.0/3.0);
    
    if (bFirst) {
      if (eel == eelGRF)
	please_cite(log,"Tironi95a");
      fprintf(log,"%s:\n"
	      "epsRF = %10g, I   = %10g, volume = %10g, kappa  = %10g\n"
	      "rc    = %10g, krf = %10g, crf    = %10g, epsfac = %10g\n",
	      eel_names[eel],eps,I,vol,*kappa,Rc,*krf,*crf,*epsfac);
      fprintf(log,
	      "The electrostatics potential has its minimum at rc = %g\n",
	      rmin);
      
      bFirst=FALSE;
    }
  }
  else {
    /* If we're not using a reaction field, set the factor to 0
     * and multiply the dielectric constant by 1/eps
     */
    *kappa  = 0.0;
    *krf    = 0.0;
    *crf    = 0.0;
    *epsfac = ONE_4PI_EPS0/eps;
  } 
}

void update_forcerec(FILE *log,t_forcerec *fr,matrix box)
{
  calc_rffac(log,fr->eeltype,
	     fr->epsilon_r,fr->rcoulomb,fr->temp,fr->zsquare,box,
	     &fr->kappa,&fr->epsfac,&fr->k_rf,&fr->c_rf);
}

static double calc_avcsix(FILE *log,real *nbfp,int atnr,
			  int natoms,int type[],bool bBHAM)
{
  int    i,j,tpi,tpj;
  double csix;

  /* Check this code: do we really need a double loop? */  
  csix = 0;
  for(i=0; (i<natoms); i++) {
    tpi = type[i];
#ifdef DEBUG
    if (tpi >= atnr)
      fatal_error(0,"Atomtype[%d] = %d, maximum = %d",i,tpi,atnr);
#endif
    for(j=0; (j<natoms); j++) {
      tpj   = type[j];
#ifdef DEBUG
      if (tpj >= atnr)
	fatal_error(0,"Atomtype[%d] = %d, maximum = %d",j,tpj,atnr);
#endif
      if (bBHAM)
	csix += BHAMC(nbfp,atnr,tpi,tpj);
      else
	csix += C6(nbfp,atnr,tpi,tpj);
    }
  }
  csix /= (natoms*natoms);
  if (debug)
    fprintf(debug,"Average C6 parameter is: %10g\n",csix);
  
  return csix;
}

void set_avcsix(FILE *log,t_forcerec *fr,t_mdatoms *mdatoms)
{
  fr->avcsix=calc_avcsix(log,fr->nbfp,fr->ntype,mdatoms->nr,
			 mdatoms->typeA,fr->bBHAM);
}

static void set_bham_b_max(FILE *log,t_forcerec *fr,t_mdatoms *mdatoms)
{
  int  i,j,tpi,tpj,ntypes,natoms,*type;
  real b,bmin;
  real *nbfp;

  fprintf(log,"Determining largest Buckingham b parameter for table\n");
  nbfp   = fr->nbfp;
  ntypes = fr->ntype;
  type   = mdatoms->typeA;
  natoms = mdatoms->nr;

  bmin           = -1;
  fr->bham_b_max = 0;
  for(i=0; (i<natoms); i++) {
    tpi = type[i];
    if (tpi >= ntypes)
      fatal_error(0,"Atomtype[%d] = %d, maximum = %d",i,tpi,ntypes);
    
    for(j=0; (j<natoms); j++) {
      tpj   = type[j];
      if (tpj >= ntypes)
	fatal_error(0,"Atomtype[%d] = %d, maximum = %d",j,tpj,ntypes);
      b = BHAMB(nbfp,ntypes,tpi,tpj);
      if (b > fr->bham_b_max)
	fr->bham_b_max = b;
      if ((b < bmin) || (bmin==-1))
	bmin = b;
    }
  }
  fprintf(log,"Buckingham b parameters, min: %g, max: %g\n",
	  bmin,fr->bham_b_max);
}

void init_forcerec(FILE *fp,
		   t_forcerec *fr,
		   t_inputrec *ir,
		   t_block    *mols,
		   t_commrec  *cr,
		   t_block    *cgs,
		   t_idef     *idef,
		   t_mdatoms  *mdatoms,
		   t_nsborder *nsb,
		   matrix     box,
		   bool bMolEpot)
{
  int    i,j,m,natoms,ngrp;
  real   q,zsq,nrdf,T;
  rvec   box_size;
  double rtab;
  
  natoms         = mdatoms->nr;

  /* Free energy */
  fr->efep       = ir->efep;
  fr->sc_alpha   = ir->sc_alpha;

  /* Neighbour searching stuff */
  fr->bGrid      = (ir->ns_type == ensGRID);
  fr->ndelta     = ir->ndelta;
  fr->eBox       = ir->eBox;
  fr->rlist      = ir->rlist;
  fr->rlistlong  = max(ir->rlist,max(ir->rcoulomb,ir->rvdw));
  fr->eeltype    = ir->coulombtype;
  fr->vdwtype    = ir->vdwtype;
  fr->bTwinRange = (fr->rlistlong > fr->rlist);
  fr->bEwald     = ((fr->eeltype == eelPME) || (fr->eeltype == eelEWALD));
  fr->bTab       = ((fr->eeltype != eelCUT) || (fr->vdwtype != evdwCUT));
  fr->bRF        = (((fr->eeltype == eelRF) || (fr->eeltype == eelGRF)) &&
		    (fr->vdwtype == evdwCUT));
  if ((fr->bRF) && (fr->vdwtype == evdwCUT))
    fr->bTab = FALSE;
  if (fp)
    fprintf(fp,"Table routines are used: %s\n",bool_names[fr->bTab]);
  
  if(fr->bEwald)
    fr->ewaldcoeff=calc_ewaldcoeff(ir->rcoulomb, ir->ewald_rtol);
  /* Tables are used for direct ewald sum */

  /* Domain decomposition parallellism... */
  fr->bDomDecomp = ir->bDomDecomp;
  fr->Dimension  = ir->decomp_dir;
  
  /* Electrostatics */
  fr->epsilon_r  = ir->epsilon_r;
  fr->fudgeQQ    = ir->fudgeQQ;
  fr->rcoulomb_switch = ir->rcoulomb_switch;
  fr->rcoulomb        = ir->rcoulomb;
  
  /* Must really support table functions with solvent_opt */
  fr->nWater     = ir->solvent_opt;
  fr->bWaterOpt  = (fr->nWater >= 0);
  
  fr->nWatMol = 0;
  if (fr->bWaterOpt) {
    for(i=START(nsb); (i<START(nsb)+HOMENR(nsb)); i++)
      if (fr->nWater == mdatoms->typeA[i])
	fr->nWatMol++;
    if (fp)
      fprintf(fp,"There are %d water molecules on cpu %d\n",
	      fr->nWatMol,nsb->pid);
  }
  
  /* Parameters for generalized RF */
  fr->zsquare = 0.0;
  fr->temp    = 0.0;
  
  if (fr->eeltype == eelGRF) {
    zsq = 0.0;
    for (i=0; (i<cgs->nr); i++) {
      q = 0;
      for(j=cgs->index[i]; (j<cgs->index[i+1]); j++)
	q+=mdatoms->chargeA[cgs->a[j]];
      if (q != 0.0)
	/* Changed from square to fabs 990314 DvdS 
	 * Does not make a difference for monovalent ions, but doe for 
	 * divalent ions (Ca2+!!)
	 */
	zsq += fabs(q);
    }
    fr->zsquare = zsq;
    
    T    = 0.0;
    nrdf = 0.0;
    for(i=0; (i<ir->opts.ngtc); i++) {
      nrdf += ir->opts.nrdf[i];
      T    += (ir->opts.nrdf[i] * ir->opts.ref_t[i]);
    }
    if (nrdf == 0) 
      fatal_error(0,"No degrees of freedom!");
    fr->temp   = T/nrdf;
  }
  else if (EEL_LR(fr->eeltype) || (fr->eeltype == eelSHIFT) || 
	   (fr->eeltype == eelUSER) || (fr->eeltype == eelSWITCH)) {
    /* We must use the long range cut-off for neighboursearching...
     * An extra range of e.g. 0.1 nm (half the size of a charge group)
     * is necessary for neighboursearching. This allows diffusion 
     * into the cut-off range (between neighborlist updates), 
     * and gives more accurate forces because all atoms within the short-range
     * cut-off rc must be taken into account, while the ns criterium takes
     * only those with the center of geometry within the cut-off.
     * (therefore we have to add half the size of a charge group, plus
     * something to account for diffusion if we have nstlist > 1)
     */
    for(m=0; (m<DIM); m++)
      box_size[m]=box[m][m];

    if (fr->phi == NULL)
      snew(fr->phi,mdatoms->nr);
    
    if ((fr->eeltype==eelPPPM) || (fr->eeltype==eelPOISSON) || 
	(fr->eeltype == eelSHIFT && fr->rcoulomb > fr->rcoulomb_switch))
	set_shift_consts(fp,fr->rcoulomb_switch,fr->rcoulomb,box_size,fr);
  }

  /* Initiate arrays */
  if (fr->bTwinRange || (EEL_LR(fr->eeltype))) {
    snew(fr->flr,natoms);
    snew(fr->fshift_lr,SHIFTS);
  }
  /* Mask that says whether or not this NBF list should be computed */
  /*  if (fr->bMask == NULL) {
    ngrp = ir->opts.ngener*ir->opts.ngener;
    snew(fr->bMask,ngrp);*/
    /* Defaults to always */
  /*    for(i=0; (i<ngrp); i++)
      fr->bMask[i] = TRUE;
      }*/

  if (fr->cg_cm == NULL)
    snew(fr->cg_cm,cgs->nr);
  if (fr->shift_vec == NULL)
    snew(fr->shift_vec,SHIFTS);
    
  if (fr->fshift == NULL)
    snew(fr->fshift,SHIFTS);
  
  if (bMolEpot && (fr->nmol==0)) {
    fr->nmol=mols->nr;
    fr->mol_nr=make_invblock(mols,natoms);
    snew(fr->mol_epot,fr->nmol);
    fr->nstcalc=ir->nstenergy;
  }
  
  if (fr->nbfp == NULL) {
    fr->ntype = idef->atnr;
    fr->bBHAM = (idef->functype[0] == F_BHAM);
    fr->nbfp  = mk_nbfp(idef,fr->bBHAM);
  }
  /* Copy the energy group exclusions */
  fr->eg_excl = ir->opts.eg_excl;

  /* Van der Waals stuff */
  fr->rvdw        = ir->rvdw;
  fr->rvdw_switch = ir->rvdw_switch;
  if ((fr->vdwtype != evdwCUT) && !fr->bBHAM) {
    if (fr->rvdw_switch >= fr->rvdw)
      fatal_error(0,"rvdw_switch (%g) must be < rvdw (%g)",
		  fr->rvdw_switch,fr->rvdw);
    if (fp)
      fprintf(fp,"Using %s Lennard-Jones, switch between %g and %g nm\n",
	      (fr->eeltype==eelSWITCH) ? "switched":"shifted",
	      fr->rvdw_switch,fr->rvdw);
  } 

  if (fp)
    fprintf(fp,"Cut-off's:   NS: %g   Coulomb: %g   %s: %g\n",
	    fr->rlist,fr->rcoulomb,fr->bBHAM ? "BHAM":"LJ",fr->rvdw);
  
  if (ir->bDispCorr)
    set_avcsix(fp,fr,mdatoms);
  if (fr->bBHAM)
    set_bham_b_max(fp,fr,mdatoms);
  
  /* Now update the rest of the vars */
  update_forcerec(fp,fr,box);
  /* if we are using LR electrostatics, and they are tabulated,
   * the tables will contain shifted coulomb interactions.
   * Since we want to use the non-shifted ones for 1-4
   * coulombic interactions, we must have an extra set of
   * tables. This should be done in tables.c, instead of this
   * ugly hack, but it works for now...
   */

#define MAX_14_DIST 1.0
  /* Shell to account for the maximum chargegroup radius (2*0.2 nm) *
   * and diffusion during nstlist steps (0.2 nm)                    */
#define TAB_EXT 0.6

  if (fr->bTab) {
    if (EEL_LR(fr->eeltype)) {
      /* generate extra tables for 1-4 interactions only
       * fake the forcerec so make_tables thinks it should
       * just create the non shifted version 
       */
      fr->bTab=FALSE;
      fr->rtab=MAX_14_DIST;
      make_tables(fp,fr,MASTER(cr));
      fr->bTab=TRUE;
      fr->VFtab14=fr->VFtab;
      fr->VFtab=NULL;
    }
    fr->rtab = max(fr->rlistlong+TAB_EXT,MAX_14_DIST);
  } 
  else if (fr->efep != efepNO) {
    if (fr->rlistlong == 0) {
      char *ptr,*envvar="FEP_TABLE_LENGTH";
      fr->rtab = 5;
      ptr = getenv(envvar);
      if (ptr) {
	sscanf(ptr,"%lf",&rtab);
	fr->rtab = rtab;
      }
      if (fp)
	fprintf(fp,"\nNote: Setting the free energy table length to %g nm\n"
		"      You can set this value with the environment variable %s"
		"\n\n",fr->rtab,envvar);
    } 
    else
      fr->rtab = max(fr->rlistlong+TAB_EXT,MAX_14_DIST);
  } 
  else
    fr->rtab = MAX_14_DIST;
  
  /* make tables for ordinary interactions */
  make_tables(fp,fr,MASTER(cr));
  if(!(EEL_LR(fr->eeltype) && fr->bTab))
    fr->VFtab14=fr->VFtab;
}

#define pr_real(fp,r) fprintf(fp,"%s: %e\n",#r,r)
#define pr_int(fp,i)  fprintf((fp),"%s: %d\n",#i,i)
#define pr_bool(fp,b) fprintf((fp),"%s: %s\n",#b,bool_names[b])

void pr_forcerec(FILE *fp,t_forcerec *fr,t_commrec *cr)
{
  pr_real(fp,fr->rlist);
  pr_real(fp,fr->rcoulomb);
  pr_real(fp,fr->fudgeQQ);
  pr_int(fp,fr->ndelta);
  pr_bool(fp,fr->bGrid);
  pr_bool(fp,fr->bTwinRange);
  /*pr_int(fp,fr->cg0);
    pr_int(fp,fr->hcg);*/
  pr_int(fp,fr->ntab);
  if (fr->ntab > 0) {
    pr_real(fp,fr->rcoulomb_switch);
    pr_real(fp,fr->rcoulomb);
  }
  
  pr_int(fp,fr->nmol);
  pr_int(fp,fr->nstcalc);
  
  fflush(fp);
}

void ns(FILE *fp,
	t_forcerec *fr,
	rvec       x[],
	rvec       f[],
	matrix     box,
	t_groups   *grps,
	t_grpopts  *opts,
	t_topology *top,
	t_mdatoms  *md,
	t_commrec  *cr,
	t_nrnb     *nrnb,
	t_nsborder *nsb,
	int        step,
	real       lambda,
	real       *dvdlambda)
{
  static bool bFirst=TRUE;
  static int  nDNL;
  char   *ptr;
  int    nsearch;
  
  if (bFirst) {
    ptr=getenv("DUMPNL");
    if (ptr) {
      nDNL=atoi(ptr);
      fprintf(fp,"nDNL = %d\n",nDNL);  
    } else
      nDNL=0;
    /* Allocate memory for the neighbor lists */
    init_neighbor_list(fp,fr,HOMENR(nsb));
      
    bFirst=FALSE;
  }
    
  if (fr->bTwinRange) 
    fr->nlr=0;

  /* Whether or not we do dynamic load balancing,
   * workload contains the proper numbers of charge groups
   * to be searched.
   */
  if (cr->pid == 0)
    fr->cg0=0;
  else
    fr->cg0=nsb->workload[cr->pid-1];
  fr->hcg=nsb->workload[cr->pid];

  nsearch = search_neighbours(fp,fr,x,box,top,grps,cr,nsb,nrnb,md,
			      lambda,dvdlambda);
  if (debug)
    fprintf(debug,"nsearch = %d\n",nsearch);
    
  /* Check whether we have to do dynamic load balancing */
  /*if ((nsb->nstDlb > 0) && (mod(step,nsb->nstDlb) == 0))
    count_nb(cr,nsb,&(top->blocks[ebCGS]),nns,fr->nlr,
    &(top->idef),opts->ngener);
  */
  if (nDNL > 0)
    dump_nblist(fp,fr,nDNL);
}

void force(FILE       *fp,     int        step,
	   t_forcerec *fr,      t_inputrec *ir,
	   t_idef     *idef,    t_nsborder *nsb,
	   t_commrec  *cr,      t_nrnb     *nrnb,
	   t_groups   *grps,    t_mdatoms  *md,
	   int        ngener,   t_grpopts  *opts,
	   rvec       x[],      rvec       f[],
	   real       epot[], 
	   bool       bVerbose, matrix     box,
	   real       lambda,   t_graph    *graph,
	   t_block    *excl,    bool       bNBFonly,
	   matrix lr_vir)
{
  int     i,nit;
  bool    bDoEpot;
  rvec    box_size;
  real    Vlr,Vself=0;
  
  /* Reset box */
  for(i=0; (i<DIM); i++)
    box_size[i]=box[i][i];
    
  bDoEpot=((fr->nmol > 0) && (fr->nstcalc > 0) && (mod(step,fr->nstcalc)==0));
  
  /* Reset epot... */
  if (bDoEpot) 
    for(i=0; (i<fr->nmol); i++)
      fr->mol_epot[i]=0.0;
  debug_gmx();
  
  /* Call the short range functions all in one go. */
  do_fnbf(fp,fr,x,f,md,
	  fr->bBHAM ? grps->estat.ee[egBHAM] : grps->estat.ee[egLJ],
	  grps->estat.ee[egCOUL],box_size,nrnb,
	  lambda,&epot[F_DVDL],FALSE,-1);
  debug_gmx();

  if (debug) 
    pr_rvecs(debug,0,"fshift after SR",fr->fshift,SHIFTS);
  
  /* Shift the coordinates. Must be done before bonded forces and PPPM, 
   * but is also necessary for SHAKE and update, therefore it can NOT 
   * go when no bonded forces have to be evaluated.
   */
  if (debug && 0)
    p_graph(debug,"DeBUGGGG",graph);
  
  /* Check whether we need to do bondeds */
  if (!bNBFonly) {
    shift_self(graph,fr->shift_vec,x);
    if (debug && 0) {
      fprintf(debug,"BBBBBBBBBBBBBBBB\n");
      fprintf(debug,"%5d\n",graph->nnodes);
      for(i=graph->start; (i<=graph->end); i++)
	fprintf(debug,"%5d%5s%5s%5d%8.3f%8.3f%8.3f\n",
		i,"A","B",i,x[i][XX],x[i][YY],x[i][ZZ]);
      fprintf(debug,"%10.5f%10.5f%10.5f\n",
	      box[XX][XX],box[YY][YY],box[ZZ][ZZ]);
    }
    inc_nrnb(nrnb,eNR_SHIFTX,graph->nnodes);
    debug_gmx();
  }
  
  if (EEL_LR(fr->eeltype)) {
    switch (fr->eeltype) {
    case eelPPPM:
      Vlr = do_pppm(fp,FALSE,x,fr->flr,md->chargeA,
		    box_size,fr->phi,cr,nsb,nrnb);
      break;
    case eelPOISSON:
      Vlr = do_poisson(fp,FALSE,ir,md->nr,x,fr->flr,md->chargeA,
		       box_size,fr->phi,cr,nrnb,&nit,TRUE);
      break;
    case eelPME:
      Vlr = do_pme(fp,FALSE,ir,x,fr->flr,md->chargeA,
		   box_size,cr,nsb,nrnb,lr_vir,fr->ewaldcoeff);
      break;
    case eelEWALD:
      Vlr = do_ewald(fp,FALSE,ir,x,fr->flr,md->chargeA,
		     box_size,cr,nsb,lr_vir,fr->ewaldcoeff);
      break;
    default:
      Vlr = 0;
      fatal_error(0,"No such electrostatics method implemented %s",
		  eel_names[fr->eeltype]);
    }
    if(fr->bEwald)
	Vself = ewald_LRcorrection(fp,nsb,cr,fr,
				   md->chargeA,excl,x,box_size,lr_vir);
    else
	Vself = shift_LRcorrection(fp,nsb,cr,fr,
			       md->chargeA,excl,x,TRUE,box_size,lr_vir);
    epot[F_LR] = Vlr - Vself;
    if (debug)
      fprintf(debug,"Vlr = %g, Vself = %g, Vlr_corr = %g\n",
	      Vlr,Vself,epot[F_LR]);
    if (debug) {
      pr_rvecs(debug,0,"lr_vir after corr",lr_vir,DIM);
      pr_rvecs(debug,0,"fshift after LR Corrections",fr->fshift,SHIFTS);
    }
  }
  debug_gmx();
  
  if (debug)    
    print_nrnb(debug,nrnb);
  debug_gmx();

  if (!bNBFonly) {
    calc_bonds(fp,idef,x,f,fr,graph,epot,nrnb,box,lambda,md,
	       opts->ngener,grps->estat.ee[egLJ14],grps->estat.ee[egCOUL14]);
    debug_gmx();
  }
  if (debug) 
    pr_rvecs(debug,0,"fshift after bondeds",fr->fshift,SHIFTS);
  
  for(i=0; (i<F_EPOT); i++)
    if (i != F_DISRES)
      epot[F_EPOT]+=epot[i];
}
