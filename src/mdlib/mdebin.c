/*
 * $Id$
 * 
 *                This source code is part of
 * 
 *                 G   R   O   M   A   C   S
 * 
 *          GROningen MAchine for Chemical Simulations
 * 
 *                        VERSION 3.2.0
 * Written by David van der Spoel, Erik Lindahl, Berk Hess, and others.
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2004, The GROMACS development team,
 * check out http://www.gromacs.org for more information.

 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 * 
 * If you want to redistribute modifications, please consider that
 * scientific software is very special. Version control is crucial -
 * bugs must be traceable. We will be happy to consider code for
 * inclusion in the official distribution, but derived work must not
 * be called official GROMACS. Details are found in the README & COPYING
 * files - if they are missing, get the official version at www.gromacs.org.
 * 
 * To help us fund GROMACS development, we humbly ask that you cite
 * the papers on the package - you can find them in the top README file.
 * 
 * For more info, check our website at http://www.gromacs.org
 * 
 * And Hey:
 * GROwing Monsters And Cloning Shrimps
 */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <string.h>
#include "typedefs.h"
#include "string2.h"
#include "mdebin.h"
#include "assert.h"
#include "smalloc.h"
#include "physics.h"
#include "enxio.h"
#include "vec.h"
#include "disre.h"
#include "main.h"
#include "network.h"
#include "names.h"
#include "orires.h"

static bool bEInd[egNR] = { TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE };

static bool bEner[F_NRE];
static char *boxs_nm[] = {
  "Box-X", "Box-Y", "Box-Z","Volume","Density (SI)",
  "pV"
};

static char *tricl_boxs_nm[] = {
  "Box-XX", "Box-YX", "Box-YY", "Box-ZX", "Box-ZY", "Box-ZZ",
  "Volume", "Density (SI)", "pV"
};

static char *boxvel_nm[] = {
  "Box-Vel-XX", "Box-Vel-YY", "Box-Vel-ZZ",
  "Box-Vel-YX", "Box-Vel-ZX", "Box-Vel-ZY"
};

static char *pcouplmu_nm[] = {
  "Pcoupl-Mu-XX", "Pcoupl-Mu-YY", "Pcoupl-Mu-ZZ",
  "Pcoupl-Mu-YX", "Pcoupl-Mu-ZX", "Pcoupl-Mu-ZY"
};

#define NBOXS asize(boxs_nm)
#define NTRICLBOXS asize(tricl_boxs_nm)

static bool bShake,bTricl,bDynBox;
static int  f_nre=0,epc,etc;

t_mdebin *init_mdebin(int fp_ene,const t_groups *grps,const t_atoms *atoms,
		      const t_idef *idef,const t_inputrec *ir,
		      t_commrec *cr)
{
  char *ener_nm[F_NRE];
  static char *vir_nm[] = {
    "Vir-XX", "Vir-XY", "Vir-XZ",
    "Vir-YX", "Vir-YY", "Vir-YZ",
    "Vir-ZX", "Vir-ZY", "Vir-ZZ"
  };
  static char *sv_nm[] = {
    "ShakeVir-XX", "ShakeVir-XY", "ShakeVir-XZ",
    "ShakeVir-YX", "ShakeVir-YY", "ShakeVir-YZ",
    "ShakeVir-ZX", "ShakeVir-ZY", "ShakeVir-ZZ"
  };
  static char *fv_nm[] = {
    "ForceVir-XX", "ForceVir-XY", "ForceVir-XZ",
    "ForceVir-YX", "ForceVir-YY", "ForceVir-YZ",
    "ForceVir-ZX", "ForceVir-ZY", "ForceVir-ZZ"
  };
  static char *pres_nm[] = {
    "Pres-XX (bar)","Pres-XY (bar)","Pres-XZ (bar)",
    "Pres-YX (bar)","Pres-YY (bar)","Pres-YZ (bar)",
    "Pres-ZX (bar)","Pres-ZY (bar)","Pres-ZZ (bar)"
  };
  static char *surft_nm[] = {
    "#Surf*SurfTen"
  };
  static char *mu_nm[] = {
    "Mu-X", "Mu-Y", "Mu-Z"
  };
  static char *vcos_nm[] = {
    "2CosZ*Vel-X"
  };
  static char *visc_nm[] = {
    "1/Viscosity (SI)"
  };
  static   char   **grpnms;
  char     **gnm;
  char     buf[256];
  t_mdebin *md;
  int      i,j,ni,nj,n,k,kk;
  bool     bBHAM,b14;
  
  bBHAM = (idef->functype[0] == F_BHAM);
  b14   = (idef->il[F_LJ14].nr > 0);

  for(i=0; i<F_NRE; i++) {
    bEner[i] = FALSE;
    if (i == F_LJ)
      bEner[i] = !bBHAM;
    else if (i == F_BHAM)
      bEner[i] = bBHAM;
    else if (i == F_COUL_LR)
      bEner[i] = (ir->rcoulomb > ir->rlist);
    else if (i == F_LJ_LR)
      bEner[i] = (!bBHAM && ir->rvdw > ir->rlist);
    else if (i == F_BHAM_LR)
      bEner[i] = (bBHAM && ir->rvdw > ir->rlist);
    else if (i == F_RF_EXCL)
      bEner[i] = (EEL_RF(ir->coulombtype) && ir->coulombtype != eelRF_OLD);
    else if (i == F_COUL_RECIP)
      bEner[i] = EEL_FULL(ir->coulombtype);
    else if (i == F_LJ14)
      bEner[i] = b14;
    else if (i == F_COUL14)
      bEner[i] = b14;
    else if ((i == F_DVDL) || (i == F_DVDLKIN))
      bEner[i] = (ir->efep != efepNO);
    else if ((strstr(interaction_function[i].name,"DUM") != NULL) ||
	     (i == F_SHAKE) || (i == F_SETTLE))
      bEner[i] = FALSE;
    else if ((i == F_COUL_SR) || (i == F_EPOT) || (i == F_ETOT) ||
	     (i == F_EKIN) || (i == F_TEMP) || (i == F_PRES))
      bEner[i] = TRUE;
    else if (i == F_DISPCORR)
      bEner[i] = (ir->eDispCorr != edispcNO);
    else if (i == F_DISRESVIOL)
      bEner[i] = (idef->il[F_DISRES].nr > 0);
    else if (i == F_ORIRESDEV)
      bEner[i] = (idef->il[F_ORIRES].nr > 0);
    else if (i == F_CONNBONDS)
      bEner[i] = FALSE;
    else
      bEner[i] = (idef->il[i].nr > 0);
  }
  if (PAR(cr))
    gmx_sumi(F_NRE,bEner,cr);

  for(i=0; i<F_NRE; i++)
    if (bEner[i]) {
      ener_nm[f_nre]=interaction_function[i].longname;
      f_nre++;
    }

  bShake = (idef->il[F_SHAKE].nr > 0) || (idef->il[F_SETTLE].nr > 0);
  if (bShake) 
    bShake = (getenv("SHAKEVIR") != NULL);
  epc = ir->epc;
  bTricl = TRICLINIC(ir->compress) || TRICLINIC(ir->deform);
  bDynBox = DYNAMIC_BOX(*ir);
  etc = ir->etc;
  
  /* Energy monitoring */
  snew(md,1);
  md->ebin  = mk_ebin();
  md->ie    = get_ebin_space(md->ebin,f_nre,ener_nm);
  if (bDynBox)
    md->ib    = get_ebin_space(md->ebin, bTricl ? NTRICLBOXS :
			       NBOXS, bTricl ? tricl_boxs_nm : boxs_nm);
  if (bShake) {
    md->isvir = get_ebin_space(md->ebin,asize(sv_nm),sv_nm);
    md->ifvir = get_ebin_space(md->ebin,asize(fv_nm),fv_nm);
  }
  md->ivir   = get_ebin_space(md->ebin,asize(vir_nm),vir_nm);
  md->ipres  = get_ebin_space(md->ebin,asize(pres_nm),pres_nm);
  md->isurft = get_ebin_space(md->ebin,asize(surft_nm),surft_nm);
  if (epc == epcPARRINELLORAHMAN) {
    md->ipc  = get_ebin_space(md->ebin,bTricl ? 6 : 3,boxvel_nm);
  } else if (epc == epcBERENDSEN || epc == epcISOTROPIC) {
    md->ipc = get_ebin_space(md->ebin,bTricl ? 6 : 3,pcouplmu_nm);
  }
  md->imu    = get_ebin_space(md->ebin,asize(mu_nm),mu_nm);
  if (grps->cosacc.cos_accel != 0) {
    md->ivcos = get_ebin_space(md->ebin,asize(vcos_nm),vcos_nm);
    md->ivisc = get_ebin_space(md->ebin,asize(visc_nm),visc_nm);
  }
  if (ir->rcoulomb > ir->rlist) 
    bEInd[egCOULLR] = TRUE;
  if (!bBHAM) {
    if (ir->rvdw > ir->rlist)
      bEInd[egLJLR]   = TRUE;
  } else {
    bEInd[egLJSR]   = FALSE;
    bEInd[egBHAMSR] = TRUE;
    if (ir->rvdw > ir->rlist)
      bEInd[egBHAMLR]   = TRUE;
  }
  if (b14) {
    bEInd[egLJ14] = TRUE;
    bEInd[egCOUL14] = TRUE;
  }
  md->nEc=0;
  for(i=0; (i<egNR); i++)
    if (bEInd[i])
      md->nEc++;
      
  n=atoms->grps[egcENER].nr;
  md->nEg=n;
  md->nE=(n*(n+1))/2;
  snew(md->igrp,md->nE);
  if (md->nE > 1) {
    n=0;
    snew(gnm,md->nEc);
    for(k=0; (k<md->nEc); k++)
      snew(gnm[k],STRLEN);
    for(i=0; (i<atoms->grps[egcENER].nr); i++) {
      ni=atoms->grps[egcENER].nm_ind[i];
      for(j=i; (j<atoms->grps[egcENER].nr); j++) {
	nj=atoms->grps[egcENER].nm_ind[j];
	for(k=kk=0; (k<egNR); k++) {
	  if (bEInd[k]) {
	    sprintf(gnm[kk],"%s:%s-%s",egrp_nm[k],
		    *(atoms->grpname[ni]),*(atoms->grpname[nj]));
	    kk++;
	  }
	}
	md->igrp[n]=get_ebin_space(md->ebin,md->nEc,gnm);
	n++;
      }
    }
    for(k=0; (k<md->nEc); k++)
      sfree(gnm[k]);
    sfree(gnm);
    
    assert(n==md->nE);
  }
  
  md->nTC=atoms->grps[egcTC].nr;
  snew(grpnms,md->nTC);
  for(i=0; (i<md->nTC); i++) {
    ni=atoms->grps[egcTC].nm_ind[i];
    sprintf(buf,"T-%s",*(atoms->grpname[ni]));
    grpnms[i]=strdup(buf);
  }
  md->itemp=get_ebin_space(md->ebin,md->nTC,grpnms);
  if (etc == etcNOSEHOOVER) {
    for(i=0; (i<md->nTC); i++) {
      ni=atoms->grps[egcTC].nm_ind[i];
      sprintf(buf,"Xi-%s",*(atoms->grpname[ni]));
      grpnms[i]=strdup(buf);
    }
    md->itc=get_ebin_space(md->ebin,md->nTC,grpnms);
  } else  if (etc == etcBERENDSEN || etc == etcYES) {
    for(i=0; (i<md->nTC); i++) {
      ni=atoms->grps[egcTC].nm_ind[i];
      sprintf(buf,"Lamb-%s",*(atoms->grpname[ni]));
      grpnms[i]=strdup(buf);
    }
    md->itc=get_ebin_space(md->ebin,md->nTC,grpnms);
  }
  sfree(grpnms);
  
  md->nU=atoms->grps[egcACC].nr;
  if (md->nU > 1) {
    snew(grpnms,3*md->nU);
    for(i=0; (i<md->nU); i++) {
      ni=atoms->grps[egcACC].nm_ind[i];
      sprintf(buf,"Ux-%s",*(atoms->grpname[ni]));
      grpnms[3*i+XX]=strdup(buf);
      sprintf(buf,"Uy-%s",*(atoms->grpname[ni]));
      grpnms[3*i+YY]=strdup(buf);
      sprintf(buf,"Uz-%s",*(atoms->grpname[ni]));
      grpnms[3*i+ZZ]=strdup(buf);
    }
    md->iu=get_ebin_space(md->ebin,3*md->nU,grpnms);
    sfree(grpnms);
  }
  
  if (fp_ene != -1)
    do_enxnms(fp_ene,&md->ebin->nener,&md->ebin->enm);
    
#ifdef DEBUG
  for(i=0; (i<md->ebin->nener); i++)
    fprintf(stdlog,"%5d  %20s\n",i,md->ebin->enm[i]);
#endif
  return md;
}

static void copy_energy(real e[],real ecpy[])
{
  int i,j;
  
  for(i=j=0; (i<F_NRE); i++)
    if (bEner[i])
      ecpy[j++] = e[i];
  assert(j == f_nre);
}

void upd_mdebin(t_mdebin *md,FILE *fp_dgdl,
		real tmass,int step,real time,
		real ener[],
		t_state *state,
		matrix  box,
		tensor svir,
		tensor fvir,
		tensor vir,
		tensor pres,
		t_groups *grps,
		rvec mu_tot)
{
  static real *ttt=NULL;
  static rvec *uuu=NULL;
  int    i,j,k,kk,m,n,gid;
  real   bs[NBOXS],tmp6[6];
  real   tricl_bs[NTRICLBOXS];
  real   eee[egNR];
  real   ecopy[F_NRE];
  real   tmp;
  
  /* Do NOT use the box in the state variable, but the separate box provided
   * as an argument. This is because we sometimes need to write the box from
   * the last timestep to match the trajectory frames.
   */
  copy_energy(ener,ecopy);
  add_ebin(md->ebin,md->ie,f_nre,ecopy,step);
  if (bDynBox || grps->cosacc.cos_accel != 0) {
    if(bTricl) {
      tricl_bs[0]=box[XX][XX];
      tricl_bs[1]=box[YY][XX];
      tricl_bs[2]=box[YY][YY];
      tricl_bs[3]=box[ZZ][XX];
      tricl_bs[4]=box[ZZ][YY];
      tricl_bs[5]=box[ZZ][ZZ];
      /* This is the volume */
      tricl_bs[6]=tricl_bs[0]*tricl_bs[2]*tricl_bs[5];
      /* This is the density */
      tricl_bs[7] = (tmass*AMU)/(tricl_bs[6]*NANO*NANO*NANO);
    } else {
      for(m=0; (m<DIM); m++) 
	bs[m]=box[m][m];
      /* This is the volume */
      bs[3] = bs[XX]*bs[YY]*bs[ZZ];      
      /* This is the density */
      bs[4] = (tmass*AMU)/(bs[3]*NANO*NANO*NANO);
    }
  }  
  if (bDynBox) {
    /* This is pV (in kJ/mol) */  
    if(bTricl) {
      tricl_bs[8] = tricl_bs[6]*ener[F_PRES]/PRESFAC;
      add_ebin(md->ebin,md->ib,NTRICLBOXS,tricl_bs,step);
    } else {
      bs[5] = bs[3]*ener[F_PRES]/PRESFAC;
      add_ebin(md->ebin,md->ib,NBOXS,bs,step);
    }
  }  
  if (bShake) {
    add_ebin(md->ebin,md->isvir,9,svir[0],step);
    add_ebin(md->ebin,md->ifvir,9,fvir[0],step);
  }
  add_ebin(md->ebin,md->ivir,9,vir[0],step);
  add_ebin(md->ebin,md->ipres,9,pres[0],step);
  tmp = (pres[ZZ][ZZ]-(pres[XX][XX]+pres[YY][YY])*0.5)*box[ZZ][ZZ];
  add_ebin(md->ebin,md->isurft,1,&tmp,step);
  if (epc == epcPARRINELLORAHMAN) {
    tmp6[0] = state->boxv[XX][XX];
    tmp6[1] = state->boxv[YY][YY];
    tmp6[2] = state->boxv[ZZ][ZZ];
    tmp6[3] = state->boxv[YY][XX];
    tmp6[4] = state->boxv[ZZ][XX];
    tmp6[5] = state->boxv[ZZ][YY];
    add_ebin(md->ebin,md->ipc,bTricl ? 6 : 3,tmp6,step);
  } else if (epc == epcBERENDSEN || epc == epcISOTROPIC) {
    tmp6[0] = state->pcoupl_mu[XX][XX];
    tmp6[1] = state->pcoupl_mu[YY][YY];
    tmp6[2] = state->pcoupl_mu[ZZ][ZZ];
    tmp6[3] = state->pcoupl_mu[YY][XX];
    tmp6[4] = state->pcoupl_mu[ZZ][XX];
    tmp6[5] = state->pcoupl_mu[ZZ][YY];
    add_ebin(md->ebin,md->ipc,bTricl ? 6 : 3,tmp6,step);
  }
  add_ebin(md->ebin,md->imu,3,mu_tot,step);

  if (grps->cosacc.cos_accel != 0) {
    add_ebin(md->ebin,md->ivcos,1,&(grps->cosacc.vcos),step);
    /* 1/viscosity, unit 1/(kg m^-1 s^-1) */
    if(bTricl) 
      tmp = 1/(grps->cosacc.cos_accel/(grps->cosacc.vcos*PICO)
	       *tricl_bs[7]*sqr(box[ZZ][ZZ]*NANO/(2*M_PI)));
    else 
      tmp = 1/(grps->cosacc.cos_accel/(grps->cosacc.vcos*PICO)
	       *bs[4]*sqr(box[ZZ][ZZ]*NANO/(2*M_PI)));
    add_ebin(md->ebin,md->ivisc,1,&tmp,step);    
  }
  if (md->nE > 1) {
    n=0;
    for(i=0; (i<md->nEg); i++) {
      for(j=i; (j<md->nEg); j++) {
	gid=GID(i,j,md->nEg);
	for(k=kk=0; (k<egNR); k++) 
	  if (bEInd[k])
	    eee[kk++]=grps->estat.ee[k][gid];
	add_ebin(md->ebin,md->igrp[n],md->nEc,eee,step);
	n++;
      }
    }
  }
  
  if(ttt == NULL) 
    snew(ttt,md->nTC);
  for(i=0; (i<md->nTC); i++)
    ttt[i]   = grps->tcstat[i].T;
  add_ebin(md->ebin,md->itemp,md->nTC,ttt,step);
  if (etc == etcNOSEHOOVER) {
    for(i=0; (i<md->nTC); i++)
      ttt[i] = state->nosehoover_xi[i];
    add_ebin(md->ebin,md->itc,md->nTC,ttt,step);
  } else if (etc == etcBERENDSEN || etc == etcYES) {
    for(i=0; (i<md->nTC); i++)
      ttt[i] = state->tcoupl_lambda[i];
    add_ebin(md->ebin,md->itc,md->nTC,ttt,step);
  }
  
  if (md->nU > 1) {
    if (uuu == NULL)
      snew(uuu,md->nU);
    for(i=0; (i<md->nU); i++)
      copy_rvec(grps->grpstat[i].u,uuu[i]);
    add_ebin(md->ebin,md->iu,3*md->nU,uuu[0],step);
  }
  if (fp_dgdl)
    fprintf(fp_dgdl,"%g %g\n",time,ener[F_DVDL]+ener[F_DVDLKIN]);
}

static void npr(FILE *log,int n,char c)
{
  for(; (n>0); n--) fprintf(log,"%c",c);
}

static void pprint(FILE *log,char *s)
{
  char   CHAR='#';
  int slen;

  slen=strlen(s);
  fprintf(log,"\t<======  "); 
  npr(log,slen,CHAR); 
  fprintf(log,"  ==>\n");
  fprintf(log,"\t<====  %s  ====>\n",s); 
  fprintf(log,"\t<==  "); 
  npr(log,slen,CHAR); 
  fprintf(log,"  ======>\n\n");
}

void print_ebin_header(FILE *log,int steps,real time,real lamb)
{
  fprintf(log,"   %12s   %12s   %12s\n"
	  "   %12d   %12.5f   %12.5f\n\n",
	  "Step","Time","Lambda",steps,time,lamb);
}

void print_ebin(int fp_ene,bool bEne,bool bDR,bool bOR,bool bDihR,
		FILE *log,int steps,real time,int mode,bool bCompact,
		t_mdebin *md,t_fcdata *fcd,t_atoms *atoms, t_grpopts *opts)
{
  static char **grpnms=NULL;
  static char *kjm="(kJ/mol)";
  char        buf[246];
  int         i,j,n,ni,nj,ndr,nor;
  int         nr[enxNR];
  real        *block[enxNR];
  t_enxframe  fr;
  
  switch (mode) {
  case eprNORMAL:
    fr.t          = time;
    fr.step       = steps;
    fr.nre        = (bEne) ? md->ebin->nener : 0;
    fr.ener       = md->ebin->e;
    fr.ndisre     = bDR ? fcd->disres.npr : 0;
    fr.rav        = fcd->disres.rav;
    fr.rt         = fcd->disres.rt;
    /* Optional additional blocks */
    for(i=0; i<enxNR; i++)
      nr[i] = 0;
    if (fcd->orires.nr > 0 && bOR) {
      diagonalize_orires_tensors(&(fcd->orires));
      nr[enxOR]     = fcd->orires.nr;
      block[enxOR]  = fcd->orires.otav;
      nr[enxORI]    = (fcd->orires.oinsl != fcd->orires.otav) ? 
	fcd->orires.nr : 0;
      block[enxORI] = fcd->orires.oinsl;
      nr[enxORT]    = fcd->orires.nex*12;
      block[enxORT] = fcd->orires.eig;
    }
    fr.nblock = 0;
    for(i=0; i<enxNR; i++)
      if (nr[i] > 0)
	fr.nblock = i+1;
    fr.nr         = nr;
    fr.block      = block;
    if (fr.nre || fr.ndisre || fr.nr[enxOR] || fr.nr[enxORI])
      do_enx(fp_ene,&fr);
    break;
  case eprAVER:
    if (log) pprint(log,"A V E R A G E S");
    break;
  case eprRMS:
    if (log) pprint(log,"R M S - F L U C T U A T I O N S");
    break;
  default:
    fatal_error(0,"Invalid print mode (%d)",mode);
  }
  
  if (log) {
    for(i=0;i<opts->ngtc;i++)
      if(opts->annealing[i]!=eannNO)
	fprintf(log,"Current ref_t for group %s: %8.1f\n",
		*(atoms->grpname[atoms->grps[egcTC].nm_ind[i]]),opts->ref_t[i]);
  
    if (mode==eprNORMAL && fcd->orires.nr>0)
      print_orires_log(log,&(fcd->orires));

    fprintf(log,"   Energies %s\n",kjm);
    pr_ebin(log,md->ebin,md->ie,f_nre,5,mode,steps,TRUE);  
    fprintf(log,"\n");

    if (!bCompact) {
      if (bDynBox) {
	pr_ebin(log,md->ebin,md->ib, bTricl ? NTRICLBOXS : NBOXS,5,mode,steps,TRUE);      
	fprintf(log,"\n");
      }
      if (bShake) {
	fprintf(log,"   Shake Virial %s\n",kjm);
	pr_ebin(log,md->ebin,md->isvir,9,3,mode,steps,FALSE);  
	fprintf(log,"\n");
	fprintf(log,"   Force Virial %s\n",kjm);
	pr_ebin(log,md->ebin,md->ifvir,9,3,mode,steps,FALSE);  
      fprintf(log,"\n");
      }
      fprintf(log,"   Total Virial %s\n",kjm);
      pr_ebin(log,md->ebin,md->ivir,9,3,mode,steps,FALSE);   
      fprintf(log,"\n");
      fprintf(log,"   Pressure (bar)\n");
      pr_ebin(log,md->ebin,md->ipres,9,3,mode,steps,FALSE);  
      fprintf(log,"\n");
      fprintf(log,"   Total Dipole (Debye)\n");
      pr_ebin(log,md->ebin,md->imu,3,3,mode,steps,FALSE);    
      fprintf(log,"\n");
      
      if (md->nE > 1) {
	if (grpnms==NULL) {
	  snew(grpnms,md->nE);
	  n=0;
	  for(i=0; (i<md->nEg); i++) {
	    ni=atoms->grps[egcENER].nm_ind[i];
	    for(j=i; (j<md->nEg); j++) {
	      nj=atoms->grps[egcENER].nm_ind[j];
	      sprintf(buf,"%s-%s",*(atoms->grpname[ni]),*(atoms->grpname[nj]));
	      grpnms[n++]=strdup(buf);
	    }
	  }
	}
	sprintf(buf,"Epot %s",kjm);
	fprintf(log,"%15s   ",buf);
	for(i=0; (i<egNR); i++) 
	  if (bEInd[i])
	    fprintf(log,"%12s   ",egrp_nm[i]);
	fprintf(log,"\n");
	for(i=0; (i<md->nE); i++) {
	  fprintf(log,"%15s",grpnms[i]);
	  pr_ebin(log,md->ebin,md->igrp[i],md->nEc,md->nEc,mode,steps,FALSE);
	}
	fprintf(log,"\n");
      }
      if (md->nTC > 1) {
	pr_ebin(log,md->ebin,md->itemp,md->nTC,4,mode,steps,TRUE);
	fprintf(log,"\n");
      }
      if (md->nU > 1) {
	fprintf(log,"%15s   %12s   %12s   %12s\n",
		"Group","Ux","Uy","Uz");
	for(i=0; (i<md->nU); i++) {
	  ni=atoms->grps[egcACC].nm_ind[i];
	  fprintf(log,"%15s",*atoms->grpname[ni]);
	  pr_ebin(log,md->ebin,md->iu+3*i,3,3,mode,steps,FALSE);
	}
	fprintf(log,"\n");
      }
    }
  }
}

