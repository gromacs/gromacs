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
 * Great Red Owns Many ACres of Sand 
 */
static char *SRCID_mdebin_c = "$Id$";

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
#include "names.h"

static bool bEInd[egNR] = { TRUE, TRUE, FALSE, FALSE, FALSE, FALSE };

t_mdebin *init_mdebin(int fp_ene,t_groups *grps,t_atoms *atoms,
		      bool bLR,bool bBHAM,bool b14)
{
  char *ener_nm[F_NRE];
  static char *boxs_nm[] = {
    "Box-X (nm)", "Box-Y (nm)", "Box-Z (nm)","Volume (nm^3)","Density (kg/l)"
  };
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
    "Pres-XX","Pres-XY","Pres-XZ",
    "Pres-YX","Pres-YY","Pres-YZ",
    "Pres-ZX","Pres-ZY","Pres-ZZ"
  };
  static char *mu_nm[] = {
    "Mu-X", "Mu-Y", "Mu-Z"
  };
  static   char   **grpnms;
  char     **gnm;
  char     buf[256];
  t_mdebin *md;
  int      i,j,ni,nj,n,k,kk;

  for(i=0; (i<F_NRE); i++)
    ener_nm[i]=interaction_function[i].longname;
    
  /* Energy monitoring */
  snew(md,1);
  md->ebin  = mk_ebin();
  md->ie    = get_ebin_space(md->ebin,F_NRE,ener_nm);
  md->ib    = get_ebin_space(md->ebin,asize(boxs_nm),boxs_nm);
  md->isvir = get_ebin_space(md->ebin,asize(sv_nm),sv_nm);
  md->ifvir = get_ebin_space(md->ebin,asize(fv_nm),fv_nm);
  md->ivir  = get_ebin_space(md->ebin,asize(vir_nm),vir_nm);
  md->ipres = get_ebin_space(md->ebin,asize(pres_nm),pres_nm);
  md->imu   = get_ebin_space(md->ebin,asize(mu_nm),mu_nm);
  if (bLR) 
    bEInd[egLR]   = TRUE;
  if (bBHAM) {
    bEInd[egLJ]   = FALSE;
    bEInd[egBHAM] = TRUE;
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
  if (md->nTC > 1) {
    snew(grpnms,2*md->nTC);
    for(i=0; (i<md->nTC); i++) {
      ni=atoms->grps[egcTC].nm_ind[i];
      sprintf(buf,"T-%s",*(atoms->grpname[ni]));
      grpnms[2*i]=strdup(buf);
      sprintf(buf,"Lamb-%s",*(atoms->grpname[ni]));
      grpnms[2*i+1]=strdup(buf);
    }
    md->itc=get_ebin_space(md->ebin,2*md->nTC,grpnms);
    sfree(grpnms);
  }
  else
    md->itc=0;
  
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

void upd_mdebin(t_mdebin *md,real tmass,int step,
		real ener[],
		matrix box,
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
  real   bs[5];
  real   eee[egNR];
  
  add_ebin(md->ebin,md->ie,F_NRE,ener,step);
  for(m=0; (m<DIM); m++) 
    bs[m]=box[m][m];
  bs[3] = bs[XX]*bs[YY]*bs[ZZ];
  bs[4] = (tmass*AMU)/(bs[3]*NANO*NANO*NANO*KILO);
  add_ebin(md->ebin,md->ib,5,bs,step);
  add_ebin(md->ebin,md->isvir,9,svir[0],step);
  add_ebin(md->ebin,md->ifvir,9,fvir[0],step);
  add_ebin(md->ebin,md->ivir,9,vir[0],step);
  add_ebin(md->ebin,md->ipres,9,pres[0],step);
  add_ebin(md->ebin,md->imu,3,mu_tot,step);
  
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
  if (md->nTC > 1) {
    if (ttt == NULL)
      snew(ttt,2*md->nTC);
    for(i=0; (i<md->nTC); i++) {
      ttt[2*i]   = grps->tcstat[i].T;
      ttt[2*i+1] = grps->tcstat[i].lambda;
    }
    add_ebin(md->ebin,md->itc,2*md->nTC,ttt,step);
  }
  
  if (md->nU > 1) {
    if (uuu == NULL)
      snew(uuu,md->nU);
    for(i=0; (i<md->nU); i++)
      copy_rvec(grps->grpstat[i].u,uuu[i]);
    add_ebin(md->ebin,md->iu,3*md->nU,uuu[0],step);
  }
}

static void npr(FILE *log,int n,char c)
{
  for(; (n>0); n--) fprintf(log,"%c",c);
}

static void pprint(FILE *log,char *s)
{
  static char *l1="==>";
  static char *l2="====>";
  static char *l3="======>";
  static char *r1="<======";
  static char *r2="<====";
  static char *r3="<==";
  char   CHAR='#';
  int slen;

  slen=strlen(s);
  fprintf(log,"\t%s  ",r1); npr(log,slen,CHAR); fprintf(log,"  %s\n",l1);
  fprintf(log,"\t%s  %s  %s\n",r2,s,l2); 
  fprintf(log,"\t%s  ",r3); npr(log,slen,CHAR); fprintf(log,"  %s\n\n",l3);
}

void print_ebin(int fp_ene,FILE *log,int steps,real time,real lamb,
		real SAfactor,int mode,bool bCompact,
		t_mdebin *md,t_groups *grps,t_atoms *atoms)
{
  static char **grpnms=NULL;
  static char *kjm="(kJ/mole)";
  t_drblock *drblock;
  char buf[246];
  int i,j,n,ni,nj;

  drblock=get_drblock();
  if (drblock->ndr == 0)
    drblock=NULL;
  switch (mode) {
  case eprNORMAL:
    if (fp_ene != -1)
      do_enx(fp_ene,&time,&steps,&md->ebin->nener,md->ebin->e,drblock);
    fprintf(log,"   %12s   %12s   %12s   %12s\n"
	    "   %12d   %12.5f   %12.5f   %12.5f\n\n",
	    "Step","Time","Lambda","Annealing",steps,time,lamb,SAfactor);
    break;
  case eprAVER:
    pprint(log,"A V E R A G E S");
    break;
  case eprRMS:
    pprint(log,"R M S - F L U C T U A T I O N S");
    break;
  default:
    fatal_error(0,"Invalid print mode (%d)",mode);
  }
#define newline() fprintf(log,"\n")

  fprintf(log,"   Energies %s\n",kjm);
  pr_ebin(log,md->ebin,md->ie,F_NRE,5,mode,steps,TRUE);  newline();
  
  if (bCompact)
    return;
    
  pr_ebin(log,md->ebin,md->ib,5,5,mode,steps,TRUE);      newline();
  fprintf(log,"   Shake Virial %s\n",kjm);
  pr_ebin(log,md->ebin,md->isvir,9,3,mode,steps,FALSE);  newline();
  fprintf(log,"   Force Virial %s\n",kjm);
  pr_ebin(log,md->ebin,md->ifvir,9,3,mode,steps,FALSE);  newline();
  fprintf(log,"   Total Virial %s\n",kjm);
  pr_ebin(log,md->ebin,md->ivir,9,3,mode,steps,FALSE);   newline();
  fprintf(log,"   Pressure (Bar)\n");
  pr_ebin(log,md->ebin,md->ipres,9,3,mode,steps,FALSE);  newline();
  fprintf(log,"   Total Dipole (Debye)\n");
  pr_ebin(log,md->ebin,md->imu,3,3,mode,steps,FALSE);    newline();
  
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
    newline();
  }
  if (md->nTC > 1) {
    pr_ebin(log,md->ebin,md->itc,2*md->nTC,4,mode,steps,TRUE);
    newline();
  }
  if (md->nU > 1) {
    fprintf(log,"%15s   %12s   %12s   %12s\n",
	    "Group","Ux","Uy","Uz");
    for(i=0; (i<md->nU); i++) {
      ni=atoms->grps[egcACC].nm_ind[i];
      fprintf(log,"%15s",*atoms->grpname[ni]);
      pr_ebin(log,md->ebin,md->iu+3*i,3,3,mode,steps,FALSE);
    }
    newline();
  }
#undef newline
}

