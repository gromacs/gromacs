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
static char *SRCID_g_chi_c = "$Id$";

#include <stdio.h>
#include <math.h>
#include "confio.h"
#include "pdbio.h"
#include "copyrite.h"
#include "fatal.h"
#include "futil.h"
#include "assert.h"
#include "gstat.h"
#include "macros.h"
#include "maths.h"
#include "physics.h"
#include "rdgroup.h"
#include "smalloc.h"
#include "statutil.h"
#include "tpxio.h"
#include "string.h"
#include "sysstuff.h"
#include "txtdump.h"
#include "typedefs.h"
#include "vec.h"
#include "strdb.h"
#include "xvgr.h"
#include "pp2shift.h"
#include "matio.h"

static bool bAllowed(real phi,real psi)
{
  static char *map[] = {
    "1100000000000000001111111000000000001111111111111111111111111",
    "1100000000000000001111110000000000011111111111111111111111111",
    "1100000000000000001111110000000000011111111111111111111111111",
    "1100000000000000001111100000000000111111111111111111111111111",
    "1100000000000000001111100000000000111111111111111111111111111",
    "1100000000000000001111100000000001111111111111111111111111111",
    "1100000000000000001111100000000001111111111111111111111111111",
    "1100000000000000001111100000000011111111111111111111111111111",
    "1110000000000000001111110000000111111111111111111111111111111",
    "1110000000000000001111110000001111111111111111111111111111111",
    "1110000000000000001111111000011111111111111111111111111111111",
    "1110000000000000001111111100111111111111111111111111111111111",
    "1110000000000000001111111111111111111111111111111111111111111",
    "1110000000000000001111111111111111111111111111111111111111111",
    "1110000000000000001111111111111111111111111111111111111111111",
    "1110000000000000001111111111111111111111111111111111111111111",
    "1110000000000000001111111111111110011111111111111111111111111",
    "1110000000000000001111111111111100000111111111111111111111111",
    "1110000000000000001111111111111000000000001111111111111111111",
    "1100000000000000001111111111110000000000000011111111111111111",
    "1100000000000000001111111111100000000000000011111111111111111",
    "1000000000000000001111111111000000000000000001111111111111110",
    "0000000000000000001111111110000000000000000000111111111111100",
    "0000000000000000000000000000000000000000000000000000000000000",
    "0000000000000000000000000000000000000000000000000000000000000",
    "0000000000000000000000000000000000000000000000000000000000000",
    "0000000000000000000000000000000000000000000000000000000000000",
    "0000000000000000000000000000000000000000000000000000000000000",
    "0000000000000000000000000000000000000000000000000000000000000",
    "0000000000000000000000000000000000000000000000000000000000000",
    "0000000000000000000000000000000000000000000000000000000000000",
    "0000000000000000000000000000000000000000000000000000000000000",
    "0000000000000000000000000000000000000000000000000000000000000",
    "0000000000000000000000000000000000000000000000000000000000000",
    "0000000000000000000000000000000000000000000000000000000000000",
    "0000000000000000000000000000000000000000000000000000000000000",
    "0000000000000000000000000000000000000000000000000000000000000",
    "0000000000000000000000000000000000000000000000000000000000000",
    "0000000000000000000000000000000000111111111111000000000000000",
    "1100000000000000000000000000000001111111111111100000000000111",
    "1100000000000000000000000000000001111111111111110000000000111",
    "0000000000000000000000000000000000000000000000000000000000000",
    "0000000000000000000000000000000000000000000000000000000000000",
    "0000000000000000000000000000000000000000000000000000000000000",
    "0000000000000000000000000000000000000000000000000000000000000",
    "0000000000000000000000000000000000000000000000000000000000000",
    "0000000000000000000000000000000000000000000000000000000000000",
    "0000000000000000000000000000000000000000000000000000000000000",
    "0000000000000000000000000000000000000000000000000000000000000",
    "0000000000000000000000000000000000000000000000000000000000000",
    "0000000000000000000000000000000000000000000000000000000000000",
    "0000000000000000000000000000000000000000000000000000000000000",
    "0000000000000000000000000000000000000000000000000000000000000",
    "0000000000000000000000000000000000000000000000000000000000000",
    "0000000000000000000000000000000000000000000000000000000000000",
    "0000000000000000000000000000000000000000000000000000000000000",
    "0000000000000000000000000000000000000000000000000000000000000",
    "0000000000000000000000000000000000000000000000000000000000000",
    "0000000000000000000000000000000000000000000000000000000000000",
    "0000000000000000000000000000000000000000000000000000000000000",
    "0000000000000000000000000000000000000000000000000000000000000"
  };
#define NPP asize(map)
  int x,y;

#define INDEX(ppp) ((((int) (360+ppp*RAD2DEG)) % 360)/6)
  x = INDEX(phi);
  y = INDEX(psi);
#undef INDEX
  return (bool) map[x][y];
}

atom_id *make_chi_ind(int nl,t_dlist dl[],int *ndih)
{
  atom_id *id;
  int     i,Xi,n;
  
  /* There are nl sidechains with max edMax dihedrals with 4 atoms each */
  snew(id,nl*edMax*4); 
  
  n=0;
  for(i=0; (i<nl); i++) { /* Phi */
    dl[i].j0[edPhi] = n/4;
    if (dl[i].atm.H == -1)
      id[n++]=dl[i].atm.minC;
    else
      id[n++]=dl[i].atm.H;
    id[n++]=dl[i].atm.N;
    id[n++]=dl[i].atm.Cn[1];
    id[n++]=dl[i].atm.C;
  }
  for(i=0; (i<nl); i++) { /* Psi */
    dl[i].j0[edPsi] = n/4;
    id[n++]=dl[i].atm.N;
    id[n++]=dl[i].atm.Cn[1];
    id[n++]=dl[i].atm.C;
    id[n++]=dl[i].atm.O;
  }
  for(i=0; (i<nl); i++) { /* Omega */
    if (has_dihedral(edOmega,&(dl[i]))) {
      dl[i].j0[edOmega] = n/4;
      id[n++]=dl[i].atm.minO;
      id[n++]=dl[i].atm.minC;
      id[n++]=dl[i].atm.N;
      id[n++]=dl[i].atm.Cn[1];
    }
  }
  for(Xi=0; (Xi<MAXCHI); Xi++) { /* Chi# */
    for(i=0; (i<nl); i++) {
      if (dl[i].atm.Cn[Xi+3] != -1) {
	dl[i].j0[edChi1+Xi] = n/4;
	id[n++]=dl[i].atm.Cn[Xi];
	id[n++]=dl[i].atm.Cn[Xi+1];
	id[n++]=dl[i].atm.Cn[Xi+2];
	id[n++]=dl[i].atm.Cn[Xi+3];
      }
    }
  }
  *ndih=n/4;
  
  return id;
}

int bin(real chi,int mult)
{
  mult=3;

  return (int) (chi*mult/360.0);
}

static void print_one(char *base,char *name,char *title,
		      int nf,real time[],real data[])
{
  FILE *fp;
  char buf[256],t2[256];
  int  k;
  
  sprintf(buf,"%s%s.xvg",base,name);
  fprintf(stderr,"\rPrinting %s  ",buf);
  sprintf(t2,"%s %s",title,name);
  fp=xvgropen(buf,t2,"Time (ps)","C(t)");
  for(k=0; (k<nf); k++)
    fprintf(fp,"%10g  %10g\n",time[k],data[k]);
  ffclose(fp);
}

static void do_dihcorr(char *fn,int nf,int ndih,real **dih,real dt,
		       int nlist,t_dlist dlist[],real time[],int maxchi,
		       bool bPhi,bool bPsi,bool bChi,bool bOmega)
{
  char name1[256],name2[256];
  int  i,j,Xi;
  
  do_autocorr(fn,"Dihedral Autocorrelation Function",
	      nf,ndih,dih,dt,eacCos,FALSE);
  /* Dump em all */
  j=0;
  for(i=0; (i<nlist); i++) {
    if (bPhi)
      print_one("corrphi",dlist[i].name,"Phi ACF for",nf/2,time,dih[j]);
    j++;
  }
  for(i=0; (i<nlist); i++) {
    if (bPsi)
      print_one("corrpsi",dlist[i].name,"Psi ACF for",nf/2,time,dih[j]);
    j++;
  }
  for(i=0; (i<nlist); i++) {
    if (has_dihedral(edOmega,&dlist[i])) {
      if (bOmega)
	print_one("corromega",dlist[i].name,"Omega ACF for",nf/2,time,dih[j]);
      j++;
    }
  }
  for(Xi=0; (Xi<maxchi); Xi++) {
    sprintf(name1, "corrchi%d", Xi+1);
    sprintf(name2, "Chi%d ACF for", Xi+1);
    for(i=0; (i<nlist); i++) {
      if (dlist[i].atm.Cn[Xi+3] != -1) {
	if (bChi)
	  print_one(name1,dlist[i].name,name2,nf/2,time,dih[j]);
	j++;
      }
    }
  }
  fprintf(stderr,"\n");
}

static void dump_em_all(int nlist,t_dlist dlist[],int nf,real time[],
			real **dih,int maxchi,
			bool bPhi,bool bPsi,bool bChi,bool bOmega)
{
  char name[256];
  int  i,j,Xi;
  
  /* Dump em all */
  j = 0;
  for(i=0; (i<nlist); i++) {
    if (bPhi)
      print_one("phi",dlist[i].name,name,nf,time,dih[j]);
      j++;
  }
  for(i=0; (i<nlist); i++) {
    if (bPsi)
      print_one("psi",dlist[i].name,name,nf,time,dih[j]);
    j++;
  }  
  for(i=0; (i<nlist); i++)
    if (has_dihedral(edOmega,&(dlist[i]))) {
      if (bOmega)
	print_one("omega",dlist[i].name,name,nf,time,dih[j]);
      j++;
    }
  
  for(Xi=0; (Xi<maxchi); Xi++)
    for(i=0; (i<nlist); i++)
      if (dlist[i].atm.Cn[Xi+3] != -1) {
	if (bChi) {
	  sprintf(name,"chi%d",Xi+1);
	  print_one(name,dlist[i].name,name,nf,time,dih[j]);
	}
	j++;
      }
  fprintf(stderr,"\n");
}

static void reset_one(real dih[],int nf,real phase)
{
  int j;
  
  for(j=0; (j<nf); j++) {
    dih[j] += phase;
    while(dih[j] < -M_PI)
      dih[j] += 2*M_PI;
    while(dih[j] >= M_PI)
      dih[j] -= 2*M_PI;
  }
}

static void reset_em_all(int nlist,t_dlist dlist[],int nf,
			 real **dih,int maxchi,bool bPhi,bool bPsi,bool bChi)
{
  int  i,j,Xi;
  
  /* Reset em all */
  j=0;
  /* Phi */
  for(i=0; (i<nlist); i++)
    if (dlist[i].atm.H == -1)
      reset_one(dih[j++],nf,0);
    else
      reset_one(dih[j++],nf,M_PI);
  /* Psi */
  for(i=0; (i<nlist); i++)
    reset_one(dih[j++],nf,M_PI);
  /* Omega */
  for(i=0; (i<nlist); i++)
    if (has_dihedral(edOmega,&dlist[i]))
      reset_one(dih[j++],nf,0);
  /* Chi 1 thru maxchi */
  for(Xi=0; (Xi<maxchi); Xi++)
    for(i=0; (i<nlist); i++)
      if (dlist[i].atm.Cn[Xi+3] != -1) {
	reset_one(dih[j],nf,0);
	j++;
      }
  fprintf(stderr,"j after resetting = %d\n",j);
}

static void histogramming(FILE *log,int naa,char **aa,
			  int nf,int maxchi,real **dih,
			  int nlist,t_dlist dlist[],
			  atom_id index[],
			  bool bPhi,bool bPsi,bool bOmega,bool bChi,
			  bool bNormalize,bool bSSHisto,char *ssdump,
			  real bfac_max,t_atoms *atoms)
{
  t_karplus kkkphi[] = {
    { "J_NHa",     6.51, -1.76,  1.6, -M_PI/3,   0.0 },
    { "J_HaC'",    4.0,   1.1,   0.1,  0.0,      0.0 },
    { "J_NHCb",    4.7,  -1.5,  -0.2,  M_PI/3,   0.0 },
    { "J_Ci-1Hai", 4.5,  -1.3,  -1.2,  2*M_PI/3, 0.0 }
  };
  t_karplus kkkpsi[] = {
    { "J_HaN",   -0.88, -0.61,-0.27,M_PI/3,  0.0 }
  };
  t_karplus kkkchi1[] = {
    { "JHaHb2",       9.5, -1.6, 1.8, -M_PI/3, 0 },
    { "JHaHb3",       9.5, -1.6, 1.8, 0, 0 }
  };
#define NKKKPHI asize(kkkphi)
#define NKKKPSI asize(kkkpsi)
#define NKKKCHI asize(kkkchi1)
#define NJC (NKKKPHI+NKKKPSI+NKKKCHI)
  
  FILE    *fp,*ssfp[3];
  char    *sss[3] = { "sheet", "helix", "coil" };
  real    S2;
  real    *normhisto;
  real    **Jc;
  int     ****his_aa_ss=NULL;
  int     ***his_aa,**his_aa1,*histmp;
  int     i,j,k,m,n,nn,Dih,nres,hindex;
  bool    bBfac,bOccup;
  char    hisfile[256],hhisfile[256],sshisfile[256],title[256],*ss_str=NULL;
  
  if (bSSHisto) {
    fp = ffopen(ssdump,"r");
    fscanf(fp,"%d",&nres);
    snew(ss_str,nres+1);
    fscanf(fp,"%s",ss_str);
    ffclose(fp);
    /* Four dimensional array... Very cool */
    snew(his_aa_ss,3);
    for(i=0; (i<3); i++) {
      snew(his_aa_ss[i],naa+1);
      for(j=0; (j<=naa); j++) {
	snew(his_aa_ss[i][j],edMax);
	for(Dih=0; (Dih<edMax); Dih++)
	  snew(his_aa_ss[i][j][Dih],NHISTO+1);
      }
    }
  }
  snew(his_aa,edMax);
  for(Dih=0; (Dih<edMax); Dih++) {
    snew(his_aa[Dih],naa+1);
    for(i=0; (i<=naa); i++) {
      snew(his_aa[Dih][i],NHISTO+1);
    }
  }
  snew(histmp,NHISTO);
  
  snew(Jc,nlist);
  for(i=0; (i<nlist); i++)
    snew(Jc[i],NJC);
  
  j=0;
  n=0;
  for (Dih=0; (Dih<NONCHI+maxchi); Dih++) {    
    for(i=0; (i<nlist); i++) {
      if (((Dih  < edOmega) ) ||
	  ((Dih == edOmega) && (has_dihedral(edOmega,&(dlist[i])))) ||
	  ((Dih  > edOmega) && (dlist[i].atm.Cn[Dih-NONCHI+3] != -1))) {
      	make_histo(log,nf,dih[j],NHISTO,histmp,-M_PI,M_PI);
	
	if (bSSHisto) {
	  /* Assume there is only one structure, the first. 
	   * Compute index in histogram.
	   */
	  /* Check the atoms to see whether their B-factors are low enough 
	   * Check atoms to see their occupancy is 1.
	   */
	  bBfac = bOccup = TRUE;
	  for(nn=0; (nn<4); nn++,n++) {
	    bBfac  = bBfac  && (atoms->pdbinfo[index[n]].bfac <= bfac_max);
	    bOccup = bOccup && (atoms->pdbinfo[index[n]].occup == 1);
	  }
	  if (bOccup && ((bfac_max <= 0) || ((bfac_max > 0) && bBfac))) {
	    hindex = ((dih[j][0]+M_PI)*NHISTO)/(2*M_PI);
	    assert(hindex >= 0);
	    assert(hindex < NHISTO);
	    /* Assign dihedral to either of the structure determined 
	     * histograms
	     */
	    switch(ss_str[dlist[i].resnr]) {
	    case 'E':
	      his_aa_ss[0][dlist[i].index][Dih][hindex]++;
	      break;
	    case 'H':
	      his_aa_ss[1][dlist[i].index][Dih][hindex]++;
	      break;
	    default:
	      his_aa_ss[2][dlist[i].index][Dih][hindex]++;
	      break;
	    }
	  }
	  else if (debug) 
	    fprintf(debug,"Res. %d has imcomplete occupancy or bfacs > %g\n",
		    dlist[i].resnr,bfac_max);
	}
	else
	  n += 4;
	  
	switch (Dih) {
	case edPhi:
	  calc_distribution_props(NHISTO,histmp,-M_PI,NKKKPHI,kkkphi,&S2);
	  
	  for(m=0; (m<NKKKPHI); m++) 
	    Jc[i][m] = kkkphi[m].Jc;
	  break;
	case edPsi:
	  calc_distribution_props(NHISTO,histmp,-M_PI,NKKKPSI,kkkpsi,&S2);
	  
	  for(m=0; (m<NKKKPSI); m++)
	    Jc[i][NKKKPHI+m] = kkkpsi[m].Jc;
	  break;
	case edOmega:
	  calc_distribution_props(NHISTO,histmp,-M_PI,0,NULL,&S2);
	  break;
	case edChi1:
	  calc_distribution_props(NHISTO,histmp,-M_PI,NKKKCHI,kkkchi1,&S2);
	  for(m=0; (m<NKKKCHI); m++)
	    Jc[i][NKKKPHI+NKKKPSI+m] = kkkchi1[m].Jc;
	  break;
	}
	dlist[i].S2[Dih]        = S2;
	
	/* Sum distribution per amino acid type as well */
	for(k=0; (k<NHISTO); k++) {
	  his_aa[Dih][dlist[i].index][k] += histmp[k];
	  histmp[k] = 0;
	}
	j++;
      }
    }
  }
  sfree(histmp);
  
  /* Print out Jcouplings */
  fprintf(log,"\n *** J-Couplings from simulation ***\n\n");
  fprintf(log,"Residue   ");
  for(i=0; (i<NKKKPHI); i++)
    fprintf(log,"%10s",kkkphi[i].name);
  for(i=0; (i<NKKKPSI); i++)
    fprintf(log,"%10s",kkkpsi[i].name);
  for(i=0; (i<NKKKCHI); i++)
    fprintf(log,"%10s",kkkchi1[i].name);
  fprintf(log,"\n");
  for(i=0; (i<NJC+1); i++)
    fprintf(log,"----------");
  fprintf(log,"\n");
  for(i=0; (i<nlist); i++) {
    fprintf(log,"%-10s",dlist[i].name);
    for(j=0; (j<NJC); j++)
      fprintf(log,"  %8.3f",Jc[i][j]);
    fprintf(log,"\n");
  }
  fprintf(log,"\n");
    
  snew(normhisto,NHISTO);
  for(i=0; (i<naa); i++) {
    for(Dih=0; (Dih<edMax); Dih++){
      /* First check whether something is in there */
      for(j=0; (j<NHISTO); j++)
	if (his_aa[Dih][i][j] != 0)
	  break;
      if ((j < NHISTO) &&
	  ((bPhi && (Dih==edPhi)) ||
	   (bPsi && (Dih==edPsi)) ||
	   (bOmega &&(Dih==edOmega)) ||
	   (bChi && (Dih>=edChi1)))) {
	if (bNormalize)
	  normalize_histo(NHISTO,his_aa[Dih][i],(360.0/NHISTO),normhisto);
	
	switch (Dih) {
	case edPhi:
	  sprintf(hisfile,"histo-phi%s",aa[i]);
	  sprintf(title,"\\8f\\4 Distribution for %s",aa[i]);
	  break;
	case edPsi:
	  sprintf(hisfile,"histo-psi%s",aa[i]);
	  sprintf(title,"\\8y\\4 Distribution for %s",aa[i]);
	  break;
	case edOmega:
	  sprintf(hisfile,"histo-omega%s",aa[i]);
	  sprintf(title,"\\8w\\4 Distribution for %s",aa[i]);
	  break;
	default:
	  sprintf(hisfile,"histo-chi%d%s",Dih-NONCHI+1,aa[i]);
	  sprintf(title,"\\8c\\4\\s%d\\N Distribution for %s",
		  Dih-NONCHI+1,aa[i]);
	}
	strcpy(hhisfile,hisfile);
	strcat(hhisfile,".xvg");
	fp=xvgropen(hhisfile,title,"Degrees","");
	fprintf(fp,"@ with g0\n");
	xvgr_world(fp,-180,0,180,0.1);
	fprintf(fp,"@ xaxis tick on\n");
	fprintf(fp,"@ xaxis tick major 90\n");
	fprintf(fp,"@ xaxis tick minor 30\n");
	fprintf(fp,"@ xaxis ticklabel prec 0\n");
	fprintf(fp,"@ yaxis tick off\n");
	fprintf(fp,"@ yaxis ticklabel off\n");
	fprintf(fp,"@ type xy\n");
	if (bSSHisto) {
	  for(k=0; (k<3); k++) {
	    sprintf(sshisfile,"%s-%s.xvg",hisfile,sss[k]);
	    ssfp[k] = ffopen(sshisfile,"w");
	  }
	}
	for(j=0; (j<NHISTO); j++) {
	  if (bNormalize)
	    fprintf(fp,"%5d  %10g\n",j-180,normhisto[j]);
	  else
	    fprintf(fp,"%5d  %10d\n",j-180,his_aa[Dih][i][j]);
	  if (bSSHisto)
	    for(k=0; (k<3); k++) 
	      fprintf(ssfp[k],"%5d  %10d\n",j-180,
		      his_aa_ss[k][i][Dih][j]);
	}
	fprintf(fp,"&\n");
	ffclose(fp);
	if (bSSHisto) {
	  for(k=0; (k<3); k++) {
	    fprintf(ssfp[k],"&\n");
	    ffclose(ssfp[k]);
	  }
	}
      }
    }
  }
  sfree(normhisto);
  
  if (bSSHisto) {
    /* Four dimensional array... Very cool */
    for(i=0; (i<3); i++) {
      for(j=0; (j<=naa); j++) {
	for(Dih=0; (Dih<edMax); Dih++)
	  sfree(his_aa_ss[i][j][Dih]);
	sfree(his_aa_ss[i][j]);
      }
      sfree(his_aa_ss[i]);
    }
    sfree(his_aa_ss);
    sfree(ss_str);
  }
}

static FILE *rama_file(char *fn,char *title,char *xaxis,char *yaxis)
{
  FILE *fp;

  fp = xvgropen(fn,title,xaxis,yaxis);  
  fprintf(fp,"@ with g0\n");
  xvgr_world(fp,-180,-180,180,180);
  fprintf(fp,"@ xaxis tick on\n");
  fprintf(fp,"@ xaxis tick major 90\n");
  fprintf(fp,"@ xaxis tick minor 30\n");
  fprintf(fp,"@ xaxis ticklabel prec 0\n");
  fprintf(fp,"@ yaxis tick on\n");
  fprintf(fp,"@ yaxis tick major 90\n");
  fprintf(fp,"@ yaxis tick minor 30\n");
  fprintf(fp,"@ yaxis ticklabel prec 0\n");
  fprintf(fp,"@    s0 type xy\n");
  fprintf(fp,"@    s0 symbol 2\n");
  fprintf(fp,"@    s0 symbol size 0.410000\n");
  fprintf(fp,"@    s0 symbol fill 1\n");
  fprintf(fp,"@    s0 symbol color 1\n");
  fprintf(fp,"@    s0 symbol linewidth 1\n");
  fprintf(fp,"@    s0 symbol linestyle 1\n");
  fprintf(fp,"@    s0 symbol center false\n");
  fprintf(fp,"@    s0 symbol char 0\n");
  fprintf(fp,"@    s0 skip 0\n");
  fprintf(fp,"@    s0 linestyle 0\n");
  fprintf(fp,"@    s0 linewidth 1\n");
  fprintf(fp,"@ type xy\n");
  
  return fp;
}

static void do_rama(int nf,int nlist,t_dlist dlist[],real **dih,
		    bool bViol,bool bRamOmega)
{
  FILE *fp,*gp=NULL;
  bool bOm;
  char fn[256];
  int  i,j,k,Xi1,Xi2,Phi,Psi,Om=0,nlevels;
#define NMAT 120
  real **mat=NULL,phi,psi,omega,axis[NMAT],lo,hi;
  t_rgb rlo = { 1.0, 0.0, 0.0 };
  t_rgb rmid= { 1.0, 1.0, 1.0 };
  t_rgb rhi = { 0.0, 0.0, 1.0 };
  
  for(i=0; (i<nlist); i++) {
    if ((has_dihedral(edPhi,&(dlist[i]))) &&
	(has_dihedral(edPsi,&(dlist[i])))) {
      sprintf(fn,"ramaPhiPsi%s.xvg",dlist[i].name);
      fp = rama_file(fn,"Ramachandran Plot",
		     "\\8f\\4 (deg)","\\8y\\4 (deg)");
      bOm = bRamOmega && has_dihedral(edOmega,&(dlist[i]));
      if (bOm) {
	Om = dlist[i].j0[edOmega];
	snew(mat,NMAT);
	for(j=0; (j<NMAT); j++) {
	  snew(mat[j],NMAT);
	  axis[j] = -180+(360*j)/NMAT;
	}
      }
      if (bViol) {
	sprintf(fn,"violPhiPsi%s.xvg",dlist[i].name);
	gp = ffopen(fn,"w");
      }
      Phi = dlist[i].j0[edPhi];   
      Psi = dlist[i].j0[edPsi];
      for(j=0; (j<nf); j++) {
	phi = RAD2DEG*dih[Phi][j];
	psi = RAD2DEG*dih[Psi][j];
	fprintf(fp,"%10g  %10g\n",phi,psi);
	if (bViol)
	  fprintf(gp,"%d\n",!bAllowed(dih[Phi][j],RAD2DEG*dih[Psi][j]));
	if (bOm) {
	  omega = RAD2DEG*dih[Om][j];
	  mat[(int)((phi*NMAT)/360)+NMAT/2][(int)((psi*NMAT)/360)+NMAT/2] 
	    += omega;
	}
      }
      if (bViol)
	fclose(gp);
      fclose(fp);
      if (bOm) {
	sprintf(fn,"ramomega%s.xpm",dlist[i].name);
	fp = ffopen(fn,"w");
	lo = hi = 0;
	for(j=0; (j<NMAT); j++)
	  for(k=0; (k<NMAT); k++) {
	    mat[j][k] /= nf;
	    lo=min(mat[j][k],lo);
	    hi=max(mat[j][k],hi);
	  }
	/* Symmetrise */
	if (fabs(lo) > fabs(hi)) 
	  hi = -lo;
	else
	  lo = -hi;
	/* Add 180 */
	for(j=0; (j<NMAT); j++)
	  for(k=0; (k<NMAT); k++)
	    mat[j][k] += 180;
	lo += 180;
	hi += 180;
	nlevels = 20;
	write_xpm3(fp,"Omega/Ramachandran Plot","Deg","Phi","Psi",
		   NMAT,NMAT,axis,axis,mat,lo,180.0,hi,rlo,rmid,rhi,&nlevels);
	fclose(fp);
	for(j=0; (j<NMAT); j++)
	  sfree(mat[j]);
	sfree(mat);
      }
    }
    if ((has_dihedral(edChi1,&(dlist[i]))) &&
	(has_dihedral(edChi2,&(dlist[i])))) {
      sprintf(fn,"ramaX1X2%s.xvg",dlist[i].name);
      fp = rama_file(fn,"\\8c\\4\\s1\\N-\\8c\\4\\s2\\N Ramachandran Plot",
		     "\\8c\\4\\s1\\N (deg)","\\8c\\4\\s2\\N (deg)");
      Xi1 = dlist[i].j0[edChi1];
      Xi2 = dlist[i].j0[edChi2];
      for(j=0; (j<nf); j++)
	fprintf(fp,"%10g  %10g\n",RAD2DEG*dih[Xi1][j],RAD2DEG*dih[Xi2][j]);
      fclose(fp);
    }
    else 
      fprintf(stderr,"No chi1 & chi2 angle for %s\n",dlist[i].name);
  }
}

static void order_params(FILE *log,
			 char *fn,int maxchi,int nlist,t_dlist dlist[],
			 char *pdbfn,real bfac_init,
			 t_atoms *atoms,rvec x[],matrix box,
			 bool bPhi,bool bPsi,bool bChi)
{
  FILE *fp;
  int  nh[edMax];
  int  i,Dih,Xi;
  real S2Max, S2Min;
  
  /* Print order parameters */
  fp=xvgropen(fn,"Dihedral Order Parameters","Residue","S2");
  
  for (Dih=0; (Dih<edMax); Dih++)
    nh[Dih]=0;
  
  fprintf(fp,"%5s  ","#Res.");
  fprintf(fp,"%10s  %10s  ","S2Min","S2Max");
  fprintf(fp,"%10s  %10s  ","Phi","Psi");
  for (Xi=1; (Xi<=maxchi); Xi++)
    fprintf(fp,"%8s%2d  ","Xi",Xi);
  fprintf(fp,"%12s\n","Res. Name");
  
  for(i=0; (i<nlist); i++) {
    S2Max=-10;
    S2Min=10;
    for (Dih=0; (Dih<NONCHI+maxchi); Dih++) {
      if (dlist[i].S2[Dih]!=0) {
	if (dlist[i].S2[Dih] > S2Max) 
	  S2Max=dlist[i].S2[Dih];
	if (dlist[i].S2[Dih] < S2Min) 
	  S2Min=dlist[i].S2[Dih];
      }
      if (dlist[i].S2[Dih] > 0.8)
	nh[Dih]++;
    }
    fprintf(fp,"%5d  ",dlist[i].resnr);
    fprintf(fp,"%10.3f  %10.3f  ",S2Min,S2Max);
    for (Dih=0; (Dih<NONCHI+maxchi); Dih++)
      fprintf(fp,"%10.3f  ",dlist[i].S2[Dih]);
    fprintf(fp,"%12s\n",dlist[i].name);
  }
  ffclose(fp);
  
  if (pdbfn) {
    int Xi;
    real x0,y0,z0;
    
    for(i=0; (i<atoms->nr); i++)
      atoms->pdbinfo[i].bfac=bfac_init;
    
    for(i=0; (i<nlist); i++) {
      atoms->pdbinfo[dlist[i].atm.N].bfac=-dlist[i].S2[0];/* Phi */
      atoms->pdbinfo[dlist[i].atm.H].bfac=-dlist[i].S2[0];/* Phi */
      atoms->pdbinfo[dlist[i].atm.C].bfac=-dlist[i].S2[1];/* Psi */
      atoms->pdbinfo[dlist[i].atm.O].bfac=-dlist[i].S2[1];/* Psi */
      for (Xi=0; (Xi<maxchi); Xi++) {           /* Chi's */
	if (dlist[i].atm.Cn[Xi+3]!=-1) {
	  atoms->pdbinfo[dlist[i].atm.Cn[Xi+1]].bfac=-dlist[i].S2[NONCHI+Xi];
	}
      }
    }
    
    fp=ffopen(pdbfn,"w");
    fprintf(fp,"REMARK generated by g_chi\n");
    fprintf(fp,"REMARK "
	    "B-factor field contains negative of dihedral order parameters\n");
    write_pdbfile(fp,NULL,atoms,x,box,0,0);
    x0=y0=z0=1000.0;
    for (i=0; (i<atoms->nr); i++) {
      x0 = min(x0, x[i][XX]);
      y0 = min(y0, x[i][YY]);
      z0 = min(z0, x[i][ZZ]);
    }
    x0*=10.0;/* nm -> angstrom */
    y0*=10.0;/* nm -> angstrom */
    z0*=10.0;/* nm -> angstrom */
    for (i=0; (i<10); i++)
      fprintf(fp,pdbformat,"ATOM  ", atoms->nr+1+i, "CA", "LEG",' ', nres+1,
	      x0, y0, z0+(1.2*i), 0.0, -0.1*i);
    ffclose(fp);
  }
  
  fprintf(log,"Dihedrals with S2 > 0.8\n");
  fprintf(log,"Dihedral: ");
  if (bPhi) 
    fprintf(log," Phi  ");
  if (bPsi) 
    fprintf(log," Psi  ");
  if (bChi)
    for(Xi=0; (Xi<maxchi); Xi++)
      fprintf(log,"Chi%d  ",Xi+1);
  fprintf(log,"\nNumber:   ");
  if (bPhi) 
    fprintf(log,"%4d  ",nh[0]);
  if (bPsi) 
    fprintf(log,"%4d  ",nh[1]);
  if (bChi)
    for(Xi=0; (Xi<maxchi); Xi++)
      fprintf(log,"%4d  ",nh[NONCHI+Xi]);
  fprintf(log,"\n");
}

int main(int argc,char *argv[])
{
  static char *desc[] = {
    "g_chi computes phi, psi, omega and chi dihedrals for all your ",
    "amino acid backbone and sidechains.",
    "It can compute dihedral angle as a function of time, and as",
    "histogram distributions.",
    "Output is in form of xvgr files, as well as a LaTeX table of the",
    "number of transitions per nanosecond.[PAR]",
    "Order parameters S2 for each of the dihedrals are calculated and",
    "output as xvgr file and optionally as a pdb file with the S2",
    "values as B-factor.[PAR]",
    "If option [TT]-c[tt] is given, the program will",
    "calculate dihedral autocorrelation functions. The function used",
    "is C(t) = < cos(chi(tau)) cos(chi(tau+t)) >. The use of cosines",
    "rather than angles themselves, resolves the problem of periodicity.",
    "(Van der Spoel & Berendsen (1997), [BB]Biophys. J. 72[bb], 2032-2041).[PAR]",
    "The option [TT]-r[tt] generates a contour plot of the average omega angle",
    "as a function of the phi and psi angles, that is, in a Ramachandran plot",
    "the average omega angle is plotted using color coding."
  };
  
  static char *bugs[] = {
    "Produces MANY output files (up to about 4 times the number of residues in the protein, twice that if autocorrelation functions are calculated). Typically several hundred files are output."
  };
  static int  r0=1,ndeg=1,maxchi=2;
  static bool bAll=FALSE;
  static bool bPhi=FALSE,bPsi=FALSE,bOmega=FALSE;
  static real bfac_init=-1.0,bfac_max=0;
  static char *maxchistr[] = { NULL, "0", "1", "2", "3",  "4", "5", "6", NULL };
  static bool bRama=FALSE,bShift=FALSE,bViol=FALSE,bRamOmega=FALSE;
  static bool bNormHisto=TRUE;
  t_pargs pa[] = {
    { "-r0",  FALSE, etINT, {&r0},
      "starting residue" },
    { "-phi",  FALSE, etBOOL, {&bPhi},
      "Output for Phi dihedral angles" },
    { "-psi",  FALSE, etBOOL, {&bPsi},
      "Output for Psi dihedral angles" },
    { "-omega",FALSE, etBOOL, {&bOmega},  
      "Output for Omega dihedrals (peptide bonds)" },
    { "-rama", FALSE, etBOOL, {&bRama},
      "Generate Phi/Psi and Chi1/Chi2 ramachandran plots" },
    { "-viol", FALSE, etBOOL, {&bViol},
      "Write a file that gives 0 or 1 for violated Ramachandran angles" },
    { "-all",  FALSE, etBOOL, {&bAll},
      "Output separate files for every dihedral." },
    { "-shift", FALSE, etBOOL, {&bShift},
	"Compute chemical shifts from Phi/Psi angles" },
    { "-run", FALSE, etINT, {&ndeg},
      "perform running average over ndeg degrees for histograms" },
    { "-maxchi", FALSE, etENUM, {maxchistr},
      "calculate first ndih Chi dihedrals" },
    { "-normhisto", FALSE, etBOOL, {&bNormHisto},
      "Normalize histograms" },
    { "-ramomega",FALSE,etBOOL, {&bRamOmega},
      "compute average omega as a function of phi/psi and plot it in an xpm plot" },
    { "-bfact", FALSE, etREAL, {&bfac_init},
      "B-factor value for pdb file for atoms with no calculated dihedral order parameter"},
    { "-bmax",  FALSE, etREAL, {&bfac_max},
      "Maximum B-factor on any of the atoms that make up a dihedral, for the dihedral angle to be considere in the statistics. Applies to database work where a number of X-Ray structures is analyzed. -bmax <= 0 means no limit." }
  };

  FILE       *log;
  int        natoms,nlist,naa,idum;
  t_atoms    atoms;
  rvec       *x;
  matrix     box;
  char       title[256];
  t_dlist    *dlist;
  char       **aa;
  bool       bChi,bCorr,bSSHisto;
  real       dt=0;

  atom_id    isize,*index;
  int        ndih,nf;
  real       **dih,*trans_frac,*aver_angle,*time;
  
  t_filenm  fnm[] = {
    { efSTX, NULL,  NULL,     ffREAD  },
    { efTRX, "-f",  NULL,     ffREAD  },
    { efXVG, "-o",  "order",  ffWRITE },
    { efPDB, "-p",  "order",  ffOPTWR },
    { efDAT, "-ss", "ssdump", ffOPTRD },
    { efXVG, "-jc", "Jcoupling", ffWRITE },
    { efXVG, "-corr",  "dihcorr",ffOPTWR },
    { efLOG, "-g",  "chi",    ffWRITE }
  };
#define NFILE asize(fnm)
  int     npargs;
  t_pargs *ppa;

  CopyRight(stderr,argv[0]);
  npargs = asize(pa);
  ppa    = add_acf_pargs(&npargs,pa);
  parse_common_args(&argc,argv,PCA_CAN_VIEW | PCA_CAN_TIME,TRUE,
		    NFILE,fnm,npargs,ppa,asize(desc),desc,asize(bugs),bugs);

  /* Handle result from enumerated type */
  sscanf(maxchistr[0],"%d",&maxchi);
  bChi = (maxchi > 0);
  
  log=ffopen(ftp2fn(efLOG,NFILE,fnm),"w");

  if (bRamOmega) {
    bOmega = TRUE;
    bPhi   = TRUE;
    bPsi   = TRUE;
  }
    
  bCorr=(opt2bSet("-corr",NFILE,fnm));
  if (bCorr) 
    fprintf(stderr,"Will calculate autocorrelation\n");
  
  if (maxchi > MAXCHI) {
    fprintf(stderr, 
	    "Will only calculate first %d Chi dihedrals in stead of %d.\n",
	    MAXCHI, maxchi);
    maxchi=MAXCHI;
  }
  bSSHisto = ftp2bSet(efDAT,NFILE,fnm);
  
  /* Find the chi angles using atoms struct and a list of amino acids */
  get_stx_coordnum(ftp2fn(efSTX,NFILE,fnm),&natoms);
  init_t_atoms(&atoms,natoms,TRUE);
  snew(x,natoms);
  read_stx_conf(ftp2fn(efSTX,NFILE,fnm),title,&atoms,x,NULL,box);
  fprintf(log,"Title: %s\n",title);
  
  naa=get_strings("aminoacids.dat",&aa);
  dlist=mk_dlist(log,&atoms,&nlist,bPhi,bPsi,bChi,maxchi,r0,naa,aa);
  fprintf(stderr,"%d residues with dihedrals found\n", nlist);
  
  if (nlist == 0) 
    fatal_error(0,"No dihedrals in your structure!\n");
  
  /* Make a linear index for reading all */
  index=make_chi_ind(nlist,dlist,&ndih);
  isize=4*ndih;
  fprintf(stderr,"%d dihedrals found\n", ndih);

  snew(dih,ndih);
    
  /* COMPUTE ALL DIHEDRALS! */
  read_ang_dih(ftp2fn(efTRX,NFILE,fnm),ftp2fn(efSTX,NFILE,fnm),
	       FALSE,TRUE,FALSE,1,&idum,
	       &nf,&time,isize,index,&trans_frac,&aver_angle,dih);

  if (bCorr) {
    if (nf < 2)
      fatal_error(0,"Need at least 2 frames for correlation");
    
    dt=(time[nf-1]-time[0])/(nf-1);
  }

  /* put angles in -M_PI to M_PI ! and correct phase factor for phi and psi */
  reset_em_all(nlist,dlist,nf,dih,maxchi,bPhi,bPsi,bChi);
  
  if (bAll)
    dump_em_all(nlist,dlist,nf,time,dih,maxchi,bPhi,bPsi,bChi,bOmega);
  
  /* Histogramming & J coupling constants */
  histogramming(log,naa,aa,nf,maxchi,dih,nlist,dlist,index,
		bPhi,bPsi,bOmega,bChi,
		bNormHisto,bSSHisto,ftp2fn(efDAT,NFILE,fnm),bfac_max,&atoms);

  /* Order parameters */  
  order_params(log,opt2fn("-o",NFILE,fnm),maxchi,nlist,dlist,
	       ftp2fn_null(efPDB,NFILE,fnm),bfac_init,
	       &atoms,x,box,bPhi,bPsi,bChi);
  
  /* Print ramachandran maps! */
  if (bRama)
    do_rama(nf,nlist,dlist,dih,bViol,bRamOmega);
  
  if (bShift)
    do_pp2shifts(log,nf,nlist,dlist,dih);
  
  pr_dlist(log,nlist,dlist,time[nf-1]-time[0]);
  ffclose(log);
  
  /* Correlation comes last because it fucks up the angles */
  if (bCorr)
    do_dihcorr(opt2fn("-corr",NFILE,fnm),nf,ndih,dih,dt,nlist,dlist,time,
	       maxchi,bPhi,bPsi,bChi,bOmega);
  
  
  xvgr_file(opt2fn("-o",NFILE,fnm),"-nxy");
  xvgr_file(opt2fn("-jc",NFILE,fnm),"-nxy");
  if (bCorr)
    xvgr_file(opt2fn("-corr",NFILE,fnm),"-nxy");
    
  thanx(stderr);
    
  return 0;
}
