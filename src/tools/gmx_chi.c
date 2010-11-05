/*
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
 * Green Red Orange Magenta Azure Cyan Skyblue
 */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include <stdio.h>
#include <math.h>

#include "confio.h"
#include "pdbio.h"
#include "copyrite.h"
#include "gmx_fatal.h"
#include "futil.h"
#include "gstat.h"
#include "macros.h"
#include "maths.h"
#include "physics.h"
#include "index.h"
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
#include "matio.h"
#include "gmx_ana.h"

static gmx_bool bAllowed(real phi,real psi)
{
  static const char *map[] = {
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
  return (gmx_bool) map[x][y];
}

atom_id *make_chi_ind(int nl,t_dlist dl[],int *ndih)
{
  atom_id *id;
  int     i,Xi,n;
  
  /* There are nl residues with max edMax dihedrals with 4 atoms each */
  snew(id,nl*edMax*4); 
  
  n=0;
  for(i=0; (i<nl); i++) 
  {
	  /* Phi, fake the first one */
	  dl[i].j0[edPhi] = n/4;
	  if(dl[i].atm.minC >= 0)
		  id[n++]=dl[i].atm.minC;
	  else
		  id[n++]=dl[i].atm.H;
	  id[n++]=dl[i].atm.N;
	  id[n++]=dl[i].atm.Cn[1];
	  id[n++]=dl[i].atm.C;
  }
  for(i=0; (i<nl); i++) 
  { 
	  /* Psi, fake the last one */
	  dl[i].j0[edPsi] = n/4;
	  id[n++]=dl[i].atm.N;
	  id[n++]=dl[i].atm.Cn[1];
	  id[n++]=dl[i].atm.C;
	  if ( i< (nl-1) )
		  id[n++]=dl[i+1].atm.N;
	  else
		  id[n++]=dl[i].atm.O;  
  }
  for(i=0; (i<nl); i++) 
  {
	  /* Omega */
	  if (has_dihedral(edOmega,&(dl[i])))
	  {
		  dl[i].j0[edOmega] = n/4;
		  id[n++]=dl[i].atm.minO;
		  id[n++]=dl[i].atm.minC;
		  id[n++]=dl[i].atm.N;
		  id[n++]=dl[i].atm.H;
	  }
  }
  for(Xi=0; (Xi<MAXCHI); Xi++)
  { 
	  /* Chi# */
	  for(i=0; (i<nl); i++) 
	  {
		  if (dl[i].atm.Cn[Xi+3] != -1) 
		  {
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


static void do_dihcorr(const char *fn,int nf,int ndih,real **dih,real dt,
		       int nlist,t_dlist dlist[],real time[],int maxchi,
		       gmx_bool bPhi,gmx_bool bPsi,gmx_bool bChi,gmx_bool bOmega,
                       const output_env_t oenv)
{
  char name1[256],name2[256];
  int  i,j,Xi;
  
  do_autocorr(fn,oenv,"Dihedral Autocorrelation Function",
	      nf,ndih,dih,dt,eacCos,FALSE);
  /* Dump em all */
  j=0;
  for(i=0; (i<nlist); i++) {
    if (bPhi)
      print_one(oenv,"corrphi",dlist[i].name,"Phi ACF for", "C(t)", nf/2,time,
                          dih[j]);
    j++;
  }
  for(i=0; (i<nlist); i++) {
    if (bPsi)
      print_one(oenv,"corrpsi",dlist[i].name,"Psi ACF for","C(t)",nf/2,time,
                dih[j]);
    j++;
  }
  for(i=0; (i<nlist); i++) {
    if (has_dihedral(edOmega,&dlist[i])) {
      if (bOmega)
	print_one(oenv,"corromega",dlist[i].name,"Omega ACF for","C(t)",
                  nf/2,time,dih[j]);
      j++;
    }
  }
  for(Xi=0; (Xi<maxchi); Xi++) {
    sprintf(name1, "corrchi%d", Xi+1);
    sprintf(name2, "Chi%d ACF for", Xi+1);
    for(i=0; (i<nlist); i++) {
      if (dlist[i].atm.Cn[Xi+3] != -1) {
	if (bChi)
	  print_one(oenv,name1,dlist[i].name,name2,"C(t)",nf/2,time,dih[j]);
	j++;
      }
    }
  }
  fprintf(stderr,"\n");
}

static void copy_dih_data(real in[], real out[], int nf, gmx_bool bLEAVE)
{
  /* if bLEAVE, do nothing to data in copying to out
   * otherwise multiply by 180/pi to convert rad to deg */ 
  int i ;
  real mult ; 
  if (bLEAVE)
    mult = 1 ; 
  else
    mult = (180.0/M_PI); 
  for (i=0;(i<nf);i++){
    out[i]=in[i]*mult ; 
  }
}

static void dump_em_all(int nlist,t_dlist dlist[],int nf,real time[],
			real **dih,int maxchi,
			gmx_bool bPhi,gmx_bool bPsi,gmx_bool bChi,gmx_bool bOmega, gmx_bool bRAD,
                        const output_env_t oenv)
{
  char name[256], titlestr[256], ystr[256]; 
  real *data ; 
  int  i,j,Xi;
  
  snew(data,nf); 
  if (bRAD) 
    strcpy(ystr,"Angle (rad)"); 
  else
    strcpy(ystr,"Angle (degrees)"); 
    
  /* Dump em all */
  j = 0;
  for(i=0; (i<nlist); i++) {
      /* grs debug  printf("OK i %d j %d\n", i, j) ; */
    if (bPhi) {
      copy_dih_data(dih[j],data,nf,bRAD); 
      print_one(oenv,"phi",dlist[i].name,"\\xf\\f{}",ystr, nf,time,data); 
    }
    j++;
  }
  for(i=0; (i<nlist); i++) {
    if (bPsi) {
      copy_dih_data(dih[j],data,nf,bRAD); 
      print_one(oenv,"psi",dlist[i].name,"\\xy\\f{}",ystr, nf,time,data);
    }
    j++;
  }  
  for(i=0; (i<nlist); i++)
    if (has_dihedral(edOmega,&(dlist[i]))) {
      if (bOmega){
	copy_dih_data(dih[j],data,nf,bRAD); 
	print_one(oenv,"omega",dlist[i].name,"\\xw\\f{}",ystr,nf,time,data);
      }
      j++;
    }
  
  for(Xi=0; (Xi<maxchi); Xi++)
    for(i=0; (i<nlist); i++)
      if (dlist[i].atm.Cn[Xi+3] != -1) {
	if (bChi) {
	  sprintf(name,"chi%d",Xi+1);
	  sprintf(titlestr,"\\xc\\f{}\\s%d\\N",Xi+1);
	  copy_dih_data(dih[j],data,nf,bRAD); 
	  print_one(oenv,name,dlist[i].name,titlestr,ystr, nf,time,data);
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

static int reset_em_all(int nlist,t_dlist dlist[],int nf,
			 real **dih,int maxchi)
{
  int  i,j,Xi;
  
  /* Reset em all */
  j=0;
  /* Phi */
  for(i=0; (i<nlist); i++)
  {
	  if (dlist[i].atm.minC == -1)
	  {  
		  reset_one(dih[j++],nf,M_PI);
	  }
	  else
      {
		  reset_one(dih[j++],nf,0);
	  }
  }
  /* Psi */
  for(i=0; (i<nlist-1); i++)
  {
	  reset_one(dih[j++],nf,0);		  
  }	  
  /* last Psi is faked from O */
  reset_one(dih[j++],nf,M_PI);		  
  
  /* Omega */
  for(i=0; (i<nlist); i++)
	  if (has_dihedral(edOmega,&dlist[i]))
		  reset_one(dih[j++],nf,0);
  /* Chi 1 thru maxchi */
  for(Xi=0; (Xi<maxchi); Xi++)
  {
	  for(i=0; (i<nlist); i++)
	  {
		  if (dlist[i].atm.Cn[Xi+3] != -1) 
		  {
			  reset_one(dih[j],nf,0);
			  j++;
		  }
	  }
  }
  fprintf(stderr,"j after resetting (nr. active dihedrals) = %d\n",j);
  return j ; 
}

static void histogramming(FILE *log,int nbin,gmx_residuetype_t rt,
			  int nf,int maxchi,real **dih,
			  int nlist,t_dlist dlist[],
			  atom_id index[],
			  gmx_bool bPhi,gmx_bool bPsi,gmx_bool bOmega,gmx_bool bChi,
			  gmx_bool bNormalize,gmx_bool bSSHisto,const char *ssdump,
			  real bfac_max,t_atoms *atoms, 
			  gmx_bool bDo_jc, const char *fn,
                          const output_env_t oenv)
{
  /* also gets 3J couplings and order parameters S2 */ 
  t_karplus kkkphi[] = {
    { "J_NHa1",    6.51, -1.76,  1.6, -M_PI/3,   0.0,  0.0 },
    { "J_NHa2",    6.51, -1.76,  1.6,  M_PI/3,   0.0,  0.0 },
    { "J_HaC'",    4.0,   1.1,   0.1,  0.0,      0.0,  0.0 },
    { "J_NHCb",    4.7,  -1.5,  -0.2,  M_PI/3,   0.0,  0.0 },
    { "J_Ci-1Hai", 4.5,  -1.3,  -1.2,  2*M_PI/3, 0.0,  0.0 }
  };
  t_karplus kkkpsi[] = {
    { "J_HaN",   -0.88, -0.61,-0.27,M_PI/3,  0.0,  0.0 }
  };
  t_karplus kkkchi1[] = {
    { "JHaHb2",       9.5, -1.6, 1.8, -M_PI/3, 0,  0.0 },
    { "JHaHb3",       9.5, -1.6, 1.8, 0, 0,  0.0 }
  };
#define NKKKPHI asize(kkkphi)
#define NKKKPSI asize(kkkpsi)
#define NKKKCHI asize(kkkchi1)
#define NJC (NKKKPHI+NKKKPSI+NKKKCHI)
  
  FILE    *fp,*ssfp[3]={NULL,NULL,NULL};
  const char *sss[3] = { "sheet", "helix", "coil" };
  real    S2;
  real    *normhisto;
  real    **Jc,**Jcsig;
  int     ****his_aa_ss=NULL;
  int     ***his_aa,**his_aa1,*histmp;
  int     i,j,k,m,n,nn,Dih,nres,hindex,angle;
  gmx_bool    bBfac,bOccup;
  char    hisfile[256],hhisfile[256],sshisfile[256],title[256],*ss_str=NULL;
  char **leg; 
  const char *residue_name;
  int     rt_size;

  rt_size = gmx_residuetype_get_size(rt);
  if (bSSHisto) {
    fp = ffopen(ssdump,"r");
    if(1 != fscanf(fp,"%d",&nres))
    {
      gmx_fatal(FARGS,"Error reading from file %s",ssdump);
    }

    snew(ss_str,nres+1);
    if( 1 != fscanf(fp,"%s",ss_str))
    {
      gmx_fatal(FARGS,"Error reading from file %s",ssdump);
    }

    ffclose(fp);
    /* Four dimensional array... Very cool */
    snew(his_aa_ss,3);
    for(i=0; (i<3); i++) {
      snew(his_aa_ss[i],rt_size+1);
      for(j=0; (j<=rt_size); j++) {
	snew(his_aa_ss[i][j],edMax);
	for(Dih=0; (Dih<edMax); Dih++)
	  snew(his_aa_ss[i][j][Dih],nbin+1);
      }
    }
  }
  snew(his_aa,edMax);
  for(Dih=0; (Dih<edMax); Dih++) {
    snew(his_aa[Dih],rt_size+1);
    for(i=0; (i<=rt_size); i++) {
      snew(his_aa[Dih][i],nbin+1);
    }
  }
  snew(histmp,nbin);
  
  snew(Jc,nlist);
  snew(Jcsig,nlist);
  for(i=0; (i<nlist); i++) {
    snew(Jc[i],NJC);
    snew(Jcsig[i],NJC);
  }
  
  j=0;
  n=0;
  for (Dih=0; (Dih<NONCHI+maxchi); Dih++) {    
    for(i=0; (i<nlist); i++) {
      if (((Dih  < edOmega) ) ||
	  ((Dih == edOmega) && (has_dihedral(edOmega,&(dlist[i])))) ||
	  ((Dih  > edOmega) && (dlist[i].atm.Cn[Dih-NONCHI+3] != -1))) {
      	make_histo(log,nf,dih[j],nbin,histmp,-M_PI,M_PI);
	
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
	    hindex = ((dih[j][0]+M_PI)*nbin)/(2*M_PI);
	    range_check(hindex,0,nbin);
	    
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
	  calc_distribution_props(nbin,histmp,-M_PI,NKKKPHI,kkkphi,&S2);
	  
	  for(m=0; (m<NKKKPHI); m++) {
	    Jc[i][m]    = kkkphi[m].Jc;
	    Jcsig[i][m] = kkkphi[m].Jcsig;
	  }
	  break;
	case edPsi:
	  calc_distribution_props(nbin,histmp,-M_PI,NKKKPSI,kkkpsi,&S2);
	  
	  for(m=0; (m<NKKKPSI); m++) {
	    Jc[i][NKKKPHI+m]    = kkkpsi[m].Jc;
	    Jcsig[i][NKKKPHI+m] = kkkpsi[m].Jcsig;
	  }
	  break;
	case edChi1:
	  calc_distribution_props(nbin,histmp,-M_PI,NKKKCHI,kkkchi1,&S2);
	  for(m=0; (m<NKKKCHI); m++) {
	    Jc[i][NKKKPHI+NKKKPSI+m]    = kkkchi1[m].Jc;
	    Jcsig[i][NKKKPHI+NKKKPSI+m] = kkkchi1[m].Jcsig;
	  }
	  break;
	default: /* covers edOmega and higher Chis than Chi1 */ 
	  calc_distribution_props(nbin,histmp,-M_PI,0,NULL,&S2);
	  break;
	}
	dlist[i].S2[Dih]        = S2;
	
	/* Sum distribution per amino acid type as well */
	for(k=0; (k<nbin); k++) {
	  his_aa[Dih][dlist[i].index][k] += histmp[k];
	  histmp[k] = 0;
	}
	j++;
      } else { /* dihed not defined */
	dlist[i].S2[Dih] = 0.0 ; 
      }
    }
  }
  sfree(histmp);
  
  /* Print out Jcouplings */
  fprintf(log,"\n *** J-Couplings from simulation (plus std. dev.) ***\n\n");
  fprintf(log,"Residue   ");
  for(i=0; (i<NKKKPHI); i++)
    fprintf(log,"%7s   SD",kkkphi[i].name);
  for(i=0; (i<NKKKPSI); i++)
    fprintf(log,"%7s   SD",kkkpsi[i].name);
  for(i=0; (i<NKKKCHI); i++)
    fprintf(log,"%7s   SD",kkkchi1[i].name);
  fprintf(log,"\n");
  for(i=0; (i<NJC+1); i++)
    fprintf(log,"------------");
  fprintf(log,"\n");
  for(i=0; (i<nlist); i++) {
    fprintf(log,"%-10s",dlist[i].name);
    for(j=0; (j<NJC); j++)
      fprintf(log,"  %5.2f %4.2f",Jc[i][j],Jcsig[i][j]);
    fprintf(log,"\n");
  }
  fprintf(log,"\n");

  /* and to -jc file... */ 
  if (bDo_jc) {
    fp=xvgropen(fn,"\\S3\\NJ-Couplings from Karplus Equation","Residue",
                "Coupling",oenv); 
    snew(leg,NJC); 
    for(i=0; (i<NKKKPHI); i++){
		leg[i] = strdup(kkkphi[i].name); 
    }
    for(i=0; (i<NKKKPSI); i++){
		leg[i+NKKKPHI]=strdup(kkkpsi[i].name); 
    }
    for(i=0; (i<NKKKCHI); i++){
      leg[i+NKKKPHI+NKKKPSI]=strdup(kkkchi1[i].name); 
    }      
    xvgr_legend(fp,NJC,(const char**)leg,oenv);
    fprintf(fp,"%5s ","#Res.");
    for(i=0; (i<NJC); i++)
      fprintf(fp,"%10s ",leg[i]); 
    fprintf(fp,"\n"); 
    for(i=0; (i<nlist); i++) {
      fprintf(fp,"%5d ",dlist[i].resnr);
      for(j=0; (j<NJC); j++)
	fprintf(fp,"  %8.3f",Jc[i][j]);
      fprintf(fp,"\n"); 
    }
    ffclose(fp);
    for(i=0; (i<NJC); i++)
      sfree(leg[i]); 
  }
  /* finished -jc stuff */ 

  snew(normhisto,nbin);
  for(i=0; (i<rt_size); i++) {
    for(Dih=0; (Dih<edMax); Dih++){
      /* First check whether something is in there */
      for(j=0; (j<nbin); j++)
	if (his_aa[Dih][i][j] != 0)
	  break;
      if ((j < nbin) &&
	  ((bPhi && (Dih==edPhi)) ||
	   (bPsi && (Dih==edPsi)) ||
	   (bOmega &&(Dih==edOmega)) ||
	   (bChi && (Dih>=edChi1)))) {
	if (bNormalize)
	  normalize_histo(nbin,his_aa[Dih][i],(360.0/nbin),normhisto);
	
	residue_name = gmx_residuetype_get_name(rt,i);
	switch (Dih) {
	case edPhi:
	  sprintf(hisfile,"histo-phi%s",residue_name);
	  sprintf(title,"\\xf\\f{} Distribution for %s",residue_name);
	  break;
	case edPsi:
	  sprintf(hisfile,"histo-psi%s",residue_name);
	  sprintf(title,"\\xy\\f{} Distribution for %s",residue_name);
	  break;
	case edOmega:
	  sprintf(hisfile,"histo-omega%s",residue_name);
	  sprintf(title,"\\xw\\f{} Distribution for %s",residue_name);
	  break;
	default:
	  sprintf(hisfile,"histo-chi%d%s",Dih-NONCHI+1,residue_name);
	  sprintf(title,"\\xc\\f{}\\s%d\\N Distribution for %s",
		  Dih-NONCHI+1,residue_name);
	}
	strcpy(hhisfile,hisfile);
	strcat(hhisfile,".xvg");
	fp=xvgropen(hhisfile,title,"Degrees","",oenv);
	fprintf(fp,"@ with g0\n");
	xvgr_world(fp,-180,0,180,0.1,oenv);
	fprintf(fp,"# this effort to set graph size fails unless you run with -autoscale none or -autoscale y flags\n"); 
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
	for(j=0; (j<nbin); j++) {
	  angle = -180 + (360/nbin)*j ; 
	  if (bNormalize)
	    fprintf(fp,"%5d  %10g\n",angle,normhisto[j]);
	  else
	    fprintf(fp,"%5d  %10d\n",angle,his_aa[Dih][i][j]);
	  if (bSSHisto)
	    for(k=0; (k<3); k++) 
	      fprintf(ssfp[k],"%5d  %10d\n",angle,
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
      for(j=0; (j<=rt_size); j++) {
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

static FILE *rama_file(const char *fn,const char *title,const char *xaxis,
                       const char *yaxis,const output_env_t oenv)
{
  FILE *fp;

  fp = xvgropen(fn,title,xaxis,yaxis,oenv);  
  fprintf(fp,"@ with g0\n");
  xvgr_world(fp,-180,-180,180,180,oenv);
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
		    gmx_bool bViol,gmx_bool bRamOmega,const output_env_t oenv)
{
  FILE *fp,*gp=NULL;
  gmx_bool bOm;
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
		     "\\8f\\4 (deg)","\\8y\\4 (deg)",oenv);
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
	ffclose(gp);
      ffclose(fp);
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
	write_xpm3(fp,0,"Omega/Ramachandran Plot","Deg","Phi","Psi",
		   NMAT,NMAT,axis,axis,mat,lo,180.0,hi,rlo,rmid,rhi,&nlevels);
	ffclose(fp);
	for(j=0; (j<NMAT); j++)
	  sfree(mat[j]);
	sfree(mat);
      }
    }
    if ((has_dihedral(edChi1,&(dlist[i]))) &&
	(has_dihedral(edChi2,&(dlist[i])))) {
      sprintf(fn,"ramaX1X2%s.xvg",dlist[i].name);
      fp = rama_file(fn,"\\8c\\4\\s1\\N-\\8c\\4\\s2\\N Ramachandran Plot",
		     "\\8c\\4\\s1\\N (deg)","\\8c\\4\\s2\\N (deg)",oenv);
      Xi1 = dlist[i].j0[edChi1];
      Xi2 = dlist[i].j0[edChi2];
      for(j=0; (j<nf); j++)
	fprintf(fp,"%10g  %10g\n",RAD2DEG*dih[Xi1][j],RAD2DEG*dih[Xi2][j]);
      ffclose(fp);
    }
    else 
      fprintf(stderr,"No chi1 & chi2 angle for %s\n",dlist[i].name);
  }
}


static void print_transitions(const char *fn,int maxchi,int nlist,
                              t_dlist dlist[], t_atoms *atoms,rvec x[],
                              matrix box, gmx_bool bPhi,gmx_bool bPsi,gmx_bool bChi,real dt,
                              const output_env_t oenv)
{
  /* based on order_params below */ 
  FILE *fp;
  int  nh[edMax];
  int  i,Dih,Xi;
	
  /*  must correspond with enum in pp2shift.h:38 */  
  char *leg[edMax];
#define NLEG asize(leg) 

  leg[0] = strdup("Phi");
  leg[1] = strdup("Psi");
  leg[2] = strdup("Omega");
  leg[3] = strdup("Chi1");
  leg[4] = strdup("Chi2");
  leg[5] = strdup("Chi3");
  leg[6] = strdup("Chi4");
  leg[7] = strdup("Chi5");
  leg[8] = strdup("Chi6");
  
  /* Print order parameters */
  fp=xvgropen(fn,"Dihedral Rotamer Transitions","Residue","Transitions/ns",
                oenv);
  xvgr_legend(fp,NONCHI+maxchi,(const char**)leg,oenv);
  
  for (Dih=0; (Dih<edMax); Dih++)
    nh[Dih]=0;
  
  fprintf(fp,"%5s ","#Res.");
  fprintf(fp,"%10s %10s %10s ",leg[edPhi],leg[edPsi],leg[edOmega]);
  for (Xi=0; Xi<maxchi; Xi++)
    fprintf(fp,"%10s ",leg[NONCHI+Xi]);
  fprintf(fp,"\n"); 
  
  for(i=0; (i<nlist); i++) {
    fprintf(fp,"%5d ",dlist[i].resnr);
    for (Dih=0; (Dih<NONCHI+maxchi); Dih++)
      fprintf(fp,"%10.3f ",dlist[i].ntr[Dih]/dt);
    /* fprintf(fp,"%12s\n",dlist[i].name);  this confuses xmgrace */ 
    fprintf(fp,"\n"); 
  }
  ffclose(fp);
}

static void order_params(FILE *log,
			 const char *fn,int maxchi,int nlist,t_dlist dlist[],
			 const char *pdbfn,real bfac_init,
			 t_atoms *atoms,rvec x[],int ePBC,matrix box,
			 gmx_bool bPhi,gmx_bool bPsi,gmx_bool bChi,const output_env_t oenv)
{
  FILE *fp;
  int  nh[edMax];
  char buf[STRLEN];
  int  i,Dih,Xi;
  real S2Max, S2Min;

  /* except for S2Min/Max, must correspond with enum in pp2shift.h:38 */  
  const char *const_leg[2+edMax]= { "S2Min","S2Max","Phi","Psi","Omega", 
                                    "Chi1", "Chi2", "Chi3", "Chi4", "Chi5", 
                                    "Chi6" };
#define NLEG asize(leg) 
  
  char *leg[2+edMax];	
	
  for(i=0;i<NLEG;i++)
    leg[i]=strdup(const_leg[i]);
	
  /* Print order parameters */
  fp=xvgropen(fn,"Dihedral Order Parameters","Residue","S2",oenv);
  xvgr_legend(fp,2+NONCHI+maxchi,const_leg,oenv);
  
  for (Dih=0; (Dih<edMax); Dih++)
    nh[Dih]=0;
  
  fprintf(fp,"%5s ","#Res.");
  fprintf(fp,"%10s %10s ",leg[0],leg[1]);
  fprintf(fp,"%10s %10s %10s ",leg[2+edPhi],leg[2+edPsi],leg[2+edOmega]);
  for (Xi=0; Xi<maxchi; Xi++)
    fprintf(fp,"%10s ",leg[2+NONCHI+Xi]);
  fprintf(fp,"\n"); 
  
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
    fprintf(fp,"%5d ",dlist[i].resnr);
    fprintf(fp,"%10.3f %10.3f ",S2Min,S2Max);
    for (Dih=0; (Dih<NONCHI+maxchi); Dih++)
      fprintf(fp,"%10.3f ",dlist[i].S2[Dih]);
    fprintf(fp,"\n"); 
    /* fprintf(fp,"%12s\n",dlist[i].name);  this confuses xmgrace */ 
  }
  ffclose(fp);
  
  if (NULL != pdbfn) {
    real x0,y0,z0;

    if (NULL == atoms->pdbinfo)
      snew(atoms->pdbinfo,atoms->nr);
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
    write_pdbfile(fp,NULL,atoms,x,ePBC,box,' ',0,NULL,TRUE);
    x0=y0=z0=1000.0;
    for (i=0; (i<atoms->nr); i++) {
      x0 = min(x0, x[i][XX]);
      y0 = min(y0, x[i][YY]);
      z0 = min(z0, x[i][ZZ]);
    }
    x0*=10.0;/* nm -> angstrom */
    y0*=10.0;/* nm -> angstrom */
    z0*=10.0;/* nm -> angstrom */
    sprintf(buf,"%s%%6.f%%6.2f\n",pdbformat);
    for (i=0; (i<10); i++) {
      fprintf(fp,buf,"ATOM  ", atoms->nr+1+i, "CA", "LEG",' ', 
	      atoms->nres+1, ' ',x0, y0, z0+(1.2*i), 0.0, -0.1*i);
    }
    ffclose(fp);
  }
  
  fprintf(log,"Dihedrals with S2 > 0.8\n");
  fprintf(log,"Dihedral: ");
  if (bPhi) fprintf(log," Phi  ");
  if (bPsi) fprintf(log," Psi ");
  if (bChi)
    for(Xi=0; (Xi<maxchi); Xi++)
      fprintf(log," %s ", leg[2+NONCHI+Xi]);
  fprintf(log,"\nNumber:   ");
  if (bPhi) fprintf(log,"%4d  ",nh[0]);
  if (bPsi) fprintf(log,"%4d  ",nh[1]);
  if (bChi)
    for(Xi=0; (Xi<maxchi); Xi++)
      fprintf(log,"%4d  ",nh[NONCHI+Xi]);
  fprintf(log,"\n");
	
  for(i=0;i<NLEG;i++)
    sfree(leg[i]);

}

int gmx_chi(int argc,char *argv[])
{
  const char *desc[] = {
    "g_chi computes phi, psi, omega and chi dihedrals for all your ",
    "amino acid backbone and sidechains.",
    "It can compute dihedral angle as a function of time, and as",
    "histogram distributions.",
    "The distributions (histo-(dihedral)(RESIDUE).xvg) are cumulative over all residues of each type.[PAR]", 
    "If option [TT]-corr[tt] is given, the program will",
    "calculate dihedral autocorrelation functions. The function used",
    "is C(t) = < cos(chi(tau)) cos(chi(tau+t)) >. The use of cosines",
    "rather than angles themselves, resolves the problem of periodicity.",
    "(Van der Spoel & Berendsen (1997), [BB]Biophys. J. 72[bb], 2032-2041).",
    "Separate files for each dihedral of each residue", 
    "(corr(dihedral)(RESIDUE)(nresnr).xvg) are output, as well as a", 
    "file containing the information for all residues (argument of [TT]-corr[tt]).[PAR]", 
    "With option [TT]-all[tt], the angles themselves as a function of time for", 
    "each residue are printed to separate files (dihedral)(RESIDUE)(nresnr).xvg.", 
    "These can be in radians or degrees.[PAR]", 
    "A log file (argument [TT]-g[tt]) is also written. This contains [BR]",
    "(a) information about the number of residues of each type.[BR]", 
    "(b) The NMR 3J coupling constants from the Karplus equation.[BR]", 
    "(c) a table for each residue of the number of transitions between ", 
    "rotamers per nanosecond,  and the order parameter S2 of each dihedral.[BR]",
    "(d) a table for each residue of the rotamer occupancy.[BR]", 
    "All rotamers are taken as 3-fold, except for omegas and chi-dihedrals",
    "to planar groups (i.e. chi2 of aromatics asp and asn, chi3 of glu", 
    "and gln, and chi4 of arg), which are 2-fold. \"rotamer 0\" means ", 
    "that the dihedral was not in the core region of each rotamer. ", 
    "The width of the core region can be set with [TT]-core_rotamer[tt][PAR]", 

    "The S2 order parameters are also output to an xvg file", 
    "(argument [TT]-o[tt] ) and optionally as a pdb file with", 
    "the S2 values as B-factor (argument [TT]-p[tt]). ", 
    "The total number of rotamer transitions per timestep", 
    "(argument [TT]-ot[tt]), the number of transitions per rotamer", 
    "(argument [TT]-rt[tt]), and the 3J couplings (argument [TT]-jc[tt]), ", 
    "can also be written to .xvg files.[PAR]", 

    "If [TT]-chi_prod[tt] is set (and maxchi > 0), cumulative rotamers, e.g.", 
    "1+9(chi1-1)+3(chi2-1)+(chi3-1) (if the residue has three 3-fold ", 
    "dihedrals and maxchi >= 3)", 
    "are calculated. As before, if any dihedral is not in the core region,", 
    "the rotamer is taken to be 0. The occupancies of these cumulative ",
    "rotamers (starting with rotamer 0) are written to the file", 
    "that is the argument of [TT]-cp[tt], and if the [TT]-all[tt] flag", 
    "is given, the rotamers as functions of time", 
    "are written to chiproduct(RESIDUE)(nresnr).xvg ", 
    "and their occupancies to histo-chiproduct(RESIDUE)(nresnr).xvg.[PAR]", 

    "The option [TT]-r[tt] generates a contour plot of the average omega angle",
    "as a function of the phi and psi angles, that is, in a Ramachandran plot",
    "the average omega angle is plotted using color coding.", 

  };
  
  const char *bugs[] = {
    "Produces MANY output files (up to about 4 times the number of residues in the protein, twice that if autocorrelation functions are calculated). Typically several hundred files are output.",
    "Phi and psi dihedrals are calculated in a non-standard way, using H-N-CA-C for phi instead of C(-)-N-CA-C, and N-CA-C-O for psi instead of N-CA-C-N(+). This causes (usually small) discrepancies with the output of other tools like g_rama.", 
    "-r0 option does not work properly", 
    "Rotamers with multiplicity 2 are printed in chi.log as if they had multiplicity 3, with the 3rd (g(+)) always having probability 0" 
  };

  /* defaults */ 
  static int  r0=1,ndeg=1,maxchi=2;
  static gmx_bool bAll=FALSE;
  static gmx_bool bPhi=FALSE,bPsi=FALSE,bOmega=FALSE;
  static real bfac_init=-1.0,bfac_max=0;
  static const char *maxchistr[] = { NULL, "0", "1", "2", "3",  "4", "5", "6", NULL };
  static gmx_bool bRama=FALSE,bShift=FALSE,bViol=FALSE,bRamOmega=FALSE;
  static gmx_bool bNormHisto=TRUE,bChiProduct=FALSE,bHChi=FALSE,bRAD=FALSE,bPBC=TRUE;
  static real core_frac=0.5 ;  
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
    { "-periodic", FALSE, etBOOL, {&bPBC},
      "Print dihedral angles modulo 360 degrees" },
    { "-all",  FALSE, etBOOL, {&bAll},
      "Output separate files for every dihedral." },
    { "-rad",  FALSE, etBOOL, {&bRAD},
      "in angle vs time files, use radians rather than degrees."}, 
    { "-shift", FALSE, etBOOL, {&bShift},
	"Compute chemical shifts from Phi/Psi angles" },
    { "-binwidth", FALSE, etINT, {&ndeg},
      "bin width for histograms (degrees)" },
    { "-core_rotamer", FALSE, etREAL, {&core_frac},
      "only the central -core_rotamer*(360/multiplicity) belongs to each rotamer (the rest is assigned to rotamer 0)" },
    { "-maxchi", FALSE, etENUM, {maxchistr},
      "calculate first ndih Chi dihedrals" },
    { "-normhisto", FALSE, etBOOL, {&bNormHisto},
      "Normalize histograms" },
    { "-ramomega",FALSE,etBOOL, {&bRamOmega},
      "compute average omega as a function of phi/psi and plot it in an xpm plot" },
    { "-bfact", FALSE, etREAL, {&bfac_init},
      "B-factor value for pdb file for atoms with no calculated dihedral order parameter"},
    { "-chi_prod",FALSE,etBOOL, {&bChiProduct},
      "compute a single cumulative rotamer for each residue"},
    { "-HChi",FALSE,etBOOL, {&bHChi},
      "Include dihedrals to sidechain hydrogens"}, 
    { "-bmax",  FALSE, etREAL, {&bfac_max},
      "Maximum B-factor on any of the atoms that make up a dihedral, for the dihedral angle to be considere in the statistics. Applies to database work where a number of X-Ray structures is analyzed. -bmax <= 0 means no limit." }
  };

  FILE       *log;
  int        natoms,nlist,idum,nbin;
  t_atoms    atoms;
  rvec       *x;
  int        ePBC;
  matrix     box;
  char       title[256],grpname[256]; 
  t_dlist    *dlist;
  gmx_bool       bChi,bCorr,bSSHisto;
  gmx_bool       bDo_rt, bDo_oh, bDo_ot, bDo_jc ; 
  real       dt=0, traj_t_ns;
  output_env_t oenv;
  gmx_residuetype_t rt;
  
  atom_id    isize,*index;
  int        ndih,nactdih,nf;
  real       **dih,*trans_frac,*aver_angle,*time;
  int        i,j,**chi_lookup,*xity; 
  
  t_filenm  fnm[] = {
    { efSTX, "-s",  NULL,     ffREAD  },
    { efTRX, "-f",  NULL,     ffREAD  },
    { efXVG, "-o",  "order",  ffWRITE },
    { efPDB, "-p",  "order",  ffOPTWR },
    { efDAT, "-ss", "ssdump", ffOPTRD },
    { efXVG, "-jc", "Jcoupling", ffWRITE },
    { efXVG, "-corr",  "dihcorr",ffOPTWR },
    { efLOG, "-g",  "chi",    ffWRITE },
    /* add two more arguments copying from g_angle */ 
    { efXVG, "-ot", "dihtrans", ffOPTWR }, 
    { efXVG, "-oh", "trhisto",  ffOPTWR },
    { efXVG, "-rt", "restrans",  ffOPTWR }, 
    { efXVG, "-cp", "chiprodhisto",  ffOPTWR }  
  };
#define NFILE asize(fnm)
  int     npargs;
  t_pargs *ppa;

  CopyRight(stderr,argv[0]);
  npargs = asize(pa);
  ppa    = add_acf_pargs(&npargs,pa);
  parse_common_args(&argc,argv,PCA_CAN_VIEW | PCA_CAN_TIME | PCA_BE_NICE,
		    NFILE,fnm,npargs,ppa,asize(desc),desc,asize(bugs),bugs,
                    &oenv);

  /* Handle result from enumerated type */
  sscanf(maxchistr[0],"%d",&maxchi);
  bChi = (maxchi > 0);
  
  log=ffopen(ftp2fn(efLOG,NFILE,fnm),"w");

  if (bRamOmega) {
    bOmega = TRUE;
    bPhi   = TRUE;
    bPsi   = TRUE;
  }
    
  /* set some options */ 
  bDo_rt=(opt2bSet("-rt",NFILE,fnm));
  bDo_oh=(opt2bSet("-oh",NFILE,fnm));
  bDo_ot=(opt2bSet("-ot",NFILE,fnm));
  bDo_jc=(opt2bSet("-jc",NFILE,fnm));
  bCorr=(opt2bSet("-corr",NFILE,fnm));
  if (bCorr) 
    fprintf(stderr,"Will calculate autocorrelation\n");
  
  if (core_frac > 1.0 ) {
    fprintf(stderr, "core_rotamer fraction > 1.0 ; will use 1.0\n"); 
    core_frac=1.0 ; 
  }
  if (core_frac < 0.0 ) {
    fprintf(stderr, "core_rotamer fraction < 0.0 ; will use 0.0\n"); 
    core_frac=0.0 ; 
  }

  if (maxchi > MAXCHI) {
    fprintf(stderr, 
	    "Will only calculate first %d Chi dihedrals in stead of %d.\n",
	    MAXCHI, maxchi);
    maxchi=MAXCHI;
  }
  bSSHisto = ftp2bSet(efDAT,NFILE,fnm);
  nbin = 360/ndeg ; 

  /* Find the chi angles using atoms struct and a list of amino acids */
  get_stx_coordnum(ftp2fn(efSTX,NFILE,fnm),&natoms);
  init_t_atoms(&atoms,natoms,TRUE);
  snew(x,natoms);
  read_stx_conf(ftp2fn(efSTX,NFILE,fnm),title,&atoms,x,NULL,&ePBC,box);
  fprintf(log,"Title: %s\n",title);
  
  gmx_residuetype_init(&rt);
  dlist=mk_dlist(log,&atoms,&nlist,bPhi,bPsi,bChi,bHChi,maxchi,r0,rt);
  fprintf(stderr,"%d residues with dihedrals found\n", nlist);
  
  if (nlist == 0) 
    gmx_fatal(FARGS,"No dihedrals in your structure!\n");
  
  /* Make a linear index for reading all. */
  index=make_chi_ind(nlist,dlist,&ndih);
  isize=4*ndih;
  fprintf(stderr,"%d dihedrals found\n", ndih);

  snew(dih,ndih);

  /* COMPUTE ALL DIHEDRALS! */
  read_ang_dih(ftp2fn(efTRX,NFILE,fnm),FALSE,TRUE,FALSE,bPBC,1,&idum,
	       &nf,&time,isize,index,&trans_frac,&aver_angle,dih,oenv);
  
  dt=(time[nf-1]-time[0])/(nf-1); /* might want this for corr or n. transit*/ 
  if (bCorr) 
  {
	  if (nf < 2)
	  {
		  gmx_fatal(FARGS,"Need at least 2 frames for correlation");
	  }
  }

  /* put angles in -M_PI to M_PI ! and correct phase factor for phi and psi 
  * pass nactdih instead of ndih to low_ana_dih_trans and get_chi_product_traj
  * to prevent accessing off end of arrays when maxchi < 5 or 6. */ 
  nactdih = reset_em_all(nlist,dlist,nf,dih,maxchi);
  
  if (bAll)
    dump_em_all(nlist,dlist,nf,time,dih,maxchi,bPhi,bPsi,bChi,bOmega,bRAD,oenv);
  
  /* Histogramming & J coupling constants & calc of S2 order params */
  histogramming(log,nbin,rt,nf,maxchi,dih,nlist,dlist,index,
		bPhi,bPsi,bOmega,bChi,
		bNormHisto,bSSHisto,ftp2fn(efDAT,NFILE,fnm),bfac_max,&atoms,
		bDo_jc,opt2fn("-jc",NFILE,fnm),oenv);

  /* transitions 
   *
   * added multiplicity */ 

  snew(xity,ndih) ;
  mk_multiplicity_lookup(xity, maxchi, dih, nlist, dlist,ndih); 
 
  strcpy(grpname, "All residues, "); 
  if(bPhi) 
    strcat(grpname, "Phi "); 
  if(bPsi) 
    strcat(grpname, "Psi "); 
  if(bOmega) 
    strcat(grpname, "Omega "); 
  if(bChi){ 
    strcat(grpname, "Chi 1-") ; 
    sprintf(grpname + strlen(grpname), "%i", maxchi); 
  }


  low_ana_dih_trans(bDo_ot, opt2fn("-ot",NFILE,fnm),
		    bDo_oh, opt2fn("-oh",NFILE,fnm),maxchi, 
		    dih, nlist, dlist, nf, nactdih, grpname, xity, 
		    *time,  dt, FALSE, core_frac,oenv) ; 

  /* Order parameters */  
  order_params(log,opt2fn("-o",NFILE,fnm),maxchi,nlist,dlist,
	       ftp2fn_null(efPDB,NFILE,fnm),bfac_init,
	       &atoms,x,ePBC,box,bPhi,bPsi,bChi,oenv);
  
  /* Print ramachandran maps! */
  if (bRama)
    do_rama(nf,nlist,dlist,dih,bViol,bRamOmega,oenv);
  
  if (bShift)
    do_pp2shifts(log,nf,nlist,dlist,dih);

  /* rprint S^2, transitions, and rotamer occupancies to log */ 
  traj_t_ns = 0.001 * (time[nf-1]-time[0]) ; 
  pr_dlist(log,nlist,dlist,traj_t_ns,edPrintST,bPhi,bPsi,bChi,bOmega,maxchi);
  pr_dlist(log,nlist,dlist,traj_t_ns,edPrintRO,bPhi,bPsi,bChi,bOmega,maxchi);  
  ffclose(log);
  /* transitions to xvg */
  if (bDo_rt)
    print_transitions(opt2fn("-rt",NFILE,fnm),maxchi,nlist,dlist,
		      &atoms,x,box,bPhi,bPsi,bChi,traj_t_ns,oenv); 
  
  /* chi_product trajectories (ie one "rotamer number" for each residue) */
  if (bChiProduct && bChi ) {
    snew(chi_lookup,nlist) ;
    for (i=0;i<nlist;i++)
      snew(chi_lookup[i], maxchi) ; 
    mk_chi_lookup(chi_lookup, maxchi, dih, nlist, dlist); 
    
    get_chi_product_traj(dih,nf,nactdih,nlist,
			 maxchi,dlist,time,chi_lookup,xity,
			 FALSE,bNormHisto, core_frac,bAll,
			 opt2fn("-cp",NFILE,fnm),oenv); 

    for (i=0;i<nlist;i++)
      sfree(chi_lookup[i]); 
  }

  /* Correlation comes last because it fucks up the angles */
  if (bCorr)
    do_dihcorr(opt2fn("-corr",NFILE,fnm),nf,ndih,dih,dt,nlist,dlist,time,
	       maxchi,bPhi,bPsi,bChi,bOmega,oenv);
  
  
  do_view(oenv,opt2fn("-o",NFILE,fnm),"-nxy");
  do_view(oenv,opt2fn("-jc",NFILE,fnm),"-nxy");
  if (bCorr)
    do_view(oenv,opt2fn("-corr",NFILE,fnm),"-nxy");
    
  gmx_residuetype_destroy(rt);

  thanx(stderr);
    
  return 0;
}
