/*
 *       $Id$
 *
 *       This source code is part of
 *
 *        G   R   O   M   A   C   S
 *
 * GROningen MAchine for Chemical Simulations
 *
 *            VERSION 1.6
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
static char *SRCID_g_chi_c = "$Id$";

#include <stdio.h>
#include <math.h>
#include "bondf.h"
#include "confio.h"
#include "pdbio.h"
#include "copyrite.h"
#include "fatal.h"
#include "futil.h"
#include "gstat.h"
#include "macros.h"
#include "maths.h"
#include "physics.h"
#include "rdgroup.h"
#include "smalloc.h"
#include "statutil.h"
#include "string.h"
#include "sysstuff.h"
#include "txtdump.h"
#include "typedefs.h"
#include "vec.h"
#include "strdb.h"
#include "xvgr.h"

enum { edPhi=0, edPsi, edChi1, edChi2, edChi3, edChi4, edChi5, edChi6, edMax };

#define NHISTO 360
#define NONCHI 2
#define MAXCHI edMax-NONCHI

typedef struct {
  int H,N,C,O,Cn[MAXCHI+3];
} t_dihatms; /* Cn[0]=N, Cn[1]=Ca, Cn[2]=Cb etc. */

typedef struct {
  char name[12];
  int  resnr;
  int  index;       /* Index for amino acids (histograms) */
  int  j0[edMax];   /* Index in dih array (phi angle is first...) */
  t_dihatms  atm;
  int  b[edMax];
  int  ntr[edMax];
  real S2[edMax];
} t_dlist;

void pr_dlist(FILE *fp,int nl,t_dlist dl[])
{
  int i, Xi;
  
  for(i=0; (i<nl); i++) {
    fprintf(fp,"%12s",dl[i].name);
    fprintf(fp," Phi[%d,%d,%d,%d]",dl[i].atm.H,dl[i].atm.N,
	    dl[i].atm.Cn[1],dl[i].atm.C);
    fprintf(fp," Psi[%d,%d,%d,%d]",dl[i].atm.N,dl[i].atm.Cn[1],
	    dl[i].atm.C,dl[i].atm.O);
    for(Xi=0; Xi<MAXCHI; Xi++)
      if (dl[i].atm.Cn[Xi+3] != -1)
	fprintf(fp," Xi%d[%d,%d,%d,%d]",Xi,dl[i].atm.Cn[Xi],
		dl[i].atm.Cn[Xi+1],dl[i].atm.Cn[Xi+2],
		dl[i].atm.Cn[Xi+3]);
    fprintf(fp,"\n");
  }
}

bool has_dihedral(int Dih,t_dlist *dl)
{
  bool b = FALSE;
  
#define BBB(x) (dl->atm.##x != -1)
  switch (Dih) {
  case edPhi:
    b = (BBB(H) && BBB(N) && BBB(Cn[1]) && BBB(C));
    break;
  case edPsi:
    b = (BBB(N) && BBB(Cn[1]) && BBB(C) && BBB(O));
    break;
  case edChi1:
  case edChi2:
  case edChi3:
  case edChi4:
  case edChi5:
  case edChi6:
    b = (BBB(Cn[Dih-2]) && BBB(Cn[Dih-1]) && BBB(Cn[Dih]) && BBB(Cn[Dih+1]));
    break;
  default:
    pr_dlist(stdout,1,dl);
    fatal_error(0,"Non existant dihedral %d in file %s, line %d",
		Dih,__FILE__,__LINE__);
  }
  return b;
}

t_dlist *mk_dlist(FILE *log, 
		  t_atoms *atoms, int *nlist,
		  bool bPhi, bool bPsi, bool bChi, int maxchi,
		  int r0,int naa,char **aa)
{
  int     ires,i,j,k;
  t_dihatms atm;
  int     nl=0,nc[edMax],ndih;
  bool    bDih;
  char    *thisres;
  t_dlist *dl;
  
  snew(dl,atoms->nres+1);
  
  for(i=0; (i<edMax); i++)
    nc[i]=0;
  ires = -1;
  i    =  0;
  while (i<atoms->nr) {
    ires=atoms->atom[i].resnr;
    
    /* Initiate all atom numbers to -1 */
    atm.H=atm.N=atm.C=atm.O=-1;
    for(j=0; (j<MAXCHI+3); j++)
      atm.Cn[j]=-1;
      
    /* Look for atoms in this residue */
    while ( (atoms->atom[i].resnr==ires) && (i<atoms->nr) ) {
      if ((strcmp(*(atoms->atomname[i]),"H") == 0) ||
	  (strcmp(*(atoms->atomname[i]),"H1") == 0) )
	atm.H=i;
      else if (strcmp(*(atoms->atomname[i]),"N") == 0)
	atm.N=i;
      else if (strcmp(*(atoms->atomname[i]),"C") == 0)
	atm.C=i;
      else if ((strcmp(*(atoms->atomname[i]),"O") == 0) ||
	       (strcmp(*(atoms->atomname[i]),"O1") == 0))
	atm.O=i;
      else if (strcmp(*(atoms->atomname[i]),"CA") == 0)
	atm.Cn[1]=i;
      else if (strcmp(*(atoms->atomname[i]),"CB") == 0)
	atm.Cn[2]=i;
      else if ((strcmp(*(atoms->atomname[i]),"CG") == 0) ||
	       (strcmp(*(atoms->atomname[i]),"CG1") == 0) ||
	       (strcmp(*(atoms->atomname[i]),"OG") == 0) ||
	       (strcmp(*(atoms->atomname[i]),"OG1") == 0) ||
	       (strcmp(*(atoms->atomname[i]),"SG") == 0))
	atm.Cn[3]=i;
      else if ((strcmp(*(atoms->atomname[i]),"CD") == 0) ||
	       (strcmp(*(atoms->atomname[i]),"CD1") == 0) ||
	       (strcmp(*(atoms->atomname[i]),"SD") == 0) ||
	       (strcmp(*(atoms->atomname[i]),"OD1") == 0) ||
	       (strcmp(*(atoms->atomname[i]),"ND1") == 0))
	atm.Cn[4]=i;
      else if ((strcmp(*(atoms->atomname[i]),"CE") == 0) ||
	       (strcmp(*(atoms->atomname[i]),"CE1") == 0) ||
	       (strcmp(*(atoms->atomname[i]),"OE1") == 0) ||
	       (strcmp(*(atoms->atomname[i]),"NE") == 0))
	atm.Cn[5]=i;
      else if ((strcmp(*(atoms->atomname[i]),"CZ") == 0) ||
	       (strcmp(*(atoms->atomname[i]),"NZ") == 0))
	atm.Cn[6]=i;
      else if (strcmp(*(atoms->atomname[i]),"NH1") == 0)
	atm.Cn[7]=i;
      i++;
    }
    
    /* Special case for Pro, has no H */
    if (strcmp(*(atoms->resname[ires]),"PRO") == 0) 
      atm.H=atm.Cn[4];
    thisres=*(atoms->resname[ires]);	
    
    /* Check how many dihedrals we have */
    if ((atm.N != -1) && (atm.Cn[1] != -1) && (atm.C != -1) &&
	(atm.O != -1) && (atm.H != -1)) {
      dl[nl].resnr     = ires+1;
      dl[nl].atm       = atm;
      dl[nl].atm.Cn[0] = atm.N;
      if ((atm.Cn[3] != -1) && (atm.Cn[2] != -1) && (atm.Cn[1] != -1)) {
	nc[0]++;
	if (atm.Cn[4] != -1) {
	  nc[1]++;
	  if (atm.Cn[5] != -1) {
	    nc[2]++;
	    if (atm.Cn[6] != -1) {
	      nc[3]++;
	      if (atm.Cn[7] != -1) {
		nc[4]++;
		if (atm.Cn[8] != -1) {
		  nc[5]++;
		}
	      }
	    }
	  }
	}
      }
      for(k=0; (k<naa); k++) {
	if (strcasecmp(aa[k],thisres) == 0)
	  break;
      }
      dl[nl].index=k;
      
      sprintf(dl[nl].name,"%s%d",thisres,ires+r0);
      nl++;
    }
    else if (debug)
      fprintf(debug,"Could not find N atom but could find other atoms"
	      " in residue %s%d\n",thisres,ires+r0);
  }
  fprintf(stderr,"\n");
  fprintf(log,"\n");
  fprintf(log,"There are %d residues with dihedrals\n",nl);
  j=0;
  if (bPhi) j+=nl;
  if (bPsi) j+=nl;
  if (bChi)
    for(i=0; (i<maxchi); i++)
      j+=nc[i];
  fprintf(log,"There are %d dihedrals\n",j);
  fprintf(log,"Dihedral: ");
  if (bPhi) 
    fprintf(log," Phi  ");
  if (bPsi) 
    fprintf(log," Psi  ");
  if (bChi)
    for(i=0; (i<maxchi); i++)
      fprintf(log,"Chi%d  ",i+1);
  fprintf(log,"\nNumber:   ");
  if (bPhi) 
    fprintf(log,"%4d  ",nl);
  if (bPsi) 
    fprintf(log,"%4d  ",nl);
  if (bChi)
    for(i=0; (i<maxchi); i++)
      fprintf(log,"%4d  ",nc[i]);
  fprintf(log,"\n");
  
  *nlist=nl;

  return dl;
}

int pr_trans(FILE *fp,int nl,t_dlist dl[],real dt,int Xi)
{
  int  i,nn,nz;
  
  nz=0;
  fprintf(fp,"\\begin{table}[h]\n");
  fprintf(fp,"\\caption{Number of dihedral transitions per nanosecond}\n");
  fprintf(fp,"\\begin{tabular}{|l|l|}\n");
  fprintf(fp,"\\hline\n");
  fprintf(fp,"Residue\t&$\\chi_%d$\t\\\\\n",Xi+1);
  for(i=0; (i<nl); i++) {
    nn=dl[i].ntr[Xi]/dt;
    
    if (nn == 0) {
      fprintf(fp,"%s\t&\\HL{%d}\t\\\\\n",dl[i].name,nn);
      nz++;
    }
    else if (nn > 0)
      fprintf(fp,"%s\t&\\%d\t\\\\\n",dl[i].name,nn);
  }
  fprintf(fp,"\\hline\n");
  fprintf(fp,"\\end{tabular}\n");
  fprintf(fp,"\\end{table}\n\n");

  return nz;  
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

void print_one(char *base,char *name,char *title,
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

int bin(real chi,int mult)
{
  mult=3;

  return (int) (chi*mult/360.0);
}

static void do_dihcorr(char *fn,int nf,int ndih,real **dih,real dt,
		       int nlist,t_dlist dlist[],real time[],int maxchi,
		       bool bPhi,bool bPsi,bool bChi)
{
  char name1[256],name2[256];
  int  i,j,Xi;
  
  do_autocorr(fn,"Dihedral Autocorrelation Function",
	      nf,ndih,dih,dt,eacCos,FALSE,NULL,NULL);
  /* Dump em all */
  j=0;
  for(i=0; (i<nlist); i++) {
    if (dlist[i].atm.Cn[Xi+3] != -1) {
      if (bPhi)
	print_one("corrphi",dlist[i].name,"Phi ACF for",nf/2,time,dih[j]);
      j++;
    }
  }
  for(i=0; (i<nlist); i++) {
    if (dlist[i].atm.Cn[Xi+3] != -1) {
      if (bPsi)
	print_one("corrpsi",dlist[i].name,"Psi ACF for",nf/2,time,dih[j]);
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
			real **dih,int maxchi,bool bPhi,bool bPsi,bool bChi)
{
  char name[256];
  int  i,j,Xi;
  
  /* Dump em all */
  for(i=0; (i<nlist); i++)
    if (bPhi)
      print_one("Phi",dlist[i].name,name,nf,time,dih[0]);
  for(i=0; (i<nlist); i++)
    if (bPsi)
      print_one("Psi",dlist[i].name,name,nf,time,dih[1]);
  j=2;
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

static void histogramming(FILE *log,int naa,char **aa,
			  int nf,int maxchi,real **dih,
			  int nlist,t_dlist dlist[],
			  bool bPhi,bool bPsi,bool bChi)
{
  t_karplus kkkphi[] = {
    { "J_NHa",     6.4, -1.4,  6.7, -M_PI/3, 0.0 },
    { "J_HaC'",    4.0,  1.1,  0.1, 0.0,     0.0 },
    { "J_NHCb",    4.7, -1.5, -0.2, M_PI/3,  0.0 },
    { "J_Ci-1Hai", 4.5, -1.3, -1.2, 2*M_PI/3,0.0 }
  };
  t_karplus kkkpsi[] = {
    { "J_HaN",   -0.88, -0.61,-0.27,M_PI/3,  0.0 }
  };
  t_karplus kkkchi1[] = {
    { "Jab",       9.5, -1.6, 1.8, 0, 0 }
  };
#define NJC (asize(kkkphi)+asize(kkkpsi)+asize(kkkchi1))
  
  FILE *fp;
  char buf[256];
  real S2;
  real *normhisto;
  real **Jc;
  int  ***his_aa,**his_aa1,*histmp;
  int  i,j,k,Dih;
  char name[256];
  
  snew(his_aa,edMax);
  for(Dih=0; (Dih<edMax); Dih++) {
    snew(his_aa[Dih],naa+1);
    for(i=0; (i<=naa); i++) {
      snew(his_aa[Dih][i],NHISTO);
    }
  }
  snew(histmp,NHISTO);
  
  snew(Jc,nlist);
  for(i=0; (i<nlist); i++)
    snew(Jc[i],NJC);
  
  j=0;
  for (Dih=0; (Dih<NONCHI+maxchi); Dih++) {    
    for(i=0; (i<nlist); i++) {
      if ((Dih<NONCHI) || (dlist[i].atm.Cn[Dih-NONCHI+3] != -1)) {
	make_histo(log,nf,dih[j],NHISTO,histmp,-M_PI,M_PI);

	switch (Dih) {
	case edPhi:
	  calc_distribution_props(NHISTO,histmp,-M_PI,
				  asize(kkkphi),kkkphi,&S2);
	  Jc[i][0] = kkkphi[0].Jc;
	  Jc[i][1] = kkkphi[1].Jc;
	  Jc[i][2] = kkkphi[2].Jc;
	  Jc[i][3] = kkkphi[3].Jc;
	  break;
	case edPsi:
	  calc_distribution_props(NHISTO,histmp,-M_PI,
				  asize(kkkpsi),kkkpsi,&S2);
	  Jc[i][4] = kkkpsi[0].Jc;
	  break;
	case edChi1:
	  calc_distribution_props(NHISTO,histmp,-M_PI,
				  asize(kkkchi1),kkkchi1,&S2);
	  Jc[i][5] = kkkchi1[0].Jc;
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
  for(i=0; (i<asize(kkkphi)); i++)
    fprintf(log,"%10s",kkkphi[i].name);
  for(i=0; (i<asize(kkkpsi)); i++)
    fprintf(log,"%10s",kkkpsi[i].name);
  for(i=0; (i<asize(kkkchi1)); i++)
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
	   (bChi && (Dih>=edChi1)))) {
	normalize_histo(NHISTO,his_aa[Dih][i],(360.0/NHISTO),normhisto);
	if (Dih==0) {
	  sprintf(buf,"phi%s.xvg",aa[i]);
	  sprintf(name,"\\8f\\4 Distribution for %s",aa[i]);
	} else if (Dih==1) {
	  sprintf(buf,"psi%s.xvg",aa[i]);
	  sprintf(name,"\\8y\\4 Distribution for %s",aa[i]);
	} else {
	  sprintf(buf,"chi%d%s.xvg",Dih-NONCHI+1,aa[i]);
	  sprintf(name,"\\8c\\4\\s%d\\N Distribution for %s",
		  Dih-NONCHI+1,aa[i]);
	}
	fp=xvgropen(buf,name,"Degrees","");
	fprintf(fp,"@ with g0\n");
	xvgr_world(fp,-180,0,180,0.1);
	fprintf(fp,"@ xaxis tick on\n");
	fprintf(fp,"@ xaxis tick major 90\n");
	fprintf(fp,"@ xaxis tick minor 30\n");
	fprintf(fp,"@ xaxis ticklabel prec 0\n");
	fprintf(fp,"@ yaxis tick off\n");
	fprintf(fp,"@ yaxis ticklabel off\n");
	fprintf(fp,"@ type xy\n");
	for(j=0; (j<NHISTO); j++)
	  fprintf(fp,"%5d  %10g\n",j-180,normhisto[j]);
	fprintf(fp,"&\n");
	ffclose(fp);
      }
    }
  }
  sfree(normhisto);
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

static void do_rama(int nf,int nlist,t_dlist dlist[],real **dih)
{
  FILE *fp;
  char fn[256];
  int  i,j,Xi1,Xi2,Phi,Psi;
  
  for(i=0; (i<nlist); i++) {
    if ((has_dihedral(edPhi,&(dlist[i]))) &&
	(has_dihedral(edPsi,&(dlist[i])))) {
      sprintf(fn,"ramaPhiPsi%s.xvg",dlist[i].name);
      fp = rama_file(fn,"Ramachandran Plot",
		     "\\8f\\4 (deg)","\\8y\\4 (deg)");
      Phi = dlist[i].j0[edPhi];   
      Psi = dlist[i].j0[edPsi];
      for(j=0; (j<nf); j++)
	fprintf(fp,"%10g  %10g\n",RAD2DEG*dih[Phi][j],RAD2DEG*dih[Psi][j]);
      fclose(fp);
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
  }
}

static void order_params(FILE *log,
			 char *fn,int maxchi,int nlist,t_dlist dlist[],
			 bool bPDB,char *pdbfn,real bfac_init,
			 t_topology *top,rvec x[],
			 bool bPhi,bool bPsi,bool bChi)
{
  FILE *fp;
  int  nh[edMax];
  int  i,natoms,Dih,Xi;
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
	if (dlist[i].S2[Dih]>S2Max) S2Max=dlist[i].S2[Dih];
	if (dlist[i].S2[Dih]<S2Min) S2Min=dlist[i].S2[Dih];
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
  
  if (bPDB) {
    t_pdbatom *pdba;
    int Xi;
    real x0,y0,z0;
    
    natoms=top->atoms.nr;
    pdba=atoms2pdba(&top->atoms,x);  
    pdba_trimnames(natoms,pdba);

    for(i=0; (i<natoms); i++)
      pdba[i].bfac=bfac_init;
    
    for(i=0; (i<nlist); i++) {
      pdba[dlist[i].atm.N].bfac=-dlist[i].S2[0];/* Phi */
      pdba[dlist[i].atm.H].bfac=-dlist[i].S2[0];/* Phi */
      pdba[dlist[i].atm.C].bfac=-dlist[i].S2[1];/* Psi */
      pdba[dlist[i].atm.O].bfac=-dlist[i].S2[1];/* Psi */
      for (Xi=0; (Xi<maxchi); Xi++) {           /* Chi's */
	if (dlist[i].atm.Cn[Xi+3]!=-1) {
	  pdba[dlist[i].atm.Cn[Xi+1]].bfac=-dlist[i].S2[NONCHI+Xi];
	}
      }
    }
    
    fp=ffopen(fn,"w");
    fprintf(fp,"REMARK generated by g_chi\n");
    fprintf(fp,"REMARK B-factor field contains negative of\n");
    fprintf(fp,"REMARK dihedral order parameters\n");
    print_pdbatoms(fp,natoms,pdba,NULL);
    x0=y0=z0=1000.0;
    for (i=0; (i<natoms); i++) {
      if (pdba[i].x[XX]<x0) x0=pdba[i].x[XX];
      if (pdba[i].x[YY]<y0) y0=pdba[i].x[YY];
      if (pdba[i].x[ZZ]<y0) z0=pdba[i].x[ZZ];
    }
    x0*=10.0;/* nm -> angstrom */
    y0*=10.0;/* nm -> angstrom */
    z0*=10.0;/* nm -> angstrom */
    for (i=0; (i<10); i++)
      fprintf(fp,"%6s%5d  %-4.4s%3.3s %c%4d    %8.3f%8.3f%8.3f%6.2f%6.2f\n",
	      "ATOM  ",10000+i,"CA","LEG",' ',1000,x0,y0,z0+(1.2*i),
	      0.0,-0.1*i);
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
    "g_chi computes phi, psi and chi dihedrals for all your sidechains.",
    "Output is in form of xvgr files, as well as a LaTeX table of the",
    "number of transitions per nanosecond.[PAR]",
    "Order parameters S2 for each of the dihedrals are calculated and",
    "output as xvgr file and optionally as a pdb file with the S2",
    "values as B-factor.[PAR]",
    "If option [TT]-c[tt] is given, the program will",
    "calculate dihedral autocorrelation functions. The function used",
    "is C(t) = < cos(chi(tau)) cos(chi(tau+t)) >. The use of cosines",
    "rather than angles themselves, resolves the problem of periodicity."
  };
  static char *bugs[] = {
    "Produces MANY output files (up to about 4 times the number of residues in the protein, twice that if autocorrelation functions are calculated). Typically several hundred files are output."
  };
  static int  r0=1,ndeg=1,maxchi=2,nf=10;
  static bool bAll=FALSE;
  static bool bPhi=FALSE,bPsi=FALSE,bChi=FALSE;
  static real bfac_init=-1.0;
  static bool bRama=FALSE;
  t_pargs pa[] = {
    { "-r0",  FALSE, etINT, &r0,
      "starting residue" },
    { "-phi",  FALSE, etBOOL, &bPhi,
      "Output for Phi dihedral angles" },
    { "-psi",  FALSE, etBOOL, &bPsi,
      "Output for Psi dihedral angles" },
    { "-chi",  FALSE, etBOOL, &bChi,
      "Output for Chi dihedral angles" },
    { "-rama", FALSE, etBOOL, &bRama,
      "Generate Phi/Psi and Chi1/Chi2 ramachandran plots" },
    { "-all",  FALSE, etBOOL, &bAll,
      "Output separate files for every dihedral." },
    { "-nframes", FALSE, etINT, &nf,
      "Number of frames in your trajectory" },
    { "-run", FALSE, etINT, &ndeg,
      "perform running average over ndeg degrees for histograms" },
    { "-maxchi", FALSE, etINT, &maxchi,
      "calculate first ndih Chi dihedrals (max 6)" },
    { "-bfact", FALSE, etREAL, &bfac_init,
      "bfactor value for pdb file for atoms with no calculated dihedral order parameter"}
  };

  FILE       *fp,*log;
  char       title[256];
  int        status,natoms,i,j,k,l;
  t_topology top;
  rvec       *x,*xref,*xav;
  real       t,t0,t1,lambda;
  matrix     box;
  t_dlist    *dlist;
  int        nlist,teller,naa,idum;
  char       **aa;
  int        th,th1,th2;
  bool       bCorr;
  real       dt,rh1,rh2,rj,invth,tdc,tds;

  atom_id    isize,*index;
  int        ndih;
  real       **dih,*trans_frac,*aver_angle,*time;
  
  t_filenm  fnm[] = {
    { efTPB, NULL,  NULL,     ffREAD  },
    { efTRX, "-f",  NULL,     ffREAD  },
    { efXVG, "-o",  "order",  ffWRITE },
    { efPDB, "-p",  "order",  ffOPTWR },
    { efXVG, "-jc", "Jcoupling", ffWRITE },
    { efXVG, "-c",  "dihcorr",ffOPTWR },
    { efTEX, "-t",  "trans",  ffWRITE },
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
  
  log=ffopen(ftp2fn(efLOG,NFILE,fnm),"w");
  
  bCorr=(opt2bSet("-c",NFILE,fnm));
  if (bCorr) 
    fprintf(stderr,"Will calculate autocorrelation\n");

  if (maxchi > MAXCHI) {
    fprintf(stderr, 
	    "Will only calculate first %d Chi dihedrals in stead of %d.\n",
	    MAXCHI, maxchi);
    maxchi=MAXCHI;
  }

  /* Find the chi angles using topology and a list of amino acids */
  {
    t_statheader sh;
    int dint;/* dummy */
    real dreal;/* dummy */
    
    fp=ftp2FILE(efTPB,NFILE,fnm,"r");
    rd_header(fp, &sh);
    natoms = sh.natoms;
    snew(x,natoms);
    rd_hstatus(fp,&sh,&dint,&dreal,&dreal,NULL,NULL,NULL,NULL,
	       &natoms,x,NULL,NULL,&dint,NULL,&top);
    fprintf(log,"topology name = %s\n",*top.name);
  }
  
  naa=get_strings("aminoacids.dat",&aa);
  dlist=mk_dlist(log,&(top.atoms),&nlist,bPhi,bPsi,bChi,maxchi,r0,naa,aa);
  fprintf(stderr,"%d residues with dihedrals found\n", nlist);
  
  if (nlist == 0) 
    fatal_error(0,"No dihedrals in your topology!\n");
  
  /* Make a linear index for reading all */
  index=make_chi_ind(nlist,dlist,&ndih);
  isize=4*ndih;
  fprintf(stderr,"%d dihedrals found\n", ndih);

  if (nf <= 0) 
    fatal_error(0,"No frames (%d) in trajectory ? DIY!\n",nf);
  snew(time,nf);
  snew(trans_frac,nf);
  snew(aver_angle,nf);
  snew(dih,ndih);
  for(i=0; (i<ndih); i++)
    snew(dih[i],nf);
    
  /* COMPUTE ALL DIHEDRALS! */
  read_ang_dih(opt2fn("-f",NFILE,fnm),ftp2fn(efTPB,NFILE,fnm),
	       FALSE,TRUE,FALSE,1,&idum,
	       &nf,time,isize,index,trans_frac,aver_angle,dih);
  
  if (nf < 2)
    fatal_error(0,"No frames in trajectory %s",opt2fn("-f",NFILE,fnm));
    
  dt=(time[nf-1]-time[0])/(nf-1);
  
  /* put angles in -M_PI to M_PI ! */
  for(i=0; (i<ndih); i++) {
    for(j=0; (j<nf); j++) {
      while (dih[i][j] < -M_PI)
	dih[i][j] += 2*M_PI;
      while (dih[i][j] > M_PI)
	dih[i][j] -= 2*M_PI;
    }
  }
  
  if (bAll)
    dump_em_all(nlist,dlist,nf,time,dih,maxchi,bPhi,bPsi,bChi);
  
  /* Histogramming & J coupling constants */
  histogramming(log,naa,aa,nf,maxchi,dih,nlist,dlist,bPhi,bPsi,bChi);

  /* Order parameters */  
  order_params(log,opt2fn("-o",NFILE,fnm),maxchi,nlist,dlist,
	       opt2bSet("-p",NFILE,fnm),ftp2fn(efPDB,NFILE,fnm),bfac_init,
	       &top,x,bPhi,bPsi,bChi);
  
  if (bCorr)
    do_dihcorr(opt2fn("-c",NFILE,fnm),nf,ndih,dih,dt,nlist,dlist,time,
	       maxchi,bPhi,bPsi,bChi);

  /* Print ramachandran maps! */
  if (bRama)
    do_rama(nf,nlist,dlist,dih);
  
  
  pr_dlist(log,nlist,dlist);  
  ffclose(log);
  
  xvgr_file(opt2fn("-o",NFILE,fnm),"-nxy");
  xvgr_file(opt2fn("-jc",NFILE,fnm),"-nxy");
  if (bCorr)
    xvgr_file(opt2fn("-c",NFILE,fnm),"-nxy");
    
  thanx(stdout);
    
  return 0;
}
