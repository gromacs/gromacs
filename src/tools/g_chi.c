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
 * GROningen Mixture of Alchemy and Childrens' Stories
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
#include "tpxio.h"
#include "string.h"
#include "sysstuff.h"
#include "txtdump.h"
#include "typedefs.h"
#include "vec.h"
#include "strdb.h"
#include "xvgr.h"
#include "pp2shift.h"

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

static void pr_props(FILE *fp,t_dlist *dl,int nDih,real dt)
{
  fprintf(fp,"  %6.2f  %6.2f\n",dl->ntr[nDih]/dt,dl->S2[nDih]);
}

void pr_dlist(FILE *fp,int nl,t_dlist dl[],real dt)
{
  int i, Xi;
  
  for(i=0; (i<nl); i++) {
    fprintf(fp,"Residue %s\n",dl[i].name);
    fprintf(fp," Angle [  AI,  AJ,  AK,  AL]  #tr/ns  S^2D  \n"
	       "--------------------------------------------\n");
    fprintf(fp,"   Phi [%4d,%4d,%4d,%4d]",
	    (dl[i].atm.H == -1) ? dl[i].atm.minC : dl[i].atm.H,
	    dl[i].atm.N,dl[i].atm.Cn[1],dl[i].atm.C);
    pr_props(fp,&dl[i],edPhi,dt);
    fprintf(fp,"   Psi [%4d,%4d,%4d,%4d]",dl[i].atm.N,dl[i].atm.Cn[1],
	    dl[i].atm.C,dl[i].atm.O);
    pr_props(fp,&dl[i],edPsi,dt);
    fprintf(fp," Omega [%4d,%4d,%4d,%4d]",dl[i].atm.minO,dl[i].atm.minC,
	    dl[i].atm.N,dl[i].atm.Cn[1]);
    pr_props(fp,&dl[i],edOmega,dt);
    for(Xi=0; Xi<MAXCHI; Xi++)
      if (dl[i].atm.Cn[Xi+3] != -1) {
	fprintf(fp,"   Chi%d[%4d,%4d,%4d,%4d]",Xi,dl[i].atm.Cn[Xi],
		dl[i].atm.Cn[Xi+1],dl[i].atm.Cn[Xi+2],
		dl[i].atm.Cn[Xi+3]);
	pr_props(fp,&dl[i],Xi+2,dt);
      }
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
  case edOmega:
    b = (BBB(minO) && BBB(minC) && BBB(N) && BBB(Cn[1]));
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
    pr_dlist(stdout,1,dl,1);
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
  t_dihatms atm,prev;
  int     nl=0,nc[edMax],ndih;
  bool    bDih;
  char    *thisres;
  t_dlist *dl;
  
  snew(dl,atoms->nres+1);
  prev.C = -1;
  for(i=0; (i<edMax); i++)
    nc[i]=0;
  ires = -1;
  i    =  0;
  while (i<atoms->nr) {
    ires=atoms->atom[i].resnr;
    
    /* Initiate all atom numbers to -1 */
    atm.minC=atm.H=atm.N=atm.C=atm.O=-1;
    for(j=0; (j<MAXCHI+3); j++)
      atm.Cn[j]=-1;
      
    /* Look for atoms in this residue */
    while ((i<atoms->nr) && (atoms->atom[i].resnr == ires)) {
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
    /* Carbon from previous residue */
    if (prev.C != -1)
      atm.minC = prev.C;
    if (prev.O != -1)
      atm.minO = prev.O;
    prev = atm;
      
    thisres=*(atoms->resname[ires]);	
    
    /* Check how many dihedrals we have */
    if ((atm.N != -1) && (atm.Cn[1] != -1) && (atm.C != -1) &&
      (atm.O != -1) && ((atm.H != -1) || (atm.minC != -1))) {
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
      if ((atm.minC != -1) && (atm.minO != -1))
	nc[6]++;
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
		       bool bPhi,bool bPsi,bool bChi)
{
  char name1[256],name2[256];
  int  i,j,Xi;
  
  do_autocorr(fn,"Dihedral Autocorrelation Function",
	      nf,ndih,dih,dt,eacCos,FALSE,NULL,NULL);
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
  for(i=0; (i<nlist); i++)
    if (bPhi)
      print_one("phi",dlist[i].name,name,nf,time,dih[edPhi]);
  for(i=0; (i<nlist); i++)
    if (bPsi)
      print_one("psi",dlist[i].name,name,nf,time,dih[edPsi]);
  for(i=0; (i<nlist); i++)
    if (bOmega && has_dihedral(edOmega,&(dlist[i])))
      print_one("omega",dlist[i].name,name,nf,time,dih[edOmega]);
    
  j = edChi1;
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
			  bool bPhi,bool bPsi,bool bOmega,bool bChi)
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
    { "Jab",       9.5, -1.6, 1.8, 0, 0 }
  };
#define NKKKPHI asize(kkkphi)
#define NKKKPSI asize(kkkpsi)
#define NKKKCHI asize(kkkchi1)
#define NJC (NKKKPHI+NKKKPSI+NKKKCHI)
  
  FILE *fp;
  real S2;
  real *normhisto;
  real **Jc;
  int  ***his_aa,**his_aa1,*histmp;
  int  i,j,k,m,Dih;
  char hisfile[256],title[256];
  
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
      if (((Dih  < edOmega) ) ||
	  ((Dih == edOmega) && (has_dihedral(edOmega,&(dlist[i])))) ||
	  ((Dih  > edOmega) && (dlist[i].atm.Cn[Dih-NONCHI+3] != -1))) {
      	make_histo(log,nf,dih[j],NHISTO,histmp,-M_PI,M_PI);
	
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
	normalize_histo(NHISTO,his_aa[Dih][i],(360.0/NHISTO),normhisto);
	
	switch (Dih) {
	case edPhi:
	  sprintf(hisfile,"histo-phi%s.xvg",aa[i]);
	  sprintf(title,"\\8f\\4 Distribution for %s",aa[i]);
	  break;
	case edPsi:
	  sprintf(hisfile,"histo-psi%s.xvg",aa[i]);
	  sprintf(title,"\\8y\\4 Distribution for %s",aa[i]);
	  break;
	case edOmega:
	  sprintf(hisfile,"histo-omega%s.xvg",aa[i]);
	  sprintf(title,"\\8w\\4 Distribution for %s",aa[i]);
	  break;
	default:
	  sprintf(hisfile,"histo-chi%d%s.xvg",Dih-NONCHI+1,aa[i]);
	  sprintf(title,"\\8c\\4\\s%d\\N Distribution for %s",
		  Dih-NONCHI+1,aa[i]);
	}
	fp=xvgropen(hisfile,title,"Degrees","");
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

static void do_rama(int nf,int nlist,t_dlist dlist[],real **dih,
		    bool bViol)
{
  FILE *fp,*gp;
  char fn[256];
  int  i,j,Xi1,Xi2,Phi,Psi;
  
  for(i=0; (i<nlist); i++) {
    if ((has_dihedral(edPhi,&(dlist[i]))) &&
	(has_dihedral(edPsi,&(dlist[i])))) {
      sprintf(fn,"ramaPhiPsi%s.xvg",dlist[i].name);
      fp = rama_file(fn,"Ramachandran Plot",
		     "\\8f\\4 (deg)","\\8y\\4 (deg)");
      if (bViol) {
	sprintf(fn,"violPhiPsi%s.xvg",dlist[i].name);
	gp = ffopen(fn,"w");
      }
      Phi = dlist[i].j0[edPhi];   
      Psi = dlist[i].j0[edPsi];
      for(j=0; (j<nf); j++) {
	fprintf(fp,"%10g  %10g\n",RAD2DEG*dih[Phi][j],RAD2DEG*dih[Psi][j]);
	if (bViol)
	  fprintf(gp,"%d\n",!bAllowed(dih[Phi][j],RAD2DEG*dih[Psi][j]));
      }
      if (bViol)
	fclose(gp);
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
  
  if (bPDB) {
    int Xi;
    real x0,y0,z0;
    
    for(i=0; (i<top->atoms.nr); i++)
      top->atoms.pdbinfo[i].bfac=bfac_init;
    
    for(i=0; (i<nlist); i++) {
      top->atoms.pdbinfo[dlist[i].atm.N].bfac=-dlist[i].S2[0];/* Phi */
      top->atoms.pdbinfo[dlist[i].atm.H].bfac=-dlist[i].S2[0];/* Phi */
      top->atoms.pdbinfo[dlist[i].atm.C].bfac=-dlist[i].S2[1];/* Psi */
      top->atoms.pdbinfo[dlist[i].atm.O].bfac=-dlist[i].S2[1];/* Psi */
      for (Xi=0; (Xi<maxchi); Xi++) {           /* Chi's */
	if (dlist[i].atm.Cn[Xi+3]!=-1) {
	  top->atoms.pdbinfo[dlist[i].atm.Cn[Xi+1]].bfac=-dlist[i].S2[NONCHI+Xi];
	}
      }
    }
    
    fp=ffopen(fn,"w");
    fprintf(fp,"REMARK generated by g_chi\n");
    fprintf(fp,"REMARK B-factor field contains negative of\n");
    fprintf(fp,"REMARK dihedral order parameters\n");
    write_pdbfile(fp,NULL,&(top->atoms),x,NULL,0,TRUE);
    x0=y0=z0=1000.0;
    for (i=0; (i<top->atoms.nr); i++) {
      if (x[i][XX]<x0) x0=x[i][XX];
      if (x[i][YY]<y0) y0=x[i][YY];
      if (x[i][ZZ]<y0) z0=x[i][ZZ];
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
    "(Van der Spoel & Berendsen (1997), [BB]Biophys. J. 72[bb], 2032-2041)."
  };
  
  static char *bugs[] = {
    "Produces MANY output files (up to about 4 times the number of residues in the protein, twice that if autocorrelation functions are calculated). Typically several hundred files are output."
  };
  static int  r0=1,ndeg=1,maxchi=2;
  static bool bAll=FALSE;
  static bool bPhi=FALSE,bPsi=FALSE,bOmega=FALSE;
  static real bfac_init=-1.0;
  static char *maxchistr[] = { "0", "1", "2", "3",  "4", "5", "6", NULL };
  static bool bRama=FALSE,bShift=FALSE,bViol=FALSE;
  t_pargs pa[] = {
    { "-r0",  FALSE, etINT, &r0,
      "starting residue" },
    { "-phi",  FALSE, etBOOL, &bPhi,
      "Output for Phi dihedral angles" },
    { "-psi",  FALSE, etBOOL, &bPsi,
      "Output for Psi dihedral angles" },
    { "-omega",FALSE, etBOOL, &bOmega,  
      "Output for Omega dihedrals (peptide bonds)" },
    { "-rama", FALSE, etBOOL, &bRama,
      "Generate Phi/Psi and Chi1/Chi2 ramachandran plots" },
    { "-viol", FALSE, etBOOL, &bViol,
      "Write a file that gives 0 or 1 for violated Ramachandran angles" },
    { "-all",  FALSE, etBOOL, &bAll,
      "Output separate files for every dihedral." },
    { "-shift", FALSE, etBOOL, &bShift,
	"Compute chemical shifts from Phi/Psi angles" },
    { "-run", FALSE, etINT, &ndeg,
      "perform running average over ndeg degrees for histograms" },
    { "-maxchi", FALSE, etENUM, maxchistr,
      "calculate first ndih Chi dihedrals" },
    { "-bfact", FALSE, etREAL, &bfac_init,
      "bfactor value for pdb file for atoms with no calculated dihedral order parameter"}
  };

  FILE       *fp,*log;
  char       title[256];
  int        status,natoms,i,j,k,l;
  bool       bChi;
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
  int        ndih,nf;
  real       **dih,*trans_frac,*aver_angle,*time;
  
  t_filenm  fnm[] = {
    { efTPX, NULL,  NULL,     ffREAD  },
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

  /* Handle result from enumerated type */
  sscanf(maxchistr[0],"%d",&maxchi);
  bChi = (maxchi > 0);
  
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
    t_tpxheader sh;
    int dint;/* dummy */
    real dreal;/* dummy */
    char *fn;
    
    fn=ftp2fn(efTPX,NFILE,fnm);
    read_tpxheader(fn, &sh);
    natoms = sh.natoms;
    snew(x,natoms);
    read_tpx(fn,&dint,&dreal,&dreal,NULL,NULL,
	     &natoms,x,NULL,NULL,&top);
    snew(top.atoms.pdbinfo,top.atoms.nr);
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

  snew(dih,ndih);
    
  /* COMPUTE ALL DIHEDRALS! */
  read_ang_dih(opt2fn("-f",NFILE,fnm),ftp2fn(efTPX,NFILE,fnm),
	       FALSE,TRUE,FALSE,1,&idum,
	       &nf,&time,isize,index,&trans_frac,&aver_angle,dih);
  
  if (nf < 2)
    fatal_error(0,"No frames in trajectory %s",opt2fn("-f",NFILE,fnm));
    
  dt=(time[nf-1]-time[0])/(nf-1);
  
  /* put angles in -M_PI to M_PI ! and correct phase factor for phi and psi */
  reset_em_all(nlist,dlist,nf,dih,maxchi,bPhi,bPsi,bChi);
  
  if (bAll)
    dump_em_all(nlist,dlist,nf,time,dih,maxchi,bPhi,bPsi,bChi,bOmega);
  
  /* Histogramming & J coupling constants */
  histogramming(log,naa,aa,nf,maxchi,dih,nlist,dlist,bPhi,bPsi,bOmega,bChi);

  /* Order parameters */  
  order_params(log,opt2fn("-o",NFILE,fnm),maxchi,nlist,dlist,
	       opt2bSet("-p",NFILE,fnm),ftp2fn(efPDB,NFILE,fnm),bfac_init,
	       &top,x,bPhi,bPsi,bChi);
  
  /* Print ramachandran maps! */
  if (bRama)
    do_rama(nf,nlist,dlist,dih,bViol);
  
  if (bShift)
    do_pp2shifts(log,nf,nlist,dlist,dih);
  
  pr_dlist(log,nlist,dlist,time[nf-1]-time[0]);
  ffclose(log);
  
  /* Correlation comes last because it fucks up the angles */
  if (bCorr)
    do_dihcorr(opt2fn("-c",NFILE,fnm),nf,ndih,dih,dt,nlist,dlist,time,
	       maxchi,bPhi,bPsi,bChi);

  
  xvgr_file(opt2fn("-o",NFILE,fnm),"-nxy");
  xvgr_file(opt2fn("-jc",NFILE,fnm),"-nxy");
  if (bCorr)
    xvgr_file(opt2fn("-c",NFILE,fnm),"-nxy");
    
  thanx(stdout);
    
  return 0;
}
