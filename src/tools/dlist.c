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
static char *SRCID_dlist_c = "$Id$";

#include <stdlib.h>
#include "string2.h"
#include "pp2shift.h"
#include "smalloc.h"
	
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
      else if ((strcmp(*(atoms->atomname[i]),"CG") == 0)  ||
	       (strcmp(*(atoms->atomname[i]),"CG1") == 0) ||
	       (strcmp(*(atoms->atomname[i]),"OG") == 0)  ||
	       (strcmp(*(atoms->atomname[i]),"OG1") == 0) ||
	       (strcmp(*(atoms->atomname[i]),"SG") == 0))
	atm.Cn[3]=i;
      else if ((strcmp(*(atoms->atomname[i]),"CD") == 0)  ||
	       (strcmp(*(atoms->atomname[i]),"CD1") == 0) ||
	       (strcmp(*(atoms->atomname[i]),"SD") == 0)  ||
	       (strcmp(*(atoms->atomname[i]),"OD1") == 0) ||
	       (strcmp(*(atoms->atomname[i]),"ND1") == 0) ||
	       (strcmp(*(atoms->atomname[i]),"HG")  == 0) ||
	       (strcmp(*(atoms->atomname[i]),"HG1")  == 0))
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

bool has_dihedral(int Dih,t_dlist *dl)
{
  bool b = FALSE;
  int  ddd;
  
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
    ddd = Dih - edChi1;
    b   = (BBB(Cn[ddd]) && BBB(Cn[ddd+1]) && BBB(Cn[ddd+2]) && BBB(Cn[ddd+3]));
    break;
  default:
    pr_dlist(stdout,1,dl,1);
    fatal_error(0,"Non existant dihedral %d in file %s, line %d",
		Dih,__FILE__,__LINE__);
  }
  return b;
}

static void pr_props(FILE *fp,t_dlist *dl,int nDih,real dt)
{
  fprintf(fp,"  %6.2f  %6.2f\n",(dt ==0 ) ?  0 : dl->ntr[nDih]/dt,
	  dl->S2[nDih]);
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
	fprintf(fp,"   Chi%d[%4d,%4d,%4d,%4d]",Xi+1,dl[i].atm.Cn[Xi],
		dl[i].atm.Cn[Xi+1],dl[i].atm.Cn[Xi+2],
		dl[i].atm.Cn[Xi+3]);
	pr_props(fp,&dl[i],Xi+2,dt);
      }
    fprintf(fp,"\n");
  }
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

