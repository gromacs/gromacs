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
 * Gnomes, ROck Monsters And Chili Sauce
 */
static char *SRCID_addconf_c = "$Id$";

#include "vec.h"
#include "macros.h"
#include "smalloc.h"
#include "addconf.h"
#include "gstat.h"

void add_conf(t_atoms *atoms_1,rvec *x_1,rvec *v_1,real *r_1,
	      int NTB,matrix box_1,
	      t_atoms *atoms_2,rvec *x_2,rvec *v_2,real *r_2,
	      bool bVerbose,int maxmol)
{
  int  i,j,m,resnr,oldres,nresadd;
  int  *add;
  real d2,vdw2;

  if (atoms_2->nr <= 0) {
    fatal_error(0,"Nothing to add");
  }
  
  if (bVerbose)
    fprintf(stderr,"Calculating Overlap...\n");
    
  /* The add array is filled with:
   * 0 for atoms that should NOT be added,
   * 1 for atoms that should be added
   * 2 for atoms that should be added and that are the start of a residue
   */
  snew(add,atoms_2->nr);
  
  /* First set the beginning of residues */
  add[0]  = 2;
  nresadd = 1;
  for (i=1;i<atoms_2->nr;i++) {
    add[i]=1;
    if (atoms_2->atom[i].resnr != atoms_2->atom[i-1].resnr) {
      add[i]  = 2;
      nresadd ++;
    }
  }
  fprintf(stderr,"There are %d residues before testing overlap\n",nresadd);
  if (atoms_1->nr == 0) {
    /* Fill an empty box with solvent */
    for (j=0; (j<atoms_2->nr); j++) {
      while ((add[j]==0) && (j<atoms_2->nr))
	j++;
      if (in_box(NTB,box_1,(x_2)[j]) != 0) {
	if (add[j] == 2)
	  nresadd--;
	add[j] = 0;
      }
      if (add[j]==0) {	
        resnr=atoms_2->atom[j].resnr;
	while((resnr==atoms_2->atom[j-1].resnr) && (j>0))
	  j--;
	if (j<0)
	  warning("j < 0 in addconf");
	while((resnr==atoms_2->atom[j].resnr) && (j<atoms_2->nr)) {
	  if (add[j] == 2)
	    nresadd--;
          add[j] = 0;
	  j++;
        }
      }
    }
  } 
  else {
    /* if an atom of configuration_2 is too close to an atom of configuration_1
     * add[j]=0
     * if an atom of configuration_2 is not inside the box add[j]=0
     */
    for (i=0; (i<atoms_1->nr); i++) {
      if (bVerbose && ((i % 10)==0))
	fprintf(stderr,"\r%d atoms out of %d",i,atoms_1->nr);
      for (j=0; (j<atoms_2->nr); j++) {
	d2=dist2((x_1)[i],(x_2)[j],box_1);
	vdw2=sqr((r_1)[i] + (r_2)[j]);
	if ((d2 < vdw2) || (in_box(NTB,box_1,(x_2)[j]) != 0)) {
	  if (add[j] == 2)
	    nresadd--;
	  add[j] = 0;
	}
	if (add[j]==0) {	
	  resnr=atoms_2->atom[j].resnr;
	  while((resnr==atoms_2->atom[j-1].resnr) && (j>0))
	    j--;
	  if (j<0)
	    warning("j < 0 in addconf");
	  while((resnr==atoms_2->atom[j].resnr) && (j<atoms_2->nr)) {
	    if (add[j] == 2)
	      nresadd--;
	    add[j]=0;
	    j++;
	  }
	}
      }
    }
    if (bVerbose)
      fprintf(stderr,"\n");
  }
  /* add[j]=0 of all the atoms of a molecule containing one or more atoms 
   * add[i]=0
   */
  if (maxmol == 0)
    maxmol = nresadd;
  else
    maxmol = min(maxmol,nresadd);
  fprintf(stderr,"There are %d molecules to add\n",maxmol);
  
  /*add the selected atoms_2 to atoms_1*/  
  oldres = atoms_1->nres;
  for (i=0; (i<atoms_2->nr); i++)
    if (add[i] != 0) {
      if (add[i] == 2) {
	if (maxmol == 0)
	  break;
	atoms_1->nres++;
	maxmol--;
      }
      atoms_1->nr++;
      atoms_1->atomname[atoms_1->nr-1] = atoms_2->atomname[i];
      copy_rvec((x_2)[i],(x_1)[atoms_1->nr-1]);
      copy_rvec((v_2)[i],(v_1)[atoms_1->nr-1]);
      (r_1)[atoms_1->nr-1]   =(r_2)[i];
      atoms_1->atom[atoms_1->nr-1].resnr=atoms_1->nres-1;
      atoms_1->resname[atoms_1->nres-1]=
	atoms_2->resname[atoms_2->atom[i].resnr];
    }
  sfree(add);
}/*add_conf()*/

