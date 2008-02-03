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
 * Gallium Rubidium Oxygen Manganese Argon Carbon Silicon
 */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <math.h>
#include "smalloc.h"
#include "sysstuff.h"
#include "macros.h"
#include "string2.h"
#include "topdirs.h"
#include "toputil.h"
#include "topdirs.h"
#include "toputil.h"
#include "symtab.h"
#include "gmx_fatal.h"
#include "txtdump.h"
#include "gpp_atomtype.h"

typedef struct {
  int           nr;		/* The number of atomtypes		*/
  t_atom	*atom;		/* Array of atoms			*/
  char          ***atomname;	/* Names of the atomtypes		*/
  t_param	*nb;		/* Nonbonded force default params	*/
  int           *bondatomtype;  /* The bond_atomtype for each atomtype  */
  real          *radius;        /* Radius for GBSA stuff                */
  real          *vol;           /* Effective volume for GBSA            */
  real          *surftens;      /* Surface tension with water, for GBSA */
  int           *atomnumber;    /* Atomic number, used for QM/MM        */
} gpp_atomtype;

int get_atomtype_type(char *str,t_atomtype at)
{
  gpp_atomtype *ga = (gpp_atomtype *) at;
  
  int i;

  for (i=0; (i<ga->nr); i++)
    if (strcasecmp(str,*(ga->atomname[i])) == 0)
      return i;
  
  return NOTSET;
}

int get_atomtype_ntypes(t_atomtype at)
{
  gpp_atomtype *ga = (gpp_atomtype *) at;
  
  return ga->nr;
}

char *get_atomtype_name(int nt, t_atomtype at)
{
  gpp_atomtype *ga = (gpp_atomtype *) at;
  
  if ((nt < 0) || (nt >= ga->nr))
    return NULL;
  
  return *(ga->atomname[nt]);
}

real get_atomtype_massA(int nt,t_atomtype at)
{
  gpp_atomtype *ga = (gpp_atomtype *) at;
  
  if ((nt < 0) || (nt >= ga->nr))
    return NOTSET;
  
  return ga->atom[nt].m;
}

real get_atomtype_massB(int nt,t_atomtype at)
{
  gpp_atomtype *ga = (gpp_atomtype *) at;
  
  if ((nt < 0) || (nt >= ga->nr))
    return NOTSET;
  
  return ga->atom[nt].mB;
}

real get_atomtype_qA(int nt,t_atomtype at)
{
  gpp_atomtype *ga = (gpp_atomtype *) at;
  
  if ((nt < 0) || (nt >= ga->nr))
    return NOTSET;
  
  return ga->atom[nt].q;
}

real get_atomtype_qB(int nt,t_atomtype at)
{
  gpp_atomtype *ga = (gpp_atomtype *) at;
  
  if ((nt < 0) || (nt >= ga->nr))
    return NOTSET;
  
  return ga->atom[nt].qB;
}

int get_atomtype_ptype(int nt,t_atomtype at)
{
  gpp_atomtype *ga = (gpp_atomtype *) at;
  
  if ((nt < 0) || (nt >= ga->nr))
    return NOTSET;
  
  return ga->atom[nt].ptype;
}

int get_atomtype_batype(int nt,t_atomtype at)
{
  gpp_atomtype *ga = (gpp_atomtype *) at;
  
  if ((nt < 0) || (nt >= ga->nr))
    return NOTSET;
  
  return ga->bondatomtype[nt];
}

int get_atomtype_atomnumber(int nt,t_atomtype at)
{
  gpp_atomtype *ga = (gpp_atomtype *) at;
  
  if ((nt < 0) || (nt >= ga->nr))
    return NOTSET;
  
  return ga->atomnumber[nt];
}

real get_atomtype_radius(int nt,t_atomtype at)
{
  gpp_atomtype *ga = (gpp_atomtype *) at;
  
  if ((nt < 0) || (nt >= ga->nr))
    return NOTSET;
  
  return ga->radius[nt];
}

real get_atomtype_vol(int nt,t_atomtype at)
{
  gpp_atomtype *ga = (gpp_atomtype *) at;
  
  if ((nt < 0) || (nt >= ga->nr))
    return NOTSET;
  
  return ga->vol[nt];
}

real get_atomtype_surftens(int nt,t_atomtype at)
{
  gpp_atomtype *ga = (gpp_atomtype *) at;
  
  if ((nt < 0) || (nt >= ga->nr))
    return NOTSET;
  
  return ga->surftens[nt];
}

real get_atomtype_nbparam(int nt,int param,t_atomtype at)
{
  gpp_atomtype *ga = (gpp_atomtype *) at;
  
  if ((nt < 0) || (nt >= ga->nr))
    return NOTSET;
  if ((param < 0) || (param >= MAXFORCEPARAM))
    return NOTSET;
  return ga->nb[nt].c[param];
}

t_atomtype init_atomtype(void)
{
  gpp_atomtype *ga;
  
  snew(ga,1);
  
  return (t_atomtype ) ga;
}

int set_atomtype(int nt,t_atomtype at,t_symtab *tab,
		 t_atom *a,char *name,t_param *nb,
		 int bondatomtype,
		 real radius,real vol,real surftens,int atomnumber)
{
  gpp_atomtype *ga = (gpp_atomtype *) at;
  
  if ((nt < 0) || (nt >= ga->nr))
    return NOTSET;
    
  ga->atom[nt] = *a;
  ga->atomname[nt] = put_symtab(tab,name);
  ga->nb[nt] = *nb;
  ga->bondatomtype[nt] = bondatomtype;
  ga->radius[nt] = radius;
  ga->vol[nt] = vol;
  ga->surftens[nt] = surftens;
  ga->atomnumber[nt] = atomnumber;
  
  return nt;
}

int add_atomtype(t_atomtype at,t_symtab *tab,
		 t_atom *a,char *name,t_param *nb,
		 int bondatomtype,
		 real radius,real vol,real surftens,real atomnumber)
{
  gpp_atomtype *ga = (gpp_atomtype *) at;
  
  ga->nr++;
  srenew(ga->atom,ga->nr);
  srenew(ga->atomname,ga->nr);
  srenew(ga->nb,ga->nr);
  srenew(ga->bondatomtype,ga->nr);
  srenew(ga->radius,ga->nr);
  srenew(ga->vol,ga->nr);
  srenew(ga->surftens,ga->nr);
  srenew(ga->atomnumber,ga->nr);
  
  return set_atomtype(ga->nr-1,at,tab,a,name,nb,bondatomtype,radius,
		      vol,surftens,atomnumber);
}

void print_at (FILE * out, t_atomtype at)
{
  gpp_atomtype *ga = (gpp_atomtype *) at;
  
  int i;
  t_atom     *atom = ga->atom;
  t_param    *nb   = ga->nb;
  
  fprintf (out,"[ %s ]\n",dir2str(d_atomtypes));
  fprintf (out,"; %6s  %8s  %8s  %8s  %12s  %12s\n",
	   "type","mass","charge","particle","c6","c12");
  for (i=0; (i<ga->nr); i++) 
    fprintf(out,"%8s  %8.3f  %8.3f  %8s  %12e  %12e\n",
	    *(ga->atomname[i]),atom[i].m,atom[i].q,"A",
	    nb[i].C0,nb[i].C1);
  
  fprintf (out,"\n");
}

void done_atomtype(t_atomtype *at)
{
  gpp_atomtype *ga = (gpp_atomtype *) *at;

  sfree(ga->atom);
  sfree(ga->atomname);
  sfree(ga->nb);
  sfree(ga->bondatomtype);
  sfree(ga->radius);
  sfree(ga->vol);
  sfree(ga->surftens);
  sfree(ga->atomnumber);
  ga->nr = 0;
  sfree(ga);
  *at = NULL;
}

static int search_atomtypes(t_atomtype at,int *n,int typelist[],
			    int thistype,
			    t_param param[],int ftype)
{
  int i,nn,nrfp,j,k,ntype,tli;
  bool bFound=FALSE;
  
  nn    = *n;
  nrfp  = NRFP(ftype);
  ntype = get_atomtype_ntypes(at);

  for(i=0; (i<nn); i++) 
  {
    if (typelist[i] == thistype)
    {
      /* This type number has already been added */
      break;
    }

    /* Otherwise, check if the parameters are identical to any previously added type */
    
    bFound=TRUE;
    for(j=0;j<ntype && bFound;j++) 
    {
      /* Check nonbonded parameters */
      for(k=0;k<nrfp && bFound;k++) 
      {
        bFound=(param[ntype*typelist[i]+j].c[k]==param[ntype*thistype+j].c[k]);
      }

      /* Check radius, volume, surftens */
      tli = typelist[i];
      bFound = bFound && 
	(get_atomtype_radius(tli,at) == get_atomtype_radius(thistype,at)) &&
	(get_atomtype_vol(tli,at) == get_atomtype_vol(thistype,at)) &&
	(get_atomtype_surftens(tli,at) == get_atomtype_surftens(thistype,at)) &&
	(get_atomtype_atomnumber(tli,at) == get_atomtype_atomnumber(thistype,at));
    }
    if (bFound)
    {
      break;
    }
  }
  
  if (i == nn) {
    if (debug)
      fprintf(debug,"Renumbering atomtype %d to %d\n",thistype,nn);
    if (nn == ntype)
      gmx_fatal(FARGS,"Atomtype horror n = %d, %s, %d",nn,__FILE__,__LINE__);
    typelist[nn]=thistype;
    nn++;
  }
  *n = nn;
  
  return i;
}

void renum_atype(t_params plist[],t_topology *top,
		int *wall_atomtype,
		t_atomtype at,bool bVerbose)
{
  gpp_atomtype *ga = (gpp_atomtype *) at;
  
  int      i,j,k,l,mi,mj,nat,nrfp,ftype,ntype;
  t_param  *nbsnew;
  int      *typelist;
  real     *new_radius;
  real     *new_vol;
  real     *new_surftens;
  int      *new_atomnumber;
  
  ntype = get_atomtype_ntypes(at);
  snew(typelist,ntype);

  if (bVerbose)
    fprintf(stderr,"renumbering atomtypes...\n");

  /* Since the bonded interactions have been assigned now,
   * we want to reduce the number of atom types by merging 
   * ones with identical nonbonded interactions, in addition
   * to removing unused ones.
   *
   * With Generalized-Born electrostatics, or implicit solvent
   * we also check that the atomtype radius, effective_volume
   * and surface tension match.
   *
   * With QM/MM we also check that the atom numbers match
   */
  
  /* Get nonbonded interaction type */
  if (plist[F_LJ].nr > 0)
    ftype=F_LJ;
  else
    ftype=F_BHAM;
   
  /* Renumber atomtypes by first making a list of which ones are actually used.
   * We provide the list of nonbonded parameters so search_atomtypes
   * can determine if two types should be merged. 
   */    
  nat=0;
  for(i=0; (i<top->atoms.nr); i++) {
    top->atoms.atom[i].type=
      search_atomtypes(at,&nat,typelist,top->atoms.atom[i].type,
		       plist[ftype].param,ftype);
    top->atoms.atom[i].typeB=
      search_atomtypes(at,&nat,typelist,top->atoms.atom[i].typeB,
		       plist[ftype].param,ftype);
  }

  for(i=0; i<2; i++) {
    if (wall_atomtype[i] >= 0)
      wall_atomtype[i] = search_atomtypes(at,&nat,typelist,wall_atomtype[i],
					  plist[ftype].param,ftype);
  }

  snew(new_radius,nat);
  snew(new_vol,nat);
  snew(new_surftens,nat);
  snew(new_atomnumber,nat);  

  /* We now have a list of unique atomtypes in typelist */

  if (debug)
    pr_ivec(debug,0,"typelist",typelist,nat,TRUE);
    
  /* Renumber nlist */ 
  nbsnew = NULL;
  snew(nbsnew,plist[ftype].nr);

  nrfp  = NRFP(ftype);
  
  for(i=k=0; (i<nat); i++)
  {
    mi=typelist[i];
    for(j=0; (j<nat); j++,k++) 
    {
      mj=typelist[j];
      for(l=0; (l<nrfp); l++)
      {
        nbsnew[k].c[l]=plist[ftype].param[ntype*mi+mj].c[l];
      }
    }  
    new_radius[i]     = get_atomtype_radius(mi,at);
    new_vol[i]        = get_atomtype_vol(mi,at);
    new_surftens[i]   = get_atomtype_surftens(mi,at);
    new_atomnumber[i] = get_atomtype_atomnumber(mi,at);
  }
  
  for(i=0; (i<nat*nat); i++) {
    for(l=0; (l<nrfp); l++)
      plist[ftype].param[i].c[l]=nbsnew[i].c[l];
  }
  plist[ftype].nr=i;
  top->idef.atnr = nat;
  
  sfree(ga->radius);
  sfree(ga->vol);
  sfree(ga->surftens);
  sfree(ga->atomnumber);
  
  ga->radius     = new_radius;
  ga->vol        = new_vol;
  ga->surftens   = new_surftens;
  ga->atomnumber = new_atomnumber;
  
  ga->nr=nat;

  sfree(nbsnew);
  sfree(typelist);
}

void copy_atomtype_atomtypes(t_atomtype atype,t_atomtypes *atomtypes)
{
  gpp_atomtype *ga = (gpp_atomtype *) atype;
  int i,ntype;
  
  /* Copy the atomtype data to the topology atomtype list */
  ntype = get_atomtype_ntypes(atype);
  atomtypes->nr=ntype;
  snew(atomtypes->radius,ntype);
  snew(atomtypes->vol,ntype);
  snew(atomtypes->surftens,ntype);
  snew(atomtypes->atomnumber,ntype);


  for(i=0; i<ntype; i++) {
    atomtypes->radius[i]     = ga->radius[i];
    atomtypes->vol[i]        = ga->vol[i];
    atomtypes->surftens[i]   = ga->surftens[i];
    atomtypes->atomnumber[i] = ga->atomnumber[i];
  }
}

