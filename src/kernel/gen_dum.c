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
 * GRoups of Organic Molecules in ACtion for Science
 */
static char *SRCID_gen_dum_c = "$Id$";

#include <stdio.h>
#include <math.h>
#include <string.h>
#include "assert.h"
#include "gen_dum.h"
#include "smalloc.h"
#include "resall.h"
#include "add_par.h"
#include "vec.h"
#include "toputil.h"
#include "physics.h"
#include "index.h"
#include "names.h"

static void count_bonds(int atom, t_params *psb, char ***atomname, 
			int *nrbonds, int *nrHatoms, int Hatoms[], int *Heavy,
			int *nrheavies, int heavies[])
{
  int i,heavy,other,nrb,nrH,nrhv;
  
  /* find heavy atom bound to this hydrogen */
  heavy=NOTSET;
  for(i=0; (i<psb->nr) && (heavy==NOTSET); i++)
    if (psb->param[i].AI==atom)
      heavy=psb->param[i].AJ;
    else if (psb->param[i].AJ==atom)
      heavy=psb->param[i].AI;
  if (heavy==NOTSET)
    fatal_error(0,"unbound hydrogen atom %d",atom+1);
  /* find all atoms bound to heavy atom */
  other=NOTSET;
  nrb=0;
  nrH=0;
  nrhv=0;
  for(i=0; i<psb->nr; i++) {
    if (psb->param[i].AI==heavy)
      other=psb->param[i].AJ;
    else if (psb->param[i].AJ==heavy)
      other=psb->param[i].AI;
    if (other!=NOTSET) {
      nrb++;
      if (is_hydrogen(*(atomname[other]))) {
	Hatoms[nrH]=other;
	nrH++;
      } else {
	heavies[nrhv]=other;
	nrhv++;
      }
      other=NOTSET;
    }
  }
  *Heavy   =heavy;
  *nrbonds =nrb;
  *nrHatoms=nrH;
  *nrheavies=nrhv;
}

static void print_bonds(FILE *fp, int o2n[],
			int nrHatoms, int Hatoms[], int Heavy,
			int nrheavies, int heavies[])
{
  int i;
  
  fprintf(fp,"Found: %d Hatoms: ",nrHatoms);
  for(i=0; i<nrHatoms; i++)
    fprintf(fp," %d",o2n[Hatoms[i]]+1);
  fprintf(fp,"; %d Heavy atoms: %d",nrheavies+1,o2n[Heavy]+1);
  for(i=0; i<nrheavies; i++)
    fprintf(fp," %d",o2n[heavies[i]]+1);
  fprintf(fp,"\n");
}

static int get_atype(int atom, t_atoms *at, int nrtp, t_restp rtp[])
{
  int type;
  bool bNterm;
  int  j;
  t_restp *rtpp;
  
  if (at->atom[atom].m)
    type=at->atom[atom].type;
  else {
    /* get type from rtp */
    rtpp = search_rtp(*(at->resname[at->atom[atom].resnr]),nrtp,rtp);
    bNterm = is_protein(*(at->resname[at->atom[atom].resnr])) && 
      (at->atom[atom].resnr == 0);
    j=search_jtype(rtpp,*(at->atomname[atom]),bNterm);
    type=rtpp->atom[j].type;
  }
  return type;
}

static int nm2type(char *name, t_atomtype *atype)
{
  int tp,j;
  
  tp=NOTSET;
  for(j=0; (j < atype->nr) && (tp==NOTSET); j++)
    if (strcasecmp(name,*(atype->atomname[j])) == 0)
      tp=j;
  if (tp==NOTSET)
    fatal_error(0,"Dummy mass type (%s) not found "
		"in atom type database",name);
  return tp;
}

static real get_amass(int atom, t_atoms *at, int nrtp, t_restp rtp[])
{
  real mass;
  bool bNterm;
  int  j;
  t_restp *rtpp;
  
  if (at->atom[atom].m)
    mass=at->atom[atom].m;
  else {
    /* get mass from rtp */
    rtpp = search_rtp(*(at->resname[at->atom[atom].resnr]),nrtp,rtp);
    bNterm = is_protein(*(at->resname[at->atom[atom].resnr])) && 
      (at->atom[atom].resnr == 0);
    j=search_jtype(rtpp,*(at->atomname[atom]),bNterm);
    mass=rtpp->atom[j].m;
  }
  return mass;
}

static void add_dum_param(t_params plist[], t_params *newbonds, 
			  int dummy_type[], 
			  int Heavy, int nrHatoms, int Hatoms[], 
			  int nrheavies, int heavies[], char ***atomname)
{
  int i,j,ftype,other,moreheavy;
  bool bSwapParity;
  
  for(i=0; i<nrHatoms; i++) {
    ftype=dummy_type[Hatoms[i]];
    bSwapParity = (ftype<0);
    ftype=abs(ftype);
    dummy_type[Hatoms[i]]=ftype;
    switch (ftype) {
    case F_BONDS:
      if ( (nrheavies != 1) && (nrHatoms != 1) )
	fatal_error(0,"cannot make bond in add_dum_param for %d heavy atoms "
		    "and %d hydrogen atoms",nrheavies,nrHatoms);
      add_param(newbonds,Hatoms[i],heavies[0],NULL,NULL);
      break;
    case F_DUMMY3:
    case F_DUMMY3FD:
    case F_DUMMY3OUT:
      if (nrheavies < 2) 
	fatal_error(0,"Not enough heavy atoms (%d) for %s (min 3)",nrheavies+1,
		    interaction_function[dummy_type[Hatoms[i]]].name);
      add_dum3_param(&plist[ftype],Hatoms[i],Heavy,heavies[0],heavies[1],
		     bSwapParity);
      break;
    case F_DUMMY3FAD: {
      if (nrheavies > 1)
	moreheavy=heavies[1];
      else {
	/* find more heavy atoms */
	other=moreheavy=NOTSET;
	for(j=0; (j<plist[F_BONDS].nr) && (moreheavy==NOTSET); j++) {
	  if (plist[F_BONDS].param[j].AI==heavies[0])
	    other=plist[F_BONDS].param[j].AJ;
	  else if (plist[F_BONDS].param[j].AJ==heavies[0])
	    other=plist[F_BONDS].param[j].AI;
	  if ( (other != NOTSET) && (other != Heavy) ) 
	    moreheavy=other;
	}
	if (moreheavy==NOTSET)
	  fatal_error(0,"Unbound molecule part %d-%d",Heavy+1,Hatoms[0]+1);
      }
      add_dum3_param(&plist[ftype],Hatoms[i],Heavy,heavies[0],moreheavy,
		     bSwapParity);
      break;
    }
    case F_DUMMY4FD: {
      if (nrheavies < 3) 
	fatal_error(0,"Not enough heavy atoms (%d) for %s (min 4)",nrheavies+1,
		    interaction_function[dummy_type[Hatoms[i]]].name);
      add_dum4_param(&plist[ftype],  
		     Hatoms[0], Heavy, heavies[0], heavies[1], heavies[2]);
      break;
    }
    default:
      fatal_error(0,"can't use add_dum_param for interaction function %s",
		  interaction_function[dummy_type[Hatoms[i]]].name);
    } /* switch ftype */
  } /* for i */
}

static bool is_dum(int dummy_type)
{
  switch ( abs(dummy_type) ) {
  case NOTSET:
    return FALSE;
  case F_DUMMY3:
  case F_DUMMY3FD:
  case F_DUMMY3OUT:
  case F_DUMMY3FAD:
  case F_DUMMY4FD:
    return TRUE;
  default:
    return FALSE;
  }
}

void clean_dum_bonds(t_params *psb, int dummy_type[])
{
  int i,j;
  
  /* remove bonds with dummy atoms */
  for(i=j=0; (i<psb->nr); i++)
    if (!is_dum(dummy_type[psb->param[i].AI]) && 
	!is_dum(dummy_type[psb->param[i].AJ])) {
      memcpy(&(psb->param[j]),&(psb->param[i]),(size_t)sizeof(psb->param[0]));
      j++;
    }
  i=psb->nr-j;
  psb->nr=j;
  
  if (i>0)
    fprintf(stderr,"Removed %d bonds to dummy atoms, now %d bonds\n",i,j);
}

/* this seems overdone, but it is FOOLPROOF!!! */
static char atomnamesuffix[] = "1234";

void do_dummies(int nrtp, t_restp rtp[], 
		t_atomtype *atype, real mHmult, 
		t_atoms *at, t_symtab *symtab, rvec *x[], 
		t_params plist[], t_params *newbonds, 
		int *dummy_type[], int *cgnr[])
{
  int  ftype,i,j,k,i0,nrbonds,nrHatoms,Heavy,nrheavies,add_shift;
  int  nral,ndum,nadd,tpM,tpHeavy;
  int  Hatoms[4],heavies[4];
  bool bWARNING,bAddDumParam,bFirstWater;
  real mHtot,mtot,fact,fact2;
  rvec rpar,rperp,temp;
  char name[10],tpname[10];
  rvec *newx;
  int  *o2n,*newdummy_type,*newcgnr;
  t_atom *newatom;
  t_params *params;
  char   ***newatomname;
  
  if (debug) {
    printf("Searching for atoms to make dumies...\n");
    fprintf(debug,"# # # DUMMIES # # #\n");
  }
  
  bFirstWater=TRUE;
  ndum=0;
  nadd=0;
  /* we need a marker for which atoms whould *not* be renumbered afterwards */
  add_shift = 10*at->nr;
  /* make arrays where masses can be inserted into */ 
  snew(newx,at->nr); 
  snew(newatom,at->nr);
  snew(newatomname,at->nr);
  snew(newdummy_type,at->nr);
  snew(newcgnr,at->nr);
  /* make index array to tell where the atoms go to when masses are inse'd */
  snew(o2n,at->nr);
  for(i=0; i<at->nr; i++)
    o2n[i]=i;
  /* loop over all atoms */
  for(i=0; (i<at->nr); i++) {
    /* only process hydrogen atoms which are not already set */
    if ( ((*dummy_type)[i]==NOTSET) && is_hydrogen(*(at->atomname[i]))) {
      /* find heavy atom, count #bonds from it and #H atoms bound to it
	 and return H atom numbers (Hatoms) and heavy atom numbers (heavies) */
      count_bonds(i, &plist[F_BONDS], at->atomname, 
		  &nrbonds, &nrHatoms, Hatoms, &Heavy, &nrheavies, heavies);
      /* get Heavy atom type */
      tpHeavy=get_atype(Heavy,at,nrtp,rtp);
      strcpy(tpname,type2nm(tpHeavy,atype));
      bWARNING=FALSE;
      bAddDumParam=TRUE;
      /* nested if's which check nrHatoms, nrbonds and tpname */
      if (nrHatoms == 1) {
	switch(nrbonds) {
	case 2: /* -O-H */
	  (*dummy_type)[i]=F_BONDS;
	  break;
	case 3: /* =CH-, -NH- or =NH+- */
	  (*dummy_type)[i]=F_DUMMY3FD;
	  break;
	case 4: /* --CH- (tert) */
	  (*dummy_type)[i]=F_DUMMY4FD;
	break;
	default: /* nrbonds != 2, 3 or 4 */
	  bWARNING=TRUE;
	}
      } else if ( (nrHatoms == 2) && (nrbonds == 2) && 
		  (strcasecmp(tpname,"OW")==0) ) {
	bAddDumParam=FALSE; /* this is water: skip these hydrogens */
	if (bFirstWater) {
	  bFirstWater=FALSE;
	  if (debug)
	    fprintf(debug,
		    "Not converting hydrogens in water to dummy atoms\n");
	}
      } else if ( (nrHatoms == 2) && (nrbonds == 3) && 
		  (strcasecmp(tpname,"NL")!=0) ) {
	/* =CH2 or =NH2 */
	(*dummy_type)[Hatoms[0]] = F_DUMMY3FAD;
	(*dummy_type)[Hatoms[1]] =-F_DUMMY3FAD;
	
      } else if ( (nrHatoms == 2) && (nrbonds == 4) ) {
	/* -CH2- */
	(*dummy_type)[Hatoms[0]] = F_DUMMY3OUT;
	(*dummy_type)[Hatoms[1]] =-F_DUMMY3OUT;
	
      } else if ( ( (nrHatoms == 2) && (nrbonds == 3) && 
		    (strcasecmp(tpname,"NL")==0) ) || 
		  ( (nrHatoms == 3) && (nrbonds == 4) ) ) {
	/* this tells how to set the hydrogen atoms */
	int  Hat_dum_tp[3]   = { F_DUMMY3, F_DUMMY3OUT, F_DUMMY3OUT };
	bool Hat_swap_par[3] = { FALSE,    TRUE,        FALSE };
	
	bAddDumParam=FALSE; /* we'll do this ourselves! */
	/* -NH2 (umbrella), -NH3+ or -CH3 */
	(*dummy_type)[Heavy]       = F_DUMMY3;
	for (j=0; j<nrHatoms; j++)
	  (*dummy_type)[Hatoms[j]] = Hat_dum_tp[j];
	/* get dummy mass type from first char of heavy atom type (N or C) */
	sprintf(name,"M%cH3",tpname[0]);
	tpM=nm2type(name,atype);
	/* make space for 2 masses: shift all atoms starting with 'Heavy' */
	nadd+=2;
	srenew(newx,at->nr+nadd);
	i0=Heavy;
	if (debug) 
	  fprintf(stderr,"Inserting 2 dummy masses at %d\n",o2n[i0]+1);
	for(j=i0; j<at->nr; j++)
	  o2n[j]=j+nadd;
	i0+=nadd-2;
	srenew(newatom,at->nr+nadd);
	srenew(newatomname,at->nr+nadd);
	srenew(newdummy_type,at->nr+nadd);
	srenew(newcgnr,at->nr+nadd);
	/* calculate starting position for the masses */
	mHtot=0;
	/* get atom masses, and set Heavy and Hatoms mass to zero */
	for(j=0; j<nrHatoms; j++) {
	  mHtot += get_amass(Hatoms[j],at,nrtp,rtp);
	  at->atom[Hatoms[j]].m = at->atom[Hatoms[j]].mB = 0;
	}
	mtot = mHtot + get_amass(Heavy,at,nrtp,rtp);
	at->atom[Heavy].m = at->atom[Heavy].mB = 0;
	if (mHmult != 1.0)
	  mHtot *= mHmult;
	fact2=mHtot/mtot;
	fact=sqrt(fact2);
	/* generate vectors parallel and perpendicular to rotational axis:
	 * rpar  = Heavy -> Hcom
	 * rperp = Hcom  -> H1   */
	clear_rvec(rpar);
	for(j=0; j<nrHatoms; j++)
	  rvec_inc(rpar,(*x)[Hatoms[j]]);
	svmul(1.0/nrHatoms,rpar,rpar); /* rpar = ( H1+H2+H3 ) / 3 */
	rvec_dec(rpar,(*x)[Heavy]);    /*        - Heavy          */
	rvec_sub((*x)[Hatoms[0]],(*x)[Heavy],rperp);
	rvec_dec(rperp,rpar);          /* rperp = H1 - Heavy - rpar */
	/* calc mass positions */
	svmul(fact2,rpar,temp);
	for (j=0; (j<2); j++) /* xM = xN + fact2 * rpar +/- fact * rperp */
	  rvec_add((*x)[Heavy],temp,newx[i0+j]);
	svmul(fact,rperp,temp);
	rvec_inc(newx[i0  ],temp);
	rvec_dec(newx[i0+1],temp);
	/* set atom parameters for the masses */
	for(j=0; (j<2); j++) {
	  /* make name: "M??#" or "M?#" (? is atomname, # is number) */
	  name[0]='M';
	  for(k=0; (*at->atomname[Heavy])[k] && ( k < 2 ); k++ )
	    name[k+1]=(*at->atomname[Heavy])[k];
	  name[k+1]=atomnamesuffix[j];
	  name[k+2]='\0';
	  newatomname[i0+j]   = put_symtab(symtab,name);
	  newatom[i0+j].m     = newatom[i0+j].mB    = mtot/2;
	  newatom[i0+j].q     = newatom[i0+j].qB    = 0.0;
	  newatom[i0+j].type  = newatom[i0+j].typeB = tpM;
	  newatom[i0+j].ptype = eptAtom;
	  newatom[i0+j].resnr = at->atom[Heavy].resnr;
	  newatom[i0+j].chain = at->atom[Heavy].chain;
	  newdummy_type[i0+j] = NOTSET;
	  newcgnr[i0+j]       = (*cgnr)[Heavy];
	}
	/* add bonds between dummy masses and to heavies[0] */
	/* 'add_shift' says which atoms won't be renumbered afterwards */
	add_param(newbonds, heavies[0],   add_shift+i0,   NULL, NULL);
	add_param(newbonds, heavies[0],   add_shift+i0+1, NULL, NULL);
	add_param(newbonds, add_shift+i0, add_shift+i0+1, NULL, NULL);
	
	/* generate Heavy, H1, H2 and H3 from M1, M2 and heavies[0] */
	/* note that dummy_type cannot be NOTSET, because we just set it */
	add_dum3_param  (&plist[(*dummy_type)[Heavy]],
			 Heavy,     heavies[0], add_shift+i0, add_shift+i0+1, 
			 FALSE);
	for(j=0; j<nrHatoms; j++)
	  add_dum3_param(&plist[(*dummy_type)[Hatoms[j]]],
			 Hatoms[j], heavies[0], add_shift+i0, add_shift+i0+1, 
			 Hat_swap_par[j]);
      } else
	bWARNING=TRUE;
      if (bWARNING)
	fprintf(stderr,
		"Warning: cannot convert atom %d %s (bound to a heavy atom "
		"type %s with \n"
		"         %d bonds and %d bound hydrogens atoms) to dummy "
		"atom\n",
		i+1,*(at->atomname[i]),tpname,nrbonds,nrHatoms);
      if (bAddDumParam) {
	/* add dummy parameters to topology, 
	   also get rid of negative dummy_types */
 	add_dum_param(plist, newbonds, 
		      (*dummy_type), Heavy, nrHatoms, Hatoms,
 		      nrheavies, heavies, at->atomname);
	/* transfer mass of dummy atom to Heavy atom */
	for(j=0; j<nrHatoms; j++) 
	  if (is_dum((*dummy_type)[Hatoms[j]])) {
	    at->atom[Heavy].m += at->atom[Hatoms[j]].m;
	    at->atom[Heavy].mB = at->atom[Heavy].m;
	    at->atom[Hatoms[j]].m = at->atom[Hatoms[j]].mB = 0;
	  }
      }
      ndum+=nrHatoms;
      if (debug) {
	fprintf(debug,"atom %d: ",o2n[i]+1);
	print_bonds(debug,o2n,nrHatoms,Hatoms,Heavy,nrheavies,heavies);
      }
    }
  }
  
  /* add all original atoms to the new arrays, using o2n index array */
  for(i=0; i<at->nr; i++) {
    newatomname[o2n[i]]   = at->atomname[i];
    newatom[o2n[i]]       = at->atom[i];
    newdummy_type[o2n[i]] = (*dummy_type)[i];
    newcgnr[o2n[i]]       = (*cgnr)[i];
    copy_rvec((*x)[i],newx[o2n[i]]);
  }
  /* throw away old atoms */
  sfree(at->atom);
  sfree(at->atomname);
  sfree(*dummy_type);
  sfree(*cgnr);
  sfree(*x);
  /* put in the new ones */
  at->nr      += nadd;
  at->atom     = newatom;
  at->atomname = newatomname;
  *dummy_type   = newdummy_type;
  *cgnr         = newcgnr;
  *x            = newx;
  if (at->nr > add_shift)
    fatal_error(0,"Added impossible amount of dummy masses "
		"(%d on a total of %d atoms)\n",nadd,at->nr-nadd);
  
  if (debug)
    for(i=0; i<at->nr; i++)
      fprintf(debug,"%4d %4s %4d %4s %6d %-10s\n",i+1,*(at->atomname[i]),
	      at->atom[i].resnr,*(at->resname[at->atom[i].resnr]),
	      (*cgnr)[i],
	      ((*dummy_type)[i]==NOTSET) ? 
	      "NOTSET" : interaction_function[(*dummy_type)[i]].name);
  
  /* now renumber all the interactions, including 'newbonds': */
  for (ftype=0; ftype<=F_NRE; ftype++) {
    /* this is so we don't have to write the same code for newbonds */
    if (ftype==F_NRE) {
      params=newbonds;
      nral=NRAL(F_BONDS);
    } else {
      params=&(plist[ftype]);
      nral=NRAL(ftype);
    }
    if (debug)
      fprintf(debug,"Renumbering %d %s\n", params->nr, 
	      (ftype==F_NRE)?"New Bonds":interaction_function[ftype].longname);
    for (i=0; i<params->nr; i++) {
      for (j=0; j<nral; j++)
	if (params->param[i].a[j]>=add_shift) {
	  if (debug) fprintf(debug," [%u -> %u]",params->param[i].a[j],
			     params->param[i].a[j]-add_shift);
	  params->param[i].a[j]=params->param[i].a[j]-add_shift;
	} else {
	  if (debug) fprintf(debug," [%u -> %d]",params->param[i].a[j],
			     o2n[params->param[i].a[j]]);
	  params->param[i].a[j]=o2n[params->param[i].a[j]];
	}
      if (debug) fprintf(debug,"\n");
    }
  }
  
  /* clean up */
  sfree(o2n);
  
  /* tell the user what we did */
  fprintf(stderr,"Marked %d dummy atoms\n",ndum);
  fprintf(stderr,"Added %d dummy masses\n",nadd);
  fprintf(stderr,"Added %d new bonds\n",newbonds->nr);
}

void do_h_mass(t_params *psb, int dummy_type[], t_atoms *at, real mHmult)
{
  int i,j,a;
  
  /* loop over all atoms */
  for (i=0; i<at->nr; i++)
    /* adjust masses if i is hydrogen and not a dummy atom */
    if ( !is_dum(dummy_type[i]) && is_hydrogen(*(at->atomname[i])) ) {
      /* find bonded heavy atom */
      a=NOTSET;
      for(j=0; (j<psb->nr) && (a==NOTSET); j++) {
	/* if other atom is not a dummy, it is the one we want */
	if ( (psb->param[j].AI==i) && 
	     !is_dum(dummy_type[psb->param[j].AJ]) )
	  a=psb->param[j].AJ;
	else if ( (psb->param[j].AJ==i) && 
		  !is_dum(dummy_type[psb->param[j].AI]) )
	  a=psb->param[j].AI;
      }
      if (a==NOTSET)
	fatal_error(0,"Unbound hydrogen atom (%d) found while adjusting mass",
		    i+1);
      
      /* adjust mass of i (hydrogen) with mHmult
	 and correct mass of a (bonded atom) with same amount */
      at->atom[a].m -= (mHmult-1.0)*at->atom[i].m;
      at->atom[a].mB-= (mHmult-1.0)*at->atom[i].m;
      at->atom[i].m *= mHmult;
      at->atom[i].mB*= mHmult;
    }
}

void clean_dum_angles(t_params *psa, int natom, 
		      t_params *plist, int dummy_type[])
{
  int      ftype,i,j,parnr,k,l,m,n,ndum,kept_i,dumnral,dumtype;
  atom_id  atom,constr,at1,at2;
  atom_id  dumatoms[MAXATOMLIST];
  bool     bKeep,bUsed,bPresent,bAll3FAD,bFirstTwo;
  struct { int ftype,parnr; } *pindex;
  
  /* make index into dummy entries of plist: */
  snew(pindex,natom);
  for(ftype=0; (ftype<F_NRE); ftype++)
    if (interaction_function[ftype].flags & IF_DUMMY)
      for (parnr=0; (parnr<plist[ftype].nr); parnr++) {
	k=plist[ftype].param[parnr].AI;
	pindex[k].ftype=ftype;
	pindex[k].parnr=parnr;
      }
  
  dumnral=0;
  kept_i=0;
  for(i=0; (i<psa->nr); i++) { /* for all angles in the plist */
    bKeep=FALSE;
    bAll3FAD=TRUE;
    /* check if all dummies are constructed from the same atoms */
    ndum=0;
    for(k=0; (k<3) && !bKeep; k++) { /* for all atoms in the angle */
      atom = psa->param[i].a[k];
      if (is_dum(dummy_type[atom])) {
	ndum++;
	bAll3FAD = bAll3FAD && (pindex[atom].ftype == F_DUMMY3FAD);
	if (ndum==1) {
	  /* store construction atoms of first dummy */
	  dumnral=NRAL(pindex[atom].ftype)-1;
	  for(m=0; (m<dumnral); m++)
	    dumatoms[m]=
	      plist[pindex[atom].ftype].param[pindex[atom].parnr].a[m+1];
	} else 
	  /* check if this dummy is constructed from the same atoms */
	  if (dumnral == NRAL(pindex[atom].ftype)-1 )
	    for(m=0; (m<dumnral) && !bKeep; m++) {
	      bPresent=FALSE;
	      constr=
		plist[pindex[atom].ftype].param[pindex[atom].parnr].a[m+1];
	      for(n=0; (n<dumnral) && !bPresent; n++)
		if (constr == dumatoms[n])
		  bPresent=TRUE;
	      if (!bPresent)
		bKeep=TRUE;
	    }
	  else
	    bKeep=TRUE;
      }
    }
    
    /* keep all angles with no dummies in them or 
       with dummies with more than 3 constr. atoms */
    if ( ndum == 0 && dumnral > 3 )
      bKeep=TRUE;
    
    /* check if all non-dummy atoms are used in construction: */
    bFirstTwo=TRUE;
    for(k=0; (k<3) && !bKeep; k++) { /* for all atoms in the angle */
      atom = psa->param[i].a[k];
      if (!is_dum(dummy_type[atom])) {
	bUsed=FALSE;
	for(m=0; (m<dumnral) && !bUsed; m++)
	  if (atom == dumatoms[m]) {
	    bUsed=TRUE;
	    bFirstTwo = bFirstTwo && m<2;
	  }
	if (!bUsed)
	  bKeep=TRUE;
      }
    }
    
    if ( ! ( bAll3FAD && bFirstTwo ) )
      /* check if all constructing atoms are bound together */
      for (m=0; m<dumnral && !bKeep; m++) { /* all constr. atoms */
	at1 = dumatoms[m];
	at2 = dumatoms[(m+1) % dumnral];
	bPresent=FALSE;
	for (parnr=0; parnr<plist[F_BONDS].nr && !bPresent; parnr++)
	  /* all bonds until one matches */
	  bPresent = ( ( (plist[F_BONDS].param[parnr].AI == at1) &&
			 (plist[F_BONDS].param[parnr].AJ == at2) ) || 
		       ( (plist[F_BONDS].param[parnr].AI == at2) &&
			 (plist[F_BONDS].param[parnr].AJ == at1) ) );
	if (!bPresent)
	  bKeep=TRUE;
      }
    
    if ( bKeep ) {
      /* now copy the angle to the new array */
      memcpy(&(psa->param[kept_i]),
	     &(psa->param[i]),(size_t)sizeof(psa->param[0]));
      kept_i++;
    }
  }
  
  fprintf(stderr,"Removed %d angles with dummy atoms, now %d angles\n",
	  psa->nr-kept_i,kept_i);
  psa->nr=kept_i;
  
  /* clean up */
  sfree(pindex);
}

void clean_dum_dihs(t_params *psdih, int natom, char dihname[], 
		    t_params *plist, int dummy_type[])
{
  int      ftype,i,parnr,k,l,m,n,ndum,kept_i;
  atom_id  atom,constr;
  atom_id  dumatoms[3];
  bool     bKeep,bUsed,bPresent;
  struct { int ftype,parnr; } *pindex;
  
  /* make index into dummy entries of plist: */
  snew(pindex,natom);
  for(ftype=0; (ftype<F_NRE); ftype++)
    if (interaction_function[ftype].flags & IF_DUMMY)
      for (parnr=0; (parnr<plist[ftype].nr); parnr++) {
	k=plist[ftype].param[parnr].AI;
	pindex[k].ftype=ftype;
	pindex[k].parnr=parnr;
      }
  
  kept_i=0;
  for(i=0; (i<psdih->nr); i++) { /* for all dihedrals in the plist */
    bKeep=FALSE;
    /* check if all dummies are constructed from the same atoms */
    ndum=0;
    for(k=0; (k<4) && !bKeep; k++) { /* for all atoms in the dihedral */
      atom = psdih->param[i].a[k];
      if (is_dum(dummy_type[atom])) {
	ndum++;
	if (ndum==1) {
	  /* store construction atoms of first dummy */
	  for(m=1; (m<4); m++)
	    dumatoms[m-1]=
	      plist[pindex[atom].ftype].param[pindex[atom].parnr].a[m];
	  if (debug) {
	    fprintf(debug,"dih w. dum: %u %u %u %u\n",
		    psdih->param[i].AI+1,psdih->param[i].AJ+1,
		    psdih->param[i].AK+1,psdih->param[i].AL+1);
	    fprintf(debug,"dum %u from: %u %u %u\n",
		    atom+1,dumatoms[0]+1,dumatoms[1]+1,dumatoms[2]+1);
	  }
	} else 
	  /* check if this dummy is constructed from the same atoms */
	  for(m=1; (m<4) && !bKeep; m++) {
	    bPresent=FALSE;
	    constr=plist[pindex[atom].ftype].param[pindex[atom].parnr].a[m];
	    for(n=0; (n<3) && !bPresent; n++)
	      if (constr == dumatoms[n])
		bPresent=TRUE;
	    if (!bPresent)
	      bKeep=TRUE;
	  }
      }
    }
    
    /* keep all dihedrals with no dummies in them */
    if (ndum==0)
      bKeep=TRUE;
    
    /* check if all atoms in dihedral are either dummies, or used in 
       construction of dummies. If so, keep it, if not throw away: */
    for(k=0; (k<4) && !bKeep; k++) { /* for all atoms in the dihedral */
      atom = psdih->param[i].a[k];
      if (!is_dum(dummy_type[atom])) {
	bUsed=FALSE;
	for(m=0; (m<3) && !bUsed; m++)
	  if (atom == dumatoms[m])
	    bUsed=TRUE;
	if (!bUsed) {
	  bKeep=TRUE;
	  if (debug) fprintf(debug,"unused atom in dih: %u\n",atom+1);
	}
      }
    }
      
    if ( bKeep ) {
      memcpy(&(psdih->param[kept_i]),
	     &(psdih->param[i]),(size_t)sizeof(psdih->param[0]));
      kept_i++;
    }
  }
  
  fprintf(stderr,"Removed %d %s dihedrals with dummy atoms, "
	  "now %d %s dihedrals\n",
	  psdih->nr-kept_i,dihname,kept_i,dihname);
  psdih->nr=kept_i;
  
  /* clean up */
  sfree(pindex);
}

void do_dum_excl(t_block *excl, int dummy_type[])
{
  int     i,j,k;
  t_block new_excl;
  
  new_excl.nr=excl->nr;
  snew(new_excl.index,new_excl.nr+1);
  new_excl.nra=0;
  new_excl.a=NULL;
  k=0;
  
  for (i=0; (i < excl->nr); i++) {
    new_excl.index[i]=k;
    if (is_dum(dummy_type[i])) {
      /* make space */
      new_excl.nra += excl->index[i+1] - excl->index[i];
      srenew(new_excl.a, new_excl.nra+1);
      /* copy */
      for (j=excl->index[i]; (j < excl->index[i+1]); j++)
	new_excl.a[k++] = excl->a[j];
    }
  }
  /* terminate the t_block */
  new_excl.index[new_excl.nr] = k;
  
  /* throw away old stuff */
  sfree(excl->index);
  sfree(excl->a);
  
  /* give back new stuff */
  *excl=new_excl;
}
