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

static void add_dum_atoms(t_params plist[], t_params *newbonds, 
			  int dummy_type[], 
			  int Heavy, int nrHatoms, int Hatoms[], 
			  int nrheavies, int heavies[], char ***atomname,
			  int fbond)
{
  int i,j,ftype,other,moreheavy,bb;
  bool bSwapParity;
  
  for(i=0; i<nrHatoms; i++) {
    ftype=dummy_type[Hatoms[i]];
    bSwapParity = (ftype<0);
    dummy_type[Hatoms[i]] = ftype = abs(ftype);
    if (ftype == fbond) {
      if ( (nrheavies != 1) && (nrHatoms != 1) )
	fatal_error(0,"cannot make bond in add_dum_atoms for %d heavy atoms "
		    "and %d hydrogen atoms",nrheavies,nrHatoms);
      add_param(newbonds,Hatoms[i],heavies[0],NULL,NULL);
    } else {
      switch (ftype) {
      case F_DUMMY3:
      case F_DUMMY3FD:
      case F_DUMMY3OUT:
	if (nrheavies < 2) 
	  fatal_error(0,"Not enough heavy atoms (%d) for %s (min 3)",
		      nrheavies+1,
		      interaction_function[dummy_type[Hatoms[i]]].name);
	add_dum3_atoms(&plist[ftype],Hatoms[i],Heavy,heavies[0],heavies[1],
		       bSwapParity);
	break;
      case F_DUMMY3FAD: {
	if (nrheavies > 1)
	  moreheavy=heavies[1];
	else {
	  /* find more heavy atoms */
	  other=moreheavy=NOTSET;
	  for(j=0; (j<plist[fbond].nr) && (moreheavy==NOTSET); j++) {
	    if (plist[fbond].param[j].AI==heavies[0])
	      other=plist[fbond].param[j].AJ;
	    else if (plist[fbond].param[j].AJ==heavies[0])
	      other=plist[fbond].param[j].AI;
	    if ( (other != NOTSET) && (other != Heavy) ) 
	      moreheavy=other;
	  }
	  if (moreheavy==NOTSET)
	    fatal_error(0,"Unbound molecule part %d-%d",Heavy+1,Hatoms[0]+1);
	}
	add_dum3_atoms(&plist[ftype],Hatoms[i],Heavy,heavies[0],moreheavy,
		       bSwapParity);
	break;
      }
      case F_DUMMY4FD: {
	if (nrheavies < 3) 
	  fatal_error(0,"Not enough heavy atoms (%d) for %s (min 4)",
		      nrheavies+1,
		      interaction_function[dummy_type[Hatoms[i]]].name);
	add_dum4_atoms(&plist[ftype],  
		       Hatoms[0], Heavy, heavies[0], heavies[1], heavies[2]);
	break;
      }
      default:
	fatal_error(0,"can't use add_dum_atoms for interaction function %s",
		    interaction_function[dummy_type[Hatoms[i]]].name);
      } /* switch ftype */
    } /* else */
  } /* for i */
}

static bool is_dum(int dummy_type)
{
  if (dummy_type == NOTSET)
    return FALSE;
  switch ( abs(dummy_type) ) {
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
		int *dummy_type[], int *cgnr[],
		int fbonds)
{
  int  ftype,i,j,k,i0,nrbonds,nrHatoms,Heavy,nrheavies,add_shift;
  int  nral,ndum,nadd,tpM,tpHeavy;
  int  Hatoms[4],heavies[4],bb;
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
      count_bonds(i, &plist[fbonds], at->atomname, 
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
	  /* Bond type is set to be either G96BONDS or BONDS in pdb2top */
	  (*dummy_type)[i]=fbonds;
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
	int  Hat_dummy_type[3] = { F_DUMMY3, F_DUMMY3OUT, F_DUMMY3OUT };
	bool Hat_SwapParity[3] = { FALSE,    TRUE,        FALSE };
	
	bAddDumParam=FALSE; /* we'll do this ourselves! */
	/* -NH2 (umbrella), -NH3+ or -CH3 */
	(*dummy_type)[Heavy]       = F_DUMMY3;
	for (j=0; j<nrHatoms; j++)
	  (*dummy_type)[Hatoms[j]] = Hat_dummy_type[j];
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
	add_dum3_atoms  (&plist[(*dummy_type)[Heavy]],
			 Heavy,     heavies[0], add_shift+i0, add_shift+i0+1, 
			 FALSE);
	for(j=0; j<nrHatoms; j++)
	  add_dum3_atoms(&plist[(*dummy_type)[Hatoms[j]]],
			 Hatoms[j], heavies[0], add_shift+i0, add_shift+i0+1, 
			 Hat_SwapParity[j]);
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
 	add_dum_atoms(plist, newbonds, 
		      (*dummy_type), Heavy, nrHatoms, Hatoms,
 		      nrheavies, heavies, at->atomname,fbonds);
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
      nral=NRAL(fbonds);
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

/* some things used by gen_dums_* ,
   these are the default GROMACS bondlengths and angles 
   for aromatic ring systems: */
#define bR6 0.139
#define bR5 0.133
#define bCC 0.139
#define bCO 0.136
#define bCH 0.108
#define bNH 0.100
#define bOH 0.100
/* it does not make sense to change these angles: */
#define aR6  (DEG2RAD*120)
#define aR6H (M_PI-0.5*aR6)
#define aR5  (DEG2RAD*108)
#define aR5H (M_PI-0.5*aR5)
#define aCOH (DEG2RAD*109.5)
/* force constant for bonds that we add 
   (doesn't matter when we constrain, but we also want to do em) */
#define BONDFORCE 500000
/* cosine rule: a^2 = b^2 + c^2 - 2 b c cos(alpha) */
/* get a^2 when a, b and alpha are given: */
#define cosrule(b,c,alpha) ( sqr(b) + sqr(c) - 2*b*c*cos(alpha) )
/* get cos(alpha) when a, b and c are given: */
#define acosrule(a,b,c) ( (sqr(b)+sqr(c)-sqr(a))/(2*b*c) )

static int gen_dums_6ring(t_atoms *at, int *dummy_type[], t_params plist[], 
			  t_params *newbonds, int nrfound, int *ats, 
			  bool bDoHZ)
{
  /* these MUST correspond to the atnms array in do_dum_aromatics! */
  enum { atCG, atCD1, atHD1, atCD2, atHD2, atCE1, atHE1, atCE2, atHE2, 
	 atCZ, atHZ, atNR };
  
  int i,ndum;
  real a,b,c[MAXFORCEPARAM],dCGCE,tmp1,tmp2,mtot;
  /* CG, CE1 and CE2 stay and each get 1/3 of the total mass, 
     rest gets dummified */
  if (bDoHZ)
    assert(atNR == nrfound);
  
  /* bonds between CG, CE1 and CE2: */
  for(i=0; i<MAXFORCEPARAM; i++)
    c[i]=NOTSET;
  dCGCE = sqrt( cosrule(bR6,bR6,aR6) );
  c[0]=dCGCE;
  c[1]=BONDFORCE;
  add_param(newbonds,ats[atCG] ,ats[atCE1],c,NULL);
  add_param(newbonds,ats[atCG] ,ats[atCE2],c,NULL);
  add_param(newbonds,ats[atCE1],ats[atCE2],c,NULL);
  
  /* rest will be dummy3 */
  mtot=0;
  ndum=0;
  for(i=0; i<atNR; i++) {
    mtot+=at->atom[ats[i]].m;
    at->atom[ats[i]].m = at->atom[ats[i]].mB = 0;
    if (i!=atCG && i!=atCE1 && i!=atCE2 && (bDoHZ || i!=atHZ) ) {
      (*dummy_type)[ats[i]]=F_DUMMY3;
      ndum++;
    }
  }
  at->atom[ats[atCG]].m = at->atom[ats[atCG]].mB = 
    at->atom[ats[atCE1]].m = at->atom[ats[atCE1]].mB = 
    at->atom[ats[atCE2]].m = at->atom[ats[atCE2]].mB = mtot / 3;
  
  /* dummy3 construction: r_d = r_i + a r_ij + b r_ik */
  tmp1 = dCGCE*sin(DEG2RAD*60);
  tmp2 = bR6*cos(0.5*aR6) + tmp1;
  tmp1 *= 2;
  a = b = - bCH / tmp1;
  /* HE1 and HE2: */
  add_dum3_param(&plist[F_DUMMY3],
		 ats[atHE1],ats[atCE1],ats[atCE2],ats[atCG], a,b);
  add_dum3_param(&plist[F_DUMMY3],
		 ats[atHE2],ats[atCE2],ats[atCE1],ats[atCG], a,b);
  /* CD1, CD2 and CZ: */
  a = b = tmp2 / tmp1;
  add_dum3_param(&plist[F_DUMMY3],
		 ats[atCD1],ats[atCE2],ats[atCE1],ats[atCG], a,b);
  add_dum3_param(&plist[F_DUMMY3],
		 ats[atCD2],ats[atCE1],ats[atCE2],ats[atCG], a,b);
  add_dum3_param(&plist[F_DUMMY3],
		 ats[atCZ], ats[atCG], ats[atCE1],ats[atCE2],a,b);
  /* HD1, HD2 and HZ: */
  a = b = ( bCH + tmp2 ) / tmp1;
  add_dum3_param(&plist[F_DUMMY3],
		 ats[atHD1],ats[atCE2],ats[atCE1],ats[atCG], a,b);
  add_dum3_param(&plist[F_DUMMY3],
		 ats[atHD2],ats[atCE1],ats[atCE2],ats[atCG], a,b);
  if (bDoHZ)
    add_dum3_param(&plist[F_DUMMY3],
		   ats[atHZ], ats[atCG], ats[atCE1],ats[atCE2],a,b);
  
  return ndum;
}

static int gen_dums_phe(t_atoms *at, int *dummy_type[], t_params plist[], 
			 t_params *newbonds, int nrfound, int *ats)
{
  return gen_dums_6ring(at, dummy_type, plist, newbonds, nrfound, ats, TRUE);
}

static int gen_dums_trp(t_atomtype *atype, rvec *newx[],
			t_atom *newatom[], char ***newatomname[], 
			int *o2n[], int *newdummy_type[], int *newcgnr[],
			t_symtab *symtab, int *nadd, rvec x[], int *cgnr[],
			t_atoms *at, int *dummy_type[], t_params plist[], 
			t_params *newbonds, int nrfound, int *ats,
			int add_shift)
{
#define NMASS 2
  /* these MUST correspond to the atnms array in do_dum_aromatics! */
  enum { 
    atCB,  atCG,  atCD1, atHD1, atCD2, atNE1, atHE1, atCE2, atCE3, atHE3, 
    atCZ2, atHZ2, atCZ3, atHZ3, atCH2, atHH2, atNR };
  /* weights for determining the COM's of both rings (M1 and M2): */
  real mw[NMASS][atNR] = {
    {   0,     1,     1,     1,   0.5,     1,     1,   0.5,     0,     0,
        0,     0,     0,     0,     0,     0 },
    {   0,     0,     0,     0,   0.5,     0,     0,   0.5,     1,     1,
        1,     1,     1,     1,     1,     1 }
  };
  /* flag which sets atoms to positive or negative y-axis: */
  int py[atNR] = {
        1,     1,     1,     1,     1,    -1,    -1,    -1,     1,     1,
       -1,    -1,     1,     1,    -1,    -1
  };
  real xi[atNR],yi[atNR],xM[NMASS];
  int  atM[NMASS],tpM,i,i0,j,ndum;
  real mwtot,mtot,mM[NMASS],xcom,dCBM1,dCBM2,dM1M2;
  real a,b,c[MAXFORCEPARAM];
  rvec rx,rcom,temp;
  char name[10];
    
  
  /* keep CB and generate two dummy masses (one in each ring) on the
     symmetry axis (line through CD1 and middle of CZ3-CH bond) such that
     center of mass (COM) remains the same and the dummy mass for the
     five-ring is perpendicularly placed to the CG atom. This results
     in a somewhat larger moment of inertia, but gives two perpendicular
     construction vectors for the dummy atoms which makes calculating
     the dummy constructions easier. */
  
  assert(atNR == nrfound);
  
  /* get dummy mass type */
  tpM=nm2type("MTYR",atype);
  /* make space for 2 masses: shift all atoms starting with CB */
  *nadd+=NMASS;
  srenew(*newx,at->nr+*nadd);
  i0=ats[atCB];
  for(j=0; j<NMASS; j++)
    atM[j]=i0+j;
  if (debug)
    fprintf(stderr,"Inserting %d dummy masses at %d\n",NMASS,(*o2n)[i0]+1);
  for(j=i0; j<at->nr; j++)
    (*o2n)[j]=j+*nadd;
  srenew(*newatom,at->nr+*nadd);
  srenew(*newatomname,at->nr+*nadd);
  srenew(*newdummy_type,at->nr+*nadd);
  srenew(*newcgnr,at->nr+*nadd);
  
  /* calculate the relative positions of all atoms (along symmetry axis): */
  /* take x=0 at CD2/CE2: */
  xi[atCD2] = xi[atCE2]   = 0;
  /* first 5-ring: */
  xi[atCG]  = xi[atNE1]   = xi[atCD2] - sin(aR5)*bR5;
  xi[atHE1]               = xi[atNE1] - sin(aR5H-aR5)*bNH;
  xi[atCB]                = xi[atCG]  - sin(aR5H-aR5)*bCC;
  xi[atCD1]               = xi[atNE1] - sin(2*aR5-M_PI)*bR5;
  xi[atHD1]               = xi[atCD1] - bCH;
  /* then 6-ring: */
  xi[atCE3] = xi[atHE3] = 
    xi[atCZ2] = xi[atHZ2] = xi[atCD2] + sin(aR6)*bR6;
  xi[atCZ3] = xi[atCH2]   = xi[atCE3] + sin(2*aR6-M_PI)*bR6;
  xi[atHZ3] = xi[atHH2]   = xi[atCZ3] + sin(aR6H+M_PI-2*aR6)*bCH;
  /* also position perp. to symmetry axis: */
  yi[atCD1] = yi[atHD1]   = 0;
  yi[atCG]  = yi[atNE1]   = sin(0.5*aR5)*bR5;
  yi[atCB]                = yi[atCG]  + cos(aR5H-aR5)*bCC;
  yi[atHE1]               = yi[atNE1] + cos(aR5H-aR5)*bNH;
  yi[atCD2] = yi[atCE2] = 
    yi[atCZ3] = yi[atCH2] = 0.5*bR6;
  yi[atCE3] = yi[atCZ2]   = yi[atCD2] - cos(aR6)*bR6;
  yi[atHE3] = yi[atHZ2]   = yi[atCE3] + cos(aR6H-aR6)*bCH;
  yi[atHZ3] = yi[atHH2]   = yi[atCZ3] - cos(aR6)*bCH;
  /* now introduce positive or negative y: */
  for(i=0; i<atNR; i++)
    yi[i]*=py[i];
  
  /* first get COM: */
  xcom=0;
  for(j=0; j<NMASS; j++)
    mM[j]=0;
  for(i=0; i<atNR; i++)
    if (i!=atCB) {
      for(j=0; j<NMASS; j++)
	mM[j] += mw[j][i] * at->atom[ats[i]].m;
      xcom += at->atom[ats[i]].m * xi[i];
    }
  mtot=0;
  for(j=0; j<NMASS; j++)
    mtot += mM[j];
  xcom /= mtot;
  /* redefine coordinates to get zero at COM: */
  for(i=0; i<atNR; i++)
    xi[i] -= xcom;
  /* now we set M1 at the same 'x' as CB: */
  xM[0] = xi[atCB];
  /* then set M2 so that com is conserved (mM1 * xM1 + mM2 * xM2 = 0): */
  xM[1] = - mM[0] * xM[0] / mM[1];
  
  /* make a unitvector that defines our symmetry axis: */
  rvec_sub(x[ats[atCZ3]],x[ats[atCD2]],rx);
  unitv(rx,rx);        /* rx = ( CZ3 - CD2 ) / | CZ3 - CD2 | */
  /* make vector that points to origin (= COM): */
  rvec_add(x[ats[atCE2]],x[ats[atCD2]],rcom);
  svmul(HALF,rcom,rcom);
  svmul(xcom,rx,temp);
  rvec_inc(rcom,temp); /* rcom = 0.5 * ( CE2 + CD2 ) + xcom * rx */
  
  /* calc initial position for dummy masses: */
  for(j=0; j<NMASS; j++) {
    svmul(xM[j],rx,temp);
    rvec_add(rcom,temp,(*newx)[atM[j]]); /* rM = rcom + xM * rx */
  }
  
  /* set parameters for the masses */
  for(j=0; j<NMASS; j++) {
    sprintf(name,"MTRP%d",j+1);
    (*newatomname)[atM[j]]= put_symtab(symtab,name);
    (*newatom)[atM[j]].m     = (*newatom)[atM[j]].mB    = mM[j];
    (*newatom)[atM[j]].q     = (*newatom)[atM[j]].qB    = 0.0;
    (*newatom)[atM[j]].type  = (*newatom)[atM[j]].typeB = tpM;
    (*newatom)[atM[j]].ptype = eptAtom;
    (*newatom)[atM[j]].resnr = at->atom[i0].resnr;
    (*newatom)[atM[j]].chain = at->atom[i0].chain;
    (*newdummy_type)[atM[j]] = NOTSET;
    (*newcgnr)[atM[j]]       = (*cgnr)[i0];
  }
  /* renumber cgnr: */
  for(i=i0; i<at->nr; i++)
    (*cgnr)[i]++;

  /* bonds between CB, M1 and M2 */
  /* 'add_shift' says which atoms won't be renumbered afterwards */
  dCBM1 = yi[atCB];
  dM1M2 = xM[1]-xM[0];
  dCBM2 = sqrt( sqr(dCBM1) + sqr(dM1M2) );
  for(i=0; i<MAXFORCEPARAM; i++)
    c[i]=NOTSET;
  c[0]=dCBM1;
  c[1]=BONDFORCE;
  add_param(newbonds, ats[atCB],        add_shift+atM[0], c, NULL);
  c[0]=dCBM2;
  add_param(newbonds, ats[atCB],        add_shift+atM[1], c, NULL);
  c[0]=dM1M2;
  add_param(newbonds, add_shift+atM[0], add_shift+atM[1], c, NULL);
  
  /* rest will be dummy3 */
  ndum=0;
  for(i=0; i<atNR; i++)
    if (i!=atCB) {
      at->atom[ats[i]].m = at->atom[ats[i]].mB = 0;
      (*dummy_type)[ats[i]] = F_DUMMY3;
      ndum++;
    }
  
  /* now define all dummies from M1, M2, CB, ie:
     r_d = r_M1 + a r_M!M2 + b r_M1_CB */
  for(i=0; i<atNR; i++)
    if ( (*dummy_type)[ats[i]] == F_DUMMY3 )
      add_dum3_param(&plist[F_DUMMY3],
		     ats[i],add_shift+atM[0],add_shift+atM[1],ats[atCB],
		     (xi[i]-xM[0])/dM1M2,yi[i]/dCBM1);
  
  return ndum;
}

static int gen_dums_tyr(t_atoms *at, int *dummy_type[], t_params plist[], 
			t_params *newbonds, int nrfound, int *ats)
{
  int ndum,i;
  real c[MAXFORCEPARAM],dCGCE,dCEOH,dCGHH,tmp1;
  
  /* these MUST correspond to the atnms array in do_dum_aromatics! */
  enum { atCG, atCD1, atHD1, atCD2, atHD2, atCE1, atHE1, atCE2, atHE2, 
	 atCZ, atOH, atHH, atNR };
  /* CG, CE1, CE2 (as in general 6-ring) and OH and HH stay, 
     rest gets dummified.
     Now we have two linked triangles with one improper keeping them flat */
  assert(atNR == nrfound);

  /* first do 6 ring as default, but without HZ (we don't have that): */
  ndum = gen_dums_6ring(at, dummy_type, plist, newbonds, nrfound, ats, FALSE);
  
  /* bonds between CE1, CE2 and OH */
  for(i=0; i<MAXFORCEPARAM; i++)
    c[i]=NOTSET;
  dCGCE = sqrt( cosrule(bR6,bR6,aR6) );
  dCEOH = sqrt( cosrule(bR6,bCO,aR6) );
  c[0]=dCEOH;
  c[1]=BONDFORCE;
  add_param(newbonds,ats[atCE1],ats[atOH],c,NULL);
  add_param(newbonds,ats[atCE2],ats[atOH],c,NULL);
  /* assume we also want the COH angle constrained: */
  tmp1 = bR6*cos(0.5*aR6) + dCGCE*sin(DEG2RAD*60) + bCO;
  dCGHH = sqrt( cosrule(tmp1,bOH,aCOH) );
  c[0]=dCGHH;
  (*dummy_type)[ats[atHH]]=F_BONDS;
  add_param(newbonds,ats[atCG],ats[atHH],c,NULL);
  return ndum;
}
      
static int gen_dums_his(t_atoms *at, int *dummy_type[], t_params plist[], 
			t_params *newbonds, int nrfound, int *ats)
{
  int ndum,i;
  real a,b,c[MAXFORCEPARAM],alpha,dGE,dCENE,mtot,tmp1;
  real xCG, yE, mE, mG;
  real xG, xD, yD, xE, xHD, yHD, xHE, yHE;
  
  /* these MUST correspond to the atnms array in do_dum_aromatics! */
  enum { atCG, atND1, atHD1, atCD2, atCE1, atNE2, atHE2, atNR };
  /* CG, CE1 and NE2 stay, each gets part of the total mass,
     rest gets dummified */
  assert(atNR == nrfound || atNR == nrfound+1);
  
  /* bonds between CG, CE1 and NE1 */
  for(i=0; i<MAXFORCEPARAM; i++)
    c[i]=NOTSET;
  c[1]=BONDFORCE;
  dGE   = sqrt( cosrule(bR5,bR5,aR5) );
  dCENE = bR5;
  c[0]=dGE;
  add_param(newbonds,ats[atCG] ,ats[atCE1],c,NULL);
  add_param(newbonds,ats[atCG] ,ats[atNE2],c,NULL);
  c[0]=dCENE;
  add_param(newbonds,ats[atCE1],ats[atNE2],c,NULL);
  
  /* rest will be dummy3 */
  mtot=0;
  ndum=0;
  for(i=0; i<atNR; i++) 
    if (ats[i]!=NOTSET) {
      mtot+=at->atom[ats[i]].m;
      if (i!=atCG && i!=atCE1 && i!=atNE2) {
	at->atom[ats[i]].m = at->atom[ats[i]].mB = 0;
	(*dummy_type)[ats[i]]=F_DUMMY3;
	ndum++;
      }
    }
  assert(ndum+3 == nrfound);
  
  /* distribute mass so that com stays the same */
  xG = bR5 * sin(0.5*aR5) / sin(aR5);
  xE = xG * sin(0.5*aR5);
  mE = mtot * xG / ( 2 * ( xG + xE ) );
  mG = mtot - 2 * mE;
  at->atom[ats[atCG]].m = at->atom[ats[atCG]].mB = mG;
  at->atom[ats[atCE1]].m = at->atom[ats[atCE1]].mB = 
    at->atom[ats[atNE2]].m = at->atom[ats[atNE2]].mB = mE;
  
  /* these we need for all the atoms: */
  alpha = acos( acosrule(dGE,dGE,dCENE) ); /* angle CG-NE2-CE1 */
  /* we define 'x' and 'y' perp to the three construction vectors: */
  xCG = sin(alpha)*dGE;   /* perp to r_CE1-NE2 */
  yE  = sin(alpha)*dCENE; /* perp to r_CG-CE1 or r_CG_NE2 */
  
  /* HE2 */
  if (ats[atHE2]!=NOTSET) {
    yHE = sin(M_2PI-alpha-aR5H)*bNH;
    xHE = sin(aR5H)*bNH;
    a = - yHE / yE;
    b = - xHE / xCG;
    add_dum3_param(&plist[F_DUMMY3],
		   ats[atHE2],ats[atNE2],ats[atCE1],ats[atCG],a,b);
  }
  
  /* ND1, CD2 */
  xD = sin(aR5)*bR5;
  yD = yE + sin(alpha+aR5-M_PI)*bR5;
  a = yD / yE;
  b = xD / xCG;
  add_dum3_param(&plist[F_DUMMY3],
		 ats[atND1],ats[atNE2],ats[atCE1],ats[atCG],a,b);
  add_dum3_param(&plist[F_DUMMY3],
		 ats[atCD2],ats[atCE1],ats[atNE2],ats[atCG],a,b);
  
  /* HD1 */
  if (ats[atHD1]!=NOTSET) {
    xHD = xD + sin(aR5H-aR5)*bNH;
    yHD = yD + sin(alpha+aR5-aR5H)*bNH;
    a = yHD / yE;
    b = xHD / xCG;
    add_dum3_param(&plist[F_DUMMY3],
		   ats[atHD1],ats[atNE2],ats[atCE1],ats[atCG],a,b); 
  }
  
  return ndum;
}

void do_dum_aromatics(int nrtp, t_restp rtp[], t_atomtype *atype, 
		      t_atoms *at, t_symtab *symtab, rvec *x[], 
		      t_params plist[], t_params *newbonds,
		      int *dummy_type[], int *cgnr[])
{
#define MAXATOMSPERRESIDUE 16
  int i,j,k,i0,whatres,ndum,nadd,nrfound,needed,resnr,add_shift,nral,ftype;
  int *o2n,*newdummy_type,*newcgnr,ats[MAXATOMSPERRESIDUE];
  rvec *newx;
  t_atom *newatom;
  char ***newatomname;
  char *resnm;
  t_params *params;
  
  /* do_dum_aromatics will specifically convert atoms in 
     PHE, TRP, TYR and HIS to a construction of dummy atoms from 3 atoms */
  enum { resPHE, resTRP, resTYR, resHIS, resNR };
  char *resnms[resNR]   = { "PHE", "TRP", "TYR", "HIS" };
  bool bPartial[resNR]  = { FALSE, FALSE, FALSE, TRUE };
  /* the atnms for every residue MUST correspond to the enums in the 
     gen_dums_* (one for each residue) routines! */
  /* also the atom names in atnms MUST be in the same order as in the .rtp! */
  char *atnms[resNR][MAXATOMSPERRESIDUE+1] = { 
    { "CG", /* PHE */
      "CD1", "HD1", "CD2", "HD2", 
      "CE1", "HE1", "CE2", "HE2", 
      "CZ", "HZ", NULL },
    { "CB", /* TRP */
      "CG",
      "CD1", "HD1", "CD2", 
      "NE1", "HE1", "CE2", "CE3", "HE3", 
      "CZ2", "HZ2", "CZ3", "HZ3", 
      "CH2", "HH2", NULL },
    { "CG", /* TYR */
      "CD1", "HD1", "CD2", "HD2", 
      "CE1", "HE1", "CE2", "HE2", 
      "CZ", "OH", "HH", NULL },
    { "CG", /* HIS */
      "ND1", "HD1", "CD2", 
      "CE1", "NE2", "HE2", NULL }
  };
  
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
  resnr=NOTSET;
  for(i=0; (i<at->nr); i++) {
    /* NOTE: at the end of this loop, i is skipped to the next residue! */
    if (at->atom[i].resnr != resnr) {
      resnr=at->atom[i].resnr;
      resnm=*(at->resname[resnr]);
    }
    /* don't waste our effort on DNA, water etc. */
    if (is_protein(*(at->resname[resnr]))) {
      whatres=NOTSET;
      for(j=0; j<resNR && whatres==NOTSET; j++) {
	if ( ( !bPartial[j] &&
	       (strcasecmp(resnm,resnms[j])==0) ) ||
	     ( bPartial[j] && 
	       (strncasecmp(resnm,resnms[j],strlen(resnms[j]))==0) ) ) {
	  whatres=j;
	  /* get atoms we will be needing */
	  nrfound=0;
	  for (k=0; atnms[j][k]; k++) {
	    ats[k]=NOTSET;
	    i0=i;
	    while (i<at->nr && at->atom[i].resnr==resnr && ats[k]==NOTSET) {
	      if (strcasecmp(*(at->atomname[i]),atnms[j][k])==0) {
		ats[k]=i;
		nrfound++;
	      }
	      i++;
	    }
	    /* if nothing found, search next atom from same point */
	    if (ats[k]==NOTSET)
	      i=i0;
	  }
	  /* now k is number of atom names in atnms[j] */
	  if (j==resHIS)
	    needed = k-1;
	  else
	    needed = k;
	  if (nrfound<needed)
	    fatal_error(0,"not enough atoms found (%d, need %d) in "
			"residue %s %d while\n"
			"             generating aromatics dummy construction",
			nrfound,needed,resnm,resnr+1);
	}
      }
      /* the enums for every residue MUST correspond to atnms[residue] */
      switch (whatres) {
      case resPHE: 
	ndum+=gen_dums_phe(at, dummy_type, plist, newbonds, nrfound, ats);
	break;
      case resTRP: 
	ndum+=gen_dums_trp(atype, &newx, &newatom, &newatomname, &o2n, 
			   &newdummy_type, &newcgnr, symtab, &nadd, *x, cgnr,
			   at, dummy_type, plist, newbonds, nrfound, ats,
			   add_shift);
	break;
      case resTYR: 
	ndum+=gen_dums_tyr(at, dummy_type, plist, newbonds, nrfound, ats);
	break;
      case resHIS: 
	ndum+=gen_dums_his(at, dummy_type, plist, newbonds, nrfound, ats);
	break;
      case NOTSET:
	/* this means this residue won't be processed */
	break;
      default:
	fatal_error(0,"DEATH HORROR in %s at %d",__FILE__,__LINE__);
      }
    }
    /* skip to next residue */
    while(i<at->nr && at->atom[i].resnr==resnr)
      i++;
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
  fprintf(stderr,"Marked %d dummy atoms for aromatics\n",ndum);
  fprintf(stderr,"Added %d dummy masses for aromatics\n",nadd);
  fprintf(stderr,"Added %d new bonds for aromatics\n",newbonds->nr);
  
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
		      t_params *plist, int dummy_type[],
		      int fbonds)
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
	for (parnr=0; (parnr<plist[fbonds].nr) && !bPresent; parnr++)
	  /* all bonds until one matches */
	  bPresent = ( ( (plist[fbonds].param[parnr].AI == at1) &&
			 (plist[fbonds].param[parnr].AJ == at2) ) || 
		       ( (plist[fbonds].param[parnr].AI == at2) &&
			 (plist[fbonds].param[parnr].AJ == at1) ) );
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
