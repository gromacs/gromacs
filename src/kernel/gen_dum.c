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
#include "gen_dum.h"
#include "smalloc.h"
#include "resall.h"
#include "add_par.h"
#include "vec.h"
#include "toputil.h"
#include "physics.h"
#include "index.h"

void do_dummies(t_atoms *at,t_atomtype *atype,t_symtab *symtab,
		int nrtp,t_restp rtp[],rvec **x,
		int nddb, t_dumblock *ddb, bool **is_dummy,real mHmult)
{
  int     i,j,k,i0,n,m,l,resnr,add,nadd,ndum;
  int     na[4];
  ushort  type,tpM;
  char    name[STRLEN],tpstr[STRLEN],bname[STRLEN],*bres;
  bool    *bProcessed,bWild,bNterm;
  t_restp *rtpp=NULL;

  if (debug) printf("Searching for atoms to make dumies...\n");
  
  bWild=FALSE;
  snew(bProcessed,at->nr);
  /* loop over all entries in the .ddb database (atomtypes) */
  add=0;
  ndum=0;
  for (n=0; (n<nddb); n++) {
    if (debug) fprintf(debug,"processing block %s\n",ddb[n].bname);
    if (bWild)
      fatal_error(0,"Found block in .ddb (%s) after wildcard block",
		  ddb[n].bname);
    strcpy(bname,ddb[n].bname);
    bWild=(strcmp(bname,"*")==0);
    bres=NULL;
    if (bres=strchr(bname,':')) {
      bres[0]='\0';
      bres++;
      if (debug) fprintf(debug,"atomtype %s residue name %s\n",bname,bres);
    }
    resnr=-1;
    /* loop over all atoms */
    for(i=0; (i<at->nr); i++) {
      /* look up this residue in rtp */
      if (at->atom[i].resnr != resnr) {
	resnr=at->atom[i].resnr;
	rtpp=search_rtp(*(at->resname[resnr]),nrtp,rtp);
	bNterm=is_protein(*(at->resname[resnr])) && (resnr == 0);
      }
      /* look up type for this atom */
      if (at->atom[i].m == 0) {
	j=search_jtype(rtpp,*(at->atomname[i]),bNterm);
	type = rtpp->atom[j].type;
      } else
	type = at->atom[i].type;
      strcpy(name,type2nm(type,atype));
      /* check if this atom has not yet been processed 
       * and if the type matches current .ddb block (bname) */
      if ( !bProcessed[i] && (bWild || strcmp(name,bname)==0)  && 
	   ( !bres || (strcmp(*(at->resname[resnr]),bres)==0) ) ) {
	if (debug) 
	  fprintf(debug,"Found atom: %s%d %d %s (%s)\n", 
		  *(at->resname[resnr]),resnr+1,i+1,*(at->atomname[i]),name);
	/* set marker for dummies */
	if ( !bWild || 
	     (strncasecmp(/**/
			  ddb[n].dum[0].na[0],*(at->atomname[i]),
			  strlen(ddb[n].dum[0].na[0]))==0)) {
	  j=i;
	  m=0;
	  if (debug) fprintf(debug,"Dummies: ");
	  while ( (j<at->nr) && (m<ddb[n].ndum) &&
		  ( !is_hydrogen(ddb[n].dum[m].na[0]) || 
		    (j==i) || is_hydrogen(*(at->atomname[j]))) ) {
	    strcpy(tpstr,ddb[n].dum[m].na[0]);
	    if (debug && bProcessed[j]) fprintf(debug,"P%d ",j+1);
	    if (!bProcessed[j] && 
		(strncasecmp(tpstr,*(at->atomname[j]),strlen(tpstr))==0) ) {
	      if (debug) fprintf(debug,"%d:'%s' ",j+1,tpstr);
	      bProcessed[j]=TRUE;
	      (*is_dummy)[j]=TRUE;
	      ndum++;
	      m++;
	    }
	    j++;
	  }  
 	  if (debug) fprintf(debug,"\n");
	  if (m<ddb[n].ndum)
	    fprintf(stderr,
		    "Skipped %d dummy atom%s on atom type %s in residue %s%d "
		    "(probably ok)\n",
		    ddb[n].ndum-m,(ddb[n].ndum-m==1)?"":"s",
		    ddb[n].bname,*(at->resname[resnr]),resnr);
	}
	
	/* add masses to this atom */
	if (!bWild)
	  for (m=0; (m<ddb[n].nmass); m++) {
	    /* find type nr for dummy mass */
	    for(tpM=0; (tpM<atype->nr); tpM++)
	      if (strcasecmp(ddb[n].mass[m].mtype,
			     *(atype->atomname[tpM])) == 0)
		break;
	    if (tpM == atype->nr)
	      fatal_error(0,"Atom type %s not found in database "
			  "while adding dummy masses\n",ddb[n].mass[m].mtype);
	    nadd = ddb[n].mass[m].nm;
	    add+=nadd;
	    at->nr+=nadd;
	    i0=i;
	    i+=nadd;
	    srenew(at->atom,at->nr);
	    srenew(at->atomname,at->nr);
	    srenew(*x,at->nr);
	    srenew(*is_dummy,at->nr);
	    srenew(bProcessed,at->nr);
	    for(j=at->nr-1; (j>=i); j--) {
	      at->atom[j]=at->atom[j-nadd];
	      at->atomname[j]=at->atomname[j-nadd];
	      copy_rvec((*x)[j-nadd],(*x)[j]);
	      (*is_dummy)[j]=(*is_dummy)[j-nadd];
	      bProcessed[j]=bProcessed[j-nadd];
	    }
	    
	    /* find Mass control atoms */
	    j=i;
	    l=0;
	    while ((j<at->nr) && (l<ddb[n].ndum)) {
	      if (strncasecmp(
			      ddb[n].dum[l].na[0],*(at->atomname[j]),
			      strlen(ddb[n].dum[l].na[0]))==0) {
		na[l]=j;  /* set control atoms */
		l++;
	      }
	      j++;
	    }
	    /* 'ndum' should vary depending on 'mass[m].tp' */
	    if (l<ddb[n].ndum)
	      fatal_error(0,"Atom %s not found while adding dummy mass %s "
			  "(%s)\n",ddb[n].dum[l].na[0],
			  ddb[n].bname,ddb[n].mass[m].mname);
	    if (l==3) { /* -NH2 (tetraedical) i.s.o. -NH3+ */
	      ddb[n].mass[m].tp++;
	      na[3]=-1;
	    }
	    /* now calculate positions of dummy masses */
	    if (debug) 
	      fprintf(debug,"Masses (%d): %d (%g) %d (%g) %d (%g) %d (%g)\n",
		      nadd,
		      na[0]+1,at->atom[na[0]].m,na[1]+1,at->atom[na[1]].m,
		      na[2]+1,at->atom[na[2]].m,na[3]+1,at->atom[na[3]].m);
	    switch (ddb[n].mass[m].tp) {
	    case 1: {
	      real mass[4],mHtot,mtot,fact,fact2;
	      rvec xN,xH[4],xM[2],rpar,rperp,temp;
	      /* note: xH[0] is not used */
	      
	      if (nadd!=2) fatal_error(0,"Can only add 2 (not %d) dummy "
				       "masses of type 1\n",nadd);
	      
	      copy_rvec((*x)[na[0]],xN);
	      for (l=1; (l<4); l++)
		copy_rvec((*x)[na[l]],xH[l]);
	      
	      /* calc masses */
	      mtot=0;
	      for(l=0; (l<4); l++) {
		if (at->atom[na[l]].m)
		  mass[l]=at->atom[na[l]].m;
		else {
		  bool bNterm;
		  int  resnr;
		  /* get mass from rtp */
		  resnr = at->atom[na[l]].resnr;
		  bNterm=is_protein(*(at->resname[resnr])) && (resnr==0);
		  j=search_jtype(rtpp,*(at->atomname[na[l]]),bNterm);
		  mass[l]=rtpp->atom[j].m;
		}
		mtot+=mass[l];
	      }
	      if (ddb[n].mass[m].c[3]==NOTSET) { /* c[3] = mH */
		mHtot=0;
		for(l=1; (l<4); l++)
		  mHtot+=mass[l];
		if (mHmult != 1.0) mHtot*=mHmult;
		ddb[n].mass[m].c[3]=mHtot/3;
	      } else {
		mHtot=3*ddb[n].mass[m].c[3];
		if (mHmult != 1.0) mHtot*=mHmult;
	      }
	      fact2=mHtot/mtot;
	      fact=sqrt(fact2);
	      if (debug) fprintf(debug,"Masses: %g %g %g %g\n",
				 mass[0],mass[1],mass[2],mass[3]);
	      
	      /* generate vectors rpar = N->Hcom and rperp = Hcom->H1 */
	      rvec_add(xH[1],xH[2],rpar); /* rpar = H1 + H2 */
	      rvec_inc(rpar,xH[3]);       /*           + H3 */
	      svmul(1.0/3.0,rpar,rpar);   /*        / 3     */
	      rvec_dec(rpar,xN);          /*        - N     */
	      rvec_sub(xH[1],xN,rperp); /* rperp = H1 - N   */
	      rvec_dec(rperp,rpar);     /*            - rpar */

	      /* calc mass positions */	      
	      svmul(fact2,rpar,temp);
	      for (k=0; (k<2); k++)
		rvec_add(xN,temp,xM[k]); /* xM = xN + fact2 * rpar */
	      svmul(fact,rperp,temp);
	      rvec_inc(xM[0],temp);      /*       +/- fact * rperp */
	      rvec_dec(xM[1],temp);
	      
	      /* copy new positions */
	      for(k=0; (k<2); k++)
		copy_rvec(xM[k],(*x)[i0+k]);
	    }
	    break;
	    default:
	      fatal_error(0,"Invalid dummy mass type (%d) in do_dummies \n",
			  ddb[n].mass[m].tp);
	    } /* end switch */
	    
	    for(j=i0; (j<i); j++) {
	      sprintf(name,"%s%c",ddb[n].mass[m].mname,'1'+(j-i0));
	      at->atomname[j]=put_symtab(symtab,name);
	      at->atom[j].m=at->atom[j].mB=-1; /* will be calc'ed in pdb2top */
	      at->atom[j].q=at->atom[j].qB=0.0;
	      at->atom[j].type=at->atom[j].typeB=tpM;
	      at->atom[j].ptype=eptAtom;
	      at->atom[j].resnr=at->atom[i0].resnr;
	      (*is_dummy)[j]=FALSE;
	      bProcessed[j]=TRUE;
	    }
	  }
	
	/* mark atoms in angle constraints as processed */
	if (!bWild) {
	  m=0;
	  j=i;
	  if (debug) fprintf(debug,"Angle:");
	  while ( (j<at->nr) && (m<ddb[n].nang) ) {
	    if (strncasecmp(
			    ddb[n].ang[m].na[0],*(at->atomname[j]),
			    strlen(ddb[n].ang[m].na[0]))==0) {
	      if (debug) fprintf(debug," %s %d",ddb[n].ang[m].na[0],j+1);
	      bProcessed[j]=TRUE;
	      m++;
	    }
	    j++;
	  }
	  if (debug) fprintf(debug,"\n");
	  if (m<ddb[n].nang)
	    fatal_error(0,"Found only %d atoms for angle constraints, "
			"expected %d\n",m,ddb[n].nang);
	}
	/* mark atom i as done */
 	bProcessed[i]=TRUE;
      }
    }
  }
  
  fprintf(stderr,"Marked %d dummy atoms\n",ndum);
  fprintf(stderr,"Added %d dummy masses\n",add);
  if (debug) {
    fprintf(debug,"%5s%3s %3s %4s %s\n","res","nr","atom","name","dummy");
    for (i=0; (i<at->nr); i++) 
      fprintf(debug,"%5s%3d %4d %4s %s\n",
	      *(at->resname[at->atom[i].resnr]),at->atom[i].resnr+1,i+1,
	      *(at->atomname[i]),(*is_dummy)[i] ? "dummy atom" : "atom");
  }
  
  /* clean up */
  sfree(bProcessed);
}

void do_h_mass(t_params *psb, bool is_dum[], t_atoms *at, real mHmult)
{
  int i,j,a;

  /* loop over all atoms */
  for (i=0; i<at->nr; i++)
    /* adjust masses if i is hydrogen and not a dummy atom */
    if ( !is_dum[i] && is_hydrogen(*(at->atomname[i])) ) {
      /* find bonded heavy atom */
      a=NOTSET;
      for(j=0; (j<psb->nr) && (a==NOTSET); j++) {
	/* if other atom is not a dummy, it is the one we want */
	if ( (psb->param[j].AI==i) && !is_dum[psb->param[j].AJ] )
	  a=psb->param[j].AJ;
	else if ( (psb->param[j].AJ==i) && !is_dum[psb->param[j].AI] )
	  a=psb->param[j].AI;
      }
      if (a==NOTSET)
	fatal_error(0,"Unbound hydrogen atom (%d) found while adjusting mass",i+1);
      
      /* adjust mass of i (hydrogen) with mHmult
	 and correct mass of a (bonded atom) with same amount */
      at->atom[a].m -= (mHmult-1.0)*at->atom[i].m;
      at->atom[a].mB-= (mHmult-1.0)*at->atom[i].m;
      at->atom[i].m *= mHmult;
      at->atom[i].mB*= mHmult;
    }
}

void clean_dum_angles(t_params *ps, t_params *plist, bool *is_dum)
{
  int i,j,k,ndum;

  j=0;
  for(i=0; (i<ps->nr); i++) {
    /* count number of dummies in this angle: */
    ndum=0;
    for(k=0; (k<3); k++)
      if (is_dum[ps->param[i].a[k]])
	ndum++;
    /* if 0 or 1 dummies, keep angle: */
    if ( ndum<2 ) {
      memcpy(&(ps->param[j]),&(ps->param[i]),(size_t)sizeof(ps->param[0]));
      j++;
    }
  }
  fprintf(stderr,"Removed %d angles with dummy atoms, now %d angles\n",
	  ps->nr-j,j);
  ps->nr=j;
}

void clean_dum_dihs(t_params *ps, int natom, char dihname[], t_params *plist, 
		    bool *is_dum)
{
  int      i,j,k,l,m,n,ndum,kept_i;
  atom_id  atom,constr;
  atom_id  dumatoms[3];
  bool     keep,used,present;
  struct { int i,j; } *pindex;
  
  /* make index into dummy entries of plist: */
  snew(pindex,natom);
  for(i=F_DUMMY1; (i<=F_DUMMY3); i++)
    for (j=0; (j<plist[i].nr); j++) {
      k=plist[i].param[j].AI;
      pindex[k].i=i;
      pindex[k].j=j;
    }
  
  kept_i=0;
  for(i=0; (i<ps->nr); i++) { /* for all dihedrals in the plist */
    keep=FALSE;
    /* check if all dummies are constructed from the same atoms */
    ndum=0;
    for(k=0; (k<4) && !keep; k++) { /* for all atoms in the dihedral */
      atom = ps->param[i].a[k];
      if (is_dum[atom]) {
	ndum++;
	if (ndum==1) {
	  /* store construction atoms of first dummy */
	  for(m=1; (m<4); m++)
	    dumatoms[m-1]=plist[pindex[atom].i].param[pindex[atom].j].a[m];
	  if (debug) {
	    fprintf(debug,"dih w. dum: %d %d %d %d\n",
		    ps->param[i].AI+1,ps->param[i].AJ+1,
		    ps->param[i].AK+1,ps->param[i].AL+1);
	    fprintf(debug,"dum %d from: %d %d %d\n",
		    atom+1,dumatoms[0]+1,dumatoms[1]+1,dumatoms[2]+1);
	  }
	} else 
	  /* check if this dummy is constructed from the same atoms */
	  for(m=1; (m<4) && !keep; m++) {
	    present=FALSE;
	    constr=plist[pindex[atom].i].param[pindex[atom].j].a[m];
	    for(n=0; (n<3) && !present; n++)
	      if (constr == dumatoms[n])
		present=TRUE;
	    if (!present)
	      keep=TRUE;
	  }
      }
    }
    
    /* keep all dihedrals with no dummies in them */
    if (ndum==0)
      keep=TRUE;
    
    /* check if all atoms in dihedral are either dummies, or used in 
       construction of dummies. If so, keep it, if not throw away: */
    for(k=0; (k<4) && !keep; k++) { /* for all atoms in the dihedral */
      atom = ps->param[i].a[k];
      if (!is_dum[atom]) {
	used=FALSE;
	for(m=0; (m<3) && !used; m++)
	  if (atom == dumatoms[m])
	    used=TRUE;
	if (!used) {
	  keep=TRUE;
	  if (debug) fprintf(debug,"unused atom in dih: %d\n",atom+1);
	}
      }
    }
      
    if ( keep ) {
      memcpy(&(ps->param[kept_i]),
	     &(ps->param[i]),(size_t)sizeof(ps->param[0]));
      kept_i++;
    }
  }
  
  fprintf(stderr,"Removed %d %s dihedrals with dummy atoms, "
	  "now %d %s dihedrals\n",
	  ps->nr-kept_i,dihname,kept_i,dihname);
  ps->nr=kept_i;
  
  /* clean up */
  sfree(pindex);
}

static int get_bonds(int atom, int excl, 
		     int na, int **a2, t_params *psb, bool is_dum[])
{
  int i,newa,*newa2;
  
  /* find bonds belonging to this atom 
   * which are not a dummy and are not atom excl 
   */
  newa2=*a2;
  for(i=0; (i<psb->nr); i++) {
    /* check if other atom is not a dummy, it is the one we want */
    newa=NOTSET;
    if ( (psb->param[i].AI==atom) && !is_dum[psb->param[i].AJ] )
      newa=psb->param[i].AJ;
    else if ( (psb->param[i].AJ==atom) && !is_dum[psb->param[i].AI] )
      newa=psb->param[i].AI;
    if ((newa != NOTSET) && (newa!=excl)) {
      na++;
      srenew(newa2,na);
      newa2[na-1]=newa;
    }
  }
  *a2=newa2;
  return na;
}

static int check_set(int atom, int ak[4])
{
  int k;
  
  for (k=0; (k<4) && (atom!=NOTSET); k++)
    if (atom==ak[k])
      atom=NOTSET;
      
  return atom;
}
  

static int get_atom(char *atnm, int incl,
		    int nms, int *ms, int na2, int *a2, int *na3, int **a3, 
		    int ak[4],char **aname[])
{
  int  i,j,k,atom,atnmlen;
  bool bWild;
  
  atnmlen=strlen(atnm);
  bWild=(strcmp(atnm,"*")==0);
  atom=NOTSET;
  
  /* check atom 'incl' */
  if ( (incl != NOTSET) && 
       (bWild || strncasecmp(atnm,*(aname[incl]),atnmlen)==0) )
    atom=check_set(incl, ak);
  
  /* search in masses */
  if (!bWild)
    for (i=0; (i<nms) && (atom==NOTSET); i++)
      if (strncasecmp(atnm,*(aname[ms[i]]),atnmlen)==0)
	atom=check_set(ms[i], ak);
  
  /* search in 1st bonded */
  for (i=0; (i<na2) && (atom==NOTSET); i++)
    if (bWild || (strncasecmp(atnm,*(aname[a2[i]]),atnmlen)==0) )
      atom=check_set(a2[i], ak);
  
  /* search in 2nd bonded */
  for (i=0; (i<na2) && (atom==NOTSET); i++)
    for (j=0; (j<na3[i]) && (atom==NOTSET); j++)
      if (bWild || (strncasecmp(atnm,*(aname[a3[i][j]]),atnmlen)==0) )
	atom=check_set(a3[i][j], ak);
  
  return atom;
}

void do_dum_top(t_params *psb, t_params *psd2, t_params *psd3, 
		t_params *psd2FD, t_params *psd2FAD, t_params *psda,
		int nddb,t_dumblock *ddb,bool is_dum[], t_atoms *at, 
		t_atomtype *atype, int nrtp, t_restp rtp[], real mHmult)
{
  int     i,j,k,l,m,n,resnr;
  int     ai,aj,ak[4],*ad;
  int     nmba,naba,nms,*ms,na2,*a2,*na3,**a3;
  ushort  type;
  bool    *bProcessed,bWild,bLookH,bNterm;
  t_restp *rp;
  char    *bres,bname[STRLEN],tpstr[STRLEN];
  t_dmbp  *c;

  bWild=FALSE;
  nmba=0;
  naba=0;
  snew(bProcessed,at->nr);
  /* add bonds to dummy masses */
  
  /* loop over all entries in the .ddb database (atomtypes) */
  for (n=0; (n<nddb); n++) {
    if (debug) fprintf(debug,"processing block %s\n",ddb[n].bname);
    if (bWild)
      fatal_error(0,"Found block in .ddb (%s) after wildcard block",
		  ddb[n].bname);
    strcpy(bname,ddb[n].bname);
    bWild=(strcmp(bname,"*")==0);
    bres=NULL;
    if (bres=strchr(bname,':')) {
      bres[0]='\0';
      bres++;
      if (debug) fprintf(debug,"atomtype %s residue name %s\n",bname,bres);
    }
    resnr=-1;
    /* loop over all atoms */
    for(i=0; (i<at->nr); i++) {
      /* look up this residue in rtp */
      if (at->atom[i].resnr != resnr) {
	resnr=at->atom[i].resnr;
	rp=search_rtp(*(at->resname[resnr]),nrtp,rtp);
	bNterm=is_protein(*(at->resname[resnr])) && (resnr == 0);
      }
      /* look up type for this atom */
      if (at->atom[i].m == 0) {
	j=search_jtype(rp,*(at->atomname[i]),bNterm);
	type = rp->atom[j].type;
      } else
	type = at->atom[i].type;
      /* check if it has not yet been processed 
       * and if the type matches current .ddb block (bname) */
      if ( !bProcessed[i] &&
	   ( bWild || (strcmp(type2nm(type,atype),bname)==0) ) && 
	   ( !bres || (strcmp(*(at->resname[resnr]),bres)==0) ) ) {
	if (debug) fprintf(stderr,"Found atom (%s%d %d %s)\n",
			   *(at->resname[at->atom[i].resnr]),
			   at->atom[i].resnr+1,i+1,*(at->atomname[i]));
	bProcessed[i]=TRUE;/* mark this atom */
	
	/* find masses (if any) added to atom i */
	nms=0;
	for (m=0; (m<ddb[n].nmass); m++)
	  nms+=ddb[n].mass[m].nm;
	snew(ms,nms);
	for (j=0; (j<nms); j++) {
	  ms[j]=i-nms+j;
	  bProcessed[ms[j]]=TRUE;/* mark dummy mass */
	}
	if (debug && nms) {
	  fprintf(debug,"Masses from %d (%s):",i+1,*(at->atomname[i]));
	  for (m=0; (m<nms); m++)
	    fprintf(debug," %d (%s)",ms[m]+1,*(at->atomname[ms[m]]));
	  fprintf(debug,"\n");
	}
	
	/* get bonds to non-dummies from atom i */
	na2=0;
	a2=NULL;
	na2=get_bonds(i, NOTSET, na2, &a2, psb, is_dum);
	if (na2==0)
	  fprintf(stderr,
		  "WARNING: atom (%s%d %d %s) has no bonds or "
		  "only bonds to dummy atoms\n"
		  "at least one bond to non-dummy atom expected\n",
		  *(at->resname[at->atom[i].resnr]),at->atom[i].resnr+1,
		  i+1,*(at->atomname[i]));
	
	/* get bonds from these atoms, excluding i */
	snew(a3,na2);
	snew(na3,na2);
	for(j=0; (j<na2); j++) {
	  na3[j]=get_bonds(a2[j], i, na3[j], &(a3[j]), psb, is_dum);
	  if (na3[j]==1)
	    na3[j]=get_bonds(a3[j][0], a2[j], na3[j], &(a3[j]), psb, is_dum);
	}
	
	if (debug) {
	  fprintf(debug,"Bound to %s%d %d %s (%s):",
		  *(at->resname[at->atom[i].resnr]),at->atom[i].resnr+1,
		  i+1,*(at->atomname[i]),ddb[n].bname);
	  for(j=0; (j<na2); j++) {
	    fprintf(debug, " %d (", a2[j]+1);
	    for (k=0; (k<na3[j]); k++) 
	      fprintf(debug, "%d ", a3[j][k]+1);
	    fprintf(debug, ")");
	  }
	  fprintf(debug, "\n");
	}
	
	j=i;
	/* find dummy atoms */
	snew(ad,ddb[n].ndum);
	for (m=0; (m<ddb[n].ndum); m++) {
	  ad[m]=NOTSET;
	  strcpy(tpstr,ddb[n].dum[m].na[0]);
	  if (!bWild || 
	      strncasecmp(tpstr,*(at->atomname[i]),strlen(tpstr))==0) {
	    bLookH=is_hydrogen(ddb[n].dum[m].na[0]);
	    while ( (j<at->nr) && (ad[m]==NOTSET) && 
		    (resnr == at->atom[j].resnr) &&
		    (!bLookH || (j==i) || is_hydrogen(*(at->atomname[j]))) ) {
	      if (strncasecmp(tpstr,*(at->atomname[j]),strlen(tpstr))==0) {
		bProcessed[j]=TRUE;/* mark dummy atom */
		ad[m]=j;
	      }
	      j++;
	    }
	    if (ad[m]==NOTSET)
	      fprintf(stderr,
		      "Skipped dummy atom %s on atom type %s in residue %s%d "
		      "(probably ok)\n",
		      tpstr,ddb[n].bname,*(at->resname[resnr]),resnr+1);
	  }
	}
	if (debug) {
	  fprintf(debug,"Dummies from %d (%s):",i+1,*(at->atomname[i]));
	  for (m=0; (m<ddb[n].ndum); m++)
	    if (ad[m]==NOTSET)
	      fprintf(debug," NOTSET");
	    else
	      fprintf(debug," %d (%s)",ad[m]+1,*(at->atomname[ad[m]]));
	  fprintf(debug,"\n");
	}
	
	/* add bonds to dummy masses and between them */
	if (!bWild) {
	  for (m=0; (m<ddb[n].nmass); m++) {
	    real mtot;
	    
	    /* first calc masses */
	    mtot=0;
	    for (k=0; (k < ddb[n].ndum); k++)
	      if (ad[k]!=NOTSET)
		mtot+=at->atom[ad[k]].m;
	    
	    /* find bondlengths */
	    snew(c,ddb[n].mass[m].nm);
	    for (k=0; (k < ddb[n].mass[m].nm); k++) {
	      c[k][1]=500000; /* bond force constant (fictitious: SHAKE) */
	      for (l=2; (l<MAXFORCEPARAM); l++)
		c[k][l]=NOTSET;
	    }
	    switch(ddb[n].mass[m].tp) {
	    case 1:{
	      real mHtot,dpar,dperp,dCN,dNH,aCNH,dCM,dMM,x,y;
	      
	      mHtot=3*ddb[n].mass[m].c[3];/* mass of H (calc'ed in pdb2gmx) */
	      
	      /* get params */
	      dCN = ddb[n].mass[m].c[0];
	      dNH = ddb[n].mass[m].c[1];
	      aCNH= ddb[n].mass[m].c[2];
	      
	      dpar  = dCN - (mHtot/mtot)*dNH*cos(DEG2RAD*aCNH);
	      dperp =   sqrt(mHtot/mtot)*dNH*sin(DEG2RAD*aCNH);
	      dCM = sqrt(sqr(dpar)+sqr(dperp));
	      dMM = 2 * dperp;
	      
	      /* loop over bonds goes: C-M, M-M, M-C */
	      c[0][0]=dCM; /* bondlength */
	      c[1][0]=dMM; /* bondlength */
	      
	      /* now generate parameters for dummy atoms */
	      ddb[n].dum[0].a=ddb[n].dum[0].b=dCN/(2*dpar);
	      x=(dCN - dNH*cos(DEG2RAD*aCNH))/(2*dpar);
	      y=0.5*sqrt(mtot/mHtot);
	      ddb[n].dum[3].a=x+y;
	      ddb[n].dum[3].b=x-y;
	      /* keep x from previous calculation */
	      y=(dNH*sin(DEG2RAD*aCNH)*sin(DEG2RAD*30))/(2*dperp);
	      ddb[n].dum[1].a=ddb[n].dum[2].a=x-y;
	      ddb[n].dum[1].b=ddb[n].dum[2].b=x+y;
	      ddb[n].dum[1].c=( (dNH*sin(DEG2RAD*aCNH)*cos(DEG2RAD*30))/
				(2*dpar*dperp) );
	      ddb[n].dum[2].c=-ddb[n].dum[1].c;
	    }
	    break;
	    default:
	      fatal_error(0,"Invalid dummy mass type (%d) in do_dum_top\n",
			  ddb[n].mass[m].tp);
	    }/* end switch */
	    
	    /* make bonds */
	    ai=a2[0]; /* start with first bound atom */
	    for (k=0; (k < ddb[n].mass[m].nm); k++) {
	      at->atom[ms[k]].m=at->atom[ms[k]].mB=mtot/ddb[n].mass[m].nm;
	      aj=ai; /* we go around the loop */
	      ai=ms[k];  /* next dummy mass */
	      add_param(psb,ai,aj,c[k]);
	      nmba++;
	    }
	    /* close loop of bonds */
	    aj=ai;
	    ai=a2[0];
	    add_param(psb,ai,aj,c[0]);
	    nmba++;
	    sfree(c);
	  }
	}
	
	/* add dummy parameters */
	for (m=0; (m<ddb[n].ndum); m++) {
	  /* check if first atom found */
	  if (ad[m]!=NOTSET) {
	    /* init ak */
	    ak[0]=ad[m];
	    for (k=1; (k<4); k++)
	      ak[k]=NOTSET;
	    /* find next atoms along bonds */
	    for (k=1; (k<4); k++) {
	      ak[k]=get_atom(ddb[n].dum[m].na[k],is_dum[i]?NOTSET:i,
			     nms,ms,na2,a2,na3,a3,
			     ak,at->atomname);
	      if (ak[k]==NOTSET)
		fatal_error(0,"Atom %d (%s) not found along bonds "
			    "while adding dummy atoms (%s)\n",
			    k,ddb[n].dum[m].na[k],ddb[n].bname);
	    }
	    
	    /* if no dummy masses have been added 
	       add mass of dummy to atom i */
	    if (ddb[n].nmass == 0) {
	      if (at->atom[i].m == 0.0)
		fprintf(stderr,"WARNING: adding mass of dummy atom (%d) to "
			"massless atom (%d)\n",ak[0],i);
	      if (ak[0]==i)
		k=ak[1];
	      else
		k=i;
	      if (debug) fprintf(debug,"adding %g to atom %d %g\n",
				 at->atom[ak[0]].m,k+1,at->atom[k].m); 
	      at->atom[k].m +=at->atom[ak[0]].m ;
	      at->atom[k].mB+=at->atom[ak[0]].mB;
	    }
	    
	    if (debug) {
	      fprintf(stderr, "Adding dummy: type: %d atom: %d "
		      "from: %d %d %d param: %g %g",ddb[n].dum[m].tp,
		      ak[0]+1,ak[1]+1,ak[2]+1,ak[3]+1,
		      ddb[n].dum[m].a,ddb[n].dum[m].b);
	      if (ddb[n].dum[m].c != NOTSET)
		fprintf(stderr," %g\n",ddb[n].dum[m].c);
	      else
		fprintf(stderr,"\n");
	    }
	    
	    for (k=0; (k<4); k++)
	      if (ak[k]==NOTSET)
		fatal_error(0,"Failed to find all atoms (atom %d missing) "
			    "in do_dum_top (%s)",k,ddb[n].bname);
	    
	    switch(ddb[n].dum[m].tp) {
	    case 2: /* type 2 */
	      add_dum2_param(psd2,ak[0],ak[1],ak[2],ak[3],
			     ddb[n].dum[m].a,ddb[n].dum[m].b);
	      break;
	    case 3: /* type 3 */
	      add_dum3_param(psd3,ak[0],ak[1],ak[2],ak[3],ddb[n].dum[m].a,
			     ddb[n].dum[m].b,ddb[n].dum[m].c);
	      break;
	    case 4: /* type 2' */
	      add_dum2_param(psd2FD,ak[0],ak[1],ak[2],ak[3],
			     ddb[n].dum[m].a,ddb[n].dum[m].b);
	      break;
	    case 5: /* type 2'' */
	      add_dum2_param(psd2FAD,ak[0],ak[1],ak[2],ak[3],
			     ddb[n].dum[m].a,ddb[n].dum[m].b);
	      break;
	    default:
	      fatal_error(0,"Dummy type %d not supported in .ddb\n",
			  ddb[n].dum[m].tp);
	    }
	  }
	}
	
	/* add angle constraint parameters */
	for (m=0; (m<ddb[n].nang); m++)
	  if (!bWild || 
	      strncasecmp(/* this is so emacs doesn't loose track */
			  ddb[n].ang[m].na[0],*(at->atomname[i]),
			  strlen(ddb[n].ang[m].na[0]))==0) {
	    for (k=0; (k<3); k++)
	      ak[k]=NOTSET;
	    /* find first atom (will be angle-constrained) */
	    while ( (j<at->nr) && (ak[0]==NOTSET) ) {
	      if (strncasecmp(
			      ddb[n].ang[m].na[0],*(at->atomname[j]),
			      strlen(ddb[n].ang[m].na[0]))==0) {
		bProcessed[j]=TRUE;/* mark constrained atom */
		ak[0]=j;
	      }
	      j++;
	    }
	    if (ak[0]==NOTSET)
	      fatal_error(0,"Atom %s not found "
			  "while adding angle constraints (%s)\n",
			  ddb[n].ang[m].na[0],ddb[n].bname);
	    /* find next atoms along bonds */
	    for (k=1; (k<3); k++) {
	      ak[k]=get_atom(ddb[n].ang[m].na[k],is_dum[i]?NOTSET:i,
			     nms,ms,na2,a2,na3,a3,ak,at->atomname);
	      if (ak[k]==NOTSET)
		fatal_error(0,"Atom %s not found along bonds "
			    "while adding angle constraints (%s)\n",
			    ddb[n].ang[m].na[k],ddb[n].bname);
	    }
	    
	    for (k=0; (k<4); k++)
	      if (ak[k]==NOTSET)
		fatal_error(0,"DEATH HORROR ERROR in do_dum_top");
	    {
	      real cAngCon[MAXFORCEPARAM];
	      add_param(psb,ak[0],ak[2],ddb[n].ang[m].c);
	      naba++;
	    }
	    
	  }
	sfree(ms);
	sfree(a2);
	for(j=0; (j<na2); j++)
	  sfree(a3[j]);
	sfree(a3);
	sfree(na3);
      }
    }
  }
  
  if (naba>0)
    fprintf(stderr,"Added %d bonds to constrain angles, now %d bonds\n",
	    naba,psb->nr-nmba);
  if (nmba>0)
    fprintf(stderr,"Added %d bonds to dummy masses, now %d bonds\n",
	    nmba,psb->nr);
  
  /* remove bonds with dummy atoms */
  for(i=j=0; (i<psb->nr); i++)
    if (!is_dum[psb->param[i].AI] && !is_dum[psb->param[i].AJ]) {
      memcpy(&(psb->param[j]),&(psb->param[i]),(size_t)sizeof(psb->param[0]));
      j++;
    }
  i=psb->nr-j;
  psb->nr=j;
  
  if (i>0)
    fprintf(stderr,"Removed %d bonds to dummy atoms, now %d bonds\n",i,j);
  
  /* set mass of dummy atoms to 0 */
  for(i=0; (i<at->nr); i++)
    if (is_dum[i])
      at->atom[i].m=at->atom[i].mB=0;
  
  /* cleaning up */
  sfree(bProcessed);
  if (debug) fflush(debug);
}

void do_dum_excl(t_block *excl, bool is_dum[])
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
    if (is_dum[i]) {
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

