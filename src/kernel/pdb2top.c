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

#include <stdio.h>
#include <math.h>
#include <ctype.h>
#include "vec.h"
#include "copyrite.h"
#include "smalloc.h"
#include "macros.h"
#include "symtab.h"
#include "futil.h"
#include "gmx_fatal.h"
#include "pdb2top.h"
#include "gpp_nextnb.h"
#include "topdirs.h"
#include "toputil.h"
#include "h_db.h"
#include "pgutil.h"
#include "resall.h"
#include "topio.h"
#include "string2.h"
#include "physics.h"
#include "pdbio.h"
#include "gen_ad.h"
#include "filenm.h"
#include "index.h"
#include "gen_vsite.h"
#include "add_par.h"
#include "toputil.h"

/* this must correspond to enum in pdb2top.h */
const char *hh[ehisNR]   = { "HISA", "HISB", "HISH", "HIS1" };

static int missing_atoms(t_restp *rp, int resind,
			 t_atoms *at, int i0, int i, bool bCTer)
{
  int  j,k,nmiss;
  char *name;
  bool bFound, bRet;

  nmiss = 0;
  for (j=0; j<rp->natom; j++) {
    name=*(rp->atomname[j]);
    /* if ((name[0]!='H') && (name[0]!='h') && (!bCTer || (name[0]!='O'))) { */
    bFound=FALSE;
    for (k=i0; k<i; k++) 
      bFound=(bFound || !strcasecmp(*(at->atomname[k]),name));
    if (!bFound) {
      nmiss++;
      fprintf(stderr,"\nWARNING: "
	      "atom %s is missing in residue %s %d in the pdb file\n",
	      name,*(at->resinfo[resind].name),at->resinfo[resind].nr);
      if (name[0]=='H' || name[0]=='h')
	fprintf(stderr,"         You might need to add atom %s to the hydrogen database of residue %s\n"
		       "         in the file ff???.hdb (see the manual)\n",
		name,*(at->resinfo[resind].name));
      fprintf(stderr,"\n");
      }
      /* } */
  }
  
  return nmiss;
}

bool is_int(double x)
{
  const double tol = 1e-4;
  int   ix;
  
  if (x < 0)
    x=-x;
  ix=gmx_nint(x);
  
  return (fabs(x-ix) < tol);
}

typedef struct { char *desc,*fn; } t_fff;

void
choose_ff(char *forcefield, int maxlen)
{
  FILE    *in;
  t_fff   *fff;
  int     i,nff,sel;
  char    *c,buf[STRLEN],fn[32];
  char    *pret;

  in=libopen("FF.dat");
  fgets2(buf,255,in);
  sscanf(buf,"%d",&nff);
  snew(fff,nff);
  for(i=0; (i<nff); i++) {
    fgets2(buf,255,in);
    sscanf(buf,"%s",fn);
    fff[i].fn=strdup(fn);
    /* Search for next non-space character, there starts description */
    c=&(buf[strlen(fn)+1]);
    while (isspace(*c)) c++;
    fff[i].desc=strdup(c);
  }
  fclose(in);

  if (nff > 1) {
    printf("\nSelect the Force Field:\n");
    for(i=0; (i<nff); i++)
      printf("%2d: %s\n",i,fff[i].desc);
    do {
      pret = fgets(buf,STRLEN,stdin);

      if(pret != NULL)
	sscanf(buf,"%d",&sel);
    } while ( pret==NULL || (sel < 0) || (sel >= nff));
  }
  else
    sel=0;

  strncpy(forcefield,fff[sel].fn,maxlen);

  for(i=0; (i<nff); i++) {
    sfree(fff[i].desc);
    sfree(fff[i].fn);
  }
  sfree(fff);
  
}


static int name2type(t_atoms *at, int **cgnr, gpp_atomtype_t atype, 
		     t_restp restp[])
{
  int     i,j,prevresind,resind,i0,prevcg,cg,curcg;
  char    *name;
  bool    bProt, bNterm;
  double  qt;
  int     nmissat;
  t_aa_names *aan;
  
  nmissat = 0;

  resind=-1;
  bProt=FALSE;
  bNterm=FALSE;
  i0=0;
  snew(*cgnr,at->nr);
  qt=0;
  prevcg=NOTSET;
  curcg=0;
  cg=-1;
  j=NOTSET;
  aan = get_aa_names();
  
  for(i=0; (i<at->nr); i++) {
    prevresind=resind;
    if (at->atom[i].resind != resind) {
      resind = at->atom[i].resind;
      bProt = is_protein(aan,*(at->resinfo[resind].name));
      bNterm=bProt && (resind == 0);
      if (resind > 0)
	nmissat += 
	  missing_atoms(&restp[prevresind],prevresind,at,i0,i,
			(!bProt && is_protein(aan,restp[prevresind].resname)));
      i0=i;
    }
    if (at->atom[i].m == 0) {
      if (debug)
	fprintf(debug,"atom %d%s: curcg=%d, prevcg=%d, cg=%d\n",
		i+1,*(at->atomname[i]),curcg,prevcg,
		j==NOTSET ? NOTSET : restp[resind].cgnr[j]);
      qt=0;
      prevcg=cg;
      name=*(at->atomname[i]);
      j=search_jtype(&restp[resind],name,bNterm);
      at->atom[i].type = restp[resind].atom[j].type;
      at->atom[i].q    = restp[resind].atom[j].q;
      at->atom[i].m    = get_atomtype_massA(restp[resind].atom[j].type,
					    atype);
      cg = restp[resind].cgnr[j];
      if ( (cg != prevcg) || (resind != prevresind) )
	curcg++;
    } else {
      if (debug)
	fprintf(debug,"atom %d%s: curcg=%d, qt=%g, is_int=%d\n",
		i+1,*(at->atomname[i]),curcg,qt,is_int(qt));
      cg=-1;
      if (is_int(qt)) {
	qt=0;
	curcg++;
      }
      qt+=at->atom[i].q;
    }
    (*cgnr)[i]=curcg;
    at->atom[i].typeB = at->atom[i].type;
    at->atom[i].qB    = at->atom[i].q;
    at->atom[i].mB    = at->atom[i].m;
  }
  nmissat += missing_atoms(&restp[resind],resind,at,i0,i,
			   (!bProt || is_protein(aan,restp[resind].resname)));
  done_aa_names(&aan);
			   
  return nmissat;
}

static void print_top_heavy_H(FILE *out, real mHmult)
{
  if (mHmult == 2.0) 
    fprintf(out,"; Using deuterium instead of hydrogen\n\n");
  else if (mHmult == 4.0)
    fprintf(out,"#define HEAVY_H\n\n");
  else if (mHmult != 1.0)
    fprintf(stderr,"WARNING: unsupported proton mass multiplier (%g) "
	    "in pdb2top\n",mHmult);
}

void print_top_comment(FILE *out,const char *filename,const char *title,bool bITP)
{
  char tmp[256]; 
  
  nice_header(out,filename);
  fprintf(out,";\tThis is your %stopology file\n",bITP ? "include " : "");
  cool_quote(tmp,255,NULL);
  fprintf(out,";\t%s\n",title[0]?title:tmp);
  fprintf(out,";\n");
}

void print_top_header(FILE *out,const char *filename, 
		      const char *title,bool bITP,const char *ff,real mHmult)
{
  print_top_comment(out,filename,title,bITP);

  print_top_heavy_H(out, mHmult);
  fprintf(out,"; Include forcefield parameters\n");
  fprintf(out,"#include \"%s.itp\"\n\n",ff);
}

static void print_top_posre(FILE *out,const char *pr)
{
  fprintf(out,"; Include Position restraint file\n");
  fprintf(out,"#ifdef POSRES\n");
  fprintf(out,"#include \"%s\"\n",pr);
  fprintf(out,"#endif\n\n");
}
  
static void print_top_water(FILE *out,const char *water)
{
  fprintf(out,"; Include water topology\n");
  fprintf(out,"#include \"%s.itp\"\n",water);
  fprintf(out,"\n");
  fprintf(out,"#ifdef POSRES_WATER\n");
  fprintf(out,"; Position restraint for each water oxygen\n");
  fprintf(out,"[ position_restraints ]\n");
  fprintf(out,";%3s %5s %9s %10s %10s\n","i","funct","fcx","fcy","fcz");
  fprintf(out,"%4d %4d %10g %10g %10g\n",1,1,1000.0,1000.0,1000.0);
  fprintf(out,"#endif\n\n");
  fprintf(out,"; Include generic topology for ions\n");
  fprintf(out,"#include \"ions.itp\"\n");
  fprintf(out,"\n");
}

static void print_top_system(FILE *out, const char *title)
{
  fprintf(out,"[ %s ]\n",dir2str(d_system));
  fprintf(out,"; Name\n");
  fprintf(out,"%s\n\n",title[0]?title:"Protein");
}

void print_top_mols(FILE *out, const char *title, const char *water,
		    int nincl, char **incls, int nmol, t_mols *mols)
{
  int i;
  
  if (nincl>0) {
    fprintf(out,"; Include chain topologies\n");
    for (i=0; (i<nincl); i++)
      fprintf(out,"#include \"%s\"\n",incls[i]);
    fprintf(out,"\n");
  }

  if (water)
    print_top_water(out,water);
  print_top_system(out, title);
  
  if (nmol) {
    fprintf(out,"[ %s ]\n",dir2str(d_molecules));
    fprintf(out,"; %-15s %5s\n","Compound","#mols");
    for (i=0; (i<nmol); i++)
      fprintf(out,"%-15s %5d\n",mols[i].name,mols[i].nr);
  }
}

void write_top(FILE *out, char *pr,char *molname,
	       t_atoms *at,int bts[],t_params plist[],t_excls excls[],
	       gpp_atomtype_t atype,int *cgnr, int nrexcl)
     /* NOTE: nrexcl is not the size of *excl! */
{
  if (at && atype && cgnr) {
    fprintf(out,"[ %s ]\n",dir2str(d_moleculetype));
    fprintf(out,"; %-15s %5s\n","Name","nrexcl");
    fprintf(out,"%-15s %5d\n\n",molname?molname:"Protein",nrexcl);
    
    print_atoms(out, atype, at, cgnr);
    print_bondeds(out,at->nr,d_bonds,      F_BONDS,    bts[ebtsBONDS], plist);
    print_bondeds(out,at->nr,d_constraints,F_CONSTR,   0,              plist);
    print_bondeds(out,at->nr,d_constraints,F_CONSTRNC, 0,              plist);
    print_bondeds(out,at->nr,d_pairs,      F_LJ14,     0,              plist);
    print_excl(out,at->nr,excls);
    print_bondeds(out,at->nr,d_angles,     F_ANGLES,   bts[ebtsANGLES],plist);
    print_bondeds(out,at->nr,d_dihedrals,  F_PDIHS,    bts[ebtsPDIHS], plist);
    print_bondeds(out,at->nr,d_dihedrals,  F_IDIHS,    bts[ebtsIDIHS], plist);
    print_bondeds(out,at->nr,d_cmap,       F_CMAP,     bts[ebtsCMAP],  plist);
    print_bondeds(out,at->nr,d_polarization,F_POLARIZATION,   0,       plist);
    print_bondeds(out,at->nr,d_thole_polarization,F_THOLE_POL,0,       plist);
    print_bondeds(out,at->nr,d_vsites2,    F_VSITE2,   0,              plist);
    print_bondeds(out,at->nr,d_vsites3,    F_VSITE3,   0,              plist);
    print_bondeds(out,at->nr,d_vsites3,    F_VSITE3FD, 0,              plist);
    print_bondeds(out,at->nr,d_vsites3,    F_VSITE3FAD,0,              plist);
    print_bondeds(out,at->nr,d_vsites3,    F_VSITE3OUT,0,              plist);
    print_bondeds(out,at->nr,d_vsites4,    F_VSITE4FD, 0,              plist);
    print_bondeds(out,at->nr,d_vsites4,    F_VSITE4FDN, 0,             plist);
    
    if (pr)
      print_top_posre(out,pr);
  }
}

static atom_id search_res_atom(const char *type,int resind,
			       int natom,t_atom at[],
			       const char * const * const aname[],
			       const char *bondtype,bool bMissing)
{
  int i;

  for(i=0; (i<natom); i++)
    if (at[i].resind == resind)
      return search_atom(type,i,natom,at,aname,bondtype,bMissing);
  
  return NO_ATID;
}

static void do_ssbonds(t_params *ps,int natoms,t_atom atom[],char **aname[],
		       int nssbonds,t_ssbond *ssbonds,bool bMissing)
{
  int     i,ri,rj;
  atom_id ai,aj;
  
  for(i=0; (i<nssbonds); i++) {
    ri = ssbonds[i].res1;
    rj = ssbonds[i].res2;
    ai = search_res_atom(ssbonds[i].a1,ri,natoms,atom,aname,
			 "special bond",bMissing);
    aj = search_res_atom(ssbonds[i].a2,rj,natoms,atom,aname,
			 "special bond",bMissing);
    if ((ai == NO_ATID) || (aj == NO_ATID))
      gmx_fatal(FARGS,"Trying to make impossible special bond (%s-%s)!",
		  ssbonds[i].a1,ssbonds[i].a2);
    add_param(ps,ai,aj,NULL,NULL);
  }
}

static void at2bonds(t_params *psb, t_hackblock *hb,
		     int natoms, t_atom atom[], char **aname[], 
		     int nres, rvec x[], 
		     real long_bond_dist, real short_bond_dist)
{
  int     resind,i,j,k;
  atom_id ai,aj;
  real    dist2, long_bond_dist2, short_bond_dist2;
  const char *ptr;

  long_bond_dist2  = sqr(long_bond_dist);
  short_bond_dist2 = sqr(short_bond_dist);

  if (debug)
    ptr = "bond";
  else
    ptr = "check";

  fprintf(stderr,"Making bonds...\n");
  i=0;
  for(resind=0; (resind < nres) && (i<natoms); resind++) {
    /* add bonds from list of bonded interactions */
    for(j=0; j < hb[resind].rb[ebtsBONDS].nb; j++) {
      /* Unfortunately we can not issue errors or warnings
       * for missing atoms in bonds, as the hydrogens and terminal atoms
       * have not been added yet.
       */
      ai=search_atom(hb[resind].rb[ebtsBONDS].b[j].AI,i,natoms,atom,aname,
		     ptr,TRUE);
      aj=search_atom(hb[resind].rb[ebtsBONDS].b[j].AJ,i,natoms,atom,aname,
		     ptr,TRUE);
      if (ai != NO_ATID && aj != NO_ATID) {
	dist2 = distance2(x[ai],x[aj]);
	if (dist2 > long_bond_dist2 )
	  fprintf(stderr,"Warning: Long Bond (%d-%d = %g nm)\n",
		  ai+1,aj+1,sqrt(dist2));
	else if (dist2 < short_bond_dist2 )
	  fprintf(stderr,"Warning: Short Bond (%d-%d = %g nm)\n",
		  ai+1,aj+1,sqrt(dist2));
	add_param(psb,ai,aj,NULL,hb[resind].rb[ebtsBONDS].b[j].s);
      }
    }
    /* add bonds from list of hacks (each added atom gets a bond) */
    while( (i<natoms) && (atom[i].resind == resind) ) {
      for(j=0; j < hb[resind].nhack; j++)
	if ( ( hb[resind].hack[j].tp > 0 ||
	       hb[resind].hack[j].oname==NULL ) &&
	     strcmp(hb[resind].hack[j].AI,*(aname[i])) == 0 ) {
	  switch(hb[resind].hack[j].tp) {
	  case 9:          /* COOH terminus */
	    add_param(psb,i,i+1,NULL,NULL);     /* C-O  */
	    add_param(psb,i,i+2,NULL,NULL);     /* C-OA */
	    add_param(psb,i+2,i+3,NULL,NULL);   /* OA-H */
	    break;
	  default:
	    for(k=0; (k<hb[resind].hack[j].nr); k++)
	      add_param(psb,i,i+k+1,NULL,NULL);
	  }
	}
      i++;
    }
    /* we're now at the start of the next residue */
  }
}

static int pcompar(const void *a, const void *b)
{
  t_param *pa,*pb;
  int     d;
  pa=(t_param *)a;
  pb=(t_param *)b;
  
  d = pa->AI - pb->AI;
  if (d == 0) 
    d = pa->AJ - pb->AJ;
  if (d == 0)
    return strlen(pb->s) - strlen(pa->s);
  else
    return d;
}

static void clean_bonds(t_params *ps)
{
  int     i,j;
  atom_id a;
  
  if (ps->nr > 0) {
    /* swap atomnumbers in bond if first larger than second: */
    for(i=0; (i<ps->nr); i++)
      if ( ps->param[i].AJ < ps->param[i].AI ) {
	a = ps->param[i].AI;
	ps->param[i].AI = ps->param[i].AJ;
	ps->param[i].AJ = a;
      }
    
    /* Sort bonds */
    qsort(ps->param,ps->nr,(size_t)sizeof(ps->param[0]),pcompar);
    
    /* remove doubles, keep the first one always. */
    j = 1;
    for(i=1; (i<ps->nr); i++) {
      if ((ps->param[i].AI != ps->param[j-1].AI) ||
	  (ps->param[i].AJ != ps->param[j-1].AJ) ) {
	cp_param(&(ps->param[j]),&(ps->param[i]));
	j++;
      } 
    }
    fprintf(stderr,"Number of bonds was %d, now %d\n",ps->nr,j);
    ps->nr=j;
  }
  else
    fprintf(stderr,"No bonds\n");
}

void print_sums(t_atoms *atoms, bool bSystem)
{
  double m,qtot;
  int    i;
  const char *where;
  
  if (bSystem)
    where=" in system";
  else
    where="";
  
  m=0;
  qtot=0;
  for(i=0; (i<atoms->nr); i++) {
    m+=atoms->atom[i].m;
    qtot+=atoms->atom[i].q;
  }
  
  fprintf(stderr,"Total mass%s %.3f a.m.u.\n",where,m);
  fprintf(stderr,"Total charge%s %.3f e\n",where,qtot);
}

static void get_hackblocks_rtp(t_hackblock **hb, t_restp **restp, 
			       int nrtp, t_restp rtp[],
			       int nres, t_resinfo *resinfo, 
			       int nterpairs,
			       t_hackblock **ntdb, t_hackblock **ctdb,
			       int *rn, int *rc)
{
  int i, j, k, l;
  t_restp *res;
  char buf[STRLEN];
  const char *Hnum="123456";
  bool bN,bC;

  snew(*hb,nres);
  snew(*restp,nres);
  /* first the termini */
  for(i=0; i<nterpairs; i++) {
    if (rn[i]>=0)
      copy_t_hackblock(ntdb[i], &(*hb)[rn[i]]);
    if (rc[i]>=0)
      merge_t_hackblock(ctdb[i], &(*hb)[rc[i]]);
  }  

  /* then the whole rtp */
  for(i=0; i < nres; i++) {
    res = search_rtp(*resinfo[i].name,nrtp,rtp);
    copy_t_restp(res, &(*restp)[i]);
    bN = FALSE;
    for(j=0; j<nterpairs && !bN; j++)
      bN = i==rn[j];
    bC = FALSE;
    for(j=0; j<nterpairs && !bC; j++)
      bC = i==rc[j];
    merge_t_bondeds(res->rb, (*hb)[i].rb,bN,bC);
  }
  
  /* now perform t_hack's on t_restp's,
     i.e. add's and deletes from termini database will be 
     added to/removed from residue topology 
     we'll do this on one big dirty loop, so it won't make easy reading! */
  for(i=0; i < nres; i++)
    for(j=0; j < (*hb)[i].nhack; j++)
      if ( (*hb)[i].hack[j].nr ) {
	/* find atom in restp */
	for(l=0; l < (*restp)[i].natom; l++)
	  if ( ( (*hb)[i].hack[j].oname==NULL && 
		 strcmp((*hb)[i].hack[j].AI, *(*restp)[i].atomname[l])==0 ) ||
	       ( (*hb)[i].hack[j].oname!=NULL &&
		 strcmp((*hb)[i].hack[j].oname,*(*restp)[i].atomname[l])==0 ) )
	    break;
	if (l == (*restp)[i].natom) 
	  gmx_fatal(FARGS,"atom %s not found in residue %d%s "
		      "while combining tdb and rtp",
		      (*hb)[i].hack[j].oname!=NULL ? 
		      (*hb)[i].hack[j].oname : (*hb)[i].hack[j].AI, 
		      i+1,*resinfo[i].name);
	if ( (*hb)[i].hack[j].oname==NULL ) {
	  /* we're adding: */
	  if (debug) 
	    fprintf(debug,"adding atom(s) %s to atom %s in res %d%s in rtp\n",
		    (*hb)[i].hack[j].nname,
		    *(*restp)[i].atomname[l], i+1, (*restp)[i].resname);
	  strcpy(buf, (*hb)[i].hack[j].nname);
	  buf[strlen(buf)+1]='\0';
	  if ( (*hb)[i].hack[j].nr > 1 )
	    buf[strlen(buf)]='-';
	  /* make space */
	  (*restp)[i].natom += (*hb)[i].hack[j].nr;
	  srenew((*restp)[i].atom,     (*restp)[i].natom);
	  srenew((*restp)[i].atomname, (*restp)[i].natom);
	  srenew((*restp)[i].cgnr,     (*restp)[i].natom);
	  /* shift the rest */
	  for(k=(*restp)[i].natom-1; k > l+(*hb)[i].hack[j].nr; k--) {
	    (*restp)[i].atom[k] =
	      (*restp)[i].atom    [k - (*hb)[i].hack[j].nr];
	    (*restp)[i].atomname[k] =
	      (*restp)[i].atomname[k - (*hb)[i].hack[j].nr];
	    (*restp)[i].cgnr[k] =
	      (*restp)[i].cgnr    [k - (*hb)[i].hack[j].nr];
	  }
	  /* now add them */
	  for(k=0; k < (*hb)[i].hack[j].nr; k++) {
	    /* set counter in atomname */
	    if ( (*hb)[i].hack[j].nr > 1 )
	      buf[strlen(buf)-1]=Hnum[k];
	    snew( (*restp)[i].atomname[l+1+k], 1);
	    (*restp)[i].atom     [l+1+k] = *(*hb)[i].hack[j].atom;
	    *(*restp)[i].atomname[l+1+k] = strdup(buf);
	    if ( (*hb)[i].hack[j].cgnr != NOTSET )
	      (*restp)[i].cgnr   [l+1+k] = (*hb)[i].hack[j].cgnr;
	    else
	      (*restp)[i].cgnr   [l+1+k] = (*restp)[i].cgnr[l];
	  }
	} else /* oname != NULL */
	  if ( (*hb)[i].hack[j].nname==NULL ) {
	    /* we're deleting */
	    if (debug) 
	      fprintf(debug, "deleting atom %s from res %d%s in rtp\n",
			       *(*restp)[i].atomname[l], 
			       i+1,(*restp)[i].resname);
	    /* shift the rest */
	    (*restp)[i].natom--;
	    for(k=l; k < (*restp)[i].natom; k++) {
	      (*restp)[i].atom    [k] = (*restp)[i].atom    [k+1];
	      (*restp)[i].atomname[k] = (*restp)[i].atomname[k+1];
	      (*restp)[i].cgnr    [k] = (*restp)[i].cgnr    [k+1];
	    }
	    /* give back space */
	    srenew((*restp)[i].atom,     (*restp)[i].natom);
	    srenew((*restp)[i].atomname, (*restp)[i].natom);
	    srenew((*restp)[i].cgnr,     (*restp)[i].natom);
	  } else { /* nname != NULL */
	    /* we're replacing */
	    if (debug) 
	      fprintf(debug, "replacing atom %s by %s in res %d%s in rtp\n",
		      *(*restp)[i].atomname[l], (*hb)[i].hack[j].nname, 
		      i+1,(*restp)[i].resname);
	    snew( (*restp)[i].atomname[l], 1);
	    (*restp)[i].atom[l]      =       *(*hb)[i].hack[j].atom;
	    *(*restp)[i].atomname[l] = strdup((*hb)[i].hack[j].nname);
	    if ( (*hb)[i].hack[j].cgnr != NOTSET )
	      (*restp)[i].cgnr   [l] = (*hb)[i].hack[j].cgnr;
	  }
      }
}

void gen_cmap(t_params *psb, t_restp *restp, int natoms, t_atom atom[], char **aname[], int nres)
{
	int     residx,i,ii,j,k;
	atom_id ai,aj,ak,al,am;
	const char *ptr;
	
	if (debug)
		ptr = "cmap";
	else
		ptr = "check";
	
	fprintf(stderr,"Making cmap torsions...");
	i=0;
	for(residx=0; residx<nres; residx++)
	{
		/* Add CMAP terms from the list of CMAP interactions */
		for(j=0;j<restp[residx].rb[ebtsCMAP].nb; j++)
		{
			ai=search_atom(restp[residx].rb[ebtsCMAP].b[j].a[0],i,natoms,atom,aname,
						   ptr,TRUE);
			aj=search_atom(restp[residx].rb[ebtsCMAP].b[j].a[1],i,natoms,atom,aname,
						   ptr,TRUE);
			ak=search_atom(restp[residx].rb[ebtsCMAP].b[j].a[2],i,natoms,atom,aname,
						   ptr,TRUE);
			al=search_atom(restp[residx].rb[ebtsCMAP].b[j].a[3],i,natoms,atom,aname,
						   ptr,TRUE);
			am=search_atom(restp[residx].rb[ebtsCMAP].b[j].a[4],i,natoms,atom,aname,
						   ptr,TRUE);
			
			/* For now, exclude the first and last residues from cmap */
			if(residx>=1 && residx<nres-1)
			{
				add_cmap_param(psb,ai,aj,ak,al,am,0,0,restp[residx].rb[ebtsCMAP].b[j].s);
			}
		}
		
		if(residx<nres-1)
		{
			while(atom[i].resind<residx+1)
			{
				i++;
			}
		}
	}
	
	/* Start the next residue */
}

static void 
scrub_charge_groups(int *cgnr, int natoms)
{
	int i;
	
	for(i=0;i<natoms;i++)
	{
		cgnr[i]=i+1;
	}
}


void pdb2top(FILE *top_file, char *posre_fn, char *molname,
	     t_atoms *atoms, rvec **x, gpp_atomtype_t atype, t_symtab *tab,
	     int bts[], int nrtp, t_restp   rtp[],
	     int nterpairs,t_hackblock **ntdb, t_hackblock **ctdb,
	     int *rn, int *rc, bool bMissing,
	     bool bH14, bool bAlldih, bool bRemoveDih,
	     bool bVsites, bool bVsiteAromatics, char *ff, real mHmult,
	     int nssbonds, t_ssbond *ssbonds, int nrexcl, 
	     real long_bond_dist, real short_bond_dist,
	     bool bDeuterate, bool bChargeGroups)
{
  t_hackblock *hb;
  t_restp  *restp;
  t_params plist[F_NRE];
  t_excls  *excls;
  t_nextnb nnb;
  int      *cgnr;
  int      *vsite_type;
  int      i,nmissat;
  
  init_plist(plist);

  /* lookup hackblocks and rtp for all residues */
  get_hackblocks_rtp(&hb, &restp, nrtp, rtp, atoms->nres, atoms->resinfo, 
		     nterpairs, ntdb, ctdb, rn, rc);
  /* ideally, now we would not need the rtp itself anymore, but do 
     everything using the hb and restp arrays. Unfortunately, that 
     requires some re-thinking of code in gen_vsite.c, which I won't 
     do now :( AF 26-7-99 */
  
  if (debug) {
    print_resall(debug, bts, atoms->nres, restp, atype,bAlldih,nrexcl,
		 bH14,bRemoveDih);
    dump_hb(debug, atoms->nres, hb);
  }
  
  /* Make bonds */
  at2bonds(&(plist[F_BONDS]), hb, 
	   atoms->nr, atoms->atom, atoms->atomname, atoms->nres, *x, 
	   long_bond_dist, short_bond_dist);
  
  /* specbonds: disulphide bonds & heme-his */
  do_ssbonds(&(plist[F_BONDS]),
	     atoms->nr, atoms->atom, atoms->atomname, nssbonds, ssbonds,
	     bMissing);
  
  nmissat = name2type(atoms, &cgnr, atype, restp);
  if (nmissat) {
    if (bMissing)
      fprintf(stderr,"There were %d missing atoms in molecule %s\n",
	      nmissat,molname);
    else
      gmx_fatal(FARGS,"There were %d missing atoms in molecule %s, if you want to use this incomplete topology anyhow, use the option -missing",
		  nmissat,molname);
  }
  
  /* Cleanup bonds (sort and rm doubles) */ 
  clean_bonds(&(plist[F_BONDS]));
  
  snew(vsite_type,atoms->nr);
  for(i=0; i<atoms->nr; i++)
    vsite_type[i]=NOTSET;
  if (bVsites)
    /* determine which atoms will be vsites and add dummy masses 
       also renumber atom numbers in plist[0..F_NRE]! */
    do_vsites(nrtp, rtp, atype, atoms, tab, x, plist, 
	       &vsite_type, &cgnr, mHmult, bVsiteAromatics, ff);
  
  /* Make Angles and Dihedrals */
  fprintf(stderr,"Generating angles, dihedrals and pairs...\n");
  snew(excls,atoms->nr);
  init_nnb(&nnb,atoms->nr,4);
  gen_nnb(&nnb,plist);
  print_nnb(&nnb,"NNB");
  gen_pad(&nnb,atoms,nrexcl,bH14,plist,excls,hb,bAlldih,bRemoveDih,bMissing);
  done_nnb(&nnb);
  
  /* Make CMAP */
  gen_cmap(&(plist[F_CMAP]), restp, atoms->nr, atoms->atom, atoms->atomname, atoms->nres);
  fprintf(stderr, "there are %4d cmap torsions\n",plist[F_CMAP].nr);
	
  /* set mass of all remaining hydrogen atoms */
  if (mHmult != 1.0)
    do_h_mass(&(plist[F_BONDS]),vsite_type,atoms,mHmult,bDeuterate);
  sfree(vsite_type);
  
  /* Cleanup bonds (sort and rm doubles) */ 
  /* clean_bonds(&(plist[F_BONDS]));*/
   
  fprintf(stderr,
	  "There are %4d dihedrals, %4d impropers, %4d angles\n"
	  "          %4d pairs,     %4d bonds and  %4d virtual sites\n",
	  plist[F_PDIHS].nr, plist[F_IDIHS].nr, plist[F_ANGLES].nr,
	  plist[F_LJ14].nr, plist[F_BONDS].nr,
	  plist[F_VSITE2].nr +
	  plist[F_VSITE3].nr +
	  plist[F_VSITE3FD].nr +
	  plist[F_VSITE3FAD].nr +
	  plist[F_VSITE3OUT].nr +
      plist[F_VSITE4FD].nr +
      plist[F_VSITE4FDN].nr );
  
  print_sums(atoms, FALSE);
  
  if (FALSE == bChargeGroups)
  {
	  scrub_charge_groups(cgnr, atoms->nr);
  }
	
  if (top_file) {
    fprintf(stderr,"Writing topology\n");
    write_top(top_file, posre_fn, molname,
	      atoms, bts, plist, excls, atype, cgnr, nrexcl);
  }
  
  /* cleaning up */
  free_t_hackblock(atoms->nres, &hb);
  free_t_restp(atoms->nres, &restp);
    
  /* we should clean up hb and restp here, but that is a *L*O*T* of work! */
  sfree(cgnr);
  for (i=0; i<F_NRE; i++)
    sfree(plist[i].param);
  for (i=0; i<atoms->nr; i++)
    sfree(excls[i].e);
  sfree(excls);
}
