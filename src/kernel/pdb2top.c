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
static char *SRCID_pdb2top_c = "$Id$";

#include <stdio.h>
#include <math.h>
#include "vec.h"
#include "copyrite.h"
#include "assert.h"
#include "smalloc.h"
#include "macros.h"
#include "symtab.h"
#include "futil.h"
#include "fatal.h"
#include "pdb2top.h"
#include "topexcl.h"
#include "topdirs.h"
#include "toputil.h"
#include "h_db.h"
#include "pgutil.h"
#include "resall.h"
#include "topio.h"
#include "physics.h"
#include "pdbio.h"
#include "gen_ad.h"
#include "filenm.h"
#include "index.h"

real distance(rvec a,rvec b) 
{
  return sqrt((real)(sqr((real)(a[XX]-b[XX])) + 
		     sqr((real)(a[YY]-b[YY])) + sqr((real)(a[ZZ]-b[ZZ]))));
}

static void add_param(t_params *ps,int ai,int aj, real *c)
{
  int i;
  
  if ((ai < 0) || (aj < 0)) 
    fatal_error(0,"Trying to add impossible atoms: ai=%d, aj=%d",ai,aj);
  srenew(ps->param,++ps->nr);
  ps->param[ps->nr-1].AI=ai;
  ps->param[ps->nr-1].AJ=aj;
  ps->param[ps->nr-1].AK=-1;
  ps->param[ps->nr-1].AL=-1;
  if (c==NULL) 
    for(i=0; (i<MAXFORCEPARAM); i++)
      ps->param[ps->nr-1].c[i]=NOTSET;
  else
     for(i=0; (i<MAXFORCEPARAM); i++)
      ps->param[ps->nr-1].c[i]=c[i];
}

static void add_imp_param(t_params *ps,int ai,int aj,int ak,int al,
			  real c0, real c1)
{
  srenew(ps->param,++ps->nr);
  ps->param[ps->nr-1].AI=ai;
  ps->param[ps->nr-1].AJ=aj;
  ps->param[ps->nr-1].AK=ak;
  ps->param[ps->nr-1].AL=al;
  ps->param[ps->nr-1].C0=c0;
  ps->param[ps->nr-1].C1=c1;
}

static int search_jtype(t_restp *rp,char *name,bool bFirstRes)
{
  int  j,k,mn,maxk,jmax;
  char *an,buf[12];
  
  strcpy(buf,name);
  
  /* Do a best match comparison */
  if ( bFirstRes && (buf[0] == 'H') && 
       ( (buf[1] == '1') || (buf[1] == '2') || (buf[1] == '3') ) )
    buf[1]='\0';
  maxk=0;
  jmax=-1;
  for(j=0; (j<rp->natom); j++) {
    an=*(rp->atomname[j]);
    if (strcasecmp(buf,an) == 0) {
      jmax=j;
      maxk=strlen(buf);
      break;
    }
    mn=min((int)strlen(buf),(int)strlen(an));
    for(k=0; (k<mn); k++) 
      if (buf[k] != an[k])
	break;
    if (k > maxk) {
      maxk=k;
      jmax=j;
    }
  }
  if (jmax == -1)
    fatal_error(0,"Atom %s not found in database in residue %s",
		buf,rp->resname);
  if (maxk != strlen(buf))
    fatal_error(0,"Atom %s not found in database in residue %s, "
		"it looks a bit like %s",
		buf,rp->resname,*(rp->atomname[jmax]));
  return jmax;
}

static bool missing_atoms(t_restp *rp, int resnr,
			  t_atoms *at, int i0, int i, bool bCTer)
{
  int  j,k;
  char *name;
  bool bFound, bRet;

  bRet=FALSE;
  for (j=0; j<rp->natom; j++) {
    name=*(rp->atomname[j]);
    if ((name[0]!='H') && (name[0]!='h') && (!bCTer || (name[0]!='O'))) {
      bFound=FALSE;
      for (k=i0; k<i; k++) 
	bFound=(bFound || !strcasecmp(*(at->atomname[k]),name));
      if (!bFound) {
	bRet=TRUE;
	fprintf(stderr,
		"\nWARNING: atom %s is missing in residue %s %d in the pdb file\n\n",
		name,*(at->resname[resnr]),resnr+1);
      }
    }
  }

  return bRet;
}

static void name2type(t_atoms *at,t_atomtype *atype,int nrtp,t_restp rtp[])
{
  int     i,j,resnr,i0;
  char    *name;
  t_restp *rp=NULL;
  bool    bFirst, bProt;

  bFirst=TRUE;
  resnr=-1;
  bProt=FALSE;
  i0=0;
  for(i=0; (i<at->nr); i++) {
    if (bFirst || (at->atom[i].resnr != resnr)) {
      if (bFirst) {
	bFirst=FALSE;
      }
      else {
	rp=search_rtp(*(at->resname[resnr]),nrtp,rtp);
	is_protein(rp->resname);
	missing_atoms(rp,resnr,at,i0,i,(!bProt || is_protein(rp->resname)));
	bProt=is_protein(*(at->resname[resnr]));
	i0=i;
      }
      resnr=at->atom[i].resnr;
      if ((rp=search_rtp(*(at->resname[resnr]),nrtp,rtp)) == NULL)
	fatal_error(0,"Residue type %s not found, jammer...",
		    *(at->resname[resnr]));
      if (strcasecmp(rp->resname,*(at->resname[resnr])) != 0)
	fatal_error(0,"Residue type %s not found",
		    *(at->resname[resnr]));
    }
    if (at->atom[i].m == 0) {
      name=*(at->atomname[i]);
      j=search_jtype(rp,name,(resnr == 0));
      at->atom[i].type  = rp->atom[j].type;
      at->atom[i].q     = rp->atom[j].q;
      at->atom[i].m     = atype->atom[rp->atom[j].type].m;
    }
    at->atom[i].typeB = at->atom[i].type;
    at->atom[i].qB    = at->atom[i].q;
    at->atom[i].mB    = at->atom[i].m;
  }
  missing_atoms(rp,resnr,at,i0,i,(!bProt || is_protein(rp->resname)));
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

static t_block *make_cgs(t_atoms *at,int nrtp,t_restp rtp[])
{
  int     *cgnr;
  t_block *cgs;
  int     i;
  int     curcg;
  double  qt,qq,qtot;
  
  fprintf(stderr,"Making charge groups...\n");
  snew(cgnr,at->nr);
  
  qtot=0;
  qt=0;
  curcg=0;
  for(i=0; (i<at->nr); i++) {
    qq=at->atom[i].q;
    qt+=qq;
    cgnr[i]=curcg;
    if (is_int(qt)) {
      qtot+=qt;
      qt=0;
      curcg++;
    }
  }
  if (qt != 0) {
    fprintf(stderr,"WARNING: Not an integer charge %g\n",qtot+qt);
    curcg++;
  }
  else
    fprintf(stderr,"Total charge on system: %g\n",qtot);
  
  snew(cgs,1);
  snew(cgs->index,curcg+1);
  snew(cgs->a,at->nr);
  for(i=0; (i<at->nr); i++)
    cgs->a[i]=i;
  
  cgs->index[0]=0;
  cgs->nr=0;
  for(i=1; (i<at->nr); i++) {
    if (cgnr[i] != cgnr[i-1]) {
      cgs->nr++;
      cgs->index[cgs->nr]=i;
    }
  }
  cgs->nr++;
  cgs->index[cgs->nr]=at->nr;
  cgs->nra=at->nr;
  
  sfree(cgnr);
  
  return cgs;
}

static void print_top_header(FILE *out, char *title, bool bITP, 
			     char *ff, int nincl, char **incls)
{
  int i;
  
  fprintf(out,"; This is your %stopology file\n",bITP ? "include " : "");
  fprintf(out,"; %s\n",title[0]?title:cool_quote());
  
  if (!bITP) {
    fprintf(out,"; Include forcefield constants\n");
    fprintf(out,"#include \"%s.itp\"\n\n",ff);
    if (nincl>0) {
      fprintf(out,"; Include other parts of system\n");
      for (i=0; (i<nincl); i++)
	fprintf(out,"#include \"%s\"\n",incls[i]);
      fprintf(out,"\n");
    }
  }
}

static void print_top_posre(FILE *out,char *pr)
{
  fprintf(out,"; Include Position restraint file\n");
  fprintf(out,"#ifdef POSRES\n");
  fprintf(out,"#include \"%s\"\n",pr);
  fprintf(out,"#endif\n\n");
}
  
static void print_top_water(FILE *out)
{
  fprintf(out,"; Include water topology\n");
  fprintf(out,"#ifdef FLEX_SPC\n");
  fprintf(out,"#include \"flexspc.itp\"\n");
  fprintf(out,"#else\n");
  fprintf(out,"#include \"spc.itp\"\n");
  fprintf(out,"#endif\n\n");
}

static void print_top_system(FILE *out)
{
  fprintf(out,"[ %s ]\n",dir2str(d_system));
  fprintf(out,"; Name\n");
  fprintf(out,"Protein in Water\n\n");
}

static void print_top_mols(FILE *out, int nmol, char **mols)
{
  int i;
  
  fprintf(out,"[ %s ]\n",dir2str(d_molecules));
  fprintf(out,"; %-15s %5s\n","Compound","#mols");
  if (nmol==0)
    fprintf(out,"Protein\t\t1\n");
  else
    for (i=0; (i<nmol); i++)
      fprintf(out,"%-15s %5d\n",mols[i],1);
  
}

void write_top(char *ff,char *fn, char *pr,char *title, char *molname,
	       int nincl, char **incls, int nmol, char **mols,
	       t_atoms *at,t_params plist[],int nrtp,t_restp rtp[],
	       t_atomtype *atype,t_block *cgs)
{
  FILE *out;
  bool bITP;
  
  bITP=(fn2ftp(fn)==efITP);
  out=ffopen(fn,"w");
  
  print_top_header(out,title,bITP,ff,nincl,incls);
  
  if (at && atype && cgs) {
    fprintf(out,"[ %s ]\n",dir2str(d_moleculetype));
    fprintf(out,"; %-15s %5s\n","Name","nrexcl");
    fprintf(out,"%-15s %5d\n\n",molname?molname:"Protein",3);
    
    print_atoms(out,atype,at,cgs);
    print_bonds(out,at->nr,d_bonds,F_BONDS,plist,FALSE);
    print_bonds(out,at->nr,d_pairs,F_LJ14,plist,FALSE);
    print_bonds(out,at->nr,d_angles,F_ANGLES,plist,FALSE);
    print_bonds(out,at->nr,d_dihedrals,F_PDIHS,plist,FALSE);
    print_bonds(out,at->nr,d_dihedrals,F_IDIHS,plist,FALSE);
    
    if (pr)
      print_top_posre(out,pr);
  }
  
  if (!bITP) {
    print_top_water(out);
    print_top_system(out);
    print_top_mols(out, nmol, mols);
  }
  
  fclose(out);
}

static int search_res_atom(char *type,int resnr,
			   int natom,t_atom at[],char **aname[])
{
  int i;

  for(i=0; (i<natom); i++)
    if (at[i].resnr == resnr)
      return search_atom(type,i,natom,aname);
  
  return -1;
}

static void do_ssbonds(t_params *ps,int natoms,t_atom atom[],char **aname[],
		       int nssbonds,t_ssbond *ssbonds)
{
  int i,ri,rj,ai,aj;
  
  for(i=0; (i<nssbonds); i++) {
    ri = ssbonds[i].res1;
    rj = ssbonds[i].res2;
    ai = search_res_atom(ssbonds[i].a1,ri,natoms,atom,aname);
    aj = search_res_atom(ssbonds[i].a2,rj,natoms,atom,aname);
    if ((ai == -1) || (aj == -1))
      fatal_error(0,"Trying to make impossible special bond (%s-%s)!",
		  ssbonds[i].a1,ssbonds[i].a2);
    add_param(ps,ai,aj,NULL);
  }
}

/*
static void heme_his_bonds(t_params *ps,
			   int natoms,t_atom atom[],char **aname[],
			   int nres,char **resname[],
			   rvec x[])
{
  int  nhis1,nheme;
  int  *his1p,*ne2p;
  int  *hemep,*fep;
  int  i,j;
  int  ai,aj;
  real **d;
  
  fprintf(stderr,"Checking for Heme-His bonds...\n");
  snew(his1p,nres);
  snew(ne2p,nres);
  snew(hemep,nres);
  snew(fep,nres);

  nhis1 = nheme = 0;
  for(i=0; (i<nres); i++) {
    if (strcasecmp(*(resname[i]),"HIS1") == 0) {
      his1p[nhis1]=i;
      if ((ne2p[nhis1]=search_res_atom("NE2",i,natoms,atom,aname)) == -1)
	HALT("No NE2 in HIS1 residue....");
      nhis1++;
    }
    else if (strcasecmp(*(resname[i]),"HEME") == 0) {
      hemep[nheme]=i;
      if ((fep[nheme]=search_res_atom("FE",i,natoms,atom,aname)) == -1)
	HALT("No FE in HEME residue....");
      nheme++;
    }
  }
  if (nheme != nhis1)
    fatal_error(0,"You have %d HEME groups and %d HIS1 residues, try again",
		nheme,nhis1);
  if (nheme > 0) {
    snew(d,nheme);
    for(i=0; (i<nheme); i++)
      snew(d[i],nheme);
    
    for(i=0; (i<nheme); i++) 
      for(j=0; (j<nheme); j++) {
	ai=ne2p[i];
	aj=fep[j];
	d[i][j]=distance(x[ai],x[aj]);
      }
    fprintf(stderr,"Heme-His1 Distance matrix\n%8s","");
    for(i=0; (i<nheme); i++) 
      fprintf(stderr,"%4d    ",his1p[i]+1);
    fprintf(stderr,"\n");
    for(i=0; (i<nheme); i++) {
      fprintf(stderr,"%4d    ",hemep[i]+1);
      for(j=0; (j<i); j++)
	fprintf(stderr,"%8s","");
      for( ; (j<nheme); j++)
	fprintf(stderr,"%6.3f  ",d[i][j]);
      fprintf(stderr,"\n");
    }
    while (min_dist(nheme,d,&i,&j,FALSE) < 0.35) {
      fprintf(stderr,"Linking Heme %d and His1 %d...\n",hemep[i]+1,his1p[j]+1);
      ai=fep[i];
      aj=ne2p[j];
      add_param(ps,ai,aj,NULL);
      d[i][j] = 10.0;
    }
    
    for(i=0; (i<nheme); i++)
      sfree(d[i]);
    sfree(d);
  }
  else
    fprintf(stderr,"No Heme\n");
    
  sfree(hemep);
  sfree(his1p);
  sfree(fep);
  sfree(ne2p);
}
*/

static void ter2bonds(t_params *ps,
		      int natoms,t_atom atom[],char **aname[],
		      int resnr,t_terblock *tdb)
{
  int i,j,k;
  
  /* Skip to right residue */
  for(i=0; ((i<natoms) && (atom[i].resnr != resnr)); i++)
    ;
  /* Now only for this residue */
  for(   ; ((i<natoms) && (atom[i].resnr == resnr)); i++) {
    for(j=0; (j<tdb->nadd); j++) {
      if (strcmp(tdb->ab[j].na[0],*(aname[i])) == 0) {
	if (tdb->ab[j].tp == 9) {     /* COOH terminus */
	  add_param(ps,i,i+1,NULL);   /* C-O  */
	  add_param(ps,i,i+2,NULL);   /* C-OA */
	  add_param(ps,i+2,i+3,NULL); /* OA-H */
	} else {
	  for(k=0; (k<tdb->ab[j].nh); k++)
	    add_param(ps,i,i+k+1,NULL);
	}
      }
    }
  }
}

static void ter2idihs(t_params *ps,
		      int natoms,t_atom atom[],char **aname[],
		      int resnr,t_terblock *tdb)
{
  int i,j,k,aa0,a0[4],start;
  char *atomnm;
  
  /* Skip to right residue */
  start=0;
  while ((start<natoms) && (atom[start].resnr != resnr)) start++;
  /* Now only for this residue */
  for(j=0; (j<tdb->nidih); j++) {
    for(k=0; (k<4); k++) {
      atomnm=tdb->idih[j].ai[k];
      if ((aa0=search_atom(atomnm,start,natoms,aname)) == -1)
	fatal_error(0,"Trying to add improper %s %s %s %s to res %d,"
		    "atom %s not found",
		    tdb->idih[j].ai[0],tdb->idih[j].ai[1],tdb->idih[j].ai[2],
		    tdb->idih[j].ai[3],resnr+1,
		    tdb->idih[j].ai[k]);
      else {
	if (atom[aa0].resnr > (resnr+1)) {
	  fprintf(stderr,"WARNING: improper spans more than 2 residues:\n");
	  for (i=0; (i<k); i++)
	    fprintf(stderr,"atom %d %s, res %d\n",
		    a0[i],*aname[a0[i]],atom[a0[i]].resnr+1);
	}
	a0[k] = aa0;
      }
    }
    add_imp_param(ps,a0[0],a0[1],a0[2],a0[3],NOTSET,NOTSET);
  }
}

static void at2bonds(t_params *ps,
		     int nrb,t_resbond rb[],int nah,t_addh ah[],
		     int natoms,t_atom atom[],char **aname[],
		     int nres,char **resname[],bool bAlldih,
		     rvec x[])
{
  t_resbond *rb0;
  t_addh    *ah0;
  t_add_block *ab0;
  int        i,j,k,l;
  int        ai,aj;
  real       ddd;
  
  /* First generate bonds along the heavy atoms */
  fprintf(stderr,"Making bonds on the heavy atoms...\n");
  for(i=j=0; ((i<nres) && (j<natoms)); i++) {
    if ((rb0=search_rb(*(resname[i]),nrb,rb)) != NULL) {
      for(k=0; (k<rb0->nb); k++) {
	if (bAlldih) {
	  /* if bAlldih is true, don't check for hydrogens */
	  if ((ai=search_atom(rb0->rbond[k].ai,j,natoms,aname)) != -1) 
	    if ((aj=search_atom(rb0->rbond[k].aj,j,natoms,aname)) != -1)
	      if ((!is_hydrogen(rb0->rbond[k].ai) && 
		  !is_hydrogen(rb0->rbond[k].aj)) ||
		  (atom[ai].resnr == atom[aj].resnr))
		add_param(ps,ai,aj,rb0->rbond[k].c);
	} else {     
	  /* else do the normal thing */
	  if (!is_hydrogen(rb0->rbond[k].ai) && 
	      !is_hydrogen(rb0->rbond[k].aj)) {
	    if ((ai=search_atom(rb0->rbond[k].ai,j,natoms,aname)) != -1) 
	      if ((aj=search_atom(rb0->rbond[k].aj,j,natoms,aname)) != -1){
		ddd = distance(x[ai],x[aj]);
		if (ddd > 0.25)
		  fprintf(stderr,"Warning: Long Bond (%d-%d = %g nm)\n",
			  ai+1,aj+1,ddd);
		add_param(ps,ai,aj,rb0->rbond[k].c);
	      }
	  }
	}
      }
    }
    else
      fprintf(stderr,"No bond information for residue %d (%s)\n",
	      i+1,*(resname[i]));
    while ((j < natoms-1) && (atom[j].resnr == atom[j+1].resnr))
      j++;
    j++;
  }
  /* And now the hydrogen atoms. */
  fprintf(stderr,"Making bonds on the hydrogen atoms...\n");
  for(i=j=0; ((i<nres) && (j<natoms)); i++) {
    if ((ah0=search_h_db(nah,ah,*(resname[i]))) != NULL) {
      /* Loop over blocks */
      for(k=0; (k<ah0->n_add); k++) {
	ab0=&(ah0->ab[k]);
	if ((ai=search_atom(ab0->na[0],j,natoms,aname)) == -1)
	  fprintf(stderr,"WARNING: Atom %s not found in residue %s %d\n",
		  ab0->na[0],*(resname[i]),i);
	else {
	  for(l=0; (l<ab0->nh); l++)
	    add_param(ps,ai,ai+l+1,NULL);
	}
      }
    }

    while ((j < natoms-1) && (atom[j].resnr == atom[j+1].resnr))
      j++;
    j++;
  }
}

static int pcompar(const void *a, const void *b)
{
  t_param *pa,*pb;
  int     id;
  pa=(t_param *)a;
  pb=(t_param *)b;
  
  id=pa->a[0]-pb->a[0];
  if (id == 0)
    return pa->a[1]-pb->a[1];
  else
    return id;
}

static void clean_bonds(t_params *ps)
{
  int     i,j;
  atom_id ai,aj;
  
  if (ps->nr > 0) {
    for(i=0; (i<ps->nr); i++) {
      ai=ps->param[i].AI;
      aj=ps->param[i].AJ;
      if (aj < ai) {
	ps->param[i].AI=aj;
	ps->param[i].AJ=ai;
      }
    }
    
    /* Sort bonds */
    qsort(ps->param,ps->nr,(size_t)sizeof(ps->param[0]),pcompar);
    
    /* remove doubles */
    for(i=j=1; (i<ps->nr); i++) {
      if (pcompar(&(ps->param[i]),&(ps->param[j-1])) != 0) {
	memcpy(&(ps->param[j]),&(ps->param[i]),(size_t)sizeof(ps->param[0]));
	j++;
      }
    }
    fprintf(stderr,"Number of bonds was %d, now %d\n",ps->nr,j);
    ps->nr=j;
  } else
    fprintf(stderr,"No bonds\n");
}

static void print_mass(t_atoms *atoms)
{
  double m;
  int    i;
  
  m=0.0;
  for(i=0; (i<atoms->nr); i++)
    m+=atoms->atom[i].m;
    
  fprintf(stderr,"Total mass in system %g a.m.u.\n",m);
}

void pdb2top(char *ff,char *fn,char *pr,char *title,char *molname,
	     int nincl, char **incls, int nmol, char **mols,
	     t_atoms *atoms,int nah,t_addh ah[],rvec x[],
	     t_atomtype *atype,t_symtab *tab,
	     int nrb, t_resbond rb[],
	     int nrtp,t_restp rtp[],
	     int nra, t_resang ra[],
	     int nid, t_idihres idi[],
	     t_terblock *ntdb,t_terblock *ctdb,
	     bool bH14,
	     int rn,int rc,bool bAlldih,
	     int nssbonds,t_ssbond *ssbonds)
{
  t_params plist[F_NRE];
  t_nextnb nnb;
  t_block  *cgs;
  int i;
  
  init_plist(plist);

  /* Make bonds */
  at2bonds(&(plist[F_BONDS]),
	   nrb,rb,nah,ah,atoms->nr,atoms->atom,atoms->atomname,
	   atoms->nres,atoms->resname,bAlldih,x);
  
  /* Terminal bonds */
  if (rn>=0)
    ter2bonds(&(plist[F_BONDS]),atoms->nr,atoms->atom,atoms->atomname,
	      rn,ntdb);
  if (rc>=0)
    ter2bonds(&(plist[F_BONDS]),atoms->nr,atoms->atom,atoms->atomname,
	      rc,ctdb);
  if (rn>=0)
    ter2idihs(&(plist[F_IDIHS]),atoms->nr,atoms->atom,atoms->atomname,
	      rn,ntdb);
  if (rc>=0)
    ter2idihs(&(plist[F_IDIHS]),atoms->nr,atoms->atom,atoms->atomname,
	      rc,ctdb);
  
  /* Last the disulphide bonds & heme-his */
  do_ssbonds(&(plist[F_BONDS]),atoms->nr,atoms->atom,
	     atoms->atomname,nssbonds,ssbonds);
  
  /*heme_his_bonds(&(plist[F_BONDS]),atoms->nr,atoms->atom,
    atoms->atomname,
    atoms->nres,atoms->resname,x);
  */
  
  name2type(atoms,atype,nrtp,rtp);
  
  /* Cleanup bonds (sort and rm doubles) */ 
  clean_bonds(&(plist[F_BONDS]));

  fprintf(stderr,"Generating angles and dihedrals...\n");
  init_nnb(&nnb,atoms->nr,4);
  gen_nnb(&nnb,plist);
  print_nnb(&nnb,"NNB");
  gen_pad(&nnb,atoms,bH14,plist,nrtp,rtp,nra,ra,nid,idi,bAlldih);
  fprintf(stderr,
	  "There are %6d dihedrals, %6d impropers, %6d angles\n"
	  "          %6d pairs and  %6d bonds\n",
	  plist[F_PDIHS].nr,plist[F_IDIHS].nr,
	  plist[F_ANGLES].nr,plist[F_LJ14].nr,plist[F_BONDS].nr);
  
  cgs=make_cgs(atoms,nrtp,rtp);
  print_mass(atoms);
  
  fprintf(stderr,"Writing topology file\n");
  write_top(ff,fn,pr,title,molname,nincl,incls,nmol,mols,
	    atoms,plist,nrtp,rtp,atype,cgs);
  
  for (i=0; (i<F_NRE); i++)
    sfree(plist[i].param);
  done_nnb(&nnb);
  done_block(cgs);
  sfree(cgs);
}
