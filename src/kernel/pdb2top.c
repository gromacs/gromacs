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
#include "pdb2gmx.h"
#include "gen_dum.h"
#include "add_par.h"

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
	fprintf(stderr,"\nWARNING: "
		"atom %s is missing in residue %s %d in the pdb file\n\n",
		name,*(at->resname[resnr]),resnr+1);
      }
    }
  }

  return bRet;
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

static void name2type(t_atoms *at,int **cgnr,
		      t_atomtype *atype,int nrtp,t_restp rtp[])
{
  int     i,j,prevresnr,resnr,i0,prevcg,cg,curcg;
  char    *name;
  t_restp *rp=NULL,*prevrp=NULL;
  bool    bProt, bNterm;
  double  qt,qq;

  resnr=-1;
  bProt=FALSE;
  i0=0;
  snew(*cgnr,at->nr);
  qt=0;
  curcg=0;
  cg=-1;
  for(i=0; (i<at->nr); i++) {
    prevresnr=resnr;
    if (at->atom[i].resnr != resnr) {
      resnr=at->atom[i].resnr;
      bProt=is_protein(*(at->resname[resnr]));
      bNterm=bProt && (resnr == 0);
      prevrp=rp;
      rp=search_rtp(*(at->resname[resnr]),nrtp,rtp);
      if (prevrp)
	missing_atoms(prevrp,resnr,at,i0,i,
		      (!bProt && is_protein(prevrp->resname)));
      i0=i;
    }
    if (at->atom[i].m == 0) {
      qt=0;
      prevcg=cg;
      name=*(at->atomname[i]);
      j=search_jtype(rp,name,bNterm);
      at->atom[i].type  = rp->atom[j].type;
      at->atom[i].q     = rp->atom[j].q;
      at->atom[i].m     = atype->atom[rp->atom[j].type].m;
      cg = rp->cgnr[j];
      if ( (cg != prevcg) || (resnr != prevresnr) )
	curcg++;
    } else {
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
  missing_atoms(rp,resnr,at,i0,i,(!bProt || is_protein(rp->resname)));
}

static void print_top_heavy_H(FILE *out, real mHmult)
{
  if (mHmult!=1.0) {
    fprintf(out,"; heavy hydrogens:\n");
    if (mHmult==4.0)
      fprintf(out,"#define HEAVY_H\n\n");
    else
      fprintf(stderr,"WARNING: unsupported proton mass multiplier (%g) "
	      "in pdb2top\n",mHmult);
  }
}

static void print_top_header(FILE *out, char *title, bool bITP, 
			     char *ff, int nincl, char **incls, real mHmult)
{
  int i;
  
  fprintf(out,"; This is your %stopology file\n",bITP ? "include " : "");
  fprintf(out,"; %s\n",title[0]?title:cool_quote());
  
  if (!bITP) {
    print_top_heavy_H(out, mHmult);
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
	       t_atoms *at,t_params plist[],t_block *excl,
	       t_atomtype *atype,int *cgnr, int nrexcl, real mHmult)
     /* NOTE: nrexcl is not the size of *excl! */
{
  FILE *out;
  bool bITP;
  
  bITP=(fn2ftp(fn)==efITP);
  out=ffopen(fn,"w");
  
  print_top_header(out,title,bITP,ff,nincl,incls,mHmult);
  
  if (at && atype && cgnr) {
    fprintf(out,"[ %s ]\n",dir2str(d_moleculetype));
    fprintf(out,"; %-15s %5s\n","Name","nrexcl");
    fprintf(out,"%-15s %5d\n\n",molname?molname:"Protein",nrexcl);
    
    print_atoms(out,atype,at,cgnr);
    print_bonds(out,at->nr,d_bonds,F_BONDS,plist,FALSE);
    print_bonds(out,at->nr,d_pairs,F_LJ14,plist,FALSE);
    print_bonds(out,at->nr,d_angles,F_ANGLES,plist,FALSE);
    print_bonds(out,at->nr,d_dihedrals,F_PDIHS,plist,FALSE);
    print_bonds(out,at->nr,d_dihedrals,F_IDIHS,plist,FALSE);
    print_bonds(out,at->nr,d_dum2,F_DUMMY2,plist,FALSE);
    print_bonds(out,at->nr,d_dum2,F_DUMMY2FD,plist,FALSE);
    print_bonds(out,at->nr,d_dum2,F_DUMMY2FAD,plist,FALSE);
    print_bonds(out,at->nr,d_dum3,F_DUMMY3,plist,FALSE);
    
    if (excl)
      if (excl->nr > 0)
	print_excl(out,excl);
    
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

static void ter2bonds(t_params *ps,
		      int natoms,t_atom atom[],char **aname[],
		      int resnr,t_hackblock *tdb)
{
  int i,j,k;
  
  /* Skip to right residue */
  for(i=0; ((i<natoms) && (atom[i].resnr != resnr)); i++)
    ;
  /* Now only for this residue */
  for(   ; ((i<natoms) && (atom[i].resnr == resnr)); i++)
    for(j=0; (j<tdb->nadd); j++)
      if (strcmp(tdb->ab[j].na[0],*(aname[i])) == 0)
	if (tdb->ab[j].tp == 9) {         /* COOH terminus */
	  add_param(ps,i,i+1,NULL);   /* C-O  */
	  add_param(ps,i,i+2,NULL);   /* C-OA */
	  add_param(ps,i+2,i+3,NULL);   /* OA-H */
	} else
	  for(k=0; (k<tdb->ab[j].nh); k++)
	    add_param(ps,i,i+k+1,NULL);
}

static void ter2idihs(t_params *ps,
		      int natoms,t_atom atom[],char **aname[],
		      int resnr,t_hackblock *tdb)
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
	    a0[i]+1,*aname[a0[i]],atom[a0[i]].resnr+1);
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
  real       dist2;
#define LONG_BOND_DIST 0.25
  
  /* First generate bonds along the heavy atoms */
  fprintf(stderr,"Making bonds on the heavy atoms...\n");
  for(i=j=0; ((i<nres) && (j<natoms)); i++) {
    if ((rb0=search_rb(*(resname[i]),nrb,rb)) != NULL) {
      for(k=0; (k<rb0->nb); k++) {
	if (bAlldih) {
	  /* if bAlldih is true, don't check for hydrogens */
	  if ((ai=search_atom(rb0->rbond[k].ai,j,natoms,aname)) != -1)
	    if ((aj=search_atom(rb0->rbond[k].aj,j,natoms,aname)) != -1)
	      if ( ( !is_hydrogen(rb0->rbond[k].ai) &&
		     !is_hydrogen(rb0->rbond[k].aj) ) ||
		   (atom[ai].resnr == atom[aj].resnr))
		add_param(ps,ai,aj,rb0->rbond[k].c);
	} else {
	  /* else do the normal thing */
	  if (!is_hydrogen(rb0->rbond[k].ai) && 
	      !is_hydrogen(rb0->rbond[k].aj)) {
	    if ((ai=search_atom(rb0->rbond[k].ai,j,natoms,aname)) != -1)
	      if ((aj=search_atom(rb0->rbond[k].aj,j,natoms,aname)) != -1) {
		dist2 = distance2(x[ai],x[aj]);
		if (dist2 > sqr(LONG_BOND_DIST) )
		  fprintf(stderr,"Warning: Long Bond (%d-%d = %g nm)\n",
			  ai+1,aj+1,sqrt(dist2));
		add_param(ps,ai,aj,rb0->rbond[k].c);
	      }
	  }
	}
      }
    } else
      fprintf(stderr,"No bond information for residue %d (%s)\n",
	      i+1,*(resname[i]));
    while ((j < natoms-1) && (atom[j].resnr == atom[j+1].resnr))
      j++;
    j++;
  }
  /* And now the hydrogen atoms. */
  fprintf(stderr,"Making bonds on the hydrogen atoms...\n");
  for(i=j=0; ((i<nres) && (j<natoms)); i++) {
    if ((ah0=search_h_db(nah,ah,*(resname[i]))) != NULL)
      /* Loop over blocks */
      for(k=0; (k<ah0->n_add); k++) {
	ab0=&(ah0->ab[k]);
	if ((ai=search_atom(ab0->na[0],j,natoms,aname)) == -1)
	  fprintf(stderr,"WARNING: Atom %s not found in residue %s %d\n",
		  ab0->na[0]+1,*(resname[i]),i+1);
	else
	  for(l=0; (l<ab0->nh); l++) {
	    aj=ai+l+1;
	    add_param(ps,ai,aj,NULL);
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

void print_sums(t_atoms *atoms, bool bSystem)
{
  double m,qtot;
  int    i;
  char   *where;
  
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

void pdb2top(char *ff,char *fn,char *pr,char *title,char *molname,
	     int nincl, char **incls, int nmol, char **mols,
	     t_atoms *atoms,int nah,t_addh ah[],rvec **x,
	     t_atomtype *atype,t_symtab *tab,
	     int nrb, t_resbond rb[],
	     int nrtp,t_restp rtp[],
	     int nra, t_resang ra[],
	     int nrd, t_resdih rd[],
	     int nid, t_idihres idi[],
	     t_hackblock *ntdb,t_hackblock *ctdb,
	     bool bH14, int rn,int rc,bool bAlldih,
	     int nddb, t_dumblock *ddb,bool bDummies, real mHmult,
	     int nssbonds,t_ssbond *ssbonds, int nrexcl)
{
  t_params plist[F_NRE];
  t_nextnb nnb;
  int      *cgnr;
  t_block  excl;
  bool     *is_dummy;
  int      i;
  
  /* determine which atoms will be dummies and make space for dummy masses */
  if (bDummies)
    do_dummies(atoms,atype,tab,nrtp,rtp,x,nddb,ddb,&is_dummy,mHmult);
    
  init_plist(plist);
  
  /* Make bonds */
  at2bonds(&(plist[F_BONDS]),
	   nrb,rb,nah,ah,atoms->nr,atoms->atom,atoms->atomname,
	   atoms->nres,atoms->resname,bAlldih,*x);
  
  /* Terminal bonds */
  if (rn>=0)
    ter2bonds(&(plist[F_BONDS]),atoms->nr,atoms->atom,atoms->atomname,rn,ntdb);
  if (rc>=0)
    ter2bonds(&(plist[F_BONDS]),atoms->nr,atoms->atom,atoms->atomname,rc,ctdb);
  if (rn>=0)
    ter2idihs(&(plist[F_IDIHS]),atoms->nr,atoms->atom,atoms->atomname,rn,ntdb);
  if (rc>=0)
    ter2idihs(&(plist[F_IDIHS]),atoms->nr,atoms->atom,atoms->atomname,rc,ctdb);
  
  /* Last the disulphide bonds & heme-his */
  do_ssbonds(&(plist[F_BONDS]),atoms->nr,atoms->atom,
	     atoms->atomname,nssbonds,ssbonds);
  
  /*heme_his_bonds(&(plist[F_BONDS]),atoms->nr,atoms->atom,
    atoms->atomname,
    atoms->nres,atoms->resname,x);
  */
  
  name2type(atoms,&cgnr,atype,nrtp,rtp);
  
  /* Cleanup bonds (sort and rm doubles) */ 
  clean_bonds(&(plist[F_BONDS]));

  fprintf(stderr,"Generating angles and dihedrals...\n");
  init_nnb(&nnb,atoms->nr,4);
  gen_nnb(&nnb,plist);
  print_nnb(&nnb,"NNB");
  gen_pad(&nnb,atoms,bH14,plist,nrtp,rtp,nra,ra,nrd,rd,nid,idi,bAlldih);
  done_nnb(&nnb);

  /* Initiate the exclusion block, must also be done when no dummies are
   * generated!
   */
  init_block(&excl);
  
  if (bDummies) {
    /* generate exclusions for dummy masses and atoms */ 
    init_nnb(&nnb,atoms->nr,nrexcl);
    gen_nnb(&nnb,plist);
    excl.nr=atoms->nr;
    nnb2excl(&nnb,&excl);
    done_nnb(&nnb);
    
    /* only keep exclusions for dummies: */
    do_dum_excl(&excl,is_dummy);
    
    /* remove bonds to dummy atoms and add bonds to dummy masses */
    /* also generate dummy parameters and set masses of dummies to 0 */
    do_dum_top(
	       &(plist[F_BONDS]),&(plist[F_DUMMY2]),&(plist[F_DUMMY3]),
	       &(plist[F_DUMMY2FD]),&(plist[F_DUMMY2FAD]),&(plist[F_ANGLES]),
	       nddb,ddb,is_dummy,atoms,atype,nrtp,rtp,mHmult);
    
    /* remove things with dummy atoms */
    clean_dum_angles(&(plist[F_ANGLES]), plist, is_dummy);
    clean_dum_dihs(&(plist[F_PDIHS ]), atoms->nr, "proper",   plist, is_dummy);
    clean_dum_dihs(&(plist[F_IDIHS ]), atoms->nr, "improper", plist, is_dummy);
    
    sfree(is_dummy);
  }
  /* Cleanup bonds (sort and rm doubles) */ 
  clean_bonds(&(plist[F_BONDS]));
  
  fprintf(stderr,
	  "There are %4d dihedrals, %4d impropers, %4d angles\n"
	  "          %4d pairs,     %4d bonds and  %4d dummies\n",
	  plist[F_PDIHS].nr,plist[F_IDIHS].nr,plist[F_ANGLES].nr,
	  plist[F_LJ14].nr,plist[F_BONDS].nr,
	  plist[F_DUMMY1].nr+plist[F_DUMMY2].nr+plist[F_DUMMY2FD].nr+
	  plist[F_DUMMY2FAD].nr+plist[F_DUMMY3].nr);
  
  print_sums(atoms, FALSE);
  
  fprintf(stderr,"Writing topology file\n");
  write_top(ff,fn,pr,title,molname,nincl,incls,nmol,mols,
	    atoms,plist,&excl,atype,cgnr,nrexcl,mHmult);
  
  for (i=0; (i<F_NRE); i++)
    sfree(plist[i].param);
}
