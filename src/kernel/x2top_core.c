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

#include <ctype.h>
#include "maths.h"
#include "macros.h"
#include "copyrite.h"
#include "bondf.h"
#include "string2.h"
#include "smalloc.h"
#include "strdb.h"
#include "sysstuff.h"
#include "confio.h"
#include "physics.h"
#include "statutil.h"
#include "vec.h"
#include "random.h"
#include "3dview.h"
#include "txtdump.h"
#include "readinp.h"
#include "names.h"
#include "toppush.h"
#include "pdb2top.h"
#include "symtab.h"
#include "gen_ad.h"
#include "topexcl.h"
#include "vec.h"
#include "x2top_nm2type.h"
#include "x2top_core.h"
#include "atomprop.h"
#include "grompp.h"
#include "add_par.h"
#include "gmx_random.h"
#include "gpp_atomtype.h"

void mk_bonds(x2top_nm2t nmt,
	      t_atoms *atoms,rvec x[],t_params *bond,int nbond[],char *ff,
	      bool bPBC,matrix box,void *atomprop,real tol)
{
  t_param b;
  int     i,j;
  t_pbc   pbc;
  rvec    dx;
  real    dx1;
  
  for(i=0; (i<MAXATOMLIST); i++)
    b.a[i] = -1;
  for(i=0; (i<MAXFORCEPARAM); i++)
    b.c[i] = 0.0;
    
  if (bPBC)
    set_pbc(&pbc,box);
  for(i=0; (i<atoms->nr); i++) {
    if ((i % 10) == 0)
      fprintf(stderr,"\ratom %d",i);
    for(j=i+1; (j<atoms->nr); j++) {
      if (bPBC)
	pbc_dx(&pbc,x[i],x[j],dx);
      else
	rvec_sub(x[i],x[j],dx);
      
      dx1 = norm(dx);
      if (is_bond(nmt,*atoms->atomname[i],*atoms->atomname[j],dx1,tol) ||
	  is_bond(nmt,*atoms->atomname[j],*atoms->atomname[i],dx1,tol)) {
	b.AI = i;
	b.AJ = j;
	b.C0 = dx1;
	push_bondnow (bond,&b);
	nbond[i]++;
	nbond[j]++;
	if (debug) 
	  fprintf(debug,"Bonding atoms %s-%d and %s-%d at distance %g\n",
		  *atoms->atomname[i],i+1,*atoms->atomname[j],j+1,dx1);
      }
    }
  }
  fprintf(stderr,"\ratom %d\n",i);
}

int *set_cgnr(t_atoms *atoms,bool bUsePDBcharge,real *qtot,real *mtot)
{
  int    i,n=1;
  int    *cgnr;
  double qt=0,mt=0;
  
  *qtot = *mtot = 0;
  snew(cgnr,atoms->nr);
  for(i=0; (i<atoms->nr); i++) {
    if (atoms->pdbinfo && bUsePDBcharge)
      atoms->atom[i].q = atoms->pdbinfo[i].bfac;
    qt += atoms->atom[i].q;
    *qtot += atoms->atom[i].q;
    *mtot += atoms->atom[i].m;
    cgnr[i] = n;
    if (is_int(qt)) {
      n++;
      qt=0;
    }
  }
  return cgnr;
}

t_atomtype set_atom_type(t_symtab *tab,t_atoms *atoms,t_params *bonds,
			 int *nbonds,x2top_nm2t nm2t,void *atomprop)
{
  t_atomtype atype;
  int nresolved;
  
  atype = init_atomtype();
  snew(atoms->atomtype,atoms->nr);
  nresolved = nm2type(nm2t,tab,atoms,atype,nbonds,bonds,atomprop);
  if (nresolved != atoms->nr)
    gmx_fatal(FARGS,"Could only find a forcefield type for %d out of %d atoms",
	      nresolved,atoms->nr);
  
  fprintf(stderr,"There are %d different atom types in your sample\n",
	  get_atomtype_ntypes(atype));
    
  return atype;
}

void lo_set_force_const(t_params *plist,real c[],int nrfp,bool bRound,
			bool bDih,bool bParam)
{
  int    i,j;
  double cc;
  char   buf[32];
  
  for(i=0; (i<plist->nr); i++) {
    if (!bParam)
      for(j=0; j<nrfp; j++)
	c[j] = NOTSET;
    else {
      if (bRound) {
	sprintf(buf,"%.2e",plist->param[i].c[0]);
	sscanf(buf,"%lf",&cc);
	c[0] = cc;
      }
      else 
	c[0] = plist->param[i].c[0];
      if (bDih) {
	c[0] *= c[2];
	c[0] = ((int)(c[0] + 3600)) % 360;
	if (c[0] > 180)
	  c[0] -= 360;
	/* To put the minimum at the current angle rather than the maximum */
	c[0] += 180; 
      }
    }
    for(j=0; (j<nrfp); j++) {
      plist->param[i].c[j]      = c[j];
      plist->param[i].c[nrfp+j] = c[j];
    }
    set_p_string(&(plist->param[i]),"");
  }
}

void set_force_const(t_params plist[],real kb,real kt,real kp,bool bRound,
		     bool bParam)
{
  int i;
  real c[MAXFORCEPARAM];
  
  c[0] = 0;
  c[1] = kb;
  lo_set_force_const(&plist[F_BONDS],c,2,bRound,FALSE,bParam);
  c[1] = kt;
  lo_set_force_const(&plist[F_ANGLES],c,2,bRound,FALSE,bParam);
  c[1] = kp;
  c[2] = 3;
  lo_set_force_const(&plist[F_PDIHS],c,3,bRound,TRUE,bParam);
}

void calc_angles_dihs(t_params *ang,t_params *dih,rvec x[],bool bPBC,
		      matrix box)
{
  int    i,ai,aj,ak,al,t1,t2,t3;
  rvec   r_ij,r_kj,r_kl,m,n;
  real   sign,th,costh,ph,cosph;
  t_pbc  pbc;

  if (bPBC)
    set_pbc(&pbc,box);
  if (debug)
    pr_rvecs(debug,0,"X2TOP",box,DIM);
  for(i=0; (i<ang->nr); i++) {
    ai = ang->param[i].AI;
    aj = ang->param[i].AJ;
    ak = ang->param[i].AK;
    th = RAD2DEG*bond_angle(x[ai],x[aj],x[ak],bPBC ? &pbc : NULL,
			    r_ij,r_kj,&costh,&t1,&t2);
    if (debug)
      fprintf(debug,"X2TOP: ai=%3d aj=%3d ak=%3d r_ij=%8.3f r_kj=%8.3f th=%8.3f\n",
	      ai,aj,ak,norm(r_ij),norm(r_kj),th);
    ang->param[i].C0 = th;
  }
  for(i=0; (i<dih->nr); i++) {
    ai = dih->param[i].AI;
    aj = dih->param[i].AJ;
    ak = dih->param[i].AK;
    al = dih->param[i].AL;
    ph = RAD2DEG*dih_angle(x[ai],x[aj],x[ak],x[al],bPBC ? & pbc : NULL,
			   r_ij,r_kj,r_kl,m,n,&cosph,&sign,&t1,&t2,&t3);
    if (debug)
      fprintf(debug,"X2TOP: ai=%3d aj=%3d ak=%3d al=%3d r_ij=%8.3f r_kj=%8.3f r_kl=%8.3f ph=%8.3f\n",
	      ai,aj,ak,al,norm(r_ij),norm(r_kj),norm(r_kl),ph);
    dih->param[i].C0 = ph;
  }
}

void dump_hybridization(FILE *fp,t_atoms *atoms,int nbonds[])
{
  int i;
  
  for(i=0; (i<atoms->nr); i++) {
    fprintf(fp,"Atom %5s has %1d bonds\n",*atoms->atomname[i],nbonds[i]);
  }
}

static void print_pl(FILE *fp,t_params plist[],int ftp,char *name,
		     char ***atomname)
{ 
  int i,j,nral,nrfp;

  if (plist[ftp].nr > 0) {
    fprintf(fp,"\n");
    fprintf(fp,"[ %s ]\n",name);
    nral = interaction_function[ftp].nratoms;
    nrfp = interaction_function[ftp].nrfpA;
    for(i=0; (i<plist[ftp].nr); i++) {
      for(j=0; (j<nral); j++) 
	fprintf(fp,"  %5s",*atomname[plist[ftp].param[i].a[j]]);
      for(j=0; (j<nrfp); j++) 
	fprintf(fp,"  %10.3e",plist[ftp].param[i].c[j]);
      fprintf(fp,"\n");
    }
  }
}

void print_rtp(char *filenm,char *title,t_atoms *atoms,
	       t_params plist[],int cgnr[],int nbts,int bts[])
{
  FILE *fp;
  int i;
  
  fp = ffopen(filenm,"w");
  fprintf(fp,"; %s\n",title);
  fprintf(fp,"\n");
  fprintf(fp,"[ %s ]\n",*atoms->resname[0]);
  fprintf(fp,"\n");
  fprintf(fp,"[ atoms ]\n");
  for(i=0; (i<atoms->nr); i++) {
    fprintf(fp,"%-8s  %12s  %8.4f  %5d\n",
	    *atoms->atomname[i],*atoms->atomtype[i],
	    atoms->atom[i].q,cgnr[i]);
  }
  print_pl(fp,plist,F_BONDS,"bonds",atoms->atomname);
  print_pl(fp,plist,F_ANGLES,"angles",atoms->atomname);
  print_pl(fp,plist,F_PDIHS,"dihedrals",atoms->atomname);
  print_pl(fp,plist,F_IDIHS,"impropers",atoms->atomname);
  
  fclose(fp);
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
    d = pa->AK - pb->AK;
  if (d == 0) 
    d = pa->AL - pb->AL;
  /*if (d == 0)
    return strlen(pb->s) - strlen(pa->s);
    else*/
    return d;
}

static int acomp(const void *a,const void *b)
{
  atom_id *aa = (atom_id *)a;
  atom_id *ab = (atom_id *)b;
  
  return (*aa - *ab);
}

void clean_excls(int nr,t_excls excls[])
{
  int i,j,k;
  
  for(i=0; (i<nr); i++) {
    if ( excls[i].nr > 0) {
      qsort(excls[i].e,excls[i].nr,sizeof(excls[i].e[0]),acomp);
      k=0;
      for(j=0; (j<excls[i].nr); j++) {
	if (excls[i].e[j] != excls[i].e[k]) {
	  excls[i].e[++k] = excls[i].e[j];
	}
      }
      excls[i].nr = ++k;
    }
  }
}

static void clean_thole(t_params *ps)
{
  int     i,j;
  atom_id a,ai,aj,ak,al;
  
  if (ps->nr > 0) {
    /* swap atomnumbers in bond if first larger than second: */
    for(i=0; (i<ps->nr); i++)
      if ( ps->param[i].AK < ps->param[i].AI ) {
	a = ps->param[i].AI;
	ps->param[i].AI = ps->param[i].AK;
	ps->param[i].AK = a;
	a = ps->param[i].AJ;
	ps->param[i].AJ = ps->param[i].AL;
	ps->param[i].AL = a;
      }
    
    /* Sort bonds */
    qsort(ps->param,ps->nr,(size_t)sizeof(ps->param[0]),pcompar);
    
    /* remove doubles, keep the first one always. */
    j = 1;
    for(i=1; (i<ps->nr); i++) {
      if ((ps->param[i].AI != ps->param[j-1].AI) ||
	  (ps->param[i].AJ != ps->param[j-1].AJ) ||
	  (ps->param[i].AK != ps->param[j-1].AK) ||
	  (ps->param[i].AL != ps->param[j-1].AL) ) {
	cp_param(&(ps->param[j]),&(ps->param[i]));
	j++;
      } 
    }
    fprintf(stderr,"Number of Tholes was %d, now %d\n",ps->nr,j);
    ps->nr=j;
  }
  else
    fprintf(stderr,"No Tholes\n");
}

void delete_shell_interactions(t_params plist[F_NRE],t_atoms *atoms,
			       t_atomtype atype,t_nextnb *nnb,
			       t_excls excls[])
{
  int atp,jtp,jid,i,j,k,l,m,ftype,nb,nra,npol=0;
  bool *bRemove,*bHaveShell,bShell;
  int  *shell_index;
  t_param *p;
  int bt[] = { F_BONDS, F_ANGLES, F_PDIHS, F_IDIHS, F_LJ14 };
  
  pr_alloc(plist[F_BONDS].nr,&(plist[F_POLARIZATION]));
  for(i=0; (i<asize(bt)); i++) {
    ftype = bt[i];
    p     = plist[ftype].param;
    nra   = interaction_function[ftype].nratoms;
    snew(bRemove,plist[ftype].nr);
    for(j=0; (j<plist[ftype].nr); j++) {
      for(k=0; (k<nra); k++) {
	atp = atoms->atom[p[j].a[k]].type;
	if (get_atomtype_ptype(atp,atype) == eptShell) {
	  bRemove[j] = TRUE;
	  if (ftype == F_BONDS) {
	    memcpy(&plist[F_POLARIZATION].param[npol],
		   &plist[F_BONDS].param[j],
		   sizeof(plist[F_BONDS].param[j]));
	    plist[F_POLARIZATION].param[npol].C0 = atoms->atom[p[j].a[k]].qB;
	    npol++;
	    fprintf(stderr,"Adding polarization\n");
	  }
	}
      }
    }
    for(j=k=0; (j<plist[ftype].nr); j++) {
      if (!bRemove[j]) 
	memcpy(&plist[ftype].param[k++],
	       &plist[ftype].param[j],
	       sizeof(plist[ftype].param[j]));
    }
    plist[ftype].nr = k;
    sfree(bRemove);
  }
  plist[F_POLARIZATION].nr = npol;

  /* now for all atoms */
  for (i=0; (i < atoms->nr); i++) {
    atp = atoms->atom[i].type;
    if (get_atomtype_ptype(atp,atype) == eptShell) {
      for(m=3; (m<=4); m++) {
	/* for all fifth bonded atoms of atom i */
	for (j=0; (j < nnb->nrexcl[i][m]); j++) {
      
	  /* store the 1st neighbour in nb */
	  nb = nnb->a[i][m][j];
	  jtp = atoms->atom[nb].type;
	  if ((i != nb) && (strcasecmp(get_atomtype_name(jtp,atype),"SHELL") == 0)) {
	    srenew(excls[i].e,excls[i].nr+1);
	    excls[i].e[excls[i].nr++] = nb;
	    fprintf(stderr,"Excluding %d from %d\n",nb+1,i+1);
	  }
	}
      }
    }
  }
  clean_excls(atoms->nr,excls);
  pr_alloc(atoms->nr,&(plist[F_THOLE_POL]));
  npol = 0;
  snew(bHaveShell,atoms->nr);
  snew(shell_index,atoms->nr);
  for (i=0; (i < atoms->nr); i++) {
    /* for all first bonded atoms of atom i */
    for (j=0; (j < nnb->nrexcl[i][1]); j++) {
      jid = nnb->a[i][1][j];
      atp = atoms->atom[jid].type;
    
      if (get_atomtype_ptype(atp,atype) == eptShell) {
	bHaveShell[i] = TRUE;
	shell_index[i] = jid;
      }
    }
  }
  
  for (i=0; (i < atoms->nr); i++) {
    if (bHaveShell[i]) {
      /* for all first bonded atoms of atom i */
      for (j=0; (j < nnb->nrexcl[i][1]); j++) {
	jid = nnb->a[i][1][j];
	if (bHaveShell[jid]) {
	  plist[F_THOLE_POL].param[npol].AI = i;
	  plist[F_THOLE_POL].param[npol].AJ = shell_index[i];
	  plist[F_THOLE_POL].param[npol].AK = jid;
	  plist[F_THOLE_POL].param[npol].AL = shell_index[jid];
	  plist[F_THOLE_POL].param[npol].C0 = 2.6;
	  plist[F_THOLE_POL].param[npol].C1 = atoms->atom[shell_index[i]].qB;
	  plist[F_THOLE_POL].param[npol].C2 = atoms->atom[shell_index[jid]].qB;
	  npol++;
	}
      }
      /* for all second bonded atoms of atom i */
      for (j=0; (j < nnb->nrexcl[i][2]); j++) {
	jid = nnb->a[i][2][j];
	if ((jid != i) && bHaveShell[jid]) {
	  plist[F_THOLE_POL].param[npol].AI = i;
	  plist[F_THOLE_POL].param[npol].AJ = shell_index[i];
	  plist[F_THOLE_POL].param[npol].AK = jid;
	  plist[F_THOLE_POL].param[npol].AL = shell_index[jid];
	  plist[F_THOLE_POL].param[npol].C0 = 2.6;
	  plist[F_THOLE_POL].param[npol].C1 = atoms->atom[shell_index[i]].qB;
	  plist[F_THOLE_POL].param[npol].C2 = atoms->atom[shell_index[jid]].qB;
	  npol++;
	}
      }
    }
  }
  plist[F_THOLE_POL].nr = npol;
  clean_thole(&plist[F_THOLE_POL]);
  /* Add shell interactions to pairs */
  for(j=0; (j<plist[F_LJ14].nr); j++) {
    atom_id ai,aj;
    ai = plist[F_LJ14].param[j].AI;
    aj = plist[F_LJ14].param[j].AJ;
    if ((bHaveShell[ai]) || (bHaveShell[aj])) {
    }
  }
}

real calc_dip(t_atoms *atoms,rvec x[])
{
  int i;
  rvec mu,mm;
  
  clear_rvec(mu);
  for(i=0; (i<atoms->nr); i++) {
    svmul(atoms->atom[i].q,x[i],mm);
    rvec_inc(mu,mm);
  }
  return norm(mu)*ENM2DEBYE;
}
  
void reset_q(t_atoms *atoms)
{
  int i;
  
  /* Use values from file */
  for(i=0; (i<atoms->nr); i++) 
    atoms->atom[i].qB = atoms->atom[i].q;
}

void add_shells(x2top_nm2t nm2t,t_atoms **atoms,
		t_atomtype atype,t_params plist[],
		rvec **x,t_symtab *symtab,t_excls **excls)
{
  int     i,j,k,iat,shell,atp,ns=0;
  int     *renum;
  char    buf[32];
  t_param p;
  t_atom  *shell_atom;
  t_atoms *newa;
  t_excls *newexcls;
  rvec    *newx;
  
  snew(shell_atom,1);
  memset(&p,0,sizeof(p));
  snew(renum,(*atoms)->nr*2+2);
  for(i=0; (i<(*atoms)->nr); i++) {
    atp = (*atoms)->atom[i].type;
    renum[i] = i+ns;
    if (get_atomtype_qB(atp,atype) != 0) {
      ns++;
      p.AI = renum[i];
      p.AJ = renum[i]+1;
      p.C0 = get_atomtype_qB(atp,atype);
      push_bondnow(&(plist[F_POLARIZATION]),&p);
    }
  }
  renum[(*atoms)->nr] = (*atoms)->nr + ns;
  shell_atom->ptype = eptShell;
  shell = add_atomtype(atype,symtab,shell_atom,"SHELL",&p,
		       0,0,0,0,0);
  
  if (ns > 0) {
    snew(newa,1);
    init_t_atoms(newa,(*atoms)->nr+ns,TRUE);
    snew(newexcls,newa->nr);
    newa->nres = (*atoms)->nres;
    snew(newx,newa->nr);

    for(i=0; (i<(*atoms)->nr); i++) {
      newa->atom[renum[i]]     = (*atoms)->atom[i];
      newa->atomname[renum[i]] = put_symtab(symtab,*(*atoms)->atomname[i]);
      copy_rvec((*x)[i],newx[renum[i]]);
    }
    for(i=0; (i<(*atoms)->nres); i++) {
      newa->resname[i] = put_symtab(symtab,*(*atoms)->resname[i]);
    }
    
    for(i=0; (i<(*atoms)->nr); i++) {
      iat = renum[i];
      for(k=0; (k<(*excls)[i].nr); k++)
	newexcls[iat].e[k] = (*excls)[i].e[k];
      for(j=iat+1; (j<renum[i+1]); j++) {
	newa->atom[j]       = (*atoms)->atom[i];
	newa->atom[iat].q   = 0;
	newa->atom[iat].qB  = 0;
	newa->atom[j].m     = 0;
	newa->atom[j].mB    = 0;
	newa->atom[j].type  = shell;
	newa->atom[j].typeB = shell;
	newa->atom[j].resnr = (*atoms)->atom[i].resnr;
	sprintf(buf,"Sh%s",*((*atoms)->atomname[i]));
	newa->atomname[j] = put_symtab(symtab,buf);
	copy_rvec((*x)[i],newx[j]);
	for(k=0; (k<(*excls)[i].nr); k++)
	  newexcls[j].e[k] = (*excls)[i].e[k];      }
    }
    *atoms = newa;
    sfree(*x);
    *x = newx;
    *excls = newexcls;
  }
  for(i=0; (i<F_NRE); i++) {
    if (i != F_POLARIZATION)
      for(j=0; (j<plist[i].nr); j++) 
	for(k=0; (k<NRAL(i)); k++) 
	  plist[i].param[j].a[k] = renum[plist[i].param[j].a[k]];
  }
  sfree(renum);
  sfree(shell_atom);
}

static void lo_symmetrize_charges(t_atoms *atoms,t_atomtype atype,
				  t_params *bonds,int atomnumber[4])
{
  int i,j,nh,ai,aj;
  int hs[4];
  double qaver;
    
  for(i=0; (i<atoms->nr); i++) {
    nh = 0;
    if (get_atomtype_atomnumber(atoms->atom[i].type,atype) == 
	atomnumber[0]) {
      hs[nh++] = i;
      for(j=0; ((j<bonds->nr) && (nh < 4)); j++) {
	ai = bonds->param[j].AI;
	aj = bonds->param[j].AJ;
	if ((ai == i) && 
	    (get_atomtype_atomnumber(atoms->atom[aj].type,atype) == atomnumber[nh])) {
	  hs[nh++] = aj;
	}
	else if ((aj == i) && 
		 (get_atomtype_atomnumber(atoms->atom[ai].type,atype) == atomnumber[nh])) {
	  hs[nh++] = ai; 
	}
      }
      if (nh == 4) {
	qaver = (atoms->atom[hs[1]].q + atoms->atom[hs[2]].q +
		 atoms->atom[hs[3]].q)/3.0;
	atoms->atom[hs[1]].q = atoms->atom[hs[2]].q =
	  atoms->atom[hs[3]].q = qaver;
      }
    }
  }
}

void symmetrize_charges(t_atoms *atoms,t_atomtype atype,
			t_params *bonds)
{
  char **strings = NULL;
  int i,nstrings;
  int at[4];
  
  nstrings = get_file("symmetric-charges.dat",&strings);
  for(i=0; (i<nstrings); i++) {
    if (sscanf(strings[i],"%d%d%d%d",
	       &at[0],&at[1],&at[2],&at[3]) == 4)
      lo_symmetrize_charges(atoms,atype,bonds,at);
    
    sfree(strings[i]);
  }
  sfree(strings);
}
