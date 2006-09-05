#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "typedefs.h"
#include "smalloc.h"
#include "domdec.h"
#include "names.h"


typedef struct gmx_reverse_top {
  int *index; /* Index for each atom into il    */
  int *il;    /* ftype|type|a0|...|an|ftype|... */
} gmx_reverse_top_t;

#define EXCLS_ALLOC_SIZE  1000
#define IATOM_ALLOC_SIZE  1000

/* Code that only works for one moleculetype consisting of one charge group */
/* #define ONE_MOLTYPE_ONE_CG */

#ifdef ONE_MOLTYPE_ONE_CG
static int natoms_global;
#endif

static int count_excls(t_block *cgs,t_block *excls,int *n_intercg_excl)
{
  int n,n_inter,cg,at0,at1,at,excl,atj;
  
  n = 0;
  *n_intercg_excl = 0;
  for(cg=0; cg<cgs->nr; cg++) {
    at0 = cgs->index[cg];
    at1 = cgs->index[cg+1];
    for(at=at0; at<at1; at++) {
      for(excl=excls->index[at]; excl<excls->index[at+1]; excl++) {
	atj = excls->a[excl];
	if (atj > at) {
	  n++;
	  if (atj < at0 || atj >= at1)
	    (*n_intercg_excl)++;
	}
      }
    }
  }  
  
  return n;
}

static int low_make_reverse_top(t_idef *idef,int *count,gmx_reverse_top_t *rt,
				bool bAssign)
{
  int ftype,nral,i,ind_at,j;
  t_ilist *il;
  t_iatom *ia;
  atom_id a;
  int nint;
  
  nint = 0;
  for(ftype=0; ftype<F_NRE; ftype++) {
    if ((interaction_function[ftype].flags & (IF_BOND | IF_VSITE))
	|| ftype == F_SETTLE) {
      nral = NRAL(ftype);
      il = &idef->il[ftype];
      ia  = il->iatoms;
      if ((interaction_function[ftype].flags & IF_BOND) && nral > 2) {
	/* Couple to the second atom in the interaction */
	ind_at = 1;
      } else {
	/* Couple to the first atom in the interaction */
	ind_at = 0;
      }
      for(i=0; i<il->nr; i+=1+nral) {
	ia = il->iatoms + i;
	a = ia[1+ind_at];
	if (bAssign) {
	  rt->il[rt->index[a]+count[a]] = ftype;
	  for(j=0; j<1+nral; j++)
	    rt->il[rt->index[a]+count[a]+1+j] = ia[j];
	}
	count[a] += 2 + nral;
	nint++;
      }
    }
  }

  return nint;
}

static gmx_reverse_top_t *make_reverse_top(int natoms,t_idef *idef,int *nint)
{
  int *count,i;
  gmx_reverse_top_t *rt;

  snew(rt,1);

  /* Count the interactions */
  snew(count,natoms);
  low_make_reverse_top(idef,count,rt,FALSE);

  snew(rt->index,natoms+1);
  rt->index[0] = 0;
  for(i=0; i<natoms; i++) {
    rt->index[i+1] = rt->index[i] + count[i];
    count[i] = 0;
  }
  snew(rt->il,rt->index[natoms]);

  /* Store the interactions */
  *nint = low_make_reverse_top(idef,count,rt,TRUE);

  sfree(count);

  return rt;
}

void dd_make_reverse_top(FILE *fplog,
			 gmx_domdec_t *dd,t_topology *top,
			 bool bDynamics,int eeltype)
{
  int natoms,nexcl,a;
  
  fprintf(fplog,"\nLinking all bonded interactions to atoms\n");

  natoms = top->atoms.nr;

  dd->reverse_top = make_reverse_top(natoms,&top->idef,&dd->nbonded_global);
  
  nexcl = count_excls(&top->blocks[ebCGS],&top->blocks[ebEXCLS],
		      &dd->n_intercg_excl);
  if (EEL_FULL(eeltype)) {
    dd->nbonded_global += nexcl;
    if (dd->n_intercg_excl)
      fprintf(fplog,"There are %d inter charge-group exclusions,\n"
	      "will use an extra communication step for exclusion forces for %s\n",
	      dd->n_intercg_excl,eel_names[eeltype]);
  }
  
  snew(dd->ga2la,natoms);
  for(a=0; a<natoms; a++)
    dd->ga2la[a].cell = -1;
  
  if (top->idef.il[F_CONSTR].nr > 0)
    init_domdec_constraints(dd,natoms,&top->idef,bDynamics);
}

static void make_local_ilist(gmx_domdec_t *dd,t_functype ftype,
			     t_ilist *il,t_ilist *lil)
{
  int nral,nlocal,i,j,n;
  t_iatom *ia,*lia;
  atom_id a;
  
  nral = NRAL(ftype);

  if (lil->iatoms == NULL)
    /* In most cases we could do with far less memory */
    snew(lil->iatoms,il->nr);

  nlocal = dd->nat_home;

  n   = 0;
  ia  = il->iatoms;
  lia = lil->iatoms;
  for(i=0; i<il->nr; i+=1+nral) {
    ia = il->iatoms + i;
    a = dd->ga2la[ia[1]].a;
    if (a >= 0 && a < nlocal) {
      lia[0] = ia[0];
      lia[1] = a;
      for(j=2; j<=nral; j++) {
	a = dd->ga2la[ia[j]].a;
	if (a == -1)
	  /* This interaction could be assigned to a different node
	   * with the current setup (as long as the distance between
	   * all particles involved is smaller than the grid spacing).
	   * But we need an algorithm to assign it.
	   */
	  gmx_fatal(FARGS,"A bonded interaction is spread over cells");
	lia[j] = a;
      }
      lia += 1 + nral;
      n++;
    }
    ia += 1 + nral;
  }

  lil->nr = n*(1 + nral);
}

static int make_local_bondeds(gmx_domdec_t *dd,t_idef *idef)
{
  int nbonded_local,i,gat,j,ftype,nral,d,k,kc;
  int *index,*rtil;
  t_iatom *iatoms,tiatoms[1+MAXATOMLIST],*liatoms;
  bool bUse;
  ivec shift_min;
  gmx_ga2la_t *ga2la;
  t_ilist *il;

  index = dd->reverse_top->index;
  rtil  = dd->reverse_top->il;

  /* Clear the counts */
  for(ftype=0; ftype<F_NRE; ftype++)
    idef->il[ftype].nr = 0;
  nbonded_local = 0;
  
  for(i=0; i<dd->nat_tot; i++) {
    /* Get the global atom number */
    gat = dd->gatindex[i];
    /* Check all interactions assigned to this atom */
    j = index[gat];
    while (j < index[gat+1]) {
      ftype  = rtil[j++];
      iatoms = rtil + j;
      nral = NRAL(ftype);
      bUse = TRUE;
      for(d=0; d<DIM; d++)
	shift_min[d] = 1;
      for(k=1; k<=nral && bUse; k++) {
	ga2la = &dd->ga2la[iatoms[k]];
	kc = ga2la->cell;
	if (kc == -1) {
	  /* We do not have this atom of this interaction locally */
	  bUse = FALSE;
	} else {
	  tiatoms[k] = ga2la->a;
	  for(d=0; d<DIM; d++)
	    if (dd->shift[kc][d] < shift_min[d])
	      shift_min[d] = dd->shift[kc][d];
	}
      }
      if (bUse &&
	  shift_min[XX]==0 && shift_min[YY]==0 && shift_min[ZZ]==0) {
	/* Add this interaction to the local topology */
	il = &idef->il[ftype];
	if (il->nr >= il->nalloc) {
	  il->nalloc += IATOM_ALLOC_SIZE*(1+nral);
	  srenew(il->iatoms,il->nalloc);
	}
	liatoms = il->iatoms + il->nr;
	liatoms[0] = iatoms[0];
	for(k=1; k<=nral; k++)
	  liatoms[k] = tiatoms[k];
	/* Sum locally so we can check in global_stat if we have everything */
	nbonded_local++;
	il->nr += 1+nral;
      }
      j += 1 + nral;
    }
  }

  return nbonded_local;
}

static int make_local_exclusions(gmx_domdec_t *dd,t_forcerec *fr,
				 t_block *excls,t_block *lexcls)
{
  int nicell,n,count,ic,jla0,jla1,jla,cg,la0,la1,la,a,i,j;
  gmx_ga2la_t *ga2la;

  /* Since for RF and PME we need to loop over the exclusions
   * we should store each exclusion only once. This is done
   * using the same cell scheme as used for neighbor searching.
   * The exclusions involving non-home atoms are stored only
   * one way: atom j is in the excl list of i only for j > i,
   * where i and j are local atom numbers.
   */

  lexcls->nr = dd->cgindex[dd->icell[dd->nicell-1].cg1];
  if (lexcls->nr+1 > lexcls->nalloc_index) {
    lexcls->nalloc_index = over_alloc(lexcls->nr)+1;
    srenew(lexcls->index,lexcls->nalloc_index);
  }

  if (dd->n_intercg_excl)
    nicell = dd->nicell;
  else
    nicell = 1;
  n = 0;
  count = 0;
  for(ic=0; ic<nicell; ic++) {
    jla0 = dd->cgindex[dd->icell[ic].jcg0];
    jla1 = dd->cgindex[dd->icell[ic].jcg1];
    for(cg=dd->ncg_cell[ic]; cg<dd->ncg_cell[ic+1]; cg++) {
      /* Here we assume the number of exclusions in one charge group
       * is never larger than EXCLS_ALLOC_SIZE
       */
      if (n+EXCLS_ALLOC_SIZE > lexcls->nalloc_a) {
	lexcls->nalloc_a += EXCLS_ALLOC_SIZE;
	srenew(lexcls->a,lexcls->nalloc_a);
      }
      la0 = dd->cgindex[cg];
      switch (fr->solvent_type[cg]) {
      case esolNO:
	/* Copy the exclusions from the global top */
	la1 = dd->cgindex[cg+1];
	for(la=la0; la<la1; la++) {
	  lexcls->index[la] = n;
	  a = dd->gatindex[la];
	  for(i=excls->index[a]; i<excls->index[a+1]; i++) {
	    ga2la = &dd->ga2la[excls->a[i]];
	    if (ga2la->cell != -1) {
	      jla = ga2la->a;
	      /* Check to avoid double counts */
	      if (jla >= jla0 && jla < jla1) {
		lexcls->a[n++] = jla;
		if (jla > la)
		  count++;
	      }
	    }
	  }
	}
	break;
      case esolSPC:
	/* Exclude all 9 atom pairs */
	for(la=la0; la<la0+3; la++) {
	  lexcls->index[la] = n;
	  if (ic == 0) {
	    for(j=la0; j<la0+3; j++)
	      lexcls->a[n++] = j;
	  }
	}
	if (ic == 0)
	  count += 3;
	break;
      case esolTIP4P:
	/* Exclude all 16 atoms pairs */
	for(la=la0; la<la0+4; la++) {
	  lexcls->index[la] = n;
	  if (ic == 0) {
	    for(j=la0; j<la0+4; j++)
	      lexcls->a[n++] = j;
	  }
	}
	if (ic == 0)
	  count += 6;
	break;
      default:
	gmx_fatal(FARGS,"Unknown solvent type %d",fr->solvent_type[cg]);
      }
    }
  }
  if (dd->n_intercg_excl == 0) {
    /* There are no exclusions involving non-home charge groups,
     * but we need to set the indices for neighborsearching.
     */
    la0 = dd->cgindex[dd->icell[0].cg1];
    for(la=la0; la<lexcls->nr; la++)
      lexcls->index[la] = n;
  }
  lexcls->index[lexcls->nr] = n;
  lexcls->nra = n;
  if (dd->n_intercg_excl == 0) {
    /* nr is only used to loop over the exclusions for Ewald and RF,
     * so we can set it to the number of home atoms for efficiency.
     */
    lexcls->nr = dd->cgindex[dd->icell[0].cg1];
  }
  if (debug)
    fprintf(debug,"We have %d exclusions\n",lexcls->nra);

  return count;
}

#ifdef ONE_MOLTYPE_ONE_CG
static void make_local_ilist_onemoltype_onecg(gmx_domdec_t *dd,
					      t_functype ftype,
					      t_ilist *il,t_ilist *lil)
{
  int nral,nhome,i,j,n;
  t_iatom *ia,*lia;
  long long int maxlli;
  int maxi;

  nral = NRAL(ftype);

  if (lil->iatoms == NULL)
    /* In most cases we could do with far less memory */
    snew(lil->iatoms,il->nr);

  nhome = dd->comm1[0].nat;

  maxlli = (long long int)il->nr*nhome/natoms_global;
  maxi   = maxlli;

  n   = 0;
  ia  = il->iatoms;
  lia = lil->iatoms;
  for(i=0; i<maxi; i+=1+nral) {
    ia = il->iatoms + i;
    lia[0] = ia[0];
    lia[1] = ia[1];
    for(j=2; j<=nral; j++) {
      lia[j] = ia[j];
    }
    lia += 1 + nral;
    n++;
    ia += 1 + nral;
  }
  
  lil->nr = n*(1 + nral);
}

static void make_local_idef(gmx_domdec_t *dd,t_idef *idef,t_idef *lidef)
{
  int f;

  lidef->ntypes   = idef->ntypes;
  lidef->nodeid   = idef->nodeid;
  lidef->atnr     = idef->atnr;
  lidef->functype = idef->functype;
  lidef->iparams  = idef->iparams;

  for(f=0; f<F_NRE; f++)
    /* make_local_ilist(dd,f,&idef->il[f],&lidef->il[f]); */
    make_local_ilist_onemoltype_onecg(dd,f,&idef->il[f],&lidef->il[f]);
}
#endif

void make_local_cgs(gmx_domdec_t *dd,t_block *lcgs)
{
  lcgs->nr    = dd->ncg_tot;
  lcgs->index = dd->cgindex;
  lcgs->nra   = dd->nat_tot;
}

void make_local_top(FILE *fplog,gmx_domdec_t *dd,
		    t_forcerec *fr,t_topology *top,t_topology *ltop)
{
  int nexcl;

  if (debug)
    fprintf(debug,"Making local topology\n");

  ltop->name  = top->name;
  make_local_cgs(dd,&ltop->blocks[ebCGS]);

#ifdef ONE_MOLTYPE_ONE_CG
  natoms_global = top->blocks[ebCGS].nra;

  make_local_idef(dd,&top->idef,&ltop->idef);
#else
  dd->nbonded_local = make_local_bondeds(dd,&ltop->idef);
#endif

  nexcl = make_local_exclusions(dd,fr,&top->blocks[ebEXCLS],
				&ltop->blocks[ebEXCLS]);
  /* With cut-off electrostatics we don't care that exclusions
   * beyond the cut-off can be missing.
   */
  if (EEL_FULL(fr->eeltype))
    dd->nbonded_local += nexcl;
 
  ltop->atoms = top->atoms;
  ltop->atomtypes = top->atomtypes;
  /*
  for(eb=0; eb<ebNR; eb++) {
    if (eb != ebCGS)
      ltop->blocks[eb] = top->blocks[eb];
  }
  */
  ltop->symtab = top->symtab;
}
