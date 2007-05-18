#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "typedefs.h"
#include "smalloc.h"
#include "domdec.h"
#include "names.h"
#include "network.h"


typedef struct gmx_reverse_top {
  int *index; /* Index for each atom into il    */
  int *il;    /* ftype|type|a0|...|an|ftype|... */
  int n_intercg_vsite;
} gmx_reverse_top_t;

#define EXCLS_ALLOC_SIZE  1000
#define IATOM_ALLOC_SIZE  1000

/* Code that only works for one moleculetype consisting of one charge group */
/* #define ONE_MOLTYPE_ONE_CG */

#ifdef ONE_MOLTYPE_ONE_CG
static int natoms_global;
#endif

/* Static pointers only used for an error message */
static t_topology *err_top_global,*err_top_local;

void dd_print_missing_interactions(FILE *fplog,t_commrec *cr,int local_count)
{
  int  ndiff_tot,sign,cl[F_NRE],n,ndiff,rest_global,rest_local;
  int  ftype,nral;
  char *ts,buf[STRLEN];
  gmx_domdec_t *dd;

  dd = cr->dd;

  fprintf(fplog,"\nNot all bonded interactions have been properly assigned to the domain decomposition cells\n");
  fflush(fplog);

  ndiff_tot = local_count - dd->nbonded_global;
  if (ndiff_tot >= 0) {
    sign = 1;
    ts = "more than once evaluated";
  } else {
    sign = -1;
    ts = "missing";
  }

  for(ftype=0; ftype<F_NRE; ftype++) {
    nral = NRAL(ftype);
    cl[ftype] = err_top_local->idef.il[ftype].nr/(1+nral);
  }

  gmx_sumi(F_NRE,cl,cr);
  
  if (DDMASTER(dd)) {
    fprintf(fplog,"\nA list of %s interactions:\n",ts);
    fprintf(stderr,"\nA list of %s interactions:\n",ts);
    rest_global = dd->nbonded_global;
    rest_local  = local_count;
    for(ftype=0; ftype<F_NRE; ftype++) {
      if ((interaction_function[ftype].flags & IF_BOND) || ftype == F_SETTLE) {
	nral = NRAL(ftype);
	n = err_top_global->idef.il[ftype].nr/(1+nral);
	ndiff = cl[ftype] - n;
	if (ndiff != 0) {
	  sprintf(buf,"%20s of %6d %s %6d",
		  interaction_function[ftype].longname,n,ts,sign*ndiff);
	  fprintf(fplog,"%s\n",buf);
	  fprintf(stderr,"%s\n",buf);
	}
	rest_global -= n;
	rest_local  -= cl[ftype];
      }
    }
    
    ndiff = rest_local - rest_global;
    if (ndiff != 0) {
      sprintf(buf,"%20s of %6d %s %6d","exclusions",rest_global,ts,sign*ndiff);
      fprintf(fplog,"%s\n",buf);
      fprintf(stderr,"%s\n",buf);
    }

    if (sign == 1) {
      gmx_fatal(FARGS,"%d of the %d bonded interactions were evaluated multiple times. This can occur when the number of domain decomposition cell in a direction is 2 and a cell is not much larger than the cut-off length. The solution is to use 1 or 3 or more cells in such a direction.",
		sign*ndiff_tot,cr->dd->nbonded_global);
    } else {
      gmx_fatal(FARGS,"%d of the %d bonded interactions could not be calculated because some atoms involved moved further apart than the cut-off distance",
		sign*ndiff_tot,cr->dd->nbonded_global);
    }
  }
}

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

static int low_make_reverse_top(t_idef *idef,t_atom *atom,int **vsite_pbc,
				int *count,gmx_reverse_top_t *rt,
				bool bAssign)
{
  int  ftype,nral,i,ind_at,j;
  t_ilist *il;
  t_iatom *ia;
  atom_id a;
  int  nint;

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
	if (interaction_function[ftype].flags & IF_VSITE) {
	  if (bAssign) {
	    /* Add an entry to iatoms for storing 
	     * which of the constructing atoms are vsites again.
	     */
	    rt->il[rt->index[a]+count[a]+2+nral] = 0;
	    for(j=2; j<1+nral; j++) {
	      if (atom[ia[j]].ptype == eptVSite)
		rt->il[rt->index[a]+count[a]+2+nral] |= (2<<j);
	    }
	    /* Store the vsite pbc atom in a second extra entry */
	    rt->il[rt->index[a]+count[a]+2+nral+1] =
	      (vsite_pbc ? vsite_pbc[ftype-F_VSITE2][i/(1+nral)] : -2);
	  }
	  count[a] += 2 + nral + 2;
	} else {
	  /* We do not count vsites since they are always uniquely assigned
	   * and can be assigned to multiple nodes with recursive vsites.
	   */
	  nint++;
	  count[a] += 2 + nral;
	}
      }
    }
  }

  return nint;
}

static gmx_reverse_top_t *make_reverse_top(int natoms,t_idef *idef,
					   t_atom *atom,int **vsite_pbc,
					   int *nint)
{
  int *count,i;
  gmx_reverse_top_t *rt;

  snew(rt,1);

  /* Count the interactions */
  snew(count,natoms);
  low_make_reverse_top(idef,atom,vsite_pbc,count,rt,FALSE);

  snew(rt->index,natoms+1);
  rt->index[0] = 0;
  for(i=0; i<natoms; i++) {
    rt->index[i+1] = rt->index[i] + count[i];
    count[i] = 0;
  }
  snew(rt->il,rt->index[natoms]);

  /* Store the interactions */
  *nint = low_make_reverse_top(idef,atom,vsite_pbc,count,rt,TRUE);

  sfree(count);

  return rt;
}

void dd_make_reverse_top(FILE *fplog,
			 gmx_domdec_t *dd,t_topology *top,
			 gmx_vsite_t *vsite,
			 bool bDynamics,int eeltype)
{
  int natoms,n_recursive_vsite,nexcl,a;

  fprintf(fplog,"\nLinking all bonded interactions to atoms\n");

  natoms = top->atoms.nr;

  dd->reverse_top = make_reverse_top(natoms,&top->idef,top->atoms.atom,
				     vsite ? vsite->vsite_pbc : NULL,
				     &dd->nbonded_global);

  if (vsite)
    dd->reverse_top->n_intercg_vsite = vsite->n_intercg_vsite;
  else
    dd->reverse_top->n_intercg_vsite = 0;
  
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
  
  if (dd->reverse_top->n_intercg_vsite > 0) {
    fprintf(fplog,"There are %d inter charge-group virtual sites,\n"
	    "will an extra communication step for selected coordinates and forces\n",
	    dd->reverse_top->n_intercg_vsite);
    init_domdec_vsites(dd,natoms);
  }

  if (top->idef.il[F_CONSTR].nr > 0) {
    init_domdec_constraints(dd,natoms,&top->idef,&top->blocks[ebCGS],
			    bDynamics);
    if (dd->constraint_comm)
      fprintf(fplog,"There are inter charge-group constraints,\n"
	      "will communicate selected coordinates each lincs iteration\n");
  }
}

static inline void add_ifunc(int type,int nral,t_iatom *tiatoms,t_ilist *il)
{
  t_iatom *liatoms;
  int     k;

  if (il->nr >= il->nalloc) {
    il->nalloc += IATOM_ALLOC_SIZE*(1 + nral);
    srenew(il->iatoms,il->nalloc);
  }
  liatoms = il->iatoms + il->nr;
  liatoms[0] = type;
  for(k=1; k<=nral; k++)
    liatoms[k] = tiatoms[k];
  il->nr += 1 + nral;
}

static void add_vsite(gmx_domdec_t *dd,
		      int ftype,int nral,int i,t_iatom *iatoms,
		      t_idef *idef,int **vsite_pbc,int *vsite_pbc_nalloc)
{
  int  k,vsi,pbc_ga;
  t_iatom tiatoms[1+MAXATOMLIST],*iatoms_r;
  gmx_ga2la_t *ga2la;
  int  *index,*rtil;
  int  j,ftype_r,nral_r;
  
  /* We know the local index of the first atom */
  tiatoms[1] = i;

  for(k=2; k<1+nral; k++) {
    ga2la = &dd->ga2la[iatoms[k]];
    if (ga2la->cell == 0) {
      tiatoms[k] = ga2la->a;
    } else {
      /* Copy the global index, convert later in make_local_vsites */
      tiatoms[k] = -(iatoms[k] + 1);
    }
  }

  /* Add this interaction to the local topology */
  add_ifunc(iatoms[0],nral,tiatoms,&idef->il[ftype]);
  if (vsite_pbc) {
    vsi = idef->il[ftype].nr/(1+nral) - 1;
    if (vsi >= vsite_pbc_nalloc[ftype-F_VSITE2]) {
      vsite_pbc_nalloc[ftype-F_VSITE2] += IATOM_ALLOC_SIZE;
      srenew(vsite_pbc[ftype-F_VSITE2],vsite_pbc_nalloc[ftype-F_VSITE2]);
    }
    if (i >= 0) {
      pbc_ga = iatoms[1+nral+1];
      if (pbc_ga < 0) {
	/* The pbc flag is one of the following two options:
	 * -2: vsite and all constructing atoms are within the same cg, no pbc
	 * -1: vsite and its first constructing atom are in the same cg, do pbc
	 */
	vsite_pbc[ftype-F_VSITE2][vsi] = pbc_ga;
      } else {
	/* Set the pbc atom for this vsite so we can make its pbc 
	 * identical to the rest of the atoms in its charge group.
	 * Since the order of the atoms does not change within a charge group,
	 * we do no need to access to global to local atom index.
	 */
	vsite_pbc[ftype-F_VSITE2][vsi] = i + pbc_ga - iatoms[1];
      }
    } else {
      /* This vsite is non-home (required for recursion),
       * and therefore there is no charge group to match pbc with.
       * But we always turn on full_pbc to assure that higher order
       * recursion works correctly.
       */
      vsite_pbc[ftype-F_VSITE2][vsi] = -1;
    }
  }
  
  if (iatoms[1+nral]) {
    /* Check for recursion */
    index = dd->reverse_top->index;
    rtil  = dd->reverse_top->il;
    for(k=2; k<1+nral; k++) {
      if ((iatoms[1+nral] & (2<<k)) && (tiatoms[k] < 0)) {
	/* This construction atoms is a vsite and not a home atom */
	if (gmx_debug_at)
	  fprintf(debug,"Constructing atom %d of vsite atom %d is a vsite and non-home\n",iatoms[k]+1,i >= 0 ? dd->gatindex[i]+1 : -i);

	/* Find the vsite construction */

	/* Check all interactions assigned to this atom */
	j = index[iatoms[k]];
	while (j < index[iatoms[k]+1]) {
	  ftype_r = rtil[j++];
	  nral_r = NRAL(ftype_r);
	  if (interaction_function[ftype_r].flags & IF_VSITE) {
	    /* Add this vsite (recursion).
	     * Signal that the vsite atom is non-home by negation.
	     */
	    add_vsite(dd,ftype_r,nral_r,-(iatoms[k]+1),rtil+j,idef,
		      vsite_pbc,vsite_pbc_nalloc);
	    j += 1 + nral_r + 2;
	  } else {
	    j += 1 + nral_r;
	  }
	}
      }
    }
  }
}

static int make_local_bondeds(gmx_domdec_t *dd,t_idef *idef,gmx_vsite_t *vsite)
{
  int nbonded_local,i,gat,j,ftype,nral,d,k,kc;
  int *index,*rtil,**vsite_pbc,*vsite_pbc_nalloc;
  t_iatom *iatoms,tiatoms[1+MAXATOMLIST];
  bool bUse;
  ivec shift_min;
  gmx_ga2la_t *ga2la;

  if (vsite && vsite->vsite_pbc) {
    vsite_pbc        = vsite->vsite_pbc_dd;
    vsite_pbc_nalloc = vsite->vsite_pbc_dd_nalloc;
  } else {
    vsite_pbc        = NULL;
    vsite_pbc_nalloc = NULL;
  }

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
      if (interaction_function[ftype].flags & IF_VSITE) {
	/* The vsite construction goes where the vsite itself is */
	if (i < dd->nat_home)
	  add_vsite(dd,ftype,nral,i,iatoms,idef,vsite_pbc,vsite_pbc_nalloc);
	j += 1 + nral + 2;
      } else {
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
	bUse = (bUse && 
		shift_min[XX]==0 && shift_min[YY]==0 && shift_min[ZZ]==0);
	if (bUse) {
	  /* Add this interaction to the local topology */
	  add_ifunc(iatoms[0],nral,tiatoms,&idef->il[ftype]);
	  /* Sum so we can check in global_stat if we have everything */
	  nbonded_local++;
	}
	j += 1 + nral;
      }
    }
  }

  if (dd->reverse_top->n_intercg_vsite) {
    dd->nat_tot_vsite = dd_make_local_vsites(dd,idef->il);
  } else {
    dd->nat_tot_vsite = dd->nat_tot;
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
      switch (GET_CGINFO_SOLOPT(fr->cginfo[cg])) {
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
	gmx_fatal(FARGS,"Unknown solvent type %d",
		  GET_CGINFO_SOLOPT(fr->cginfo[cg]));
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
    fprintf(debug,"We have %d exclusions, check count %d\n",
	    lexcls->nra,count);

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
    make_local_ilist_onemoltype_onecg(dd,f,&idef->il[f],&lidef->il[f]);
}
#endif

void dd_make_local_cgs(gmx_domdec_t *dd,t_block *lcgs)
{
  lcgs->nr    = dd->ncg_tot;
  lcgs->index = dd->cgindex;
  lcgs->nra   = dd->nat_tot;
}

void dd_make_local_top(FILE *fplog,gmx_domdec_t *dd,
		       t_forcerec *fr,gmx_vsite_t *vsite,
		       t_topology *top,t_topology *ltop)
{
  int nexcl;

  if (debug)
    fprintf(debug,"Making local topology\n");

  ltop->name  = top->name;
  dd_make_local_cgs(dd,&ltop->blocks[ebCGS]);

#ifdef ONE_MOLTYPE_ONE_CG
  natoms_global = top->blocks[ebCGS].nra;

  make_local_idef(dd,&top->idef,&ltop->idef);
#else
  dd->nbonded_local = make_local_bondeds(dd,&ltop->idef,vsite);
#endif

  nexcl = make_local_exclusions(dd,fr,&top->blocks[ebEXCLS],
				&ltop->blocks[ebEXCLS]);
  /* With cut-off electrostatics we don't care that exclusions
   * beyond the cut-off can be missing.
   */
  if (EEL_FULL(fr->eeltype))
    dd->nbonded_local += nexcl;
  
  ltop->atoms     = top->atoms;
  ltop->atomtypes = top->atomtypes;
  ltop->symtab    = top->symtab;

  /* For an error message only */
  err_top_global = top;
  err_top_local = ltop;
}

t_topology *dd_init_local_top(t_topology *top_global)
{
  t_topology *top;
  int i;
  
  snew(top,1);
  top->idef = top_global->idef;
  for(i=0; i<F_NRE; i++) {
    top->idef.il[i].iatoms = NULL;
    top->idef.il[i].nalloc = 0;
  }

  return top;
}
