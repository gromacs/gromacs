/*
 *       @(#) copyrgt.c 1.12 9/30/97
 *
 *       This source code is part of
 *
 *        G   R   O   M   A   C   S
 *
 * GROningen MAchine for Chemical Simulations
 *
 *            VERSION 2.0b
 * 
 * Copyright (c) 1990-1997,
 * BIOSON Research Institute, Dept. of Biophysical Chemistry,
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
 * GRowing Old MAkes el Chrono Sweat
 */

#ifdef HAVE_IDENT
#ident "@(#) toppush.c 1.72 9/30/97"
#endif

#include <math.h>

#include "sysstuff.h"
#include "smalloc.h"
#include "macros.h"
#include "assert.h"
#include "topdef.h"
#include "string2.h"
#include "names.h"
#include "toputil.h"
#include "toppush.h"
#include "topdirs.h"
#include "symtab.h"
#include "fatal.h"

static char  errbuf[256];

void generate_nbparams(int comb,int ftype,t_params *plist,t_atomtype *atype)
{
  int   i,j,k,nf;
  int   nr,nrfp;
  real  c,tmp6,tmp12,eps4,sig6;
  real  *c6,*c12;

  /* Lean mean shortcuts */
  nr   = atype->nr;
  nrfp = NRFP(ftype);
  snew(plist->param,nr*nr);
  plist->nr=nr*nr;
  
  /* Fill the matrix with force parameters */
  switch(ftype) {
  case F_LJ:
    switch (comb) {
    case COMB_GROMOS:
      /* Gromos rules */
      for(i=k=0; (i<nr); i++)
	for(j=0; (j<nr); j++,k++)
	  for(nf=0; (nf<nrfp); nf++) {
	    c = sqrt(atype->nb[i].c[nf] * atype->nb[j].c[nf]);
	    plist->param[k].c[nf]      = c;
	  }
      break;
    
    case COMB_EPSSIG:
      /* c0 and c1 are epsilon and sigma */
      
      /* Create temp storage arrays */
      snew(c6,nr);
      snew(c12,nr);

      /* Fill temp arrays */
      for (i=0; (i < nr); i++) {
	eps4    = 4.0 * atype->nb[i].C1; 
	sig6    = pow (atype->nb[i].C0, 6.0);
	c6[i]   = eps4  * sig6;
	c12[i]  = c6[i] * sig6;
      }

      /* Fill the structure */
      for (i=k=0; (i < nr); i++)
	for (j=0; (j < nr); j++,k++) {
	  plist->param[k].c[0] = sqrt (c6[i]  * c6[j]);
	  plist->param[k].c[1] = sqrt (c12[i] * c12[j]);
	}
      
      /* Free temp storage */
      sfree(c6);
      sfree(c12);

      break;
    }
    assert(plist->nr == k);
    break;
    
  case F_BHAM:
    sprintf(errbuf,
	    "Buckingham parameters can not be generated automatically\n");
    warning(errbuf);
    
    break;
  default:
    sprintf(errbuf,"Invalid nonbonded type %s",
	    interaction_function[ftype].longname);
    warning(errbuf);
  }
}

void push_at (t_symtab *symtab, t_atomtype *at, char *line)
{
  typedef struct {
    char *entry;
    int  ptype;
  } t_xlate;
  t_xlate xl[eptNR] = {
    { "A",   eptAtom }, 
    { "N",   eptNucleus }, 
    { "S",   eptShell }, 
    { "B",   eptBond },
    { "D",   eptDummy }, 
  };
  
  int    nr,i,b,j,pt;
  char   type[STRLEN],ptype[STRLEN];
  double m,q;
  double c[MAXFORCEPARAM];
  char   f1[STRLEN],f2[STRLEN];
  
  nr = at->nr;

  if (sscanf (line,"%s%lf%lf%s%lf%lf",type,&m,&q,ptype,&c[0],&c[1]) != 6) {
    too_few();
    return;
  }
  for(j=2; (j<MAXFORCEPARAM); j++)
    c[j]=0.0;
    
  for(j=0; (j<eptNR); j++)
    if (strcasecmp(ptype,xl[j].entry) == 0)
      break;
  if (j == eptNR)
    fatal_error(0,"Invalid particle type %s on line %s",
		ptype,line);
  pt=xl[j].ptype;
  if (debug)
    fprintf(debug,"ptype: %s\n",ptype_str[pt]);
  
  /* Test if this type overwrites another */
  for (i=0; ((i<nr) && (strcasecmp(*(at->atomname[i]),type) != 0)); i++)
    ;

  if ((i==nr) || (nr==0)) {
    /* New atomtype, get new space for arrays */
    srenew(at->atom,nr+1);
    srenew(at->atomname,nr+1);
    srenew(at->nb,nr+1);
    at->nr++;
  }
  else {
    sprintf(errbuf,"Overriding atomtype %s",type);
    warning(errbuf);
    nr = i;
  }
  /* fill the arrays */
  at->atomname[nr] = put_symtab(symtab,type);
  at->atom[nr].ptype = pt;
  at->atom[nr].m = m;
  at->atom[nr].q = q;
  
  for (i=0; (i<2); i++)
    at->nb[nr].c[i] = c[i];
  for( ; (i<MAXFORCEPARAM); i++)
    at->nb[nr].c[i] = 0.0;
}

static void push_bondtype(t_params *bt,t_param *b,int nral,int ftype)
{
  int  i,j;
  bool bTest,bFound;  
  int  nr   = bt->nr;
  int  nrfp = NRFP(ftype);
  
  /* Check if this entry overwrites another */
  bFound=FALSE;
  for (i=0; (i < nr); i++) {
    bTest = TRUE;
    for (j=0; (j < nral); j++)
      bTest=(bTest && (b->a[j] == bt->param[i].a[j]));
    if (!bTest) {
      bTest=TRUE;
      for(j=0; (j<nral); j++)
	bTest=(bTest && (b->a[nral-1-j] == bt->param[i].a[j]));
    }
    if (bTest) {
      /* Overwrite it! */
      for (j=0; (j < nrfp); j++)
	bt->param[i].c[j] = b->c[j];
      bFound=TRUE;
    }
  }
  if (!bFound) {
    /* alloc */
    pr_alloc (2,bt);
    
    /* fill the arrays up and down */
    memcpy((char *) bt->param[nr].c,	 (char *) b->c,sizeof(b->c));
    memcpy((char *) bt->param[nr].a,	 (char *) b->a,sizeof(b->a));
    memcpy((char *) bt->param[nr+1].c,	 (char *) b->c,sizeof(b->c));
    for (j=0; (j < nral); j++) 
      bt->param[nr+1].a[j] = b->a[nral-1-j];
    
    bt->nr += 2;
  }
}

void push_bt(directive d,t_params bt[],int nral,t_atomtype *at,char *line)
{
  static   char *formal[5] = {
    "%s",
    "%s%s",
    "%s%s%s",
    "%s%s%s%s",
    "%s%s%s%s%s"
  };
  static   char *formnl[5] = {
    "%*s",
    "%*s%*s",
    "%*s%*s%*s",
    "%*s%*s%*s%*s",
    "%*s%*s%*s%*s%*s"
  };
  static   char *formlf[MAXFORCEPARAM] = {
    "%lf",
    "%lf%lf",
    "%lf%lf%lf",
    "%lf%lf%lf%lf",
    "%lf%lf%lf%lf%lf",
    "%lf%lf%lf%lf%lf%lf",
  };
  int      i,j,ft,ftype,nn,nrfp;
  char     f1[STRLEN];
  char     alc[5][20];
  double   c[MAXFORCEPARAM];
  t_param  p;
  
  /* Make format string (nral ints+functype) */
  if ((nn=sscanf(line,formal[nral],
		 alc[0],alc[1],alc[2],alc[3],alc[4])) != nral+1) {
    sprintf(errbuf,"Not enough atomtypes (%d instead of %d)",nn-1,nral);
    warning(errbuf);
    return;
  }
  
  ft    = atoi(alc[nral]);
  ftype = ifunc_index(d,ft);
  nrfp  = NRFP(ftype);
  strcpy(f1,formnl[nral]);
  strcat(f1,formlf[nrfp-1]);
  if ((nn=sscanf(line,f1,&c[0],&c[1],&c[2],&c[3],&c[4],&c[5])) 
      != nrfp) {
    for( ; (nn<nrfp); nn++)
      c[nn] = 0.0;
  }
  for(i=0; (i<nral); i++)
    p.a[i]=at2type(alc[i],at);
  for(i=0; (i<nrfp); i++)
    p.c[i]=c[i];
  push_bondtype (&(bt[ftype]),&p,nral,ftype);
}

void push_nbt(directive d,t_params nbt[],t_atomtype *atype,
	      char *pline,int nb_funct)
{
  /* swap the atoms */
  static char *form2="%*s%*s%*s%lf%lf";
  static char *form3="%*s%*s%*s%lf%lf%lf";
  char    a0[80],a1[80];
  int     i,f,k,ftype,atnr,nrfp;
  double  c[3];
  atom_id ai,aj;
  
  if (sscanf (pline,"%s%s%d",a0,a1,&f) != 3) {
    too_few();
    return;
  }
  
  ftype=ifunc_index(d,f);
  
  if (ftype != nb_funct) {
    sprintf(errbuf,"Trying to add %s while the default nonbond type is %s",
	    interaction_function[ftype].longname,
	    interaction_function[nb_funct].longname);
    warning(errbuf);
    return;
  }
  
  /* Get the force parameters */
  if (NRFP(ftype) == 2) {
    if (sscanf(pline,form2,&c[0],&c[1]) != 2) {
      too_few();
      return;
    }
  }
  else {
    if (sscanf(pline,form3,&c[0],&c[1],&c[2]) != 3) {
      too_few();
      return;
    }
  }
  /* Put the parameters in the matrix */
  ai = at2type (a0,atype);
  aj = at2type (a1,atype);

  atnr = atype->nr;
  nrfp = NRFP(ftype);
  for (i=0; (i < nrfp); i++) {
    k = atnr*ai+aj;
    nbt[ftype].param[k].c[i] = c[i];
    k = atnr*aj+ai;
    nbt[ftype].param[k].c[i] = c[i];
  }
}

static void push_atom_now(t_symtab *symtab,t_atoms *at,int atomnr,
			  int type,int ptype,int resnumber,int cgnumber,
			  char *resname,char *name,real m0,real q0,
			  int typeB,real mB,real qB)
{
  int j;
  int nr = at->nr;
  /* New atom instance
   * get new space for arrays 
   */
  srenew(at->atom,nr+1);
  srenew(at->atomname,nr+1);

  if (resnumber > at->nres) {
    at->nres=resnumber;
    srenew(at->resname,resnumber);
    at->resname[resnumber-1] = put_symtab(symtab,resname);
  }
  /* fill the list */
  at->atom[nr].type  = type;
  at->atom[nr].ptype = ptype;
  at->atom[nr].q     = q0;
  at->atom[nr].m     = m0;
  at->atom[nr].typeB = typeB;
  at->atom[nr].qB    = qB;
  at->atom[nr].mB    = mB;
  
  for(j=0; (j<egcNR); j++)
    at->atom[nr].grpnr[j] = -1;
  at->atom[nr].resnr = resnumber-1;
  at->atomname[nr] = put_symtab(symtab,name);
  at->nr++;
}

void push_cg(t_block *block, int *lastindex, int index, int a)
{
  if (debug)
    fprintf (debug,"Index %d, Atom %d\n",index,a);

  if (((block->nr) && (*lastindex != index)) || (!block->nr)) {
    /* add a new block */
    block->nr++;
    srenew(block->index,block->nr+1);
  }
  srenew(block->a,block->nra+1);
  block->a[block->nra++]=a;
  block->index[block->nr]=block->nra;
  *lastindex = index;
}

void push_atom(t_symtab *symtab,t_block *cgs,
	       t_atoms *at,t_atomtype *atype,char *line)
{
  int 		nr,ptype;
  static int    lastcg;
  int 		resnumber,cgnumber,atomnr,type,typeB,typeb,nscan;
  char 		id[STRLEN],ctype[STRLEN],
       		resname[STRLEN],name[STRLEN];
  double        m,q,mb,qb;
  real          m0,q0,mB,qB;
  char 		f1[STRLEN];

  /* Make a shortcut for writing in this molecule  */
  nr = at->nr;

  /* Fixed parameters */
  if (sscanf(line,"%s%s%d%s%s%d",
	     id,ctype,&resnumber,resname,name,&cgnumber) != 6) {
    too_few();
    return;
  }
  sscanf(id,"%d",&atomnr);
  type  = at2type(ctype,atype);
  ptype = atype->atom[type].ptype;

  /* Set default from type */  
  q0    = atype->atom[type].q;
  m0    = atype->atom[type].m;
  typeB = type;
  qB    = q0;
  mB    = m0;
  
  /* Optional parameters */
  nscan=sscanf(line,"%*s%*s%*s%*s%*s%*s%lf%lf%s%lf%lf",
	       &q,&m,ctype,&qb,&mb);
  
  /* Nasty switch that falls thru all the way down! */
  if (nscan > 0) {
    q0 = qB = q;
    if (nscan > 1) {
      m0 = mB = m;
      if (nscan > 2) {
	typeB=at2type(ctype,atype);
	qB = atype->atom[typeB].q;
	mB = atype->atom[typeB].m;
	if (nscan > 3) {
	  qB = qb;
	  if (nscan > 4)
	    mB = mb;
	}
      }
    }
  }
  if (debug) 
    fprintf(debug,"mB=%g, qB=%g, typeB=%d\n",mB,qB,typeB);
  
  push_cg(cgs,&lastcg,cgnumber,nr);

  push_atom_now(symtab,at,atomnr,type,ptype,resnumber,cgnumber,
		resname,name,m0,q0,typeB,mB,qB);
}

void push_molt(t_symtab *symtab,t_molinfo *newmol,char *line)
{
  char type[STRLEN];
  int nrexcl;

  if ((sscanf(line,"%s%d",type,&nrexcl)) != 2) {
    too_few();
    return;
  }
  
  /* Fill in the values */
  newmol->name     = put_symtab(symtab,type);
  newmol->nrexcl   = nrexcl;
  newmol->excl_set = FALSE;
}

static bool default_params(int ftype,t_params bt[],t_atoms *at,t_param *p,
			   bool bB)
{
  int      i,j;
  bool     bFound;
  t_param  *pi=NULL;
  int      nr   = bt[ftype].nr;
  int      nral = NRAL(ftype);
  int      nrfp = interaction_function[ftype].nrfpA;
  int      nrfpB= interaction_function[ftype].nrfpB;

  bFound=FALSE;
  for (i=0; ((i < nr) && !bFound); i++) {
    pi=&(bt[ftype].param[i]);
    if ((ftype == F_PDIHS)  || (ftype == F_RBDIHS)) {
      /* The j and k atoms are decisive about which dihedral type
       * we should take.
       */
      if (bB)
	bFound=((at->atom[p->AJ].typeB==pi->AI) &&
		(at->atom[p->AK].typeB==pi->AJ));
      else
	bFound=((at->atom[p->AJ].type==pi->AI) &&
		(at->atom[p->AK].type==pi->AJ));
    }
    else if (ftype == F_IDIHS) {
      /* The i and l atoms are decisive about which dihedral type
       * we should take.
       */
      if (bB)
	bFound=((at->atom[p->AI].typeB==pi->AI) &&
		(at->atom[p->AL].typeB==pi->AJ));
      else
	bFound=((at->atom[p->AI].type==pi->AI) &&
		(at->atom[p->AL].type==pi->AJ));
    }
    else {
      if (bB)
	for (j=0; ((j < nral) && 
		   (at->atom[p->a[j]].typeB == pi->a[j])); j++);
      else
	for (j=0; ((j < nral) && 
		   (at->atom[p->a[j]].type == pi->a[j])); j++);
      bFound=(j==nral);
    }
  }
  if (bFound) {
    if (bB) {
      assert(nrfp+nrfpB <= MAXFORCEPARAM);
      for (j=0; (j < nrfpB); j++)
	p->c[nrfp+j] = pi->c[j];
    }
    else
      for (j=0; (j < nrfp); j++)
	p->c[j] = pi->c[j];
  }
  else {
    for (j=0; (j < nrfp); j++)
      p->c[j] = 0.0;
  }
  return bFound;
}

void push_bondnow(t_params *bond, t_param *b)
{
  int nr   = bond->nr;

  if (debug) {
    fprintf(debug,"push_bondnow: nr = %d\n",nr);
    fflush(debug);
  }
  /* allocate one position extra */
  pr_alloc (1,bond);

  /* fill the arrays */
  memcpy ((char *)&(bond->param[nr]),(char *)b,sizeof(t_param));

  bond->nr++;
}

void push_bond(directive d,t_params bondtype[],t_params bond[],
	       t_atoms *at,char *line)
{
  static char *aaformat[4]= {
    "%d%d",
    "%d%d%d",
    "%d%d%d%d",
    "%d%d%d%d%d"
  };
  static char *asformat[4]= {
    "%*s%*s",
    "%*s%*s%*s",
    "%*s%*s%*s%*s",
    "%*s%*s%*s%*s%*s"
  };
  static char *ccformat[6]= {
    "%lf",
    "%lf%lf",
    "%lf%lf%lf",
    "%lf%lf%lf%lf",
    "%lf%lf%lf%lf%lf",
    "%lf%lf%lf%lf%lf%lf"
  };
  int      nr,j,nrfp,nrfpA,nral,nread,ftype;
  char     format[STRLEN];
  double   cc[MAXFORCEPARAM];
  int      aa[MAXATOMLIST+1];
  t_param  param,paramB;
  bool     bFoundA,bFoundB,bDef,bPert;
  
  ftype = ifunc_index(d,1);
  nral  = NRAL(ftype);
  nread = sscanf(line,aaformat[nral-1],&aa[0],&aa[1],&aa[2],&aa[3],&aa[4]);
  if (nread < nral) {
    too_few();
    return;
  }
  else if (nread == nral) 
    ftype = ifunc_index(d,1);
  else
    ftype = ifunc_index(d,aa[nral]);

  /* default force parameters  */
  for(j=0; (j<MAXATOMLIST); j++)
    param.a[j] = aa[j]-1;
  for(j=0; (j<MAXFORCEPARAM); j++)
    param.c[j] = 0.0;
  
  /* Get force params for normal and free energy perturbation
   * studies, as determined by types!
   */
  bFoundA = default_params(ftype,bondtype,at,&param,FALSE);
  bFoundB = default_params(ftype,bondtype,at,&param,TRUE);
  bDef    = TRUE;
  
  nrfp  = NRFP(ftype);
  if (nread > nral) {  
    strcpy(format,asformat[nral-1]);
    strcat(format,ccformat[nrfp-1]);
    
    nread = sscanf(line,format,&cc[0],&cc[1],&cc[2],&cc[3],&cc[4],&cc[5]);
    if (nread > nrfp) {
      warning("Too many parameters");
      nread = nrfp;
    }
    
    for(j=0; (j<nread); j++)
      param.c[j]=cc[j];
      
    /* Check whether we have to use the defaults */
    if (nread == nrfp)
      bDef = FALSE;
  }
  else {
    nread = 0;
  }
  /* nread now holds the number of force parameters read! */
  
  if (bDef) {
    /* Use defaults */
    nrfpA=interaction_function[ftype].nrfpA;
  
    if (nread < nrfpA) {
      if (!bFoundA) {
	sprintf(errbuf,"No default %s types, using zeroes",
		interaction_function[ftype].longname);
	warning(errbuf);
      }
    }
    else {
      /* This implies nread < nrfp, and so perturbed values have not
       * been read 
       */
      if (!bFoundB) {
	/* implies !bFoundB */
	/* We only have to issue a warning if these atoms are perturbed! */
	bPert = FALSE;
	for(j=0; (j<nral); j++)
	  bPert = bPert || PERTURBED(at->atom[param.a[j]]);
	
	if (bPert) {
	  sprintf(errbuf,"No default %s types for perturbed atoms, using normal values",interaction_function[ftype].longname);
	  warning(errbuf);
	}
	for(j=nrfpA; (j<nrfp); j++)
	  param.c[j]=param.c[j-nrfpA];
      }
    }
  }
  /* Put the values in the appropriate arrays */
  push_bondnow (&bond[ftype],&param);
}

void push_mol(int nrmols,t_molinfo mols[],char *pline,int *whichmol,
		  int *nrcopies)
{
  int  i,copies;
  char type[STRLEN];

  *nrcopies=0;
  if (sscanf(pline,"%s%d",type,&copies) != 2) {
    too_few();
    return;
  }
  
  /* search moleculename */
  for (i=0; ((i<nrmols) && strcasecmp(type,*(mols[i].name))); i++)
    ;

  if (i<nrmols) {
    *nrcopies        = copies;
    *whichmol        = i;
  }
  else {
    sprintf(errbuf,"No such moleculetype %s",type);
    warning(errbuf);
  }
}

void init_block2(t_block2 *b2, int natom)
{
  int i;

  b2->nr=natom;
  snew(b2->nra,b2->nr);
  snew(b2->a,b2->nr);
  for(i=0; (i<b2->nr); i++)
    b2->a[i]=NULL;
}

void done_block2(t_block2 *b2)
{
  int i;
  
  if (b2->nr) {
    for(i=0; (i<b2->nr); i++)
      sfree(b2->a[i]);
    sfree(b2->a);
    sfree(b2->nra);
    b2->nr=0;
  }
}

void push_excl(char *line, t_block2 *b2)
{
  int  i,j;
  int  n;
  char base[STRLEN],format[STRLEN];

  if (sscanf(line,"%d",&i)==0)
    return;
    
  if ((1 <= i) && (i <= b2->nr))
    i--;
  else {
    if (debug)
      fprintf(debug,"Unbound atom %d\n",i-1);
    return;
  }
  strcpy(base,"%*d");
  do {
    strcpy(format,base);
    strcat(format,"%d");
    n=sscanf(line,format,&j);
    if (n == 1) {
      if ((1 <= j) && (j <= b2->nr)) {
	srenew(b2->a[i],++(b2->nra[i]));
	b2->a[i][b2->nra[i]-1]=j-1;
	strcat(base,"%*d");
      }
      else 
	fatal_error(0,"Invalid Atomnr j: %d, b2->nr: %d\n",j,b2->nr);
    }
  } while (n == 1);
}

void b_to_b2(t_block *b, t_block2 *b2)
{
  int     i;
  atom_id j,a;

  for(i=0; (i<b->nr); i++)
    for(j=b->index[i]; (j<b->index[i+1]); j++) {
      a=b->a[j];
      srenew(b2->a[i],++b2->nra[i]);
      b2->a[i][b2->nra[i]-1]=a;
    }
}

void b2_to_b(t_block2 *b2, t_block *b)
{
  int     i,nra;
  atom_id j;

  nra=0;
  for(i=0; (i<b2->nr); i++) {
    b->index[i]=nra;
    for(j=0; (j<b2->nra[i]); j++)
      b->a[nra+j]=b2->a[i][j];
    nra+=b2->nra[i];
  }
  /* terminate list */
  b->index[i]=nra;
}

static int icomp(const void *v1, const void *v2)
{
  return (*((atom_id *) v1))-(*((atom_id *) v2));
}

void merge_excl(t_block *excl, t_block2 *b2)
{
  int     i,k;
  atom_id j;
  int     nra;

  if (!b2->nr)
    return;
  else if (b2->nr != excl->nr) {
    fatal_error(0,"DEATH HORROR: b2->nr = %d, while excl->nr = %d",
		b2->nr,excl->nr);
  }
  else
    printf("Entering merge_excl\n");

  /* First copy all entries from excl to b2 */
  b_to_b2(excl,b2);

  /* Count and sort the exclusions */
  nra=0;
  for(i=0; (i<b2->nr); i++) {
    if (b2->nra[i] > 0) {
      /* remove double entries */
      qsort(b2->a[i],b2->nra[i],sizeof(b2->a[i][0]),icomp);
      k=1;
      for(j=1; (j<b2->nra[i]); j++)
	if (b2->a[i][j]!=b2->a[i][k-1]) {
	  b2->a[i][k]=b2->a[i][j];
	  k++;
	}
      b2->nra[i]=k;
      nra+=b2->nra[i];
    }
  }
  excl->nra=nra;
  srenew(excl->a,excl->nra);

  b2_to_b(b2,excl);
}

