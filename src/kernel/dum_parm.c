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
static char *SRCID_dum_parm_c = "$Id$";

#include <stdio.h>
#include <math.h>
#include <string.h>
#include "assert.h"
#include "dum_parm.h"
#include "smalloc.h"
#include "resall.h"
#include "add_par.h"
#include "vec.h"
#include "toputil.h"
#include "physics.h"
#include "index.h"
#include "names.h"
#include "fatal.h"

static void default_params(int ftype,t_params bt[],t_atoms *at,t_param *p)
{
  int      i,j;
  bool     bFound;
  t_param  *pi=NULL;

  bFound=FALSE;
  for (i=0; ((i < bt[ftype].nr) && !bFound); i++) {
    pi=&(bt[ftype].param[i]);
    switch(ftype) {
    case F_PDIHS:
    case F_RBDIHS:
      bFound=((at->atom[p->AJ].type==pi->AI) &&
	      (at->atom[p->AK].type==pi->AJ));
      break;
    case F_IDIHS:
      bFound=((at->atom[p->AI].type==pi->AI) &&
	      (at->atom[p->AL].type==pi->AJ));
      break;
    default:
      bFound=TRUE;
      for (j=0; j<NRAL(ftype) && bFound; j++)
	bFound = (at->atom[p->a[j]].type == pi->a[j]);
    } /* switch */
  }
  if (bFound)
    for (j=0; (j < interaction_function[ftype].nrfpA); j++)
      p->c[j] = pi->c[j];
  else
    /* if not found, return all NOTSET */
    for (j=0; (j < interaction_function[ftype].nrfpA); j++)
      p->c[j] = NOTSET;
}

typedef struct {
  t_iatom a[2];
  real    c;
} t_mybond;

typedef struct {
  t_iatom a[3];
  real    c;
} t_myang;

typedef struct {
  t_iatom a[4];
  real    c;
} t_myidih;

static void enter_bond(int *nrbond, t_mybond **bonds, t_param param, 
		       int ftype, t_params ptype[], t_atoms *atoms)
{
  int j;
  
  (*nrbond)++;
  srenew(*bonds, *nrbond);
  
  /* copy atom numbers */
  for(j=0; j<2; j++)
    (*bonds)[*nrbond-1].a[j] = param.a[j];
  
  /* but if ftype is dummy we don't have the parameters yet */
  if (interaction_function[ftype].flags & IF_DUMMY)
    default_params(F_BONDS,ptype,atoms,&param);
  
  /* copy bondlength */
  (*bonds)[*nrbond-1].c = param.C0;
}

static void enter_angle(int *nrang, t_myang **angles, t_param param)
{
  int j;
  
  (*nrang)++;
  srenew(*angles, *nrang);
  
  for(j=0; j<3; j++)
    (*angles)[*nrang-1].a[j] = param.a[j];
  (*angles)[*nrang-1].c = param.C0; /* angle */
}

static void enter_idih(int *nridih, t_myidih **idihs, t_param param)
{
  int j;
  
  (*nridih)++;
  srenew(*idihs, *nridih);
  for(j=0; j<4; j++)
    (*idihs)[*nridih-1].a[j] = param.a[j];
  (*idihs)[*nridih-1].c = param.C0; /* angle */
}

static void get_bondeds(int nrat, t_iatom atoms[], 
			t_params ptype[], t_params plist[], t_atoms *at,
			int *nrbond, t_mybond **bonds,
			int *nrang,  t_myang  **angles,
			int *nridih, t_myidih **idihs )
{
  int     i,j,k,ftype;
  int     nra,nrd,tp,nrcheck;
  t_iatom *ia,adum;
  bool    bCheck;
  t_param param;
  
  if (debug) {
    fprintf(debug,"getting bondeds for %d (nr=%d):",atoms[0]+1,nrat);
    for(k=1; k<nrat; k++)
      fprintf(debug," %d",atoms[k]+1);
    fprintf(debug,"\n");
  }
  for(ftype=0; (ftype<F_NRE); ftype++) {
    if (interaction_function[ftype].flags & IF_DUMMY)
      /* this is a dummy, we treat it like a bond */
      nrcheck = 2;
    else
      switch(ftype) {
      case F_BONDS:
      case F_SHAKE:
	nrcheck = 2;
	break;
      case F_ANGLES: 
	nrcheck = 3;
	break;
      case F_IDIHS: 
	nrcheck = 4;
	break;
      default:
	nrcheck = 0;
      } /* case */
    
    if (nrcheck)
      for(i=0; i<plist[ftype].nr; i++) {
	/* now we have something, check if it is between atoms[] */
	bCheck=TRUE;
	for(j=0; j<nrcheck && bCheck; j++) {
	  bCheck=FALSE;
	  for(k=0; k<nrat; k++)
	    bCheck = bCheck || (plist[ftype].param[i].a[j]==atoms[k]);
	}
	
	if (bCheck)
	  /* use nrcheck to see if we're adding bond, angle or idih */
	  switch (nrcheck) {
	  case 2: enter_bond (nrbond,bonds, plist[ftype].param[i],
			      ftype,ptype,at);
	  break;
	  case 3: enter_angle(nrang, angles,plist[ftype].param[i]);
	    break;
	  case 4: enter_idih (nridih,idihs, plist[ftype].param[i]);
	    break;
	  } /* case */
      } /* for i */
  } /* for ftype */
}

/* for debug */
static void print_bad(FILE *fp, 
		      int nrbond, t_mybond *bonds,
		      int nrang,  t_myang  *angles,
		      int nridih, t_myidih *idihs )
{
  int i;
  
  if (nrbond) {
    fprintf(fp,"bonds:");
    for(i=0; i<nrbond; i++)
      fprintf(fp," %u-%u (%g)", 
	      bonds[i].AI+1, bonds[i].AJ+1, bonds[i].c);
    fprintf(fp,"\n");
  }
  if (nrang) {
    fprintf(fp,"angles:");
    for(i=0; i<nrang; i++)
      fprintf(fp," %u-%u-%u (%g)", 
	      angles[i].AI+1, angles[i].AJ+1, 
	      angles[i].AK+1, angles[i].c);
    fprintf(fp,"\n");
  }
  if (nridih) {
    fprintf(fp,"idihs:");
    for(i=0; i<nridih; i++)
      fprintf(fp," %u-%u-%u-%u (%g)", 
	      idihs[i].AI+1, idihs[i].AJ+1, 
	      idihs[i].AK+1, idihs[i].AL+1, idihs[i].c);
    fprintf(fp,"\n");
  }
}

static void print_param(FILE *fp, int ftype, int i, t_param *param)
{
  static int pass = 0;
  static int prev_ftype= NOTSET;
  static int prev_i    = NOTSET;
  int j;
  
  if ( (ftype!=prev_ftype) || (i!=prev_i) ) {
    pass = 0;
    prev_ftype= ftype;
    prev_i    = i;
  }
  fprintf(fp,"(%d) plist[%s].param[%d]",
	  pass,interaction_function[ftype].name,i);
  for(j=0; j<NRFP(ftype); j++)
    fprintf(fp,".c[%d]=%g ",j,param->c[j]);
  fprintf(fp,"\n");
  pass++;
}

static real get_def_bond_length(t_params *bt, t_atoms *at,
				t_iatom ai, t_iatom aj)
{
  int  i;
  real bondlen;

  bondlen=NOTSET;
  for (i=0; i < bt->nr && (bondlen==NOTSET); i++) {
    /* check one way only, bt is filled with both */
    if ( (at->atom[ai].type == bt->param[i].AI) &&
	 (at->atom[aj].type == bt->param[i].AJ) )
      bondlen = bt->param[i].C0;
  }
  if (bondlen==NOTSET)
    fatal_error(0,"No default bondlength for atoms %d-%d while "
		"calculating dummy parameters",ai+1,aj+1);
  return bondlen;
}

static real get_def_angle(t_params *bt, t_atoms *at,
			  t_iatom ai, t_iatom aj, t_iatom ak)
{
  int  i;
  real angle;
  
  angle=NOTSET;
  for (i=0; i < bt->nr && (angle==NOTSET); i++) {
    /* check one way only, bt is filled with both */
    if ( (at->atom[ai].type == bt->param[i].AI) &&
	 (at->atom[aj].type == bt->param[i].AJ) &&
	 (at->atom[ak].type == bt->param[i].AK) )
      angle = DEG2RAD*bt->param[i].C0;
  }
  if (angle==NOTSET)
    fatal_error(0,"No default angle for atoms %d-%d-%d while "
		"calculating dummy parameters",ai+1,aj+1,ak+1);
  return angle;
}

static real get_bond_length(int nrbond, t_mybond bonds[], 
			    t_iatom ai, t_iatom aj)
{
  int  i;
  real bondlen;
  
  bondlen=NOTSET;
  for (i=0; i < nrbond && (bondlen==NOTSET); i++) {
    /* check both ways */
    if ( ( (ai == bonds[i].AI) && (aj == bonds[i].AJ) ) || 
	 ( (ai == bonds[i].AJ) && (aj == bonds[i].AI) ) )
      bondlen = bonds[i].c; /* note: bonds[i].c might be NOTSET */
  }
  /* no fatal_error here: just return NOTSET if not found */
  return bondlen;
}

static real get_angle(int nrang, t_myang angles[], 
		      t_iatom ai, t_iatom aj, t_iatom ak)
{
  int  i;
  real angle;
  
  angle=NOTSET;
  for (i=0; i < nrang && (angle==NOTSET); i++) {
    /* check both ways */
    if ( ( (ai==angles[i].AI) && (aj==angles[i].AJ) && (ak==angles[i].AK) ) || 
	 ( (ai==angles[i].AK) && (aj==angles[i].AJ) && (ak==angles[i].AI) ) )
      angle = DEG2RAD*angles[i].c;
  }
  /* no fatal_error here: just return NOTSET if not found */
  return angle;
}

static real get_an_angle(int nrang, t_myang angles[], 
			 t_params *bt, t_atoms *at,
			 t_iatom ai, t_iatom aj, t_iatom ak)
{
  real angle;
  
  angle = get_angle(nrang, angles, ai, aj, ak);
  if (angle==NOTSET)
    angle = get_def_angle(bt, at, ai, aj, ak);
  
  return angle;
}

static bool calc_dum3_param(t_params ptype[], t_atomtype *atype,
			       t_param *param, t_atoms *at,
			       int nrbond, t_mybond *bonds,
			       int nrang,  t_myang  *angles,
			       int nridih, t_myidih *idihs )
{
  /* i = dummy atom            |    ,k
   * j = 1st bonded heavy atom | i-j
   * k,l = 2nd bonded atoms    |    `l
   */
  
  bool bXH3,bError;
  real bij,bjk,bjl,aijk,aijl,akjl,pijk,pijl,a=-1,b=-1,c;
  /* check if this is part of a NH3 or CH3 group,
   * i.e. if atom k and l are dummy masses (MNH3 or MCH3) */
  if (debug) {
    int i;
    for (i=0; i<4; i++)
      fprintf(debug,"atom %u type %s ",
	      param->a[i]+1,type2nm(at->atom[param->a[i]].type,atype));
    fprintf(debug,"\n");
  }
  bXH3 = 
    ( (strcasecmp(type2nm(at->atom[param->AK].type,atype),"MNH3")==0) &&
      (strcasecmp(type2nm(at->atom[param->AL].type,atype),"MNH3")==0) ) ||
    ( (strcasecmp(type2nm(at->atom[param->AK].type,atype),"MCH3")==0) &&
      (strcasecmp(type2nm(at->atom[param->AL].type,atype),"MCH3")==0) );
  
  bjk = get_bond_length(nrbond, bonds, param->AJ, param->AK);
  bjl = get_bond_length(nrbond, bonds, param->AJ, param->AL);
  bError = (bjk==NOTSET) || (bjl==NOTSET);
  if (bXH3) {
    /* now we get some XH3 group specific construction */
    /* note: we call the heavy atom 'C' and the X atom 'N' */
    real bMM,bCM,bCN,bNH,aCNH,dH,rH,dM,rM;
    t_param xparam;
    int aN, xatoms[3];
    
    /* check if bonds from heavy atom (j) to dummy masses (k,l) are equal: */
    bError = bError || (bjk!=bjl);
    
    /* the X atom (C or N) in the XH3 group is the first after the masses: */
    aN = max(param->AK,param->AL)+1;
    
    /* get common bonds */
    bMM = get_bond_length(nrbond, bonds, param->AK, param->AL);
    bCM = bjk;
    bCN = get_def_bond_length(&ptype[F_BONDS], at, param->AJ, aN);
    
    /* calculate common things */
    rM  = 0.5*bMM;
    dM  = sqrt( sqr(bCM) - sqr(rM) );
    
    /* are we dealing with the X atom? */
    if ( param->AI == aN ) {
      /* this is trivial */
      a = b = 0.5 * bCN/dM;
      
    } else {
      /* get other bondlengths and angles: */
      bNH = get_def_bond_length(&ptype[F_BONDS], at, aN, param->AI);
      aCNH= get_def_angle(&ptype[F_ANGLES], at, param->AJ, aN, param->AI);
      
      /* calculate */
      dH  = bCN - bNH * cos(aCNH);
      rH  = bNH * sin(aCNH);
      
      a = 0.5 * ( dH/dM + rH/rM );
      b = 0.5 * ( dH/dM - rH/rM );
    }
  } else
    fatal_error(0,"calc_dum3_param not implemented for the general case "
		"(atom %d)",param->AI+1);
  
  param->C0 = a;
  param->C1 = b;
  
  if (debug)
    fprintf(debug,"params for dummy3 %u: %g %g\n",
	    param->AI+1,param->C0,param->C1);
  
  return bError;
}

static bool calc_dum3fd_param(t_params ptype[], t_param *param, t_atoms *at,
			      int nrbond, t_mybond *bonds,
			      int nrang,  t_myang  *angles,
			      int nridih, t_myidih *idihs )
{
  /* i = dummy atom            |    ,k
   * j = 1st bonded heavy atom | i-j
   * k,l = 2nd bonded atoms    |    `l
   */

  bool bError;
  real bij,bjk,bjl,aijk,aijl,rk,rl;
  
  bij = get_def_bond_length(&ptype[F_BONDS], at, param->AI, param->AJ);
  bjk = get_bond_length(nrbond, bonds, param->AJ, param->AK);
  bjl = get_bond_length(nrbond, bonds, param->AJ, param->AL);
  aijk= get_an_angle(nrang, angles, &ptype[F_ANGLES], at, 
		     param->AI, param->AJ, param->AK);
  aijl= get_an_angle(nrang, angles, &ptype[F_ANGLES], at,
		     param->AI, param->AJ, param->AL);
  bError = (bjk==NOTSET) || (bjl==NOTSET);
  
  rk = bjk * sin(aijk);
  rl = bjl * sin(aijl);
  param->C0 = rk / (rk + rl);
  param->C1 = -bij; /* 'bond'-length for fixed distance dummy */
  
  if (debug)
    fprintf(debug,"params for dummy3fd %u: %g %g\n",
	    param->AI+1,param->C0,param->C1);
  return bError;
}

static bool calc_dum3fad_param(t_params ptype[], t_param *param, t_atoms *at,
			       int nrbond, t_mybond *bonds,
			       int nrang,  t_myang  *angles,
			       int nridih, t_myidih *idihs )
{
  /* i = dummy atom            |
   * j = 1st bonded heavy atom | i-j
   * k = 2nd bonded heavy atom |    `k-l
   * l = 3d bonded heavy atom  |
   */

  bool bSwapParity;
  real bij,aijk;
  
  bSwapParity = ( param->C1 == -1 );
  
  bij = get_def_bond_length(&ptype[F_BONDS], at, param->AI, param->AJ);
  aijk = get_def_angle(&ptype[F_ANGLES], at, param->AI, param->AJ, param->AK);
  
  param->C1 = bij;          /* 'bond'-length for fixed distance dummy */
  param->C0 = RAD2DEG*aijk; /* 'bond'-angle for fixed angle dummy */
  
  if (bSwapParity)
    param->C0 = 360 - param->C0;
  
  if (debug)
    fprintf(debug,"params for dummy3fad %u: %g %g\n",
	    param->AI+1,param->C0,param->C1);
  /* if something goes wrong, get_def_angle gives fatal_error */
  return FALSE;
}

static bool calc_dum3out_param(t_params ptype[], t_atomtype *atype,
			       t_param *param, t_atoms *at,
			       int nrbond, t_mybond *bonds,
			       int nrang,  t_myang  *angles,
			       int nridih, t_myidih *idihs )
{
  /* i = dummy atom            |    ,k
   * j = 1st bonded heavy atom | i-j
   * k,l = 2nd bonded atoms    |    `l
   * NOTE: i is out of the j-k-l plane!
   */
  
  bool bXH3,bError,bSwapParity;
  real bij,bjk,bjl,aijk,aijl,akjl,pijk,pijl,a,b,c;
  
  /* check if this is part of a NH3 or CH3 group,
   * i.e. if atom k and l are dummy masses (MNH3 or MCH3) */
  if (debug) {
    int i;
    for (i=0; i<4; i++)
      fprintf(debug,"atom %u type %s ",
	      param->a[i]+1,type2nm(at->atom[param->a[i]].type,atype));
    fprintf(debug,"\n");
  }
  bXH3 = 
    ( (strcasecmp(type2nm(at->atom[param->AK].type,atype),"MNH3")==0) &&
      (strcasecmp(type2nm(at->atom[param->AL].type,atype),"MNH3")==0) ) ||
    ( (strcasecmp(type2nm(at->atom[param->AK].type,atype),"MCH3")==0) &&
      (strcasecmp(type2nm(at->atom[param->AL].type,atype),"MCH3")==0) );
  
  /* check if construction parity must be swapped */  
  bSwapParity = ( param->C1 == -1 );
  
  bjk = get_bond_length(nrbond, bonds, param->AJ, param->AK);
  bjl = get_bond_length(nrbond, bonds, param->AJ, param->AL);
  bError = (bjk==NOTSET) || (bjl==NOTSET);
  if (bXH3) {
    /* now we get some XH3 group specific construction */
    /* note: we call the heavy atom 'C' and the X atom 'N' */
    real bMM,bCM,bCN,bNH,aCNH,dH,rH,rHx,rHy,dM,rM;
    t_param xparam;
    int aN, xatoms[3];
    
    /* check if bonds from heavy atom (j) to dummy masses (k,l) are equal: */
    bError = bError || (bjk!=bjl);
    
    /* the X atom (C or N) in the XH3 group is the first after the masses: */
    aN = max(param->AK,param->AL)+1;
    
    /* get all bondlengths and angles: */
    bMM = get_bond_length(nrbond, bonds, param->AK, param->AL);
    bCM = bjk;
    bCN = get_def_bond_length(&ptype[F_BONDS], at, param->AJ, aN);
    bNH = get_def_bond_length(&ptype[F_BONDS], at, aN, param->AI);
    aCNH= get_def_angle(&ptype[F_ANGLES], at, param->AJ, aN, param->AI);
    
    /* calculate */
    dH  = bCN - bNH * cos(aCNH);
    rH  = bNH * sin(aCNH);
    /* we assume the H's are symmetrically distributed */
    rHx = rH*cos(DEG2RAD*30);
    rHy = rH*sin(DEG2RAD*30);
    rM  = 0.5*bMM;
    dM  = sqrt( sqr(bCM) - sqr(rM) );
    a   = 0.5*( (dH/dM) - (rHy/rM) );
    b   = 0.5*( (dH/dM) + (rHy/rM) );
    c   = rHx / (2*dM*rM);
    
  } else {
    /* this is the general construction */
    
    bij = get_def_bond_length(&ptype[F_BONDS], at, param->AI, param->AJ);
    aijk= get_an_angle(nrang, angles, &ptype[F_ANGLES], at, 
		       param->AI, param->AJ, param->AK);
    aijl= get_an_angle(nrang, angles, &ptype[F_ANGLES], at,
		       param->AI, param->AJ, param->AL);
    akjl= get_angle(nrang, angles, param->AK, param->AJ, param->AL);
  
    pijk = cos(aijk)*bij;
    pijl = cos(aijl)*bij;
    a = ( pijk + (pijk*cos(akjl)-pijl) * cos(akjl) / sqr(sin(akjl)) ) / bjk;
    b = ( pijl + (pijl*cos(akjl)-pijk) * cos(akjl) / sqr(sin(akjl)) ) / bjl;
    c = - sqrt( sqr(bij) - 
		( sqr(pijk) - 2*pijk*pijl*cos(akjl) + sqr(pijl) ) 
		/ sqr(sin(akjl)) )
      / ( bjk*bjl*sin(akjl) );
  }
  
  param->C0 = a;
  param->C1 = b;
  if (bSwapParity)
    param->C2 = -c;
  else
    param->C2 =  c;
  if (debug)
    fprintf(debug,"params for dummy3out %u: %g %g %g\n",
	    param->AI+1,param->C0,param->C1,param->C2);
  return bError;
}

static bool calc_dum4fd_param(t_params ptype[], t_param *param, t_atoms *at,
			      int nrbond, t_mybond *bonds,
			      int nrang,  t_myang  *angles,
			      int nridih, t_myidih *idihs )
{
  /* i = dummy atom            |    ,k
   * j = 1st bonded heavy atom | i-j-m
   * k,l,m = 2nd bonded atoms  |    `l
   */
  
  bool bError;
  real bij,bjk,bjl,bjm,aijk,aijl,aijm,akjm,akjl;
  real pk,pl,pm,cosakl,cosakm,sinakl,sinakm,cl,cm;
  
  bij = get_def_bond_length(&ptype[F_BONDS], at, param->AI, param->AJ);
  bjk = get_bond_length(nrbond, bonds, param->AJ, param->AK);
  bjl = get_bond_length(nrbond, bonds, param->AJ, param->AL);
  bjm = get_bond_length(nrbond, bonds, param->AJ, param->AM);
  aijk= get_an_angle(nrang, angles, &ptype[F_ANGLES], at, 
		     param->AI, param->AJ, param->AK);
  aijl= get_an_angle(nrang, angles, &ptype[F_ANGLES], at,
		     param->AI, param->AJ, param->AL);
  aijm= get_an_angle(nrang, angles, &ptype[F_ANGLES], at,
		     param->AI, param->AJ, param->AM);
  akjm= get_angle(nrang, angles, param->AK, param->AJ, param->AM);
  akjl= get_angle(nrang, angles, param->AK, param->AJ, param->AL);
  bError = (bjk==NOTSET) || (bjl==NOTSET) || (bjm==NOTSET);
  
  pk = bjk*sin(aijk);
  pl = bjl*sin(aijl);
  pm = bjm*sin(aijm);
  cosakl = (cos(akjl) - cos(aijk)*cos(aijl)) / (sin(aijk)*sin(aijl));
  cosakm = (cos(akjm) - cos(aijk)*cos(aijm)) / (sin(aijk)*sin(aijm));
  if ( cosakl < -1 || cosakl > 1 || cosakm < -1 || cosakm > 1 )
    fatal_error(0,"invalid construction in calc_dum4fd for atom %u: "
		"cosakl=%g cosakm=%g\n",param->AI+1,cosakl,cosakm);
  sinakl = sqrt(1-sqr(cosakl));
  sinakm = sqrt(1-sqr(cosakm));
  
  /* note: there is a '+' because of the way the sines are calculated */
  cl = -pk / ( pl*cosakl - pk + pl*sinakl*(pm*cosakm-pk)/(pm*sinakm) );
  cm = -pk / ( pm*cosakm - pk + pm*sinakm*(pl*cosakl-pk)/(pl*sinakl) );
  
  param->C0 = cl;
  param->C1 = cm;
  param->C2 = -bij;
  if (debug)
    fprintf(debug,"params for dummy4fd %u: %g %g %g\n",
	    param->AI+1,param->C0,param->C1,param->C2);
  
  return bError;
}

void set_dummies(bool bVerbose, t_atoms *atoms, t_atomtype atype,
		 t_params ptype[], t_params plist[])
{
  int i,j,ftype;
  int nrbond,nrang,nridih,nrset;
  bool bSet,bERROR=TRUE;
  t_mybond *bonds;
  t_myang  *angles;
  t_myidih *idihs;
  
  if (bVerbose)
    fprintf(stderr,"Calculating parameters for dummy atoms\n");
  if (debug)
    fprintf(debug, "\nCalculating parameters for dummy atoms\n");  
  for(ftype=0; (ftype<F_NRE); ftype++)
    if (interaction_function[ftype].flags & IF_DUMMY) {
      nrset=0;
      for(i=0; (i<plist[ftype].nr); i++) {
	/* check if all parameters are set */
	bSet=TRUE;
	for(j=0; j<NRFP(ftype) && bSet; j++)
	  bSet = plist[ftype].param[i].c[j]!=NOTSET;

	
	if (debug) {
	  fprintf(debug,"bSet=%s ",bool_names[bSet]);
	  print_param(debug,ftype,i,&plist[ftype].param[i]);
	}
	if (!bSet) {
	  nrbond=nrang=nridih=0;
	  bonds = NULL;
	  angles= NULL;
	  idihs = NULL;
	  nrset++;
	  /* now set the dummy parameters: */
	  get_bondeds(NRAL(ftype),plist[ftype].param[i].a, 
		      ptype, plist, atoms, 
		      &nrbond, &bonds, 
		      &nrang,  &angles, 
		      &nridih, &idihs);
	  if (debug) {
	    fprintf(debug, "Found %d bonds, %d angles and %d idihs "
		    "in dummy atom %u (%s)\n",nrbond,nrang,nridih,
		    plist[ftype].param[i].AI+1,
		    interaction_function[ftype].longname);
	    print_bad(debug, nrbond, bonds, nrang, angles, nridih, idihs);
	  } /* debug */
	  switch(ftype) {
	  case F_DUMMY3: 
	    bERROR = 
	      calc_dum3_param(ptype, &atype, 
			      &(plist[ftype].param[i]), atoms,
			      nrbond, bonds, nrang, angles, nridih, idihs);
	    break;
	  case F_DUMMY3FD:
	    bERROR = 
	      calc_dum3fd_param(ptype, &(plist[ftype].param[i]), atoms,
				nrbond, bonds, nrang, angles, nridih, idihs);
	    break;
	  case F_DUMMY3FAD:
	    bERROR = 
	      calc_dum3fad_param(ptype, &(plist[ftype].param[i]), atoms,
				 nrbond, bonds, nrang, angles, nridih, idihs);
	    break;
	  case F_DUMMY3OUT:
	    bERROR = 
	      calc_dum3out_param(ptype, &atype, 
				 &(plist[ftype].param[i]), atoms,
				 nrbond, bonds, nrang, angles, nridih, idihs);
	    break;
	  case F_DUMMY4FD:
	    bERROR = 
	      calc_dum4fd_param(ptype, &(plist[ftype].param[i]), atoms,
				nrbond, bonds, nrang, angles, nridih, idihs);
	    break;
	  default:
	    fatal_error(0,"Automatic parameter generation not supported "
			"for %s atom %d",
			interaction_function[ftype].longname,
			plist[ftype].param[i].AI+1);
	  } /* switch */
	  if (bERROR)
	    fatal_error(0,"Automatic parameter generation not supported "
			"for %s atom %d for this bonding configuration",
			interaction_function[ftype].longname,
			plist[ftype].param[i].AI+1);
	  sfree(bonds);
	  sfree(angles);
	  sfree(idihs);
	} /* if bSet */
      } /* for i */
      if (debug && plist[ftype].nr)
	fprintf(stderr,"Calculated parameters for %d out of %d %s atoms\n",
		nrset,plist[ftype].nr,interaction_function[ftype].longname);
    } /* if IF_DUMMY */
}

void set_dummies_ptype(bool bVerbose, t_idef *idef, t_atoms *atoms)
{
  int i,ftype;
  int nra,nrd,tp;
  t_iatom *ia,adum;
  
  if (bVerbose)
    fprintf(stderr,"Setting particle type to Dummy for dummy atoms\n");
  if (debug)
    fprintf(stderr,"checking %d functypes\n",F_NRE);
  for(ftype=0; (ftype<F_NRE); ftype++) {
    if (interaction_function[ftype].flags & IF_DUMMY) {
      nra    = interaction_function[ftype].nratoms;
      nrd    = idef->il[ftype].nr;
      ia     = idef->il[ftype].iatoms;
      
      if (debug && nrd)
	fprintf(stderr,"doing %d %s dummies\n",
		(nrd / (nra+1)),interaction_function[ftype].longname);
      
      for(i=0; (i<nrd); ) {
	tp   = ia[0];
	assert(ftype == idef->functype[tp]);
	
	/* The dummy atom */
	adum = ia[1];
	atoms->atom[adum].ptype=eptDummy;
	
	i  += nra+1;
	ia += nra+1;
      }
    }
  }
  
}
