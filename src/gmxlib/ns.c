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
 * Good ROcking Metal Altar for Chronical Sinners */
static char *SRCID_ns_c = "$Id$";

#include <math.h>
#include <string.h>
#include "sysstuff.h"
#include "assert.h"
#include "smalloc.h"
#include "macros.h"
#include "maths.h"
#include "vec.h"
#include "network.h"
#include "nsgrid.h"
#include "force.h"
#include "ns.h"
#include "pbc.h"
#include "names.h"
#include "fatal.h"
#include "nrnb.h"
#include "ns.h"
#include "fnbf.h"

#define MAX_CG 1024

typedef struct {
  int     ncg;
  int     nj;
  atom_id jcg[MAX_CG];
} t_ns_buf;

/* 
 *    E X C L U S I O N   H A N D L I N G
 */
typedef unsigned long t_excl;

#ifdef DEBUG
static void SETEXCL_(t_excl e[],atom_id i,atom_id j)
{   e[j] = e[j] | (1<<i); }
static void RMEXCL_(t_excl e[],atom_id i,atom_id j) 
{ e[j]=e[j] & ~(1<<i); }
static bool ISEXCL_(t_excl e[],atom_id i,atom_id j) 
{ return (bool)(e[j] & (1<<i)); }
static bool NOTEXCL_(t_excl e[],atom_id i,atom_id j)
{  return !(ISEXCL(e,i,j)); }
#else
#define SETEXCL(e,i,j) (e)[((atom_id) j)] |= (1<<((atom_id) i))
#define RMEXCL(e,i,j)  (e)[((atom_id) j)] &= (~(1<<((atom_id) i)))
#define ISEXCL(e,i,j)  (bool) ((e)[((atom_id) j)] & (1<<((atom_id) i)))
#define NOTEXCL(e,i,j) !(ISEXCL(e,i,j))
#endif

/************************************************
 *
 *  U T I L I T I E S    F O R    N S
 *
 ************************************************/

static int NLJ_INC = 16384;
static int NLI_INC = 2048;

void init_neighbor_list(FILE *log,t_forcerec *fr,int nn)
{
  fr->nn = nn;
  snew(fr->vdw,nn);
  snew(fr->free,nn);
  snew(fr->coul,nn);
  if (fr->bTwinRange) {
    snew(fr->vdw_lr,nn);
    snew(fr->free_lr,nn);
    snew(fr->coul_lr,nn);
  }
}

static void reset_neighbor_list(FILE *log,t_forcerec *fr,bool bLR)
{
  int i;
  
  for(i=0; (i < fr->nn); i++) {
    if (bLR) {
      if (fr->bTwinRange) {
	fr->vdw_lr[i].nri  = fr->vdw_lr[i].nrj  = 0;
	fr->coul_lr[i].nri = fr->coul_lr[i].nrj = 0;
	fr->free_lr[i].nri = fr->free_lr[i].nrj = 0;
      }
    }
    else {
      fr->vdw[i].nri  = fr->vdw[i].nrj  = 0;
      fr->coul[i].nri = fr->coul[i].nrj = 0;
      fr->free[i].nri = fr->free[i].nrj = 0;
    }
  }
}

static int count_neighbours(t_forcerec *fr,int *n_c,int *n_ljc,int *n_free)
{
  int i,nrj_c=0,nrj_ljc=0,nrj_free=0;
  
  for(i=0; (i < fr->nn); i++) {
    nrj_c    += fr->coul[i].nrj;
    nrj_ljc  += fr->vdw[i].nrj;
    nrj_free += fr->free[i].nrj;
    if (fr->bTwinRange) {
      nrj_c    += fr->coul_lr[i].nrj;
      nrj_ljc  += fr->vdw_lr[i].nrj;
      nrj_free += fr->free_lr[i].nrj;
    }
  }
  *n_c    = nrj_c;
  *n_ljc  = nrj_ljc;
  *n_free = nrj_free;
  
  return nrj_c+nrj_ljc+nrj_free;
}

static void new_i_nblist(FILE *log,t_nblist *nlist,int ftype,int i_atom,
			 bool bWater,int shift)
{
  t_nl_i *nli;
#ifdef _amb_
  int    ma;
#endif
    
  if (nlist->maxnri <= nlist->nri) {
#ifdef _amb_
    ma=memavail();
    fprintf(log,"Adding %5d I particles for nblist %s. Memavail = %d bytes\n",
	    NLI_INC,interaction_function[ftype].longname,ma);
#else
    fprintf(log,"Adding %5d I particles for nblist %s\n",NLI_INC,
	    interaction_function[ftype].longname);
#endif
    nlist->maxnri += NLI_INC;
    srenew(nlist->nl_i,nlist->maxnri);
  }
  if (nlist->maxnrj <= nlist->nrj + NLJ_INC-1) {
#ifdef _amb_
    ma=memavail();
    fprintf(log,"Adding %5d J particles for nblist %s. Memavail = %d bytes\n",
	    NLJ_INC,interaction_function[ftype].longname,ma);
#else
    fprintf(log,"Adding %5d J particles for nblist %s\n",NLJ_INC,
	    interaction_function[ftype].longname);
#endif
    nlist->maxnrj += NLJ_INC;
    srenew(nlist->nl_j,nlist->maxnrj);
  }
  
  /* Shortcut, copy the pointer... */
  nli          = &(nlist->nl_i[nlist->nri]);
  nli->i_atom  = i_atom;
  nli->shift   = shift;
  nli->nj      = 0;
  nli->j_index = nlist->nrj;
  nli->bWater  = bWater;
}

static void close_i_nblist(t_nblist *nlist) 
{
  if ((nlist->nl_i[nlist->nri].nj = 
       nlist->nrj - nlist->nl_i[nlist->nri].j_index) > 0)
    nlist->nri ++;
}

static void add_j_to_nblist(t_nblist nlist[],int gid,int j_atom)
{
  int nrj=nlist[gid].nrj;
  
  nlist[gid].nl_j[nrj] = j_atom;
  nlist[gid].nrj ++;
}

static void put_in_list(FILE *log,t_iparams ip[],int atnr,int nWater,
			int ngener,t_mdatoms *md,
			int icg,int nj,atom_id jjcg[],
			atom_id index[],atom_id a[],
			t_excl bExcl[],int shift,
			t_forcerec *fr,bool bLR,bool bCoulOnly)
{
  t_nblist  *vdw,*coul,*free;
  
  int 	    i,j,jcg,nit,igid,jgid,gid,ind_ij;
  atom_id   jj,jj0,jj1,i_atom,j_atom;
  int       i0,nicg;
  
  int       *type;
  ushort    *cENER;
  real      *charge;
  real      qi,qq,rlj;
  bool      bWater,bFree,bFreeJ,bNotEx,*bPert;
  
  /* Copy some pointers */
  charge  = md->chargeA;
  type    = md->typeA;
  cENER   = md->cENER;
  bPert   = md->bPerturbed;
  if (bLR) {
    vdw   = fr->vdw_lr;
    coul  = fr->coul_lr;
    free  = fr->free_lr;
  }
  else {
    vdw   = fr->vdw;
    coul  = fr->coul;
    free  = fr->free;
  }
  i0   = index[icg];
  nicg = index[icg+1]-i0;
  
  /* Check whether this molecule is a water molecule */
  bWater = (((type[a[i0]] == nWater) && (nicg == 3)) && 
	    (!bPert[a[i0]]) && (!bPert[a[i0+1]]) && (!bPert[a[i0+2]]));
  if (bWater)
    nicg = 1;

  /* Loop over the atoms in the i charge group */    
  for(i=0; (i<nicg); i++) {
    i_atom = a[i0+i];
    nit    = type[i_atom]*atnr;
    qi     = charge[i_atom];
    igid   = cENER[i_atom];
    bFree  = bPert[i_atom];
    
    /* Create new i_atom for each energy group */
    for(gid=0; (gid<fr->nn); gid++) {
      new_i_nblist(log,&vdw[gid], F_LJ,  i_atom,bWater,shift);
      new_i_nblist(log,&coul[gid],F_SR,  i_atom,bWater,shift);
      new_i_nblist(log,&free[gid],F_DVDL,i_atom,bWater,shift);
    }

    /* Loop over the j charge groups */
    for(j=0; (j<nj); j++) {
      jcg=jjcg[j];
      
      if (bWater && (jcg==icg))
	continue;

      /* The atoms in this j charge group run from jj0 to jj1 */
      jj0=index[jcg];
      jj1=index[jcg+1];
      
      /* Patch for the same charge group, to prevent double interactions within
       * one charge group.
       */
      if (icg == jcg)
	jj0 += i+1;
	
      /* Finally loop over the atoms in the j-charge group */	
      for(jj=jj0; (jj<jj1); jj++) {
	j_atom = a[jj];
	jgid   = cENER[j_atom];
	gid    = GID(igid,jgid,ngener);
	bFreeJ = bFree || bPert[j_atom];
	bNotEx = NOTEXCL(bExcl,i,j_atom);
	
	if (bNotEx) {
	  if (bFreeJ) 
	    add_j_to_nblist(free,gid,j_atom);
	  else if (bCoulOnly) 
	    /* This is done whether or  not bWater is set */
	    add_j_to_nblist(coul,gid,j_atom);
	  else if (bWater)
	    add_j_to_nblist(vdw,gid,j_atom);
	  else {
	    ind_ij = nit+type[j_atom];
	    rlj    = ip[ind_ij].lj.c6+ip[ind_ij].lj.c12;
	    
	    if (rlj) 
	      add_j_to_nblist(vdw,gid,j_atom);
	    else if (qi*charge[j_atom] != 0)
	      add_j_to_nblist(coul,gid,j_atom);
	  }
	}
      }
    }
    for(gid=0; (gid<fr->nn); gid++) {
      close_i_nblist(&coul[gid]);
      close_i_nblist(&vdw[gid]);
      close_i_nblist(&free[gid]);
    }
  }
}  

void setexcl(FILE *log,int nri,atom_id ia[],t_block *excl,bool b,
	     t_excl bexcl[])
{
  atom_id k;
  int     i,inr;
  
  if (b) {
    for(i=0; (i<nri); i++) {
      inr = ia[i];
      for(k=excl->index[inr]; (k<excl->index[inr+1]); k++) {
	SETEXCL(bexcl,i,excl->a[k]);
      }
    }
  }
  else {
    for(i=0; (i<nri); i++) {
      inr = ia[i];
      for(k=excl->index[inr]; (k<excl->index[inr+1]); k++) {
	RMEXCL(bexcl,i,excl->a[k]);
      }
    }
  }
}

int calc_naaj(FILE *log,int icg,int cgtot)
{
  int naaj;
  
  if ((cgtot % 2) == 1) {
    /* Odd number of charge groups, easy */
    naaj = 1+(cgtot/2);
  }
  else if ((cgtot % 4) == 0) {
    /* Multiple of four is hard */
    if (icg < cgtot/2) {
      if ((icg % 2) == 0)
	naaj=1+(cgtot/2);
      else
	naaj=cgtot/2;
    }
    else {
      if ((icg % 2) == 1)
	naaj=1+(cgtot/2);
      else
	naaj=cgtot/2;
    }
  }
  else {
    /* cgtot/2 = odd */
    if ((icg % 2) == 0)
      naaj=1+(cgtot/2);
    else
      naaj=cgtot/2;
  }
#ifdef DEBUG
  fprintf(log,"naaj=%d\n",naaj);
#endif
  return naaj;
}

/************************************************
 *
 *  S I M P L E      C O R E     S T U F F
 *
 ************************************************/
real calc_image_rect_test(rvec xi,rvec xj,rvec box_size,
			  rvec b_inv,t_ishift *shift)
{
  /* This code sucks! */
  real dd,hb,r2;
  ivec tt;
  int  m;
  
  r2 = 0;
  for(m=0; (m<DIM); m++) {
    hb = 0.5*box_size[m];
    dd = xj[m]-xi[m];
    if (dd < -hb)
      tt[m]=1;
    else if (dd >= hb)
      tt[m]=-1;
    else
      tt[m]=0;
    r2+=sqr(dd+tt[m]*box_size[m]);
  }
  
  *shift=XYZ2IS(tt[XX],tt[YY],tt[ZZ]);
  
  return r2;
}

real calc_image_rect(rvec xi,rvec xj,rvec box_size,
		     rvec b_inv,t_ishift *shift)
{
  const real h15=1.5;
  real ddx,ddy,ddz;
  real dx,dy,dz;
  real r2;
  int  tx,ty,tz;
  
  /* Compute diff vector */
  dx=xj[XX]-xi[XX];
  dy=xj[YY]-xi[YY];
  dz=xj[ZZ]-xi[ZZ];
  
  /* Perform NINT operation, using trunc operation, therefore
   * we first add 1.5 then subtract 1 again
   */
  tx=dx*b_inv[XX]+h15;
  ty=dy*b_inv[YY]+h15;
  tz=dz*b_inv[ZZ]+h15;
  tx--;
  ty--;
  tz--;
  
  /* Correct diff vector for translation */
  ddx=tx*box_size[XX]-dx;
  ddy=ty*box_size[YY]-dy;
  ddz=tz*box_size[ZZ]-dz;

  /* Distance squared */
  r2=(ddx*ddx)+(ddy*ddy)+(ddz*ddz);
  
  *shift=XYZ2IS(tx,ty,tz);
  
  return r2;
}

static void ns_inner_rect(FILE *log,rvec x[],int icg,int njcg,atom_id jcg[],
			  rvec box_size,rvec b_inv,real rcut2,
			  t_block *cgs,t_ns_buf ns_buf[])
{
  t_ishift shift;
  int      j,nrj;
  atom_id  cg_j;

  shift = CENTRAL;
  for(j=0; (j<njcg); j++) {
    cg_j   = jcg[j];
    nrj    = cgs->index[cg_j+1]-cgs->index[cg_j];
    
    if ((rcut2 == 0) || 
	(calc_image_rect(x[icg],x[cg_j],box_size,b_inv,&shift) < rcut2)) {
      if (ns_buf[shift].ncg >= MAX_CG) 
	fatal_error(0,"Too many charge groups (%d) in buffer",
		    ns_buf[shift].ncg);
      ns_buf[shift].jcg[ns_buf[shift].ncg++]=cg_j;
      ns_buf[shift].nj += nrj;
    }
  }
}

static int ns_simple_core(FILE *log,t_forcerec *fr,
			  rvec x[],matrix box,int ngener,
			  t_topology *top,t_groups *grps,
			  t_mdatoms *md,
			  rvec box_size,t_excl bexcl[],
			  t_ns_buf ns_buf[])
{
  static   atom_id  *aaj=NULL;
  int      naaj,k;
  real     rlist2;
  int      nsearch;
  int      icg,jcg;
  int      i0,nri;
  atom_id  *i_atoms;
  t_block  *cgs=&(top->blocks[ebCGS]);
  t_block  *excl=&(top->atoms.excl);
  rvec     b_inv;
  int      m;
  
  if (aaj==NULL) {
    snew(aaj,2*cgs->nr);
    for(jcg=0; (jcg<cgs->nr); jcg++) {
      aaj[jcg]=jcg;
      aaj[jcg+cgs->nr]=jcg;
    }
  }
  rlist2 = sqr(fr->rlist);

  for(m=0; (m<DIM); m++)
    b_inv[m]=divide(1.0,box_size[m]);

  nsearch=0;
  for (icg=fr->cg0; (icg<fr->hcg); icg++) {
    i0      = cgs->index[icg];
    nri     = cgs->index[icg+1]-i0;
    i_atoms = &(cgs->a[i0]);
    setexcl(log,nri,i_atoms,excl,TRUE,bexcl);
    
    naaj=calc_naaj(log,icg,cgs->nr);
    ns_inner_rect(log,fr->cg_cm,icg,naaj,&(aaj[icg]),
		  box_size,b_inv,rlist2,cgs,ns_buf);
    nsearch += naaj;
    
    for(k=0; (k<SHIFTS); k++) {
      if (ns_buf[k].ncg > 0) { 
	put_in_list(log,top->idef.iparams,top->idef.atnr,fr->nWater,
		    ngener,md,icg,ns_buf[k].ncg,ns_buf[k].jcg,
		    cgs->index,cgs->a,bexcl,k,fr,FALSE,FALSE);
	ns_buf[k].ncg=ns_buf[k].nj=0;
      }
    }
    setexcl(log,nri,i_atoms,excl,FALSE,bexcl);
  }

  return nsearch;
}

/************************************************
 *
 *    N S 5     G R I D     S T U F F
 *
 ************************************************/
static bool get_dx(int cx,int Nx,int tx,int delta,int *dx0,int *dx1)
{
  if (tx == 1) {
    if (cx >= delta)
      return TRUE;
    else {
      *dx0=Nx+cx-delta;
      *dx1=Nx-1;
    }
  } else if (tx == 0) {
    if (cx < delta)
      *dx0=0;
    else
      *dx0=cx-delta;
    if (cx < Nx-delta)
      *dx1=cx+delta;
    else
      *dx1=Nx-1;
  }
  else {
    if (cx < Nx-delta)
      return TRUE;
    else {
      *dx0=0;
      *dx1=cx+delta-Nx;
    }
  }
 
  return FALSE; 
}

#define sqr(x) ((x)*(x))
#define calc_dx2(XI,YI,ZI,y) (sqr(XI-y[XX])+sqr(YI-y[YY])+sqr(ZI-y[ZZ]))
#define calc_cyl_dx2(XI,YI,y) (sqr(XI-y[XX])+sqr(YI-y[YY]))
/****************************************************
 *
 *    F A S T   N E I G H B O R  S E A R C H I N G
 *
 *    Optimized neighboursearching routine using grid 
 *    of at least 5x5x5, see GROMACS manual
 *
 ****************************************************/

static void do_longrange(FILE *log,t_topology *top,t_forcerec *fr,
			 int ngener,t_mdatoms *md,int icg,int nlr,
			 atom_id lr[],t_excl bexcl[],int shift,
			 rvec x[],rvec box_size,t_nrnb *nrnb,
			 real lambda,real *dvdlambda,
			 t_groups *grps,bool bCoulOnly)
{
  int ftypes[] = { F_SR, F_DVDL, F_LJ };
#define NFTYPES asize(ftypes)
  int i;
  
  /* Put the long range particles in a list */
  put_in_list(log,top->idef.iparams,top->idef.atnr,fr->nWater,
	      ngener,md,icg,nlr,lr,top->blocks[ebCGS].index,
	      top->blocks[ebCGS].a,bexcl,shift,fr,TRUE,bCoulOnly);

  /* Evaluate the forces */
  for(i=0; (i<NFTYPES); i++) {
    do_fnbf(log,ftypes[i],fr,x,fr->flr,md,
	    grps->estat.ee[egLJLR],grps->estat.ee[egLR],box_size,
	    nrnb,lambda,dvdlambda,TRUE);
  }
  reset_neighbor_list(log,fr,TRUE);
}
  
static int ns5_core(FILE *log,t_forcerec *fr,int cg_index[],
		    matrix box,rvec box_size,int ngener,
		    t_topology *top,t_groups *grps,
		    t_grid *grid,rvec x[],t_excl bexcl[],
		    t_nrnb *nrnb,t_mdatoms *md,
		    real lambda,real *dvdlambda)
{
  static atom_id *nl_lr_ljc, *nl_lr_coul, *nl_sr=NULL;
  t_block *cgs=&(top->blocks[ebCGS]);
  int  tx,ty,tz,cx,cy,cz,dx,dy,dz;
  int  cj;
  int  dx0,dx1,dy0,dy1,dz0,dz1;
  int  Nx,Ny,Nz,delta,shift;
  real r2;
  int  nlr_ljc,nlr_coul,nsr,nns;
  int  j,nrj;
  int  icg,iicg,total_cg,i0,nri,naaj,min_icg,icg_naaj,jjcg,cgj0;
  atom_id *i_atoms;
  int  *grida,*gridnra,*gridind;
  rvec xi,*cgcm,*svec;
  real rs2,rvdw2,rcoul2,XI,YI,ZI;
  bool bCyl;
  
  total_cg = cgs->nr;
  rs2      = sqr(fr->rlist);
  rvdw2    = sqr(fr->rvdw);
  rcoul2   = sqr(fr->rcoulomb);
  
  if (nl_sr == NULL) {
    snew(nl_sr,MAX_CG);
    if (rvdw2 > rs2)
      snew(nl_lr_ljc,MAX_CG);
    if (rcoul2 > rvdw2)
      snew(nl_lr_coul,MAX_CG);
    if (debug)
      fprintf(debug,"ns5_core: rs2 = %g, rvdw2 = %g, rcoul2 = %g (nm^2)\n",
	      rs2,rvdw2,rcoul2);
  }
  
  /* Unpack arrays */
  cgcm    = fr->cg_cm;
  svec    = fr->shift_vec;
  Nx      = grid->nrx;
  Ny      = grid->nry;
  Nz      = grid->nrz;
  grida   = grid->a;
  gridind = grid->index;
  gridnra = grid->nra;
  delta   = grid->delta;
  nns     = 0;
  
  /* Loop over charge groups */
  for(iicg=fr->cg0; (iicg < fr->hcg); iicg++) {
    icg      = cg_index[iicg];
    /* Consistency check */
    if (icg != iicg)
      fatal_error(0,"icg = %d, iicg = %d, file %s, line %d",icg,iicg,__FILE__,
		  __LINE__);
    i0       = cgs->index[icg];
    nri      = cgs->index[icg+1]-i0;
    i_atoms  = &(cgs->a[i0]);

    /* Set the exclusions for the atoms in charge group icg using
     * a bitmask
     */    
    setexcl(log,nri,i_atoms,&top->atoms.excl,TRUE,bexcl);
    
    /* Compute the number of charge groups that fall within the control
     * of this one (icg)
     */
    naaj     = calc_naaj(log,icg,cgs->nr);
    icg_naaj = icg+naaj;
    min_icg  = icg_naaj-total_cg;
    
    /* Changed iicg to icg, DvdS 990115 
     * (but see consistency check above, DvdS 990330) 
     */
    ci2xyz(grid,icg,&cx,&cy,&cz);
#ifdef NS5DB
    fprintf(log,"icg=%5d, naaj=%5d, cx=%2d, cy=%2d, cz=%2d\n",
	    icg,naaj,cx,cy,cz);
#endif
    /* Loop over shift vectors in three dimensions */
    for (tx=-1; (tx<=1); tx++) {
      /* Calculate range of cells in X direction that have the shift tx */
      if (get_dx(cx,Nx,tx,delta,&dx0,&dx1))
	continue;
      for (ty=-1; (ty<=1); ty++) {
	/* Calculate range of cells in Y direction that have the shift ty */
	if (get_dx(cy,Ny,ty,delta,&dy0,&dy1))
	  continue;
	for (tz=-1; (tz<=1); tz++) {
	  /* Calculate range of cells in Z direction that have the shift tz */
	  if (get_dx(cz,Nz,tz,delta,&dz0,&dz1))
	    continue;

	  /* Get shift vector */	  
	  shift=XYZ2IS(tx,ty,tz);
#ifdef NS5DB
	  assert(shift >= 0);
	  assert(shift < SHIFTS);
#endif
	  XI       = cgcm[icg][XX]+svec[shift][XX];
	  YI       = cgcm[icg][YY]+svec[shift][YY];
	  ZI       = cgcm[icg][ZZ]+svec[shift][ZZ];
	  nsr      = 0;
	  nlr_ljc  = 0;
	  nlr_coul = 0;
#ifdef NS5DB
	  fprintf(log,"shift: %2d, dx0,1: %2d,%2d, dy0,1: %2d,%2d, dz0,1: %2d,%2d\n",
		  shift,dx0,dx1,dy0,dy1,dz0,dz1);
	  fprintf(log,"cgcm: %8.3f  %8.3f  %8.3f\n",cgcm[icg][XX],
		  cgcm[icg][YY],cgcm[icg][ZZ]);
	  fprintf(log,"xi:   %8.3f  %8.3f  %8.3f\n",XI,YI,ZI);
#endif
	  for (dx=dx0; (dx<=dx1); dx++) {
	    for (dy=dy0; (dy<=dy1); dy++) {
	      for (dz=dz0; (dz<=dz1); dz++) {
		/* Find the grid-cell cj in which possible neighbours are */
		cj   = xyz2ci(Ny,Nz,dx,dy,dz);
		
		/* Check out how many cgs (nrj) there in this cell */
		nrj  = gridnra[cj];
		
		/* Find the offset in the cg list */
		cgj0 = gridind[cj];
		
		/* Loop over cgs */
		for (j=0; (j<nrj); j++) {
		  jjcg = grida[cgj0+j];
		  
		  /* check whether this guy is in range! */
		  if (((jjcg >= icg) && (jjcg < icg_naaj)) ||
		      ((jjcg < min_icg))) {
		    r2=calc_dx2(XI,YI,ZI,cgcm[jjcg]);
		    if (r2 < rs2) {
		      if (nsr >= MAX_CG) {
			put_in_list(log,top->idef.iparams,
				    top->idef.atnr,fr->nWater,
				    ngener,md,icg,nsr,nl_sr,
				    cgs->index,cgs->a,bexcl,
				    shift,fr,FALSE,FALSE);
			nsr=0;
		      }
		      nl_sr[nsr++]=jjcg;
		    } else if (r2 < rvdw2) {
		      if (nlr_ljc >= MAX_CG) {
			do_longrange(log,top,fr,ngener,md,icg,nlr_ljc,
				     nl_lr_ljc,bexcl,shift,x,box_size,nrnb,
				     lambda,dvdlambda,grps,FALSE);
			nlr_ljc=0;
		      }
		      nl_lr_ljc[nlr_ljc++]=jjcg;
		    } else if (r2 < rcoul2) {
		      if (nlr_coul >= MAX_CG) {
			do_longrange(log,top,fr,ngener,md,icg,nlr_coul,
				     nl_lr_coul,bexcl,shift,x,box_size,nrnb,
				     lambda,dvdlambda,grps,TRUE);
			nlr_coul=0;
		      }
		      nl_lr_coul[nlr_coul++]=jjcg;
		    }
		    nns++;
		  }
		}
	      }
	    }
	  }
	  /* CHECK whether there is anything left in the buffers */
	  if (nsr > 0)
	    put_in_list(log,top->idef.iparams,top->idef.atnr,fr->nWater,
			ngener,md,icg,nsr,nl_sr,cgs->index,cgs->a,bexcl,
			shift,fr,FALSE,FALSE);
			
	  if (nlr_ljc > 0) 
	    do_longrange(log,top,fr,ngener,md,icg,nlr_ljc,
			 nl_lr_ljc,bexcl,shift,x,box_size,nrnb,
			 lambda,dvdlambda,grps,FALSE);
	  
	  if (nlr_coul > 0) 
	    do_longrange(log,top,fr,ngener,md,icg,nlr_coul,
			 nl_lr_coul,bexcl,shift,x,box_size,nrnb,
			 lambda,dvdlambda,grps,TRUE);
	  
	}
      }
    }
    setexcl(log,nri,i_atoms,&top->atoms.excl,FALSE,bexcl);
  }
  return nns;
}

static rvec *sptr;
static int  sdim;

static int  rv_comp(const void *a,const void *b)
{
  int  ia = *(int *)a;
  int  ib = *(int *)b;
  real diff;
  
  diff = sptr[ia][sdim] - sptr[ib][sdim];
  if (diff < 0)
    return -1;
  else if (diff == 0)
    return 0;
  else
    return 1;
}

static void sort_charge_groups(FILE *log,
			       t_commrec *cr,int cg_index[],int slab_index[],
			       rvec cg_cm[],int Dimension)
{
  int i,nrcg,cgind;
  
  nrcg = slab_index[cr->nprocs];
  sptr = cg_cm;
  sdim = Dimension;
  qsort(cg_index,nrcg,sizeof(cg_index[0]),rv_comp);
  
  if (debug) {
    fprintf(debug,"Just sorted the cg_cm array on dimension %d\n",Dimension);
    fprintf(debug,"Index:  Coordinates of cg_cm\n");
    for(i=0; (i<nrcg); i++) {
      cgind = cg_index[i];
      fprintf(debug,"%8d%10.3f%10.3f%10.3f\n",cgind,
	      cg_cm[cgind][XX],cg_cm[cgind][YY],cg_cm[cgind][ZZ]);
    }
  }
  sptr = NULL;
  sdim = -1;
}

int search_neighbours(FILE *log,t_forcerec *fr,
		      rvec x[],matrix box,
		      t_topology *top,t_groups *grps,
		      t_commrec *cr,t_nsborder *nsb,
		      t_nrnb *nrnb,t_mdatoms *md,
		      real lambda,real *dvdlambda)
{
  static   bool        bFirst=TRUE;
  static   t_grid      *grid=NULL;
  static   t_excl      *bexcl;
  static   t_ns_buf    *ns_buf=NULL;
  static   int         *cg_index=NULL,*slab_index=NULL;
  static   bool        bSwitched=FALSE;
  
  t_block  *cgs=&(top->blocks[ebCGS]);
  rvec     box_size;
  int      i,m,ngener;

  int      nsearch;
  bool     bGrid;
  char     *ptr;
  
  /* Set some local variables */
  bGrid=fr->bGrid;
  ngener=top->atoms.grps[egcENER].nr;
  
  for(m=0; (m<DIM); m++)
    box_size[m]=box[m][m];

  /* First time initiation of arrays etc. */  
  if (bFirst) {
    int icg,nr_in_cg,maxcg;
    
    /* Compute largest charge groups size (# atoms) */
    nr_in_cg=1;
    for (icg=0; (icg < cgs->nr); icg++)
      nr_in_cg=max(nr_in_cg,(int)(cgs->index[icg+1]-cgs->index[icg]));

    /* Verify whether largest charge group is <= max cg.
     * This is determined by the type of the local exclusion type 
     * Exclusions are stored in bits. (If the type is not large
     * enough, enlarge it, unsigned char -> unsigned short -> unsigned long)
     */
    maxcg=sizeof(t_excl)*8;
    if (nr_in_cg > maxcg)
      fatal_error(0,"Max #atoms in a charge group: %d > %d\n",
		  nr_in_cg,maxcg);
      
    snew(bexcl,cgs->nra);
    debug_gmx();

    if ((ptr=getenv("NLIST")) != NULL) {
      sscanf(ptr,"%d",&NLJ_INC);
      NLI_INC = NLJ_INC/8;
      
      fprintf(log,"%s: I will increment I-lists by %d and J-lists by %d\n",
	      __FILE__,NLI_INC,NLJ_INC);
    }

    /* Check whether we have to do domain decomposition,
     * if so set local variables for the charge group index and the
     * slab index.
     */
    if (fr->bDomDecomp) {
      snew(slab_index,cr->nprocs+1);
      for(i=0; (i<=cr->nprocs); i++)
	slab_index[i] = i*((real) cgs->nr/((real) cr->nprocs));
      fr->cg0 = slab_index[cr->pid];
      fr->hcg = slab_index[cr->pid+1] - fr->cg0;
      if (debug)
	fprintf(debug,"Will use DOMAIN DECOMPOSITION, "
		"from charge group index %d to %d on cpu %d\n",
		fr->cg0,fr->cg0+fr->hcg,cr->pid);
    }
    snew(cg_index,cgs->nr+1);
    for(i=0; (i<=cgs->nr);  i++)
      cg_index[i] = i;
    
    if (bGrid) {
      snew(grid,1);
      init_grid(log,grid,fr->ndelta,box,fr->rlistlong,cgs->nr);
    }
    bFirst=FALSE;
  }
  debug_gmx();
  
  /* Reset the neighbourlists */
  reset_neighbor_list(log,fr,FALSE);
  
  if (bGrid) {
    grid_first(log,grid,box,fr->rlistlong);
    /* Check if box is big enough to do grid searching... */
    if ( !( (grid->nrx >= 2*grid->delta+1) && 
	    (grid->nry >= 2*grid->delta+1) && 
	    (grid->nrz >= 2*grid->delta+1) ) ) {
      if (!bSwitched)
	fprintf(log,"WARNING: Box too small for grid-search, "
		"switching to simple neighboursearch.\n");
      if (fr->bTwinRange)
	fatal_error(0,"TWIN-RANGE cut-off with Simple "
		    "neighboursearching not implemented.\n"
		    "Use grid neighboursearching, and make (rlong < 0.4 box)");
      bGrid=FALSE;
      bSwitched=TRUE;
    } else {
      if (bSwitched)
	fprintf(log,"WARNING: Box large enough again for grid-search\n");
      bSwitched=FALSE;
    }
  }
  debug_gmx();
  
  if (bGrid) {
    /* Don't know why this all is... (DvdS 3/99) */
#ifndef SEGV
    int start = 0;
    int end   = cgs->nr;
#else
    int start = fr->cg0;
    int end   = (cgs->nr+1)/2;
#endif

    if (fr->bDomDecomp) {
      sort_charge_groups(log,cr,cg_index,slab_index,fr->cg_cm,fr->Dimension);
    }
    debug_gmx();
    
    fill_grid(log,fr->bDomDecomp,cg_index,
	      grid,box,cgs->nr,fr->cg0,fr->hcg,fr->cg_cm);
    debug_gmx();

    if (PAR(cr))
      mv_grid(cr,fr->bDomDecomp,cg_index,grid,nsb->workload);
    debug_gmx();
      
    calc_elemnr(log,fr->bDomDecomp,cg_index,grid,start,end,cgs->nr);
    calc_ptrs(grid);
    grid_last(log,fr->bDomDecomp,cg_index,grid,start,end,cgs->nr);

    if (debug) {
      check_grid(debug,grid);
      print_grid(debug,grid,fr->bDomDecomp,cg_index);
    }
  }
  debug_gmx();
  
  /* Do the core! */
  if (bGrid)
    nsearch = ns5_core(log,fr,cg_index,box,box_size,ngener,top,grps,
		       grid,x,bexcl,nrnb,md,lambda,dvdlambda);
  else {
    /* Only allocate this when necessary, saves 100 kb */
    if (ns_buf == NULL)
      snew(ns_buf,SHIFTS);
    nsearch = ns_simple_core(log,fr,x,box,ngener,top,grps,md,box_size,
			     bexcl,ns_buf);
  }
  debug_gmx();
  
#ifdef DEBUG
  pr_nsblock(log);
#endif

  inc_nrnb(nrnb,eNR_NS,nsearch);
  inc_nrnb(nrnb,eNR_LR,fr->nlr);

  if (debug) {
    int n_c,n_ljc,n_free,ntot;
    
    ntot = count_neighbours(fr,&n_c,&n_ljc,&n_free);
    fprintf(debug,"%d neighbours\n",ntot);
  }
  return nsearch;
}

static void put_in_nblr(FILE *log,int nrj,atom_id j_atoms[],
			t_nblist_lr *nlr)
{
  int j,n;
  
  n=nlr->nj;
  if (n+nrj >= MAXNB_LR) 
    fatal_error(0,"Too many long range neighbours. Increase MAXNB_LR");
  
  for(j=0; (j<nrj); j++) {
    nlr->nlj[n++]=j_atoms[j];
  }
  nlr->nj+=nrj;    
}

