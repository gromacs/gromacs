#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <stdio.h>
#include <time.h>
#include <math.h>
#include <string.h>
#include "typedefs.h"
#include "smalloc.h"
#include "vec.h"
#include "domdec.h"
#include "nrnb.h"
#include "pbc.h"
#include "constr.h"
#include "mdatoms.h"
#include "names.h"

#ifdef GMX_MPI
#include <mpi.h>
#endif


/* Code only works for one moleculetype consisting of one charge group */
/* #define ONE_MOLTYPE_ONE_CG */

#ifdef ONE_MOLTYPE_ONE_CG
static int natoms_global;
#endif


#define CG_ALLOC_SIZE     1000
#define EXCLS_ALLOC_SIZE  1000
#define IATOM_ALLOC_SIZE  1000

static const int cell_perm[3][4] = { {0,0,0,0},{1,0,0,0},{3,0,1,2} };

#define dd_c3n 8
static const ivec dd_c3[dd_c3n] =
  {{0,0,0},{0,0,1},{0,1,1},{0,1,0},{1,1,0},{1,0,0},{1,0,1},{1,1,1}};
#define dd_cp3n 4
static const ivec dd_cp3[dd_cp3n] =
  {{0,0,8},{1,3,6},{2,5,6},{3,5,7}};

#define dd_c2n 4
static const ivec dd_c2[dd_c2n] =
  {{0,0,0},{0,1,0},{1,1,0},{1,0,0}};
#define dd_cp2n 2
static const ivec dd_cp2[dd_cp2n] =
  {{0,0,4},{1,3,4}};

#define dd_c1n 2
static const ivec dd_c1[dd_c1n] =
  {{0,0,0},{1,0,0}};
#define dd_cp1n 1
static const ivec dd_cp1[dd_cp1n] =
  {{0,0,2}};

/*
#define dd_index(n,i) ((((i)[ZZ]*(n)[YY] + (i)[YY])*(n)[XX]) + (i)[XX])

static void index2xyz(ivec nc,int ind,ivec xyz)
{
  xyz[XX] = ind % nc[XX];
  xyz[YY] = (ind / nc[XX]) % nc[YY];
  xyz[ZZ] = ind / (nc[YY]*nc[XX]);
}
*/

/* This order is required to minimize the coordinate communication in PME
 * which uses decomposition in the x direction.
 */
#define dd_index(n,i) ((((i)[XX]*(n)[YY] + (i)[YY])*(n)[ZZ]) + (i)[ZZ])

void gmx_ddindex2xyz(ivec nc,int ind,ivec xyz)
{
  xyz[XX] = ind / (nc[YY]*nc[ZZ]);
  xyz[YY] = (ind / nc[ZZ]) % nc[YY];
  xyz[ZZ] = ind % nc[ZZ];
}

int glatnr(gmx_domdec_t *dd,int i)
{
  int atnr;

  if (dd == NULL) {
    atnr = i + 1;
  } else {
    if (i >= dd->nat_tot_con)
      gmx_fatal(FARGS,"glatnr called with %d, which is larger than the local number of atoms (%d)",i,dd->nat_tot_con);
    atnr = dd->gatindex[i] + 1;
  }

  return atnr;
}

void dd_get_ns_ranges(gmx_domdec_t *dd,int icg,
		      int *jcg0,int *jcg1,ivec shift0,ivec shift1)
{
  int i;

  i = 0;
  while (icg >= dd->icell[i].cg1)
    i++;

  if (i == 0)
    *jcg0 = icg;
  else if (i < dd->nicell)
    *jcg0 = dd->icell[i].jcg0;
  else
    gmx_fatal(FARGS,"DD icg %d out of range: icell (%d) >= nicell (%d)",
	      icg,i,dd->nicell);

  *jcg1 = dd->icell[i].jcg1;
  copy_ivec(dd->icell[i].shift0,shift0);
  copy_ivec(dd->icell[i].shift1,shift1);
}

void dd_sendrecv_int(const gmx_domdec_t *dd,
		     int dim,int direction,
		     int *buf_s,int n_s,
		     int *buf_r,int n_r)
{
#ifdef GMX_MPI
  int rank_s,rank_r;
  MPI_Status stat;

  rank_s = dd->neighbor[dim][direction==ddForward ? 0 : 1];
  rank_r = dd->neighbor[dim][direction==ddForward ? 1 : 0];

  MPI_Sendrecv(buf_s,n_s*sizeof(int),MPI_BYTE,rank_s,0,
	       buf_r,n_r*sizeof(int),MPI_BYTE,rank_r,0,
	       dd->all,&stat);
#endif
}

void dd_sendrecv_rvec(const gmx_domdec_t *dd,
		      int dim,int direction,
		      rvec *buf_s,int n_s,
		      rvec *buf_r,int n_r)
{
#ifdef GMX_MPI
  int rank_s,rank_r;
  MPI_Status stat;

  rank_s = dd->neighbor[dim][direction==ddForward ? 0 : 1];
  rank_r = dd->neighbor[dim][direction==ddForward ? 1 : 0];

  MPI_Sendrecv(buf_s[0],n_s*sizeof(rvec),MPI_BYTE,rank_s,0,
	       buf_r[0],n_r*sizeof(rvec),MPI_BYTE,rank_r,0,
	       dd->all,&stat);
#endif
}

void dd_move_x(gmx_domdec_t *dd,rvec x[],rvec buf[])
{
  int  ncell,nat_tot,n,dim,i,j;
  int  *index,*cgindex;
  gmx_domdec_comm_t *cc;
  
  cgindex = dd->cgindex;

  ncell = 1;
  nat_tot = dd->nat_home;
  for(dim=0; dim<dd->ndim; dim++) {
    cc = &dd->comm[dim];
    index = cc->index;
    n = 0;
    for(i=0; i<cc->nsend[ncell]; i++) {
      for(j=cgindex[index[i]]; j<cgindex[index[i]+1]; j++) {
	copy_rvec(x[j],buf[n]);
	n++;
      }
    }
    /* Send and receive the coordinates */
    dd_sendrecv_rvec(dd, dim, ddBackward,
		     buf,       cc->nsend[ncell+1],
		     x+nat_tot, cc->nrecv[ncell+1]);
    nat_tot += cc->nrecv[ncell+1];
    ncell += ncell;
  }
}

void dd_move_f(gmx_domdec_t *dd,rvec f[],rvec buf[])
{
  int  ncell,nat_tot,n,dim,i,j;
  int  *index,*cgindex;
  gmx_domdec_comm_t *cc;
  
  cgindex = dd->cgindex;

  ncell = 1;
  nat_tot = dd->nat_home;
  n = 0;
  ncell = dd->ncell/2;
  nat_tot = dd->nat_tot;
  for(dim=dd->ndim-1; dim>=0; dim--) {
    cc = &dd->comm[dim];
    nat_tot -= cc->nrecv[ncell+1];
    /* Communicate the forces */
    dd_sendrecv_rvec(dd,dim,ddForward,
		     f+nat_tot, cc->nrecv[ncell+1],
		     buf,       cc->nsend[ncell+1]);
    index = cc->index;
    /* Add the received forces */
    n = 0;
    for(i=0; i<cc->nsend[ncell]; i++) {
      for(j=cgindex[index[i]]; j<cgindex[index[i]+1]; j++) {
	rvec_inc(f[j],buf[n]);
	n++;
      }
    }
    ncell /= 2;
  }
}

static void dd_collect_cg(gmx_domdec_t *dd)
{
  int buf2[2],i;

  buf2[0] = dd->ncg_home;
  buf2[1] = dd->nat_home;

  /* Collect the charge group and atom counts on the master */
#ifdef GMX_MPI
  MPI_Gather(buf2,2*sizeof(int),MPI_BYTE,
	     dd->ma.ibuf,2*sizeof(int),MPI_BYTE,
	     DDMASTERRANK(dd),dd->all);
#endif
  
  if (DDMASTER(dd)) {
    dd->ma.index[0] = 0;
    for(i=0; i<dd->nnodes; i++) {
      dd->ma.ncg[i] = dd->ma.ibuf[2*i];
      dd->ma.nat[i] = dd->ma.ibuf[2*i+1];
      dd->ma.index[i+1] = dd->ma.index[i] + dd->ma.ncg[i];
      
    }
    /* Make byte counts and indices */
    for(i=0; i<dd->nnodes; i++) {
      dd->ma.ibuf[i] = dd->ma.ncg[i]*sizeof(int);
      dd->ma.ibuf[dd->nnodes+i] = dd->ma.index[i]*sizeof(int);
    }
    if (debug) {
      fprintf(debug,"Initial charge group distribution: ");
      for(i=0; i<dd->nnodes; i++)
	fprintf(debug," %d",dd->ma.ncg[i]);
      fprintf(debug,"\n");
    }
  }
  
  /* Collect the charge group indices the master */
#ifdef GMX_MPI
  MPI_Gatherv(dd->index_gl,dd->ncg_home*sizeof(int),MPI_BYTE,
	      dd->ma.cg,dd->ma.ibuf,dd->ma.ibuf+dd->nnodes,MPI_BYTE,
	      DDMASTERRANK(dd),dd->all);
#endif

  dd->bMasterHasAllCG = TRUE;
}

void dd_collect_vec(gmx_domdec_t *dd,t_block *cgs,rvec *lv,rvec *v)
{
  int  n,i,c,a;
  rvec *buf;

  if (!dd->bMasterHasAllCG)
    dd_collect_cg(dd);

  if (!DDMASTER(dd)) {
#ifdef GMX_MPI
    MPI_Send(lv,dd->nat_home*sizeof(rvec),MPI_BYTE,DDMASTERRANK(dd),
	     dd->nodeid,dd->all);
#endif
  } else {
    /* Copy the master coordinates to the global array */
    n = 0;
    a = 0;
    for(i=dd->ma.index[n]; i<dd->ma.index[n+1]; i++)
      for(c=cgs->index[dd->ma.cg[i]]; c<cgs->index[dd->ma.cg[i]+1]; c++)
	copy_rvec(lv[a++],v[c]);

    /* Use the unused part of lv as a temporary buffer */
    buf = lv + dd->ma.nat[n];

    for(n=0; n<dd->nnodes; n++) {
      if (n != dd->nodeid) {
#ifdef GMX_MPI
	MPI_Recv(buf,dd->ma.nat[n]*sizeof(rvec),MPI_BYTE,DDRANK(dd,n),
		 n,dd->all,MPI_STATUS_IGNORE);
#endif
	a = 0;
	for(i=dd->ma.index[n]; i<dd->ma.index[n+1]; i++)
	  for(c=cgs->index[dd->ma.cg[i]]; c<cgs->index[dd->ma.cg[i]+1]; c++)
	    copy_rvec(buf[a++],v[c]);
      }
    }
  }
}

void dd_collect_state(gmx_domdec_t *dd,t_block *cgs,
		      t_state *state_local,t_state *state)
{
  dd_collect_vec(dd,cgs,state_local->x,state->x);
  if (state->v)
    dd_collect_vec(dd,cgs,state_local->v,state->v);
  if (state->sd_X)
    dd_collect_vec(dd,cgs,state_local->sd_X,state->sd_X);
}

static void dd_distribute_vec(gmx_domdec_t *dd,t_block *cgs,rvec *v,rvec *lv)
{
  int n,i,c,a;

  if (DDMASTER(dd)) {
    for(n=0; n<dd->nnodes; n++) {
      if (n != dd->nodeid) {
	/* Use lv as a temporary buffer */
	a = 0;
	for(i=dd->ma.index[n]; i<dd->ma.index[n+1]; i++)
	  for(c=cgs->index[dd->ma.cg[i]]; c<cgs->index[dd->ma.cg[i]+1]; c++)
	    copy_rvec(v[c],lv[a++]);
	if (a != dd->ma.nat[n])
	  gmx_fatal(FARGS,"Internal error a (%d) != nat (%d)",a,dd->ma.nat[n]);
#ifdef GMX_MPI
	MPI_Send(lv,dd->ma.nat[n]*sizeof(rvec),MPI_BYTE,
		 DDRANK(dd,n),n,dd->all);
#endif
      }
    }
    n = 0;
    a = 0;
    for(i=dd->ma.index[n]; i<dd->ma.index[n+1]; i++)
      for(c=cgs->index[dd->ma.cg[i]]; c<cgs->index[dd->ma.cg[i]+1]; c++)
	copy_rvec(v[c],lv[a++]);
  } else {
#ifdef GMX_MPI
    MPI_Recv(lv,dd->nat_home*sizeof(rvec),MPI_BYTE,DDMASTERRANK(dd),
	     MPI_ANY_TAG,dd->all,MPI_STATUS_IGNORE);
#endif
  }
}

void dd_distribute_state(gmx_domdec_t *dd,t_block *cgs,
			 t_state *state,t_state *state_local)
{
  dd_distribute_vec(dd,cgs,state->x,state_local->x);
  if (state->v)
    dd_distribute_vec(dd,cgs,state->v,state_local->v);
  if (state->sd_X)
    dd_distribute_vec(dd,cgs,state->sd_X,state_local->sd_X);
}

static void make_dd_indices(gmx_domdec_t *dd,t_block *gcgs,int cg_start,
			    t_forcerec *fr)
{
  int cell,cg0,cg,cg_gl,a,a_gl;
  int *cell_ncg,*index_gl,*cgindex,*gatindex;
  gmx_ga2la_t *ga2la;
  bool bMakeSolventType;

  if (dd->nat_tot > dd->gatindex_nalloc) {
    dd->gatindex_nalloc = over_alloc(dd->nat_tot);
    srenew(dd->gatindex,dd->gatindex_nalloc);
  }
  
  if (fr->solvent_opt == esolNO) {
    /* Since all entries are identical, we can use the global array */
    fr->solvent_type = fr->solvent_type_global;
    bMakeSolventType = FALSE;
  } else {
    if (dd->ncg_tot > fr->solvent_type_nalloc) {
      fr->solvent_type_nalloc = over_alloc(dd->ncg_tot);
      srenew(fr->solvent_type,fr->solvent_type_nalloc);
    }
    bMakeSolventType = TRUE;
  }
  

  cell_ncg   = dd->ncg_cell;
  index_gl   = dd->index_gl;
  cgindex    = gcgs->index;
  gatindex   = dd->gatindex;

  /* Make the local to global and global to local atom index */
  a = dd->cgindex[cg_start];
  for(cell=0; cell<dd->ncell; cell++) {
    if (cell == 0)
      cg0 = cg_start;
    else
      cg0 = cell_ncg[cell];
    for(cg=cg0; cg<cell_ncg[cell+1]; cg++) {
      cg_gl = index_gl[cg];
      for(a_gl=cgindex[cg_gl]; a_gl<cgindex[cg_gl+1]; a_gl++) {
	gatindex[a] = a_gl;
	ga2la = &dd->ga2la[a_gl];
	ga2la->cell = cell;
	ga2la->a    = a++;
      }
      if (bMakeSolventType)
	fr->solvent_type[cg] = fr->solvent_type_global[cg_gl];
    }
  }
}

static void clear_dd_indices(gmx_domdec_t *dd,int a_start)
{
  int i;

  /* Clear the indices without have to go over all the atoms in the system */
  for(i=a_start; i<dd->nat_tot; i++)
    dd->ga2la[dd->gatindex[i]].cell = -1;
  if (dd->constraints)
    clear_local_constraint_indices(dd);
}

static void distribute_cg(FILE *fplog,matrix box,t_block *cgs,rvec pos[],
			  gmx_domdec_t *dd)
{
  int **tmp_ind=NULL,*tmp_nalloc=NULL;
  int  i,icg,k,k0,k1,d;
  rvec g_inv,cg_cm;
  ivec ind;
  real nrcg,inv_ncg;
  atom_id *cgindex;

  if (tmp_ind == NULL) {
    snew(tmp_nalloc,dd->nnodes);
    snew(tmp_ind,dd->nnodes);
    for(i=0; i<dd->nnodes; i++) {
      tmp_nalloc[i] = (cgs->nr/CG_ALLOC_SIZE + 2)*CG_ALLOC_SIZE;
      snew(tmp_ind[i],tmp_nalloc[i]);
    }
  }

  /* Clear the count */
  for(i=0; i<dd->nnodes; i++) {
    dd->ma.ncg[i] = 0;
    dd->ma.nat[i] = 0;
  }
  
  for(d=0; (d<DIM); d++)
    g_inv[d] = divide(dd->nc[d],box[d][d]);

  cgindex = cgs->index;
  
  /* Compute the center of geometry for all charge groups */
  for(icg=0; icg<cgs->nr; icg++) {
    k0      = cgindex[icg];
    k1      = cgindex[icg+1];
    nrcg    = k1 - k0;
    if (nrcg == 1) {
      copy_rvec(pos[k0],cg_cm);
    }
    else {
      inv_ncg = 1.0/nrcg;
      
      clear_rvec(cg_cm);
      for(k=k0; (k<k1); k++)
	rvec_inc(cg_cm,pos[k]);
      for(d=0; (d<DIM); d++)
	cg_cm[d] = inv_ncg*cg_cm[d];
    }
    /* Put the charge group in the box and determine the cell index */
    for(d=DIM-1; d>=0; d--) {
      while(cg_cm[d] >= box[d][d]) {
	rvec_dec(cg_cm,box[d]);
	for(k=k0; (k<k1); k++)
	  rvec_dec(pos[k],box[d]);
      }
      while(cg_cm[d] < 0) {
	rvec_inc(cg_cm,box[d]);
	for(k=k0; (k<k1); k++)
	  rvec_inc(pos[k],box[d]);
      }
      ind[d] = cg_cm[d]*g_inv[d];
      if (ind[d] >= dd->nc[d])
	ind[d] = dd->nc[d] - 1;
    }
    i = dd_index(dd->nc,ind);
    if (dd->ma.ncg[i] == tmp_nalloc[i]) {
      tmp_nalloc[i] += CG_ALLOC_SIZE;
      srenew(tmp_ind[i],tmp_nalloc[i]);
    }
    tmp_ind[i][dd->ma.ncg[i]] = icg;
    dd->ma.ncg[i]++;
    dd->ma.nat[i] += cgindex[icg+1] - cgindex[icg];
  }
  
  k1 = 0;
  for(i=0; i<dd->nnodes; i++) {
    dd->ma.index[i] = k1;
    for(k=0; k<dd->ma.ncg[i]; k++)
      dd->ma.cg[k1++] = tmp_ind[i][k];
  }
  dd->ma.index[dd->nnodes] = k1;

  sfree(tmp_nalloc);
  for(i=0; i<dd->nnodes; i++)
    sfree(tmp_ind[i]);
  sfree(tmp_ind);

  fprintf(fplog,"Charge group distribution:");
  for(i=0; i<dd->nnodes; i++)
    fprintf(fplog," %d",dd->ma.ncg[i]);
  fprintf(fplog,"\n");
}

static void get_cg_distribution(FILE *fplog,gmx_domdec_t *dd,
				t_block *cgs,matrix box,rvec pos[])
{
  int i,buf2[2],cg_gl;

  clear_dd_indices(dd,0);

  if (DDMASTER(dd)) {
    if (dd->ma.ncg == NULL) {
      snew(dd->ma.ncg,dd->nnodes);
      snew(dd->ma.index,dd->nnodes+1);
      snew(dd->ma.cg,cgs->nr);
      snew(dd->ma.nat,dd->nnodes);
      snew(dd->ma.ibuf,dd->nnodes*2);
    }      
    
    distribute_cg(fplog,box,cgs,pos,dd);
    for(i=0; i<dd->nnodes; i++) {
      dd->ma.ibuf[2*i]   = dd->ma.ncg[i];
      dd->ma.ibuf[2*i+1] = dd->ma.nat[i];
    }
  }
#ifdef GMX_MPI
  MPI_Scatter(dd->ma.ibuf,2*sizeof(int),MPI_BYTE,
	      buf2,2*sizeof(int),MPI_BYTE,
	      DDMASTERRANK(dd),dd->all);
#endif

  dd->ncg_home = buf2[0];
  dd->nat_home = buf2[1];
  if (dd->ncg_home > dd->cg_nalloc || dd->cg_nalloc == 0) {
    dd->cg_nalloc = over_alloc(dd->ncg_home);
    srenew(dd->index_gl,dd->cg_nalloc);
    srenew(dd->cgindex,dd->cg_nalloc+1);
  }
  if (DDMASTER(dd)) {
    for(i=0; i<dd->nnodes; i++) {
      dd->ma.ibuf[i] = dd->ma.ncg[i]*sizeof(int);
      dd->ma.ibuf[dd->nnodes+i] = dd->ma.index[i]*sizeof(int);
    }
  }

#ifdef GMX_MPI
  MPI_Scatterv(dd->ma.cg,dd->ma.ibuf,dd->ma.ibuf+dd->nnodes,MPI_BYTE,
	       dd->index_gl,dd->ncg_home*sizeof(int),MPI_BYTE,
	       DDMASTERRANK(dd),dd->all);
#endif

  /* Determine the home charge group sizes */
  dd->cgindex[0] = 0;
  for(i=0; i<dd->ncg_home; i++) {
    cg_gl = dd->index_gl[i];
    dd->cgindex[i+1] =
      dd->cgindex[i] + cgs->index[cg_gl+1] - cgs->index[cg_gl];
  }

  if (debug) {
    fprintf(debug,"Home charge groups:\n");
    for(i=0; i<dd->ncg_home; i++) {
      fprintf(debug," %d",dd->index_gl[i]);
      if (i % 10 == 9) 
	fprintf(debug,"\n");
    }
    fprintf(debug,"\n");
  }

  dd->bMasterHasAllCG = TRUE;
}

#define DD_MAXNBCELL 27

static int compact_and_copy_vec(int ncg,int *cell,int local,int *cgindex,
				int nnbc,int *buf_pos,rvec *vec,rvec *buf)
{
  int c,icg,i;
  int local_pos;
  

  local_pos = 0;

  for(icg=0; icg<ncg; icg++) {
    c = cell[icg];
    if (c == local) {
      /* Locally compact the array */
      for(i=cgindex[icg]; i<cgindex[icg+1]; i++)
	copy_rvec(vec[i],vec[local_pos++]);
    } else {
      /* Copy to the communication buffer */
      for(i=cgindex[icg]; i<cgindex[icg+1]; i++)
	copy_rvec(vec[i],buf[buf_pos[c]++]);
    }
  }

  return local_pos;
}

static int compact_and_copy_cg(int ncg,int *cell,int local,
			       int nnbc,int *ind,int *index_gl,
			       int *cgindex,int *gatindex,gmx_ga2la_t *ga2la,
			       int *buf,rvec *cg_cm,int *solvent_type)
{
  int c,cg,nat,a0,a1,a,a_gl;
  int local_pos,buf_pos[DD_MAXNBCELL];
  
  for(c=0; c<nnbc; c++)
    buf_pos[c] = ind[c];
  
  local_pos = 0;
  nat = 0;
  for(cg=0; cg<ncg; cg++) {
    c = cell[cg];
    a0 = cgindex[cg];
    a1 = cgindex[cg+1];
    if (c == local) {
      /* Locally compact the arrays.
       * Anything that can be done here avoids access to global arrays.
       */
      cgindex[local_pos] = nat;
      for(a=a0; a<a1; a++) {
	a_gl = gatindex[a];
	gatindex[nat] = a_gl;
	/* The cell number stays 0, so we don't need to set it */
	ga2la[a_gl].a = nat;
	nat++;
      }
      index_gl[local_pos] = index_gl[cg];
      copy_rvec(cg_cm[cg],cg_cm[local_pos]);
      if (solvent_type)
	solvent_type[local_pos] = solvent_type[cg];
      local_pos++;
    } else {
      /* Copy to the communication buffer */
      buf[buf_pos[c]++] = index_gl[cg];
      /* Clear the global indices */
      for(a=a0; a<a1; a++) {
	a_gl = gatindex[a];
	ga2la[a_gl].cell = -1;
      }
    }
  }
  cgindex[local_pos] = nat;

  return local_pos;
}

static int dd_redistribute_cg(FILE *fplog,
			      gmx_domdec_t *dd,t_block *gcgs,
			      t_state *state,rvec cg_cm[],
			      t_forcerec *fr,t_mdatoms *md,
			      t_nrnb *nrnb)
{
#define dd_nbindex(ms0,ns,i) ((((i)[ZZ] - ms0[ZZ])*ns[YY] + (i)[YY] - ms0[YY])*ns[XX] + (i)[XX] - ms0[XX])

  int  *cell,*buf_cg;
  int  ncg[DD_MAXNBCELL],nat[DD_MAXNBCELL];
  int  ncg_r[DD_MAXNBCELL],nat_r[DD_MAXNBCELL],nat_r_max;
  int  nnbc,nvec,c,i,cg,k,k0,k1,d,x,y,z,nreq,ncg_new;
  int  local=-1,dest[DD_MAXNBCELL],src[DD_MAXNBCELL];
  int  sbuf[DD_MAXNBCELL*2],rbuf[DD_MAXNBCELL*2];
  int  buf_cg_ind[DD_MAXNBCELL+1],buf_vs_ind[DD_MAXNBCELL+1];
  int  local_pos,buf_pos[DD_MAXNBCELL],buf_pos_r;
  rvec inv_box;
  ivec ms0,ms1,ns,ind,dev,vlocal,vdest,vsrc;
  int  cg_gl,nrcg;
  real inv_ncg;
  atom_id *cgindex;
#ifdef GMX_MPI
  MPI_Request mpi_req[DD_MAXNBCELL*2];
#endif

  /* Set the required grid shift ranges */
  /* This could be done in setup_dd_grid */
  for(d=0; d<DIM; d++) {
    if (dd->nc[d] == 1) {
      ms0[d] = 0;
      ms1[d] = 0;
    } else if (dd->nc[d] == 2) {
      ms0[d] = 0;
      ms1[d] = 1;
    } else {
      ms0[d] = -1;
      ms1[d] = 1;
    }
    ns[d] = 1 + ms1[d] - ms0[d];
  }

  /* Make a local (possibly) 3x3x3 communication grid */
  /* This could be done in setup_dd_grid */
  /* This should be corrected for triclinic pbc !!! */
  gmx_ddindex2xyz(dd->nc,dd->nodeid,vlocal);
  nnbc = 0;
  for(z=ms0[ZZ]; z<=ms1[ZZ]; z++) {
    vdest[ZZ] = (vlocal[ZZ] + z + dd->nc[ZZ]) % dd->nc[ZZ];
    vsrc[ZZ]  = (vlocal[ZZ] - z + dd->nc[ZZ]) % dd->nc[ZZ];
    for(y=ms0[YY]; y<=ms1[YY]; y++) {
      vdest[YY] = (vlocal[YY] + y + dd->nc[YY]) % dd->nc[YY];
      vsrc[YY]  = (vlocal[YY] - y + dd->nc[YY]) % dd->nc[YY];
      for(x=ms0[XX]; x<=ms1[XX]; x++) {
	vdest[XX] = (vlocal[XX] + x + dd->nc[XX]) % dd->nc[XX];
	vsrc[XX]  = (vlocal[XX] - x + dd->nc[XX]) % dd->nc[XX];
	dest[nnbc] = dd_index(dd->nc,vdest);
	src[nnbc]  = dd_index(dd->nc,vsrc);
	if (x==0 && y==0 && z==0)
	  local = nnbc;
	nnbc++;
      }
    }
  }

  nvec = 1;
  if (state->v)
    nvec++;
  if (state->sd_X)
    nvec++;

  if (dd->ncg_tot > dd->nalloc_i1) {
    dd->nalloc_i1 = over_alloc(dd->ncg_tot);
    srenew(dd->buf_i1,dd->nalloc_i1);
  }
  cell = dd->buf_i1;

  /* Clear the count */
  for(c=0; c<DD_MAXNBCELL; c++) {
    ncg[c] = 0;
    nat[c] = 0;
  }
  
  for(d=0; d<DIM; d++)
    inv_box[d] = divide(dd->nc[d],state->box[d][d]);

  cgindex = dd->cgindex;

  /* Compute the center of geometry for all home charge groups
   * and put them in the box.
   */
  for(cg=0; cg<dd->ncg_home; cg++) {
    k0      = cgindex[cg];
    k1      = cgindex[cg+1];
    nrcg    = k1 - k0;
    if (nrcg == 1) {
      copy_rvec(state->x[k0],cg_cm[cg]);
    }
    else {
      inv_ncg = 1.0/nrcg;
      
      clear_rvec(cg_cm[cg]);
      for(k=k0; (k<k1); k++)
	rvec_inc(cg_cm[cg],state->x[k]);
      for(d=0; (d<DIM); d++)
	cg_cm[cg][d] = inv_ncg*cg_cm[cg][d];
    }
    for(d=DIM-1; d>=0; d--) {
      /* Determine the cell shift.
       * This does not involve pbc shifts as charge groups
       * should not have been put in the box between repartitioning.
       */
      /* We have to do this in two steps to obtain the correct rounding */
      ind[d] = cg_cm[cg][d]*inv_box[d] + 2;
      ind[d] -= 2;
      dev[d] = ind[d] - vlocal[d];
      if (dev[d] < -1 || dev[d] > 1)
	gmx_fatal(FARGS,"The charge group starting at atom %d moved more than one cell: %d %d %d, coords %f %f %f",
		  gcgs->index[dd->index_gl[cg]]+1,
		  dev[XX],dev[YY],dev[ZZ],
		  cg_cm[cg][XX],cg_cm[cg][YY],cg_cm[cg][ZZ]);
      if (ind[d] < 0) {
	ind[d] += dd->nc[d];
	for(k=k0; (k<k1); k++)  {
	  rvec_inc(state->x[k],state->box[d]);
	}
	rvec_inc(cg_cm[cg],state->box[d]);
      } else if (ind[d] >= dd->nc[d]) {
	ind[d] -= dd->nc[d];
	for(k=k0; (k<k1); k++)  {
	  rvec_dec(state->x[k],state->box[d]);
	}
	rvec_dec(cg_cm[cg],state->box[d]);
      }
      if (ns[d] == 1) {
	/* There is only one cell in this dimension, the cg can stay here */
	dev[d] = 0;
      } else if (ns[d] == 2 && dev[d] == -1) {
	/* Two cells in this dimension: shift -1 = shift +1 */
	/* Needs to be corrected for triclinic unit-cells !!! */
	dev[d] = 1;
      }
    }
    i = dd_nbindex(ms0,ns,dev);
    /*
    if (dest[i] != dd_index(dd->nc,ind))
      gmx_fatal(FARGS,"Charge group %d (local index %d) dest[%d] (%d) != ind (%d),shift %d %d %d, coords %f %f %f",
		cc->index_gl[cg],cg,
		i,dest[c],dd_index(dd->nc,ind),
		dev[XX],dev[YY],dev[ZZ],
		cg_cm[cg][XX],cg_cm[cg][YY],cg_cm[cg][ZZ]);
    */
    cell[cg] = i;
    ncg[i] ++;
    nat[i] += nrcg;
  }

  inc_nrnb(nrnb,eNR_RESETX,dd->nat_home);
  
  if (debug) {
    fprintf(debug,"Sending:");
    for(c=0; c<nnbc; c++)
      fprintf(debug," %d %d (%d)",dest[c],ncg[c],nat[c]);
    fprintf(debug,"\n");
  }

  buf_cg_ind[0] = 0;
  buf_vs_ind[0] = 0;
  nreq = 0;
  for(c=0; c<nnbc; c++) {
    buf_cg_ind[c+1] = buf_cg_ind[c];
    buf_vs_ind[c+1] = buf_vs_ind[c];
    if (c != local) {
      sbuf[c*2]   = ncg[c];
      sbuf[c*2+1] = nat[c];
#ifdef GMX_MPI
      MPI_Isend(sbuf+c*2,2*sizeof(int),MPI_BYTE,
		DDRANK(dd,dest[c]),c,dd->all,&mpi_req[nreq++]);
      MPI_Irecv(rbuf+c*2,2*sizeof(int),MPI_BYTE,
		DDRANK(dd,src[c]),c,dd->all,&mpi_req[nreq++]);
#endif
      /* Make index into communication buffers,
       * but without the charge groups that stay on this node.
       */
      buf_cg_ind[c+1] += ncg[c];
      buf_vs_ind[c+1] += nvec*nat[c];
    }
    buf_pos[c] = buf_vs_ind[c];
  }

  if (buf_cg_ind[nnbc] > dd->nalloc_i2) {
    dd->nalloc_i2 = over_alloc(buf_cg_ind[nnbc]);
    srenew(dd->buf_i2,dd->nalloc_i2);
  }
  buf_cg = dd->buf_i2;

  if (buf_vs_ind[nnbc] > dd->nalloc_vs) {
    dd->nalloc_vs = over_alloc(buf_vs_ind[nnbc]);
    srenew(dd->buf_vs,dd->nalloc_vs);
  }

  local_pos = compact_and_copy_vec(dd->ncg_home,cell,local,cgindex,
				   nnbc,buf_pos,state->x,dd->buf_vs);
  if (state->v)
    compact_and_copy_vec(dd->ncg_home,cell,local,cgindex,
			 nnbc,buf_pos,state->v,dd->buf_vs);
  if (state->sd_X) {
    compact_and_copy_vec(dd->ncg_home,cell,local,cgindex,
			 nnbc,buf_pos,state->sd_X,dd->buf_vs);
  }

#ifdef GMX_MPI
  MPI_Waitall(nreq,mpi_req,MPI_STATUSES_IGNORE);
#endif

 
  ncg_new = 0;
  nat_r_max = 0;
  for(c=0; c<nnbc; c++){
    if (c == local) {
      ncg_r[c] = ncg[c];
      nat_r[c] = nat[c];
    } else {
      ncg_r[c] = rbuf[c*2];
      nat_r[c] = rbuf[c*2+1];
      if (nat_r[c] > nat_r_max)
	nat_r_max = nat_r[c];
    }
    ncg_new += ncg_r[c];
  }

  if (debug) {
    fprintf(debug,"Receiving:");
    for(c=0; c<nnbc; c++)
      fprintf(debug," %d %d (%d)",src[c],ncg_r[c],nat_r[c]);
    fprintf(debug,"\n");
  }
 
  if (nvec*nat_r_max > dd->nalloc_vr) {
    dd->nalloc_vr = over_alloc(nvec*nat_r_max);
    srenew(dd->buf_vr,dd->nalloc_vr);
  }

  nreq = 0;
  for(c=0; c<nnbc; c++) {
    if (c != local && nat[c] > 0) {
#ifdef GMX_MPI
      MPI_Isend(dd->buf_vs+buf_vs_ind[c],nvec*nat[c]*sizeof(rvec),MPI_BYTE,
		DDRANK(dd,dest[c]),c,dd->all,&mpi_req[nreq++]);
#endif
    }
  }
  
  for(c=0; c<nnbc; c++) {
    if (c != local && nat_r[c] > 0) {
#ifdef GMX_MPI
      MPI_Recv(dd->buf_vr,nvec*nat_r[c]*sizeof(rvec),MPI_BYTE,
	       DDRANK(dd,src[c]),c,dd->all,MPI_STATUS_IGNORE);
#endif
      memcpy(state->x+local_pos,dd->buf_vr,nat_r[c]*sizeof(rvec));
      buf_pos_r = nat_r[c];
      if (state->v) {
	memcpy(state->v+local_pos,dd->buf_vr+buf_pos_r,nat_r[c]*sizeof(rvec));
	buf_pos_r += nat_r[c];
      }
      if (state->sd_X) {
	memcpy(state->sd_X+local_pos,dd->buf_vr+buf_pos_r,
		 nat_r[c]*sizeof(rvec));
	buf_pos_r += nat_r[c];
      }
      local_pos += nat_r[c];
    }
  }

  if (ncg_new > dd->cg_nalloc || dd->cg_nalloc == 0) {
    dd->cg_nalloc = over_alloc(ncg_new);
    srenew(dd->index_gl,dd->cg_nalloc);
    srenew(dd->cgindex,dd->cg_nalloc+1);
  }

  /* Clear the local indices, except for the home cell.
   * The home cell indices are updated in compact_and_copy_cg.
   */
  clear_dd_indices(dd,dd->nat_home);

  local_pos =
    compact_and_copy_cg(dd->ncg_home,cell,local,
			nnbc,buf_cg_ind,dd->index_gl,
			dd->cgindex,dd->gatindex,dd->ga2la,
			buf_cg,cg_cm,
			fr->solvent_opt==esolNO ? NULL : fr->solvent_type);
  
#ifdef GMX_MPI
  MPI_Waitall(nreq,mpi_req,MPI_STATUSES_IGNORE);
#endif
  
  nreq = 0;
  for(c=0; c<nnbc; c++) {
    if (c != local) {
#ifdef GMX_MPI
      if (ncg[c])
	MPI_Isend(buf_cg+buf_cg_ind[c],ncg[c]*sizeof(int),MPI_BYTE,
		  DDRANK(dd,dest[c]),c,dd->all,&mpi_req[nreq++]);
      if (ncg_r[c])
	/* Receive the cg's in place */
	MPI_Irecv(dd->index_gl+local_pos,ncg_r[c]*sizeof(int),MPI_BYTE,
		  DDRANK(dd,src[c]),c,dd->all,&mpi_req[nreq++]);
#endif
      local_pos += ncg_r[c];
    }
  }
  
  dd->ncg_home = ncg[local];
  dd->nat_home = nat[local];
  for(c=0; c<nnbc; c++){
    if (c != local) {
      dd->ncg_home += ncg_r[c];
      dd->nat_home += nat_r[c];
    }
  }
  
#ifdef GMX_MPI
  MPI_Waitall(nreq,mpi_req,MPI_STATUSES_IGNORE);
#endif
  
  dd->bMasterHasAllCG = FALSE;

  /* Compute the center of mass for the new home charge groups.
   * We could also communicate this, but as we will only have
   * a few new ones after nstlist steps it will be faster to recompute.
   */
  k0 = dd->cgindex[ncg[local]];
  for(cg=ncg[local]; cg<dd->ncg_home; cg++) {
    cg_gl = dd->index_gl[cg];
    nrcg = gcgs->index[cg_gl+1] - gcgs->index[cg_gl];
    k1 = k0 + nrcg;
    clear_rvec(cg_cm[cg]);
    for(k=k0; (k<k1); k++) {
      rvec_inc(cg_cm[cg],state->x[k]);
    }
    if (nrcg > 1) {
      inv_ncg = 1.0/nrcg;
      for(d=0; (d<DIM); d++)
	cg_cm[cg][d] = inv_ncg*cg_cm[cg][d];
    }
    /* Set the charge group size */
    dd->cgindex[cg+1] = k1;
    k0 = k1;
  }
  
  inc_nrnb(nrnb,eNR_CGCM,dd->nat_home-nat[local]);

  if (debug)
    fprintf(debug,"Finished repartitioning\n");

  return ncg[local];
}

void setup_dd_grid(FILE *fplog,matrix box,gmx_domdec_t *dd)
{
  int  d,i,j,m;
  ivec xyz,tmp,h,s;
  int  ncell,ncellp;
  ivec *dd_c,dd_cp[DD_MAXICELL];
  gmx_domdec_ns_t *icell;

  gmx_ddindex2xyz(dd->nc,dd->nodeid,xyz);

  /* The partition order is z, y, x */
  dd->ndim = 0;
  for(d=DIM-1; d>=0; d--) {
    if (dd->nc[d] > 1) {
      dd->dim[dd->ndim] = d;
      copy_ivec(xyz,tmp);
      tmp[d] = (tmp[d] + 1) % dd->nc[d];
      dd->neighbor[dd->ndim][0] = dd_index(dd->nc,tmp);
      copy_ivec(xyz,tmp);
      tmp[d] = (tmp[d] - 1 + dd->nc[d]) % dd->nc[d];
      dd->neighbor[dd->ndim][1] = dd_index(dd->nc,tmp);
      dd->ndim++;
    }
  }
  
  dd_c = dd->shift;
  
  fprintf(fplog,"Making %dD domain decomposition %d x %d x %d\n",
	  dd->ndim,dd->nc[XX],dd->nc[YY],dd->nc[ZZ]);
  if (DDMASTER(dd))
    fprintf(stderr,"Making %dD domain decomposition %d x %d x %d\n",
	    dd->ndim,dd->nc[XX],dd->nc[YY],dd->nc[ZZ]);
  switch (dd->ndim) {
  case 3:
    ncell = dd_c3n;
    for(i=0; i<ncell; i++)
      copy_ivec(dd_c3[i],dd_c[i]);
    ncellp = dd_cp3n;
    for(i=0; i<ncellp; i++)
      copy_ivec(dd_cp3[i],dd_cp[i]);
    break;
  case 2:
    ncell = dd_c2n;
    for(i=0; i<ncell; i++) {
      m = 0;
      for(d=0; d<DIM; d++) {
	if (dd->nc[d] > 1)
	  dd_c[i][d] = dd_c2[i][m++];
	else 
	  dd_c[i][d] = 0;
      }
    }
    ncellp = dd_cp2n;
    for(i=0; i<ncellp; i++)
      copy_ivec(dd_cp2[i],dd_cp[i]);
    break;
  case 1:
    ncell = dd_c1n;
    for(i=0; i<ncell; i++) {
      m = 0;
      for(d=0; d<DIM; d++) {
	if (dd->nc[d] > 1)
	  dd_c[i][d] = dd_c1[i][m++];
	else 
	  dd_c[i][d] = 0;
      }
    }
    ncellp = dd_cp1n;
    for(i=0; i<ncellp; i++)
      copy_ivec(dd_cp1[i],dd_cp[i]);
    break;
  default:
    gmx_fatal(FARGS,"Can only do 1, 2 or 3D domain decomposition");
    ncell = 0;
    ncellp = 0;
  }

  dd->ncell  = ncell;
  gmx_ddindex2xyz(dd->nc,dd->nodeid,h);
  for(i=0; i<ncell; i++) {
    for(d=0; d<DIM; d++) {
      s[d] = h[d] + dd_c[i][d];
      /* Currently only works for rectangular boxes.
       * For even numbers of grid cells it is easy to extend
       * to things like rhombic dodecahedrons.
       */
      if (s[d] < 0)
	s[d] += dd->nc[d];
      else if (s[d] >= dd->nc[d])
	s[d] -= dd->nc[d];
    }
    for(d=0; d<DIM; d++) {
      s[d] = h[d] - dd_c[i][d];
      /* Currently only works for rectangular boxes.
       * For even numbers of grid cells it is easy to extend
       * to things like rhombic dodecahedrons.
       */
      if (s[d] < 0)
	s[d] += dd->nc[d];
      else if (s[d] >= dd->nc[d])
	s[d] -= dd->nc[d];
    }
  }
  dd->nicell = ncellp;
  for(i=0; i<dd->nicell; i++) {
    if (dd_cp[i][0] != i)
      gmx_fatal(FARGS,"Internal inconsistency in the dd grid setup");
    icell = &dd->icell[i];
    icell->j0 = dd_cp[i][1];
    icell->j1 = dd_cp[i][2];
    for(d=0; d<DIM; d++) {
      if (dd->nc[d] == 1) {
	icell->shift0[d] = -1;
	icell->shift1[d] = 1;
      } else {
	/*
	icell->shift0[d] = 0;
	icell->shift1[d] = 0;
	for(j=icell->j0; j<icell->j1; j++) {
	  if (dd->shift[j][d] > dd->shift[i][d])
	    icell->shift0[d] = -1;
	  if (dd->shift[j][d] < dd->shift[i][d])
	    icell->shift1[d] = 1;
	}
	*/
	
	int shift_diff;
	
	/* Assume the shift are not more than 1 cell */
	icell->shift0[d] = 1;
	icell->shift1[d] = -1;
	for(j=icell->j0; j<icell->j1; j++) {
	  shift_diff = dd->shift[j][d] - dd->shift[i][d];
	  if (shift_diff < icell->shift0[d])
	    icell->shift0[d] = shift_diff;
	  if (shift_diff > icell->shift1[d])
	    icell->shift1[d] = shift_diff;
	}
      }
    }
  }
}

static int count_intercg_excls(t_block *cgs,t_block *excls)
{
  int n,cg,at0,at1,at,excl,atj;
  
  n = 0;
  for(cg=0; cg<cgs->nr; cg++) {
    at0 = cgs->index[cg];
    at1 = cgs->index[cg+1];
    for(at=at0; at<at1; at++) {
      for(excl=excls->index[at]; excl<excls->index[at+1]; excl++) {
	atj = excls->a[excl];
	if (atj > at && (atj < at0 || atj >= at1))
	  n++;
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

static gmx_reverse_top_t make_reverse_top(int natoms,t_idef *idef,int *nint)
{
  int *count,i;
  gmx_reverse_top_t rt;

  /* Count the interactions */
  snew(count,natoms);
  low_make_reverse_top(idef,count,&rt,FALSE);

  snew(rt.index,natoms+1);
  rt.index[0] = 0;
  for(i=0; i<natoms; i++) {
    rt.index[i+1] = rt.index[i] + count[i];
    count[i] = 0;
  }
  snew(rt.il,rt.index[natoms]);

  /* Store the interactions */
  *nint = low_make_reverse_top(idef,count,&rt,TRUE);

  sfree(count);

  return rt;
}

void dd_make_reverse_top(FILE *fplog,
			 gmx_domdec_t *dd,t_topology *top,
			 bool bDynamics,int eeltype)
{
  int natoms,a;
  
  natoms = top->atoms.nr;

  dd->reverse_top = make_reverse_top(natoms,&top->idef,&dd->nbonded_global);
  
  dd->n_intercg_excl = count_intercg_excls(&top->blocks[ebCGS],
					   &top->blocks[ebEXCLS]);
  if (EEL_FULL(eeltype) && dd->n_intercg_excl)
    fprintf(fplog,"There are %d inter charge-group exclusions,\n"
	    "will use an extra communication step for exclusion forces for %s\n",
	    dd->n_intercg_excl,eel_names[eeltype]);
 
  snew(dd->ga2la,natoms*sizeof(gmx_ga2la_t));
  for(a=0; a<natoms; a++)
    dd->ga2la[a].cell = -1;
  
  if (top->idef.il[F_CONSTR].nr > 0)
    dd->constraints = init_domdec_constraints(natoms,&top->idef,bDynamics);
}

static int *dd_pmenodes(t_commrec *cr)
{
  int *pmenodes;
  int n,i,p0,p1;

  snew(pmenodes,cr->npmenodes);
  n = 0;
  for(i=0; i<cr->dd->nnodes; i++) {
    p0 = ( i   *cr->npmenodes + cr->npmenodes/2)/cr->dd->nnodes;
    p1 = ((i+1)*cr->npmenodes + cr->npmenodes/2)/cr->dd->nnodes;
    if (i+1 == cr->dd->nnodes || p1 > p0) {
      if (debug)
	fprintf(debug,"pmenode[%d] = %d\n",n,i+1+n);
      pmenodes[n] = i + 1 + n;
      n++;
    }
  }

  return pmenodes;
}

static void dd_coords2pmecoords(gmx_domdec_t *dd,ivec coords,ivec coords_pme)
{
  int nc,ntot;

  nc   = dd->nc[dd->pmedim];
  ntot = dd->ntot[dd->pmedim];
  copy_ivec(coords,coords_pme);
  coords_pme[dd->pmedim] =
    nc + (coords[dd->pmedim]*(ntot - nc) + (ntot - nc)/2)/nc;
}

int gmx_ddindex2pmeslab(t_commrec *cr,int ddindex)
{
  gmx_domdec_t *dd;
  ivec coords,coords_pme,nc;
  int  slab;

  dd = cr->dd;
  if (dd->bCartesian) {
    gmx_ddindex2xyz(dd->nc,ddindex,coords);
    dd_coords2pmecoords(dd,coords,coords_pme);
    copy_ivec(dd->ntot,nc);
    nc[dd->pmedim]         -= dd->nc[dd->pmedim];
    coords_pme[dd->pmedim] -= dd->nc[dd->pmedim];

    slab = (coords_pme[XX]*nc[YY] + coords_pme[YY])*nc[ZZ] + coords_pme[ZZ];
  } else {
    slab = (ddindex*cr->npmenodes + cr->npmenodes/2)/dd->nnodes;
  }

  return slab;
}

int gmx_ddindex2nodeid(t_commrec *cr,int ddindex)
{
  ivec coords;
  int  nodeid;

  if (cr->dd->bCartesian) {
    gmx_ddindex2xyz(cr->dd->nc,ddindex,coords);
#ifdef GMX_MPI
    MPI_Cart_rank(cr->mpi_comm_mysim,coords,&nodeid);
#endif
  } else {
    nodeid = ddindex;
    if (cr->dd->pmenodes)
      nodeid += gmx_ddindex2pmeslab(cr,ddindex);
  }
  
  return nodeid;
}

static int dd_node2pmenode(t_commrec *cr,int nodeid)
{
  gmx_domdec_t *dd;
  ivec coords,coords_pme;
  int  i;
  int  pmenode;
  
  dd = cr->dd;

  pmenode = -1;
  /* This assumes a uniform x domain decomposition grid cell size */
  if (dd->bCartesian) {
#ifdef GMX_MPI
    MPI_Cart_coords(cr->mpi_comm_mysim,nodeid,DIM,coords);
    if (coords[dd->pmedim] < dd->nc[dd->pmedim]) {
      /* This is a PP node */
      dd_coords2pmecoords(cr->dd,coords,coords_pme);
      MPI_Cart_rank(cr->mpi_comm_mysim,coords_pme,&pmenode);
    }
#endif
  } else {
    /* This assumes DD cells with identical x coordinates
     * are numbered sequentially.
     */
    if (dd->pmenodes == NULL) {
      if (nodeid < dd->nnodes) {
	pmenode = dd->nnodes +
	  (nodeid*cr->npmenodes + cr->npmenodes/2)/dd->nnodes;
      }
    } else {
      i = 0;
      while (nodeid > dd->pmenodes[i])
	i++;
      if (nodeid < dd->pmenodes[i])
	pmenode = dd->pmenodes[i];
    }
  }

  return pmenode;
}

bool gmx_pmeonlynode(t_commrec *cr,int nodeid)
{
  bool bPMEOnlyNode;

  if (DOMAINDECOMP(cr)) {
    bPMEOnlyNode = (dd_node2pmenode(cr,nodeid) == -1);
  } else {
    bPMEOnlyNode = FALSE;
  }

  return bPMEOnlyNode;
}

static bool receive_vir_ener(t_commrec *cr)
{
  int  pmenode,coords[DIM],rank;
  bool bReceive;

  pmenode = dd_node2pmenode(cr,cr->nodeid);

  bReceive = TRUE;
  if (cr->npmenodes < cr->dd->nnodes) {
    if (cr->dd->bCartesian) {
#ifdef GMX_MPI
      MPI_Cart_coords(cr->mpi_comm_mysim,cr->nodeid,DIM,coords);
      coords[cr->dd->pmedim]++;
      if (coords[cr->dd->pmedim] < cr->dd->nc[cr->dd->pmedim]) {
	MPI_Cart_rank(cr->mpi_comm_mysim,coords,&rank);
	if (dd_node2pmenode(cr,rank) == pmenode) {
	  /* This is not the last PP node for pmenode */
	  bReceive = FALSE;
	}
      }
#endif  
    } else {
      if (cr->nodeid+1 < cr->nnodes &&
	  dd_node2pmenode(cr,cr->nodeid+1) == pmenode) {
	/* This is not the last PP node for pmenode */
	bReceive = FALSE;
      }
    }
  }
  
  return bReceive;
}

void make_dd_communicators(FILE *fplog,t_commrec *cr,bool bCartesian)
{
  gmx_domdec_t *dd;
  bool bDiv[DIM];
  int  i;
  ivec periods,coords;
#ifdef GMX_MPI
  MPI_Comm comm_cart;
#endif

  dd = cr->dd;

  copy_ivec(dd->nc,dd->ntot);
  
  dd->bCartesian = bCartesian;
  if (dd->bCartesian) {
    if (cr->npmenodes > 0) {
      for(i=1; i<DIM; i++)
	bDiv[i] = ((cr->npmenodes*dd->nc[i]) % (dd->nnodes) == 0);
      if (bDiv[YY] || bDiv[ZZ]) {
	/* We choose the direction that provides the thinnest slab
	 * of PME only nodes as this will have the least effect
	 * on the PP communication.
	 * But for the PME communication the opposite might be better.
	 */
	if (bDiv[YY] && (!bDiv[ZZ] || dd->nc[YY] <= dd->nc[ZZ])) {
	  dd->pmedim = YY;
	} else {
	  dd->pmedim = ZZ;
	}
	dd->ntot[dd->pmedim] += (cr->npmenodes*dd->nc[dd->pmedim])/dd->nnodes;
      } else {
	dd->bCartesian = FALSE;
      }
    }
  }

#ifdef GMX_MPI
  if (dd->bCartesian) {
    fprintf(fplog,"Will use a Cartesian communicator: %d x %d x %d\n",
	    dd->ntot[XX],dd->ntot[YY],dd->ntot[ZZ]);

    for(i=0; i<DIM; i++) {
      periods[i] = TRUE;
    }

    MPI_Cart_create(cr->mpi_comm_mysim,DIM,dd->ntot,periods,TRUE,&comm_cart);
    
    /* With this assigment we loose the link to the original communicator
     * which will usually be MPI_COMM_WORLD, unless have multisim.
     */
    cr->mpi_comm_mysim = comm_cart;
    
    MPI_Comm_rank(cr->mpi_comm_mysim,&cr->nodeid);

    MPI_Cart_coords(cr->mpi_comm_mysim,cr->nodeid,DIM,coords);

    fprintf(fplog,"Cartesian nodeid %d, coordinates %d %d %d\n\n",
	    cr->nodeid,coords[0],coords[1],coords[2]);
    
    if (coords[dd->pmedim] < dd->nc[dd->pmedim])
      cr->duty |= DUTY_PP;
    if (cr->npmenodes == 0 || coords[dd->pmedim] >= dd->nc[dd->pmedim])
      cr->duty |= DUTY_PME;

    /* We always split here, even when we do not have separate PME nodes,
     * as MPI_Comm_split is the only MPI call that allows us to reorder
     * the nodeids, so we do not have to keep an extra index ourselves.
     */
    MPI_Comm_split(cr->mpi_comm_mysim,
		   cr->duty,
		   dd_index(dd->ntot,coords),
		   &cr->mpi_comm_mygroup);
  } else {
    fprintf(fplog,"Will not use a Cartesian communicator\n\n");
    
    if (cr->npmenodes == 0) {
      cr->duty |= (DUTY_PP | DUTY_PME);

      cr->mpi_comm_mygroup = cr->mpi_comm_mysim;
    } else {
      if (getenv("GMX_ORDER_PP_PME") == NULL) {
	/* Interleave the PP-only and PME-only nodes,
	 * as on clusters with dual-core machines this will double
	 * the communication bandwidth of the PME processes
	 * and thus speed up the PP <-> PME and inter PME communication.
	 */
	cr->dd->pmenodes = dd_pmenodes(cr);
      }

      if (dd_node2pmenode(cr,cr->nodeid) == -1)
	cr->duty |= DUTY_PME;
      else
	cr->duty |= DUTY_PP;

      MPI_Comm_split(cr->mpi_comm_mysim,
		     cr->duty,
		     cr->nodeid,
		     &cr->mpi_comm_mygroup);
    }
  }
  
  dd->all = cr->mpi_comm_mygroup;
  
  MPI_Comm_rank(dd->all,&dd->nodeid);
  
  fprintf(fplog,"Domain decomposition nodeid %d\n\n",dd->nodeid);
  
  if (cr->npmenodes > 0) {
    if (cr->duty & DUTY_PP) {
      /* Make the ring smaller */
      cr->left  = (dd->nodeid - 1 + dd->nnodes) % dd->nnodes;
      cr->right = (dd->nodeid + 1) % dd->nnodes;
    }

    fprintf(fplog,"This is a %s only node\n\n",
	    (cr->duty & DUTY_PP) ? "particle-particle" : "PME-mesh");
  }
#endif
  
  if (!(cr->duty & DUTY_PME)) {
    dd->pme_nodeid = dd_node2pmenode(cr,cr->nodeid);
    dd->pme_receive_vir_ener = receive_vir_ener(cr);
  }
}

gmx_domdec_t *init_domain_decomposition(FILE *fplog,t_commrec *cr,ivec nc)
{
  gmx_domdec_t *dd;
  
  fprintf(fplog,
	  "Domain decomposition grid %d x %d x %d, separate PME nodes %d\n",
	  nc[XX],nc[YY],nc[ZZ],cr->npmenodes);

  snew(dd,1);

  copy_ivec(nc,dd->nc);
  dd->nnodes = dd->nc[XX]*dd->nc[YY]*dd->nc[ZZ];
  if (dd->nnodes != cr->nnodes - cr->npmenodes)
    gmx_fatal(FARGS,"The size of the domain decomposition grid (%d) does not match the number of nodes (%d)\n",dd->nnodes,cr->nnodes - cr->npmenodes);

  return dd;
}

static void setup_dd_communication(FILE *fplog,gmx_domdec_t *dd,t_block *gcgs,
				   rvec x[],rvec buf[],matrix box,rvec cg_cm[],
				   real r_comm)
{
  int dim_ind,dim,nat_tot,ncell,cell,celli,c,i,cg,cg_gl,nrcg,d;
  int *ncg_cell,*index_gl,*cgindex;
  gmx_domdec_comm_t *cc;
  ivec xyz;
  rvec corner;
  real corner0=0,corner1=0,r_comm2,r2,inv_ncg;
  int  nsend,nat;

  if (debug)
    fprintf(debug,"Setting up DD communication\n");

  gmx_ddindex2xyz(dd->nc,dd->nodeid,xyz);
  /* This should be changed for non-uniform grids */
  for(d=0; d<DIM; d++)
    corner[d] = xyz[d]*box[d][d]/dd->nc[d];
  if (dd->ndim >= 2) {
    /* Set the upper-right corner for rounding */
    d = dd->dim[0];
    corner0 = ((xyz[d] + 1) % dd->nc[d])*box[d][d]/dd->nc[d];
    d = dd->dim[1];
    corner1 = ((xyz[d] + 1) % dd->nc[d])*box[d][d]/dd->nc[d];
  }

  r_comm2 = sqr(r_comm);

  for(dim=0; dim<dd->ndim; dim++) {
    d = dd->dim[dim];
    /* Should be fixed for triclinic boxes */
    if (box[d][d]/dd->nc[d] < r_comm)
      gmx_fatal(FARGS,"One of the domain decomposition grid cell sizes (%f %f %f) is smaller than the cut-off (%f)",
		box[XX][XX]/dd->nc[XX],
		box[YY][YY]/dd->nc[YY],
		box[ZZ][ZZ]/dd->nc[ZZ],r_comm);
  }

  ncg_cell = dd->ncg_cell;
  index_gl = dd->index_gl;
  cgindex  = dd->cgindex;
  
  ncg_cell[0] = 0;
  ncg_cell[1] = dd->ncg_home;

  nat_tot = dd->nat_home;
  ncell = 1;
  for(dim_ind=0; dim_ind<dd->ndim; dim_ind++) {
    dim = dd->dim[dim_ind];
    cc = &dd->comm[dim_ind];
    nsend = 0;
    nat = 0;
    for(cell=0; cell<ncell; cell++) {
      cc->nsend[cell] = 0;
      celli = cell_perm[dim_ind][cell];
      for(cg=ncg_cell[celli]; cg<ncg_cell[celli+1]; cg++) {
	r2 = sqr(cg_cm[cg][dim] - corner[dim]);
	/* Rounding gives at most a 16% reduction in communicated atoms */
	if (dim_ind >= 1 && (celli == 1 || celli == 2))
	  r2 += sqr(cg_cm[cg][dd->dim[0]] - corner0);
	if (dim_ind == 2 && (celli == 2 || celli == 3))
	  r2 += sqr(cg_cm[cg][dd->dim[1]] - corner1);
	if (r2 < r_comm2) {
	  /* Make an index to the local charge groups */
	  if (nsend >= cc->nalloc) {
	    cc->nalloc += CG_ALLOC_SIZE;
	    srenew(cc->index,cc->nalloc);
	  }
	  if (nsend >= dd->nalloc_i1) {
	    dd->nalloc_i1 += CG_ALLOC_SIZE;
	    srenew(dd->buf_i1,dd->nalloc_i1);
	  }
	  cc->index[nsend] = cg;
	  dd->buf_i1[nsend] = index_gl[cg];
	  cc->nsend[cell]++;
	  nsend++;
	  for(i=cgindex[cg]; i<cgindex[cg+1]; i++) {
	    copy_rvec(x[i],buf[nat]);
	    nat++;
	  }
	}
      }
    }
    cc->nsend[ncell]   = nsend;
    cc->nsend[ncell+1] = nat;
    /* Communicate the number of cg's and atoms to receive */
    dd_sendrecv_int(dd, dim_ind, ddBackward,
		    cc->nsend, ncell+2,
		    cc->nrecv, ncell+2);
    /* Communicate the global cg indices, receive in place */
    if (ncg_cell[ncell] + cc->nrecv[ncell] > dd->cg_nalloc
	|| dd->cg_nalloc == 0) {
      dd->cg_nalloc = over_alloc(ncg_cell[ncell] + cc->nrecv[ncell]);
      srenew(index_gl,dd->cg_nalloc);
      srenew(cgindex,dd->cg_nalloc+1);
    }
    dd_sendrecv_int(dd, dim_ind, ddBackward,
		    dd->buf_i1,               nsend,
		    index_gl+ncg_cell[ncell], cc->nrecv[ncell]);
    /* Communicate the coordinate, receive in place */
    dd_sendrecv_rvec(dd, dim_ind, ddBackward,
		     buf,       cc->nsend[ncell+1],
		     x+nat_tot, cc->nrecv[ncell+1]);
    /* Make the charge group index and determine cgcm */
    for(cell=ncell; cell<2*ncell; cell++) {
      ncg_cell[cell+1] = ncg_cell[cell] + cc->nrecv[cell-ncell];
      for(cg=ncg_cell[cell]; cg<ncg_cell[cell+1]; cg++) {
	cg_gl = index_gl[cg];
	nrcg = gcgs->index[cg_gl+1] - gcgs->index[cg_gl];
	clear_rvec(cg_cm[cg]);
	for(i=0; i<nrcg; i++)
	  rvec_inc(cg_cm[cg],x[nat_tot++]);
	if (nrcg > 1) {
	  inv_ncg = 1.0/nrcg;
	  for(d=0; (d<DIM); d++)
	    cg_cm[cg][d] *= inv_ncg;
	}
	cgindex[cg+1] = cgindex[cg] + nrcg;
      }
    }
    ncell += ncell;
  }
  dd->index_gl = index_gl;
  dd->cgindex  = cgindex;

  dd->ncg_tot = ncg_cell[dd->ncell];
  dd->nat_tot = nat_tot;
  dd->nat_tot_con = nat_tot;

  if (debug) {
    fprintf(debug,"Finished setting up DD communication, cells:");
    for(c=0; c<dd->ncell; c++)
      fprintf(debug," %d",dd->ncg_cell[c+1]-dd->ncg_cell[c]);
    fprintf(debug,"\n");
  }
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

static void make_local_bondeds(gmx_domdec_t *dd,t_idef *idef)
{
  int i,gat,j,ftype,nral,d,k,kc;
  int *index,*rtil;
  t_iatom *iatoms,tiatoms[1+MAXATOMLIST],*liatoms;
  bool bUse;
  ivec shift_min;
  gmx_ga2la_t *ga2la;
  t_ilist *il;

  index = dd->reverse_top.index;
  rtil  = dd->reverse_top.il;

  /* Clear the counts */
  for(ftype=0; ftype<F_NRE; ftype++)
    idef->il[ftype].nr = 0;
  dd->nbonded_local = 0;
  
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
	  srenew(il->iatoms,il->nalloc*sizeof(t_iatom));
	}
	liatoms = il->iatoms + il->nr;
	liatoms[0] = iatoms[0];
	for(k=1; k<=nral; k++)
	  liatoms[k] = tiatoms[k];
	/* Sum locally so we can check in global_stat if we have everything */
	dd->nbonded_local++;
	il->nr += 1+nral;
      }
      j += 1 + nral;
    }
  }
}

static void make_local_exclusions(gmx_domdec_t *dd,t_forcerec *fr,
				  t_block *excls,t_block *lexcls)
{
  int nicell,n,ic,jlat0,jlat1,jlat,cg,la0,la1,la,a,i,j;
  gmx_ga2la_t *ga2la;

  /* Since for RF and PME we need to loop over the exclusions
   * we should store each exclusion only once. This is done
   * using the same cell scheme as used for neighbor searching.
   * The exclusions involving non-home atoms are stored only
   * one way: atom j is in the excl list of i only for j > i,
   * where i and j are local atom numbers.
   */

  /* We do not check for missing exclusions. For cut-off simulations
   * this is fine, but it could be a problem with full electrostatics.
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
  for(ic=0; ic<nicell; ic++) {
    jlat0 = dd->cgindex[dd->icell[ic].jcg0];
    jlat1 = dd->cgindex[dd->icell[ic].jcg1];
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
	      jlat = ga2la->a;
	      /* Check to avoid double counts */
	      if (jlat >= jlat0 && jlat < jlat1) 
		lexcls->a[n++] = jlat;
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
	break;
      default:
	gmx_fatal(FARGS,"Unkown solvent type %d",fr->solvent_type[cg]);
      }
    }
  }
  if (dd->n_intercg_excl == 0) {
    /* There are no exclusions involving non-home charge groups */
    la0 = dd->cgindex[dd->icell[nicell].cg1];
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

static void make_local_cgs(gmx_domdec_t *dd,t_block *lcgs)
{
  lcgs->nr    = dd->ncg_tot;
  lcgs->index = dd->cgindex;
  lcgs->nra   = dd->nat_tot;
}

static void set_cg_boundaries(gmx_domdec_t *dd)
{
  int c;

  for(c=0; c<dd->nicell; c++) {
    dd->icell[c].cg1  = dd->ncg_cell[c+1];
    dd->icell[c].jcg0 = dd->ncg_cell[dd->icell[c].j0];
    dd->icell[c].jcg1 = dd->ncg_cell[dd->icell[c].j1];
  }
}

static void dd_update_ns_border(gmx_domdec_t *dd,t_nsborder *nsb)
{
  nsb->nodeid    = 0 /* dd->nodeid */;
  nsb->cgtotal   = dd->ncg_tot;
  nsb->index[nsb->nodeid]  = 0;
  nsb->homenr[nsb->nodeid] = dd->nat_home;
  nsb->cgload[nsb->nodeid] = nsb->cgtotal;
}

static void make_local_top(FILE *fplog,gmx_domdec_t *dd,
			   t_forcerec *fr,
			   t_topology *top,t_topology *ltop,
			   t_nsborder *nsb)
{
  if (debug)
    fprintf(debug,"Making local topology\n");

  ltop->name  = top->name;
  make_local_cgs(dd,&ltop->blocks[ebCGS]);

  make_local_exclusions(dd,fr,&top->blocks[ebEXCLS],&ltop->blocks[ebEXCLS]);

#ifdef ONE_MOLTYPE_ONE_CG
  natoms_global = top->blocks[ebCGS].nra;

  make_local_idef(dd,&top->idef,&ltop->idef);
#else
  make_local_bondeds(dd,&ltop->idef);
#endif
 
  ltop->atoms = top->atoms;
  ltop->atomtypes = top->atomtypes;
  /*
  for(eb=0; eb<ebNR; eb++) {
    if (eb != ebCGS)
      ltop->blocks[eb] = top->blocks[eb];
  }
  */
  ltop->symtab = top->symtab;

  dd_update_ns_border(dd,nsb);
}

static void dump_conf(gmx_domdec_t *dd,
		      char *name,rvec *x,matrix box)
{
  char str[STRLEN];
  FILE *fp;
  int c,i,j,a;
  
  sprintf(str,"%s%d.pdb",name,dd->nodeid);
  fp = fopen(str,"w");
  fprintf(fp,"CRYST1%9.3f%9.3f%9.3f%7.2f%7.2f%7.2f P 1           1\n",
	  10*norm(box[XX]),10*norm(box[YY]),10*norm(box[ZZ]),
	  90.0,90.0,90.0);
  a = 0;
  for(c=0; c<dd->ncell; c++) {
    for(i=dd->ncg_cell[c]; i<dd->ncg_cell[c+1]; i++)
      for(j=dd->cgindex[i]; j<dd->cgindex[i+1]; j++) {
	fprintf(fp,"%-6s%5u  %-4.4s%3.3s %c%4d    %8.3f%8.3f%8.3f%6.2f%6.2f\n",
		"ATOM",dd->gatindex[i]+1,"C","ALA",' ',a+1,
		10*x[a][XX],10*x[a][YY],10*x[a][ZZ],
		1.0,1.0*c/(dd->ncell-1.0));
	a++;
    }
  }
  fclose(fp);
}

void dd_partition_system(FILE         *fplog,
			 gmx_domdec_t *dd,
			 bool         bMasterState,
			 t_state      *state_global,
			 t_topology   *top_global,
			 t_inputrec   *ir,
			 t_state      *state_local,
			 rvec         *buf,
			 t_mdatoms    *mdatoms,
			 t_topology   *top_local,
			 t_nsborder   *nsb,
			 t_forcerec   *fr,
			 t_nrnb       *nrnb)
{
  int cg0;

  if (ir->ePBC == epbcXYZ)
    gmx_fatal(FARGS,"pbc=%s is not supported (yet), use %s",
	      epbc_names[epbcXYZ],epbc_names[epbcFULL]);
  
  if (ir->eConstrAlg == estSHAKE)
    gmx_fatal(FARGS,"%s is not supported (yet), use %s",
	      eshake_names[estSHAKE],eshake_names[estLINCS]);

  if (bMasterState) {
    get_cg_distribution(fplog,dd,
			&top_global->blocks[ebCGS],
			state_global->box,state_global->x);
    
    dd_distribute_state(dd,&top_global->blocks[ebCGS],
			state_global,state_local);
    
    make_local_cgs(dd,&top_local->blocks[ebCGS]);

    calc_cgcm(fplog,0,dd->ncg_home,
	      &top_local->blocks[ebCGS],state_local->x,fr->cg_cm);
    
    inc_nrnb(nrnb,eNR_CGCM,dd->nat_home);

    cg0 = 0;
  } else {
    cg0 = dd_redistribute_cg(fplog,dd,&top_global->blocks[ebCGS],
			     state_local,fr->cg_cm,fr,mdatoms,nrnb);
  }
  
  /* Setup up the communication and communicate the coordinates */
  setup_dd_communication(fplog,dd,&top_global->blocks[ebCGS],
			 state_local->x,buf,state_local->box,
			 fr->cg_cm,fr->rlistlong);
  
  /* Set the indices */
  make_dd_indices(dd,&top_global->blocks[ebCGS],cg0,fr);

  /* Set the charge group boundaries for neighbor searching */
  set_cg_boundaries(dd);

  /* Update the rest of the forcerec */
  fr->cg0 = 0;
  fr->hcg = dd->ncg_tot;
  /* The normal force array should also be dynamic and reallocate here */
  if (fr->bTwinRange) {
    fr->f_twin_n = dd->nat_tot;
    if (fr->f_twin_n > fr->f_twin_nalloc) {
      fr->f_twin_nalloc = over_alloc(fr->f_twin_n);
      srenew(fr->f_twin,fr->f_twin_nalloc);
    }
  }
  if (EEL_FULL(fr->eeltype)) {
    fr->f_el_recip_n = (dd->n_intercg_excl ? dd->nat_tot : dd->nat_home);
    if (fr->f_el_recip_n > fr->f_el_recip_nalloc) {
      fr->f_el_recip_nalloc = over_alloc(fr->f_el_recip_n);
      srenew(fr->f_el_recip,fr->f_el_recip_nalloc);
    }
  }
  
  /* Extract a local topology from the global topology */
  make_local_top(fplog,dd,fr,top_global,top_local,nsb);

  if (top_global->idef.il[F_CONSTR].nr > 0) {
    make_local_constraints(dd,top_global->idef.il[F_CONSTR].iatoms,
			   ir->nProjOrder,dd->constraints);
    /* Make space for the extra coordinates for constraint communication */
    /* This is not a nice solution.
     * state->natoms is always equal to the global number of atoms.
     * Reallocation will only happen for very small systems
     * with cell sizes close to the cut-off.
     */
    if (dd->nat_tot_con > state_local->natoms)
      srenew(state_local->x,dd->nat_tot_con);
  } else {
    dd->nat_tot_con = dd->nat_tot;
  }

  /* We make the all mdatoms up to nat_tot_con.
   * We could save some work by only setting invmass
   * between nat_tot and nat_tot_con.
   */
  atoms2md(&top_global->atoms,ir,top_global->idef.il[F_ORIRES].nr,
	   dd->nat_tot_con,dd->gatindex,mdatoms);

  if (dd->constraints || top_global->idef.il[F_SETTLE].nr>0)
    init_constraints(fplog,top_global,&top_local->idef.il[F_SETTLE],ir,mdatoms,
		     START(nsb),HOMENR(nsb),
		     ir->eI!=eiSteep,NULL,dd);
  
  /*dump_conf(dd,"ap",state_local->x,state_local->box);*/
  
  /* Calculate the centers of mass of the communicated charge groups.
   * The home charge groups have been done in dd_redistribute_cg.
   */
  calc_cgcm(fplog,dd->ncg_home,dd->ncg_tot,
	    &top_local->blocks[ebCGS],state_local->x,fr->cg_cm);

  inc_nrnb(nrnb,eNR_CGCM,dd->nat_tot-dd->nat_home);
}
