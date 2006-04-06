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

#define dd_c3n 8
static const ivec dd_c3[dd_c3n] =
  {{0,0,0},{1,0,0},{1,1,0},{0,1,0},{0,1,1},{0,0,1},{1,0,1},{1,1,1}};
#define dd_cp3n 4
static const ivec dd_cp3[dd_cp3n] =
  {{0,0,8},{1,3,6},{2,5,6},{3,5,7}};

#define dd_c2n 4
static const ivec dd_c2[dd_c2n] =
  {{0,0,0},{1,0,0},{0,1,0},{1,1,0}};
#define dd_cp2n 2
static const ivec dd_cp2[dd_cp2n] =
  {{0,0,4},{1,2,3}};

#define dd_c1n 2
static const ivec dd_c1[dd_c1n] =
  {{0,0,0},{1,0,0}};
#define dd_cp1n 1
static const ivec dd_cp1[dd_cp1n] =
  {{0,0,2}};

#define dd_index(n,i) ((((i)[ZZ]*(n)[YY] + (i)[YY])*(n)[XX]) + (i)[XX])

static void index2xyz(ivec nc,int ind,ivec xyz)
{
  xyz[XX] = ind % nc[XX];
  xyz[YY] = (ind / nc[XX]) % nc[YY];
  xyz[ZZ] = ind / (nc[YY]*nc[XX]);
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

int dd_nicg(gmx_domdec_t *dd)
{
  return dd->icell[dd->nicell-1].cg1;
}

int dd_ncg_tot(gmx_domdec_t *dd)
{
  return dd->ncg_tot;
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

static void low_dd_move_x(gmx_domdec_t *dd,rvec x[],rvec buf[],
			  int *nreq
#ifdef GMX_MPI
			  ,MPI_Request *req
#endif
			  )
{
  int buf_pos,c,i,cg,j,k,at_index;
  atom_id *cgindex;
  gmx_domdec_comm_t *cc;

  cgindex = dd->cgindex;

  at_index = dd->comm1[0].nat;
  buf_pos = 0;
  for(c=1; c<dd->ncell; c++) {
    cc = &dd->comm0[c];
    k = buf_pos;
    for(i=0; i<cc->ncg; i++) {
      cg = cc->index[i];
      for(j=cgindex[cg]; j<cgindex[cg+1]; j++) {
	if (debug)
	  fprintf(debug,"  %d %d %d %d %d %d %d\n",
		  c,i,cg,cgindex[cg],cgindex[cg+1],j,k);
	copy_rvec(x[j],buf[k]);
	k++;
      }
    }
#ifdef GMX_MPI
    MPI_Isend(buf[buf_pos],cc->nat*sizeof(rvec),MPI_BYTE,
	      DDRANK(dd,cc->cell),c,dd->all,&req[(*nreq)++]);
#endif
    buf_pos += cc->nat;
#ifdef GMX_MPI
    MPI_Irecv(x[at_index],
	      dd->comm1[c].nat*sizeof(rvec),MPI_BYTE,
	      DDRANK(dd,dd->comm1[c].cell),c,dd->all,&req[(*nreq)++]);
#endif
    at_index += dd->comm1[c].nat;
  }
}

void dd_move_x(gmx_domdec_t *dd,rvec x[],rvec buf[])
{
  int nmpi_req;
#ifdef GMX_MPI
  MPI_Request mpi_req[(DD_MAXCELL-1)*2];
  
  nmpi_req = 0;
  low_dd_move_x(dd,x,buf,&nmpi_req,mpi_req);
  
  MPI_Waitall(nmpi_req,mpi_req,MPI_STATUSES_IGNORE);
#endif
}

void dd_start_move_x(gmx_domdec_t *dd,rvec x[],rvec buf[])
{
  if (dd->nmpi_req > 0)
    gmx_fatal(FARGS,"dd_start_move_x called while there are %d outstanding MPI requests",dd->nmpi_req);

#ifdef GMX_MPI  
  low_dd_move_x(dd,x,buf,&dd->nmpi_req,dd->mpi_req);
#endif
}

void dd_finish_move_x(gmx_domdec_t *dd)
{
  if (dd->nmpi_req > 0) {
#ifdef GMX_MPI
    MPI_Waitall(dd->nmpi_req,dd->mpi_req,MPI_STATUSES_IGNORE);
#endif
    dd->nmpi_req = 0;
  }
}

void dd_move_f(gmx_domdec_t *dd,rvec f[],rvec buf[])
{
  int c,i,j,cg,at_index,buf_pos,nreq;
  atom_id *cgindex;
  gmx_domdec_comm_t *cc;
#ifdef GMX_MPI
  MPI_Request mpi_req[DD_MAXCELL*2];

  cgindex = dd->cgindex;

  at_index = dd->comm1[0].nat;
  buf_pos = 0;
  nreq = 0;
  for(c=1; c<dd->ncell; c++) {
    cc = &dd->comm1[c];
    MPI_Isend(f[at_index],
	      cc->nat*sizeof(rvec),MPI_BYTE,
	      DDRANK(dd,cc->cell),c,dd->all,&mpi_req[nreq++]);
    at_index += cc->nat;

    cc = &dd->comm0[c];
    MPI_Irecv(buf[buf_pos],cc->nat*sizeof(rvec),MPI_BYTE,
	      DDRANK(dd,cc->cell),c,dd->all,&mpi_req[nreq++]);
    buf_pos += cc->nat;
  }

  /* Also tried interleaving the Waits and the summing,
   * but that was a lot slower (on Infiniband) B.H.
   */
  MPI_Waitall(nreq,mpi_req,MPI_STATUSES_IGNORE);
#endif

  buf_pos = 0;
  for(c=1; c<dd->ncell; c++) {
    cc = &dd->comm0[c];
    for(i=0; i<cc->ncg; i++) {
      cg = cc->index[i];
      for(j=cgindex[cg]; j<cgindex[cg+1]; j++)
	rvec_inc(f[j],buf[buf_pos++]);
    }
  }
}

static void dd_collect_cg(gmx_domdec_t *dd)
{
  int buf2[2],i;
  gmx_domdec_comm_t *cc;

  cc = &dd->comm1[0];

  buf2[0] = cc->ncg;
  buf2[1] = cc->nat;

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
      fprintf(debug,"Charge group distribution: ");
      for(i=0; i<dd->nnodes; i++)
	fprintf(debug," %d",dd->ma.ncg[i]);
      fprintf(debug,"\n");
    }
  }
  
  /* Collect the charge group indices the master */
#ifdef GMX_MPI
  MPI_Gatherv(cc->index_gl,cc->ncg*sizeof(int),MPI_BYTE,
	      dd->ma.cg,dd->ma.ibuf,dd->ma.ibuf+dd->nnodes,MPI_BYTE,
	      DDMASTERRANK(dd),dd->all);
#endif

  dd->bMasterHasAllCG = TRUE;
}

void dd_collect_vec(gmx_domdec_t *dd,t_block *cgs,rvec *lv,rvec *v)
{
  int  n,i,c,a;
  rvec *buf;
#ifdef GMX_MPI
  MPI_Request mpi_req;
#endif

  if (!dd->bMasterHasAllCG)
    dd_collect_cg(dd);

  if (!DDMASTER(dd)) {
#ifdef GMX_MPI
    MPI_Send(lv,dd->comm1[0].nat*sizeof(rvec),MPI_BYTE,DDMASTERRANK(dd),
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
    MPI_Recv(lv,dd->comm1[0].nat*sizeof(rvec),MPI_BYTE,DDMASTERRANK(dd),
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

static void make_local_indices(gmx_domdec_t *dd,int *cgindex,bool bNonHome)
{
  int c0,c1,c,i,ncg,nat,cg,a_gl;
  gmx_domdec_comm_t *cc;

  if (bNonHome) {
    c0 = 1;
    c1 = dd->ncell;
  } else {
    c0 = 0;
    c1 = 1;
  }

  dd->ncg_tot = 0;
  for(c=0; c<c1; c++)
    dd->ncg_tot += dd->comm1[c].ncg;
  
  /* Make local charge group index */
  if (dd->ncg_tot+1 > dd->cgindex_nalloc) {
    dd->cgindex_nalloc = over_alloc(dd->ncg_tot) + 1;
    srenew(dd->cgindex,dd->cgindex_nalloc);
  }

  if (bNonHome) {
    ncg = dd->comm1[0].ncg;
  } else {
    ncg = 0;
    dd->cgindex[ncg] = 0;
    dd->nat_tot = 0;
  }
  nat = dd->nat_tot;

  /* Count the total numbers of known atoms on this node for allocation */
  for(c=c0; c<c1; c++)
    dd->nat_tot += dd->comm1[c].nat;
  if (dd->nat_tot > dd->gatindex_nalloc) {
    dd->gatindex_nalloc = over_alloc(dd->nat_tot);
    srenew(dd->gatindex,dd->gatindex_nalloc);
  }

  for(c=c0; c<c1; c++) {
    cc = &dd->comm1[c];
    for(i=0; i<cc->ncg; i++) {
      cg = cc->index_gl[i];
      dd->cgindex[ncg+1] = dd->cgindex[ncg] + cgindex[cg+1] - cgindex[cg];
      if (debug)
	fprintf(debug," %d %d %d %d\n",
		ncg,cg,dd->cgindex[ncg],dd->cgindex[ncg+1]);
      for(a_gl=cgindex[cg]; a_gl<cgindex[cg+1]; a_gl++) {
	/* Make local to global atom index */
	dd->gatindex[nat] = a_gl;
	/* Make global to local atom index */
	dd->ga2la[a_gl].cell = c;
	dd->ga2la[a_gl].a    = nat;
	nat++;
      }
      ncg++;
    }
  }
}

static void clear_local_indices(gmx_domdec_t *dd)
{
  int i;

  for(i=0; i<dd->nat_tot; i++)
    dd->ga2la[dd->gatindex[i]].cell = -1;
  if (dd->constraints)
    clear_local_constraint_indices(dd);
}

static void distribute_cg(FILE *fplog,matrix box,t_block *cgs,rvec pos[],
			  gmx_domdec_t *dd)
{
  int **tmp_ind=NULL,*tmp_nalloc=NULL;
  int  i,icg,ai,k,k0,k1,d;
  rvec g_inv,cg_cm;
  ivec ind;
  real nrcg,inv_ncg;
  atom_id *cga,*cgindex;

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
  int i,buf2[2];
  gmx_domdec_comm_t *cc;

  clear_local_indices(dd);

  if (DDMASTER(dd)) {
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

  cc = &dd->comm1[0];
  cc->ncg = buf2[0];
  cc->nat = buf2[1];
  if (cc->ncg > cc->nalloc) {
    cc->nalloc = over_alloc(cc->ncg);
    srenew(cc->index_gl,cc->nalloc);
  }
  if (DDMASTER(dd)) {
    for(i=0; i<dd->nnodes; i++) {
      dd->ma.ibuf[i] = dd->ma.ncg[i]*sizeof(int);
      dd->ma.ibuf[dd->nnodes+i] = dd->ma.index[i]*sizeof(int);
    }
  }

#ifdef GMX_MPI
  MPI_Scatterv(dd->ma.cg,dd->ma.ibuf,dd->ma.ibuf+dd->nnodes,MPI_BYTE,
	       cc->index_gl,cc->ncg*sizeof(int),MPI_BYTE,
	       DDMASTERRANK(dd),dd->all);
#endif

  /* Make local charge group index for the home charge groups */
  make_local_indices(dd,cgs->index,FALSE);
  
  if (debug) {
    fprintf(debug,"Home charge groups:\n");
    for(i=0; i<cc->ncg; i++) {
      fprintf(debug," %d",cc->index_gl[i]);
      if (i % 10 == 9) 
	fprintf(debug,"\n");
    }
    fprintf(debug,"\n");
  }

  dd->bMasterHasAllCG = TRUE;
}

static void dd_calc_cgcm_home(gmx_domdec_t *dd,t_block *cgs,rvec x[],
			      rvec cg_cm[])
{
  int  ncg,*ind_gl,*cgindex,icg,nrcg,i,j,k,d;
  real inv_ncg;
  
  ncg    = dd->comm1[0].ncg;
  ind_gl = dd->comm1[0].index_gl;
  
  cgindex = cgs->index;

  k = 0;
  for(icg=0; icg<ncg; icg++) {
    i = ind_gl[icg];
    nrcg = cgindex[i+1] - cgindex[i];
    if (nrcg == 1) {
      copy_rvec(x[k++],cg_cm[icg]);
    } else {
      inv_ncg = 1.0/nrcg;
      
      clear_rvec(cg_cm[icg]);
      for(j=0; j<nrcg; j++)
	rvec_inc(cg_cm[icg],x[k++]);
      for(d=0; (d<DIM); d++)
	cg_cm[icg][d] = inv_ncg*cg_cm[icg][d];
    }
  }
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
			       int nnbc,int *ind,int *cg,int *buf,
			       rvec *cg_cm)
{
  int c,icg;
  int local_pos,buf_pos[DD_MAXNBCELL];
  
  for(c=0; c<nnbc; c++)
    buf_pos[c] = ind[c];

  local_pos = 0;

  for(icg=0; icg<ncg; icg++) {
    c = cell[icg];
    if (c == local) {
      /* Locally compact the arrays */
      cg[local_pos] = cg[icg];
      copy_rvec(cg_cm[icg],cg_cm[local_pos]);
      local_pos++;
    } else {
      /* Copy to the communication buffer */
      buf[buf_pos[c]++] = cg[icg];
    }
  }

  return local_pos;
}

static void dd_redistribute_cg(FILE *fplog,
			       gmx_domdec_t *dd,t_block *gcgs,
			       t_state *state,rvec cg_cm[],
			       t_nrnb *nrnb)
{
#define dd_nbindex(ms0,ns,i) ((((i)[ZZ] - ms0[ZZ])*ns[YY] + (i)[YY] - ms0[YY])*ns[XX] + (i)[XX] - ms0[XX])

  int  *cell,*buf_cg;
  int  ncg[DD_MAXNBCELL],nat[DD_MAXNBCELL];
  int  ncg_r[DD_MAXNBCELL],nat_r[DD_MAXNBCELL],nat_r_max;
  int  nnbc,nvec,c,i,icg,ai,k,k0,k1,d,x,y,z,nreq,nreq2,ncg_new;
  int  local=-1,dest[DD_MAXNBCELL],src[DD_MAXNBCELL];
  int  sbuf[DD_MAXNBCELL*2],rbuf[DD_MAXNBCELL*2];
  int  buf_cg_ind[DD_MAXNBCELL+1],buf_vs_ind[DD_MAXNBCELL+1],buf_vr_size;
  int  local_pos,buf_pos[DD_MAXNBCELL],buf_pos_r;
  rvec inv_box;
  ivec ms0,ms1,ns,ind,dev,vlocal,vdest,vsrc;
  real nrcg,inv_ncg;
  atom_id *cga,*cgindex;
  gmx_domdec_comm_t *cc;
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
  index2xyz(dd->nc,dd->nodeid,vlocal);
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

  cc = &dd->comm1[0];

  cgindex = dd->cgindex;

  /* Compute the center of geometry for all home charge groups
   * and put them in the box.
   */
  for(icg=0; icg<cc->ncg; icg++) {
    k0      = cgindex[icg];
    k1      = cgindex[icg+1];
    nrcg    = k1 - k0;
    if (nrcg == 1) {
      copy_rvec(state->x[k0],cg_cm[icg]);
    }
    else {
      inv_ncg = 1.0/nrcg;
      
      clear_rvec(cg_cm[icg]);
      for(k=k0; (k<k1); k++)
	rvec_inc(cg_cm[icg],state->x[k]);
      for(d=0; (d<DIM); d++)
	cg_cm[icg][d] = inv_ncg*cg_cm[icg][d];
    }
    for(d=DIM-1; d>=0; d--) {
      /* Determine the cell shift.
       * This does not involve pbc shifts as charge groups
       * should not have been put in the box between repartitioning.
       */
      /* We have to do this in two steps to obtain the correct rounding */
      ind[d] = cg_cm[icg][d]*inv_box[d] + 2;
      ind[d] -= 2;
      dev[d] = ind[d] - vlocal[d];
      if (dev[d] < -1 || dev[d] > 1)
	gmx_fatal(FARGS,"Charge group %d (local index %d) moved more than one cell: %d %d %d, coords %f %f %f",
		  cc->index_gl[icg],icg,dev[XX],dev[YY],dev[ZZ],
		  cg_cm[icg][XX],cg_cm[icg][YY],cg_cm[icg][ZZ]);
      if (ind[d] < 0) {
	ind[d] += dd->nc[d];
	for(k=k0; (k<k1); k++)  {
	  rvec_inc(state->x[k],state->box[d]);
	}
	rvec_inc(cg_cm[icg],state->box[d]);
      } else if (ind[d] >= dd->nc[d]) {
	ind[d] -= dd->nc[d];
	for(k=k0; (k<k1); k++)  {
	  rvec_dec(state->x[k],state->box[d]);
	}
	rvec_dec(cg_cm[icg],state->box[d]);
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
		cc->index_gl[icg],icg,
		i,dest[c],dd_index(dd->nc,ind),
		dev[XX],dev[YY],dev[ZZ],
		cg_cm[icg][XX],cg_cm[icg][YY],cg_cm[icg][ZZ]);
    */
    cell[icg] = i;
    ncg[i] ++;
    nat[i] += nrcg;
  }

  inc_nrnb(nrnb,eNR_RESETX,cc->nat);
  
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

  local_pos = compact_and_copy_vec(cc->ncg,cell,local,cgindex,
				   nnbc,buf_pos,state->x,dd->buf_vs);
  if (state->v)
    compact_and_copy_vec(cc->ncg,cell,local,cgindex,
			 nnbc,buf_pos,state->v,dd->buf_vs);
  if (state->sd_X) {
    compact_and_copy_vec(cc->ncg,cell,local,cgindex,
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

  if (ncg_new > cc->nalloc) {
    cc->nalloc = over_alloc(ncg_new);
    srenew(cc->index_gl,cc->nalloc);
  }

  clear_local_indices(dd);

  local_pos = compact_and_copy_cg(cc->ncg,cell,local,
				  nnbc,buf_cg_ind,cc->index_gl,buf_cg,
				  cg_cm);
  
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
	MPI_Irecv(cc->index_gl+local_pos,ncg_r[c]*sizeof(int),MPI_BYTE,
		  DDRANK(dd,src[c]),c,dd->all,&mpi_req[nreq++]);
#endif
      local_pos += ncg_r[c];
    }
  }
  
  cc->ncg = ncg[local];
  cc->nat = nat[local];
  for(c=0; c<nnbc; c++){
    if (c != local) {
      cc->ncg += ncg_r[c];
      cc->nat += nat_r[c];
    }
  }
  
#ifdef GMX_MPI
  MPI_Waitall(nreq,mpi_req,MPI_STATUSES_IGNORE);
#endif
  
  dd->bMasterHasAllCG = FALSE;

  make_local_indices(dd,gcgs->index,FALSE);

  /* Compute the center of mass for the new home charge groups.
   * We could also communicate this, but as we will only have
   * a few new ones after nstlist steps it will be faster to recompute.
   */
  for(icg=ncg[local]; icg<cc->ncg; icg++) {
    k0      = cgindex[icg];
    k1      = cgindex[icg+1];
    nrcg    = k1 - k0;
    if (nrcg == 1) {
      copy_rvec(state->x[k0],cg_cm[icg]);
    }
    else {
      inv_ncg = 1.0/nrcg;
      
      clear_rvec(cg_cm[icg]);
      for(k=k0; (k<k1); k++)
	rvec_inc(cg_cm[icg],state->x[k]);
      for(d=0; (d<DIM); d++)
	cg_cm[icg][d] = inv_ncg*cg_cm[icg][d];
    }
  }

  inc_nrnb(nrnb,eNR_CGCM,cc->nat-nat[local]);

  if (debug)
    fprintf(debug,"Finished repartitioning\n");
}

void setup_dd_grid(FILE *fplog,matrix box,gmx_domdec_t *dd)
{
  int  d,i,j,m;
  ivec xyz,tmp,h,s;
  int  ncell,ncellp;
  ivec *dd_c,dd_cp[DD_MAXICELL];
  gmx_domdec_ns_t *icell;

  index2xyz(dd->nc,dd->nodeid,xyz);

  dd->ndim = 0;
  for(d=0; d<DIM; d++) {
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
  index2xyz(dd->nc,dd->nodeid,h);
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
    dd->comm1[i].cell = dd_index(dd->nc,s);
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
    dd->comm0[i].cell = dd_index(dd->nc,s);
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

static void low_make_at2iatoms(t_idef *idef,
			       gmx_at2iatoms_t *ga2iatoms,bool bAssign)
{
  int ftype,nral,i,ind_at;
  t_ilist *il;
  t_iatom *ia;
  atom_id a;
  
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
	  ga2iatoms[a].ftype [ga2iatoms[a].n] = ftype;
	  ga2iatoms[a].iatoms[ga2iatoms[a].n] = ia;
	}
	ga2iatoms[a].n++;
      }
    }
  }
}

static gmx_at2iatoms_t *make_at2iatoms(int natoms,t_idef *idef)
{
  int i;
  gmx_at2iatoms_t *ga2iatoms;

  snew(ga2iatoms,natoms);
  
  /* Count the interactions */
  low_make_at2iatoms(idef,ga2iatoms,FALSE);

  for(i=0; i<natoms; i++) {
    snew(ga2iatoms[i].ftype, ga2iatoms[i].n);
    snew(ga2iatoms[i].iatoms,ga2iatoms[i].n);
    ga2iatoms[i].n = 0;
  }

  /* Store the interactions */
  low_make_at2iatoms(idef,ga2iatoms,TRUE);

  return ga2iatoms;
}

gmx_domdec_t *init_domain_decomposition(FILE *fplog,t_commrec *cr,
					ivec nc,int ncg,int natoms,
					matrix box,
					t_idef *idef,
					bool bDynamics)
{
  gmx_domdec_t *dd;
  int a;

  fprintf(fplog,"Will do domain decomposition\n");
  
  snew(dd,1);

  copy_ivec(nc,dd->nc);
  dd->nnodes = cr->nnodes - cr->npmenodes;
  dd->nodeid = cr->nodeid;
#ifdef GMX_MPI
  dd->all    = cr->mpi_comm_mygroup;
#endif
  if (DDMASTER(dd)) {
    snew(dd->ma.ncg,dd->nnodes);
    snew(dd->ma.index,dd->nnodes+1);
    snew(dd->ma.cg,ncg);
    snew(dd->ma.nat,dd->nnodes);
    snew(dd->ma.ibuf,dd->nnodes*2);
  }
  dd->comm1[0].cell = dd->nodeid;

  dd->ga2iatoms = make_at2iatoms(natoms,idef);
  
  snew(dd->ga2la,natoms*sizeof(gmx_ga2la_t));
  for(a=0; a<natoms; a++)
    dd->ga2la[a].cell = -1;
  
  if (idef->il[F_CONSTR].nr > 0)
    dd->constraints = init_domdec_constraints(natoms,idef,bDynamics);

  return dd;
}

static void setup_dd_communication(FILE *fplog,gmx_domdec_t *dd,t_block *gcgs,
				   rvec cg_cm[],matrix box,real r_comm)
{
  int *cgindex,*index_gl,ncg,c,d,i,cg;
  gmx_domdec_comm_t *cc;
  ivec xyz,shift;
  rvec corner[DD_MAXCELL],corner_c;
  real r_comm2,r2;
  bool bUse;
  int  buf_cs[DD_MAXCELL*2],buf_cr[DD_MAXCELL*2],nreq;
#ifdef GMX_MPI
  MPI_Request mpi_req[DD_MAXCELL*2];
#endif

  if (debug)
    fprintf(debug,"Setting up DD communication\n");

  cgindex = dd->cgindex;

  for(c=1; c<dd->ncell; c++) {
    cc = &dd->comm0[c];
    index2xyz(dd->nc,dd->nodeid,xyz);
    /* This should be changed for non-uniform grids */
    for(d=0; d<DIM; d++)
      corner[c][d] = xyz[d]*box[d][d]/dd->nc[d];
    cc->ncg = 0;
    cc->nat = 0;
  }

  if (debug)
    fprintf(debug,"Cell %d, ncg %d",dd->nodeid,dd->comm1[0].ncg);

  r_comm2 = sqr(r_comm);
  ncg      = dd->comm1[0].ncg;
  index_gl = dd->comm1[0].index_gl;
  nreq = 0;
  for(c=1; c<dd->ncell; c++) {
    cc = &dd->comm0[c];
    copy_ivec(dd->shift[c],shift);
    copy_rvec(corner[c],corner_c);
    for(i=0; i<ncg; i++) {
      /* This code currently only works for shifts 0 and +1 */
      /*
      bUse = TRUE;
      for(d=0; d<DIM; d++) {
	if (dd->shift[c][d] == 1 &&
	    cg_cm[i][d] - corner[c][d] >= r_comm)
	  bUse = FALSE;
      }
      */
      /* This code currently only works for shifts 0 and +1 */
      r2 = 0;
      for(d=0; d<DIM; d++)
	r2 += shift[d]*sqr(cg_cm[i][d] - corner_c[d]);
      if (r2 < r_comm2) {
	/* Make an index to the local charge groups */
	if (cc->ncg >= cc->nalloc) {
	  cc->nalloc += CG_ALLOC_SIZE;
	  srenew(cc->index_gl,cc->nalloc);
	  srenew(cc->index,cc->nalloc);
	}
	cg = index_gl[i];
	cc->index_gl[cc->ncg] = cg;
	cc->index[cc->ncg] = i;
	cc->ncg++;
	cc->nat += cgindex[i+1] - cgindex[i];
      }
    }

    if (debug)
      fprintf(debug,", %d > %d",cc->cell,cc->ncg);

    /* Send the number of charge groups and atom to communicate */
    buf_cs[c*2]   = cc->ncg;
    buf_cs[c*2+1] = cc->nat;
#ifdef GMX_MPI
    MPI_Isend(buf_cs+c*2,2*sizeof(int),MPI_BYTE,
	      DDRANK(dd,cc->cell),c,dd->all,&mpi_req[nreq++]);
    MPI_Irecv(buf_cr+c*2,2*sizeof(int),MPI_BYTE,
	      DDRANK(dd,dd->comm1[c].cell),c,dd->all,&mpi_req[nreq++]);
#endif
  }

#ifdef GMX_MPI
  MPI_Waitall(nreq,mpi_req,MPI_STATUSES_IGNORE);
#endif

  /* We could send out the charge group indices before,
   * but apparently some MPI implementations do not allow
   * multiple messages in flight between the same pairs of nodes.
   */
  nreq = 0;
  for(c=1; c<dd->ncell; c++) {
    cc = &dd->comm0[c];
#ifdef GMX_MPI
    if (cc->ncg > 0)
      MPI_Isend(cc->index_gl,cc->ncg*sizeof(int),MPI_BYTE,
		DDRANK(dd,cc->cell),c,dd->all,&mpi_req[nreq++]);
#endif

    /* Receive the number of charge groups and atoms to communicate */
    cc = &dd->comm1[c];
    cc->ncg = buf_cr[c*2];
    cc->nat = buf_cr[c*2+1];
    if (debug)
      fprintf(debug," %d < %d",cc->cell,cc->ncg);
    if (cc->ncg > cc->nalloc) {
      cc->nalloc = over_alloc(cc->ncg);
      srenew(cc->index_gl,cc->nalloc);
    }
#ifdef GMX_MPI
    if (cc->ncg > 0)
      MPI_Irecv(cc->index_gl,cc->ncg*sizeof(int),MPI_BYTE,
		DDRANK(dd,cc->cell),c,dd->all,&mpi_req[nreq++]);
#endif
  }
  if (debug)
    fprintf(debug,"\n");
  
#ifdef GMX_MPI
  MPI_Waitall(nreq,mpi_req,MPI_STATUSES_IGNORE);
#endif

  /* Make local indices for the non-home charge groups */
  make_local_indices(dd,gcgs->index,TRUE);

  if (debug) {
    fprintf(debug,"Finished setting up DD communication, cells:");
    for(c=0; c<dd->ncell; c++)
      fprintf(debug," %d",dd->comm1[c].ncg);
    fprintf(debug,"\n");
  }
}

static void make_local_ilist(gmx_domdec_t *dd,t_functype ftype,
			     t_ilist *il,t_ilist *lil)
{
  int nral,nhome,i,j,n;
  t_iatom *ia,*lia;
  atom_id a;
  
  nral = NRAL(ftype);

  if (lil->iatoms == NULL)
    /* In most cases we could do with far less memory */
    snew(lil->iatoms,il->nr);

  nhome = dd->comm1[0].nat;

  n   = 0;
  ia  = il->iatoms;
  lia = lil->iatoms;
  for(i=0; i<il->nr; i+=1+nral) {
    ia = il->iatoms + i;
    a = dd->ga2la[ia[1]].a;
    if (a >= 0 && a < nhome) {
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
  gmx_domdec_comm_t *cc;
  gmx_at2iatoms_t *ga2ia;
  t_iatom *iatoms,tiatoms[1+MAXATOMLIST],*liatoms;
  bool bUse;
  ivec shift_min;
  gmx_ga2la_t *ga2la;
  t_ilist *il;

  /* Clear the counts */
  for(ftype=0; ftype<F_NRE; ftype++)
    idef->il[ftype].nr = 0;

  for(i=0; i<dd->nat_tot; i++) {
    /* Get the global atom number */
    gat = dd->gatindex[i];
    ga2ia = &dd->ga2iatoms[gat];
    /* Check all interactions assigned to this atom */
    for(j=0; j<ga2ia->n; j++) {
      ftype  = ga2ia->ftype[j];
      iatoms = ga2ia->iatoms[j];
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
	if (debug) {
	  fprintf(debug,"ftype %d ip %d at",ftype,liatoms[0]);
	  for(k=1; k<=nral; k++)
	    fprintf(debug," %d (%d)",liatoms[k],iatoms[k]);
	  fprintf(debug,"\n");
	}
	il->nr += 1+nral;
      }
    }
  }
}

static void make_local_exclusions(gmx_domdec_t *dd,
				  t_block *excls,t_block *lexcls)
{
  int n,la,a,i;
  gmx_ga2la_t *ga2la;
  
  if (dd->nat_tot+1 > lexcls->nalloc_index) {
    lexcls->nalloc_index = over_alloc(dd->nat_tot)+1;
    srenew(lexcls->index,lexcls->nalloc_index);
  }
  
  n = 0;
  for(la=0; la<dd->nat_tot; la++) {
    lexcls->index[la] = n;
    a = dd->gatindex[la];
    for(i=excls->index[a]; i<excls->index[a+1]; i++) {
      ga2la = &dd->ga2la[excls->a[i]];
      if (ga2la->cell != -1) {
	if (n >= lexcls->nalloc_a) {
	  lexcls->nalloc_a += EXCLS_ALLOC_SIZE;
	  srenew(lexcls->a,lexcls->nalloc_a);
	}
	lexcls->a[n++] = ga2la->a;
      }
    }
  }
  lexcls->index[la] = n;
  lexcls->nr = la;
  lexcls->nra = n;
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
  int i,c;

  lcgs->nr    = dd->ncg_tot;
  lcgs->index = dd->cgindex;
  lcgs->nra   = dd->nat_tot;
}

static void set_cg_boundaries(gmx_domdec_t *dd,t_block *lcgs)
{
  int cgend[DD_MAXCELL+1],c;

  dd->ncg_tot = 0;
  cgend[0] = 0;
  for(c=0; c<dd->ncell; c++)
    cgend[c+1] = cgend[c] + dd->comm1[c].ncg;

  dd->ncg_tot = cgend[dd->ncell];

  for(c=0; c<dd->nicell; c++) {
    dd->icell[c].cg1  = cgend[c+1];
    dd->icell[c].jcg0 = cgend[dd->icell[c].j0];
    dd->icell[c].jcg1 = cgend[dd->icell[c].j1];
  }
}

static void dd_update_ns_border(gmx_domdec_t *dd,t_nsborder *nsb)
{
  nsb->nodeid    = 0 /* dd->nodeid */;
  nsb->cgtotal   = dd_ncg_tot(dd);
  //nsb->natoms    = dd->gl.nat[dd->nodeid];
  nsb->index[nsb->nodeid]  = 0;
  nsb->homenr[nsb->nodeid] = dd->comm1[0].nat;
  nsb->cgload[nsb->nodeid] = nsb->cgtotal;
}

static void make_local_top(FILE *fplog,gmx_domdec_t *dd,
			   t_topology *top,t_topology *ltop,
			   t_nsborder *nsb)
{
  int a,i,cg,eb;

  if (debug)
    fprintf(debug,"Making local topology\n");

  ltop->name  = top->name;
  make_local_cgs(dd,&ltop->blocks[ebCGS]);

  make_local_exclusions(dd,&top->blocks[ebEXCLS],&ltop->blocks[ebEXCLS]);

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

  set_cg_boundaries(dd,&ltop->blocks[ebCGS]);

  dd_update_ns_border(dd,nsb);
}

static void make_local_forcerec(gmx_domdec_t *dd,t_forcerec *fr)
{
  int c,cg,i;
  gmx_domdec_comm_t *cc;

  fr->cg0 = 0;
  fr->hcg = dd_ncg_tot(dd);

  if (fr->solvent_opt == esolNO) {
    /* Since all entries are identical, we can use the global array */
    fr->solvent_type = fr->solvent_type_global;
  } else {
    if (dd->ncg_tot > fr->solvent_type_nalloc) {
      fr->solvent_type_nalloc = over_alloc(dd->ncg_tot);
      srenew(fr->solvent_type,fr->solvent_type_nalloc);
    }
    cg = 0;
    for(c=0; c<dd->ncell; c++) {
      cc = &dd->comm1[c];
      for(i=0; i<cc->ncg; i++)
	fr->solvent_type[cg++] = fr->solvent_type_global[cc->index_gl[i]];
    }
  }
}

static void dump_conf(gmx_domdec_t *dd,
		      char *name,rvec *x,matrix box)
{
  char str[STRLEN];
  FILE *fp;
  int c,i,j;
  
  sprintf(str,"%s%d.pdb",name,dd->nodeid);
  fp = fopen(str,"w");
  fprintf(fp,"CRYST1%9.3f%9.3f%9.3f%7.2f%7.2f%7.2f P 1           1\n",
	  10*norm(box[XX]),10*norm(box[YY]),10*norm(box[ZZ]),
	  90.0,90.0,90.0);
  i = 0;
  for(c=0; c<dd->ncell; c++) {
    for(j=0; j<dd->comm1[c].nat; j++) {
      fprintf(fp,"%-6s%5u  %-4.4s%3.3s %c%4d    %8.3f%8.3f%8.3f%6.2f%6.2f\n",
	    "ATOM",dd->gatindex[i]+1,"C","ALA",' ',i+1,
	    10*x[i][XX],10*x[i][YY],10*x[i][ZZ],
	    1.0,1.0*c/(dd->ncell-1.0));
      i++;
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
  if (bMasterState) {
    get_cg_distribution(fplog,dd,
			&top_global->blocks[ebCGS],
			state_global->box,state_global->x);
    
    dd_distribute_state(dd,&top_global->blocks[ebCGS],
			state_global,state_local);
    
    dd_calc_cgcm_home(dd,&top_global->blocks[ebCGS],
		      state_local->x,fr->cg_cm);
  } else {
    dd_redistribute_cg(fplog,dd,&top_global->blocks[ebCGS],
		       state_local,fr->cg_cm,nrnb);
  }
  
  setup_dd_communication(fplog,dd,&top_global->blocks[ebCGS],
			 fr->cg_cm,state_local->box,fr->rlistlong);
  
  dd_start_move_x(dd,state_local->x,buf);
  
  make_local_top(fplog,dd,top_global,top_local,nsb);

  dd_finish_move_x(dd);

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
  atoms2md(fplog,NULL,&top_global->atoms,ir->opts.nFreeze,
	   ir->eI,ir->delta_t,ir->bd_fric,ir->opts.tau_t,
	   ir->efep!=efepNO,dd->nat_tot_con,dd->gatindex,mdatoms,FALSE);

  make_local_forcerec(dd,fr);

  if (dd->constraints)
    init_constraints(fplog,top_global,ir,mdatoms,
		     START(nsb),HOMENR(nsb),
		     ir->eI!=eiSteep,NULL,dd);
  
  //dump_conf(dd,"ap",state_local->x,state_local->box);
  
  /* Calculate the centers of mass of the communicated charge groups.
   * The home charge groups have been done in dd_redistribute_cg.
   */
  calc_cgcm(fplog,dd->comm1[0].ncg,dd->ncg_tot,
	    &top_local->blocks[ebCGS],state_local->x,fr->cg_cm);

  inc_nrnb(nrnb,eNR_CGCM,dd->nat_tot-dd->comm1[0].nat);
}
