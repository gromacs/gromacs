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
#include "pdbio.h"
#include "futil.h"
#include "pme.h"

#ifdef GMX_MPI
#include <mpi.h>
#endif

typedef struct gmx_domdec_master {
  /* The cell boundaries */
  real **cell_x;
  /* The global charge group division */
  int  *ncg;     /* Number of home charge groups for each node */
  int  *index;   /* Index of nnodes+1 into cg */
  int  *cg;      /* Global charge group index */
  int  *nat;     /* Number of home atoms for each node. */
  int  *ibuf;    /* Buffer for communication */
} gmx_domdec_master_t;

typedef struct {
  /* The numbers of charge groups and atoms to send and receive */
  int nsend[DD_MAXICELL+2];
  int nrecv[DD_MAXICELL+2];
  /* The charge groups to send */
  int *index;
  int nalloc;
} gmx_domdec_ind_t;

typedef struct gmx_domdec_comm {
  /* The width of the communicated boundaries */
  real distance;

  /* The indices to communicate */
  gmx_domdec_ind_t ind[DIM];

  /* Communication buffer for general use */
  int  *buf_int;
  int  nalloc_int;

  /* Communication buffers for local redistribution */
  int  **cggl_flag;
  int  cggl_flag_nalloc[DIM*2];
  rvec **cgcm_state;
  int  cgcm_state_nalloc[DIM*2];
  rvec *buf_vr;
  int  nalloc_vr;
} gmx_domdec_comm_t;

/* The size per charge group of the cggl_flag buffer in gmx_domdec_comm_t */
#define DD_CGIBS 2

#define CG_ALLOC_SIZE     1000

/* Cell permutation required to obtain consecutive charge groups
 * for neighbor searching.
 */
static const int cell_perm[3][4] = { {0,0,0,0},{1,0,0,0},{3,0,1,2} };

/* The DD cell order */
static const ivec dd_co[DD_MAXCELL] =
  {{0,0,0},{1,0,0},{1,1,0},{0,1,0},{0,1,1},{0,0,1},{1,0,1},{1,1,1}};

/* The 3D setup */
#define dd_c3n  8
#define dd_cp3n 4
static const ivec dd_cp3[dd_cp3n] = {{0,0,8},{1,3,6},{2,5,6},{3,5,7}};

/* The 2D setup */
#define dd_c2n  4
#define dd_cp2n 2
static const ivec dd_cp2[dd_cp2n] = {{0,0,4},{1,3,4}};

/* The 1D setup */
#define dd_c1n  2
#define dd_cp1n 1
static const ivec dd_cp1[dd_cp1n] = {{0,0,2}};

static bool bDDDump;

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
  int icell,d,dim;

  icell = 0;
  while (icg >= dd->icell[icell].cg1)
    icell++;

  if (icell == 0)
    *jcg0 = icg;
  else if (icell < dd->nicell)
    *jcg0 = dd->icell[icell].jcg0;
  else
    gmx_fatal(FARGS,"DD icg %d out of range: icell (%d) >= nicell (%d)",
	      icg,icell,dd->nicell);

  *jcg1 = dd->icell[icell].jcg1;

  for(d=0; d<dd->ndim; d++) {
    dim = dd->dim[d];
    shift0[dim] = dd->icell[icell].shift0[dim];
    shift1[dim] = dd->icell[icell].shift1[dim];
    if (dd->tric_dir[dim] || (dd->bDynLoadBal && d > 0)) {
      /* A conservative approach, this can be optimized */
      shift0[dim] -= 1;
      shift1[dim] += 1;
    }
  }
}

void dd_sendrecv_int(const gmx_domdec_t *dd,
		     int ddim,int direction,
		     int *buf_s,int n_s,
		     int *buf_r,int n_r)
{
#ifdef GMX_MPI
  int rank_s,rank_r;
  MPI_Status stat;

  rank_s = dd->neighbor[ddim][direction==ddForward ? 0 : 1];
  rank_r = dd->neighbor[ddim][direction==ddForward ? 1 : 0];

  MPI_Sendrecv(buf_s,n_s*sizeof(int),MPI_BYTE,rank_s,0,
	       buf_r,n_r*sizeof(int),MPI_BYTE,rank_r,0,
	       dd->all,&stat);
#endif
}

void dd_sendrecv_rvec(const gmx_domdec_t *dd,
		      int ddim,int direction,
		      rvec *buf_s,int n_s,
		      rvec *buf_r,int n_r)
{
#ifdef GMX_MPI
  int rank_s,rank_r;
  MPI_Status stat;

  rank_s = dd->neighbor[ddim][direction==ddForward ? 0 : 1];
  rank_r = dd->neighbor[ddim][direction==ddForward ? 1 : 0];

  MPI_Sendrecv(buf_s[0],n_s*sizeof(rvec),MPI_BYTE,rank_s,0,
	       buf_r[0],n_r*sizeof(rvec),MPI_BYTE,rank_r,0,
	       dd->all,&stat);
#endif
}

void dd_move_x(gmx_domdec_t *dd,matrix box,rvec x[],rvec buf[])
{
  int  ncell,nat_tot,n,dim,i,j;
  int  *index,*cgindex;
  gmx_domdec_comm_t *comm;
  gmx_domdec_ind_t *ind;
  rvec shift;

  comm = dd->comm;
  
  cgindex = dd->cgindex;

  ncell = 1;
  nat_tot = dd->nat_home;
  for(dim=0; dim<dd->ndim; dim++) {
    ind = &comm->ind[dim];
    index = ind->index;
    n = 0;
    for(i=0; i<ind->nsend[ncell]; i++) {
      if (dd->ci[dd->dim[dim]] == 0) {
	/* We need to shift the coordinates */
	copy_rvec(box[dd->dim[dim]],shift);
	for(j=cgindex[index[i]]; j<cgindex[index[i]+1]; j++) {
	  rvec_add(x[j],shift,buf[n]);
	  n++;
	}
      } else {
	for(j=cgindex[index[i]]; j<cgindex[index[i]+1]; j++) {
	  copy_rvec(x[j],buf[n]);
	  n++;
	}
      }
    }
    /* Send and receive the coordinates */
    dd_sendrecv_rvec(dd, dim, ddBackward,
		     buf,       ind->nsend[ncell+1],
		     x+nat_tot, ind->nrecv[ncell+1]);
    nat_tot += ind->nrecv[ncell+1];
    ncell += ncell;
  }
}

void dd_move_f(gmx_domdec_t *dd,rvec f[],rvec buf[],rvec *fshift)
{
  int  ncell,nat_tot,n,dim,i,j;
  int  *index,*cgindex;
  gmx_domdec_comm_t *comm;
  gmx_domdec_ind_t *ind;
  ivec vis;
  int  is;

  comm = dd->comm;
  
  cgindex = dd->cgindex;

  ncell = 1;
  nat_tot = dd->nat_home;
  n = 0;
  ncell = dd->ncell/2;
  nat_tot = dd->nat_tot;
  for(dim=dd->ndim-1; dim>=0; dim--) {
    ind = &comm->ind[dim];
    nat_tot -= ind->nrecv[ncell+1];
    /* Communicate the forces */
    dd_sendrecv_rvec(dd,dim,ddForward,
		     f+nat_tot, ind->nrecv[ncell+1],
		     buf,       ind->nsend[ncell+1]);
    index = ind->index;
    /* Add the received forces */
    n = 0;
    for(i=0; i<ind->nsend[ncell]; i++) {
      if (fshift && dd->ci[dd->dim[dim]] == 0) {
	clear_ivec(vis);
	vis[dd->dim[dim]] = 1;
	is = IVEC2IS(vis);
	for(j=cgindex[index[i]]; j<cgindex[index[i]+1]; j++) {
	  rvec_inc(f[j],buf[n]);
	  /* Add this force to the shift force */
	  rvec_inc(fshift[is],buf[n]);
	  n++;
	}
      } else {
	for(j=cgindex[index[i]]; j<cgindex[index[i]+1]; j++) {
	  rvec_inc(f[j],buf[n]);
	  n++;
	}
      }
    }
    ncell /= 2;
  }
}

static void dd_bcast(gmx_domdec_t *dd,int nbytes,void *data)
{
#ifdef GMX_MPI
  MPI_Bcast(data,nbytes,MPI_BYTE,
	    DDMASTERRANK(dd),dd->all);
#endif
}

static void dd_scatter(gmx_domdec_t *dd,int nbytes,void *src,void *dest)
{
#ifdef GMX_MPI
  MPI_Scatter(src,nbytes,MPI_BYTE,
	      dest,nbytes,MPI_BYTE,
	      DDMASTERRANK(dd),dd->all);
#endif
}

static void dd_gather(gmx_domdec_t *dd,int nbytes,void *src,void *dest)
{
#ifdef GMX_MPI
  MPI_Gather(src,nbytes,MPI_BYTE,
	     dest,nbytes,MPI_BYTE,
	     DDMASTERRANK(dd),dd->all);
#endif
}

static void dd_scatterv(gmx_domdec_t *dd,
			int *scounts,int *disps,void *sbuf,
			int rcount,void *rbuf)
{
#ifdef GMX_MPI
  MPI_Scatterv(sbuf,scounts,disps,MPI_BYTE,
	       rbuf,rcount,MPI_BYTE,
	       DDMASTERRANK(dd),dd->all);
#endif
}

static void dd_gatherv(gmx_domdec_t *dd,
		       int scount,void *sbuf,
		       int *rcounts,int *disps,void *rbuf)
{
#ifdef GMX_MPI
  MPI_Gatherv(sbuf,scount,MPI_BYTE,
	      rbuf,rcounts,disps,MPI_BYTE,
	      DDMASTERRANK(dd),dd->all);
#endif
}

static void dd_collect_cg(gmx_domdec_t *dd)
{
  gmx_domdec_master_t *ma;
  int buf2[2],*ibuf,i;

  ma = dd->ma;

  buf2[0] = dd->ncg_home;
  buf2[1] = dd->nat_home;
  if (DDMASTER(dd)) {
    ibuf = ma->ibuf;
  } else {
    ibuf = NULL;
  }
  /* Collect the charge group and atom counts on the master */
  dd_gather(dd,2*sizeof(int),buf2,ibuf);

  if (DDMASTER(dd)) {
    ma->index[0] = 0;
    for(i=0; i<dd->nnodes; i++) {
      ma->ncg[i] = ma->ibuf[2*i];
      ma->nat[i] = ma->ibuf[2*i+1];
      ma->index[i+1] = ma->index[i] + ma->ncg[i];
      
    }
    /* Make byte counts and indices */
    for(i=0; i<dd->nnodes; i++) {
      ma->ibuf[i] = ma->ncg[i]*sizeof(int);
      ma->ibuf[dd->nnodes+i] = ma->index[i]*sizeof(int);
    }
    if (debug) {
      fprintf(debug,"Initial charge group distribution: ");
      for(i=0; i<dd->nnodes; i++)
	fprintf(debug," %d",ma->ncg[i]);
      fprintf(debug,"\n");
    }
  }
  
  /* Collect the charge group indices on the master */
  dd_gatherv(dd,
	     dd->ncg_home*sizeof(int),dd->index_gl,
	     DDMASTER(dd) ? ma->ibuf : NULL,
	     DDMASTER(dd) ? ma->ibuf+dd->nnodes : NULL,
	     DDMASTER(dd) ? ma->cg : NULL);

  dd->bMasterHasAllCG = TRUE;
}

void dd_collect_vec(gmx_domdec_t *dd,t_block *cgs,rvec *lv,rvec *v)
{
  gmx_domdec_master_t *ma;
  int  n,i,c,a;
  rvec *buf;

  ma = dd->ma;

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
    for(i=ma->index[n]; i<ma->index[n+1]; i++)
      for(c=cgs->index[ma->cg[i]]; c<cgs->index[ma->cg[i]+1]; c++)
	copy_rvec(lv[a++],v[c]);

    /* Use the unused part of lv as a temporary buffer */
    buf = lv + ma->nat[n];

    for(n=0; n<dd->nnodes; n++) {
      if (n != dd->nodeid) {
#ifdef GMX_MPI
	MPI_Recv(buf,ma->nat[n]*sizeof(rvec),MPI_BYTE,DDRANK(dd,n),
		 n,dd->all,MPI_STATUS_IGNORE);
#endif
	a = 0;
	for(i=ma->index[n]; i<ma->index[n+1]; i++)
	  for(c=cgs->index[ma->cg[i]]; c<cgs->index[ma->cg[i]+1]; c++)
	    copy_rvec(buf[a++],v[c]);
      }
    }
  }
}

void dd_collect_state(gmx_domdec_t *dd,t_block *cgs,
		      t_state *state_local,t_state *state)
{
  if (DDMASTER(dd)) {
    copy_mat(state_local->box,state->box);
  }
  dd_collect_vec(dd,cgs,state_local->x,state->x);
  if (state->v)
    dd_collect_vec(dd,cgs,state_local->v,state->v);
  if (state->sd_X)
    dd_collect_vec(dd,cgs,state_local->sd_X,state->sd_X);
}

static void dd_distribute_vec(gmx_domdec_t *dd,t_block *cgs,rvec *v,rvec *lv)
{
  gmx_domdec_master_t *ma;
  int n,i,c,a;

  if (DDMASTER(dd)) {
    ma  = dd->ma;

    for(n=0; n<dd->nnodes; n++) {
      if (n != dd->nodeid) {
	/* Use lv as a temporary buffer */
	a = 0;
	for(i=ma->index[n]; i<ma->index[n+1]; i++)
	  for(c=cgs->index[ma->cg[i]]; c<cgs->index[ma->cg[i]+1]; c++)
	    copy_rvec(v[c],lv[a++]);
	if (a != ma->nat[n])
	  gmx_fatal(FARGS,"Internal error a (%d) != nat (%d)",a,ma->nat[n]);

#ifdef GMX_MPI
	MPI_Send(lv,ma->nat[n]*sizeof(rvec),MPI_BYTE,
		 DDRANK(dd,n),n,dd->all);
#endif
      }
    }
    n = 0;
    a = 0;
    for(i=ma->index[n]; i<ma->index[n+1]; i++)
      for(c=cgs->index[ma->cg[i]]; c<cgs->index[ma->cg[i]+1]; c++)
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
  if (DDMASTER(dd)) {
    copy_mat(state->box,state_local->box);
  }
  dd_bcast(dd,sizeof(state_local->box),state_local->box);
  dd_distribute_vec(dd,cgs,state->x,state_local->x);
  if (state->v)
    dd_distribute_vec(dd,cgs,state->v,state_local->v);
  if (state->sd_X)
    dd_distribute_vec(dd,cgs,state->sd_X,state_local->sd_X);
}

static char dim2char(int dim)
{
  char c='?';

  switch (dim) {
  case XX: c = 'x'; break;
  case YY: c = 'y'; break;
  case ZZ: c = 'z'; break;
  default: gmx_fatal(FARGS,"Unknown dim %d",dim);
  }

  return c;
}

static void write_dd_pdb(char *fn,int step,char *title,t_atoms *atoms,
			 gmx_domdec_t *dd,int natoms,
			 rvec x[],matrix box)
{
  char fname[STRLEN],format[STRLEN];
  FILE *out;
  int  i,ii,resnr,c;
  real b;

  sprintf(fname,"%s_%d_n%d.pdb",fn,step,dd->nodeid);

  sprintf(format,"%s%s\n",pdbformat,"%6.2f%6.2f");

  out = ffopen(fname,"w");

  fprintf(out,"TITLE     %s\n",title);
  gmx_write_pdb_box(out,box);
  for(i=0; i<natoms; i++) {
    ii = dd->gatindex[i];
    resnr = atoms->atom[ii].resnr;
    if (i < dd->nat_tot) {
      c = 0;
      while (i >= dd->cgindex[dd->ncg_cell[c+1]]) {
	c++;
      }
      b = c;
    } else if (i < dd->nat_tot_vsite) {
      b = dd->ncell;
    } else {
      b = dd->ncell + 1;
    }
    fprintf(out,format,"ATOM",(ii+1)%100000,
	    *atoms->atomname[ii],*atoms->resname[resnr],' ',(resnr+1)%10000,
	    10*x[i][XX],10*x[i][YY],10*x[i][ZZ],1.0,b);
  }
  fprintf(out,"TER\n");
  
  fclose(out);
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

  /* Clear the indices without looping over all the atoms in the system */
  for(i=a_start; i<dd->nat_tot; i++)
    dd->ga2la[dd->gatindex[i]].cell = -1;

  clear_local_vsite_indices(dd);

  if (dd->constraints)
    clear_local_constraint_indices(dd);
}

static void set_dd_cell_sizes(gmx_domdec_t *dd,matrix box,bool bMaster)
{
  int  d,j;
  real tot;

  for(d=0; d<DIM; d++) {
    dd->tric_dir[d] = 0;
    for(j=d+1; j<DIM; j++)
      if (box[j][d] != 0) {
	dd->tric_dir[d] = 1;
	if (dd->nc[j] > 1 && dd->nc[d] == 1)
	  gmx_fatal(FARGS,"Domain decomposition has not been implemented for box vectors that have non-zero components in directions that do not use domain decomposition: ncells = %d %d %d, box vector[%d] = %f %f %f",
		    dd->nc[XX],dd->nc[YY],dd->nc[ZZ],
		    j+1,box[j][XX],box[j][YY],box[j][ZZ]);
      }
    if (dd->cell_load[d] == NULL || dd->nc[d] == 1) {
      /* Uniform grid */
      if (bMaster) {
       for(j=0; j<dd->nc[d]+1; j++)
	 dd->ma->cell_x[d][j] =      j*box[d][d]/dd->nc[d];
      } else {
	dd->cell_x0[d] = (dd->ci[d]  )*box[d][d]/dd->nc[d];
	dd->cell_x1[d] = (dd->ci[d]+1)*box[d][d]/dd->nc[d];
      }
    } else {
      /* Load balanced grid */
      tot = 0;
      for(j=0; j<dd->nc[d]; j++)
	tot += 1/dd->cell_load[d][j];
      
      if (bMaster) {
	dd->ma->cell_x[d][0] = 0;
	for(j=0; j<dd->nc[d]; j++)
	  dd->ma->cell_x[d][j+1] =
	    dd->ma->cell_x[d][j] + box[d][d]/(dd->cell_load[d][j]*tot);
      } else {
	dd->cell_x0[d] = 0;
	for(j=0; j<dd->ci[d]; j++)
	  dd->cell_x0[d] += box[d][d]/(dd->cell_load[d][j]*tot);
	dd->cell_x1[d] =
	  dd->cell_x0[d] +  box[d][d]/(dd->cell_load[d][dd->ci[d]]*tot);
      }
    }
  }
}

static void distribute_cg(FILE *fplog,matrix box,t_block *cgs,rvec pos[],
			  gmx_domdec_t *dd)
{
  gmx_domdec_master_t *ma;
  int **tmp_ind=NULL,*tmp_nalloc=NULL;
  int  i,icg,j,k,k0,k1,d;
  rvec invbox,cg_cm;
  ivec ind;
  real nrcg,inv_ncg,pos_d;
  atom_id *cgindex;

  ma = dd->ma;

  /* Set the cell boundaries */
  set_dd_cell_sizes(dd,box,TRUE);

  if (tmp_ind == NULL) {
    snew(tmp_nalloc,dd->nnodes);
    snew(tmp_ind,dd->nnodes);
    for(i=0; i<dd->nnodes; i++) {
      tmp_nalloc[i] = (cgs->nr/(dd->nnodes*CG_ALLOC_SIZE) + 2)*CG_ALLOC_SIZE;
      snew(tmp_ind[i],tmp_nalloc[i]);
    }
  }

  /* Clear the count */
  for(i=0; i<dd->nnodes; i++) {
    ma->ncg[i] = 0;
    ma->nat[i] = 0;
  }

  for(d=0; (d<DIM); d++)
    invbox[d] = divide(1,box[d][d]);

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
      pos_d = cg_cm[d];
      if (dd->tric_dir[d] && dd->nc[d] > 1)
	for(j=d+1; j<DIM; j++)
	  pos_d -= cg_cm[j]*box[j][d]*invbox[j];
      while(pos_d >= box[d][d]) {
	pos_d -= box[d][d];
	rvec_dec(cg_cm,box[d]);
	for(k=k0; (k<k1); k++)
	  rvec_dec(pos[k],box[d]);
      }
      while(pos_d < 0) {
	pos_d += box[d][d];
	rvec_inc(cg_cm,box[d]);
	for(k=k0; (k<k1); k++)
	  rvec_inc(pos[k],box[d]);
      }
      /* This could be done more efficiently */
      ind[d] = 0;
      while(ind[d]+1 < dd->nc[d] && pos_d >= ma->cell_x[d][ind[d]+1])
	ind[d]++;
    }
    i = dd_index(dd->nc,ind);
    if (ma->ncg[i] == tmp_nalloc[i]) {
      tmp_nalloc[i] += CG_ALLOC_SIZE;
      srenew(tmp_ind[i],tmp_nalloc[i]);
    }
    tmp_ind[i][ma->ncg[i]] = icg;
    ma->ncg[i]++;
    ma->nat[i] += cgindex[icg+1] - cgindex[icg];
  }
  
  k1 = 0;
  for(i=0; i<dd->nnodes; i++) {
    ma->index[i] = k1;
    for(k=0; k<ma->ncg[i]; k++)
      ma->cg[k1++] = tmp_ind[i][k];
  }
  ma->index[dd->nnodes] = k1;

  for(i=0; i<dd->nnodes; i++)
    sfree(tmp_ind[i]);
  sfree(tmp_ind);
  sfree(tmp_nalloc);

  fprintf(fplog,"Charge group distribution:");
  for(i=0; i<dd->nnodes; i++)
    fprintf(fplog," %d",ma->ncg[i]);
  fprintf(fplog,"\n");
}

static void get_cg_distribution(FILE *fplog,gmx_domdec_t *dd,
				t_block *cgs,matrix box,rvec pos[])
{
  gmx_domdec_master_t *ma=NULL;
  int i,cg_gl;
  int *ibuf,buf2[2] = { 0, 0 };

  clear_dd_indices(dd,0);

  if (DDMASTER(dd)) {
    ma = dd->ma;
    if (ma->ncg == NULL) {
      snew(ma->ncg,dd->nnodes);
      snew(ma->index,dd->nnodes+1);
      snew(ma->cg,cgs->nr);
      snew(ma->nat,dd->nnodes);
      snew(ma->ibuf,dd->nnodes*2);
      snew(ma->cell_x,DIM);
      for(i=0; i<DIM; i++)
	snew(ma->cell_x[i],dd->nc[i]+1);
    }      
    
    distribute_cg(fplog,box,cgs,pos,dd);
    for(i=0; i<dd->nnodes; i++) {
      ma->ibuf[2*i]   = ma->ncg[i];
      ma->ibuf[2*i+1] = ma->nat[i];
    }
    ibuf = ma->ibuf;
  } else {
    ibuf = NULL;
  }
  dd_scatter(dd,2*sizeof(int),ibuf,buf2);

  dd->ncg_home = buf2[0];
  dd->nat_home = buf2[1];
  if (dd->ncg_home > dd->cg_nalloc || dd->cg_nalloc == 0) {
    dd->cg_nalloc = over_alloc(dd->ncg_home);
    srenew(dd->index_gl,dd->cg_nalloc);
    srenew(dd->cgindex,dd->cg_nalloc+1);
  }
  if (DDMASTER(dd)) {
    for(i=0; i<dd->nnodes; i++) {
      ma->ibuf[i] = ma->ncg[i]*sizeof(int);
      ma->ibuf[dd->nnodes+i] = ma->index[i]*sizeof(int);
    }
  }

  dd_scatterv(dd,
	      DDMASTER(dd) ? ma->ibuf : NULL,
	      DDMASTER(dd) ? ma->ibuf+dd->nnodes : NULL,
	      DDMASTER(dd) ? ma->cg : NULL,
	      dd->ncg_home*sizeof(int),dd->index_gl);

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

static int compact_and_copy_vec_at(int ncg,int *move,
				   int *cgindex,
				   int nvec,int vec,
				   rvec *src,gmx_domdec_comm_t *comm)
{
  int m,icg,i,i0,i1,nrcg;
  int home_pos;
  int pos_vec[DIM*2];

  home_pos = 0;

  for(m=0; m<DIM*2; m++)
    pos_vec[m] = 0;

  i0 = 0;
  for(icg=0; icg<ncg; icg++) {
    i1 = cgindex[icg+1];
    m = move[icg];
    if (m == -1) {
      /* Compact the home array in place */
      for(i=i0; i<i1; i++)
	copy_rvec(src[i],src[home_pos++]);
    } else {
      /* Copy to the communication buffer */
      nrcg = i1 - i0;
      pos_vec[m] += 1 + vec*nrcg;
      for(i=i0; i<i1; i++)
	copy_rvec(src[i],comm->cgcm_state[m][pos_vec[m]++]);
      pos_vec[m] += (nvec - vec - 1)*nrcg;
    }
    i0 = i1;
  }

  return home_pos;
}

static int compact_and_copy_vec_cg(int ncg,int *move,
				   int *cgindex,
				   int nvec,rvec *src,gmx_domdec_comm_t *comm)
{
  int m,icg,i0,i1,nrcg;
  int home_pos;
  int pos_vec[DIM*2];

  home_pos = 0;

  for(m=0; m<DIM*2; m++)
    pos_vec[m] = 0;

  i0 = 0;
  for(icg=0; icg<ncg; icg++) {
    i1 = cgindex[icg+1];
    m = move[icg];
    if (m == -1) {
      /* Compact the home array in place */
      copy_rvec(src[icg],src[home_pos++]);
    } else {
      nrcg = i1 - i0;
      /* Copy to the communication buffer */
      copy_rvec(src[icg],comm->cgcm_state[m][pos_vec[m]]);
      pos_vec[m] += 1 + nrcg*nvec;
    }
    i0 = i1;
  }

  return home_pos;
}

static int compact_ind(int ncg,int *move,
		       int *index_gl,int *cgindex,
		       int *gatindex,gmx_ga2la_t *ga2la,
		       int *solvent_type)
{
  int cg,nat,a0,a1,a,a_gl;
  int home_pos;
  
  home_pos = 0;
  nat = 0;
  for(cg=0; cg<ncg; cg++) {
    a0 = cgindex[cg];
    a1 = cgindex[cg+1];
    if (move[cg] == -1) {
      /* Compact the home arrays in place.
       * Anything that can be done here avoids access to global arrays.
       */
      cgindex[home_pos] = nat;
      for(a=a0; a<a1; a++) {
	a_gl = gatindex[a];
	gatindex[nat] = a_gl;
	/* The cell number stays 0, so we don't need to set it */
	ga2la[a_gl].a = nat;
	nat++;
      }
      index_gl[home_pos] = index_gl[cg];
      if (solvent_type)
	solvent_type[home_pos] = solvent_type[cg];
      home_pos++;
    } else {
      /* Clear the global indices */
      for(a=a0; a<a1; a++) {
	a_gl = gatindex[a];
	ga2la[a_gl].cell = -1;
      }
    }
  }
  cgindex[home_pos] = nat;
  
  return home_pos;
}

static int dd_redistribute_cg(FILE *fplog,
			      gmx_domdec_t *dd,t_block *gcgs,
			      t_state *state,rvec cg_cm[],
			      t_forcerec *fr,t_mdatoms *md,
			      t_nrnb *nrnb)
{
  int  *move;
  int  ncg[DIM*2],nat[DIM*2];
  int  c,i,cg,k,k0,k1,d,dim,dim2,dir,d2,d3,d4,cell_d;
  int  mc,cdd,nrcg,ncg_recv,nat_recv,nvs,nvr,nvec,vec;
  int  sbuf[2],rbuf[2];
  int  home_pos_cg,home_pos_at,ncg_stay_home,buf_pos;
  int  flag;
  ivec tric_dir,dev;
  real inv_ncg,pos_d;
  rvec invbox,cell_x0,cell_x1,limit0,limit1;
  atom_id *cgindex;
  gmx_domdec_comm_t *comm;

  comm = dd->comm;

  if (dd->ncg_tot > comm->nalloc_int) {
    comm->nalloc_int = over_alloc(dd->ncg_tot);
    srenew(comm->buf_int,comm->nalloc_int);
  }
  move = comm->buf_int;

  /* Clear the count */
  for(c=0; c<dd->ndim*2; c++) {
    ncg[c] = 0;
    nat[c] = 0;
  }

  for(d=0; (d<DIM); d++) {
    invbox[d] = divide(1,state->box[d][d]);
    cell_x0[d] = dd->cell_x0[d];
    cell_x1[d] = dd->cell_x1[d];
    c = dd->ci[d] - 1;
    if (c < 0)
      c = dd->nc[d] - 1;
    limit0[d] = cell_x0[d] - comm->distance;
    c = dd->ci[d] + 1;
    if (c >= dd->nc[d])
      c = 0;
    limit1[d] = cell_x1[d] + comm->distance;
    if (dd->tric_dir[d] && dd->nc[d] > 1)
      tric_dir[d] = 1;
    else
      tric_dir[d] = 0;
  }

  cgindex = dd->cgindex;

  /* Compute the center of geometry for all home charge groups
   * and put them in the box and determine where they should go.
   */
  for(cg=0; cg<dd->ncg_home; cg++) {
    k0   = cgindex[cg];
    k1   = cgindex[cg+1];
    nrcg = k1 - k0;
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
      if (dd->nc[d] > 1) {
	pos_d = cg_cm[cg][d];
	if (tric_dir[d])
	  for(d2=d+1; d2<DIM; d2++)
	    pos_d -= cg_cm[cg][d2]*state->box[d2][d]*invbox[d2];
	/* Put the charge group in the triclinic unit-cell */
	if (pos_d >= cell_x1[d]) {
	  if (pos_d >= limit1[d])
	    gmx_fatal(FARGS,"The charge group starting at atom %d moved more than the cut-off (%f) in direction %c: distance out of cell %f, coords %f %f %f",
		      glatnr(dd,dd->cgindex[cg]),
		      comm->distance,dim2char(d),pos_d - cell_x1[d],
		      cg_cm[cg][XX],cg_cm[cg][YY],cg_cm[cg][ZZ]);
	  dev[d] = 1;
	  if (dd->ci[d] == dd->nc[d] - 1) {
	    rvec_dec(cg_cm[cg],state->box[d]);
	    for(k=k0; (k<k1); k++)
	      rvec_dec(state->x[k],state->box[d]);
	  }
	} else if (pos_d < cell_x0[d]) {
	  if (pos_d < limit0[d])
	    gmx_fatal(FARGS,"The charge group starting at atom %d moved more than the cut-off (%f) in direction %c: distance out of cell %f, coords %f %f %f",
		      glatnr(dd,dd->cgindex[cg]),
		      comm->distance,dim2char(d),pos_d - cell_x0[d],
		      cg_cm[cg][XX],cg_cm[cg][YY],cg_cm[cg][ZZ]);
	  dev[d] = -1;
	  if (dd->ci[d] == 0) {
	    rvec_inc(cg_cm[cg],state->box[d]);
	    for(k=k0; (k<k1); k++)
	      rvec_inc(state->x[k],state->box[d]);
	  }
	} else {
	  dev[d] = 0;
	}
      } else {
	/* Put the charge group in the rectangular unit-cell */
	while (cg_cm[cg][d] >= state->box[d][d]) {
	  rvec_dec(cg_cm[cg],state->box[d]);
	  for(k=k0; (k<k1); k++)
	    rvec_dec(state->x[k],state->box[d]);
	}
	while (cg_cm[cg][d] < 0) {
	  rvec_inc(cg_cm[cg],state->box[d]);
	  for(k=k0; (k<k1); k++)
	    rvec_inc(state->x[k],state->box[d]);
	}
      }
    }

    /* Determine where this cg should go */
    flag = 0;
    mc = -1;
    for(d=0; d<dd->ndim; d++) {
      dim = dd->dim[d];
      if (dev[dim] < -1 || dev[dim] > 1)
	gmx_fatal(FARGS,"The charge group starting at atom %d moved more than one cell: %d %d %d, coords %f %f %f",
		  glatnr(dd,dd->cgindex[cg]),
		  dev[XX],dev[YY],dev[ZZ],
		  cg_cm[cg][XX],cg_cm[cg][YY],cg_cm[cg][ZZ]);
      if (dev[dim] != 0) {
	if (dd->nc[dim] == 2 || dev[dim] == 1) {
	  flag |= 1<<(16+d*2);
	  if (mc == -1)
	    mc = d*2;
	} else {
	  flag |= 1<<(16+d*2+1);
	  if (mc == -1)
	    mc = d*2 + 1;
	}
      }
    }
    move[cg] = mc;
    if (mc >= 0) {
      if (ncg[mc]+1 > comm->cggl_flag_nalloc[mc]) {
	comm->cggl_flag_nalloc[mc] = over_alloc(ncg[mc]+1);
	srenew(comm->cggl_flag[mc],comm->cggl_flag_nalloc[mc]*DD_CGIBS);
      }
      comm->cggl_flag[mc][ncg[mc]*DD_CGIBS  ] = dd->index_gl[cg];
      /* We store the cg size in the lower 16 bits
       * and the place where the charge group should go
       * in the next 6 bits. This saves some communication volume.
       */
      comm->cggl_flag[mc][ncg[mc]*DD_CGIBS+1] = nrcg | flag;
      ncg[mc] += 1;
      nat[mc] += nrcg;
    }
  }

  inc_nrnb(nrnb,eNR_CGCM,dd->nat_home);
  inc_nrnb(nrnb,eNR_RESETX,dd->ncg_home);

  nvec = 1;
  if (state->v)
    nvec++;
  if (state->sd_X)
    nvec++;

  /* Make sure the communication buffers are large enough */
  for(mc=0; mc<dd->ndim*2; mc++) {
    nvr = ncg[mc] + nat[mc]*nvec;
    if (nvr>comm->cgcm_state_nalloc[mc]) {
      comm->cgcm_state_nalloc[mc] = over_alloc(nvr);
      srenew(comm->cgcm_state[mc],comm->cgcm_state_nalloc[mc]);
    }
  }

  /* Recalculating cg_cm might be cheaper than communicating,
   * but that could give rise to rounding issues.
   */
  home_pos_cg =
    compact_and_copy_vec_cg(dd->ncg_home,move,cgindex,
			    nvec,cg_cm,comm);

  vec = 0;
  home_pos_at =
    compact_and_copy_vec_at(dd->ncg_home,move,cgindex,
			    nvec,vec++,state->x,comm);
  if (state->v)
    compact_and_copy_vec_at(dd->ncg_home,move,cgindex,
			    nvec,vec++,state->v,comm);
  if (state->sd_X)
    compact_and_copy_vec_at(dd->ncg_home,move,cgindex,
			    nvec,vec++,state->sd_X,comm);
  
  compact_ind(dd->ncg_home,move,
	      dd->index_gl,dd->cgindex,
	      dd->gatindex,dd->ga2la,
	      fr->solvent_opt==esolNO ? NULL : fr->solvent_type);

  ncg_stay_home = home_pos_cg;
  for(d=0; d<dd->ndim; d++) {
    dim = dd->dim[d];
    ncg_recv = 0;
    nat_recv = 0;
    nvr      = 0;
    for(dir=0; dir<(dd->nc[dim]==2 ? 1 : 2); dir++) {
      cdd = d*2 + dir;
      /* Communicate the cg and atom counts */
      sbuf[0] = ncg[cdd];
      sbuf[1] = nat[cdd];
      if (debug)
	fprintf(debug,"Sending ddim %d dir %d: ncg %d nat %d\n",
		d,dir,sbuf[0],sbuf[1]);
      dd_sendrecv_int(dd,d,dir,sbuf,2,rbuf,2);

      if ((ncg_recv+rbuf[0])*DD_CGIBS > comm->nalloc_int) {
	comm->nalloc_int = over_alloc((ncg_recv+rbuf[0])*DD_CGIBS);
	srenew(comm->buf_int,comm->nalloc_int);
      }

      /* Communicate the charge group indices, sizes and flags */
      dd_sendrecv_int(dd,d,dir,
		      comm->cggl_flag[cdd],sbuf[0]*DD_CGIBS,
		      comm->buf_int+ncg_recv*DD_CGIBS,rbuf[0]*DD_CGIBS);
      
      nvs = ncg[cdd] + nat[cdd]*nvec;
      i   = rbuf[0]  + rbuf[1] *nvec;
      if (nvr+i > comm->nalloc_vr) {
	comm->nalloc_vr = over_alloc(nvr+i);
	srenew(comm->buf_vr,comm->nalloc_vr);
      }

      /* Communicate cgcm and state */
      dd_sendrecv_rvec(dd,d,dir,
		       comm->cgcm_state[cdd],nvs,
		       comm->buf_vr+nvr,i);
      ncg_recv += rbuf[0];
      nat_recv += rbuf[1];
      nvr      += i;
    }

    /* Process the received charge groups */
    buf_pos = 0;
    for(cg=0; cg<ncg_recv; cg++) {
      flag = comm->buf_int[cg*DD_CGIBS+1];
      mc = -1;
      if (d < dd->ndim-1) {
	if (dd->bDynLoadBal) {
	  /* Clear the move flags */
	  flag &= ~(63<<16);
	  /* Check where this cg should go */
	  for(d2=d+1; (d2<dd->ndim && mc==-1); d2++) {
	    dim2 = dd->dim[d2];
	    pos_d = comm->buf_vr[buf_pos][dim2];
	    if (tric_dir[dim2])
	      for(d3=dim2+1; d3<DIM; d3++)
		pos_d -= comm->buf_vr[buf_pos][d3]*state->box[d3][dim2]*invbox[d3];
	    if (pos_d >= cell_x1[dim2]) {
	      flag |= (1<<(16+d2*2));
	      mc = d2*2;
	    } else if (pos_d < cell_x0[dim2]) {
	      if (dd->nc[dim2] == 2) {
		flag |= (1<<(16+d2*2));
		mc = d2*2;
	      } else {
		flag |= (1<<(16+d2*2+1));
		mc = d2*2+1;
	      }
	    }
	  }
	  comm->buf_int[cg*DD_CGIBS+1] = flag;
	} else {
	  /* The communicated flag tells where this cg should go */
	  for(d2=d+1; (d2<dd->ndim && mc==-1); d2++) {
	    if (flag & (1<<(16+d2*2))) {
	      mc = d2*2;
	    } else if (flag & (1<<(16+d2*2+1))) {
	      mc = d2*2+1;
	    }
	  }
	}
      }
      
      nrcg = flag & 65535;
      if (mc == -1) {
	if (home_pos_cg+1 > dd->cg_nalloc) {
	  dd->cg_nalloc = over_alloc(home_pos_cg+1);
	  srenew(dd->index_gl,dd->cg_nalloc);
	  srenew(dd->cgindex,dd->cg_nalloc+1);
	}
	/* If the other home arrays were dynamic, we would need
	 * to reallocate here.
	 */
	/* Set the global charge group index and size */
	dd->index_gl[home_pos_cg] = comm->buf_int[cg*DD_CGIBS];
	dd->cgindex[home_pos_cg+1] = dd->cgindex[home_pos_cg] + nrcg;
	/* Copy the state from the buffer */
	copy_rvec(comm->buf_vr[buf_pos++],cg_cm[home_pos_cg]);
	for(i=0; i<nrcg; i++)
	  copy_rvec(comm->buf_vr[buf_pos++],state->x[home_pos_at+i]);
	if (state->v)
	  for(i=0; i<nrcg; i++)
	    copy_rvec(comm->buf_vr[buf_pos++],state->v[home_pos_at+i]);
	if (state->sd_X)
	  for(i=0; i<nrcg; i++)
	    copy_rvec(comm->buf_vr[buf_pos++],state->sd_X[home_pos_at+i]);
	home_pos_cg += 1;
	home_pos_at += nrcg;
      } else {
	/* Reallocate the buffers if necessary  */
	if (ncg[mc]+1 > comm->cggl_flag_nalloc[mc]) {
	  comm->cggl_flag_nalloc[mc] = over_alloc(ncg[mc]+1);
	  srenew(comm->cggl_flag[mc],comm->cggl_flag_nalloc[mc]*DD_CGIBS);
	}
	nvr = ncg[mc] + nat[mc]*nvec;
	if (nvr + 1 + nrcg*nvec > comm->cgcm_state_nalloc[mc]) {
	  comm->cgcm_state_nalloc[mc] = over_alloc(nvr + 1 + nrcg*nvec);
	  srenew(comm->cgcm_state[mc],comm->cgcm_state_nalloc[mc]);
	}
	/* Copy from the receive to the send buffers */
	memcpy(comm->cggl_flag[mc] + ncg[mc]*DD_CGIBS,
	       comm->buf_int + cg*DD_CGIBS,
	       DD_CGIBS*sizeof(int));
	memcpy(comm->cgcm_state[mc][nvr],
	       comm->buf_vr[buf_pos],
	       (1+nrcg*nvec)*sizeof(rvec));
	buf_pos += 1 + nrcg*nvec;
	ncg[mc] += 1;
	nat[mc] += nrcg;
      }
    }
  }

  /* Clear the local indices, except for the home cell.
   * The home cell indices were updated and cleaned in compact_ind.
   */
  clear_dd_indices(dd,dd->nat_home);

  dd->ncg_home = home_pos_cg;
  dd->nat_home = home_pos_at;

  dd->bMasterHasAllCG = FALSE;

  if (debug)
    fprintf(debug,"Finished repartitioning\n");

  return ncg_stay_home;
}

void setup_dd_grid(FILE *fplog,matrix box,gmx_domdec_t *dd)
{
  bool bZYX;
  int  start,inc;
  int  d,i,j,m;
  ivec tmp,s;
  int  ncell,ncellp;
  ivec dd_cp[DD_MAXICELL];
  gmx_domdec_ns_ranges_t *icell;

  gmx_ddindex2xyz(dd->nc,dd->nodeid,dd->ci);

  if (getenv("GMX_DD_ORDER_ZYX")) {
    /* Decomposition order z,y,x */
    fprintf(fplog,"Using domain decomposition order z, y, x\n");
    start = DIM-1;
    inc   = -1;
  } else {
    /* Decomposition order x,y,z */
    start = 0;
    inc   = 1;
  }

  /* The x/y/z communication order is set by this loop */
  dd->ndim = 0;
  for(d=start; (d>=0 && d<DIM); d+=inc) {
    if (dd->nc[d] > 1) {
      dd->dim[dd->ndim] = d;
      copy_ivec(dd->ci,tmp);
      tmp[d] = (tmp[d] + 1) % dd->nc[d];
      dd->neighbor[dd->ndim][0] = dd_index(dd->nc,tmp);
      copy_ivec(dd->ci,tmp);
      tmp[d] = (tmp[d] - 1 + dd->nc[d]) % dd->nc[d];
      dd->neighbor[dd->ndim][1] = dd_index(dd->nc,tmp);
      dd->ndim++;
    }
  }
  
  fprintf(fplog,"Making %dD domain decomposition %d x %d x %d, home cell index %d %d %d\n",
	  dd->ndim,
	  dd->nc[XX],dd->nc[YY],dd->nc[ZZ],
	  dd->ci[XX],dd->ci[YY],dd->ci[ZZ]);
  if (DDMASTER(dd))
    fprintf(stderr,"Making %dD domain decomposition %d x %d x %d\n",
	    dd->ndim,dd->nc[XX],dd->nc[YY],dd->nc[ZZ]);
  switch (dd->ndim) {
  case 3:
    ncell  = dd_c3n;
    ncellp = dd_cp3n;
    for(i=0; i<ncellp; i++)
      copy_ivec(dd_cp3[i],dd_cp[i]);
    break;
  case 2:
    ncell  = dd_c2n;
    ncellp = dd_cp2n;
    for(i=0; i<ncellp; i++)
      copy_ivec(dd_cp2[i],dd_cp[i]);
    break;
  case 1:
    ncell  = dd_c1n;
    ncellp = dd_cp1n;
    for(i=0; i<ncellp; i++)
      copy_ivec(dd_cp1[i],dd_cp[i]);
    break;
  default:
    gmx_fatal(FARGS,"Can only do 1, 2 or 3D domain decomposition");
    ncell = 0;
    ncellp = 0;
  }
    
  for(i=0; i<ncell; i++) {
    m = 0;
    for(d=start; (d>=0 && d<DIM); d+=inc) {
      if (dd->nc[d] > 1)
	dd->shift[i][d] = dd_co[i][m++];
      else 
	dd->shift[i][d] = 0;
    }
  }

  dd->ncell  = ncell;
  for(i=0; i<ncell; i++) {
    for(d=0; d<DIM; d++) {
      s[d] = dd->ci[d] - dd->shift[i][d];
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
	/* All shifts should be allowed */
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
  int  nodeid=-1;

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
  int  pmenode=-1;
 
  dd = cr->dd;

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

  if (DDMASTER(dd) && dd->ma == NULL) {
    snew(dd->ma,1);
  }
}

static real *get_cell_load(FILE *fplog,char *dir,int nc,char *load_string)
{
  real *cell_load;
  int  i,n;
  double dbl;

  cell_load = NULL;
  if (nc > 1 && load_string != NULL) {
    fprintf(fplog,"Using static load balancing for the %s direction\n",dir);
    snew(cell_load,nc);
    for (i=0; i<nc; i++) {
      dbl = 0;
      sscanf(load_string,"%lf%n",&dbl,&n);
      if (dbl == 0)
	gmx_fatal(FARGS,
		  "Incorrect or not enough load entries for direction %s",
		  dir);
      cell_load[i] = dbl;
      load_string += n;
    }
  }
  
  return cell_load;
}

gmx_domdec_t *init_domain_decomposition(FILE *fplog,t_commrec *cr,ivec nc,
					bool bDynLoadBal,
					char *loadx,char *loady,char *loadz)
{
  gmx_domdec_t *dd;
  int  d;
  
  fprintf(fplog,
	  "Domain decomposition grid %d x %d x %d, separate PME nodes %d\n",
	  nc[XX],nc[YY],nc[ZZ],cr->npmenodes);

  snew(dd,1);
  snew(dd->comm,1);
  snew(dd->comm->cggl_flag,DIM*2);
  snew(dd->comm->cgcm_state,DIM*2);

  copy_ivec(nc,dd->nc);
  dd->nnodes = dd->nc[XX]*dd->nc[YY]*dd->nc[ZZ];
  if (dd->nnodes != cr->nnodes - cr->npmenodes)
    gmx_fatal(FARGS,"The size of the domain decomposition grid (%d) does not match the number of nodes (%d)\n",dd->nnodes,cr->nnodes - cr->npmenodes);

  dd->bDynLoadBal = bDynLoadBal;
  snew(dd->cell_load,DIM);
  if (!dd->bDynLoadBal) {
    dd->cell_load[XX] = get_cell_load(fplog,"x",dd->nc[XX],loadx);
    dd->cell_load[YY] = get_cell_load(fplog,"y",dd->nc[YY],loady);
    dd->cell_load[ZZ] = get_cell_load(fplog,"z",dd->nc[ZZ],loadz);
  }

  bDDDump = (getenv("GMX_DD_DUMP") != NULL);

  return dd;
}

static void setup_dd_communication(FILE *fplog,int step,
				   gmx_domdec_t *dd,t_block *gcgs,
				   rvec buf[],matrix box,rvec cg_cm[],
				   real comm_distance)
{
  int dim_ind,dim,nat_tot,ncell,cell,celli,c,i,cg,cg_gl,nrcg,d,d1,d2;
  int *ncg_cell,*index_gl,*cgindex;
  gmx_domdec_comm_t *comm;
  gmx_domdec_ind_t *ind;
  rvec corner,v[DIM][DIM],skew_fac;
  real corner0=0,corner1=0,r_comm2,r,r2,inv_ncg,dep;
  int  nsend,nat;
  bool bTric;

  if (debug)
    fprintf(debug,"Setting up DD communication\n");

  comm = dd->comm;
  
  comm->distance = comm_distance;

  bTric = TRICLINIC(box);

  /* This should be changed for non-uniform grids */
  for(d=0; d<DIM; d++)
    corner[d] = dd->cell_x0[d];
  if (dd->ndim >= 2) {
    /* Set the upper-right corner for rounding */
    d = dd->dim[0];
    corner0 = dd->cell_x1[d];
    d = dd->dim[1];
    corner1 = dd->cell_x1[d];
  }

  r_comm2 = sqr(comm->distance);

  for(dim=0; dim<dd->ndim; dim++) {
    d = dd->dim[dim];
    /* Determine the cell thickness in this direction */
    r2 = sqr(box[d][d]);
    if (bTric) {
      if (d == XX || d == YY) {
	/* Normalize such that the "diagonal" is 1 */
	svmul(1/box[d+1][d+1],box[d+1],v[d][d+1]);
	for(i=0; i<d; i++)
	  v[d][d+1][i] = 0;
	r2 -= sqr(iprod(box[d],v[d][d+1]));
	if (d == XX) {
	  /* Normalize such that the "diagonal" is 1 */
	  svmul(1/box[d+2][d+2],box[d+2],v[d][d+2]);
	  for(i=0; i<d; i++)
	    v[d][d+2][i] = 0;
	  /* Make vector [d+2] perpendicular to vector [d+1],
	   * this does not affect the normalization.
	   */
	  dep = iprod(v[d][d+1],v[d][d+2])/norm2(v[d][d+1]);
	  for(i=0; i<DIM; i++)
	    v[d][d+2][i] -= dep*v[d][d+1][i];
	  r2 -= sqr(iprod(box[d],v[d][d+2]));
	}
      }
      if (debug) {
	fprintf(debug,"box[%d]  %.3f %.3f %.3f",
		d,box[d][XX],box[d][YY],box[d][ZZ]);
	for(i=d+1; i<DIM; i++)
	  fprintf(debug,"  v[%d] %.3f %.3f %.3f",
		  i,v[d][i][XX],v[d][i][YY],v[d][i][ZZ]);
	fprintf(debug,"\n");
      }
    }

    r2 *= sqr((dd->cell_x1[d] - dd->cell_x0[d])/box[d][d]);

    if (r2 < r_comm2)
      gmx_fatal(FARGS,"Step %d: One of the domain decomposition grid cell sizes (%f %f %f: %c size %f) is smaller than the cut-off (%f)",
		step,
		dd->cell_x1[XX] - dd->cell_x0[XX],
		dd->cell_x1[YY] - dd->cell_x0[YY],
		dd->cell_x1[ZZ] - dd->cell_x0[ZZ],
		dim2char(d),sqrt(r2),comm->distance);
  }

  if (bTric) {
    /* Determine correction factors for transforming
     * x or y distances^2 to distances^2 between triclinic cell planes.
     */
    for(d=0; d<DIM; d++) {
      skew_fac[d] = 1;
      for(i=d+1; i<DIM; i++)
	skew_fac[d] -= sqr(v[d][i][d]);
    }
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
    ind = &comm->ind[dim_ind];
    nsend = 0;
    nat = 0;
    for(cell=0; cell<ncell; cell++) {
      ind->nsend[cell] = 0;
      celli = cell_perm[dim_ind][cell];
      for(cg=ncg_cell[celli]; cg<ncg_cell[celli+1]; cg++) {
	if (!bTric) {
	  /* Rectangular box, easy */
	  r2 = sqr(cg_cm[cg][dim] - corner[dim]);
	  /* Rounding gives at most a 16% reduction in communicated atoms */
	  if (dim_ind >= 1 && (celli == 1 || celli == 2))
	    r2 += sqr(cg_cm[cg][dd->dim[0]] - corner0);
	  if (dim_ind == 2 && (celli == 2 || celli == 3))
	    r2 += sqr(cg_cm[cg][dd->dim[1]] - corner1);
	} else {
	  /* Triclinic box, complicated */
	  r = cg_cm[cg][dim] - corner[dim];
	  for(i=dim+1; i<DIM; i++)
	    r -= cg_cm[cg][i]*v[dim][i][dim];
	  r2 = sqr(r)*skew_fac[dim];
	  /* Rounding, conservative as the skew_fac multiplication
	   * will slighlty underestimate the distance.
	   */
	  if (dim_ind >= 1 && (celli == 1 || celli == 2)) {
	    r = cg_cm[cg][dd->dim[0]] - corner0;
	    for(i=dd->dim[0]+1; i<DIM; i++)
	      r -= cg_cm[cg][i]*v[dd->dim[0]][i][dd->dim[0]];
	    r2 += sqr(r)*skew_fac[dd->dim[0]];
	  }
	  if (dim_ind == 2 && (celli == 2 || celli == 3)) {
	    r = cg_cm[cg][dd->dim[1]] - corner1;
	    for(i=dd->dim[1]+1; i<DIM; i++)
	      r -= cg_cm[cg][i]*v[dd->dim[1]][i][dd->dim[1]];
	    r2 += sqr(r)*skew_fac[dd->dim[1]];
	  }
	}
	if (r2 < r_comm2) {
	  /* Make an index to the local charge groups */
	  if (nsend >= ind->nalloc) {
	    ind->nalloc += CG_ALLOC_SIZE;
	    srenew(ind->index,ind->nalloc);
	  }
	  if (nsend >= comm->nalloc_int) {
	    comm->nalloc_int += CG_ALLOC_SIZE;
	    srenew(comm->buf_int,comm->nalloc_int);
	  }
	  ind->index[nsend] = cg;
	  comm->buf_int[nsend] = index_gl[cg];
	  ind->nsend[cell]++;
	  if (dd->ci[dim] == 0) {
	    /* Correct cg_cm for pbc */
	    rvec_add(cg_cm[cg],box[dim],buf[nsend]);
	  } else {
	    copy_rvec(cg_cm[cg],buf[nsend]);
	  }
	  nsend++;
	  nat += cgindex[cg+1] - cgindex[cg];
	}
      }
    }
    ind->nsend[ncell]   = nsend;
    ind->nsend[ncell+1] = nat;
    /* Communicate the number of cg's and atoms to receive */
    dd_sendrecv_int(dd, dim_ind, ddBackward,
		    ind->nsend, ncell+2,
		    ind->nrecv, ncell+2);
    /* Communicate the global cg indices, receive in place */
    if (ncg_cell[ncell] + ind->nrecv[ncell] > dd->cg_nalloc
	|| dd->cg_nalloc == 0) {
      dd->cg_nalloc = over_alloc(ncg_cell[ncell] + ind->nrecv[ncell]);
      srenew(index_gl,dd->cg_nalloc);
      srenew(cgindex,dd->cg_nalloc+1);
    }
    dd_sendrecv_int(dd, dim_ind, ddBackward,
		    comm->buf_int,            nsend,
		    index_gl+ncg_cell[ncell], ind->nrecv[ncell]);
    /* Communicate cg_cm, receive in place */
    dd_sendrecv_rvec(dd, dim_ind, ddBackward,
		     buf,                   nsend,
		     cg_cm+ncg_cell[ncell], ind->nrecv[ncell]);
    /* Make the charge group index */
    for(cell=ncell; cell<2*ncell; cell++) {
      ncg_cell[cell+1] = ncg_cell[cell] + ind->nrecv[cell-ncell];
      for(cg=ncg_cell[cell]; cg<ncg_cell[cell+1]; cg++) {
	cg_gl = index_gl[cg];
	nrcg = gcgs->index[cg_gl+1] - gcgs->index[cg_gl];
	cgindex[cg+1] = cgindex[cg] + nrcg;
	nat_tot += nrcg;
      }
    }
    ncell += ncell;
  }
  dd->index_gl = index_gl;
  dd->cgindex  = cgindex;

  dd->ncg_tot = ncg_cell[dd->ncell];
  dd->nat_tot       = nat_tot;
  dd->nat_tot_vsite = nat_tot;
  dd->nat_tot_con   = nat_tot;

  if (debug) {
    fprintf(debug,"Finished setting up DD communication, cells:");
    for(c=0; c<dd->ncell; c++)
      fprintf(debug," %d",dd->ncg_cell[c+1]-dd->ncg_cell[c]);
    fprintf(debug,"\n");
  }
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

void dd_partition_system(FILE         *fplog,
			 int          step,
			 t_commrec    *cr,
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
  gmx_domdec_t *dd;
  int i,j,cg0;

  dd = cr->dd;

  set_dd_cell_sizes(dd,bMasterState ? state_global->box : state_local->box,
		    FALSE);

  if (ir->ePBC == epbcNONE)
    gmx_fatal(FARGS,"pbc type %s is not supported with domain decomposition",
	      epbc_names[epbcNONE]);

  if (dd->ndim == DIM)
    fr->ePBC = epbcXYZ;
  else
    fr->ePBC = epbcFULL;

  if (ir->ns_type == ensSIMPLE)
    gmx_fatal(FARGS,"ns type %s is not supported with domain decomposition",
	      ens_names[ensSIMPLE]);
  
  if (ir->eConstrAlg == estSHAKE)
    gmx_fatal(FARGS,
	      "%s is not supported (yet) with domain decomposition, use %s",
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
  setup_dd_communication(fplog,step,dd,&top_global->blocks[ebCGS],
			 buf,state_local->box,
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
  make_local_top(fplog,dd,fr,top_global,top_local);

  dd_update_ns_border(dd,nsb);

  if (top_global->idef.il[F_CONSTR].nr > 0) {
    make_local_constraints(dd,top_global->idef.il[F_CONSTR].iatoms,
			   ir->nProjOrder);
  } else {
    dd->nat_tot_con = dd->nat_tot_vsite;
  }
  /* Make space for the extra coordinates for virtual site
   * or constraint communication.
   */
  /* This is not a nice solution.
   * state->natoms is always equal to the global number of atoms.
   * Reallocation will only happen for very small systems
   * with cell sizes close to the cut-off.
   */
  if (dd->nat_tot_con > state_local->natoms) {
    state_local->natoms = dd->nat_tot_con;
    srenew(state_local->x,state_local->natoms);
    srenew(state_local->v,state_local->natoms);
  }

  /* We make the all mdatoms up to nat_tot_con.
   * We could save some work by only setting invmass
   * between nat_tot and nat_tot_con.
   */
  atoms2md(&top_global->atoms,ir,top_global->idef.il[F_ORIRES].nr,
	   dd->nat_tot_con,dd->gatindex,mdatoms);

  if (!(cr->duty & DUTY_PME))
    /* Send the charges to our PME only node */
    gmx_pme_send_x_q(cr,state_local->box,NULL,
		     mdatoms->chargeA,mdatoms->chargeB,
		     mdatoms->nChargePerturbed,0,FALSE);

  if (dd->constraints || top_global->idef.il[F_SETTLE].nr>0)
    init_constraints(fplog,top_global,&top_local->idef.il[F_SETTLE],ir,mdatoms,
		     START(nsb),HOMENR(nsb),
		     ir->eI!=eiSteep,NULL,dd);

  /* We need the constructing atom coordinates of the virtual sites
   * when spreading the forces.
   */
  dd_move_x_vsites(dd,state_local->box,state_local->x);

  if (bDDDump) {
    dd_move_x(dd,state_local->box,state_local->x,buf);
    write_dd_pdb("dd_dump",step,"dump",&top_global->atoms,dd,dd->nat_tot,
		 state_local->x,state_local->box);
  }
}
