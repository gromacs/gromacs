#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <stdio.h>
#include <time.h>
#include <math.h>
#include "typedefs.h"
#include "smalloc.h"
#include "vec.h"
#include "domdec.h"

#ifdef GMX_MPI
#include <mpi.h>
#endif


#define dd_c3n 8
static const ivec dd_c3[dd_c3n] =
  {{0,0,0},{1,0,0},{1,1,0},{0,1,0},{0,1,1},{0,0,1},{1,0,1},{1,1,1}};
#define dd_cp3n 4
static const ivec dd_cp3[dd_cp3n] =
  {{0,0,8},{1,3,6},{2,5,7},{3,6,7}};

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

#define dd_index(n,i) ((((i)[ZZ]*(n)[YY] + (i)[YY])*n[XX]) + (i)[XX])

static void index2xyz(ivec nc,int ind,ivec xyz)
{
  xyz[XX] = ind % nc[XX];
  xyz[YY] = (ind / nc[XX]) % nc[YY];
  xyz[ZZ] = ind / (nc[YY]*nc[XX]);
}

int dd_nicg(gmx_domdec_t *dd)
{
  return dd->icellcg1[dd->nicell-1];
}

int dd_ncg_tot(gmx_domdec_t *dd)
{
  return dd->ncg_tot;
}

int dd_cg_j0(gmx_domdec_t *dd,int icg)
{
  int icell,jcg0=-1;

  icell = 0;
  while (icg >= dd->icellcg1[icell])
    icell++;

  if (icell == 0)
    jcg0 = icg;
  else if (icell < dd->nicell)
    jcg0 = dd->icelljcg0[icell];
  else
    gmx_fatal(FARGS,"DD icg %d out of range: icell (%d) >= nicell (%d)",
	      icg,icell,dd->nicell);

  return jcg0;
}

int dd_cg_j1(gmx_domdec_t *dd,int icg)
{
  int icell;

  icell = 0;
  while (icg >= dd->icellcg1[icell])
    icell++;

  return dd->icelljcg1[icell];
}

void dd_move_x(gmx_domdec_t *dd,rvec x[])
{
  int c;
#ifdef GMX_MPI
  MPI_Request reqs[DD_MAXCELL],reqr[DD_MAXCELL];
#endif

  /*
  fprintf(stderr,"%d",dd->gl.nodeid);
  for(c=1; c<dd->ncell; c++) {
    fprintf(stderr," %d %3d %d %3d",
	    dd->fcell[c],dd->at_index[0],
	    dd->cell[c],dd->at_index[c]-dd->at_index[c-1]);
  }
  fprintf(stderr,"\n");
  */

#ifdef GMX_MPI
  for(c=1; c<dd->ncell; c++) {
    if (MPI_Isend(x[0],dd->at_index[0]*sizeof(rvec),MPI_BYTE,
		 dd->fcell[c],c,dd->gl.all,&reqs[c-1]) != 0)
      gmx_fatal(FARGS,"MPI_Send failed from %d to %d",
		dd->nodeid,dd->fcell[c]);
    //fprintf(stderr,"Recving %4d to %4d from %d to %d\n",
    //    dd->at_index[c-1],dd->at_index[c],dd->cell[c],dd->nodeid);
    if (MPI_Irecv(x[dd->at_index[c-1]],
		 (dd->at_index[c]-dd->at_index[c-1])*sizeof(rvec),MPI_BYTE,
		 dd->cell[c],c,dd->gl.all,&reqr[c-1]) != 0)
      gmx_fatal(FARGS,"MPI_Send failed from %d to %d",
		dd->cell[c],dd->nodeid);
  }

  MPI_Waitall(dd->ncell-1,reqs,MPI_STATUSES_IGNORE);
  MPI_Waitall(dd->ncell-1,reqr,MPI_STATUSES_IGNORE);
#endif
}

void dd_move_f(gmx_domdec_t *dd,rvec f[],rvec buf[])
{
  int c,i;
#ifdef GMX_MPI
  MPI_Request reqs[DD_MAXCELL];

  for(c=1; c<dd->ncell; c++) {

    /*
    fprintf(stderr,"node %d send %d %d recv %d %d\n",
	    dd->nodeid,
	    dd->cell[c],dd->at_index[c]-dd->at_index[c-1],
	    dd->fcell[c],dd->at_index[0]);
    */
    if (MPI_Isend(f[dd->at_index[c-1]],
		  (dd->at_index[c]-dd->at_index[c-1])*sizeof(rvec),MPI_BYTE,
		  dd->cell[c],c,dd->gl.all,&reqs[c-1]) != 0)
      gmx_fatal(FARGS,"MPI_Isend failed from %d to %d",
		dd->nodeid,dd->fcell[c]);

    if (MPI_Recv(buf[0],dd->at_index[0]*sizeof(rvec),MPI_BYTE,
		 dd->fcell[c],c,dd->gl.all,MPI_STATUS_IGNORE) != 0)
      gmx_fatal(FARGS,"MPI_Irecv failed from %d to %d",
		dd->cell[c],dd->nodeid);
    for(i=0; i<dd->at_index[0]; i++)
      rvec_inc(f[i],buf[i]);
  }

  MPI_Waitall(dd->ncell-1,reqs,MPI_STATUSES_IGNORE);
#endif
}

void dd_collect_vec(gmx_domdec_t *dd,t_block *cgs,rvec *lv,rvec *v)
{
  int  n,i,c,a;
  rvec *buf;
#ifdef GMX_MPI
  MPI_Request mpi_req;
#endif

  if (!DDMASTER(dd)) {
#ifdef GMX_MPI
    if (MPI_Send(lv,dd->at_index[0]*sizeof(rvec),MPI_BYTE,DDMASTERNODE,
		 dd->nodeid,dd->gl.all) != 0)
      gmx_fatal(FARGS,"MPI_Send from %d to %d failed",
		dd->nodeid,DDMASTERNODE);
#endif
  } else {
    /* Copy the master coordinates to the global array */
    n = DDMASTERNODE;
    a = 0;
    for(i=dd->gl.index[n]; i<dd->gl.index[n+1]; i++)
	for(c=cgs->index[dd->gl.cg[i]]; c<cgs->index[dd->gl.cg[i]+1]; c++)
	  copy_rvec(lv[a++],v[c]);

    /* Use the unused part of lv as a temporary buffer */
    buf = lv + dd->gl.nat[DDMASTERNODE];

    for(n=0; n<dd->gl.nnodes; n++) {
      if (n != DDMASTERNODE) {
#ifdef GMX_MPI
	if (MPI_Recv(buf,dd->gl.nat[n]*sizeof(rvec),MPI_BYTE,n,
		     MPI_ANY_TAG,dd->gl.all,MPI_STATUS_IGNORE) != 0)
	  gmx_fatal(FARGS,"MPI_Recv from %d to %d failed",DDMASTERNODE);
#endif
	a = 0;
	for(i=dd->gl.index[n]; i<dd->gl.index[n+1]; i++)
	  for(c=cgs->index[dd->gl.cg[i]]; c<cgs->index[dd->gl.cg[i]+1]; c++)
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
#ifdef GMX_MPI
  MPI_Request mpi_req;
#endif

  if (DDMASTER(dd)) {
    for(n=0; n<dd->gl.nnodes; n++) {
      if (n != DDMASTERNODE) {
	/* Use lv as a temporary buffer */
	a = 0;
	for(i=dd->gl.index[n]; i<dd->gl.index[n+1]; i++)
	  for(c=cgs->index[dd->gl.cg[i]]; c<cgs->index[dd->gl.cg[i]+1]; c++)
	    copy_rvec(v[c],lv[a++]);
	i = cgs->index[dd->gl.cg[dd->gl.index[n]]];
#ifdef GMX_MPI
	//if (MPI_Isend(lx,a*sizeof(rvec),MPI_BYTE,n,n,dd->gl.all,&mpi_req) != 0)
	//  gmx_fatal(FARGS,"MPI_Isend from %d to %d failed",dd->gl.master,n);
	if (MPI_Send(lv,a*sizeof(rvec),MPI_BYTE,n,n,dd->gl.all) != 0)
	  gmx_fatal(FARGS,"MPI_Send from %d to %d failed",DDMASTERNODE);
#endif
      }
    }
    n = DDMASTERNODE;
    a = 0;
    for(i=dd->gl.index[n]; i<dd->gl.index[n+1]; i++)
      for(c=cgs->index[dd->gl.cg[i]]; c<cgs->index[dd->gl.cg[i]+1]; c++)
	copy_rvec(v[c],lv[a++]);
  } else {
#ifdef GMX_MPI
    //fprintf(stderr,"%d nat %d\n",dd->nodeid,dd->at_index[0]);
    MPI_Recv(lv,dd->at_index[0]*sizeof(rvec),MPI_BYTE,DDMASTERNODE,
	     MPI_ANY_TAG,dd->gl.all,MPI_STATUS_IGNORE);
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

static void distribute_cg(FILE *fplog,matrix box,t_block *cgs,rvec pos[],
			  gmx_domdec_t *dd)
{
  static int **tmp_ind=NULL,*tmp_nalloc=NULL;
#define TMP_BLOCK_SIZE 1000

  int  i,icg,ai,k,k0,k1,d;
  rvec g_inv,cg_cm;
  ivec ind;
  real nrcg,inv_ncg;
  atom_id *cga,*cgindex;

  if (tmp_ind == NULL) {
    snew(tmp_nalloc,dd->gl.nnodes);
    snew(tmp_ind,dd->gl.nnodes);
    for(i=0; i<dd->gl.nnodes; i++) {
      tmp_nalloc[i] = (cgs->nr/TMP_BLOCK_SIZE + 2)*TMP_BLOCK_SIZE;
      snew(tmp_ind[i],tmp_nalloc[i]);
    }
  }

  /* Clear the count */
  for(i=0; i<dd->gl.nnodes; i++) {
    dd->gl.ncg[i] = 0;
    dd->gl.nat[i] = 0;
  }
  
  for(d=0; (d<DIM); d++)
    g_inv[d] = divide(dd->gl.nc[d],box[d][d]);

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
      for(k=k0; (k<k1); k++)  {
	for(d=0; (d<DIM); d++)
	  cg_cm[d] += pos[k][d];
      }
      for(d=0; (d<DIM); d++)
	cg_cm[d] = inv_ncg*cg_cm[d];
    }
    for(d=DIM-1; d>=0; d--) {
      if (cg_cm[d] < 0) 
	rvec_inc(cg_cm,box[d]);
      ind[d] = cg_cm[d]*g_inv[d];
      if (ind[d] >= dd->gl.nc[d])
	ind[d] = dd->gl.nc[d] - 1;
    }
    i = dd_index(dd->gl.nc,ind);
    if (dd->gl.ncg[i] == tmp_nalloc[i]) {
      tmp_nalloc[i] += TMP_BLOCK_SIZE;
      srenew(tmp_ind[i],tmp_nalloc[i]);
    }
    tmp_ind[i][dd->gl.ncg[i]] = icg;
    dd->gl.ncg[i]++;
    dd->gl.nat[i] += cgindex[icg+1] - cgindex[icg];
  }
  
  dd->gl.index[0] = 0;
  k1 = 0;
  for(i=0; i<dd->gl.nnodes; i++) {
    dd->gl.index[i] = k1;
    for(k=0; k<dd->gl.ncg[i]; k++)
      dd->gl.cg[k1++] = tmp_ind[i][k];
  }
  dd->gl.index[dd->gl.nnodes] = k1;

  fprintf(fplog,"Charge group distribution:");
  for(i=0; i<dd->gl.nnodes; i++)
    fprintf(fplog," %d",dd->gl.ncg[i]);
  fprintf(fplog,"\n");
}

void get_cg_distribution(FILE *fplog,gmx_domdec_t *dd,
			 matrix box,t_block *cgs,rvec pos[])
{
  int i;

  if (DDMASTER(dd)) {
    distribute_cg(fplog,box,cgs,pos,dd);
  }
#ifdef GMX_MPI
  MPI_Bcast(dd->gl.index,(dd->gl.nnodes+1+cgs->nr)*sizeof(int),MPI_BYTE,
	    DDMASTERNODE,dd->gl.all);
#endif
  for(i=0; i<dd->gl.nnodes; i++)
    dd->gl.ncg[i] = dd->gl.index[i+1] - dd->gl.index[i];
}

void setup_dd_grid(FILE *fplog,matrix box,gmx_domdec_t *dd)
{
  int  ndiv,d,i,m;
  ivec h,s;
  int  ncell,ncellp;
  ivec dd_c[8],dd_cp[4];
  
  ndiv = 0;
  for(d=0; d<DIM; d++)
    if (dd->gl.nc[d] > 1)
      ndiv++;
  
  fprintf(stderr,"Making %dD domain decomposition %d x %d x %d\n",
	  ndiv,dd->gl.nc[XX],dd->gl.nc[YY],dd->gl.nc[ZZ]);
  if (ndiv == 3) {
    ncell = dd_c3n;
    for(i=0; i<ncell; i++)
      copy_ivec(dd_c3[i],dd_c[i]);
    ncellp = dd_cp3n;
    for(i=0; i<ncellp; i++)
      copy_ivec(dd_cp3[i],dd_cp[i]);
  } else if (ndiv == 2) {
    ncell = dd_c2n;
    for(i=0; i<ncell; i++) {
      m = 0;
      for(d=0; d<DIM; d++) {
	if (dd->gl.nc[d] > 1)
	  dd_c[i][d] = dd_c2[i][m++];
	else 
	  dd_c[i][d] = 0;
      }
    }
    ncellp = dd_cp2n;
    for(i=0; i<ncellp; i++)
      copy_ivec(dd_cp2[i],dd_cp[i]);
  } else if (ndiv == 1) {
    ncell = dd_c1n;
    for(i=0; i<ncell; i++) {
      m = 0;
      for(d=0; d<DIM; d++) {
	if (dd->gl.nc[d] > 1)
	  dd_c[i][d] = dd_c1[i][m++];
	else 
	  dd_c[i][d] = 0;
      }
    }
    ncellp = dd_cp1n;
    for(i=0; i<ncellp; i++)
      copy_ivec(dd_cp1[i],dd_cp[i]);
  } else {
    gmx_fatal(FARGS,"Can only do 1, 2 or 3D domain decomposition");
    ncell = 0;
    ncellp = 0;
  }

  dd->ncell  = ncell;
  index2xyz(dd->gl.nc,dd->nodeid,h);
  for(i=0; i<ncell; i++) {
    for(d=0; d<DIM; d++) {
      s[d] = h[d] + dd_c[i][d];
      /* Currently only works for rectangular boxes.
       * For even numbers of grid cells it is easy to extend
       * to things like rhombic dodecahedrons.
       */
      if (s[d] < 0)
	s[d] += dd->gl.nc[d];
      else if (s[d] >= dd->gl.nc[d])
	s[d] -= dd->gl.nc[d];
    }
    dd->cell[i] = dd_index(dd->gl.nc,s);
    for(d=0; d<DIM; d++) {
      s[d] = h[d] - dd_c[i][d];
      /* Currently only works for rectangular boxes.
       * For even numbers of grid cells it is easy to extend
       * to things like rhombic dodecahedrons.
       */
      if (s[d] < 0)
	s[d] += dd->gl.nc[d];
      else if (s[d] >= dd->gl.nc[d])
	s[d] -= dd->gl.nc[d];
    }
    dd->fcell[i] = dd_index(dd->gl.nc,s);
  }
  dd->nicell = ncellp;
  for(i=0; i<dd->nicell; i++) {
    if (dd_cp[i][0] != i)
      gmx_fatal(FARGS,"Internal inconsistency in the dd grid setuo");
    dd->icellj0[i] = dd_cp[i][1];
    dd->icellj1[i] = dd_cp[i][2];
  }
  //fprintf(stderr,"%d dd->icellj1 %d\n",dd->nodeid,dd->icellj1[dd->nicell-1]);
}

gmx_domdec_t *init_domain_decomposition(FILE *fplog,t_commrec *cr,
					ivec nc,int ncg,int natoms,
					matrix box)
{
  gmx_domdec_t *dd;

  fprintf(fplog,"Will do domain decomposition\n");
  
  snew(dd,1);

  copy_ivec(nc,dd->gl.nc);
  dd->nodeid    = cr->nodeid;
  dd->gl.nnodes = cr->nnodes - cr->npmenodes;
#ifdef GMX_MPI
  dd->gl.all    = MPI_COMM_WORLD;
#endif
  snew(dd->gl.ncg,dd->gl.nnodes);
  snew(dd->gl.nat,dd->gl.nnodes);
  snew(dd->gl.index,dd->gl.nnodes+1+ncg);
  dd->gl.cg = dd->gl.index + dd->gl.nnodes + 1;
  snew(dd->ga2la,natoms);

  return dd;
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

  nhome = dd->at_index[0];

  n   = 0;
  ia  = il->iatoms;
  lia = lil->iatoms;
  for(i=0; i<il->nr; i+=1+nral) {
    ia = il->iatoms + i;
    a = dd->ga2la[ia[1]];
    //fprintf(stderr,"%d ftype %d nral %d i %d a[0] %d\n",
    //	    dd->nodeid,ftype,nral,i,a);
    if (a >= 0 && a < nhome) {
      lia[0] = ia[0];
      lia[1] = a;
      for(j=2; j<=nral; j++) {
	a = dd->ga2la[ia[j]];
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
  for(i=0; i<MAXNODES; i++)
    if (i < dd->nodeid)
      lil->multinr[i] = 0;
    else
      lil->multinr[i] = lil->nr;

  if (debug && il->nr)
    fprintf(debug,"ftype %d nr %d of %d\n",ftype,lil->nr,il->nr);
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
    make_local_ilist(dd,f,&idef->il[f],&lidef->il[f]);
}

static void make_local_cgs(gmx_domdec_t *dd,t_block *cgs,t_block *lcgs)
{
  int i,j,natoms_orig,orig_cg;
  atom_id a,a_orig;

  if (lcgs->index == NULL) {
    /* In most cases we could do with far less memory */
    snew(lcgs->index,cgs->nr+1);
  }

  for(i=0; i<MAXNODES; i++) {
    if (i <  dd->nodeid)
      lcgs->multinr[i] = 0;
    else
      lcgs->multinr[i] = dd->gl.ncg[dd->nodeid];
  }

  natoms_orig = cgs->index[cgs->nr];
  for(i=0; i<natoms_orig; i++)
    dd->ga2la[i] = -1;

  lcgs->nr = 0;
  a = 0;
  for(i=0; i<dd->ncell; i++) {
    for(j=dd->gl.index[dd->cell[i]]; j<dd->gl.index[dd->cell[i]+1]; j++) {
      lcgs->index[lcgs->nr] = a;
      orig_cg = dd->gl.cg[j];
      for(a_orig=cgs->index[orig_cg]; a_orig<cgs->index[orig_cg+1]; a_orig++) {
	dd->ga2la[a_orig] = a;
	a++;
      }
      lcgs->nr++;
    }
    dd->at_index[i] = a;
  }
  lcgs->index[lcgs->nr] = a;
  lcgs->nra = a;
}

static void set_cg_boundaries(gmx_domdec_t *dd,t_block *lcgs)
{
  int cgend[DD_MAXCELL+1],c;

  dd->ncg_tot = 0;
  cgend[0] = 0;
  for(c=0; c<dd->ncell; c++)
    cgend[c+1] = cgend[c] + dd->gl.ncg[dd->cell[c]];

  dd->ncg_tot = cgend[dd->ncell];

  for(c=0; c<dd->nicell; c++) {
    dd->icellcg1[c] = cgend[c+1];
    dd->icelljcg0[c] = cgend[dd->icellj0[c]];
    dd->icelljcg1[c] = cgend[dd->icellj1[c]];
  }
  /*
  fprintf(stderr,"%d ncg %d\n",dd->nodeid,dd->gl.ncg[dd->nodeid]);
  fprintf(stderr,"%d icellj1 %d\n",dd->nodeid,dd->icellj1[0]);
  fprintf(stderr,"%d icelljcg1 %d\n",dd->nodeid,dd->icelljcg1[0]);
  */
}

static void dd_update_ns_border(gmx_domdec_t *dd,t_nsborder *nsb)
{
  nsb->nodeid    = 0 /* dd->nodeid */;
  nsb->cgtotal   = dd_ncg_tot(dd);
  //nsb->natoms    = dd->gl.nat[dd->nodeid];
  nsb->index[nsb->nodeid]  = 0;
  nsb->homenr[nsb->nodeid] = dd->at_index[0];
  nsb->cgload[nsb->nodeid] = dd->icelljcg1[0];
}

void make_local_top(FILE *fplog,gmx_domdec_t *dd,
		    t_topology *top,t_topology *ltop,
		    t_nsborder *nsb)
{
  int a,i,cg,eb;

  ltop->name  = top->name;
  make_local_cgs(dd,&top->blocks[ebCGS],&ltop->blocks[ebCGS]);

  make_local_idef(dd,&top->idef,&ltop->idef);
 
  ltop->atoms = top->atoms;
  ltop->atomtypes = top->atomtypes;
  for(eb=0; eb<ebNR; eb++) {
    if (eb != ebCGS)
      ltop->blocks[eb] = top->blocks[eb];
  }
  ltop->symtab = top->symtab;

  set_cg_boundaries(dd,&ltop->blocks[ebCGS]);

  dd_update_ns_border(dd,nsb);
}
