#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "smalloc.h"
#include "vec.h"
#include "constr.h"
#include "domdec.h"

typedef struct {
  int nsend;
  int *a;
  int a_nalloc;
  int nrecv;
} gmx_specatsend_t;

typedef struct gmx_domdec_specat_comm {
  /* The atom indices we need from the surrounding cells */
  int  nind_req;
  int  *ind_req;
  int  ind_req_nalloc;
  /* The number of indices to receive during the setup */
  int  nreq[DIM][2][2];
  /* The atoms to send */
  gmx_specatsend_t spas[DIM][2];
  bool *bSendAtom;
  int   bSendAtom_nalloc;
  /* Send buffers */
  int  *ibuf;
  int  ibuf_nalloc;
  rvec *vbuf;
  int  vbuf_nalloc;
  rvec *vbuf2;
  int  vbuf2_nalloc;
} gmx_domdec_specat_comm_t;

#define SPECAT_ALLOC_SIZE     1000
#define IATOM_ALLOC_SIZE      1000
#define CONSTRAINT_ALLOC_SIZE 1000

static void dd_move_f_specat(gmx_domdec_t *dd,gmx_domdec_specat_comm_t *spac,
			     int end,rvec *f,rvec *fshift)
{
  gmx_specatsend_t *spas;
  rvec *vbuf;
  int  n,n0,n1,d,dir,i;
  ivec vis;
  int  is;
  
  n = end;
  for(d=dd->ndim-1; d>=0; d--) {
    if (dd->nc[dd->dim[d]] > 2) {
      /* Pulse the grid forward and backward */
      spas = spac->spas[d];
      n0 = spas[0].nrecv;
      n1 = spas[1].nrecv;
      n -= n1 + n0;
      vbuf = spac->vbuf;
      /* Send and receive the coordinates */
      dd_sendrecv2_rvec(dd,d,
			f+n   ,n0,vbuf              ,spas[0].nsend,
			f+n+n0,n1,vbuf+spas[0].nsend,spas[1].nsend);
      for(dir=0; dir<2; dir++) {
	/* Sum the buffer into the required forces */
	if (fshift &&
	    ((dir == 0 && dd->ci[dd->dim[d]] == 0) || 
	     (dir == 1 && dd->ci[dd->dim[d]] == dd->nc[dd->dim[d]]-1))) {
	  clear_ivec(vis);
	  vis[dd->dim[d]] = (dir==0 ? 1 : -1);
	  is = IVEC2IS(vis);
	  for(i=0; i<spas->nsend; i++) {
	    rvec_inc(f[spas->a[i]],*vbuf);
	    rvec_inc(fshift[is],*vbuf);
	    vbuf++;
	  }
	} else {
	  for(i=0; i<spas->nsend; i++) {
	    rvec_inc(f[spas->a[i]],*vbuf);
	    vbuf++;
	  }
	}
      }
    } else {
      /* Two cells, so we only need to communicate one way */
      spas = &spac->spas[d][0];
      n -= spas->nrecv;
      /* Send and receive the coordinates */
      dd_sendrecv_rvec(dd,d,ddForward,
		       f+n,spas->nrecv,spac->vbuf,spas->nsend);
      /* Sum the buffer into the required forces */
      for(i=0; i<spas->nsend; i++)
	rvec_inc(f[spas->a[i]],spac->vbuf[i]);
    }
  }
}

void dd_move_f_vsites(gmx_domdec_t *dd,rvec *f,rvec *fshift)
{
  if (dd->vsite_comm)
    dd_move_f_specat(dd,dd->vsite_comm,dd->nat_tot_vsite,f,fshift);
}

static void dd_move_x_specat(gmx_domdec_t *dd,gmx_domdec_specat_comm_t *spac,
			     matrix box,int start,rvec *x0,rvec *x1)
{
  gmx_specatsend_t *spas;
  rvec *x,*vbuf,*rbuf;
  int  nvec,v,n,nn,ns0,ns1,nr0,nr1,nr,d,dir,i;
  rvec shift;
  
  nvec = 1;
  if (x1)
    nvec++;

  n = start;
  for(d=0; d<dd->ndim; d++) {
    if (dd->nc[dd->dim[d]] > 2) {
      /* Pulse the grid forward and backward */
      vbuf = spac->vbuf;
      for(dir=0; dir<2; dir++) {
	spas = &spac->spas[d][dir];
	for(v=0; v<nvec; v++) {
	  x = (v == 0 ? x0 : x1);
	  /* Copy the required coordinates to the send buffer */
	  if (dir == 0 &&
	      dd->ci[dd->dim[d]] == 0) {
	    copy_rvec(box[dd->dim[d]],shift);
	    for(i=0; i<spas->nsend; i++) {
	      rvec_add(x[spas->a[i]],shift,*vbuf);
	      vbuf++;
	    }
	  } else if (dir == 1 &&
		     dd->ci[dd->dim[d]] == dd->nc[dd->dim[d]]-1) {
	    copy_rvec(box[dd->dim[d]],shift);
	    for(i=0; i<spas->nsend; i++) {
	      rvec_sub(x[spas->a[i]],shift,*vbuf);
	      vbuf++;
	    }
	  } else {
	    for(i=0; i<spas->nsend; i++) {
	      copy_rvec(x[spas->a[i]],*vbuf);
	      vbuf++;
	    }
	  }
	}
      }
      /* Send and receive the coordinates */
      spas = spac->spas[d];
      ns0 = spas[0].nsend;
      nr0 = spas[0].nrecv;
      ns1 = spas[1].nsend;
      nr1 = spas[1].nrecv;
      if (nvec == 1) {
	dd_sendrecv2_rvec(dd,d,
			  spac->vbuf+ns0,ns1,x0+n    ,nr1,
			  spac->vbuf    ,ns0,x0+n+nr1,nr0);
      } else {
	/* Communicate both vectors in one buffer */
	rbuf = spac->vbuf2;
	dd_sendrecv2_rvec(dd,d,
			  spac->vbuf+2*ns0,2*ns1,rbuf      ,2*nr1,
			  spac->vbuf      ,2*ns0,rbuf+2*nr1,2*nr0);
	/* Split the buffer into the two vectors */
	nn = n;
	for(dir=1; dir>=0; dir--) {
	  nr = spas[dir].nrecv;
	  for(v=0; v<2; v++) {
	    x = (v == 0 ? x0 : x1);
	    for(i=0; i<nr; i++) {
	      copy_rvec(*rbuf,x[nn+i]);
	      rbuf++;
	    }
	  }
	  nn += nr;
	}
      }
      n += nr0 + nr1;
    } else {
      spas = &spac->spas[d][0];
      /* Copy the required coordinates to the send buffer */
      vbuf = spac->vbuf;
      for(v=0; v<nvec; v++) {
	x = (v == 0 ? x0 : x1);
	for(i=0; i<spas->nsend; i++) {
	  copy_rvec(x[spas->a[i]],*vbuf);
	  vbuf++;
	}
      }
      /* Send and receive the coordinates */
      if (nvec == 1) {
	dd_sendrecv_rvec(dd,d,ddBackward,
			 spac->vbuf,spas->nsend,x0+n,spas->nrecv);
      } else {
	/* Communicate both vectors in one buffer */
	rbuf = spac->vbuf2;
	dd_sendrecv_rvec(dd,d,ddBackward,
			 spac->vbuf,2*spas->nsend,rbuf,2*spas->nrecv);
	/* Split the buffer into the two vectors */
	nr = spas[0].nrecv;
	for(v=0; v<2; v++) {
	  x = (v == 0 ? x0 : x1);
	  for(i=0; i<nr; i++) {
	    copy_rvec(*rbuf,x[n+i]);
	    rbuf++;
	  }
	}
      }
      n += spas->nrecv;
    }
  }
}

void dd_move_x_constraints(gmx_domdec_t *dd,matrix box,rvec *x0,rvec *x1)
{
  if (dd->constraint_comm)
    dd_move_x_specat(dd,dd->constraint_comm,box,dd->nat_tot_vsite,x0,x1);
}

void dd_move_x_vsites(gmx_domdec_t *dd,matrix box,rvec *x)
{
  if (dd->vsite_comm)
    dd_move_x_specat(dd,dd->vsite_comm,box,dd->nat_tot,x,NULL);
}

void dd_clear_local_constraint_indices(gmx_domdec_t *dd)
{
  gmx_domdec_constraints_t *dc;
  int i;
  
  dc = dd->constraints;
  
  for(i=0; i<dc->ncon; i++)
    dc->gc2lc[dc->con[i]] = -1;
  
  for(i=dd->nat_tot_vsite; i<dd->nat_tot_con; i++)
    dc->ga2la[dd->gatindex[i]] = -1;
}

void dd_clear_local_vsite_indices(gmx_domdec_t *dd)
{
  int i;

  for(i=dd->nat_tot; i<dd->nat_tot_vsite; i++)
    dd->ga2la_vsite[dd->gatindex[i]] = -1;
}

static int setup_specat_communication(gmx_domdec_t *dd,
				      gmx_domdec_specat_comm_t *spac,
				      int *ga2la_specat,
				      int at_start,
				      int vbuf_fac,
				      char *specat_type,char *add_err)
{
  int  d,dir,nsend[2],nlast,ndir,nr,ns,i,nrecv_local,n0,start,ireq,ind,buf[2];
  int  nat_tot_specat,nat_tot_prev,nalloc_old;
  bool bFirst;
  gmx_specatsend_t *spas;

  if (debug)
    fprintf(debug,"Begin setup_specat_communication for %s\n",specat_type);

  /* nsend[0]: the number of atoms requested by this node only,
   *           we communicate this for more efficients checks
   * nsend[1]: the total number of requested atoms
   */
  nsend[0] = spac->nind_req;
  nsend[1] = nsend[0];
  nlast    = nsend[1];
  for(d=dd->ndim-1; d>=0; d--) {
    /* Pulse the grid forward and backward */
    if (dd->nc[dd->dim[d]] > 2)
      ndir = 2;
    else
      ndir = 1;
    for(dir=0; dir<ndir; dir++) {
      /* Communicate the number of indices */
      dd_sendrecv_int(dd,d,dir==0 ? ddForward : ddBackward,
		      nsend,2,spac->nreq[d][dir],2);
      nr = spac->nreq[d][dir][1];
      if (nlast+nr > spac->ind_req_nalloc) {
	spac->ind_req_nalloc = over_alloc(nlast+nr);
	srenew(spac->ind_req,spac->ind_req_nalloc);
      }
      /* Communicate the indices */
      dd_sendrecv_int(dd,d,dir==0 ? ddForward : ddBackward,
		      spac->ind_req,nsend[1],spac->ind_req+nlast,nr);
      nlast += nr;
    }
    nsend[1] = nlast;
  }
  if (debug)
    fprintf(debug,"Communicated the counts\n");

  /* Search for the requested atoms and communicate the indices we have */
  nat_tot_specat = at_start;
  nrecv_local = 0;
  for(d=0; d<dd->ndim; d++) {
    bFirst = (d == 0);
    /* Pulse the grid forward and backward */
    if (dd->nc[dd->dim[d]] > 2)
      ndir = 2;
    else
      ndir = 1;
    nat_tot_prev = nat_tot_specat;
    for(dir=ndir-1; dir>=0; dir--) {
      if (nat_tot_specat > spac->bSendAtom_nalloc) {
	nalloc_old = spac->bSendAtom_nalloc;
	spac->bSendAtom_nalloc = over_alloc(nat_tot_specat);
	srenew(spac->bSendAtom,spac->bSendAtom_nalloc);
	for(i=nalloc_old; i<spac->bSendAtom_nalloc; i++)
	  spac->bSendAtom[i] = FALSE;
      }
      spas = &spac->spas[d][dir];
      n0 = spac->nreq[d][dir][0];
      nr = spac->nreq[d][dir][1];
      if (debug)
	fprintf(debug,"dim=%d, dir=%d, searching for %d atoms\n",d,dir,nr);
      start = nlast - nr;
      spas->nsend = 0;
      nsend[0] = 0;
      for(i=0; i<nr; i++) {
	ireq = spac->ind_req[start+i];
	ind = -1;
	if (dd->ga2la[ireq].cell == 0) {
	  /* We have this atom locally */
	  ind = dd->ga2la[ireq].a;
	} else if (!bFirst) {
	  /* Search in the communicated atoms */
	  ind = ga2la_specat[ireq];
	}
	if (ind >= 0) {
	  if (i < n0 || !spac->bSendAtom[ind]) {
	    if (spas->nsend >= spas->a_nalloc) {
	      spas->a_nalloc += SPECAT_ALLOC_SIZE;
	      srenew(spas->a,spas->a_nalloc);
	    }
	    /* Store the local index so we know which coordinates
	     * to send out later.
	     */
	    spas->a[spas->nsend] = ind;
	    spac->bSendAtom[ind] = TRUE;
	    if (spas->nsend >= spac->ibuf_nalloc) {
	      spac->ibuf_nalloc += SPECAT_ALLOC_SIZE;
	      srenew(spac->ibuf,spac->ibuf_nalloc);
	    }
	    /* Store the global index so we can send it now */
	    spac->ibuf[spas->nsend] = ireq;
	    if (i < n0)
	      nsend[0]++;
	    spas->nsend++;
	  }
	}
      }
      nlast = start;
      /* Clear the local flags */
      for(i=0; i<spas->nsend; i++)
	spac->bSendAtom[spas->a[i]] = FALSE;
      /* Send and receive the number of indices to communicate */
      nsend[1] = spas->nsend;
      dd_sendrecv_int(dd,d,dir==0 ? ddBackward : ddForward,
		      nsend,2,buf,2);
      if (debug) {
	fprintf(debug,"Send to node %d, %d (%d) indices, "
		"receive from node %d, %d (%d) indices\n",
		dd->neighbor[d][1-dir],nsend[1],nsend[0],
		dd->neighbor[d][dir],buf[1],buf[0]);
	if (gmx_debug_at) {
	  for(i=0; i<spas->nsend; i++)
	    fprintf(debug," %d",spac->ibuf[i]+1);
	  fprintf(debug,"\n");
	}
      }
      nrecv_local += buf[0];
      spas->nrecv  = buf[1];
      if (nat_tot_specat + spas->nrecv > dd->gatindex_nalloc) {
	dd->gatindex_nalloc = over_alloc(nat_tot_specat + spas->nrecv);
	srenew(dd->gatindex,dd->gatindex_nalloc);
      }
      /* Send and receive the indices */
      dd_sendrecv_int(dd,d,dir==0 ? ddBackward : ddForward,
		      spac->ibuf,spas->nsend,
		      dd->gatindex+nat_tot_specat,spas->nrecv);
      nat_tot_specat += spas->nrecv;
    }

    /* Allocate the x/f communication buffers */
    ns = spac->spas[d][0].nsend;
    nr = spac->spas[d][0].nrecv;
    if (ndir == 2) {
      ns += spac->spas[d][1].nsend;
      nr += spac->spas[d][1].nrecv;
    }
    if (vbuf_fac*ns > spac->vbuf_nalloc) {
      spac->vbuf_nalloc = over_alloc(vbuf_fac*ns);
      srenew(spac->vbuf,spac->vbuf_nalloc);
    }
    if (vbuf_fac == 2 && vbuf_fac*nr > spac->vbuf2_nalloc) {
      spac->vbuf2_nalloc = over_alloc(vbuf_fac*nr);
      srenew(spac->vbuf2,spac->vbuf2_nalloc);
    }

    /* Make a global to local index for the communication atoms */
    for(i=nat_tot_prev; i<nat_tot_specat; i++)
      ga2la_specat[dd->gatindex[i]] = i;
  }
  
  /* Check that in the end we got the number of atoms we asked for */
  if (nrecv_local != spac->nind_req) {
    if (debug) {
      fprintf(debug,"Requested %d, received %d (tot recv %d)\n",
	      spac->nind_req,nrecv_local,nat_tot_specat-at_start);
      if (gmx_debug_at) {
	for(i=0; i<spac->nind_req; i++)
	  fprintf(debug," %s%d",
		  ga2la_specat[spac->ind_req[i]]>=0 ? "" : "!",
		  spac->ind_req[i]+1);
	fprintf(debug,"\n");
      }
    }
    gmx_fatal(FARGS,"Node %d could only obtain %d of the %d atoms that are connected via %ss from the neighboring cells. This probably means you %s lengths are too long compared to the domain decomposition cell size. Decrease the number of domain decomposition grid cells%s.",
	      dd->sim_nodeid,nrecv_local,spac->nind_req,specat_type,
	      specat_type,add_err);
  }

  if (debug)
    fprintf(debug,"Done setup_specat_communication\n");

  return nat_tot_specat;
}

static void walk_out(int con,int a,int nrec,const t_iatom *ia,
		     const gmx_ga2la_t *ga2la,bool bHomeConnect,
		     gmx_domdec_constraints_t *dc,
		     gmx_domdec_specat_comm_t *dcc)
{
  int i,coni,b;

  if (dc->gc2lc[con] == -1) {
    /* Add this non-home constraint to the list */
    if (dc->ncon >= dc->con_nalloc) {
      dc->con_nalloc += CONSTRAINT_ALLOC_SIZE;
      srenew(dc->con,dc->con_nalloc);
      srenew(dc->con_nlocat,dc->con_nalloc);
    }
    dc->con[dc->ncon] = con;
    dc->con_nlocat[dc->ncon] = (bHomeConnect ? 1 : 0);
    dc->gc2lc[con] = dc->ncon;
    dc->ncon++;
  }
  /* Check to not ask for the same atom more than once */
  if (dc->ga2la[a] == -1) {
    /* Add this non-home atom to the list */
    if (dcc->nind_req >= dcc->ind_req_nalloc) {
      dcc->ind_req_nalloc += SPECAT_ALLOC_SIZE;
      srenew(dcc->ind_req,dcc->ind_req_nalloc);
    }
    dcc->ind_req[dcc->nind_req++] = a;
    /* Temporarily mark with -2, we get the index later */
    dc->ga2la[a] = -2;
  }

  if (nrec > 0) {
    for(i=dc->at2con.index[a]; i<dc->at2con.index[a+1]; i++) {
      coni = dc->at2con.a[i];
      if (coni != con) {
	/* Walk further */
	if (a == ia[coni*3+1])
	  b = ia[coni*3+2];
	else
	  b = ia[coni*3+1];
	if (ga2la[b].cell != 0)
	  walk_out(coni,b,nrec-1,ia,ga2la,FALSE,dc,dcc);
      }
    }
  }
}

void dd_make_local_constraints(gmx_domdec_t *dd,t_iatom *ia,int nrec)
{
  t_block at2con;
  gmx_ga2la_t *ga2la;
  t_iatom *iap;
  int nhome,a,ag,bg,i,con;
  gmx_domdec_constraints_t *dc;

  dc = dd->constraints;

  at2con = dc->at2con;
  ga2la  = dd->ga2la;

  dc->ncon = 0;
  nhome = 0;
  if (dd->constraint_comm)
    dd->constraint_comm->nind_req = 0;
  for(a=0; a<dd->nat_home; a++) {
    ag = dd->gatindex[a];
    for(i=at2con.index[ag]; i<at2con.index[ag+1]; i++) {
      con = at2con.a[i];
      iap = ia + con*3;
      if (ag == iap[1]) {
	bg = iap[2];
      } else {
	bg = iap[1];
      }
      if (ga2la[bg].cell == 0) {
	/* Add this fully home constraint at the first atom */
	if (ag < bg) {
	  if (dc->ncon >= dc->con_nalloc) {
	    dc->con_nalloc += CONSTRAINT_ALLOC_SIZE;
	    srenew(dc->con,dc->con_nalloc);
	    srenew(dc->con_nlocat,dc->con_nalloc);
	  }
	  dc->con[dc->ncon] = con;
	  dc->con_nlocat[dc->ncon] = 2;
	  dc->gc2lc[con] = dc->ncon;
	  dc->ncon++;
	  nhome++;
	}
      } else {
	/* We need to walk out of the home cell by nrec constraints */
	walk_out(con,bg,nrec,ia,dd->ga2la,TRUE,dc,dd->constraint_comm);
      }
    }
  }

  if (debug)
    fprintf(debug,
	    "Constraints: home %3d border %3d atoms: %3d\n",
	    nhome,dc->ncon-nhome,
	    dd->constraint_comm ? dd->constraint_comm->nind_req : 0);

  if (dd->constraint_comm)
    dd->nat_tot_con =
      setup_specat_communication(dd,dd->constraint_comm,dd->constraints->ga2la,
				 dd->nat_tot_vsite,2,"constraint",
				 " or lincs-order");
}

void dd_make_local_vsites(gmx_domdec_t *dd,t_ilist *lil)
{
  gmx_domdec_specat_comm_t *spac;
  int  *ga2la_specat;
  int  ftype,nral,i,j,gat,a;
  t_ilist *lilf;
  t_iatom *iatoms;

  spac         = dd->vsite_comm;
  ga2la_specat = dd->ga2la_vsite;

  spac->nind_req = 0;
  /* Loop over all the home vsites */
  for(ftype=0; ftype<F_NRE; ftype++) {
    if (interaction_function[ftype].flags & IF_VSITE) {
      nral = NRAL(ftype);
      lilf = &lil[ftype];
      for(i=0; i<lilf->nr; i+=1+nral) {
	iatoms = lilf->iatoms + i;
	/* Check if we have the other atoms */
	for(j=1; j<1+nral; j++) {
	  if (iatoms[j] < 0) {
	    /* This is not a home atom, we need to ask our neighbors */
	    a = -iatoms[j] - 1;
	    /* Check to not ask for the same atom more than once */
	    if (ga2la_specat[a] == -1) {
	      /* Add this non-home atom to the list */
	      if (spac->nind_req >= spac->ind_req_nalloc) {
		spac->ind_req_nalloc += SPECAT_ALLOC_SIZE;
		srenew(spac->ind_req,spac->ind_req_nalloc);
	      }
	      spac->ind_req[spac->nind_req++] = a;
	      /* Temporarily mark with -2, we get the index later */
	      ga2la_specat[a] = -2;
	    }
	  }
	}
      }
    }
  }

  dd->nat_tot_vsite =
    setup_specat_communication(dd,dd->vsite_comm,ga2la_specat,
			       dd->nat_tot,1,"vsite","");

  /* Fill in the missing indices */
  for(ftype=0; ftype<F_NRE; ftype++) {
    if (interaction_function[ftype].flags & IF_VSITE) {
      nral = NRAL(ftype);
      lilf = &lil[ftype];
      for(i=0; i<lilf->nr; i+=1+nral) {
	iatoms = lilf->iatoms + i;
	for(j=1; j<1+nral; j++) {
	  if (iatoms[j] < 0)
	    iatoms[j] = ga2la_specat[-iatoms[j]-1];
	}
      }
    }
  }
}

void init_domdec_constraints(gmx_domdec_t *dd,
			     int natoms,t_idef *idef,t_block *cgs,
			     bool bDynamics)
{
  int c,cg,cg0,cg1,a,a2;
  gmx_domdec_constraints_t *dc;
  bool bInterCGConstraints;

  if (debug)
    fprintf(debug,"Begin init_domdec_constraints\n");

  snew(dd->constraints,1);
  dc = dd->constraints;

  dc->at2con = make_at2con(0,natoms,idef,bDynamics,
			   &dc->ncon_global,&dc->nflexcon_global);
  dc->iatoms = idef->il[F_CONSTR].iatoms;

  snew(dc->gc2lc,idef->il[F_CONSTR].nr/3);
  for(c=0; c<idef->il[F_CONSTR].nr/3; c++)
    dc->gc2lc[c] = -1;

  snew(dc->ga2la,natoms);
  for(a=0; a<natoms; a++)
    dc->ga2la[a] = -1;

  bInterCGConstraints = FALSE;
  for(cg=0; cg<cgs->nr && !bInterCGConstraints; cg++) {
    cg0 = cgs->index[cg];
    cg1 = cgs->index[cg+1];
    for(a=cg0; a<cg1; a++)
      for(c=dc->at2con.index[a]; c<dc->at2con.index[a+1]; c++) {
	/* Single all constraints are stored for both atoms
	 * we only need to check the second atom.
	 */
	a2 = dc->iatoms[3*dc->at2con.a[c]+2];
	if (a2 < cg0 || a2 >= cg1)
	  bInterCGConstraints = TRUE;
      }
  }
  
  if (bInterCGConstraints)
    snew(dd->constraint_comm,1);
}

void init_domdec_vsites(gmx_domdec_t *dd,int natoms)
{
  int i;
  gmx_domdec_constraints_t *dc;
  
  if (debug)
    fprintf(debug,"Begin init_domdec_vsites\n");

  snew(dd->ga2la_vsite,natoms);
  for(i=0; i<natoms; i++)
    dd->ga2la_vsite[i] = -1;

  snew(dd->vsite_comm,1);
}
