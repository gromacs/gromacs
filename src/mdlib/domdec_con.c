#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "smalloc.h"
#include "vec.h"
#include "constr.h"
#include "domdec.h"

#define CONSTRAINTS_ALLOC_SIZE 1000

void dd_move_x_constraints(gmx_domdec_t *dd,rvec *x)
{
  gmx_domdec_constraints_t *dc;
  gmx_conatomsend_t *cas;
  int n,d,ndir,dir,i;

  dc = dd->constraints;
  
  n = dd->nat_tot;
  for(d=0; d<dd->ndim; d++) {
    /* Pulse the grid forward and backward */
    if (dd->nc[dd->dim[d]] > 2)
      ndir = 2;
    else
      ndir = 1;
    for(dir=ndir-1; dir>=0; dir--) {
      cas = &dc->cas[d][dir];
      /* Copy the required coordinates to the send buffer */
      for(i=0; i<cas->nsend; i++)
	copy_rvec(x[cas->a[i]],dc->vbuf[i]);
      /* Send and receive the coordinates */
      dd_sendrecv_rvec(dd,d,dir==0 ? ddBackward : ddForward,
		       dc->vbuf,cas->nsend,x+n,cas->nrecv);
      n += cas->nrecv;
    }
  }
}

void clear_local_constraint_indices(gmx_domdec_t *dd)
{
  gmx_domdec_constraints_t *dc;
  int i;
  
  dc = dd->constraints;
  
  for(i=0; i<dc->ncon; i++)
    dc->gc2lc[dc->con[i]] = -1;
  
  for(i=dd->nat_tot; i<dd->nat_tot_con; i++)
    dc->ga2la[dd->gatindex[i]] = -1;
}

static void setup_constraint_communication(gmx_domdec_t *dd)
{
  gmx_domdec_constraints_t *dc;
  int  d,dir,nsend[2],nlast,ndir,nr,i,nrecv_local,n0,start,ireq,ind,buf[2];
  int  nat_tot_con_prev,nalloc_old;
  bool bFirst;
  gmx_conatomsend_t *cas;

  if (debug)
    fprintf(debug,"Begin setup_constraint_communication\n");

  dc = dd->constraints;
  
  /* nsend[0]: the number of atoms requested by this node only,
   *           we communicate this for more efficients checks
   * nsend[1]: the total number of requested atoms
   */
  nsend[0] = dc->nind_req;
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
		      nsend,2,dc->nreq[d][dir],2);
      nr = dc->nreq[d][dir][1];
      if (nlast+nr > dc->ind_req_nalloc) {
	dc->ind_req_nalloc = over_alloc(nlast+nr);
	srenew(dc->ind_req,dc->ind_req_nalloc);
      }
      /* Communicate the indices */
      dd_sendrecv_int(dd,d,dir==0 ? ddForward : ddBackward,
		      dc->ind_req,nsend[1],dc->ind_req+nlast,nr);
      nlast += nr;
    }
    nsend[1] = nlast;
  }
  if (debug)
    fprintf(debug,"Communicated the counts\n");

  /* Search for the requested atoms and communicate the indices we have */
  dd->nat_tot_con = dd->nat_tot;
  nrecv_local = 0;
  for(d=0; d<dd->ndim; d++) {
    bFirst = (d == 0);
    /* Pulse the grid forward and backward */
    if (dd->nc[dd->dim[d]] > 2)
      ndir = 2;
    else
      ndir = 1;
    nat_tot_con_prev = dd->nat_tot_con;
    for(dir=ndir-1; dir>=0; dir--) {
      if (dd->nat_tot_con > dc->bSendAtom_nalloc) {
	nalloc_old = dc->bSendAtom_nalloc;
	dc->bSendAtom_nalloc = over_alloc(dd->nat_tot_con);
	srenew(dc->bSendAtom,dc->bSendAtom_nalloc);
	for(i=nalloc_old; i<dc->bSendAtom_nalloc; i++)
	  dc->bSendAtom[i] = FALSE;
      }
      cas = &dc->cas[d][dir];
      n0 = dc->nreq[d][dir][0];
      nr = dc->nreq[d][dir][1];
      if (debug)
	fprintf(debug,"dim=%d, dir=%d, searching for %d atoms\n",d,dir,nr);
      start = nlast - nr;
      cas->nsend = 0;
      nsend[0] = 0;
      for(i=0; i<nr; i++) {
	ireq = dc->ind_req[start+i];
	ind = -1;
	if (dd->ga2la[ireq].cell == 0) {
	  /* We have this atom locally */
	  ind = dd->ga2la[ireq].a;
	} else if (!bFirst) {
	  /* Search in the communicated atoms */
	  ind = dc->ga2la[ireq];
	}
	if (ind >= 0) {
	  if (i < n0 || !dc->bSendAtom[ind]) {
	    if (cas->nsend >= cas->a_nalloc) {
	      cas->a_nalloc += CONSTRAINTS_ALLOC_SIZE;
	      srenew(cas->a,cas->a_nalloc);
	    }
	    /* Store the local index so we know which coordinates
	     * to send out later.
	     */
	    cas->a[cas->nsend] = ind;
	    dc->bSendAtom[ind] = TRUE;
	    if (cas->nsend >= dc->buf_nalloc) {
	      dc->buf_nalloc += CONSTRAINTS_ALLOC_SIZE;
	      srenew(dc->ibuf,dc->buf_nalloc);
	      srenew(dc->vbuf,dc->buf_nalloc);
	    }
	    /* Store the global index so we can send it now */
	    dc->ibuf[cas->nsend] = ireq;
	    if (i < n0)
	      nsend[0]++;
	    cas->nsend++;
	  }
	}
      }
      nlast = start;
      /* Clear the local flags */
      for(i=0; i<cas->nsend; i++)
	dc->bSendAtom[cas->a[i]] = FALSE;
      /* Send and receive the number of indices to communicate */
      nsend[1] = cas->nsend;
      dd_sendrecv_int(dd,d,dir==0 ? ddBackward : ddForward,
		      nsend,2,buf,2);
      if (debug) {
	fprintf(debug,"Send to node %d, %d (%d) indices, "
		"receive from node %d, %d (%d) indices\n",
		dd->neighbor[d][1-dir],nsend[1],nsend[0],
		dd->neighbor[d][dir],buf[1],buf[0]);
	for(i=0; i<cas->nsend; i++)
	  fprintf(debug," %d",dc->ibuf[i]);
	fprintf(debug,"\n");
      }
      nrecv_local += buf[0];
      cas->nrecv   = buf[1];
      if (dd->nat_tot_con + cas->nrecv > dd->gatindex_nalloc) {
	dd->gatindex_nalloc = over_alloc(dd->nat_tot_con + cas->nrecv);
	srenew(dd->gatindex,dd->gatindex_nalloc);
      }
      /* Send and receive the indices */
      dd_sendrecv_int(dd,d,dir==0 ? ddBackward : ddForward,
		      dc->ibuf,cas->nsend,
		      dd->gatindex+dd->nat_tot_con,cas->nrecv);
      dd->nat_tot_con += cas->nrecv;
    }

    /* Make a global to local index for the communication atoms */
    for(i=nat_tot_con_prev; i<dd->nat_tot_con; i++)
      dc->ga2la[dd->gatindex[i]] = i;
  }
  
  /* Check that in the end we got the number of atoms we asked for */
  if (nrecv_local != dc->nind_req) {
    if (debug) {
      fprintf(debug,"Requested %d, received %d (tot recv %d)\n",
	      dc->nind_req,nrecv_local,dd->nat_tot_con-dd->nat_tot);
      for(i=0; i<dc->nind_req; i++)
	fprintf(debug," %s%d",
		(dc->ga2la[dc->ind_req[i]]>=0 ? "" : "!"),dc->ind_req[i]);
      fprintf(debug,"\n");
    }
    gmx_fatal(FARGS,"Node %d could only obtain %d of the %d atoms that are connected via constraints from the neighboring cells. This probably means you constraint lengths are too long compared to the domain decomposition cell size. Decrease lincs-order or decrease the number of domain decomposition grid cells.",dd->nodeid,nrecv_local,dc->nind_req);
  }

  if (debug)
    fprintf(debug,"Done setup_constraint_communication\n");
}

static void walk_out(int con,int a,int nrec,const t_iatom *ia,
		     const gmx_ga2la_t *ga2la,bool bLocalConnect,
		     gmx_domdec_constraints_t *dc)
{
  int i,coni,b;

  if (dc->gc2lc[con] == -1) {
    /* Add this non-local constraint to the list */
    if (dc->ncon >= dc->con_nalloc) {
      dc->con_nalloc += CONSTRAINTS_ALLOC_SIZE;
      srenew(dc->con,dc->con_nalloc);
      srenew(dc->con_nlocat,dc->con_nalloc);
    }
    dc->con[dc->ncon] = con;
    dc->con_nlocat[dc->ncon] = (bLocalConnect ? 1 : 0);
    dc->gc2lc[con] = dc->ncon;
    dc->ncon++;
  }
  /* Check to not ask for the same atom more than once */
  if (dc->ga2la[a] == -1) {
    /* Add this non-local atom to the list */
    if (dc->nind_req >= dc->ind_req_nalloc) {
      dc->ind_req_nalloc += CONSTRAINTS_ALLOC_SIZE;
      srenew(dc->ind_req,dc->ind_req_nalloc);
    }
    dc->ind_req[dc->nind_req++] = a;
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
	  walk_out(coni,b,nrec-1,ia,ga2la,FALSE,dc);
      }
    }
  }
}

void make_local_constraints(gmx_domdec_t *dd,t_iatom *ia,int nrec,
			    gmx_domdec_constraints_t *dc)
{
  t_block at2con;
  gmx_ga2la_t *ga2la;
  t_iatom *iap;
  int nlocal,a,ag,bg,i,con;

  at2con = dc->at2con;
  ga2la  = dd->ga2la;

  dc->ncon = 0;
  nlocal = 0;
  dc->nind_req = 0;
  for(a=0; a<dd->nat_local; a++) {
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
	/* Add this fully local constraint at the first atom */
	if (ag < bg) {
	  if (dc->ncon >= dc->con_nalloc) {
	    dc->con_nalloc += CONSTRAINTS_ALLOC_SIZE;
	    srenew(dc->con,dc->con_nalloc);
	    srenew(dc->con_nlocat,dc->con_nalloc);
	  }
	  dc->con[dc->ncon] = con;
	  dc->con_nlocat[dc->ncon] = 2;
	  dc->gc2lc[con] = dc->ncon;
	  dc->ncon++;
	  nlocal++;
	}
      } else {
	/* We need to walk out of the local cell by nrec constraints */
	walk_out(con,bg,nrec,ia,dd->ga2la,TRUE,dc);
      }
    }
  }
  if (debug)
    fprintf(debug,
	    "Constraints: local %3d border %3d atoms: %3d\n",
	    nlocal,dc->ncon-nlocal,dc->nind_req);

  setup_constraint_communication(dd);
}

gmx_domdec_constraints_t *init_domdec_constraints(int natoms,t_idef *idef,
						  bool bDynamics)
{
  int i;
  gmx_domdec_constraints_t *dc;
  
  if (debug)
    fprintf(debug,"Begin init_domdec_constraints\n");

  snew(dc,1);
  dc->at2con = make_at2con(0,natoms,idef,bDynamics,
			   &dc->ncon_global,&dc->nflexcon_global);
  dc->iatoms = idef->il[F_CONSTR].iatoms;

  snew(dc->gc2lc,idef->il[F_CONSTR].nr/3);
  for(i=0; i<idef->il[F_CONSTR].nr/3; i++)
    dc->gc2lc[i] = -1;

  snew(dc->ga2la,natoms);
  for(i=0; i<natoms; i++)
    dc->ga2la[i] = -1;
  
  return dc;
}
