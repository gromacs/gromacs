/* -*- mode: c; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; c-file-style: "stroustrup"; -*-
 *
 * 
 * This file is part of Gromacs        Copyright (c) 1991-2008
 * David van der Spoel, Erik Lindahl, Berk Hess, University of Groningen.
 *
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 *
 * To help us fund GROMACS development, we humbly ask that you cite
 * the research papers on the package. Check out http://www.gromacs.org
 * 
 * And Hey:
 * Gnomes, ROck Monsters And Chili Sauce
 */

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "smalloc.h"
#include "vec.h"
#include "constr.h"
#include "domdec.h"
#include "domdec_network.h"
#include "mtop_util.h"
#include "gmx_ga2la.h"

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
    gmx_bool *bSendAtom;
    int   bSendAtom_nalloc;
    /* Send buffers */
    int  *ibuf;
    int  ibuf_nalloc;
    rvec *vbuf;
    int  vbuf_nalloc;
    rvec *vbuf2;
    int  vbuf2_nalloc;
    /* The range in the local buffer(s) for received atoms */
    int  at_start;
    int  at_end;
} gmx_domdec_specat_comm_t;

typedef struct gmx_domdec_constraints {
    int  *molb_con_offset;
    int  *molb_ncon_mol;
    /* The fully local and connected constraints */
    int  ncon;
    /* The global constraint number, only required for clearing gc_req */
    int  *con_gl;
    int  *con_nlocat;
    int  con_nalloc;
    /* Boolean that tells if a global constraint index has been requested */
    char *gc_req;
    /* Global to local communicated constraint atom only index */
    int  *ga2la;
} gmx_domdec_constraints_t;


static void dd_move_f_specat(gmx_domdec_t *dd,gmx_domdec_specat_comm_t *spac,
                             rvec *f,rvec *fshift)
{
    gmx_specatsend_t *spas;
    rvec *vbuf;
    int  n,n0,n1,d,dim,dir,i;
    ivec vis;
    int  is;
    gmx_bool bPBC,bScrew;
    
    n = spac->at_end;
    for(d=dd->ndim-1; d>=0; d--)
    {
        dim = dd->dim[d];
        if (dd->nc[dim] > 2)
        {
            /* Pulse the grid forward and backward */
            spas = spac->spas[d];
            n0 = spas[0].nrecv;
            n1 = spas[1].nrecv;
            n -= n1 + n0;
            vbuf = spac->vbuf;
            /* Send and receive the coordinates */
            dd_sendrecv2_rvec(dd,d,
                              f+n+n1,n0,vbuf              ,spas[0].nsend,
                              f+n   ,n1,vbuf+spas[0].nsend,spas[1].nsend);
            for(dir=0; dir<2; dir++)
            {
                bPBC   = ((dir == 0 && dd->ci[dim] == 0) || 
                          (dir == 1 && dd->ci[dim] == dd->nc[dim]-1));
                bScrew = (bPBC && dd->bScrewPBC && dim == XX);
                
                spas = &spac->spas[d][dir];
                /* Sum the buffer into the required forces */
                if (!bPBC || (!bScrew && fshift == NULL))
                {
                    for(i=0; i<spas->nsend; i++)
                    {
                        rvec_inc(f[spas->a[i]],*vbuf);
                        vbuf++;
                    }
                }
                else
                {
                    clear_ivec(vis);
                    vis[dim] = (dir==0 ? 1 : -1);
                    is = IVEC2IS(vis);
                    if (!bScrew)
                    {
                        /* Sum and add to shift forces */
                        for(i=0; i<spas->nsend; i++)
                        {
                            rvec_inc(f[spas->a[i]],*vbuf);
                            rvec_inc(fshift[is],*vbuf);
                            vbuf++;
                        }
                    }
                    else
                    {	    
                        /* Rotate the forces */
                        for(i=0; i<spas->nsend; i++)
                        {
                            f[spas->a[i]][XX] += (*vbuf)[XX];
                            f[spas->a[i]][YY] -= (*vbuf)[YY];
                            f[spas->a[i]][ZZ] -= (*vbuf)[ZZ];
                            if (fshift)
                            {
                                rvec_inc(fshift[is],*vbuf);
                            }
                            vbuf++;
                        }
                    }
                }
            }
        }
        else
        {
            /* Two cells, so we only need to communicate one way */
            spas = &spac->spas[d][0];
            n -= spas->nrecv;
            /* Send and receive the coordinates */
            dd_sendrecv_rvec(dd,d,dddirForward,
                             f+n,spas->nrecv,spac->vbuf,spas->nsend);
            /* Sum the buffer into the required forces */
            if (dd->bScrewPBC && dim == XX &&
                (dd->ci[dim] == 0 ||
                 dd->ci[dim] == dd->nc[dim]-1))
            {
                for(i=0; i<spas->nsend; i++)
                {
                    /* Rotate the force */
                    f[spas->a[i]][XX] += spac->vbuf[i][XX];
                    f[spas->a[i]][YY] -= spac->vbuf[i][YY];
                    f[spas->a[i]][ZZ] -= spac->vbuf[i][ZZ];
                }
            }
            else
            {
                for(i=0; i<spas->nsend; i++)
                {
                    rvec_inc(f[spas->a[i]],spac->vbuf[i]);
                }
            }
        }
    }
}

void dd_move_f_vsites(gmx_domdec_t *dd,rvec *f,rvec *fshift)
{
    if (dd->vsite_comm)
    {
        dd_move_f_specat(dd,dd->vsite_comm,f,fshift);
    }
}

void dd_clear_f_vsites(gmx_domdec_t *dd,rvec *f)
{
    int i;
    
    if (dd->vsite_comm)
    {
        for(i=dd->vsite_comm->at_start; i<dd->vsite_comm->at_end; i++)
        {
            clear_rvec(f[i]);
        }
    }
}

static void dd_move_x_specat(gmx_domdec_t *dd,gmx_domdec_specat_comm_t *spac,
                             matrix box,rvec *x0,rvec *x1)
{
    gmx_specatsend_t *spas;
    rvec *x,*vbuf,*rbuf;
    int  nvec,v,n,nn,ns0,ns1,nr0,nr1,nr,d,dim,dir,i;
    gmx_bool bPBC,bScrew=FALSE;
    rvec shift={0,0,0};
    
    nvec = 1;
    if (x1)
    {
        nvec++;
    }
    
    n = spac->at_start;
    for(d=0; d<dd->ndim; d++)
    {
        dim = dd->dim[d];
        if (dd->nc[dim] > 2)
        {
            /* Pulse the grid forward and backward */
            vbuf = spac->vbuf;
            for(dir=0; dir<2; dir++)
            {
                if (dir == 0 && dd->ci[dim] == 0)
                {
                    bPBC   = TRUE;
                    bScrew = (dd->bScrewPBC && dim == XX);
                    copy_rvec(box[dim],shift);
                }
                else if (dir == 1 && dd->ci[dim] == dd->nc[dim]-1)
                {
                    bPBC = TRUE;
                    bScrew = (dd->bScrewPBC && dim == XX);
                    for(i=0; i<DIM; i++)
                    {
                        shift[i] = -box[dim][i];
                    }
                }
                else
                {
                    bPBC = FALSE;
                    bScrew = FALSE;
                }
                spas = &spac->spas[d][dir];
                for(v=0; v<nvec; v++)
                {
                    x = (v == 0 ? x0 : x1);
                    /* Copy the required coordinates to the send buffer */
                    if (!bPBC)
                    {
                        /* Only copy */
                        for(i=0; i<spas->nsend; i++)
                        {
                            copy_rvec(x[spas->a[i]],*vbuf);
                            vbuf++;
                        }
                    }
                    else if (!bScrew)
                    {
                        /* Shift coordinates */
                        for(i=0; i<spas->nsend; i++)
                        {
                            rvec_add(x[spas->a[i]],shift,*vbuf);
                            vbuf++;
                        }
                    }
                    else
                    {
                        /* Shift and rotate coordinates */
                        for(i=0; i<spas->nsend; i++)
                        {
                            (*vbuf)[XX] =               x[spas->a[i]][XX] + shift[XX];
                            (*vbuf)[YY] = box[YY][YY] - x[spas->a[i]][YY] + shift[YY];
                            (*vbuf)[ZZ] = box[ZZ][ZZ] - x[spas->a[i]][ZZ] + shift[ZZ];
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
            if (nvec == 1)
            {
                dd_sendrecv2_rvec(dd,d,
                                  spac->vbuf+ns0,ns1,x0+n    ,nr1,
                                  spac->vbuf    ,ns0,x0+n+nr1,nr0);
            }
            else
            {
                /* Communicate both vectors in one buffer */
                rbuf = spac->vbuf2;
                dd_sendrecv2_rvec(dd,d,
                                  spac->vbuf+2*ns0,2*ns1,rbuf      ,2*nr1,
                                  spac->vbuf      ,2*ns0,rbuf+2*nr1,2*nr0);
                /* Split the buffer into the two vectors */
                nn = n;
                for(dir=1; dir>=0; dir--)
                {
                    nr = spas[dir].nrecv;
                    for(v=0; v<2; v++)
                    {
                        x = (v == 0 ? x0 : x1);
                        for(i=0; i<nr; i++)
                        {
                            copy_rvec(*rbuf,x[nn+i]);
                            rbuf++;
                        }
                    }
                    nn += nr;
                }
            }
            n += nr0 + nr1;
        }
        else
        {
            spas = &spac->spas[d][0];
            /* Copy the required coordinates to the send buffer */
            vbuf = spac->vbuf;
            for(v=0; v<nvec; v++)
            {
                x = (v == 0 ? x0 : x1);
                if (dd->bScrewPBC && dim == XX &&
                    (dd->ci[XX] == 0 || dd->ci[XX] == dd->nc[XX]-1))
                {
                    /* Here we only perform the rotation, the rest of the pbc
                     * is handled in the constraint or viste routines.
                     */
                    for(i=0; i<spas->nsend; i++)
                    {
                        (*vbuf)[XX] =               x[spas->a[i]][XX];
                        (*vbuf)[YY] = box[YY][YY] - x[spas->a[i]][YY];
                        (*vbuf)[ZZ] = box[ZZ][ZZ] - x[spas->a[i]][ZZ];
                        vbuf++;
                    }	  
                }
                else
                {
                    for(i=0; i<spas->nsend; i++)
                    {
                        copy_rvec(x[spas->a[i]],*vbuf);
                        vbuf++;
                    }
                }
            }
            /* Send and receive the coordinates */
            if (nvec == 1)
            {
                dd_sendrecv_rvec(dd,d,dddirBackward,
                                 spac->vbuf,spas->nsend,x0+n,spas->nrecv);
            }
            else
            {
                /* Communicate both vectors in one buffer */
                rbuf = spac->vbuf2;
                dd_sendrecv_rvec(dd,d,dddirBackward,
                                 spac->vbuf,2*spas->nsend,rbuf,2*spas->nrecv);
                /* Split the buffer into the two vectors */
                nr = spas[0].nrecv;
                for(v=0; v<2; v++)
                {
                    x = (v == 0 ? x0 : x1);
                    for(i=0; i<nr; i++)
                    {
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
    {
        dd_move_x_specat(dd,dd->constraint_comm,box,x0,x1);
    }
}

void dd_move_x_vsites(gmx_domdec_t *dd,matrix box,rvec *x)
{
    if (dd->vsite_comm)
    {
        dd_move_x_specat(dd,dd->vsite_comm,box,x,NULL);
    }
}

int *dd_constraints_nlocalatoms(gmx_domdec_t *dd)
{
    if (dd->constraints)
    {
        return dd->constraints->con_nlocat;
    }
    else
    {
        return NULL;
    }
}

void dd_clear_local_constraint_indices(gmx_domdec_t *dd)
{
    gmx_domdec_constraints_t *dc;
    int i;
    
    dc = dd->constraints;
    
    for(i=0; i<dc->ncon; i++)
    {
        dc->gc_req[dc->con_gl[i]] = 0;
    }
  
    if (dd->constraint_comm)
    {
        for(i=dd->constraint_comm->at_start; i<dd->constraint_comm->at_end; i++)
        {
            dc->ga2la[dd->gatindex[i]] = -1;
        }
    }
}

void dd_clear_local_vsite_indices(gmx_domdec_t *dd)
{
    int i;
    
    if (dd->vsite_comm)
    {
        for(i=dd->vsite_comm->at_start; i<dd->vsite_comm->at_end; i++)
        {
            dd->ga2la_vsite[dd->gatindex[i]] = -1;
        }
    }
}

static int setup_specat_communication(gmx_domdec_t *dd,
                                      gmx_domdec_specat_comm_t *spac,
                                      int *ga2la_specat,
                                      int at_start,
                                      int vbuf_fac,
                                      const char *specat_type,
                                      const char *add_err)
{
    int  nsend[2],nlast,nsend_zero[2]={0,0},*nsend_ptr;
    int  d,dim,ndir,dir,nr,ns,i,nrecv_local,n0,start,ireq,ind,buf[2];
    int  nat_tot_specat,nat_tot_prev,nalloc_old;
    gmx_bool bPBC,bFirst;
    gmx_specatsend_t *spas;
    
    if (debug)
    {
        fprintf(debug,"Begin setup_specat_communication for %s\n",specat_type);
    }
    
    /* nsend[0]: the number of atoms requested by this node only,
     *           we communicate this for more efficients checks
     * nsend[1]: the total number of requested atoms
     */
    nsend[0] = spac->nind_req;
    nsend[1] = nsend[0];
    nlast    = nsend[1];
    for(d=dd->ndim-1; d>=0; d--)
    {
        /* Pulse the grid forward and backward */
        dim = dd->dim[d];
        bPBC = (dim < dd->npbcdim);
        if (dd->nc[dim] == 2)
        {
            /* Only 2 cells, so we only need to communicate once */
            ndir = 1;
        }
        else
        {
            ndir = 2;
        }
        for(dir=0; dir<ndir; dir++)
        {
            if (!bPBC && 
                dd->nc[dim] > 2 &&
                ((dir == 0 && dd->ci[dim] == dd->nc[dim] - 1) ||
                 (dir == 1 && dd->ci[dim] == 0)))
            {
                /* No pbc: the fist/last cell should not request atoms */
                nsend_ptr = nsend_zero;
            }
            else
            {
                nsend_ptr = nsend;
            }
            /* Communicate the number of indices */
            dd_sendrecv_int(dd,d,dir==0 ? dddirForward : dddirBackward,
                            nsend_ptr,2,spac->nreq[d][dir],2);
            nr = spac->nreq[d][dir][1];
            if (nlast+nr > spac->ind_req_nalloc)
            {
                spac->ind_req_nalloc = over_alloc_dd(nlast+nr);
                srenew(spac->ind_req,spac->ind_req_nalloc);
            }
            /* Communicate the indices */
            dd_sendrecv_int(dd,d,dir==0 ? dddirForward : dddirBackward,
                            spac->ind_req,nsend_ptr[1],spac->ind_req+nlast,nr);
            nlast += nr;
        }
        nsend[1] = nlast;
    }
    if (debug)
    {
        fprintf(debug,"Communicated the counts\n");
    }
    
    /* Search for the requested atoms and communicate the indices we have */
    nat_tot_specat = at_start;
    nrecv_local = 0;
    for(d=0; d<dd->ndim; d++)
    {
        bFirst = (d == 0);
        /* Pulse the grid forward and backward */
        if (dd->dim[d] >= dd->npbcdim || dd->nc[dd->dim[d]] > 2)
        {
            ndir = 2;
        }
        else
        {
            ndir = 1;
        }
        nat_tot_prev = nat_tot_specat;
        for(dir=ndir-1; dir>=0; dir--)
        {
            if (nat_tot_specat > spac->bSendAtom_nalloc)
            {
                nalloc_old = spac->bSendAtom_nalloc;
                spac->bSendAtom_nalloc = over_alloc_dd(nat_tot_specat);
                srenew(spac->bSendAtom,spac->bSendAtom_nalloc);
                for(i=nalloc_old; i<spac->bSendAtom_nalloc; i++)
                {
                    spac->bSendAtom[i] = FALSE;
                }
            }
            spas = &spac->spas[d][dir];
            n0 = spac->nreq[d][dir][0];
            nr = spac->nreq[d][dir][1];
            if (debug)
            {
                fprintf(debug,"dim=%d, dir=%d, searching for %d atoms\n",
                        d,dir,nr);
            }
            start = nlast - nr;
            spas->nsend = 0;
            nsend[0] = 0;
            for(i=0; i<nr; i++)
            {
                ireq = spac->ind_req[start+i];
                ind = -1;
                /* Check if this is a home atom and if so ind will be set */
                if (!ga2la_get_home(dd->ga2la,ireq,&ind))
                {
                    /* Search in the communicated atoms */
                    ind = ga2la_specat[ireq];
                }
                if (ind >= 0)
                {
                    if (i < n0 || !spac->bSendAtom[ind])
                    {
                        if (spas->nsend+1 > spas->a_nalloc)
                        {
                            spas->a_nalloc = over_alloc_large(spas->nsend+1);
                            srenew(spas->a,spas->a_nalloc);
                        }
                        /* Store the local index so we know which coordinates
                         * to send out later.
                         */
                        spas->a[spas->nsend] = ind;
                        spac->bSendAtom[ind] = TRUE;
                        if (spas->nsend+1 > spac->ibuf_nalloc)
                        {
                            spac->ibuf_nalloc = over_alloc_large(spas->nsend+1);
                            srenew(spac->ibuf,spac->ibuf_nalloc);
                        }
                        /* Store the global index so we can send it now */
                        spac->ibuf[spas->nsend] = ireq;
                        if (i < n0)
                        {
                            nsend[0]++;
                        }
                        spas->nsend++;
                    }
                }
            }
            nlast = start;
            /* Clear the local flags */
            for(i=0; i<spas->nsend; i++)
            {
                spac->bSendAtom[spas->a[i]] = FALSE;
            }
            /* Send and receive the number of indices to communicate */
            nsend[1] = spas->nsend;
            dd_sendrecv_int(dd,d,dir==0 ? dddirBackward : dddirForward,
                            nsend,2,buf,2);
            if (debug)
            {
                fprintf(debug,"Send to node %d, %d (%d) indices, "
                        "receive from node %d, %d (%d) indices\n",
                        dd->neighbor[d][1-dir],nsend[1],nsend[0],
                        dd->neighbor[d][dir],buf[1],buf[0]);
                if (gmx_debug_at)
                {
                    for(i=0; i<spas->nsend; i++)
                    {
                        fprintf(debug," %d",spac->ibuf[i]+1);
                    }
                    fprintf(debug,"\n");
                }
            }
            nrecv_local += buf[0];
            spas->nrecv  = buf[1];
            if (nat_tot_specat + spas->nrecv > dd->gatindex_nalloc)
            {
                dd->gatindex_nalloc =
                    over_alloc_dd(nat_tot_specat + spas->nrecv);
                srenew(dd->gatindex,dd->gatindex_nalloc);
            }
            /* Send and receive the indices */
            dd_sendrecv_int(dd,d,dir==0 ? dddirBackward : dddirForward,
                            spac->ibuf,spas->nsend,
                            dd->gatindex+nat_tot_specat,spas->nrecv);
            nat_tot_specat += spas->nrecv;
        }
        
        /* Allocate the x/f communication buffers */
        ns = spac->spas[d][0].nsend;
        nr = spac->spas[d][0].nrecv;
        if (ndir == 2)
        {
            ns += spac->spas[d][1].nsend;
            nr += spac->spas[d][1].nrecv;
        }
        if (vbuf_fac*ns > spac->vbuf_nalloc)
        {
            spac->vbuf_nalloc = over_alloc_dd(vbuf_fac*ns);
            srenew(spac->vbuf,spac->vbuf_nalloc);
        }
        if (vbuf_fac == 2 && vbuf_fac*nr > spac->vbuf2_nalloc)
        {
            spac->vbuf2_nalloc = over_alloc_dd(vbuf_fac*nr);
            srenew(spac->vbuf2,spac->vbuf2_nalloc);
        }
        
        /* Make a global to local index for the communication atoms */
        for(i=nat_tot_prev; i<nat_tot_specat; i++)
        {
            ga2la_specat[dd->gatindex[i]] = i;
        }
    }
    
    /* Check that in the end we got the number of atoms we asked for */
    if (nrecv_local != spac->nind_req)
    {
        if (debug)
        {
            fprintf(debug,"Requested %d, received %d (tot recv %d)\n",
                    spac->nind_req,nrecv_local,nat_tot_specat-at_start);
            if (gmx_debug_at)
            {
                for(i=0; i<spac->nind_req; i++)
                {
                    fprintf(debug," %s%d",
                            ga2la_specat[spac->ind_req[i]]>=0 ? "" : "!",
                            spac->ind_req[i]+1);
                }
                fprintf(debug,"\n");
            }
        }
        fprintf(stderr,"\nDD cell %d %d %d: Neighboring cells do not have atoms:",
                dd->ci[XX],dd->ci[YY],dd->ci[ZZ]);
        for(i=0; i<spac->nind_req; i++)
        {
            if (ga2la_specat[spac->ind_req[i]] < 0)
            {
                fprintf(stderr," %d",spac->ind_req[i]+1);
            }
        }
        fprintf(stderr,"\n");
        gmx_fatal(FARGS,"DD cell %d %d %d could only obtain %d of the %d atoms that are connected via %ss from the neighboring cells. This probably means your %s lengths are too long compared to the domain decomposition cell size. Decrease the number of domain decomposition grid cells%s%s.",
                  dd->ci[XX],dd->ci[YY],dd->ci[ZZ],
                  nrecv_local,spac->nind_req,specat_type,
                  specat_type,add_err,
                  dd->bGridJump ? " or use the -rcon option of mdrun" : "");
    }
    
    spac->at_start = at_start;
    spac->at_end   = nat_tot_specat;
    
    if (debug)
    {
        fprintf(debug,"Done setup_specat_communication\n");
    }
    
    return nat_tot_specat;
}

static void walk_out(int con,int con_offset,int a,int offset,int nrec,
                     int ncon1,const t_iatom *ia1,const t_iatom *ia2,
                     const t_blocka *at2con,
                     const gmx_ga2la_t ga2la,gmx_bool bHomeConnect,
                     gmx_domdec_constraints_t *dc,
                     gmx_domdec_specat_comm_t *dcc,
                     t_ilist *il_local)
{
    int a1_gl,a2_gl,a_loc,i,coni,b;
    const t_iatom *iap;
  
    if (dc->gc_req[con_offset+con] == 0)
    {
        /* Add this non-home constraint to the list */
        if (dc->ncon+1 > dc->con_nalloc)
        {
            dc->con_nalloc = over_alloc_large(dc->ncon+1);
            srenew(dc->con_gl,dc->con_nalloc);
            srenew(dc->con_nlocat,dc->con_nalloc);
        }
        dc->con_gl[dc->ncon] = con_offset + con;
        dc->con_nlocat[dc->ncon] = (bHomeConnect ? 1 : 0);
        dc->gc_req[con_offset+con] = 1;
        if (il_local->nr + 3 > il_local->nalloc)
        {
            il_local->nalloc = over_alloc_dd(il_local->nr+3);
            srenew(il_local->iatoms,il_local->nalloc);
        }
        iap = constr_iatomptr(ncon1,ia1,ia2,con);
        il_local->iatoms[il_local->nr++] = iap[0];
        a1_gl = offset + iap[1];
        a2_gl = offset + iap[2];
        /* The following indexing code can probably be optizimed */
        if (ga2la_get_home(ga2la,a1_gl,&a_loc))
        {
            il_local->iatoms[il_local->nr++] = a_loc;
        }
        else
        {
            /* We set this index later */
            il_local->iatoms[il_local->nr++] = -a1_gl - 1;
        }
        if (ga2la_get_home(ga2la,a2_gl,&a_loc))
        {
            il_local->iatoms[il_local->nr++] = a_loc;
        }
        else
        {
            /* We set this index later */
            il_local->iatoms[il_local->nr++] = -a2_gl - 1;
        }
        dc->ncon++;
    }
    /* Check to not ask for the same atom more than once */
    if (dc->ga2la[offset+a] == -1)
    {
        /* Add this non-home atom to the list */
        if (dcc->nind_req+1 > dcc->ind_req_nalloc)
        {
            dcc->ind_req_nalloc = over_alloc_large(dcc->nind_req+1);
            srenew(dcc->ind_req,dcc->ind_req_nalloc);
        }
        dcc->ind_req[dcc->nind_req++] = offset + a;
        /* Temporarily mark with -2, we get the index later */
        dc->ga2la[offset+a] = -2;
    }
    
    if (nrec > 0)
    {
        for(i=at2con->index[a]; i<at2con->index[a+1]; i++)
        {
            coni = at2con->a[i];
            if (coni != con)
            {
                /* Walk further */
                iap = constr_iatomptr(ncon1,ia1,ia2,coni);
                if (a == iap[1])
                {
                    b = iap[2];
                }
                else
                {
                    b = iap[1];
                }
                if (!ga2la_get_home(ga2la,offset+b,&a_loc))
                {
                    walk_out(coni,con_offset,b,offset,nrec-1,
                             ncon1,ia1,ia2,at2con,
                             ga2la,FALSE,dc,dcc,il_local);
                }
            }
        }
    }
}

int dd_make_local_constraints(gmx_domdec_t *dd,int at_start,
                              gmx_mtop_t *mtop,
                              gmx_constr_t constr,int nrec,
                              t_ilist *il_local)
{
    t_blocka *at2con_mt,*at2con;
    gmx_ga2la_t ga2la;
    int ncon1,ncon2;
    gmx_molblock_t *molb;
    t_iatom *ia1,*ia2,*iap;
    int nhome,a,a_gl,a_mol,a_loc,b_lo,offset,mb,molnr,b_mol,i,con,con_offset;
    gmx_domdec_constraints_t *dc;
    int at_end,*ga2la_specat,j;
    
    dc = dd->constraints;
    
    at2con_mt = atom2constraints_moltype(constr);
    ga2la  = dd->ga2la;
    
    dc->ncon     = 0;
    il_local->nr = 0;
    nhome = 0;
    if (dd->constraint_comm)
    {
        dd->constraint_comm->nind_req = 0;
    }
    for(a=0; a<dd->nat_home; a++)
    {
        a_gl = dd->gatindex[a];
        
        gmx_mtop_atomnr_to_molblock_ind(mtop,a_gl,&mb,&molnr,&a_mol);
        molb = &mtop->molblock[mb];
        
        ncon1 = mtop->moltype[molb->type].ilist[F_CONSTR].nr/3;
        ncon2 = mtop->moltype[molb->type].ilist[F_CONSTRNC].nr/3;
        if (ncon1 > 0 || ncon2 > 0)
        {
            ia1 = mtop->moltype[molb->type].ilist[F_CONSTR].iatoms;
            ia2 = mtop->moltype[molb->type].ilist[F_CONSTRNC].iatoms;

            /* Calculate the global constraint number offset for the molecule.
             * This is only required for the global index to make sure
             * that we use each constraint only once.
             */
            con_offset = dc->molb_con_offset[mb] + molnr*dc->molb_ncon_mol[mb];
            
            /* The global atom number offset for this molecule */
            offset = a_gl - a_mol;
            at2con = &at2con_mt[molb->type];
            for(i=at2con->index[a_mol]; i<at2con->index[a_mol+1]; i++)
            {
                con = at2con->a[i];
                iap = constr_iatomptr(ncon1,ia1,ia2,con);
                if (a_mol == iap[1])
                {
                    b_mol = iap[2];
                }
                else
                {
                    b_mol = iap[1];
                }
                if (ga2la_get_home(ga2la,offset+b_mol,&a_loc))
                {
                    /* Add this fully home constraint at the first atom */
                    if (a_mol < b_mol)
                    {
                        if (dc->ncon+1 > dc->con_nalloc)
                        {
                            dc->con_nalloc = over_alloc_large(dc->ncon+1);
                            srenew(dc->con_gl,dc->con_nalloc);
                            srenew(dc->con_nlocat,dc->con_nalloc);
                        }
                        dc->con_gl[dc->ncon] = con_offset + con;
                        dc->con_nlocat[dc->ncon] = 2;
                        if (il_local->nr + 3 > il_local->nalloc)
                        {
                            il_local->nalloc = over_alloc_dd(il_local->nr + 3);
                            srenew(il_local->iatoms,il_local->nalloc);
                        }
                        b_lo = a_loc;
                        il_local->iatoms[il_local->nr++] = iap[0];
                        il_local->iatoms[il_local->nr++] = (a_gl == iap[1] ? a    : b_lo);
                        il_local->iatoms[il_local->nr++] = (a_gl == iap[1] ? b_lo : a   );
                        dc->ncon++;
                        nhome++;
                    }
                }
                else
                {
                    /* We need the nrec constraints coupled to this constraint,
                     * so we need to walk out of the home cell by nrec+1 atoms,
                     * since already atom bg is not locally present.
                     * Therefore we call walk_out with nrec recursions to go
                     * after this first call.
                     */
                    walk_out(con,con_offset,b_mol,offset,nrec,
                             ncon1,ia1,ia2,at2con,
                             dd->ga2la,TRUE,dc,dd->constraint_comm,il_local);
                }
            }
        }
    }
    
    if (debug)
    {
        fprintf(debug,
                "Constraints: home %3d border %3d atoms: %3d\n",
                nhome,dc->ncon-nhome,
                dd->constraint_comm ? dd->constraint_comm->nind_req : 0);
    }

    if (dd->constraint_comm) {
        at_end =
            setup_specat_communication(dd,dd->constraint_comm,
                                       dd->constraints->ga2la,
                                       at_start,2,
                                       "constraint"," or lincs-order");
        
        /* Fill in the missing indices */
        ga2la_specat = dd->constraints->ga2la;
        for(i=0; i<il_local->nr; i+=3)
        {
            iap = il_local->iatoms + i;
            for(j=1; j<3; j++)
            {
                if (iap[j] < 0)
                {
                    iap[j] = ga2la_specat[-iap[j]-1];
                }
            }
        }
    }
    else
    {
        at_end = at_start;
    }
    
    return at_end;
}

int dd_make_local_vsites(gmx_domdec_t *dd,int at_start,t_ilist *lil)
{
    gmx_domdec_specat_comm_t *spac;
    int  *ga2la_specat;
    int  ftype,nral,i,j,gat,a;
    t_ilist *lilf;
    t_iatom *iatoms;
    int  at_end;
    
    spac         = dd->vsite_comm;
    ga2la_specat = dd->ga2la_vsite;
    
    spac->nind_req = 0;
    /* Loop over all the home vsites */
    for(ftype=0; ftype<F_NRE; ftype++)
    {
        if (interaction_function[ftype].flags & IF_VSITE)
        {
            nral = NRAL(ftype);
            lilf = &lil[ftype];
            for(i=0; i<lilf->nr; i+=1+nral)
            {
                iatoms = lilf->iatoms + i;
                /* Check if we have the other atoms */
                for(j=1; j<1+nral; j++)
                {
                    if (iatoms[j] < 0) {
                        /* This is not a home atom,
                         * we need to ask our neighbors.
                         */
                        a = -iatoms[j] - 1;
                        /* Check to not ask for the same atom more than once */
                        if (ga2la_specat[a] == -1)
                        {
                            /* Add this non-home atom to the list */
                            if (spac->nind_req+1 > spac->ind_req_nalloc)
                            {
                                spac->ind_req_nalloc =
                                    over_alloc_small(spac->nind_req+1);
                                srenew(spac->ind_req,spac->ind_req_nalloc);
                            }
                            spac->ind_req[spac->nind_req++] = a;
                            /* Temporarily mark with -2,
                             * we get the index later.
                             */
                            ga2la_specat[a] = -2;
                        }
                    }
                }
            }
        }
    }
    
    at_end = setup_specat_communication(dd,dd->vsite_comm,ga2la_specat,
                                        at_start,1,"vsite","");
    
    /* Fill in the missing indices */
    for(ftype=0; ftype<F_NRE; ftype++)
    {
        if (interaction_function[ftype].flags & IF_VSITE)
        {
            nral = NRAL(ftype);
            lilf = &lil[ftype];
            for(i=0; i<lilf->nr; i+=1+nral)
            {
                iatoms = lilf->iatoms + i;
                for(j=1; j<1+nral; j++)
                {
                    if (iatoms[j] < 0)
                    {
                        iatoms[j] = ga2la_specat[-iatoms[j]-1];
                    }
                }
            }
        }
    }
    
    return at_end;
}

void init_domdec_constraints(gmx_domdec_t *dd,
                             int natoms,gmx_mtop_t *mtop,
                             gmx_constr_t constr)
{
    gmx_domdec_constraints_t *dc;
    gmx_molblock_t *molb;
    int mb,ncon,c,a;
    
    if (debug)
    {
        fprintf(debug,"Begin init_domdec_constraints\n");
    }
    
    snew(dd->constraints,1);
    dc = dd->constraints;
    
    snew(dc->molb_con_offset,mtop->nmolblock);
    snew(dc->molb_ncon_mol,mtop->nmolblock);
    
    ncon = 0;
    for(mb=0; mb<mtop->nmolblock; mb++)
    {
        molb = &mtop->molblock[mb];
        dc->molb_con_offset[mb] = ncon;
        dc->molb_ncon_mol[mb] =
            mtop->moltype[molb->type].ilist[F_CONSTR].nr/3 +
            mtop->moltype[molb->type].ilist[F_CONSTRNC].nr/3;
        ncon += molb->nmol*dc->molb_ncon_mol[mb];
    }
    
    snew(dc->gc_req,ncon);
    for(c=0; c<ncon; c++)
    {
        dc->gc_req[c] = 0;
    }
    
    snew(dc->ga2la,natoms);
    for(a=0; a<natoms; a++)
    {
        dc->ga2la[a] = -1;
    }
    
    snew(dd->constraint_comm,1);
}

void init_domdec_vsites(gmx_domdec_t *dd,int natoms)
{
    int i;
    gmx_domdec_constraints_t *dc;
    
    if (debug)
    {
        fprintf(debug,"Begin init_domdec_vsites\n");
    }
    
    snew(dd->ga2la_vsite,natoms);
    for(i=0; i<natoms; i++)
    {
        dd->ga2la_vsite[i] = -1;
    }
    
    snew(dd->vsite_comm,1);
}
