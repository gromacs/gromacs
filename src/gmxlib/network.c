/*
 * $Id$
 * 
 *                This source code is part of
 * 
 *                 G   R   O   M   A   C   S
 * 
 *          GROningen MAchine for Chemical Simulations
 * 
 *                        VERSION 3.2.0
 * Written by David van der Spoel, Erik Lindahl, Berk Hess, and others.
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2004, The GROMACS development team,
 * check out http://www.gromacs.org for more information.

 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 * 
 * If you want to redistribute modifications, please consider that
 * scientific software is very special. Version control is crucial -
 * bugs must be traceable. We will be happy to consider code for
 * inclusion in the official distribution, but derived work must not
 * be called official GROMACS. Details are found in the README & COPYING
 * files - if they are missing, get the official version at www.gromacs.org.
 * 
 * To help us fund GROMACS development, we humbly ask that you cite
 * the papers on the package - you can find them in the top README file.
 * 
 * For more info, check our website at http://www.gromacs.org
 * 
 * And Hey:
 * GROningen Mixture of Alchemy and Childrens' Stories
 */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <string.h>
#include "gmx_fatal.h"
#include "main.h"
#include "smalloc.h"
#include "network.h"
#include "copyrite.h"
#include "statutil.h"
#include "ctype.h"

#ifdef GMX_MPI
#include <mpi.h>
#include "mpelogging.h"
#endif

#ifdef GMX_MPI
static MPI_Request mpi_req_tx=MPI_REQUEST_NULL,mpi_req_rx;
#endif

/* Try setting MPI_TEST when you experience unexplainable crashes, *
 * up til now these crashes have only occured with IRIX 6.5        */
/* #define MPI_TEST */

void gmx_tx(const t_commrec *cr,int nodeid,void *buf,int bufsize)
{
#ifndef GMX_MPI
  gmx_call("gmx_tx"); 
#else
  int        tag,flag;
  MPI_Status status;
  
#ifdef DEBUG
  fprintf(stderr,"gmx_tx: nodeid=%d, buf=%x, bufsize=%d\n",
	  nodeid,buf,bufsize);
#endif
#ifdef MPI_TEST
  /* workaround for crashes encountered with MPI on IRIX 6.5 */
  if (mpi_req_tx != MPI_REQUEST_NULL) {
    MPI_Test(&mpi_req_tx,&flag,&status);
    if (flag==FALSE) {
      fprintf(stdlog,"gmx_tx called before previous send was complete: nodeid=%d, buf=%x, bufsize=%d\n",
	      nodeid,buf,bufsize);
      gmx_tx_wait(nodeid);
    }
  }
#endif
  tag = 0;
  if (MPI_Isend(buf,bufsize,MPI_BYTE,RANK(cr,nodeid),tag,cr->mpi_comm_mygroup,&mpi_req_tx) != 0)
    gmx_comm("MPI_Isend Failed");
#endif
}

void gmx_tx_wait(int nodeid)
{
#ifndef GMX_MPI
  gmx_call("gmx_tx_wait");
#else
  MPI_Status  status;
  int mpi_result;
  
  if ((mpi_result=MPI_Wait(&mpi_req_tx,&status)) != 0)
    gmx_fatal(FARGS,"MPI_Wait: result=%d",mpi_result);
#endif
}

void gmx_txs(const t_commrec *cr,int nodeid,void *buf,int bufsize)
{
#ifndef GMX_MPI
  gmx_call("gmx_txs");
#else
  int tag;

#ifdef DEBUG
  fprintf(stderr,"gmx_txs: nodeid=%d, buf=%x, bufsize=%d\n",
	  nodeid,buf,bufsize);
#endif
  tag = 0;
  if (MPI_Send(buf,bufsize,MPI_BYTE,RANK(cr,nodeid),tag,cr->mpi_comm_mygroup) != 0)
    gmx_comm("MPI_Send Failed");
#endif
}

void gmx_rx(const t_commrec *cr,int nodeid,void *buf,int bufsize)
{
#ifndef GMX_MPI
  gmx_call("gmx_rx");
#else
  int        tag;

#ifdef DEBUG
  fprintf(stderr,"gmx_rx: nodeid=%d, buf=%x, bufsize=%d\n",
	  nodeid,buf,bufsize);
#endif
  tag = 0;
  if (MPI_Irecv( buf, bufsize, MPI_BYTE, RANK(cr,nodeid), tag, cr->mpi_comm_mygroup, &mpi_req_rx) != 0 )
    gmx_comm("MPI_Recv Failed");
#endif
}

void gmx_rx_wait(int nodeid)
{
#ifndef GMX_MPI
  gmx_call("gmx_rx_wait");
#else
  MPI_Status  status;
  int mpi_result;
  
  if ((mpi_result=MPI_Wait(&mpi_req_rx,&status)) != 0)
    gmx_fatal(FARGS,"MPI_Wait: result=%d",mpi_result);
#endif
}

int gmx_rx_probe(int nodeid)
{
#ifndef GMX_MPI
  gmx_call("gmx_rx_probe");
  return 0;
#else
  MPI_Status  status;
  int mpi_result,flag=0;
  
  if ((mpi_result = MPI_Test(&mpi_req_rx,&flag,&status)) != MPI_SUCCESS)
    gmx_fatal(FARGS,"MPI_Test: result=%d",mpi_result);
    
  return flag;
#endif
}

void gmx_rxs(const t_commrec *cr,int nodeid,void *buf,int bufsize)
{
#ifndef GMX_MPI
  gmx_call("gmx_rxs");
#else
  MPI_Status stat;
  int        tag;

#ifdef DEBUG
  fprintf(stderr,"gmx_rxs: nodeid=%d, buf=%x, bufsize=%d\n",
	  nodeid,buf,bufsize);
#endif
  tag = 0;
  if (MPI_Recv( buf, bufsize, MPI_BYTE, RANK(cr,nodeid), tag, cr->mpi_comm_mygroup, &stat) != 0 )
    gmx_comm("MPI_Recv Failed");
#endif
}

bool gmx_mpi_initialized(void)
{
  int n;
#ifndef GMX_MPI
  return 0;
#else
  MPI_Initialized(&n);
  
  return n;
#endif
}

int gmx_setup(int *argc,char **argv,int *nnodes)
{
#ifndef GMX_MPI
  gmx_call("gmx_setup");
  return 0;
#else
  char   buf[256];
  int    resultlen;               /* actual length of node name      */
  int    i,flag;
  int  mpi_num_nodes;
  int  mpi_my_rank;
  char mpi_hostname[MPI_MAX_PROCESSOR_NAME];

  /* Call the MPI routines */
  (void) MPI_Init(argc,&argv);
  (void) MPI_Comm_size( MPI_COMM_WORLD, &mpi_num_nodes );
  (void) MPI_Comm_rank( MPI_COMM_WORLD, &mpi_my_rank );
  (void) MPI_Get_processor_name( mpi_hostname, &resultlen );


#ifdef USE_MPE
  /* MPE logging routines */
  /* get event IDs from MPE: */  
  ev_pme_start               = MPE_Log_get_event_number( );
  ev_pme_finish              = MPE_Log_get_event_number( );
  ev_spread_on_grid_start    = MPE_Log_get_event_number( );
  ev_spread_on_grid_finish   = MPE_Log_get_event_number( );
  ev_sum_qgrid_start         = MPE_Log_get_event_number( );
  ev_sum_qgrid_finish        = MPE_Log_get_event_number( );
  ev_gmxfft3d_start          = MPE_Log_get_event_number( );
  ev_gmxfft3d_finish         = MPE_Log_get_event_number( );
  ev_solve_pme_start         = MPE_Log_get_event_number( );
  ev_solve_pme_finish        = MPE_Log_get_event_number( );
  ev_gather_f_bsplines_start = MPE_Log_get_event_number( );
  ev_gather_f_bsplines_finish= MPE_Log_get_event_number( );
  ev_timestep1               = MPE_Log_get_event_number( );
  ev_timestep2               = MPE_Log_get_event_number( );
  ev_force_start             = MPE_Log_get_event_number( );
  ev_force_finish            = MPE_Log_get_event_number( );
  ev_do_fnbf_start           = MPE_Log_get_event_number( );
  ev_do_fnbf_finish          = MPE_Log_get_event_number( );
  ev_reduce_start            = MPE_Log_get_event_number( );
  ev_reduce_finish           = MPE_Log_get_event_number( );
  ev_rscatter_start          = MPE_Log_get_event_number( );
  ev_rscatter_finish         = MPE_Log_get_event_number( );
  ev_alltoall_start          = MPE_Log_get_event_number( );
  ev_alltoall_finish         = MPE_Log_get_event_number( );
  ev_ns_start                = MPE_Log_get_event_number( );
  ev_ns_finish               = MPE_Log_get_event_number( );
  ev_calc_bonds_start        = MPE_Log_get_event_number( );
  ev_calc_bonds_finish       = MPE_Log_get_event_number( );
  ev_pmeredist_start         = MPE_Log_get_event_number( );
  ev_pmeredist_finish        = MPE_Log_get_event_number( );
  ev_init_pme_start          = MPE_Log_get_event_number( );      
  ev_init_pme_finish         = MPE_Log_get_event_number( );
  ev_global_stat_start       = MPE_Log_get_event_number( );
  ev_global_stat_finish      = MPE_Log_get_event_number( );
  ev_shift_start             = MPE_Log_get_event_number( );
  ev_shift_finish            = MPE_Log_get_event_number( );
  ev_send_coordinates_start  = MPE_Log_get_event_number( );
  ev_send_coordinates_finish = MPE_Log_get_event_number( );
  ev_update_fr_start         = MPE_Log_get_event_number( );
  ev_update_fr_finish        = MPE_Log_get_event_number( );
  ev_clear_rvecs_start       = MPE_Log_get_event_number( );
  ev_clear_rvecs_finish      = MPE_Log_get_event_number( ); 
  ev_update_start            = MPE_Log_get_event_number( ); 
  ev_update_finish           = MPE_Log_get_event_number( ); 
  ev_output_start            = MPE_Log_get_event_number( ); 
  ev_output_finish           = MPE_Log_get_event_number( ); 
  ev_sum_lrforces_start      = MPE_Log_get_event_number( ); 
  ev_sum_lrforces_finish     = MPE_Log_get_event_number( ); 
  ev_virial_start            = MPE_Log_get_event_number( );
  ev_virial_finish           = MPE_Log_get_event_number( );
  ev_sort_start              = MPE_Log_get_event_number( );
  ev_sort_finish             = MPE_Log_get_event_number( );
  ev_sum_qgrid_start         = MPE_Log_get_event_number( );
  ev_sum_qgrid_finish        = MPE_Log_get_event_number( );
  /* describe events: */
  if ( mpi_my_rank == 0 ) 
  {
    MPE_Describe_state(ev_pme_start,               ev_pme_finish,               "doing PME",       "grey" );
    MPE_Describe_state(ev_spread_on_grid_start,    ev_spread_on_grid_finish,    "spread",          "dark orange" );   
    MPE_Describe_state(ev_sum_qgrid_start,         ev_sum_qgrid_finish,         "sum qgrid",       "slate blue");
    MPE_Describe_state(ev_gmxfft3d_start,          ev_gmxfft3d_finish,          "fft3d",           "snow2" );   
    MPE_Describe_state(ev_solve_pme_start,         ev_solve_pme_finish,         "solve PME",       "indian red" );   
    MPE_Describe_state(ev_gather_f_bsplines_start, ev_gather_f_bsplines_finish, "bsplines",        "light sea green" );   
    MPE_Describe_state(ev_timestep1,               ev_timestep2,                "timestep START",  "magenta" );
    MPE_Describe_state(ev_force_start,             ev_force_finish,             "force",           "cornflower blue" );
    MPE_Describe_state(ev_do_fnbf_start,           ev_do_fnbf_finish,           "do_fnbf",         "navy" );
    MPE_Describe_state(ev_reduce_start,            ev_reduce_finish,            "reduce",          "cyan1" );
    MPE_Describe_state(ev_rscatter_start,          ev_rscatter_finish,          "rscatter",        "cyan3" );
    MPE_Describe_state(ev_alltoall_start,          ev_alltoall_finish,          "alltoall",        "LightCyan4" );
    MPE_Describe_state(ev_ns_start,                ev_ns_finish,                "neighbor search", "tomato" );
    MPE_Describe_state(ev_calc_bonds_start,        ev_calc_bonds_finish,        "bonded forces",   "slate blue" );
    MPE_Describe_state(ev_pmeredist_start,         ev_pmeredist_finish,         "pmeredist",       "thistle" );
    MPE_Describe_state(ev_init_pme_start,          ev_init_pme_finish,          "init PME",        "snow4");
    MPE_Describe_state(ev_global_stat_start,       ev_global_stat_finish,       "global stat",     "firebrick3");
    MPE_Describe_state(ev_shift_start,             ev_shift_finish,             "shift",           "orange");
    MPE_Describe_state(ev_send_coordinates_start,  ev_send_coordinates_finish,  "send_coordinates","blue");
    MPE_Describe_state(ev_update_fr_start,         ev_update_fr_finish,         "update forcerec", "goldenrod");
    MPE_Describe_state(ev_clear_rvecs_start,       ev_clear_rvecs_finish,       "clear rvecs",     "bisque");
    MPE_Describe_state(ev_update_start,            ev_update_finish,            "update",          "cornsilk");
    MPE_Describe_state(ev_output_start,            ev_output_finish,            "output",          "black");
    MPE_Describe_state(ev_sum_lrforces_start,      ev_sum_lrforces_finish,      "sum_LRforces",    "lime green");
    MPE_Describe_state(ev_virial_start,            ev_virial_finish,            "calc_virial",     "thistle4");
    MPE_Describe_state(ev_sort_start,              ev_sort_finish,              "sort pme atoms",  "brown");
    MPE_Describe_state(ev_sum_qgrid_start,         ev_sum_qgrid_finish,         "sum charge grid", "medium orchid");
    }
  MPE_Init_log();
#endif
  
  fprintf(stderr,"NNODES=%d, MYRANK=%d, HOSTNAME=%s\n",
	  mpi_num_nodes,mpi_my_rank,mpi_hostname);
  
  *nnodes=mpi_num_nodes;
  
  return mpi_my_rank;
#endif
}

int  gmx_node_num(void)
{
#ifndef GMX_MPI
  return 1;
#else
  int i;
  (void) MPI_Comm_size(MPI_COMM_WORLD, &i);
  return i;
#endif
}

int gmx_node_rank(void)
{
#ifndef GMX_MPI
  return 0;
#else
  int i;
  (void) MPI_Comm_rank(MPI_COMM_WORLD, &i);
  return i;
#endif
}

void gmx_setup_nodecomm(FILE *fplog,t_commrec *cr)
{
  gmx_nodecomm_t *nc;
  int  n,rank,resultlen,hostnum,i,ng;
#ifdef GMX_MPI
  char mpi_hostname[MPI_MAX_PROCESSOR_NAME];
#endif

  /* Many MPI implementations do not optimize MPI_Allreduce
   * (and probably also other global communication calls)
   * for multi-core nodes connected by a network.
   * We can optimize such communication by using one MPI call
   * within each node and one between the nodes.
   * For MVAPICH2 and Intel MPI this reduces the time for
   * the global_stat communication by 25%
   * for 2x2-core 3 GHz Woodcrest connected by mixed DDR/SDR Infiniband.
   * B. Hess, November 2007
   */

  nc = &cr->nc;

  nc->bUse = FALSE;
  if (getenv("GMX_NO_NODECOMM") == NULL) {
#ifdef GMX_MPI
    MPI_Comm_size(cr->mpi_comm_mygroup,&n);
    MPI_Comm_rank(cr->mpi_comm_mygroup,&rank);
    MPI_Get_processor_name(mpi_hostname,&resultlen);
    /* This procedure can only differentiate nodes with host names
     * that end on unique numbers.
     */
    i = resultlen - 1;
    while(i > 0 && isdigit(mpi_hostname[i-1]))
      i--;
    if (isdigit(mpi_hostname[i])) {
      hostnum = atoi(mpi_hostname+i);
    } else {
      hostnum = 0;
    }

    /* The intra-node communicator, split on node number */
    MPI_Comm_split(cr->mpi_comm_mygroup,hostnum,rank,&nc->comm_intra);
    MPI_Comm_rank(nc->comm_intra,&nc->rank_intra);
    /* The inter-node communicator, split on rank_intra.
     * We actually only need the one for rank=0,
     * but it is easier to create them all.
     */
    MPI_Comm_split(cr->mpi_comm_mygroup,nc->rank_intra,rank,&nc->comm_inter);
    /* Check if this really created two step communication */
    MPI_Comm_size(nc->comm_inter,&ng);
    if (ng > 1 && ng < n) {
      nc->bUse = TRUE;
      if (fplog)
	fprintf(fplog,"Using two step summing over %d groups of on average %.1f processes\n\n",ng,(real)n/(real)ng);
      if (nc->comm_intra > 0)
	MPI_Comm_free(&nc->comm_inter);
    } else {
      /* One group or all processes in a separate group, use normal summing */
      MPI_Comm_free(&nc->comm_inter);
      MPI_Comm_free(&nc->comm_intra);
    }
#endif
    
  }
}

int gmx_idle_send(void)
{
  return 0;
}

int gmx_idle_rec(void)
{
  return 0;
}

void gmx_left_right(int nnodes,int nodeid,int *left,int *right)
{
  *left  = (nnodes+nodeid-1) % nnodes;
  *right = (nodeid + 1) % nnodes;
}

void gmx_tx_rx(const t_commrec *cr,
	       int send_nodeid,void *send_buf,int send_bufsize,
	       int rec_nodeid,void *rec_buf,int rec_bufsize)
{
#ifndef GMX_MPI
  gmx_call("gmx_tx_rx");
#else
  int tx_tag = 0,rx_tag = 0;
  MPI_Status stat;
  
  MPI_Sendrecv(send_buf,send_bufsize,MPI_BYTE,RANK(cr,send_nodeid),tx_tag,
	       rec_buf,rec_bufsize,MPI_BYTE,RANK(cr,rec_nodeid),rx_tag,
	       cr->mpi_comm_mygroup,&stat);
#endif
}
		 
void gmx_tx_rx_real(const t_commrec *cr,
		    int send_nodeid,real *send_buf,int send_bufsize,
		    int rec_nodeid,real *rec_buf,int rec_bufsize)
{
#ifndef GMX_MPI
  gmx_call("gmx_tx_rx_real");
#else
  int tx_tag = 0,rx_tag = 0;
  MPI_Status stat;
#ifdef GMX_DOUBLE
#define mpi_type MPI_DOUBLE
#else
#define mpi_type MPI_FLOAT
#endif
  
  MPI_Sendrecv(send_buf,send_bufsize,mpi_type,RANK(cr,send_nodeid),tx_tag,
	       rec_buf,rec_bufsize,mpi_type,RANK(cr,rec_nodeid),rx_tag,
	       cr->mpi_comm_mygroup,&stat);
#undef mpi_type
#endif
}
		 
void gmx_wait(int left,int right)
{
#ifndef GMX_MPI
  gmx_call("gmx_wait");
#else
  gmx_tx_wait(left);
  gmx_rx_wait(right);
#endif
}

void gmx_sync_ring(const t_commrec *cr,int nodeid,int nnodes,int left,int right)
{
#ifndef GMX_MPI
  gmx_call("gmx_sync_ring");
#else
  int i;
  int tag=0;

  for (i=0; (i<nnodes); i++) {
    if (nodeid == 0) {
      gmx_txs(cr,right,&tag,sizeof(tag));
      gmx_rxs(cr,left,&tag,sizeof(tag));
    }
    else {
      gmx_rxs(cr,left,&tag,sizeof(tag));
      gmx_txs(cr,right,&tag,sizeof(tag));
    }
  }
#endif
}

void gmx_stat(FILE *fp,char *msg)
{
  fprintf(fp,"def_stat: %s (from %s, %d)\n",msg,__FILE__,__LINE__);
}

void gmx_reset_idle(void)
{
  ;
}

void gmx_abort(int noderank,int nnodes,int errorno)
{
#ifndef GMX_MPI
  gmx_call("gmx_abort");
#else
  if (nnodes > 1)
    fprintf(stderr,"Halting parallel program %s on CPU %d out of %d\n",
	    ShortProgram(),noderank,nnodes);
  else
    fprintf(stderr,"Halting program %s\n",ShortProgram());
  if (stdlog)
    fclose(stdlog);
  thanx(stderr);
  MPI_Abort(MPI_COMM_WORLD,errorno);
  exit(1);
#endif
}

void gmx_sumd(int nr,double r[],const t_commrec *cr)
{
#ifndef GMX_MPI
  gmx_call("gmx_sumd");
#else
#define GMX_MPI_SUM
#ifdef GMX_MPI_SUM
  static double *buf=NULL;
  static int nalloc=0;
  int i;
  
  if (nr > nalloc) {
    nalloc = nr;
    srenew(buf,nalloc);
  }
  if (cr->nc.bUse) {
    /* Use two step summing */
    MPI_Allreduce(r,buf,nr,MPI_DOUBLE,MPI_SUM,cr->nc.comm_intra);
    if (cr->nc.rank_intra == 0) {
      /* Sum with the buffers reversed */
      MPI_Allreduce(buf,r,nr,MPI_DOUBLE,MPI_SUM,cr->nc.comm_inter);
    }
    MPI_Bcast(r,nr,MPI_DOUBLE,0,cr->nc.comm_intra);
  } else {
    MPI_Allreduce(r,buf,nr,MPI_DOUBLE,MPI_SUM,cr->mpi_comm_mygroup);
    for(i=0; i<nr; i++)
      r[i] = buf[i];
  }
#else
  double *buf[2];
  int  NR,bufs,j,i,cur=0;
#define next (1-cur)

  bufs=nr*sizeof(buf[0][0]);
  NR=nr;

  snew(buf[0],NR);
  snew(buf[1],NR);

  for(i=0; (i<nr); i++)
    buf[cur][i]=r[i];
  for(j=0; (j<(cr->nnodes-cr->npmenodes-1)); j++) {
    gmx_tx(cr->left,buf[cur],bufs);
    gmx_rx(cr->right,buf[next],bufs);
    gmx_tx_wait(cr->left);
    gmx_rx_wait(cr->right);
    for(i=0; (i<nr); i++)
      r[i]+=buf[next][i];

    cur=next;
  }
  sfree(buf[1]);
  sfree(buf[0]);
#endif
#endif
}

void gmx_sumf(int nr,float r[],const t_commrec *cr)
{
#ifndef GMX_MPI
  gmx_call("gmx_sumf");
#else
#ifdef GMX_MPI_SUM
  static float *buf=NULL;
  static int nalloc=0;
  int i;
  
  if (nr > nalloc) {
    nalloc = nr;
    srenew(buf,nalloc);
  }
  if (cr->nc.bUse) {
    /* Use two step summing */
    MPI_Allreduce(r,buf,nr,MPI_FLOAT,MPI_SUM,cr->nc.comm_intra);
    if (cr->nc.rank_intra == 0) {
      /* Sum with the buffers reversed */
      MPI_Allreduce(buf,r,nr,MPI_FLOAT,MPI_SUM,cr->nc.comm_inter);
    }
    MPI_Bcast(r,nr,MPI_FLOAT,0,cr->nc.comm_intra);
  } else {
    MPI_Allreduce(r,buf,nr,MPI_FLOAT,MPI_SUM,cr->mpi_comm_mygroup);
    for(i=0; i<nr; i++)
      r[i] = buf[i];
  }
#else
  float *buf[2];
  int  NR,bufs,j,i,cur=0;
#define next (1-cur)

  bufs=nr*sizeof(float);
  NR=nr;

  snew(buf[0],NR);
  snew(buf[1],NR);

  for(i=0; (i<nr); i++)
    buf[cur][i]=r[i];
  for(j=0; (j<(cr->nnodes-cr->npmenodes-1)); j++) {
    gmx_tx(cr,cr->left,buf[cur],bufs);
    gmx_rx(cr,cr->right,buf[next],bufs);
    gmx_wait(cr->left,cr->right);
    for(i=0; (i<nr); i++)
      r[i]+=buf[next][i];

    cur=next;
  }
  sfree(buf[1]);
  sfree(buf[0]);
#endif
#endif
}

void gmx_sumi(int nr,int r[],const t_commrec *cr)
{
#ifndef GMX_MPI
  gmx_call("gmx_sumi");
#else
#define GMX_MPI_SUM
#ifdef GMX_MPI_SUM
  static int *buf=NULL;
  static int nalloc=0;
  int i;
  
  if (nr > nalloc) {
    nalloc = nr;
    srenew(buf,nalloc);
  }
  if (cr->nc.bUse) {
    /* Use two step summing */
    MPI_Allreduce(r,buf,nr,MPI_INT,MPI_SUM,cr->nc.comm_intra);
    if (cr->nc.rank_intra == 0) {
      /* Sum with the buffers reversed */
      MPI_Allreduce(buf,r,nr,MPI_INT,MPI_SUM,cr->nc.comm_inter);
    }
    MPI_Bcast(r,nr,MPI_INT,0,cr->nc.comm_intra);
  } else {
    MPI_Allreduce(r,buf,nr,MPI_INT,MPI_SUM,cr->mpi_comm_mygroup);
    for(i=0; i<nr; i++)
      r[i] = buf[i];
  }
#else
  int *buf[2];
  int  NR,bufs,j,i,cur=0;
#define next (1-cur)

  bufs=nr*sizeof(int);
  NR=nr;

  snew(buf[0],NR);
  snew(buf[1],NR);

  for(i=0; (i<nr); i++)
    buf[cur][i]=r[i];
  for(j=0; (j<(cr->nnodes-cr->npmenodes-1)); j++) {
    gmx_tx(cr,cr->left,buf[cur],bufs);
    gmx_rx(cr,cr->right,buf[next],bufs);
    gmx_wait(cr->left,cr->right);
    for(i=0; (i<nr); i++)
      r[i]+=buf[next][i];

    cur=next;
  }
  sfree(buf[1]);
  sfree(buf[0]);
#endif
#endif
}

void gmx_sumd_sim(int nr,double r[],const gmx_multisim_t *ms)
{
#ifndef GMX_MPI
  gmx_call("gmx_sumd");
#else
  static double *buf=NULL;
  static int nalloc=0;
  int i;
  
  if (nr > nalloc) {
    nalloc = nr;
    srenew(buf,nalloc);
  }
  MPI_Allreduce(r,buf,nr,MPI_DOUBLE,MPI_SUM,ms->mpi_comm_masters);
  for(i=0; i<nr; i++)
    r[i] = buf[i];
#endif
}

void gmx_sumf_sim(int nr,float r[],const gmx_multisim_t *ms)
{
#ifndef GMX_MPI
  gmx_call("gmx_sumd");
#else
  static float *buf=NULL;
  static int nalloc=0;
  int i;
  
  if (nr > nalloc) {
    nalloc = nr;
    srenew(buf,nalloc);
  }
  MPI_Allreduce(r,buf,nr,MPI_FLOAT,MPI_SUM,ms->mpi_comm_masters);
  for(i=0; i<nr; i++)
    r[i] = buf[i];
#endif
}

void gmx_sumi_sim(int nr,int r[],const gmx_multisim_t *ms)
{
#ifndef GMX_MPI
  gmx_call("gmx_sumd");
#else
  static int *buf=NULL;
  static int nalloc=0;
  int i;
  
  if (nr > nalloc) {
    nalloc = nr;
    srenew(buf,nalloc);
  }
  MPI_Allreduce(r,buf,nr,MPI_INT,MPI_SUM,ms->mpi_comm_masters);
  for(i=0; i<nr; i++)
    r[i] = buf[i];
#endif
}

void gmx_finalize(const t_commrec *cr)
{
  int ret;
#ifndef GMX_MPI
  gmx_call("gmx_finalize");
#else
  /* We sync the processes here to try to avoid problems
   * with buggy MPI implementations that could cause
   * unfinished processes to terminate.
   */
  MPI_Barrier(MPI_COMM_WORLD);

  /*
  if (DOMAINDECOMP(cr)) {
    if (cr->npmenodes > 0 || cr->dd->bCartesian) 
      MPI_Comm_free(&cr->mpi_comm_mygroup);
    if (cr->dd->bCartesian)
      MPI_Comm_free(&cr->mpi_comm_mysim);
  }
  */

  /* Apparently certain mpich implementations cause problems
   * with MPI_Finalize. In that case comment out MPI_Finalize.
   */
  if (debug)
    fprintf(debug,"Will call MPI_Finalize now\n");
  ret = MPI_Finalize();
  if (debug)
    fprintf(debug,"Return code from MPI_Finalize = %d\n",ret);
#endif
}

