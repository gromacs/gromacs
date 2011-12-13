/*
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
#include "macros.h"

#ifdef GMX_LIB_MPI
#include <mpi.h>
#endif

#ifdef GMX_THREADS
#include "tmpi.h"
#endif

#include "mpelogging.h"

/* The source code in this file should be thread-safe. 
      Please keep it that way. */

gmx_bool gmx_mpi_initialized(void)
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
#ifdef GMX_LIB_MPI
#ifdef GMX_FAHCORE
  (void) fah_MPI_Init(argc,&argv);
#else
  (void) MPI_Init(argc,&argv);
#endif
#endif
  (void) MPI_Comm_size( MPI_COMM_WORLD, &mpi_num_nodes );
  (void) MPI_Comm_rank( MPI_COMM_WORLD, &mpi_my_rank );
  (void) MPI_Get_processor_name( mpi_hostname, &resultlen );


#ifdef USE_MPE
  /* MPE logging routines. Get event IDs from MPE: */
  /* General events */
  ev_timestep1               = MPE_Log_get_event_number( );
  ev_timestep2               = MPE_Log_get_event_number( );
  ev_force_start             = MPE_Log_get_event_number( );
  ev_force_finish            = MPE_Log_get_event_number( );
  ev_do_fnbf_start           = MPE_Log_get_event_number( );
  ev_do_fnbf_finish          = MPE_Log_get_event_number( );
  ev_ns_start                = MPE_Log_get_event_number( );
  ev_ns_finish               = MPE_Log_get_event_number( );
  ev_calc_bonds_start        = MPE_Log_get_event_number( );
  ev_calc_bonds_finish       = MPE_Log_get_event_number( );
  ev_global_stat_start       = MPE_Log_get_event_number( );
  ev_global_stat_finish      = MPE_Log_get_event_number( );
  ev_virial_start            = MPE_Log_get_event_number( );
  ev_virial_finish           = MPE_Log_get_event_number( );
  
  /* Shift related events */
  ev_shift_start             = MPE_Log_get_event_number( );
  ev_shift_finish            = MPE_Log_get_event_number( );
  ev_unshift_start           = MPE_Log_get_event_number( );
  ev_unshift_finish          = MPE_Log_get_event_number( );
  ev_mk_mshift_start         = MPE_Log_get_event_number( );
  ev_mk_mshift_finish        = MPE_Log_get_event_number( );
  
  /* PME related events */
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
  ev_reduce_start            = MPE_Log_get_event_number( );
  ev_reduce_finish           = MPE_Log_get_event_number( );
  ev_rscatter_start          = MPE_Log_get_event_number( );
  ev_rscatter_finish         = MPE_Log_get_event_number( );
  ev_alltoall_start          = MPE_Log_get_event_number( );
  ev_alltoall_finish         = MPE_Log_get_event_number( );
  ev_pmeredist_start         = MPE_Log_get_event_number( );
  ev_pmeredist_finish        = MPE_Log_get_event_number( );
  ev_init_pme_start          = MPE_Log_get_event_number( );      
  ev_init_pme_finish         = MPE_Log_get_event_number( );
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
  ev_sort_start              = MPE_Log_get_event_number( );
  ev_sort_finish             = MPE_Log_get_event_number( );
  ev_sum_qgrid_start         = MPE_Log_get_event_number( );
  ev_sum_qgrid_finish        = MPE_Log_get_event_number( );
  
  /* Essential dynamics related events */
  ev_edsam_start             = MPE_Log_get_event_number( );
  ev_edsam_finish            = MPE_Log_get_event_number( );
  ev_get_coords_start        = MPE_Log_get_event_number( );
  ev_get_coords_finish       = MPE_Log_get_event_number( );
  ev_ed_apply_cons_start     = MPE_Log_get_event_number( );
  ev_ed_apply_cons_finish    = MPE_Log_get_event_number( );
  ev_fit_to_reference_start  = MPE_Log_get_event_number( );
  ev_fit_to_reference_finish = MPE_Log_get_event_number( );
  
  /* describe events: */
  if ( mpi_my_rank == 0 ) 
  {
    /* General events */
    MPE_Describe_state(ev_timestep1,               ev_timestep2,                "timestep START",  "magenta" );
    MPE_Describe_state(ev_force_start,             ev_force_finish,             "force",           "cornflower blue" );
    MPE_Describe_state(ev_do_fnbf_start,           ev_do_fnbf_finish,           "do_fnbf",         "navy" );
    MPE_Describe_state(ev_ns_start,                ev_ns_finish,                "neighbor search", "tomato" );
    MPE_Describe_state(ev_calc_bonds_start,        ev_calc_bonds_finish,        "bonded forces",   "slate blue" );
    MPE_Describe_state(ev_global_stat_start,       ev_global_stat_finish,       "global stat",     "firebrick3");
    MPE_Describe_state(ev_update_fr_start,         ev_update_fr_finish,         "update forcerec", "goldenrod");
    MPE_Describe_state(ev_clear_rvecs_start,       ev_clear_rvecs_finish,       "clear rvecs",     "bisque");
    MPE_Describe_state(ev_update_start,            ev_update_finish,            "update",          "cornsilk");
    MPE_Describe_state(ev_output_start,            ev_output_finish,            "output",          "black");
    MPE_Describe_state(ev_virial_start,            ev_virial_finish,            "calc_virial",     "thistle4");
    
    /* PME related events */
    MPE_Describe_state(ev_pme_start,               ev_pme_finish,               "doing PME",       "grey" );
    MPE_Describe_state(ev_spread_on_grid_start,    ev_spread_on_grid_finish,    "spread",          "dark orange" );   
    MPE_Describe_state(ev_sum_qgrid_start,         ev_sum_qgrid_finish,         "sum qgrid",       "slate blue");
    MPE_Describe_state(ev_gmxfft3d_start,          ev_gmxfft3d_finish,          "fft3d",           "snow2" );   
    MPE_Describe_state(ev_solve_pme_start,         ev_solve_pme_finish,         "solve PME",       "indian red" );   
    MPE_Describe_state(ev_gather_f_bsplines_start, ev_gather_f_bsplines_finish, "bsplines",        "light sea green" );   
    MPE_Describe_state(ev_reduce_start,            ev_reduce_finish,            "reduce",          "cyan1" );
    MPE_Describe_state(ev_rscatter_start,          ev_rscatter_finish,          "rscatter",        "cyan3" );
    MPE_Describe_state(ev_alltoall_start,          ev_alltoall_finish,          "alltoall",        "LightCyan4" );
    MPE_Describe_state(ev_pmeredist_start,         ev_pmeredist_finish,         "pmeredist",       "thistle" );
    MPE_Describe_state(ev_init_pme_start,          ev_init_pme_finish,          "init PME",        "snow4");
    MPE_Describe_state(ev_send_coordinates_start,  ev_send_coordinates_finish,  "send_coordinates","blue");
    MPE_Describe_state(ev_sum_lrforces_start,      ev_sum_lrforces_finish,      "sum_LRforces",    "lime green");
    MPE_Describe_state(ev_sort_start,              ev_sort_finish,              "sort pme atoms",  "brown");
    MPE_Describe_state(ev_sum_qgrid_start,         ev_sum_qgrid_finish,         "sum charge grid", "medium orchid");
    
    /* Shift related events */
    MPE_Describe_state(ev_shift_start,             ev_shift_finish,             "shift",           "orange");
    MPE_Describe_state(ev_unshift_start,           ev_unshift_finish,           "unshift",         "dark orange");    
    MPE_Describe_state(ev_mk_mshift_start,         ev_mk_mshift_finish,         "mk_mshift",       "maroon");
        
    /* Essential dynamics related events */
    MPE_Describe_state(ev_edsam_start,             ev_edsam_finish,             "EDSAM",           "deep sky blue");
    MPE_Describe_state(ev_get_coords_start,        ev_get_coords_finish,        "ED get coords",   "steel blue");
    MPE_Describe_state(ev_ed_apply_cons_start,     ev_ed_apply_cons_finish,     "ED apply constr", "forest green");
    MPE_Describe_state(ev_fit_to_reference_start,  ev_fit_to_reference_finish,  "ED fit to ref",   "lavender");
       
  }
  MPE_Init_log();
#endif
 
#ifdef GMX_LIB_MPI 
  fprintf(stderr,"NNODES=%d, MYRANK=%d, HOSTNAME=%s\n",
	  mpi_num_nodes,mpi_my_rank,mpi_hostname);
#endif
  
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
  int  n,rank,resultlen,hostnum,i,j,ng,ni;
#ifdef GMX_MPI
  char mpi_hostname[MPI_MAX_PROCESSOR_NAME],num[MPI_MAX_PROCESSOR_NAME];
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
#ifndef GMX_THREADS
  if (getenv("GMX_NO_NODECOMM") == NULL) {
#ifdef GMX_MPI
    MPI_Comm_size(cr->mpi_comm_mygroup,&n);
    MPI_Comm_rank(cr->mpi_comm_mygroup,&rank);
    MPI_Get_processor_name(mpi_hostname,&resultlen);
    /* This procedure can only differentiate nodes with host names
     * that end on unique numbers.
     */
    i = 0;
    j = 0;
    /* Only parse the host name up to the first dot */
    while(i < resultlen && mpi_hostname[i] != '.') {
      if (isdigit(mpi_hostname[i])) {
	num[j++] = mpi_hostname[i];
      }
      i++;
    }
    num[j] = '\0';
    if (j == 0) {
      hostnum = 0;
    } else {
      /* Use only the last 9 decimals, so we don't overflow an int */
      hostnum = strtol(num + max(0,j-9), NULL, 10); 
    }

    if (debug) {
      fprintf(debug,
	      "In gmx_setup_nodecomm: splitting communicator of size %d\n",
	      n);
      fprintf(debug,"In gmx_setup_nodecomm: hostname '%s', hostnum %d\n",
	      mpi_hostname,hostnum);
    }

    /* The intra-node communicator, split on node number */
    MPI_Comm_split(cr->mpi_comm_mygroup,hostnum,rank,&nc->comm_intra);
    MPI_Comm_rank(nc->comm_intra,&nc->rank_intra);
    if (debug) {
      fprintf(debug,"In gmx_setup_nodecomm: node rank %d rank_intra %d\n",
	      rank,nc->rank_intra);
    }
    /* The inter-node communicator, split on rank_intra.
     * We actually only need the one for rank=0,
     * but it is easier to create them all.
     */
    MPI_Comm_split(cr->mpi_comm_mygroup,nc->rank_intra,rank,&nc->comm_inter);
    /* Check if this really created two step communication */
    MPI_Comm_size(nc->comm_inter,&ng);
    MPI_Comm_size(nc->comm_intra,&ni);
    if (debug) {
      fprintf(debug,"In gmx_setup_nodecomm: groups %d, my group size %d\n",
	      ng,ni);
    }
    if ((ng > 1 && ng < n) || (ni > 1 && ni < n)) {
      nc->bUse = TRUE;
      if (fplog)
	fprintf(fplog,"Using two step summing over %d groups of on average %.1f processes\n\n",ng,(real)n/(real)ng);
      if (nc->rank_intra > 0)
	MPI_Comm_free(&nc->comm_inter);
    } else {
      /* One group or all processes in a separate group, use normal summing */
      MPI_Comm_free(&nc->comm_inter);
      MPI_Comm_free(&nc->comm_intra);
    }
#endif
  }
#endif
}

void gmx_barrier(const t_commrec *cr)
{
#ifndef GMX_MPI
  gmx_call("gmx_barrier");
#else
  MPI_Barrier(cr->mpi_comm_mygroup);
#endif
}

void gmx_abort(int noderank,int nnodes,int errorno)
{
#ifndef GMX_MPI
  gmx_call("gmx_abort");
#else
#ifdef GMX_THREADS
  fprintf(stderr,"Halting program %s\n",ShortProgram());
  thanx(stderr);
  exit(1);
#else
  if (nnodes > 1)
  {
      fprintf(stderr,"Halting parallel program %s on CPU %d out of %d\n",
              ShortProgram(),noderank,nnodes);
  }
  else
  {
      fprintf(stderr,"Halting program %s\n",ShortProgram());
  }

  thanx(stderr);
  MPI_Abort(MPI_COMM_WORLD,errorno);
  exit(1);
#endif
#endif
}

void gmx_bcast(int nbytes,void *b,const t_commrec *cr)
{
#ifndef GMX_MPI
  gmx_call("gmx_bast");
#else
  MPI_Bcast(b,nbytes,MPI_BYTE,MASTERRANK(cr),cr->mpi_comm_mygroup);
#endif
}

void gmx_bcast_sim(int nbytes,void *b,const t_commrec *cr)
{
#ifndef GMX_MPI
  gmx_call("gmx_bast");
#else
  MPI_Bcast(b,nbytes,MPI_BYTE,MASTERRANK(cr),cr->mpi_comm_mysim);
#endif
}

void gmx_sumd(int nr,double r[],const t_commrec *cr)
{
#ifndef GMX_MPI
    gmx_call("gmx_sumd");
#else
#if defined(MPI_IN_PLACE_EXISTS) || defined(GMX_THREADS)
    if (cr->nc.bUse) {
        if (cr->nc.rank_intra == 0)
        {
            /* Use two step summing. */
            MPI_Reduce(MPI_IN_PLACE,r,nr,MPI_DOUBLE,MPI_SUM,0,
                       cr->nc.comm_intra);
            /* Sum the roots of the internal (intra) buffers. */
            MPI_Allreduce(MPI_IN_PLACE,r,nr,MPI_DOUBLE,MPI_SUM,
                          cr->nc.comm_inter);
        }
        else
        {
            /* This is here because of the silly MPI specification
                that MPI_IN_PLACE should be put in sendbuf instead of recvbuf */
            MPI_Reduce(r,NULL,nr,MPI_DOUBLE,MPI_SUM,0,cr->nc.comm_intra);
        }
        MPI_Bcast(r,nr,MPI_DOUBLE,0,cr->nc.comm_intra);
    } 
    else 
    {
        MPI_Allreduce(MPI_IN_PLACE,r,nr,MPI_DOUBLE,MPI_SUM, 
                      cr->mpi_comm_mygroup);
    }
#else
    int i;

    if (nr > cr->mpb->dbuf_alloc) {
        cr->mpb->dbuf_alloc = nr;
        srenew(cr->mpb->dbuf,cr->mpb->dbuf_alloc);
    }
    if (cr->nc.bUse) {
        /* Use two step summing */
        MPI_Allreduce(r,cr->mpb->dbuf,nr,MPI_DOUBLE,MPI_SUM,cr->nc.comm_intra);
        if (cr->nc.rank_intra == 0) {
            /* Sum with the buffers reversed */
            MPI_Allreduce(cr->mpb->dbuf,r,nr,MPI_DOUBLE,MPI_SUM, 
                          cr->nc.comm_inter);
        }
        MPI_Bcast(r,nr,MPI_DOUBLE,0,cr->nc.comm_intra);
    } else {
        MPI_Allreduce(r,cr->mpb->dbuf,nr,MPI_DOUBLE,MPI_SUM,
                      cr->mpi_comm_mygroup);
        for(i=0; i<nr; i++)
            r[i] = cr->mpb->dbuf[i];
    }
#endif
#endif
}

void gmx_sumf(int nr,float r[],const t_commrec *cr)
{
#ifndef GMX_MPI
    gmx_call("gmx_sumf");
#else
#if defined(MPI_IN_PLACE_EXISTS) || defined(GMX_THREADS)
    if (cr->nc.bUse) {
        /* Use two step summing.  */
        if (cr->nc.rank_intra == 0)
        {
            MPI_Reduce(MPI_IN_PLACE,r,nr,MPI_FLOAT,MPI_SUM,0,
                       cr->nc.comm_intra);
            /* Sum the roots of the internal (intra) buffers */
            MPI_Allreduce(MPI_IN_PLACE,r,nr,MPI_FLOAT,MPI_SUM,
                          cr->nc.comm_inter);
        }
        else
        {
            /* This is here because of the silly MPI specification
                that MPI_IN_PLACE should be put in sendbuf instead of recvbuf */
            MPI_Reduce(r,NULL,nr,MPI_FLOAT,MPI_SUM,0,cr->nc.comm_intra);
        }
        MPI_Bcast(r,nr,MPI_FLOAT,0,cr->nc.comm_intra);
    } 
    else 
    {
        MPI_Allreduce(MPI_IN_PLACE,r,nr,MPI_FLOAT,MPI_SUM,cr->mpi_comm_mygroup);
    }
#else
    int i;

    if (nr > cr->mpb->fbuf_alloc) {
        cr->mpb->fbuf_alloc = nr;
        srenew(cr->mpb->fbuf,cr->mpb->fbuf_alloc);
    }
    if (cr->nc.bUse) {
        /* Use two step summing */
        MPI_Allreduce(r,cr->mpb->fbuf,nr,MPI_FLOAT,MPI_SUM,cr->nc.comm_intra);
        if (cr->nc.rank_intra == 0) {
            /* Sum with the buffers reversed */
            MPI_Allreduce(cr->mpb->fbuf,r,nr,MPI_FLOAT,MPI_SUM, 
                          cr->nc.comm_inter);
        }
        MPI_Bcast(r,nr,MPI_FLOAT,0,cr->nc.comm_intra);
    } else {
        MPI_Allreduce(r,cr->mpb->fbuf,nr,MPI_FLOAT,MPI_SUM,
                      cr->mpi_comm_mygroup);
        for(i=0; i<nr; i++)
            r[i] = cr->mpb->fbuf[i];
    }
#endif
#endif
}

void gmx_sumi(int nr,int r[],const t_commrec *cr)
{
#ifndef GMX_MPI
    gmx_call("gmx_sumi");
#else
#if defined(MPI_IN_PLACE_EXISTS) || defined(GMX_THREADS)
    if (cr->nc.bUse) {
        /* Use two step summing */
        if (cr->nc.rank_intra == 0) 
        {
            MPI_Reduce(MPI_IN_PLACE,r,nr,MPI_INT,MPI_SUM,0,cr->nc.comm_intra);
            /* Sum with the buffers reversed */
            MPI_Allreduce(MPI_IN_PLACE,r,nr,MPI_INT,MPI_SUM,cr->nc.comm_inter);
        }
        else
        {
            /* This is here because of the silly MPI specification
                that MPI_IN_PLACE should be put in sendbuf instead of recvbuf */
            MPI_Reduce(r,NULL,nr,MPI_INT,MPI_SUM,0,cr->nc.comm_intra);
        }
        MPI_Bcast(r,nr,MPI_INT,0,cr->nc.comm_intra);
    } 
    else 
    {
        MPI_Allreduce(MPI_IN_PLACE,r,nr,MPI_INT,MPI_SUM,cr->mpi_comm_mygroup);
    }
#else
    int i;

    if (nr > cr->mpb->ibuf_alloc) {
        cr->mpb->ibuf_alloc = nr;
        srenew(cr->mpb->ibuf,cr->mpb->ibuf_alloc);
    }
    if (cr->nc.bUse) {
        /* Use two step summing */
        MPI_Allreduce(r,cr->mpb->ibuf,nr,MPI_INT,MPI_SUM,cr->nc.comm_intra);
        if (cr->nc.rank_intra == 0) {
            /* Sum with the buffers reversed */
            MPI_Allreduce(cr->mpb->ibuf,r,nr,MPI_INT,MPI_SUM,cr->nc.comm_inter);
        }
        MPI_Bcast(r,nr,MPI_INT,0,cr->nc.comm_intra);
    } else {
        MPI_Allreduce(r,cr->mpb->ibuf,nr,MPI_INT,MPI_SUM,cr->mpi_comm_mygroup);
        for(i=0; i<nr; i++)
            r[i] = cr->mpb->ibuf[i];
    }
#endif
#endif
}

void gmx_sumli(int nr,gmx_large_int_t r[],const t_commrec *cr)
{
#ifndef GMX_MPI
    gmx_call("gmx_sumli");
#else
#if defined(MPI_IN_PLACE_EXISTS) || defined(GMX_THREADS)
    if (cr->nc.bUse) {
        /* Use two step summing */
        if (cr->nc.rank_intra == 0) 
        {
            MPI_Reduce(MPI_IN_PLACE,r,nr,GMX_MPI_LARGE_INT,MPI_SUM,0,
                       cr->nc.comm_intra);
            /* Sum with the buffers reversed */
            MPI_Allreduce(MPI_IN_PLACE,r,nr,GMX_MPI_LARGE_INT,MPI_SUM,
                          cr->nc.comm_inter);
        }
        else
        {
            /* This is here because of the silly MPI specification
                that MPI_IN_PLACE should be put in sendbuf instead of recvbuf */
            MPI_Reduce(r,NULL,nr,GMX_MPI_LARGE_INT,MPI_SUM,0,cr->nc.comm_intra);
        }
        MPI_Bcast(r,nr,GMX_MPI_LARGE_INT,0,cr->nc.comm_intra);
    } 
    else 
    {
        MPI_Allreduce(MPI_IN_PLACE,r,nr,GMX_MPI_LARGE_INT,MPI_SUM,cr->mpi_comm_mygroup);
    }
#else
    int i;

    if (nr > cr->mpb->libuf_alloc) {
        cr->mpb->libuf_alloc = nr;
        srenew(cr->mpb->libuf,cr->mpb->libuf_alloc);
    }
    if (cr->nc.bUse) {
        /* Use two step summing */
        MPI_Allreduce(r,cr->mpb->libuf,nr,GMX_MPI_LARGE_INT,MPI_SUM,
                      cr->nc.comm_intra);
        if (cr->nc.rank_intra == 0) {
            /* Sum with the buffers reversed */
            MPI_Allreduce(cr->mpb->libuf,r,nr,GMX_MPI_LARGE_INT,MPI_SUM,
                          cr->nc.comm_inter);
        }
        MPI_Bcast(r,nr,GMX_MPI_LARGE_INT,0,cr->nc.comm_intra);
    } else {
        MPI_Allreduce(r,cr->mpb->libuf,nr,GMX_MPI_LARGE_INT,MPI_SUM,
                      cr->mpi_comm_mygroup);
        for(i=0; i<nr; i++)
            r[i] = cr->mpb->libuf[i];
    }
#endif
#endif
}



#ifdef GMX_MPI
void gmx_sumd_comm(int nr,double r[],MPI_Comm mpi_comm)
{
#if defined(MPI_IN_PLACE_EXISTS) || defined(GMX_THREADS)
    MPI_Allreduce(MPI_IN_PLACE,r,nr,MPI_DOUBLE,MPI_SUM,mpi_comm);
#else
    /* this function is only used in code that is not performance critical,
       (during setup, when comm_rec is not the appropriate communication  
       structure), so this isn't as bad as it looks. */
    double *buf;
    int i;

    snew(buf, nr);
    MPI_Allreduce(r,buf,nr,MPI_DOUBLE,MPI_SUM,mpi_comm);
    for(i=0; i<nr; i++)
        r[i] = buf[i];
    sfree(buf);
#endif
}
#endif

#ifdef GMX_MPI
void gmx_sumf_comm(int nr,float r[],MPI_Comm mpi_comm)
{
#if defined(MPI_IN_PLACE_EXISTS) || defined(GMX_THREADS)
    MPI_Allreduce(MPI_IN_PLACE,r,nr,MPI_FLOAT,MPI_SUM,mpi_comm);
#else
    /* this function is only used in code that is not performance critical,
       (during setup, when comm_rec is not the appropriate communication  
       structure), so this isn't as bad as it looks. */
    float *buf;
    int i;

    snew(buf, nr);
    MPI_Allreduce(r,buf,nr,MPI_FLOAT,MPI_SUM,mpi_comm);
    for(i=0; i<nr; i++)
        r[i] = buf[i];
    sfree(buf);
#endif
}
#endif

void gmx_sumd_sim(int nr,double r[],const gmx_multisim_t *ms)
{
#ifndef GMX_MPI
  gmx_call("gmx_sumd_sim");
#else
  gmx_sumd_comm(nr,r,ms->mpi_comm_masters);
#endif
}

void gmx_sumf_sim(int nr,float r[],const gmx_multisim_t *ms)
{
#ifndef GMX_MPI
  gmx_call("gmx_sumf_sim");
#else
  gmx_sumf_comm(nr,r,ms->mpi_comm_masters);
#endif
}

void gmx_sumi_sim(int nr,int r[], const gmx_multisim_t *ms)
{
#ifndef GMX_MPI
    gmx_call("gmx_sumi_sim");
#else
#if defined(MPI_IN_PLACE_EXISTS) || defined(GMX_THREADS)
    MPI_Allreduce(MPI_IN_PLACE,r,nr,MPI_INT,MPI_SUM,ms->mpi_comm_masters);
#else
    /* this is thread-unsafe, but it will do for now: */
    int i;

    if (nr > ms->mpb->ibuf_alloc) {
        ms->mpb->ibuf_alloc = nr;
        srenew(ms->mpb->ibuf,ms->mpb->ibuf_alloc);
    }
    MPI_Allreduce(r,ms->mpb->ibuf,nr,MPI_INT,MPI_SUM,ms->mpi_comm_masters);
    for(i=0; i<nr; i++)
        r[i] = ms->mpb->ibuf[i];
#endif
#endif
}

void gmx_sumli_sim(int nr,gmx_large_int_t r[], const gmx_multisim_t *ms)
{
#ifndef GMX_MPI
    gmx_call("gmx_sumli_sim");
#else
#if defined(MPI_IN_PLACE_EXISTS) || defined(GMX_THREADS)
    MPI_Allreduce(MPI_IN_PLACE,r,nr,GMX_MPI_LARGE_INT,MPI_SUM,
                  ms->mpi_comm_masters);
#else
    /* this is thread-unsafe, but it will do for now: */
    int i;

    if (nr > ms->mpb->libuf_alloc) {
        ms->mpb->libuf_alloc = nr;
        srenew(ms->mpb->libuf,ms->mpb->libuf_alloc);
    }
    MPI_Allreduce(r,ms->mpb->libuf,nr,GMX_MPI_LARGE_INT,MPI_SUM,
                  ms->mpi_comm_masters);
    for(i=0; i<nr; i++)
        r[i] = ms->mpb->libuf[i];
#endif
#endif
}


void gmx_finalize(void)
{
#ifndef GMX_MPI
  gmx_call("gmx_finalize");
#else
  int ret;

  /* just as a check; we don't want to finalize twice */
  int finalized;
  MPI_Finalized(&finalized);
  if (finalized)
      return;

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

