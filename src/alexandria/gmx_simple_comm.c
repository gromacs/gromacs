/* -*- mode: c; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; c-file-style: "stroustrup"; -*-
 * $Id: tune_dip.c,v 1.39 2009/04/12 21:24:26 spoel Exp $
 * 
 *                This source code is part of
 * 
 *                 G   R   O   M   A   C   S
 * 
 *          GROningen MAchine for Chemical Simulations
 * 
 *                        VERSION 4.0.99
 * Written by David van der Spoel, Erik Lindahl, Berk Hess, and others.
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2008, The GROMACS development team,
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
 * Groningen Machine for Chemical Simulation
 */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <string.h>
#ifdef GMX_LIB_MPI
#include <mpi.h>
#endif

#ifdef GMX_THREADS
#include "tmpi.h"
#endif

#include "typedefs.h"
#include "smalloc.h"
#include "network.h"
#include "gmx_simple_comm.h"
	
void gmx_send(const t_commrec *cr,int dest,void *buf,int bufsize)
{
#ifndef GMX_MPI
    gmx_call("gmx_send"); 
#else
    int         tag=0;
    MPI_Status  status;
    MPI_Request req;
    
    if (MPI_Isend(buf,bufsize,MPI_BYTE,RANK(cr,dest),tag,
		  cr->mpi_comm_mygroup,&req) != MPI_SUCCESS)
      gmx_comm("MPI_Isend Failed");
    if (MPI_Wait(&req,&status) != MPI_SUCCESS)
      gmx_comm("MPI_Wait failed");
#endif
}

void gmx_recv(const t_commrec *cr,int src,void *buf,int bufsize)
{
#ifndef GMX_MPI
    gmx_call("gmx_recv");
#else
    int         tag=0;
    MPI_Status  status;
    MPI_Request req;
    
    if (MPI_Irecv(buf, bufsize, MPI_BYTE, src, tag, 
		  cr->mpi_comm_mygroup,&req) != MPI_SUCCESS)
      gmx_comm("MPI_Irecv Failed");
    if (MPI_Wait(&req,&status) != MPI_SUCCESS)
      gmx_comm("MPI_Wait failed");
#endif
}

void gmx_send_str(t_commrec *cr,int dest,const char *ptr)
{
    int len;
    
    if (NULL == ptr)
        len = 0;
    else
        len = strlen(ptr)+1;
    if (NULL != debug)
        fprintf(debug,"Sending string '%s' to %d\n",ptr,dest);
    gmx_send(cr,dest,&len,sizeof(len));
    if (NULL != ptr)
        gmx_send(cr,dest,(void *)ptr,len); 
}

char *gmx_recv_str(t_commrec *cr,int src)
{
    int  len;
    char *buf;
    
    gmx_recv(cr,src,&len,sizeof(len));
    if (len == 0)
    {
        return NULL;
    }
    else
    {
        snew(buf,len);
        gmx_recv(cr,src,buf,len);
        
        return buf;
    }
}

void gmx_send_double(t_commrec *cr,int dest,double d)
{
    if (NULL != debug)
        fprintf(debug,"Sending double '%g' to %d\n",d,dest);
    gmx_send(cr,dest,&d,sizeof(d)); 
}

double gmx_recv_double(t_commrec *cr,int src)
{
    double d;
    
    gmx_recv(cr,src,&d,sizeof(d));
    
    return d;
}

void gmx_send_int(t_commrec *cr,int dest,int d)
{
    if (NULL != debug)
        fprintf(debug,"Sending int '%d' to %d\n",d,dest);
    gmx_send(cr,dest,&d,sizeof(d)); 
}

int gmx_recv_int(t_commrec *cr,int src)
{
    int d;
    
    gmx_recv(cr,src,&d,sizeof(d));
    
    return d;
}

