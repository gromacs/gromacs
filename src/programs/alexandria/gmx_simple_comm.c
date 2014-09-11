/*
 * This source file is part of the Alexandria project.
 *
 * Copyright (C) 2014 David van der Spoel and Paul J. van Maaren
 *
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.
 */
/*! \internal \brief
 * Implements part of the alexandria program.
 * \author David van der Spoel <david.vanderspoel@icm.uu.se>
 */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <string.h>
#ifdef GMX_LIB_MPI
#include <mpi.h>
#endif

#include "gromacs/utility/smalloc.h"
#include "gmx_simple_comm.h"

void gmx_send(const t_commrec *cr, int dest, void *buf, int bufsize)
{
#ifndef GMX_MPI
    gmx_call("gmx_send");
#else
    int         tag = 0;
    MPI_Status  status;
    MPI_Request req;

    if (MPI_Isend(buf, bufsize, MPI_BYTE, RANK(cr, dest), tag,
                  cr->mpi_comm_mygroup, &req) != MPI_SUCCESS)
    {
        gmx_comm("MPI_Isend Failed");
    }
    if (MPI_Wait(&req, &status) != MPI_SUCCESS)
    {
        gmx_comm("MPI_Wait failed");
    }
#endif
}

void gmx_recv(const t_commrec *cr, int src, void *buf, int bufsize)
{
#ifndef GMX_MPI
    gmx_call("gmx_recv");
#else
    int         tag = 0;
    MPI_Status  status;
    MPI_Request req;

    if (MPI_Irecv(buf, bufsize, MPI_BYTE, src, tag,
                  cr->mpi_comm_mygroup, &req) != MPI_SUCCESS)
    {
        gmx_comm("MPI_Irecv Failed");
    }
    if (MPI_Wait(&req, &status) != MPI_SUCCESS)
    {
        gmx_comm("MPI_Wait failed");
    }
#endif
}

void gmx_send_str(t_commrec *cr, int dest, const char *ptr)
{
    int len;

    if (NULL == ptr)
    {
        len = 0;
    }
    else
    {
        len = strlen(ptr)+1;
    }
    if (NULL != debug)
    {
        fprintf(debug, "Sending string '%s' to %d\n", ptr, dest);
    }
    gmx_send(cr, dest, &len, sizeof(len));
    if (NULL != ptr)
    {
        gmx_send(cr, dest, (void *)ptr, len);
    }
}

char *gmx_recv_str(t_commrec *cr, int src)
{
    int   len;
    char *buf;

    gmx_recv(cr, src, &len, sizeof(len));
    if (len == 0)
    {
        return NULL;
    }
    else
    {
        snew(buf, len);
        gmx_recv(cr, src, buf, len);

        return buf;
    }
}

void gmx_send_double(t_commrec *cr, int dest, double d)
{
    if (NULL != debug)
    {
        fprintf(debug, "Sending double '%g' to %d\n", d, dest);
    }
    gmx_send(cr, dest, &d, sizeof(d));
}

double gmx_recv_double(t_commrec *cr, int src)
{
    double d;

    gmx_recv(cr, src, &d, sizeof(d));

    return d;
}

void gmx_send_int(t_commrec *cr, int dest, int d)
{
    if (NULL != debug)
    {
        fprintf(debug, "Sending int '%d' to %d\n", d, dest);
    }
    gmx_send(cr, dest, &d, sizeof(d));
}

int gmx_recv_int(t_commrec *cr, int src)
{
    int d;

    gmx_recv(cr, src, &d, sizeof(d));

    return d;
}
