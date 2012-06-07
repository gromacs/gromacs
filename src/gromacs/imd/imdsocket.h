/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2014, by the GROMACS development team, led by
 * Mark Abraham, David van der Spoel, Berk Hess, and Erik Lindahl,
 * and including many others, as listed in the AUTHORS file in the
 * top-level source directory and at http://www.gromacs.org.
 *
 * GROMACS is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public License
 * as published by the Free Software Foundation; either version 2.1
 * of the License, or (at your option) any later version.
 *
 * GROMACS is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with GROMACS; if not, see
 * http://www.gnu.org/licenses, or write to the Free Software Foundation,
 * Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA.
 *
 * If you want to redistribute modifications to GROMACS, please
 * consider that scientific software is very special. Version
 * control is crucial - bugs must be traceable. We will be happy to
 * consider code for inclusion in the official distribution, but
 * derived work must not be called official GROMACS. Details are found
 * in the README & COPYING files - if they are missing, get the
 * official version at http://www.gromacs.org.
 *
 * To help us fund GROMACS development, we humbly ask that you cite
 * the research papers on the package. Check out http://www.gromacs.org.
 */

/* \libinternal \file
 *
 * \brief
 * Implements parts of vmdsock.h interface from the original IMD API from scratch.
 *
 * \authors Martin Hoefling, Carsten Kutzner <ckutzne@gwdg.de>
 *
 * For more information, see:
 * IMD general information: https://www-s.ks.uiuc.edu/Research/vmd/imd/
 * IMD reference implementation/API https://www-s.ks.uiuc.edu/Research/vmd/imd/code/imdapi.tar.gz
 *
 * \inlibraryapi
 * \ingroup group_mdrun
 */

#ifndef GMX_IMD_IMDSOCKET_H
#define GMX_IMD_IMDSOCKET_H

/* Check if we should use winsock or standard unix sockets. */
#ifdef GMX_NATIVE_WINDOWS
#include <Winsock.h>
#define NOFLAGS 0
#else
#include <sys/socket.h>
#include <netinet/in.h>
#endif


#ifdef __cplusplus
extern "C" {
#endif

/* The IMD socket structure: */
typedef struct
{
    struct sockaddr_in address;
    int                addresslen;
    int                sockfd;
} IMDSocket;


/* public functions, see imdsocket.c */

#ifdef GMX_NATIVE_WINDOWS
gmx_bool imdsock_winsockinit();
#endif

extern IMDSocket *imdsock_create();

extern int imdsock_bind(IMDSocket *sock, int port);

extern int imd_sock_listen(IMDSocket *sock);

extern IMDSocket *imdsock_accept(IMDSocket *sock);

extern int imdsock_write(IMDSocket *sock, const char *, int len);

extern int imdsock_read(IMDSocket *sock, char *, int len);

extern void imdsock_shutdown(IMDSocket *sock);

extern int imdsock_destroy(IMDSocket *sock);

extern int imdsock_tryread(IMDSocket *sock, int timeoutsec, int timeoutusec);

#ifdef __cplusplus
}
#endif

#endif
