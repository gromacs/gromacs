/* This file is part of Gromacs, a molecular simulation package.
 * Copyright (c) 2012- The Gromacs development team
 * (see AUTHORS file in the distribution)
 *
 * Gromacs is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 *
 * See COPYING file in the distribution and http://www.gromacs.org
 * for more information.
 */

/*This file implements parts of vmdsock.h interface from original IMD API from scratch,
 *For more information, see:
 *IMD general information: https://www-s.ks.uiuc.edu/Research/vmd/imd/
 *IMD reference implementation/API https://www-s.ks.uiuc.edu/Research/vmd/imd/code/imdapi.tar.gz
 */

#ifndef IMDSOCK_H__
#define IMDSOCK_H__

/*Currently only 1 client connection is supported.*/
#define MAXIMDCONNECTIONS 1

/*Check if we should use winsock or standard unix sockets.*/
#if ((defined WIN32 || defined _WIN32 || defined WIN64 || defined _WIN64) && !defined __CYGWIN__ && !defined __CYGWIN32__)
#include <Winsock.h>
#define NOFLAGS 0
#else
#include <sys/socket.h>
#include <netinet/in.h>
#endif

/*gromacs includes*/
#include "typedefs.h"

/*the IMD socket structure*/
typedef struct
{
    struct sockaddr_in address;
    int addresslen;
    int sockfd;
} IMDSocket;

/*public functions, see imdsocket.c*/

#if ((defined WIN32 || defined _WIN32 || defined WIN64 || defined _WIN64) && !defined __CYGWIN__ && !defined __CYGWIN32__)
gmx_bool imdsock_winsockinit();
#endif

IMDSocket *imdsock_create();
int imdsock_bind(IMDSocket *sock, int port);
int imd_sock_listen(IMDSocket *sock);
IMDSocket *imdsock_accept(IMDSocket *sock);
int imdsock_write(IMDSocket *sock, const char *, int len);
int imdsock_read(IMDSocket *sock, char *, int len);
void imdsock_shutdown(IMDSocket *sock);
int imdsock_destroy(IMDSocket *sock);
int imdsock_tryread(IMDSocket *sock, int timeoutsec, int timeoutusec);

#endif

