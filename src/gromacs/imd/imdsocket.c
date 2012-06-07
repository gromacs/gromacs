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

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif


/* This file re-implements vmdsock.c functions from original IMD API from scratch,
 * see imdsocket.h for references to the IMD API */

#include <string.h>
#include "smalloc.h"
#include "imdsocket.h"
#include "imd.h"

#ifdef GMX_NATIVE_WINDOWS

/* define socklen type on windows */
typedef int socklen_t;


/* and a function to initialize winsock. */
extern int imdsock_winsockinit()
{
    WSADATA wsd;


    /* We use winsock 1.1 compatibility for now. Though I guess no one will try on Windows 95. */
    return (WSAStartup(MAKEWORD(1, 1), &wsd));
}
#else
/* On UNIX, we can use nice errors from errno.h */
#include <errno.h>
#include <unistd.h>
#endif


/* Error handling */
#ifdef GMX_NATIVE_WINDOWS
#define ERR_ARGS __FILE__, __LINE__, NULL
#else
#define ERR_ARGS __FILE__, __LINE__, strerror(errno)
#endif


/* Currently only 1 client connection is supported. */
#define MAXIMDCONNECTIONS 1


/* Print a nice error message on UNIX systems, using errno.h */
static void print_IMD_error(char *file, int line, char *msg)
{
    fprintf(stderr, "%s Error in file %s on line %d.\n", IMDstr, file, line);

    if (NULL != msg)
    {
        fprintf(stderr, "%s\n", msg);
    }
}


/* Create an IMD master socket. */
extern IMDSocket* imdsock_create()
{
    IMDSocket *sock;


    snew(sock, 1);
    /* Try to create socket: */
    if ((sock->sockfd = socket(PF_INET, SOCK_STREAM, 0)) == -1)
    {
        print_IMD_error(ERR_ARGS);
        sfree(sock);

        return NULL;
    }
    else
    {
        return sock;
    }
}


/* Bind the socket to address and port. */
extern int imdsock_bind(IMDSocket *sock, int port)
{
    int ret;


    memset(&(sock->address), 0, sizeof(sock->address));
    sock->address.sin_family = PF_INET;
    sock->address.sin_port   = htons(port);

    /* Try to bind to address and port ...*/
    ret = bind(sock->sockfd, (struct sockaddr *) &sock->address, sizeof(sock->address));

    if (ret)
    {
        print_IMD_error(ERR_ARGS);
    }

    return ret;
}


/* Set socket to listening state. */
extern int imd_sock_listen(IMDSocket *sock)
{
    int ret;


    /* Try to set to listening state */
    ret = listen(sock->sockfd, MAXIMDCONNECTIONS);

    if (ret)
    {
        print_IMD_error(ERR_ARGS);
    }

    return ret;
}


/* Accept incoming connection and redirect to client socket. */
extern IMDSocket* imdsock_accept(IMDSocket *sock)
{
    int       ret;
    socklen_t length;


    length = sizeof(sock->address);
    ret    = accept(sock->sockfd, (struct sockaddr *) &sock->address, &length);

    /* successful, redirect to distinct clientsocket */
    if (ret >= 0)
    {
        IMDSocket *newsock;

        snew(newsock, 1);
        newsock->address    = sock->address;
        newsock->addresslen = sock->addresslen;
        newsock->sockfd     = ret;

        return newsock;
    }
    else
    {
        print_IMD_error(ERR_ARGS);

        return NULL;
    }
}


/* Write to socket */
extern int imdsock_write(IMDSocket *sock, const char *buffer, int length)
{
    /* No read and write on windows, we have to use send and recv instead... */
#ifdef GMX_NATIVE_WINDOWS
#ifdef GMX_HAVE_WINSOCK
    return send(sock->sockfd, (const char *) buffer, length, NOFLAGS);
#endif
#else
    return write(sock->sockfd, buffer, length);
#endif
}

/* Read from socket */
extern int imdsock_read(IMDSocket *sock, char *buffer, int length)
{
    /* See above... */
#ifdef GMX_NATIVE_WINDOWS
#ifdef GMX_HAVE_WINSOCK
    return recv(sock->sockfd, (char *) buffer, length, NOFLAGS);
#endif
#else
    return read(sock->sockfd, buffer, length);
#endif
}

/* Shutdown the socket */
extern void imdsock_shutdown(IMDSocket *sock)
{
    /* is the socket already NULL? */
    if (sock == NULL)
    {
        return;
    }

    /* If not, try to properly shut down. */
    if (shutdown(sock->sockfd, 1) == -1)
    {
        fprintf(stderr, "%s Failed to shutdown socket. Did the client already disconnect?\n", IMDstr);
        print_IMD_error(ERR_ARGS);
    }
}


/* Close the socket and free the sock struct memory. */
extern int imdsock_destroy(IMDSocket *sock)
{
    int ret = -1;


    if (sock == NULL)
    {
        return 1;
    }

#ifdef GMX_NATIVE_WINDOWS
    /* On Windows, this function is called closesocket */
#ifdef GMX_HAVE_WINSOCK
    ret = closesocket(sock->sockfd);
#endif
#else
    ret = close(sock->sockfd);
#endif

    if (ret == -1)
    {
        sfree(sock);
        print_IMD_error(ERR_ARGS);

        return 0;
    }

    return 1;
}


/* Try to read from the socket */
extern int imdsock_tryread(IMDSocket *sock, int timeoutsec, int timeoutusec)
{
    int             ret;
    fd_set          readfds;
    /* Create new time structure with sec and usec. */
    struct timeval *tval;


    snew(tval, 1);

    /* clear the set */
    FD_ZERO(&readfds);
    /* add the socket to the read set */
    FD_SET(sock->sockfd, &readfds);

    /* set the timeout */
    tval->tv_sec  = timeoutsec;
    tval->tv_usec = timeoutusec;
    do
    {
        /* check the set for read readiness. */
        ret = select(sock->sockfd + 1, &readfds, NULL, NULL, tval);
        /* redo on system interrupt */
    }
    while (ret < 0 && errno == EINTR);

    sfree(tval);

    if (ret < 0)
    {
        print_IMD_error(ERR_ARGS);
    }

    return ret;
}
