/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2014- The GROMACS Authors
 * and the project initiators Erik Lindahl, Berk Hess and David van der Spoel.
 * Consult the AUTHORS/COPYING files and https://www.gromacs.org for details.
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
 * https://www.gnu.org/licenses, or write to the Free Software Foundation,
 * Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA.
 *
 * If you want to redistribute modifications to GROMACS, please
 * consider that scientific software is very special. Version
 * control is crucial - bugs must be traceable. We will be happy to
 * consider code for inclusion in the official distribution, but
 * derived work must not be called official GROMACS. Details are found
 * in the README & COPYING files - if they are missing, get the
 * official version at https://www.gromacs.org.
 *
 * To help us fund GROMACS development, we humbly ask that you cite
 * the research papers on the package. Check out https://www.gromacs.org.
 */

/*! \internal \file
 *
 * \brief
 * Implements functions of imdsocket.h.
 *
 * This file re-implements vmdsock.c functions from original IMD API from scratch,
 * see imdsocket.h for references to the IMD API.
 *
 * \author Martin Hoefling, Carsten Kutzner <ckutzne@gwdg.de>
 *
 * \ingroup module_imd
 */

#include "gmxpre.h"

#include "imdsocket.h"

#include "config.h"

#include <cerrno>
#include <cstdio>
#include <cstring>
#include <ctime>

#include "gromacs/imd/imd.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/smalloc.h"

#if GMX_IMD

#    if GMX_NATIVE_WINDOWS

#        include <Windows.h>
#        include <Winsock.h>

//! Constant for passing no flags
constexpr int c_noFlags = 0;
/*! \brief Define socklen type on Windows. */
typedef int socklen_t;

#    else

#        include <unistd.h>

#        include <netinet/in.h>
#        include <sys/select.h>
#        include <sys/socket.h>
#        include <sys/time.h>

#    endif

#endif

namespace gmx
{

/*! \internal
 *
 * \brief
 * IMD (interactive molecular dynamics) socket structure
 *
 */
struct IMDSocket
{
#if GMX_IMD
    struct sockaddr_in address; /**< address of socket                   */
    int                sockfd;  /**< socket file descriptor              */
#endif
};

/*! \brief Define a function to initialize winsock. */
int imdsock_winsockinit()
{
#if GMX_IMD && GMX_NATIVE_WINDOWS
    fprintf(stderr, "%s Initializing winsock.\n", IMDstr);
    int ret = -1;

    WSADATA wsd;

    /* We use winsock 1.1 compatibility for now. Though I guess no one will try on Windows 95. */
    ret = WSAStartup(MAKEWORD(1, 1), &wsd);
    return ret;
#else
    return 0;
#endif
}


/*! \brief Simple error handling. */
#if GMX_NATIVE_WINDOWS
#    define ERR_ARGS __FILE__, __LINE__, NULL
#else
#    define ERR_ARGS __FILE__, __LINE__, strerror(errno)
#endif


/*! \brief Currently only 1 client connection is supported. */
constexpr int c_maxConnections = 1;


/*! \brief Print a nice error message on UNIX systems, using errno.h. */
static void print_IMD_error(const char* file, int line, char* msg)
{
    fprintf(stderr, "%s Error in file %s on line %d.\n", IMDstr, file, line);

    if (nullptr != msg)
    {
        fprintf(stderr, "%s\n", msg);
    }
}


IMDSocket* imdsock_create()
{
    IMDSocket* sock = nullptr;


#if GMX_IMD
    snew(sock, 1);
    /* Try to create socket: */
    if ((sock->sockfd = socket(PF_INET, SOCK_STREAM, 0)) == -1)
    {
        print_IMD_error(ERR_ARGS);
        sfree(sock);

        return nullptr;
    }
    else
#endif
    {
        return sock;
    }
}

void imd_sleep(unsigned int seconds)
{
#if GMX_IMD
#    if GMX_NATIVE_WINDOWS
    Sleep(seconds);
#    else
    sleep(seconds);
#    endif
#else
    GMX_UNUSED_VALUE(seconds);
#endif
}


int imdsock_bind(IMDSocket* sock, int port)
{
    int ret;


#if GMX_IMD
    memset(&(sock->address), 0, sizeof(sock->address));
    sock->address.sin_family = PF_INET;
    sock->address.sin_port   = htons(port);

    /* Try to bind to address and port ...*/
    ret = bind(sock->sockfd, reinterpret_cast<struct sockaddr*>(&sock->address), sizeof(sock->address));
#else
    GMX_UNUSED_VALUE(port);
    GMX_UNUSED_VALUE(sock);
    ret = -1;
#endif

    if (ret)
    {
        print_IMD_error(ERR_ARGS);
    }

    return ret;
}


int imd_sock_listen(IMDSocket* sock)
{
    int ret;


#if GMX_IMD
    /* Try to set to listening state */
    ret = listen(sock->sockfd, c_maxConnections);
#else
    GMX_UNUSED_VALUE(c_maxConnections);
    GMX_UNUSED_VALUE(sock);
    ret = -1;
#endif

    if (ret)
    {
        print_IMD_error(ERR_ARGS);
    }

    return ret;
}


IMDSocket* imdsock_accept(IMDSocket* sock)
{

#if GMX_IMD
    socklen_t length = sizeof(sock->address);
    int ret = accept(sock->sockfd, reinterpret_cast<struct sockaddr*>(&sock->address), &length);

    /* successful, redirect to distinct clientsocket */
    if (ret >= 0)
    {
        IMDSocket* newsock;

        snew(newsock, 1);
        newsock->address = sock->address;
        newsock->sockfd  = ret;

        return newsock;
    }
    else
#else
    GMX_UNUSED_VALUE(sock);
#endif
    {
        print_IMD_error(ERR_ARGS);

        return nullptr;
    }
}


int imdsock_getport(IMDSocket* sock, int* port)
{
    int ret;
#if GMX_IMD
    socklen_t len;


    len = sizeof(struct sockaddr_in);
    ret = getsockname(sock->sockfd, reinterpret_cast<struct sockaddr*>(&(sock->address)), &len);
    if (ret)
    {
        fprintf(stderr, "%s getsockname failed with error %d.\n", IMDstr, ret);
        print_IMD_error(ERR_ARGS);
    }
    else
    {
        *port = ntohs(sock->address.sin_port);
    }
#else
    GMX_UNUSED_VALUE(port);
    GMX_UNUSED_VALUE(sock);
    gmx_incons("imdsock_getport called without IMD support.");
    ret = -1;
#endif

    return ret;
}


int imd_htonl(int src)
{
#if GMX_IMD
    return htonl(src);
#else
    return src;
#endif
}

int imd_ntohl(int src)
{
#if GMX_IMD
    return ntohl(src);
#else
    return src;
#endif
}

int imdsock_write(IMDSocket* sock, const char* buffer, int length)
{
#if GMX_IMD
    /* No read and write on windows, we have to use send and recv instead... */
#    if GMX_NATIVE_WINDOWS
    return send(sock->sockfd, (const char*)buffer, length, c_noFlags);
#    else
    return write(sock->sockfd, buffer, length);
#    endif
#else
    GMX_UNUSED_VALUE(buffer);
    GMX_UNUSED_VALUE(length);
    GMX_UNUSED_VALUE(sock);
    return 0;
#endif
}


int imdsock_read(IMDSocket* sock, char* buffer, int length)
{
#if GMX_IMD
    /* See above... */
#    if GMX_NATIVE_WINDOWS
    return recv(sock->sockfd, (char*)buffer, length, c_noFlags);
#    else
    return read(sock->sockfd, buffer, length);
#    endif
#else
    GMX_UNUSED_VALUE(buffer);
    GMX_UNUSED_VALUE(length);
    GMX_UNUSED_VALUE(sock);
    return 0;
#endif
}


void imdsock_shutdown(IMDSocket* sock)
{
    int ret = -1;


    /* is the socket already NULL? */
    if (sock == nullptr)
    {
        return;
    }

#if GMX_IMD
    /* If not, try to properly shut down. */
    ret = shutdown(sock->sockfd, 1);
#endif

    if (ret == -1)
    {
        fprintf(stderr, "%s Failed to shutdown socket. Did the client already disconnect?\n", IMDstr);
        print_IMD_error(ERR_ARGS);
    }
}


int imdsock_destroy(IMDSocket* sock)
{
    int ret = -1;


    if (sock == nullptr)
    {
        return 1;
    }

#if GMX_IMD
#    if GMX_NATIVE_WINDOWS
    ret = closesocket(sock->sockfd);
#    else
    ret = close(sock->sockfd);
#    endif
#endif

    if (ret == -1)
    {
        sfree(sock);
        print_IMD_error(ERR_ARGS);

        return 0;
    }

    return 1;
}


int imdsock_tryread(IMDSocket* sock, int timeoutsec, int timeoutusec)
{
    int ret = -1;


#if GMX_IMD
    fd_set readfds;
    /* Create new time structure with sec and usec. */
    struct timeval* tval;


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
        ret = select(sock->sockfd + 1, &readfds, nullptr, nullptr, tval);
        /* redo on system interrupt */
    } while (ret < 0 && errno == EINTR);

    sfree(tval);
#else
    GMX_UNUSED_VALUE(sock);
    GMX_UNUSED_VALUE(timeoutsec);
    GMX_UNUSED_VALUE(timeoutusec);
#endif

    if (ret < 0)
    {
        print_IMD_error(ERR_ARGS);
    }

    return ret;
}

} // namespace gmx
