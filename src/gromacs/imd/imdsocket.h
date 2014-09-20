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

/*! \libinternal \file
 *
 * \brief
 * Implements the parts of the vmdsock.h interface needed for IMD communication.
 *
 * \author Martin Hoefling, Carsten Kutzner <ckutzne@gwdg.de>
 *
 * For more information, see https://www-s.ks.uiuc.edu/Research/vmd/imd/ for general
 * IMD information and https://www-s.ks.uiuc.edu/Research/vmd/imd/code/imdapi.tar.gz
 * for the IMD reference implementation and API.
 *
 * \inlibraryapi
 * \ingroup module_imd
 */

#ifndef GMX_IMD_IMDSOCKET_H
#define GMX_IMD_IMDSOCKET_H

#include "config.h"

/* Check if we can/should use winsock or standard UNIX sockets. */
#ifdef GMX_NATIVE_WINDOWS
  #ifdef GMX_HAVE_WINSOCK
  #include <Winsock.h>
  #define GMX_IMD
  #endif
#else
#include <netinet/in.h>
#include <sys/socket.h>
#define GMX_IMD
#endif


#ifdef __cplusplus
extern "C" {
#endif



/*! \internal
 *
 * \brief
 * IMD (interactive molecular dynamics) socket structure
 *
 */
typedef struct
{
#ifdef GMX_IMD
    struct sockaddr_in address;      /**< address of socket                   */
#endif
    int                sockfd;       /**< socket file descriptor              */
} IMDSocket;



#if defined(GMX_NATIVE_WINDOWS) && defined(GMX_HAVE_WINSOCK)
/*! \internal
 *
 * \brief Function to initialize winsock
 *
 * \returns 0 if successful.
 */
extern int imdsock_winsockinit();
#endif


/*! \brief Create an IMD master socket.
 *
 * \returns  The IMD socket if successful. Otherwise prints an error message and returns NULL.
 */
extern IMDSocket *imdsock_create();


/*! \brief Bind the IMD socket to address and port.
 *
 * Prints out an error message if unsuccessful.
 * If port == 0, bind() assigns a free port automatically.
 *
 *
 * \param sock      The IMD socket.
 * \param port      The port to bind to.
 *
 * \returns 0 if successful.
 */
extern int imdsock_bind(IMDSocket *sock, int port);


/*! \brief Set socket to listening state.
 *
 * Prints out an error message if unsuccessful.
 *
 * \param sock      The IMD socket.
 *
 * \returns 0 if successful.
 *
 */
extern int imd_sock_listen(IMDSocket *sock);


/*! \brief Accept incoming connection and redirect to client socket.
 *
 * Prints out an error message if unsuccessful.
 *
 * \param sock      The IMD socket.
 *
 * \returns IMD socket if successful, NULL otherwise.
 */
extern IMDSocket *imdsock_accept(IMDSocket *sock);


/*! \brief Get the port number used for IMD connection.
 *
 * Prints out an error message if unsuccessful.
 *
 * \param sock      The IMD socket.
 * \param port      The assigned port number.
 *
 * \returns 0 if successful, an error code otherwise.
 */
extern int imdsock_getport(IMDSocket *sock, int *port);


/*! \brief  Write to socket.
 *
 *
 * \param sock      The IMD socket.
 * \param buffer    The data to write.
 * \param length    Number of bytes to write.
 *
 * \returns The number of bytes written, or -1.
 */
extern int imdsock_write(IMDSocket *sock, const char *buffer, int length);


/*! \brief  Read from socket.
 *
 * \param sock      The IMD socket.
 * \param buffer    Buffer to put the read data.
 * \param length    Number of bytes to read.
 *
 * \returns The number of bytes read, or -1 for errors.
 */
extern int imdsock_read(IMDSocket *sock, char *buffer, int length);


/*! \brief Shutdown the socket.
 *
 * \param sock      The IMD socket.
 *
 */
extern void imdsock_shutdown(IMDSocket *sock);


/*! \brief Close the socket and free the sock struct memory.
 *
 * Writes an error message if unsuccessful.
 *
 * \param sock      The IMD socket.
 *
 * \returns 1 on success, or 0 if unsuccessful.
 */
extern int imdsock_destroy(IMDSocket *sock);


/*! \brief Try to read from the socket.
 *
 * Time out after waiting the interval specified.
 * Print an error message if unsuccessful.
 *
 * \param sock         The IMD socket.
 * \param timeoutsec   Time out seconds
 * \param timeoutusec  Time out microseconds
 *
 */
extern int imdsock_tryread(IMDSocket *sock, int timeoutsec, int timeoutusec);


#ifdef __cplusplus
}
#endif


#endif
