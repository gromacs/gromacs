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

/*! \libinternal \file
 *
 * \brief
 * Implements the parts of the vmdsock.h interface needed for IMD communication.
 *
 * \author Martin Hoefling, Carsten Kutzner <ckutzne@gwdg.de>
 *
 * For more information, see https://www.ks.uiuc.edu/Research/vmd/imd/ for general
 * IMD information and https://www.ks.uiuc.edu/Research/vmd/imd/code/imdapi.tar.gz
 * for the IMD reference implementation and API.
 *
 * \inlibraryapi
 * \ingroup module_imd
 */

#ifndef GMX_IMD_IMDSOCKET_H
#define GMX_IMD_IMDSOCKET_H

namespace gmx
{

struct IMDSocket;


/*! \internal
 *
 * \brief Function to initialize winsock
 *
 * \returns 0 if successful.
 */
int imdsock_winsockinit();

/*! \brief Create an IMD main socket.
 *
 * \returns  The IMD socket if successful. Otherwise prints an error message and returns NULL.
 */
IMDSocket* imdsock_create();

//! Portability wrapper around sleep function
void imd_sleep(unsigned int seconds);


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
int imdsock_bind(IMDSocket* sock, int port);


/*! \brief Set socket to listening state.
 *
 * Prints out an error message if unsuccessful.
 *
 * \param sock      The IMD socket.
 *
 * \returns 0 if successful.
 *
 */
int imd_sock_listen(IMDSocket* sock);


/*! \brief Accept incoming connection and redirect to client socket.
 *
 * Prints out an error message if unsuccessful.
 *
 * \param sock      The IMD socket.
 *
 * \returns IMD socket if successful, NULL otherwise.
 */
IMDSocket* imdsock_accept(IMDSocket* sock);


/*! \brief Get the port number used for IMD connection.
 *
 * Prints out an error message if unsuccessful.
 *
 * \param sock      The IMD socket.
 * \param port      The assigned port number.
 *
 * \returns 0 if successful, an error code otherwise.
 */
int imdsock_getport(IMDSocket* sock, int* port);

//! Portability wrapper around system htonl function.
int imd_htonl(int src);

//! Portability wrapper around system ntohl function.
int imd_ntohl(int src);


/*! \brief  Write to socket.
 *
 *
 * \param sock      The IMD socket.
 * \param buffer    The data to write.
 * \param length    Number of bytes to write.
 *
 * \returns The number of bytes written, or -1.
 */
int imdsock_write(IMDSocket* sock, const char* buffer, int length);


/*! \brief  Read from socket.
 *
 * \param sock      The IMD socket.
 * \param buffer    Buffer to put the read data.
 * \param length    Number of bytes to read.
 *
 * \returns The number of bytes read, or -1 for errors.
 */
int imdsock_read(IMDSocket* sock, char* buffer, int length);


/*! \brief Shutdown the socket.
 *
 * \param sock      The IMD socket.
 *
 */
void imdsock_shutdown(IMDSocket* sock);


/*! \brief Close the socket and free the sock struct memory.
 *
 * Writes an error message if unsuccessful.
 *
 * \param sock      The IMD socket.
 *
 * \returns 1 on success, or 0 if unsuccessful.
 */
int imdsock_destroy(IMDSocket* sock);


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
int imdsock_tryread(IMDSocket* sock, int timeoutsec, int timeoutusec);

} // namespace gmx

#endif
