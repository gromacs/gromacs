/*This file re-implements vmdsock.h from original IMD API*/
#ifdef GMX_IMD

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

#endif
