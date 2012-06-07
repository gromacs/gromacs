#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#ifdef GMX_IMD

/*This file re-implements vmdsock.c functions from original IMD API*/

#include <string.h>

/*gromacs includes*/
#include "smalloc.h"
#include "imdsocket.h"
#include "gmx_fatal.h"
#include "imd.h"

#if ((defined WIN32 || defined _WIN32 || defined WIN64 || defined _WIN64) && !defined __CYGWIN__ && !defined __CYGWIN32__)
/*define socklen type on windows*/
typedef int socklen_t;

/*and a function to initialize winsock.*/
extern gmx_bool imdsock_winsockinit()
{
    WSADATA wsd;
    /*We use winsock 1.1 compatibility for now. Though I guess no one will try on Windows 95.*/
    return (0==WSAStartup(MAKEWORD(1,1),&wsd));
}
#else
/*On unix, we can use nice errors from errno.h*/
#define UNIXERR
#include <errno.h>
#include <string.h>
#include <unistd.h>
#endif

/*creates an IMD master socket.*/
extern IMDSocket* imdsock_create()
{
    IMDSocket *sock;
    snew(sock, 1);
    /*try to create socket*/
    if ((sock->sockfd = socket(PF_INET, SOCK_STREAM, 0)) == -1)
    {
#ifdef UNIXERR
        /*print actual error on unix*/
        fprintf(stderr, "%s %s\n", IMDstr, strerror(errno));
#endif
        sfree(sock);
        return NULL;
    } 
    else
    {
        return sock;
    }
}

/*binds the socket...*/
extern int imdsock_bind(IMDSocket *sock, int port)
{
    int ret;
    memset(&(sock->address), 0, sizeof(sock->address));
    sock->address.sin_family = PF_INET;
    sock->address.sin_port = htons(port);
    /*try to bind to address and port...*/
    ret = bind(sock->sockfd, (struct sockaddr *) &sock->address, sizeof(sock->address));
#ifdef UNIXERR
    if (ret != 0)
    {
        /*print actual error on unix.*/
        fprintf(stderr, "%s %s\n", IMDstr, strerror(errno));
        ret = errno;
    }
#endif
    return ret;
}

/*sets socket to listen state*/
extern int imd_sock_listen(IMDSocket *sock)
{
    int ret;
    /*try to set to listen state*/
    ret = listen(sock->sockfd, MAXIMDCONNECTIONS);
#ifdef UNIXERR
    if (ret != 0)
    {
        /*print actual error on unix.*/
        fprintf(stderr, "%s %s\n", IMDstr, strerror(errno));
        ret = errno;
    }
#endif
    return ret;
}

/*accept incoming connection and redirect to client socket.*/
extern IMDSocket* imdsock_accept(IMDSocket *sock)
{
    int ret;
    socklen_t length;

    length = sizeof(sock->address);
    ret = accept(sock->sockfd, (struct sockaddr *) &sock->address, &length);

    /*successful, redirect to distinct clientsocket*/
    if (ret >= 0)
    {
        IMDSocket *newsock;
        snew(newsock, 1);
        newsock->address = sock->address;
        newsock->addresslen = sock->addresslen;
        newsock->sockfd = ret;
        return newsock;
    }
    /*error otherwise*/
    else
    {
#ifdef UNIXERR
        if (ret != 0)
        {
            /*print unix error.*/
            fprintf(stderr, "%s %s\n", IMDstr, strerror(errno));
        }
#endif
        return NULL;
    }
}

/*write to socket*/
extern int imdsock_write(IMDSocket *sock, const char *buffer, int length)
{
    /*No read and write on windows, we have to use send and recv instead...*/
#if ((defined WIN32 || defined _WIN32 || defined WIN64 || defined _WIN64) && !defined __CYGWIN__ && !defined __CYGWIN32__)
    return send(sock->sockfd,(const char *) buffer, length, NOFLAGS);
#else
    return write(sock->sockfd, buffer, length);
#endif
}

/*read from socket*/
extern int imdsock_read(IMDSocket *sock, char *buffer, int length)
{
    /*see above...*/
#if ((defined WIN32 || defined _WIN32 || defined WIN64 || defined _WIN64) && !defined __CYGWIN__ && !defined __CYGWIN32__)
    return recv(sock->sockfd,(char *) buffer, length, NOFLAGS);
#else
    return read(sock->sockfd, buffer, length);
#endif
}

/*shutdown the socket*/
extern void imdsock_shutdown(IMDSocket *sock)
{
    /*is the socket already NULL?*/
    if (sock == NULL)
    {
        return;
    }
    /*If not try to properly shut down.*/
    if (shutdown(sock->sockfd, 1) == -1)
    {
        fprintf(stderr, "%s Failed to shutdown socket. Maybe the client has already disconnected?\n", IMDstr);
#ifdef UNIXERR
        /*print unix error*/
        fprintf(stderr, "%s %s\n", IMDstr, strerror(errno));
#endif
    }
}

/*closes the socket and frees the sock struct memory.*/
extern int imdsock_destroy(IMDSocket *sock)
{
    if (sock == NULL)
    {
        return TRUE;
    }
#if ((defined WIN32 || defined _WIN32 || defined WIN64 || defined _WIN64) && !defined __CYGWIN__ && !defined __CYGWIN32__)
    /*on windows, this function is called closesocket*/
    if (closesocket(sock->sockfd) == -1)
#else
    if (close(sock->sockfd) == -1)
#endif
    {
        sfree(sock);
#ifdef UNIXERR
        /*print unix error*/
        fprintf(stderr, "%s %s\n", IMDstr, strerror(errno));
#endif
        return FALSE;
    }
    return TRUE;
}

/*try to read from the socket*/
extern int imdsock_tryread(IMDSocket *sock, int timeoutsec, int timeoutusec)
{
    int ret;
    fd_set readfds;

    /*create new time structure with sec and usec*/
    struct timeval *tval;
    snew(tval, 1);

    /*clear the set*/
    FD_ZERO(&readfds);
    /*add the socket to the read set.*/
    FD_SET(sock->sockfd, &readfds);

    /*set the timeout*/
    tval->tv_sec = timeoutsec;
    tval->tv_usec = timeoutusec;
    do
    {
        /*check the set for read readiness.*/
        ret = select(sock->sockfd + 1, &readfds, NULL, NULL, tval);
        /*redo on system interrupt*/
    } while (ret < 0 && errno == EINTR);

    sfree(tval);
#ifdef UNIXERR
    if (ret < 0)
        /*otherwise print unix error*/
        fprintf(stderr, "%s %s\n", IMDstr, strerror(errno));
#endif
    return ret;
}

#endif
