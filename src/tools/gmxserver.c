/*
 *       @(#) copyrgt.c 1.12 9/30/97
 *
 *       This source code is part of
 *
 *        G   R   O   M   A   C   S
 *
 * GROningen MAchine for Chemical Simulations
 *
 *            VERSION 2.0b
 * 
 * Copyright (c) 1990-1997,
 * BIOSON Research Institute, Dept. of Biophysical Chemistry,
 * University of Groningen, The Netherlands
 *
 * Please refer to:
 * GROMACS: A message-passing parallel molecular dynamics implementation
 * H.J.C. Berendsen, D. van der Spoel and R. van Drunen
 * Comp. Phys. Comm. 91, 43-56 (1995)
 *
 * Also check out our WWW page:
 * http://rugmd0.chem.rug.nl/~gmx
 * or e-mail to:
 * gromacs@chem.rug.nl
 *
 * And Hey:
 * Gyas ROwers Mature At Cryogenic Speed
 */
#include <stdio.h>
#include <unistd.h>
#include <signal.h>
#include <sys/types.h>
#include <string.h>
#include <sys/wait.h>
#include <ctype.h>
#include <signal.h>

#include <stdlib.h>
#include <string.h>
#include <sys/types.h>
#include <sys/socket.h>
#include <errno.h>
#include <netinet/in.h>
#include <limits.h>
#include <netdb.h>
#include <arpa/inet.h>


/* Take a service name, and a service type, and return a port number.  If the
   service name is not found, it tries it as a decimal number.  The number
   returned is byte ordered for the network. */
int atoport(char *service, char *proto) {
  int port;
  long int lport;
  struct servent *serv;
  char *errpos;

  /* First try to read it from /etc/services */
  serv = getservbyname(service, proto);
  if (serv != NULL)
    port = serv->s_port;
  else { /* Not in services, maybe a number? */
    lport = strtol(service,&errpos,0);
    if ( (errpos[0] != 0) || (lport < 1) || (lport > 65535) )
      return -1; /* Invalid port address */
    port = htons(lport);
  }
  return port;
}

/* Converts ascii text to in_addr struct.  NULL is returned if the address
   can not be found. */
struct in_addr *atoaddr(char *address) {
  struct hostent *host;
  static struct in_addr saddr;

  /* First try it as aaa.bbb.ccc.ddd. */
  saddr.s_addr = inet_addr(address);
  if (saddr.s_addr != -1) {
    return &saddr;
  }
  host = gethostbyname(address);
  if (host != NULL) {
    return (struct in_addr *) *host->h_addr_list;
  }
  return NULL;
}

int get_connection(int socket_type, u_short port, int *listener) 
{
  struct sockaddr_in address;
  int listening_socket;
  int connected_socket = -1;
  int new_process;
  int reuse_addr = 1;

  /* Setup internet address information.  
     This is used with the bind() call */
  memset((char *) &address, 0, sizeof(address));
  address.sin_family = AF_INET;
  address.sin_port = port;
  address.sin_addr.s_addr = htonl(INADDR_ANY);

  listening_socket = socket(AF_INET, socket_type, 0);
  if (listening_socket < 0) {
    perror("socket");
    exit(-1);
  }

  if (listener != NULL)
    *listener = listening_socket;

  setsockopt(listening_socket, SOL_SOCKET, SO_REUSEADDR,(void *)&reuse_addr, sizeof(reuse_addr));

  if (bind(listening_socket, (struct sockaddr *) &address, 
    sizeof(address)) < 0) {
    perror("bind");
    close(listening_socket);
    exit(-1);
  }

  if (socket_type == SOCK_STREAM) {
    listen(listening_socket, 5); /* Queue up to five connections before
                                  having them automatically rejected. */

    while(connected_socket < 0) {
      connected_socket = accept(listening_socket, NULL, NULL);
      if (connected_socket < 0) {
        /* Either a real error occured, or blocking was interrupted for
           some reason.  Only abort execution if a real error occured. */
        if (errno != EINTR) {
          perror("accept");
          close(listening_socket);
          exit(-1);
        }
      }

      new_process = fork();
      if (new_process < 0) {
        perror("fork");
        close(connected_socket);
        connected_socket = -1;
      }
      else { /* We have a new process... */
        if (new_process == 0) {
          /* This is the new process. */
          close(listening_socket); /* Close our copy of this socket */
          *listener = -1; /* Closed in this process.  We are not responsible
                             for it. */
        }
        else {
          /* This is the main loop.  Close copy of connected socket, and
             continue loop. */
          close(connected_socket);
          connected_socket = -1;
        }
      }
    }
    return connected_socket;
  }
  else
    return listening_socket;
}

/* This is a generic function to make a connection to a given server/port.
   service is the port name/number,
   type is either SOCK_STREAM or SOCK_DGRAM, and
   netaddress is the host name to connect to.
   The function returns the socket, ready for action.*/
int make_connection(char *service, int type, char *netaddress) {
  /* First convert service from a string, to a number... */
  int port = -1;
  struct in_addr *addr;
  int sock, connected;
  struct sockaddr_in address;

  if (type == SOCK_STREAM) 
    port = atoport(service, "tcp");
  if (type == SOCK_DGRAM)
    port = atoport(service, "udp");
  if (port == -1) {
    return -1;
  }
  addr = atoaddr(netaddress);
  if (addr == NULL) {
    return -1;
  }
 
  memset((char *) &address, 0, sizeof(address));
  address.sin_family = AF_INET;
  address.sin_port = (port);
  address.sin_addr.s_addr = addr->s_addr;

  sock = socket(AF_INET, type, 0);

  if (type == SOCK_STREAM) {
    connected = connect(sock, (struct sockaddr *) &address, 
      sizeof(address));
    if (connected < 0) {
      return -1;
    }
    return sock;
  }
  /* Otherwise, must be for udp, so bind to address. */
  if (bind(sock, (struct sockaddr *) &address, sizeof(address)) < 0) {
    return -1;
  }
  return sock;
}

/* This is just like the read() system call, accept that it will make
   sure that all your data goes through the socket. */
int sock_read(int sockfd, char *buf, size_t count) {
  size_t bytes_read = 0;
  int this_read;

  while (bytes_read < count) {
    do
      this_read = read(sockfd, buf, count - bytes_read);
    while ( (this_read < 0) && (errno == EINTR) );
    if (this_read <= 0)
      return this_read;
    bytes_read += this_read;
    buf += this_read;
  }
  return count;
}

/* This is just like the write() system call, accept that it will
   make sure that all data is transmitted. */
int sock_write(int sockfd, const char *buf, size_t count) {
  size_t bytes_sent = 0;
  int this_write;

  while (bytes_sent < count) {
    do
      this_write = write(sockfd, buf, count - bytes_sent);
    while ( (this_write < 0) && (errno == EINTR) );
    if (this_write <= 0)
      return this_write;
    bytes_sent += this_write;
    buf += this_write;
  }
  return count;
}

/* This function reads from a socket, until it recieves a linefeed
   character.  It fills the buffer "str" up to the maximum size "count".

   This function will return -1 if the socket is closed during the read
   operation.

   Note that if a single line exceeds the length of count, the extra data
   will be read and discarded!  You have been warned. */
int sock_gets(int sockfd, char *str, size_t count) {
  int bytes_read;
  int total_count = 0;
  char *current_position;
  char last_read = 0;

  current_position = str;
  while (last_read != 10) {
    bytes_read = read(sockfd, &last_read, 1);
    if (bytes_read <= 0) {
      /* The other side may have closed unexpectedly */
      return -1; /* Is this effective on other platforms than linux? */
    }
    if ( (total_count < count) && (last_read != 10) && (last_read !=13) ) {
      current_position[0] = last_read;
      current_position++;
      total_count++;
    }
  }
  if (count > 0)
    current_position[0] = 0;
  return total_count;
}

/* This function writes a character string out to a socket.  It will 
   return -1 if the connection is closed while it is trying to write. */
int sock_puts(int sockfd, const char *str) {
  return sock_write(sockfd, str, strlen(str));
}















int listensock = -1; /* So that we can close sockets on ctrl-c */
int connectsock = -1;

/* This waits for all children, so that they don't become zombies. */
void sig_cld(int signal_type) {
#ifdef BSD
  int pid;
  int status;

  while ( (pid = wait3(&status, WNOHANG, NULL)) > 0);
#endif
}

typedef struct {
  char *name;
  int score;
} t_entr;

typedef int bool;

#define STRLEN 255
#define TRUE 1
#define FALSE 0

int main(int argc, char *argv[]) 
{
  int sock;
  int connected = 1;
  char buffer[1024];
  char *progname;
  int port = -1;
  struct sigaction act, oldact;

  FILE *fp;
  int nentr;
  t_entr *entr;
  int l,c;
  char name[255];
  bool bDone;
  char line[STRLEN];
  char *logfile;

  if (argc != 2) {
    fprintf(stderr,"Usage:  gmxserver logfile\n");
    exit(-1);
  }
  logfile = argv[1];
  sigemptyset(&act.sa_mask);
  act.sa_flags = 0;
#ifdef BSD
  act.sa_handler = sig_cld;
  sigaction(SIGCLD, &act, &oldact);
#else
  act.sa_handler = SIG_IGN;
  sigaction(SIGCLD, &act, &oldact);
#endif

  port = atoport("8123", "tcp");
  if (port == -1) {
    fprintf(stderr,"Unable to find service: %s\n","8123");
    exit(-1);
  }

  sock = get_connection(SOCK_STREAM, port, &listensock);

  connectsock = sock;
  sock_puts(sock,"Welcome to gromacs counter server.\n");
  while (connected) {
    /* Read input */
    if ( sock_gets(sock, buffer, 1024) < 0) {
      connected = 0;
    }
    else {

      /* delete leading slashes */
      progname = buffer;
      for(c=0;(c<(int)strlen(buffer));c++)
	if ( buffer[c] == '/' )
	  progname= buffer + c + 1;

      bDone = FALSE;
      fp = fopen(logfile,"r");
      fgets(line,STRLEN,fp);
      sscanf(line,"%d",&nentr);
      entr = (t_entr *)malloc(sizeof(t_entr) * nentr);
      for(l=0;(l<nentr);l++) {
	fgets(line,STRLEN,fp);
	sscanf(line,"%s %d",name,&(entr[l].score));
	entr[l].name = (char *)malloc(strlen(name) + 2);
	strcpy(entr[l].name,name);
	if ( strcmp(entr[l].name,progname)==0) {
	  entr[l].score++;
	  bDone=TRUE;
	}
      }
      fclose(fp);
      if ( bDone == FALSE ) {
	nentr++;
	entr = (t_entr *)realloc(entr,sizeof(t_entr) * nentr);
	entr[nentr-1].name = (char *)malloc(strlen(progname)+2);
	strcpy(entr[nentr-1].name,progname);
	entr[nentr-1].score=1;
      }

      fp = fopen(logfile,"w");
      fprintf(fp,"%d\n",nentr);
      for(l=0;(l<nentr);l++)
	fprintf(fp,"%s %d\n",entr[l].name,entr[l].score);
      fclose(fp);
      
      
      

      if (sock_puts(sock,"ok") < 0) {
        connected = 0;
      }
    }
  }
  close(sock);
  return 0;
}
