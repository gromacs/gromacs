#ifndef _GMX_SIMPLE_COMM_H
#define _GMX_SIMPLE_COMM_H

#include "network.h"

#ifdef __cplusplus
extern "C"
#endif
void gmx_send(const t_commrec *cr,int dest,void *buf,int bufsize);

#ifdef __cplusplus
extern "C"
#endif
void gmx_recv(const t_commrec *cr,int src,void *buf,int bufsize);

#ifdef __cplusplus
extern "C"
#endif
void gmx_send_str(t_commrec *cr,int dest,const char *ptr);

#ifdef __cplusplus
extern "C"
#endif
char *gmx_recv_str(t_commrec *cr,int src);

#ifdef __cplusplus
extern "C"
#endif
void gmx_send_double(t_commrec *cr,int dest,double d);

#ifdef __cplusplus
extern "C"
#endif
double gmx_recv_double(t_commrec *cr,int src);

#ifdef __cplusplus
extern "C"
#endif
void gmx_send_int(t_commrec *cr,int dest,int d);

#ifdef __cplusplus
extern "C"
#endif
int gmx_recv_int(t_commrec *cr,int src);

#endif
