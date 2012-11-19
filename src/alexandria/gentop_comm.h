#ifndef GENTOP_COMM_H
#define GENTOP_COMM_H

extern void gmx_send(const t_commrec *cr,int dest,const void *buf,int bufsize);

extern void gmx_recv(const t_commrec *cr,int src,void *buf,int bufsize);

extern void gmx_send_str(t_commrec *cr,int dest,const char *ptr);

extern char *gmx_recv_str(t_commrec *cr,int src);

extern void gmx_send_double(t_commrec *cr,int dest,double d);

extern double gmx_recv_double(t_commrec *cr,int src);

extern void gmx_send_int(t_commrec *cr,int dest,int d);

extern int gmx_recv_int(t_commrec *cr,int src);

#endif
