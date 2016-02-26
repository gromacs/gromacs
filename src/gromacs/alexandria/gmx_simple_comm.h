/*! \internal \brief
 * Implements part of the alexandria program.
 * \author David van der Spoel <david.vanderspoel@icm.uu.se>
 */
#ifndef GMX_SIMPLE_COMM_H
#define GMX_SIMPLE_COMM_H

#include "gromacs/gmxlib/network.h"
#include "gromacs/mdtypes/commrec.h"

void gmx_send(const t_commrec *cr, int dest, void *buf, int bufsize);

void gmx_recv(const t_commrec *cr, int src, void *buf, int bufsize);

void gmx_send_str(t_commrec *cr, int dest, const char *ptr);

char *gmx_recv_str(t_commrec *cr, int src);

void gmx_send_double(t_commrec *cr, int dest, double d);

double gmx_recv_double(t_commrec *cr, int src);

void gmx_send_int(t_commrec *cr, int dest, int d);

int gmx_recv_int(t_commrec *cr, int src);

#endif
