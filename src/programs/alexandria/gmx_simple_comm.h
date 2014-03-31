/*
 * This source file is part of the Alexandria project.
 *
 * Copyright (C) 2014 David van der Spoel and Paul J. van Maaren
 *
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.
 */
/*! \internal \brief
 * Implements part of the alexandria program.
 * \author David van der Spoel <david.vanderspoel@icm.uu.se>
 */
#ifndef GMX_SIMPLE_COMM_H
#define GMX_SIMPLE_COMM_H

#include "gromacs/legacyheaders/types/commrec.h"
#include "gromacs/legacyheaders/typedefs.h"
#include "gromacs/legacyheaders/network.h"

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
