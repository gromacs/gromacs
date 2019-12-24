/*
 * This source file is part of the Alexandria Chemistry Toolkit.
 *
 * Copyright (C) 2014-2020 
 *
 * Developers:
 *             Mohammad Mehdi Ghahremanpour, 
 *             Paul J. van Maaren, 
 *             David van der Spoel (Project leader)
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
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, 
 * Boston, MA  02110-1301, USA.
 */
 
/*! \internal \brief
 * Implements part of the alexandria program.
 * \author Mohammad Mehdi Ghahremanpour <mohammad.ghahremanpour@icm.uu.se>
 * \author David van der Spoel <david.vanderspoel@icm.uu.se>
 */
 
 
#ifndef GMX_SIMPLE_COMM_H
#define GMX_SIMPLE_COMM_H

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>

#include <algorithm>
#include <string>
#include <vector>

#include "gromacs/gmxlib/network.h"
#include "gromacs/mdtypes/commrec.h"

void gmx_send(const t_commrec *cr, int dest, void *buf, int bufsize);

void gmx_recv(const t_commrec *cr, int src, void *buf, int bufsize);

void gmx_send_str(const t_commrec *cr, int dest, const std::string *str);

void gmx_recv_str(const t_commrec *cr, int src, std::string *str);

void gmx_send_double(const t_commrec *cr, int dest, double d);

double gmx_recv_double(const t_commrec *cr, int src);

void gmx_send_int(const t_commrec *cr, int dest, int d);

int gmx_recv_int(const t_commrec *cr, int src);

#endif
