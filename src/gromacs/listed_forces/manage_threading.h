/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2004, The GROMACS development team.
 * Copyright (c) 2013,2014,2015,2016,2018,2019, by the GROMACS development team, led by
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
 * \brief Declares functions for managing threading of listed forces
 *
 * \author Mark Abraham <mark.j.abraham@gmail.com>
 * \inlibraryapi
 * \ingroup module_listed_forces
 */
#ifndef GMX_LISTED_FORCES_MANAGE_THREADING_H
#define GMX_LISTED_FORCES_MANAGE_THREADING_H

#include <cstdio>

struct bonded_threading_t;
struct t_idef;

/*! \brief Divide the listed interactions over the threads and GPU
 *
 * Uses fr->nthreads for the number of threads, and sets up the
 * thread-force buffer reduction.
 * This should be called each time the bonded setup changes;
 * i.e. at start-up without domain decomposition and at DD.
 */
void setup_bonded_threading(bonded_threading_t* bt, int numAtoms, bool useGpuForBondes, const t_idef& idef);

//! Destructor.
void tear_down_bonded_threading(bonded_threading_t* bt);

/*! \brief Initialize the bonded threading data structures
 *
 * Allocates and initializes a bonded threading data structure.
 * A pointer to this struct is returned as \p *bb_ptr.
 *
 * \todo Avoid explicit pointers by using Impl
 */
bonded_threading_t* init_bonded_threading(FILE* fplog, int nenergrp);

#endif
