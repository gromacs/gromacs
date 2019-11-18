/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2019, by the GROMACS development team, led by
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
 * \brief Declares object that provides buffers for the copies
 * necessary when MPI_IN_PLACE is not supported by the MPI library.
 *
 * The struct is declared in the header so that functionality it
 * supports can be inlined.
 *
 * \author Mark Abraham <mark.j.abraham@gmail.com>
 * \inlibraryapi
 * \ingroup module_utility
 */
#ifndef GMX_UTILTIY_MPIINPLACEBUFFERS_H
#define GMX_UTILTIY_MPIINPLACEBUFFERS_H

#include <cstdint>

struct mpi_in_place_buf_t
{
    /* these buffers are used as destination buffers if MPI_IN_PLACE isn't
       supported.*/
    int* ibuf; /* for ints */
    int  ibuf_alloc;

    int64_t* libuf;
    int      libuf_alloc;

    float* fbuf; /* for floats */
    int    fbuf_alloc;

    double* dbuf; /* for doubles */
    int     dbuf_alloc;
};

//! Cleans up the buffers
void done_mpi_in_place_buf(mpi_in_place_buf_t* buf);

#endif
