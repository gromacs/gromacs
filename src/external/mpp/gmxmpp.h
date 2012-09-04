/*
 *
 *                This source code is part of
 *
 *                 G   R   O   M   A   C   S
 *
 *          GROningen MAchine for Chemical Simulations
 *
 * Written by David van der Spoel, Erik Lindahl, Berk Hess, and others.
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2009, The GROMACS development team,
 * check out http://www.gromacs.org for more information.

 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 *
 * If you want to redistribute modifications, please consider that
 * scientific software is very special. Version control is crucial -
 * bugs must be traceable. We will be happy to consider code for
 * inclusion in the official distribution, but derived work must not
 * be called official GROMACS. Details are found in the README & COPYING
 * files - if they are missing, get the official version at www.gromacs.org.
 *
 * To help us fund GROMACS development, we humbly ask that you cite
 * the papers on the package - you can find them in the top README file.
 *
 * For more info, check our website at http://www.gromacs.org
 */
#pragma once

#include "gromacs/utility/gmx_header_config.h"
// For GMX_LIB_MPI
#ifdef GMX_CXX11
#define MPP_CXX11_RVALREF
#endif
#ifdef GMX_LIB_MPI
#include "mpp.h"
#else
#ifdef GMX_THREAD_MPI
#include <tmpi.h>
#else
typedef void* MPI_Datatype;
extern const MPI_Datatype MPI_DOUBLE;
extern const MPI_Datatype MPI_INT;
extern const MPI_Datatype MPI_CHAR;
extern const MPI_Datatype MPI_FLOAT;
extern const MPI_Datatype MPI_LONG;
extern const MPI_Datatype MPI_BYTE;
int MPI_Type_commit(MPI_Datatype *datatype);
#endif //GMX_THREAD_MPI
#include <cstddef>
typedef std::size_t MPI_Aint;
extern const MPI_Datatype MPI_UNSIGNED_LONG;
extern const MPI_Datatype MPI_DATATYPE_NULL;
int MPI_Get_address(void*,MPI_Aint*);
int MPI_Type_indexed(int, int*, int*, MPI_Datatype, MPI_Datatype*);
int MPI_Type_free(MPI_Datatype *datatype);
int MPI_Type_dup(MPI_Datatype type, MPI_Datatype *newtype);
#define MPP_NO_MPI_INCL
#include "type_traits.h"
#endif //GMX_LIB_MPI
