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
/*! \internal \file
 * \brief
 * An MPI CPP Interface
 * \authors Ryan Johnson <ryanphjohnson@gmail.com>
 */

// For GMX_LIB_MPI
#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include "mpp.h"

namespace mpi {
comm comm::world = comm(MPI_COMM_NULL);
std::vector<MPI_Datatype> cached_types;
// returns MPI_Datatype containing all datatypes that were added to the
//     mpi_type_builder. If none were added, it makes an empty MPI_Datatype.
MPI_Datatype mpi_type_builder::build() {
    //TODO: another possible optimization would be for is_pod types, which should be completely send
    //only OK on non-heterogenous hardware (so that MPI_BYTE copy is OK)
       MPI_Datatype dt;
       if (type.empty()) { MPI_Type_contiguous (0, MPI_INT, &dt); }
       else {
           std::vector<int> size(type.size(),1);
           MPI_Type_create_struct (type.size(), &size.front(), &addr.front(), &type.front(), &dt);
           if (base!=0 && curr_extent!=real_extent) {
               MPI_Datatype dt_resized;
               MPI_Type_commit(&dt);
               MPI_Type_create_resized(dt,0,real_extent,&dt_resized);
               MPI_Type_free(&dt);
               dt=dt_resized;
           }
       }
       MPI_Type_commit (&dt);
       std::vector<MPI_Datatype>::iterator it;
       std::vector<bool>::iterator itf;
       for (it = type.begin(), itf = bFree.begin(); it!=type.end(); ++it, ++itf) {
           if (*itf) {
               MPI_Type_free(&*it);
           }
       }
       return dt;
}
} // mpi

