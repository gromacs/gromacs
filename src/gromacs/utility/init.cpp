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
 *
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
 * Implements gmx::init
 *
 * \author Ryan Johnson <ryanphjohnson@gmail.com>
 */
#ifdef HAVE_CONFIG_H
#include "config.h"
#endif
#include "external/mpp/gmxmpp.h"

#include "gromacs/utility/init.h"

namespace gmx
{
namespace
{
/*! \brief
 * Does nothing but call constructor and destructor functions.
 */
class GlobalInitializer
{
public:
    GlobalInitializer (int &argc, char** &argv)
    {
#ifdef GMX_LIB_MPI
        mpi::init (&argc, &argv);
#endif
    }
    ~GlobalInitializer()
    {
#ifdef GMX_LIB_MPI
        mpi::finalize();
#endif
    }
};
} // end namespace
const ProgramInfo &init(const char *realBinaryName, int& argc, char** &argv)
{
    static GlobalInitializer initSingleton(argc, argv);
    return ProgramInfo::init(realBinaryName, argc, argv);
}
} // end gmx namespace
