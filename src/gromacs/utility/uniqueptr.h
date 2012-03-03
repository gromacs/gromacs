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
/*! \libinternal \file
 * \brief
 * Declares functions for OS-independent path handling.
 *
 * \author Roland Schulz <roland@utk.edu>
 * \author John Eblen <jeblen@acm.org>
 * \inlibraryapi
 */

#ifndef GMX_UTILITY_UNIQUEPTR_H
#define GMX_UTILITY_UNIQUEPTR_H

#ifdef HAVE_CXX11 // C++11 Compiler
#include <memory>
#include <algorithm>
#else      // C++03 Compiler
#include <boost/shared_ptr.hpp>
#endif

namespace gmx
{

#ifdef HAVE_CXX11 // C++11 Compiler
using std::move;
template<typename T>
struct gmx_unique_ptr {
        typedef std::unique_ptr<T> type;
};
#else // C++03 Compiler
template<typename T>
const boost::shared_ptr<T> &move(const boost::shared_ptr<T> &ptr) {
        return ptr;
}
template<typename T>
boost::shared_ptr<T> &move(boost::shared_ptr<T> &ptr) {
                return ptr;
}
template<typename T>
struct gmx_unique_ptr {
        typedef boost::shared_ptr<T> type;
};
#endif

}

#endif
