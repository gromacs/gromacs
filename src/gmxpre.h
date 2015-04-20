/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2014, by the GROMACS development team, led by
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
 * \brief
 * Prerequisite header file for \Gromacs build.
 *
 * This header should be included as the first header in all source files, but
 * not in header files.  It is intended to contain definitions that must appear
 * before any other code to work properly (e.g., macro definitions that
 * influence behavior of system headers).  This frees other code from include
 * order dependencies that may raise from requirements of getting these
 * definitions from some header.
 *
 * The definitions here should be kept to a minimum, and should be as static as
 * possible (typically not change as a result of user choices in the build
 * system), as any change will trigger a full rebuild.  Avoid including any
 * actual headers to not hide problems with include-what-you-use, and to keep
 * build times to a minimum.  Also, installer headers should avoid relying on
 * the definitions from here (if possible), as this header will not be
 * available to the user.
 *
 * \inlibraryapi
 */
//! \cond
#ifdef HAVE_CONFIG_H
#include "gmxpre-config.h"
#endif

/* We use a few GNU functions for thread affinity and other low-level stuff.
 * However, all such uses should be accompanied by #ifdefs and a feature test
 * at CMake level, so that the actual uses will be compiled only when available.
 * But since the define affects system headers, it should be defined before
 * including any system headers, and this is a robust location to do that.
 * If this were defined only in source files that needed it, it would clutter
 * the list of includes somewhere close to the beginning and make automatic
 * sorting of the includes more difficult.
 */
#ifndef _GNU_SOURCE
#define _GNU_SOURCE 1
#endif

/* Some C++(?) compilers require these to be defined to get the integer limits
 * and format specifier macros from stdint.h and inttypes.h, respectively.
 * The macros are in C99 and C++11, but not in C++98...
 * As with _GNU_SOURCE, these need to be defined before these headers get first
 * included.  Unlike _GNU_SOURCE, these headers are included indirectly in most
 * header and source files (even though the macros are not used that often), so
 * there is no easy alternative to defining them here, either.
 * If someone happens to use such a compiler to compile against the installed
 * Gromacs headers, they need for now take care to define the macros themselves
 * (as there is no way Gromacs can do that consistently).
 */
#define __STDC_LIMIT_MACROS
#define __STDC_FORMAT_MACROS

#ifdef GMX_FAHCORE
#define FULLINDIRECT 1
#define USE_FAH_XDR  1
#include "swindirect.h"
#endif
//! \endcond
