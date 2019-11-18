/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2014,2015,2018,2019, by the GROMACS development team, led by
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
#include "gmxpre.h"

#include "baseversion.h"

#include "config.h"

#include "baseversion_gen.h"

const char* gmx_version()
{
    return _gmx_ver_string;
}

const char* gmx_version_git_full_hash()
{
    return _gmx_full_git_hash;
}

const char* gmx_version_git_central_base_hash()
{
    return _gmx_central_base_hash;
}

const char* gmxDOI()
{
    return gmxSourceDoiString;
}

const char* gmxReleaseSourceChecksum()
{
    return gmxReleaseSourceFileChecksum;
}

const char* gmxCurrentSourceChecksum()
{
    return gmxCurrentSourceFileChecksum;
}

#if GMX_DOUBLE
void gmx_is_double_precision() {}
#else
void gmx_is_single_precision() {}
#endif

/* Note that this array (and some which follow) must match the "GPU
 * support enumeration" in src/config.h.cmakein */
static const char* const gpuImplementationStrings[] = { "disabled", "CUDA", "OpenCL" };

const char* getGpuImplementationString()
{
    return gpuImplementationStrings[GMX_GPU];
}
