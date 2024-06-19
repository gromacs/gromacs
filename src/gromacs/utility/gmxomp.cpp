/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2012- The GROMACS Authors
 * and the project initiators Erik Lindahl, Berk Hess and David van der Spoel.
 * Consult the AUTHORS/COPYING files and https://www.gromacs.org for details.
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
 * https://www.gnu.org/licenses, or write to the Free Software Foundation,
 * Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA.
 *
 * If you want to redistribute modifications to GROMACS, please
 * consider that scientific software is very special. Version
 * control is crucial - bugs must be traceable. We will be happy to
 * consider code for inclusion in the official distribution, but
 * derived work must not be called official GROMACS. Details are found
 * in the README & COPYING files - if they are missing, get the
 * official version at https://www.gromacs.org.
 *
 * To help us fund GROMACS development, we humbly ask that you cite
 * the research papers on the package. Check out https://www.gromacs.org.
 */
/*! \internal \file
 * \brief
 * Implements functions from gmxomp.h.
 *
 * \ingroup module_utility
 */
#include "gmxpre.h"

#include "gromacs/utility/gmxomp.h"

#include "config.h"

#include <cstdio>
#include <cstdlib>

#include <string>

#if GMX_OPENMP
#    include <omp.h>
#endif

#include "gromacs/utility/basedefinitions.h"
#include "gromacs/utility/cstringutil.h"
#include "gromacs/utility/exceptions.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/programcontext.h"
#include "gromacs/utility/stringutil.h"

int gmx_omp_get_max_threads()
{
#if GMX_OPENMP
    return omp_get_max_threads();
#else
    return 1;
#endif
}

int gmx_omp_get_num_procs()
{
#if GMX_OPENMP
    return omp_get_num_procs();
#else
    return 1;
#endif
}

int gmx_omp_get_thread_num()
{
#if GMX_OPENMP
    return omp_get_thread_num();
#else
    return 0;
#endif
}

void gmx_omp_set_num_threads(int num_threads)
{
#if GMX_OPENMP
    omp_set_num_threads(num_threads);
#else
    GMX_UNUSED_VALUE(num_threads);
#endif
}

bool gmx_omp_check_thread_affinity(char** message)
{
    bool shouldSetAffinity = true;

    *message = nullptr;
#if GMX_OPENMP
    /* We assume that the affinity setting is available on all platforms
     * gcc supports. Even if this is not the case (e.g. Mac OS) the user
     * will only get a warning. */
#    if defined(__GNUC__)
    const char* programName;
    try
    {
        programName = gmx::getProgramContext().displayName();
    }
    GMX_CATCH_ALL_AND_EXIT_WITH_FATAL_ERROR

    const char* const gomp_env            = getenv("GOMP_CPU_AFFINITY");
    const bool        bGompCpuAffinitySet = (gomp_env != nullptr);

    /* turn off internal pinning if GOMP_CPU_AFFINITY is set & non-empty */
    if (bGompCpuAffinitySet && *gomp_env != '\0')
    {
        try
        {
            std::string buf = gmx::formatString(
                    "NOTE: GOMP_CPU_AFFINITY set, will turn off %s internal affinity\n"
                    "      setting as the two can conflict and cause performance degradation.\n"
                    "      To keep using the %s internal affinity setting, unset the\n"
                    "      GOMP_CPU_AFFINITY environment variable.",
                    programName,
                    programName);
            *message = gmx_strdup(buf.c_str());
        }
        GMX_CATCH_ALL_AND_EXIT_WITH_FATAL_ERROR
        shouldSetAffinity = false;
    }
#    endif /* __GNUC__ */

#endif /* GMX_OPENMP */
    return shouldSetAffinity;
}
