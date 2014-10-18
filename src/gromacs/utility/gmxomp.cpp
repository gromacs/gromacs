/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2012,2013,2014, by the GROMACS development team, led by
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
/*! \internal \file
 * \brief
 * Implements functions from gmxomp.h.
 *
 * \ingroup module_utility
 */
#include "gmxpre.h"

#include "gmxomp.h"

#include "config.h"

#include <stdio.h>
#include <stdlib.h>

#ifdef GMX_OPENMP
#include <omp.h>
#endif

#include "gromacs/utility/basedefinitions.h"
#include "gromacs/utility/cstringutil.h"
#include "gromacs/utility/exceptions.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/programcontext.h"
#include "gromacs/utility/stringutil.h"

int gmx_omp_get_max_threads(void)
{
#ifdef GMX_OPENMP
    return omp_get_max_threads();
#else
    return 1;
#endif
}

int gmx_omp_get_num_procs(void)
{
#ifdef GMX_OPENMP
    return omp_get_num_procs();
#else
    return 1;
#endif
}

int gmx_omp_get_thread_num(void)
{
#ifdef GMX_OPENMP
    return omp_get_thread_num();
#else
    return 0;
#endif
}

void gmx_omp_set_num_threads(int num_threads)
{
#ifdef GMX_OPENMP
    omp_set_num_threads(num_threads);
#else
    GMX_UNUSED_VALUE(num_threads);
#endif
}

gmx_bool gmx_omp_check_thread_affinity(char **message)
{
    bool shouldSetAffinity = true;

    *message = NULL;
#ifdef GMX_OPENMP
    /* We assume that the affinity setting is available on all platforms
     * gcc supports. Even if this is not the case (e.g. Mac OS) the user
     * will only get a warning. */
#if defined(__GNUC__) || defined(__INTEL_COMPILER)
    const char *programName;
    try
    {
        programName = gmx::getProgramContext().displayName();
    }
    GMX_CATCH_ALL_AND_EXIT_WITH_FATAL_ERROR;

    const char *const gomp_env            = getenv("GOMP_CPU_AFFINITY");
    const bool        bGompCpuAffinitySet = (gomp_env != NULL);

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
                        programName, programName);
            *message = gmx_strdup(buf.c_str());
        }
        GMX_CATCH_ALL_AND_EXIT_WITH_FATAL_ERROR;
        shouldSetAffinity = false;
    }
#endif /* __GNUC__ || __INTEL_COMPILER */

#if defined(__INTEL_COMPILER)
    const char *const kmp_env         = getenv("KMP_AFFINITY");
    const bool        bKmpAffinitySet = (kmp_env != NULL);

    /* disable Intel OpenMP affinity if neither KMP_AFFINITY nor
     * GOMP_CPU_AFFINITY is set (Intel uses the GNU env. var as well) */
    if (!bKmpAffinitySet && !bGompCpuAffinitySet)
    {
        int retval;

#ifdef _MSC_VER
        /* Windows not POSIX */
        retval = _putenv_s("KMP_AFFINITY", "disabled");
#else
        /* POSIX */
        retval = setenv("KMP_AFFINITY", "disabled", 0);
#endif  /* _MSC_VER */

        if (debug)
        {
            fprintf(debug, "Disabling Intel OpenMP affinity by setting the KMP_AFFINITY=disabled env. var.\n");
        }

        if (retval != 0)
        {
            gmx_warning("Disabling Intel OpenMp affinity setting failed!");
        }
    }

    /* turn off internal pinning KMP_AFFINITY != "disabled" */
    if (bKmpAffinitySet && (gmx_strncasecmp(kmp_env, "disabled", 8) != 0))
    {
        try
        {
            std::string buf = gmx::formatString(
                        "NOTE: KMP_AFFINITY set, will turn off %s internal affinity\n"
                        "      setting as the two can conflict and cause performance degradation.\n"
                        "      To keep using the %s internal affinity setting, set the\n"
                        "      KMP_AFFINITY=disabled environment variable.",
                        programName, programName);
            *message = gmx_strdup(buf.c_str());
        }
        GMX_CATCH_ALL_AND_EXIT_WITH_FATAL_ERROR;
        shouldSetAffinity = false;
    }
#endif /* __INTEL_COMPILER */

#endif /* GMX_OPENMP */
    return shouldSetAffinity;
}
