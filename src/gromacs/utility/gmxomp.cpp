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

std::optional<std::string> messageWhenOpenMPLibraryWillSetAffinity()
{
    std::optional<std::string> message;
#if GMX_OPENMP
    // The OMP_PROC_BIND is the standard environment variable used by
    // OpenMP implementations. The KMP_AFFINITY environment variable
    // is used by Intel's compiler, and GOMP_CPU_AFFINITY by the GNU
    // compilers (Intel also honors it as well).
    std::string environmentVariablesFoundToBeSet;
    for (const char* nameOfEnvironmentVariable :
         { "OMP_PROC_BIND", "GOMP_CPU_AFFINITY", "KMP_AFFINITY" })
    {
        const char* const environmentVariableValue = std::getenv(nameOfEnvironmentVariable);
        if (environmentVariableValue != nullptr && *environmentVariableValue != '\0')
        {
            if (!environmentVariablesFoundToBeSet.empty())
            {
                environmentVariablesFoundToBeSet += ", ";
            }
            environmentVariablesFoundToBeSet += nameOfEnvironmentVariable;
        }
    }
    if (!environmentVariablesFoundToBeSet.empty())
    {
        try
        {
            const char* programName = gmx::getProgramContext().displayName();
            message                 = gmx::formatString(
                    "NOTE: Thread-affinity variable(s) set, will turn off %s internal affinity\n"
                                    "      setting as the two can conflict and cause performance degradation.\n"
                                    "      To keep using the %s internal affinity setting, unset the\n"
                                    "      following environment variable(s):\n"
                                    "      %s",
                    programName,
                    programName,
                    environmentVariablesFoundToBeSet.c_str());
        }
        GMX_CATCH_ALL_AND_EXIT_WITH_FATAL_ERROR
    }
#endif /* GMX_OPENMP */
    return message;
}

namespace gmx
{

std::string openmpDescription()
{
    return GMX_OPENMP ? formatString("enabled (GMX_OPENMP_MAX_THREADS = %d)", GMX_OPENMP_MAX_THREADS)
                      : "disabled";
}

} // namespace gmx
