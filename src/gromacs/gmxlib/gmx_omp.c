/* -*- mode: c; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; c-file-style: "stroustrup"; -*-
 *
 *
 *                This source code is part of
 *
 *                 G   R   O   M   A   C   S
 *
 *          GROningen MAchine for Chemical Simulations
 *
 * Written by the Gromacs development team under coordination of
 * David van der Spoel, Berk Hess, and Erik Lindahl.
 *
 * This library is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public License
 * as published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 *
 * To help us fund GROMACS development, we humbly ask that you cite
 * the research papers on the package. Check out http://www.gromacs.org
 *
 * And Hey:
 * Gnomes, ROck Monsters And Chili Sauce
 */

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#ifdef GMX_OPENMP
#include <omp.h>
#endif

#include <stdio.h>

#include "md_logging.h"
#include "gmx_fatal.h"
#include "statutil.h"
#include "string2.h"
#include "gmx_omp.h"

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
    return;
#endif
}

/*!
 * Thread affinity set by the OpenMP library can conflict with the GROMACS
 * internal affinity setting.
 *
 * While GNU OpenMP does not set affinity by default, the Intel OpenMP library
 * does. This conflicts with the internal affinity (especially thread-MPI)
 * setting, results in incorrectly locked threads, and causes dreadful performance.
 *
 * The KMP_AFFINITY environment variable is used by Intel, GOMP_CPU_AFFINITY
 * by the GNU compilers (Intel also honors it well). If any of the variables
 * is set, we honor it, disable the internal pinning, and warn the user.
 * When using Intel OpenMP, we will disable affinity if the user did not set it
 * anually through one of the aforementioned environment variables.
 *
 * Note that the Intel OpenMP affinity disabling iwll only take effect if this
 * function is called before the OpenMP library gets initialized which happens
 * when the first call is made into a compilation unit that contains OpenMP
 * pragmas.
 */
void gmx_omp_check_thread_affinity(FILE *fplog, const t_commrec *cr,
                                   gmx_hw_opt_t *hw_opt)
{
    gmx_bool bKmpAffinitySet, bGompCpuAffinitySet;
    char    *kmp_env, *gomp_env;

    /* no need to worry if internal thread pinning is turned off */
    if (hw_opt->thread_affinity == threadaffOFF)
    {
        return;
    }

#if defined(GMX_OPENMP)

    /* We assume that the affinity setting is available on all platforms
     * gcc supports. Even if this is not the case (e.g. Mac OS) the user
     * will only get a warning.*/
    bGompCpuAffinitySet = FALSE;
    gomp_env            = NULL;
#if defined(__GNUC__)
    gomp_env            = getenv("GOMP_CPU_AFFINITY");
    bGompCpuAffinitySet = (gomp_env != NULL);
#endif /* __GNUC__ */

    bKmpAffinitySet = FALSE;
#if defined(__INTEL_COMPILER)
    kmp_env         = getenv("KMP_AFFINITY");
    bKmpAffinitySet = (kmp_env != NULL);

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
        /* TODO: with -pin auto we should only warn when using all cores */
        md_print_warn(cr, fplog,
                      "NOTE: KMP_AFFINITY set, will turn off %s internal affinity\n"
                      "      setting as the two can conflict and cause performance degradation.\n"
                      "      To keep using the %s internal affinity setting, set the\n"
                      "      KMP_AFFINITY=disabled environment variable.",
                      ShortProgram(), ShortProgram());

        hw_opt->thread_affinity = threadaffOFF;
    }
#endif /* __INTEL_COMPILER */

#if defined(__INTEL_COMPILER) || defined(__GNUC__)
    /* turn off internal pinning f GOMP_CPU_AFFINITY is set & non-empty */
    if (bGompCpuAffinitySet && gomp_env != NULL && gomp_env != '\0')
    {
        /* TODO: with -pin auto we should only warn when using all cores */
        md_print_warn(cr, fplog,
                      "NOTE: GOMP_CPU_AFFINITY set, will turn off %s internal affinity\n"
                      "      setting as the two can conflict and cause performance degradation.\n"
                      "      To keep using the %s internal affinity setting, unset the\n"
                      "      GOMP_CPU_AFFINITY environment variable.",
                      ShortProgram(), ShortProgram());

        hw_opt->thread_affinity = threadaffOFF;
    }
#endif /* __INTEL_COMPILER || __GNUC__ */

#endif /* GMX_OPENMP */
}
