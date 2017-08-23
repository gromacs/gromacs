/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2009,2010,2011,2012,2013,2014,2015,2016,2017, by the GROMACS development team, led by
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

#include "config.h"

#include <cmath>
#include <cstdlib>
#include <cstring>
#include <ctime>

#include <algorithm>
#include <string>

#ifdef HAVE_SYS_TIME_H
#include <sys/time.h>
#endif

#include "gromacs/commandline/pargs.h"
#include "gromacs/ewald/pme.h"
#include "gromacs/fft/calcgrid.h"
#include "gromacs/fileio/checkpoint.h"
#include "gromacs/fileio/tpxio.h"
#include "gromacs/gmxana/gmx_ana.h"
#include "gromacs/math/utilities.h"
#include "gromacs/math/vec.h"
#include "gromacs/mdlib/perf_est.h"
#include "gromacs/mdtypes/commrec.h"
#include "gromacs/mdtypes/inputrec.h"
#include "gromacs/mdtypes/md_enums.h"
#include "gromacs/mdtypes/state.h"
#include "gromacs/pbcutil/pbc.h"
#include "gromacs/taskassignment/usergpuids.h"
#include "gromacs/timing/walltime_accounting.h"
#include "gromacs/topology/topology.h"
#include "gromacs/utility/arraysize.h"
#include "gromacs/utility/baseversion.h"
#include "gromacs/utility/cstringutil.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/futil.h"
#include "gromacs/utility/gmxassert.h"
#include "gromacs/utility/smalloc.h"
#include "gromacs/utility/stringutil.h"

/* Enum for situations that can occur during log file parsing, the
 * corresponding string entries can be found in do_the_tests() in
 * const char* ParseLog[] */
/* TODO clean up CamelCasing of these enum names */
enum {
    eParselogOK,
    eParselogNotFound,
    eParselogNoPerfData,
    eParselogTerm,
    eParselogResetProblem,
    eParselogNoDDGrid,
    eParselogTPXVersion,
    eParselogNotParallel,
    eParselogLargePrimeFactor,
    eParselogMismatchOfNumberOfPPRanksAndAvailableGPUs,
    eParselogGpuProblem,
    eParselogFatal,
    eParselogNr
};


typedef struct
{
    int     nPMEnodes;    /* number of PME-only nodes used in this test */
    int     nx, ny, nz;   /* DD grid */
    int     guessPME;     /* if nPMEnodes == -1, this is the guessed number of PME nodes */
    double *Gcycles;      /* This can contain more than one value if doing multiple tests */
    double  Gcycles_Av;
    float  *ns_per_day;
    float   ns_per_day_Av;
    float  *PME_f_load;     /* PME mesh/force load average*/
    float   PME_f_load_Av;  /* Average average ;) ... */
    char   *mdrun_cmd_line; /* Mdrun command line used for this test */
} t_perf;


typedef struct
{
    int             nr_inputfiles;   /* The number of tpr and mdp input files */
    gmx_int64_t     orig_sim_steps;  /* Number of steps to be done in the real simulation */
    gmx_int64_t     orig_init_step;  /* Init step for the real simulation */
    real           *rcoulomb;        /* The coulomb radii [0...nr_inputfiles] */
    real           *rvdw;            /* The vdW radii */
    real           *rlist;           /* Neighbourlist cutoff radius */
    int            *nkx, *nky, *nkz;
    real           *fsx, *fsy, *fsz; /* Fourierspacing in x,y,z dimension */
} t_inputinfo;


static void sep_line(FILE *fp)
{
    fprintf(fp, "\n------------------------------------------------------------\n");
}


/* Wrapper for system calls */
static int gmx_system_call(char *command)
{
    return ( system(command) );
}


/* Check if string starts with substring */
static gmx_bool str_starts(const char *string, const char *substring)
{
    return ( std::strncmp(string, substring, std::strlen(substring)) == 0);
}


static void cleandata(t_perf *perfdata, int test_nr)
{
    perfdata->Gcycles[test_nr]    = 0.0;
    perfdata->ns_per_day[test_nr] = 0.0;
    perfdata->PME_f_load[test_nr] = 0.0;

    return;
}


static void remove_if_exists(const char *fn)
{
    if (gmx_fexist(fn))
    {
        fprintf(stdout, "Deleting %s\n", fn);
        remove(fn);
    }
}


static void finalize(const char *fn_out)
{
    char  buf[STRLEN];
    FILE *fp;


    fp = fopen(fn_out, "r");
    fprintf(stdout, "\n\n");

    while (fgets(buf, STRLEN-1, fp) != nullptr)
    {
        fprintf(stdout, "%s", buf);
    }
    fclose(fp);
    fprintf(stdout, "\n\n");
}


enum {
    eFoundNothing, eFoundDDStr, eFoundAccountingStr, eFoundCycleStr
};

static int parse_logfile(const char *logfile, const char *errfile,
                         t_perf *perfdata, int test_nr, int presteps, gmx_int64_t cpt_steps,
                         int nnodes)
{
    FILE           *fp;
    char            line[STRLEN], dumstring[STRLEN], dumstring2[STRLEN];
    const char      matchstrdd[]  = "Domain decomposition grid";
    const char      matchstrcr[]  = "resetting all time and cycle counters";
    const char      matchstrbal[] = "Average PME mesh/force load:";
    const char      matchstring[] = "R E A L   C Y C L E   A N D   T I M E   A C C O U N T I N G";
    const char      errSIG[]      = "signal, stopping at the next";
    int             iFound;
    float           dum1, dum2, dum3, dum4;
    int             ndum;
    int             npme;
    gmx_int64_t     resetsteps     = -1;
    gmx_bool        bFoundResetStr = FALSE;
    gmx_bool        bResetChecked  = FALSE;


    if (!gmx_fexist(logfile))
    {
        fprintf(stderr, "WARNING: Could not find logfile %s.\n", logfile);
        cleandata(perfdata, test_nr);
        return eParselogNotFound;
    }

    fp = fopen(logfile, "r");
    perfdata->PME_f_load[test_nr] = -1.0;
    perfdata->guessPME            = -1;

    iFound = eFoundNothing;
    if (1 == nnodes)
    {
        iFound = eFoundDDStr; /* Skip some case statements */
    }

    while (fgets(line, STRLEN, fp) != nullptr)
    {
        /* Remove leading spaces */
        ltrim(line);

        /* Check for TERM and INT signals from user: */
        if (std::strstr(line, errSIG) != nullptr)
        {
            fclose(fp);
            cleandata(perfdata, test_nr);
            return eParselogTerm;
        }

        /* Check whether cycle resetting  worked */
        if (presteps > 0 && !bFoundResetStr)
        {
            if (std::strstr(line, matchstrcr) != nullptr)
            {
                sprintf(dumstring, "step %s", "%" GMX_SCNd64);
                sscanf(line, dumstring, &resetsteps);
                bFoundResetStr = TRUE;
                if (resetsteps == presteps+cpt_steps)
                {
                    bResetChecked = TRUE;
                }
                else
                {
                    sprintf(dumstring, "%" GMX_PRId64, resetsteps);
                    sprintf(dumstring2, "%" GMX_PRId64, presteps+cpt_steps);
                    fprintf(stderr, "WARNING: Time step counters were reset at step %s,\n"
                            "         though they were supposed to be reset at step %s!\n",
                            dumstring, dumstring2);
                }
            }
        }

        /* Look for strings that appear in a certain order in the log file: */
        switch (iFound)
        {
            case eFoundNothing:
                /* Look for domain decomp grid and separate PME nodes: */
                if (str_starts(line, matchstrdd))
                {
                    sscanf(line, "Domain decomposition grid %d x %d x %d, separate PME ranks %d",
                           &(perfdata->nx), &(perfdata->ny), &(perfdata->nz), &npme);
                    if (perfdata->nPMEnodes == -1)
                    {
                        perfdata->guessPME = npme;
                    }
                    else if (perfdata->nPMEnodes != npme)
                    {
                        gmx_fatal(FARGS, "PME ranks from command line and output file are not identical");
                    }
                    iFound = eFoundDDStr;
                }
                /* Catch a few errors that might have occurred: */
                else if (str_starts(line, "There is no domain decomposition for"))
                {
                    fclose(fp);
                    return eParselogNoDDGrid;
                }
                else if (str_starts(line, "The number of ranks you selected"))
                {
                    fclose(fp);
                    return eParselogLargePrimeFactor;
                }
                else if (str_starts(line, "reading tpx file"))
                {
                    fclose(fp);
                    return eParselogTPXVersion;
                }
                else if (str_starts(line, "The -dd or -npme option request a parallel simulation"))
                {
                    fclose(fp);
                    return eParselogNotParallel;
                }
                break;
            case eFoundDDStr:
                /* Even after the "Domain decomposition grid" string was found,
                 * it could be that mdrun had to quit due to some error. */
                if (str_starts(line, "Incorrect launch configuration: mismatching number of"))
                {
                    fclose(fp);
                    return eParselogMismatchOfNumberOfPPRanksAndAvailableGPUs;
                }
                else if (str_starts(line, "Some of the requested GPUs do not exist"))
                {
                    fclose(fp);
                    return eParselogGpuProblem;
                }
                /* Look for PME mesh/force balance (not necessarily present, though) */
                else if (str_starts(line, matchstrbal))
                {
                    sscanf(&line[std::strlen(matchstrbal)], "%f", &(perfdata->PME_f_load[test_nr]));
                }
                /* Look for matchstring */
                else if (str_starts(line, matchstring))
                {
                    iFound = eFoundAccountingStr;
                }
                break;
            case eFoundAccountingStr:
                /* Already found matchstring - look for cycle data */
                if (str_starts(line, "Total  "))
                {
                    sscanf(line, "Total %*f %lf", &(perfdata->Gcycles[test_nr]));
                    iFound = eFoundCycleStr;
                }
                break;
            case eFoundCycleStr:
                /* Already found cycle data - look for remaining performance info and return */
                if (str_starts(line, "Performance:"))
                {
                    ndum = sscanf(line, "%s %f %f %f %f", dumstring, &dum1, &dum2, &dum3, &dum4);
                    /* (ns/day) is the second last entry, depending on whether GMX_DETAILED_PERF_STATS was set in print_perf(), nrnb.c */
                    perfdata->ns_per_day[test_nr] = (ndum == 5) ? dum3 : dum1;
                    fclose(fp);
                    if (bResetChecked || presteps == 0)
                    {
                        return eParselogOK;
                    }
                    else
                    {
                        return eParselogResetProblem;
                    }
                }
                break;
        }
    } /* while */

    /* Close the log file */
    fclose(fp);

    /* Check why there is no performance data in the log file.
     * Did a fatal errors occur? */
    if (gmx_fexist(errfile))
    {
        fp = fopen(errfile, "r");
        while (fgets(line, STRLEN, fp) != nullptr)
        {
            if (str_starts(line, "Fatal error:") )
            {
                if (fgets(line, STRLEN, fp) != nullptr)
                {
                    fprintf(stderr, "\nWARNING: An error occurred during this benchmark:\n"
                            "%s\n", line);
                }
                fclose(fp);
                cleandata(perfdata, test_nr);
                return eParselogFatal;
            }
        }
        fclose(fp);
    }
    else
    {
        fprintf(stderr, "WARNING: Could not find stderr file %s.\n", errfile);
    }

    /* Giving up ... we could not find out why there is no performance data in
     * the log file. */
    fprintf(stdout, "No performance data in log file.\n");
    cleandata(perfdata, test_nr);

    return eParselogNoPerfData;
}


static gmx_bool analyze_data(
        FILE         *fp,
        const char   *fn,
        t_perf      **perfdata,
        int           nnodes,
        int           ntprs,
        int           ntests,
        int           nrepeats,
        t_inputinfo  *info,
        int          *index_tpr,    /* OUT: Nr of mdp file with best settings */
        int          *npme_optimal) /* OUT: Optimal number of PME nodes */
{
    int      i, j, k;
    int      line  = 0, line_win = -1;
    int      k_win = -1, i_win = -1, winPME;
    double   s     = 0.0; /* standard deviation */
    t_perf  *pd;
    char     strbuf[STRLEN];
    char     str_PME_f_load[13];
    gmx_bool bCanUseOrigTPR;
    gmx_bool bRefinedCoul, bRefinedVdW, bRefinedGrid;


    if (nrepeats > 1)
    {
        sep_line(fp);
        fprintf(fp, "Summary of successful runs:\n");
        fprintf(fp, "Line tpr PME ranks  Gcycles Av.     Std.dev.       ns/day        PME/f");
        if (nnodes > 1)
        {
            fprintf(fp, "    DD grid");
        }
        fprintf(fp, "\n");
    }


    for (k = 0; k < ntprs; k++)
    {
        for (i = 0; i < ntests; i++)
        {
            /* Select the right dataset: */
            pd = &(perfdata[k][i]);

            pd->Gcycles_Av    = 0.0;
            pd->PME_f_load_Av = 0.0;
            pd->ns_per_day_Av = 0.0;

            if (pd->nPMEnodes == -1)
            {
                sprintf(strbuf, "(%3d)", pd->guessPME);
            }
            else
            {
                sprintf(strbuf, "     ");
            }

            /* Get the average run time of a setting */
            for (j = 0; j < nrepeats; j++)
            {
                pd->Gcycles_Av    += pd->Gcycles[j];
                pd->PME_f_load_Av += pd->PME_f_load[j];
            }
            pd->Gcycles_Av    /= nrepeats;
            pd->PME_f_load_Av /= nrepeats;

            for (j = 0; j < nrepeats; j++)
            {
                if (pd->ns_per_day[j] > 0.0)
                {
                    pd->ns_per_day_Av += pd->ns_per_day[j];
                }
                else
                {
                    /* Somehow the performance number was not aquired for this run,
                     * therefor set the average to some negative value: */
                    pd->ns_per_day_Av = -1.0f*nrepeats;
                    break;
                }
            }
            pd->ns_per_day_Av /= nrepeats;

            /* Nicer output: */
            if (pd->PME_f_load_Av > 0.0)
            {
                sprintf(str_PME_f_load, "%12.3f", pd->PME_f_load_Av);
            }
            else
            {
                sprintf(str_PME_f_load, "%s", "         -  ");
            }


            /* We assume we had a successful run if both averages are positive */
            if (pd->Gcycles_Av > 0.0 && pd->ns_per_day_Av > 0.0)
            {
                /* Output statistics if repeats were done */
                if (nrepeats > 1)
                {
                    /* Calculate the standard deviation */
                    s = 0.0;
                    for (j = 0; j < nrepeats; j++)
                    {
                        s += gmx::square( pd->Gcycles[j] - pd->Gcycles_Av );
                    }
                    s /= (nrepeats - 1);
                    s  = std::sqrt(s);

                    fprintf(fp, "%4d %3d %4d%s %12.3f %12.3f %12.3f %s",
                            line, k, pd->nPMEnodes, strbuf, pd->Gcycles_Av, s,
                            pd->ns_per_day_Av, str_PME_f_load);
                    if (nnodes > 1)
                    {
                        fprintf(fp, "  %3d %3d %3d", pd->nx, pd->ny, pd->nz);
                    }
                    fprintf(fp, "\n");
                }
                /* Store the index of the best run found so far in 'winner': */
                if ( (k_win == -1) || (pd->Gcycles_Av < perfdata[k_win][i_win].Gcycles_Av) )
                {
                    k_win    = k;
                    i_win    = i;
                    line_win = line;
                }
                line++;
            }
        }
    }

    if (k_win == -1)
    {
        gmx_fatal(FARGS, "None of the runs was successful! Check %s for problems.", fn);
    }

    sep_line(fp);

    winPME = perfdata[k_win][i_win].nPMEnodes;

    if (1 == ntests)
    {
        /* We stuck to a fixed number of PME-only nodes */
        sprintf(strbuf, "settings No. %d", k_win);
    }
    else
    {
        /* We have optimized the number of PME-only nodes */
        if (winPME == -1)
        {
            sprintf(strbuf, "%s", "the automatic number of PME ranks");
        }
        else
        {
            sprintf(strbuf, "%d PME ranks", winPME);
        }
    }
    fprintf(fp, "Best performance was achieved with %s", strbuf);
    if ((nrepeats > 1) && (ntests > 1))
    {
        fprintf(fp, " (see line %d)", line_win);
    }
    fprintf(fp, "\n");

    /* Only mention settings if they were modified: */
    bRefinedCoul = !gmx_within_tol(info->rcoulomb[k_win], info->rcoulomb[0], GMX_REAL_EPS);
    bRefinedVdW  = !gmx_within_tol(info->rvdw[k_win], info->rvdw[0], GMX_REAL_EPS);
    bRefinedGrid = !(info->nkx[k_win] == info->nkx[0] &&
                     info->nky[k_win] == info->nky[0] &&
                     info->nkz[k_win] == info->nkz[0]);

    if (bRefinedCoul || bRefinedVdW || bRefinedGrid)
    {
        fprintf(fp, "Optimized PME settings:\n");
        bCanUseOrigTPR = FALSE;
    }
    else
    {
        bCanUseOrigTPR = TRUE;
    }

    if (bRefinedCoul)
    {
        fprintf(fp, "   New Coulomb radius: %f nm (was %f nm)\n", info->rcoulomb[k_win], info->rcoulomb[0]);
    }

    if (bRefinedVdW)
    {
        fprintf(fp, "   New Van der Waals radius: %f nm (was %f nm)\n", info->rvdw[k_win], info->rvdw[0]);
    }

    if (bRefinedGrid)
    {
        fprintf(fp, "   New Fourier grid xyz: %d %d %d (was %d %d %d)\n", info->nkx[k_win], info->nky[k_win], info->nkz[k_win],
                info->nkx[0], info->nky[0], info->nkz[0]);
    }

    if (bCanUseOrigTPR && ntprs > 1)
    {
        fprintf(fp, "and original PME settings.\n");
    }

    fflush(fp);

    /* Return the index of the mdp file that showed the highest performance
     * and the optimal number of PME nodes */
    *index_tpr    = k_win;
    *npme_optimal = winPME;

    return bCanUseOrigTPR;
}


/* Get the commands we need to set up the runs from environment variables */
static void get_program_paths(gmx_bool bThreads, char *cmd_mpirun[], char *cmd_mdrun[])
{
    char      *cp;
    const char def_mpirun[]   = "mpirun";

    const char empty_mpirun[] = "";

    /* Get the commands we need to set up the runs from environment variables */
    if (!bThreads)
    {
        if ( (cp = getenv("MPIRUN")) != nullptr)
        {
            *cmd_mpirun = gmx_strdup(cp);
        }
        else
        {
            *cmd_mpirun = gmx_strdup(def_mpirun);
        }
    }
    else
    {
        *cmd_mpirun = gmx_strdup(empty_mpirun);
    }

    if (*cmd_mdrun == nullptr)
    {
        /* The use of MDRUN is deprecated, but made available in 5.1
           for backward compatibility. It may be removed in a future
           version. */
        if ( (cp = getenv("MDRUN" )) != nullptr)
        {
            *cmd_mdrun = gmx_strdup(cp);
        }
        else
        {
            gmx_fatal(FARGS, "The way to call mdrun must be set in the -mdrun command-line flag.");
        }
    }
}

/* Check that the commands will run mdrun (perhaps via mpirun) by
 * running a very quick test simulation. Requires MPI environment or
 * GPU support to be available if applicable. */
/* TODO implement feature to parse the log file to get the list of
   compatible GPUs from mdrun, if the user of gmx tune-pme has not
   given one. */
static void check_mdrun_works(gmx_bool    bThreads,
                              const char *cmd_mpirun,
                              const char *cmd_np,
                              const char *cmd_mdrun,
                              gmx_bool    bNeedGpuSupport)
{
    char      *command = nullptr;
    char      *cp;
    char       line[STRLEN];
    FILE      *fp;
    const char filename[]     = "benchtest.log";

    /* This string should always be identical to the one in copyrite.c,
     * gmx_print_version_info() in the GMX_MPI section */
    const char match_mpi[]     = "MPI library:        MPI";
    const char match_mdrun[]   = "Executable: ";
    const char match_nogpu[]   = "GPU support:        disabled";
    gmx_bool   bMdrun          = FALSE;
    gmx_bool   bMPI            = FALSE;
    gmx_bool   bHaveGpuSupport = TRUE;

    /* Run a small test to see whether mpirun + mdrun work  */
    fprintf(stdout, "Making sure that mdrun can be executed. ");
    if (bThreads)
    {
        snew(command, std::strlen(cmd_mdrun) + std::strlen(cmd_np) + std::strlen(filename) + 50);
        sprintf(command, "%s%s -version -maxh 0.001 1> %s 2>&1", cmd_mdrun, cmd_np, filename);
    }
    else
    {
        snew(command, std::strlen(cmd_mpirun) + std::strlen(cmd_np) + std::strlen(cmd_mdrun) + std::strlen(filename) + 50);
        sprintf(command, "%s%s%s -version -maxh 0.001 1> %s 2>&1", cmd_mpirun, cmd_np, cmd_mdrun, filename);
    }
    fprintf(stdout, "Trying '%s' ... ", command);
    make_backup(filename);
    gmx_system_call(command);

    /* Check if we find the characteristic string in the output: */
    if (!gmx_fexist(filename))
    {
        gmx_fatal(FARGS, "Output from test run could not be found.");
    }

    fp = fopen(filename, "r");
    /* We need to scan the whole output file, since sometimes the queuing system
     * also writes stuff to stdout/err */
    while (!feof(fp) )
    {
        cp = fgets(line, STRLEN, fp);
        if (cp != nullptr)
        {
            if (str_starts(line, match_mdrun) )
            {
                bMdrun = TRUE;
            }
            if (str_starts(line, match_mpi) )
            {
                bMPI = TRUE;
            }
            if (str_starts(line, match_nogpu) )
            {
                bHaveGpuSupport = FALSE;
            }
        }
    }
    fclose(fp);

    if (bThreads)
    {
        if (bMPI)
        {
            gmx_fatal(FARGS, "Need a threaded version of mdrun. This one\n"
                      "(%s)\n"
                      "seems to have been compiled with MPI instead.",
                      cmd_mdrun);
        }
    }
    else
    {
        if (bMdrun && !bMPI)
        {
            gmx_fatal(FARGS, "Need an MPI-enabled version of mdrun. This one\n"
                      "(%s)\n"
                      "seems to have been compiled without MPI support.",
                      cmd_mdrun);
        }
    }

    if (!bMdrun)
    {
        gmx_fatal(FARGS, "Cannot execute mdrun. Please check %s for problems!",
                  filename);
    }

    if (bNeedGpuSupport && !bHaveGpuSupport)
    {
        gmx_fatal(FARGS, "The mdrun executable did not have the expected GPU support.");
    }

    fprintf(stdout, "passed.\n");

    /* Clean up ... */
    remove(filename);
    sfree(command);
}

/* Handles the no-GPU case by emitting an empty string. */
static std::string make_gpu_id_command_line(const char *eligible_gpu_ids)
{
    /* If the user has given no eligible GPU IDs, or we're trying the
     * default behaviour, then there is nothing for tune_pme to give
     * to mdrun -gpu_id */
    if (eligible_gpu_ids != nullptr)
    {
        return gmx::formatString("-gpu_id %s", eligible_gpu_ids);
    }


    return std::string();
}

static void launch_simulation(
        gmx_bool    bLaunch,          /* Should the simulation be launched? */
        FILE       *fp,               /* General log file */
        gmx_bool    bThreads,         /* whether to use threads */
        char       *cmd_mpirun,       /* Command for mpirun */
        char       *cmd_np,           /* Switch for -np or -ntmpi or empty */
        char       *cmd_mdrun,        /* Command for mdrun */
        char       *args_for_mdrun,   /* Arguments for mdrun */
        const char *simulation_tpr,   /* This tpr will be simulated */
        int         nPMEnodes,        /* Number of PME ranks to use */
        const char *eligible_gpu_ids) /* Available GPU IDs for
                                       * constructing mdrun command lines */
{
    char  *command;


    /* Make enough space for the system call command,
     * (200 extra chars for -npme ... etc. options should suffice): */
    snew(command, std::strlen(cmd_mpirun)+std::strlen(cmd_mdrun)+std::strlen(cmd_np)+std::strlen(args_for_mdrun)+std::strlen(simulation_tpr)+200);

    auto cmd_gpu_ids = make_gpu_id_command_line(eligible_gpu_ids);

    /* Note that the -passall options requires args_for_mdrun to be at the end
     * of the command line string */
    if (bThreads)
    {
        sprintf(command, "%s%s-npme %d -s %s %s %s",
                cmd_mdrun, cmd_np, nPMEnodes, simulation_tpr, args_for_mdrun, cmd_gpu_ids.c_str());
    }
    else
    {
        sprintf(command, "%s%s%s -npme %d -s %s %s %s",
                cmd_mpirun, cmd_np, cmd_mdrun, nPMEnodes, simulation_tpr, args_for_mdrun, cmd_gpu_ids.c_str());
    }

    fprintf(fp, "%s this command line to launch the simulation:\n\n%s", bLaunch ? "Using" : "Please use", command);
    sep_line(fp);
    fflush(fp);

    /* Now the real thing! */
    if (bLaunch)
    {
        fprintf(stdout, "\nLaunching simulation with best parameters now.\nExecuting '%s'", command);
        sep_line(stdout);
        fflush(stdout);
        gmx_system_call(command);
    }
}


static void modify_PMEsettings(
        gmx_int64_t     simsteps,    /* Set this value as number of time steps */
        gmx_int64_t     init_step,   /* Set this value as init_step */
        const char     *fn_best_tpr, /* tpr file with the best performance */
        const char     *fn_sim_tpr)  /* name of tpr file to be launched */
{
    t_state        state;
    gmx_mtop_t     mtop;
    char           buf[200];

    t_inputrec     irInstance;
    t_inputrec    *ir = &irInstance;
    read_tpx_state(fn_best_tpr, ir, &state, &mtop);

    /* Reset nsteps and init_step to the value of the input .tpr file */
    ir->nsteps    = simsteps;
    ir->init_step = init_step;

    /* Write the tpr file which will be launched */
    sprintf(buf, "Writing optimized simulation file %s with nsteps=%s.\n", fn_sim_tpr, "%" GMX_PRId64);
    fprintf(stdout, buf, ir->nsteps);
    fflush(stdout);
    write_tpx_state(fn_sim_tpr, ir, &state, &mtop);
}

static gmx_bool can_scale_rvdw(int vdwtype)
{
    return (evdwCUT == vdwtype ||
            evdwPME == vdwtype);
}

#define EPME_SWITCHED(e) ((e) == eelPMESWITCH || (e) == eelPMEUSERSWITCH)

/* Make additional TPR files with more computational load for the
 * direct space processors: */
static void make_benchmark_tprs(
        const char     *fn_sim_tpr,      /* READ : User-provided tpr file                 */
        char           *fn_bench_tprs[], /* WRITE: Names of benchmark tpr files           */
        gmx_int64_t     benchsteps,      /* Number of time steps for benchmark runs       */
        gmx_int64_t     statesteps,      /* Step counter in checkpoint file               */
        real            rmin,            /* Minimal Coulomb radius                        */
        real            rmax,            /* Maximal Coulomb radius                        */
        real            bScaleRvdw,      /* Scale rvdw along with rcoulomb                */
        int            *ntprs,           /* No. of TPRs to write, each with a different
                                            rcoulomb and fourierspacing                   */
        t_inputinfo    *info,            /* Contains information about mdp file options   */
        FILE           *fp)              /* Write the output here                         */
{
    int           i, j, d;
    t_state       state;
    gmx_mtop_t    mtop;
    real          nlist_buffer;     /* Thickness of the buffer regions for PME-switch potentials */
    char          buf[200];
    rvec          box_size;
    gmx_bool      bNote = FALSE;
    real          add;              /* Add this to rcoul for the next test    */
    real          fac = 1.0;        /* Scaling factor for Coulomb radius      */
    real          fourierspacing;   /* Basic fourierspacing from tpr          */


    sprintf(buf, "Making benchmark tpr file%s with %s time step%s",
            *ntprs > 1 ? "s" : "", "%" GMX_PRId64, benchsteps > 1 ? "s" : "");
    fprintf(stdout, buf, benchsteps);
    if (statesteps > 0)
    {
        sprintf(buf, " (adding %s steps from checkpoint file)", "%" GMX_PRId64);
        fprintf(stdout, buf, statesteps);
        benchsteps += statesteps;
    }
    fprintf(stdout, ".\n");

    t_inputrec  irInstance;
    t_inputrec *ir = &irInstance;
    read_tpx_state(fn_sim_tpr, ir, &state, &mtop);

    /* Check if some kind of PME was chosen */
    if (EEL_PME(ir->coulombtype) == FALSE)
    {
        gmx_fatal(FARGS, "Can only do optimizations for simulations with %s electrostatics.",
                  EELTYPE(eelPME));
    }

    /* Check if rcoulomb == rlist, which is necessary for plain PME. */
    if (  (ir->cutoff_scheme != ecutsVERLET) &&
          (eelPME == ir->coulombtype) && !(ir->rcoulomb == ir->rlist))
    {
        gmx_fatal(FARGS, "%s requires rcoulomb (%f) to be equal to rlist (%f).",
                  EELTYPE(eelPME), ir->rcoulomb, ir->rlist);
    }
    /* For other PME types, rcoulomb is allowed to be smaller than rlist */
    else if (ir->rcoulomb > ir->rlist)
    {
        gmx_fatal(FARGS, "%s requires rcoulomb (%f) to be equal to or smaller than rlist (%f)",
                  EELTYPE(ir->coulombtype), ir->rcoulomb, ir->rlist);
    }

    if (bScaleRvdw && ir->rvdw != ir->rcoulomb)
    {
        fprintf(stdout, "NOTE: input rvdw != rcoulomb, will not scale rvdw\n");
        bScaleRvdw = FALSE;
    }

    /* Reduce the number of steps for the benchmarks */
    info->orig_sim_steps = ir->nsteps;
    ir->nsteps           = benchsteps;
    /* We must not use init_step from the input tpr file for the benchmarks */
    info->orig_init_step = ir->init_step;
    ir->init_step        = 0;

    /* For PME-switch potentials, keep the radial distance of the buffer region */
    nlist_buffer   = ir->rlist - ir->rcoulomb;

    /* Determine length of triclinic box vectors */
    for (d = 0; d < DIM; d++)
    {
        box_size[d] = 0;
        for (i = 0; i < DIM; i++)
        {
            box_size[d] += state.box[d][i]*state.box[d][i];
        }
        box_size[d] = std::sqrt(box_size[d]);
    }

    if (ir->fourier_spacing > 0)
    {
        info->fsx[0] = ir->fourier_spacing;
        info->fsy[0] = ir->fourier_spacing;
        info->fsz[0] = ir->fourier_spacing;
    }
    else
    {
        /* Reconstruct fourierspacing per dimension from the number of grid points and box size */
        info->fsx[0] = box_size[XX]/ir->nkx;
        info->fsy[0] = box_size[YY]/ir->nky;
        info->fsz[0] = box_size[ZZ]/ir->nkz;
    }

    /* If no value for the fourierspacing was provided on the command line, we
     * use the reconstruction from the tpr file */
    if (ir->fourier_spacing > 0)
    {
        /* Use the spacing from the tpr */
        fourierspacing = ir->fourier_spacing;
    }
    else
    {
        /* Use the maximum observed spacing */
        fourierspacing = std::max(std::max(info->fsx[0], info->fsy[0]), info->fsz[0]);
    }

    fprintf(stdout, "Calculating PME grid points on the basis of a fourierspacing of %f nm\n", fourierspacing);

    /* For performance comparisons the number of particles is useful to have */
    fprintf(fp, "   Number of particles  : %d\n", mtop.natoms);

    /* Print information about settings of which some are potentially modified: */
    fprintf(fp, "   Coulomb type         : %s\n", EELTYPE(ir->coulombtype));
    fprintf(fp, "   Grid spacing x y z   : %f %f %f\n",
            box_size[XX]/ir->nkx, box_size[YY]/ir->nky, box_size[ZZ]/ir->nkz);
    fprintf(fp, "   Van der Waals type   : %s\n", EVDWTYPE(ir->vdwtype));
    if (ir_vdw_switched(ir))
    {
        fprintf(fp, "   rvdw_switch          : %f nm\n", ir->rvdw_switch);
    }
    if (EPME_SWITCHED(ir->coulombtype))
    {
        fprintf(fp, "   rlist                : %f nm\n", ir->rlist);
    }

    /* Print a descriptive line about the tpr settings tested */
    fprintf(fp, "\nWill try these real/reciprocal workload settings:\n");
    fprintf(fp, " No.   scaling  rcoulomb");
    fprintf(fp, "  nkx  nky  nkz");
    fprintf(fp, "   spacing");
    if (can_scale_rvdw(ir->vdwtype))
    {
        fprintf(fp, "      rvdw");
    }
    if (EPME_SWITCHED(ir->coulombtype))
    {
        fprintf(fp, "     rlist");
    }
    fprintf(fp, "  tpr file\n");

    /* Loop to create the requested number of tpr input files */
    for (j = 0; j < *ntprs; j++)
    {
        /* The first .tpr is the provided one, just need to modify nsteps,
         * so skip the following block */
        if (j != 0)
        {
            /* Determine which Coulomb radii rc to use in the benchmarks */
            add = (rmax-rmin)/(*ntprs-1);
            if (gmx_within_tol(rmin, info->rcoulomb[0], GMX_REAL_EPS))
            {
                ir->rcoulomb = rmin + j*add;
            }
            else if (gmx_within_tol(rmax, info->rcoulomb[0], GMX_REAL_EPS))
            {
                ir->rcoulomb = rmin + (j-1)*add;
            }
            else
            {
                /* rmin != rcoul != rmax, ergo test between rmin and rmax */
                add          = (rmax-rmin)/(*ntprs-2);
                ir->rcoulomb = rmin + (j-1)*add;
            }

            /* Determine the scaling factor fac */
            fac = ir->rcoulomb/info->rcoulomb[0];

            /* Scale the Fourier grid spacing */
            ir->nkx = ir->nky = ir->nkz = 0;
            calcFftGrid(nullptr, state.box, fourierspacing*fac, minimalPmeGridSize(ir->pme_order),
                        &ir->nkx, &ir->nky, &ir->nkz);

            /* Adjust other radii since various conditions need to be fulfilled */
            if (eelPME == ir->coulombtype)
            {
                /* plain PME, rcoulomb must be equal to rlist TODO only in the group scheme? */
                ir->rlist = ir->rcoulomb;
            }
            else
            {
                /* rlist must be >= rcoulomb, we keep the size of the buffer region */
                ir->rlist = ir->rcoulomb + nlist_buffer;
            }

            if (bScaleRvdw && can_scale_rvdw(ir->vdwtype))
            {
                if (ecutsVERLET == ir->cutoff_scheme ||
                    evdwPME == ir->vdwtype)
                {
                    /* With either the Verlet cutoff-scheme or LJ-PME,
                       the van der Waals radius must always equal the
                       Coulomb radius */
                    ir->rvdw = ir->rcoulomb;
                }
                else
                {
                    /* For vdw cutoff, rvdw >= rlist */
                    ir->rvdw = std::max(info->rvdw[0], ir->rlist);
                }
            }
        } /* end of "if (j != 0)" */

        /* for j==0: Save the original settings
         * for j >0: Save modified radii and Fourier grids */
        info->rcoulomb[j]  = ir->rcoulomb;
        info->rvdw[j]      = ir->rvdw;
        info->nkx[j]       = ir->nkx;
        info->nky[j]       = ir->nky;
        info->nkz[j]       = ir->nkz;
        info->rlist[j]     = ir->rlist;
        info->fsx[j]       = fac*fourierspacing;
        info->fsy[j]       = fac*fourierspacing;
        info->fsz[j]       = fac*fourierspacing;

        /* Write the benchmark tpr file */
        std::strncpy(fn_bench_tprs[j], fn_sim_tpr, std::strlen(fn_sim_tpr)-std::strlen(".tpr"));
        sprintf(buf, "_bench%.2d.tpr", j);
        std::strcat(fn_bench_tprs[j], buf);
        fprintf(stdout, "Writing benchmark tpr %s with nsteps=", fn_bench_tprs[j]);
        fprintf(stdout, "%" GMX_PRId64, ir->nsteps);
        if (j > 0)
        {
            fprintf(stdout, ", scaling factor %f\n", fac);
        }
        else
        {
            fprintf(stdout, ", unmodified settings\n");
        }

        write_tpx_state(fn_bench_tprs[j], ir, &state, &mtop);

        /* Write information about modified tpr settings to log file */
        fprintf(fp, "%4d%10f%10f", j, fac, ir->rcoulomb);
        fprintf(fp, "%5d%5d%5d", ir->nkx, ir->nky, ir->nkz);
        fprintf(fp, " %9f ", info->fsx[j]);
        if (can_scale_rvdw(ir->vdwtype))
        {
            fprintf(fp, "%10f", ir->rvdw);
        }
        if (EPME_SWITCHED(ir->coulombtype))
        {
            fprintf(fp, "%10f", ir->rlist);
        }
        fprintf(fp, "  %-14s\n", fn_bench_tprs[j]);

        /* Make it clear to the user that some additional settings were modified */
        if (!gmx_within_tol(ir->rvdw, info->rvdw[0], GMX_REAL_EPS)
            || !gmx_within_tol(ir->rlist, info->rlist[0], GMX_REAL_EPS) )
        {
            bNote = TRUE;
        }
    }
    if (bNote)
    {
        fprintf(fp, "\nNote that in addition to the Coulomb radius and the Fourier grid\n"
                "other input settings were also changed (see table above).\n"
                "Please check if the modified settings are appropriate.\n");
    }
    fflush(stdout);
    fflush(fp);
}


/* Rename the files we want to keep to some meaningful filename and
 * delete the rest */
static void cleanup(const t_filenm *fnm, int nfile, int k, int nnodes,
                    int nPMEnodes, int nr, gmx_bool bKeepStderr)
{
    char        numstring[STRLEN];
    char        newfilename[STRLEN];
    const char *fn = nullptr;
    int         i;
    const char *opt;


    fprintf(stdout, "Cleaning up, deleting benchmark temp files ...\n");

    for (i = 0; i < nfile; i++)
    {
        opt = (char *)fnm[i].opt;
        if (std::strcmp(opt, "-p") == 0)
        {
            /* do nothing; keep this file */
            ;
        }
        else if (std::strcmp(opt, "-bg") == 0)
        {
            /* Give the log file a nice name so one can later see which parameters were used */
            numstring[0] = '\0';
            if (nr > 0)
            {
                sprintf(numstring, "_%d", nr);
            }
            sprintf(newfilename, "%s_no%d_np%d_npme%d%s", opt2fn("-bg", nfile, fnm), k, nnodes, nPMEnodes, numstring);
            if (gmx_fexist(opt2fn("-bg", nfile, fnm)))
            {
                fprintf(stdout, "renaming log file to %s\n", newfilename);
                make_backup(newfilename);
                rename(opt2fn("-bg", nfile, fnm), newfilename);
            }
        }
        else if (std::strcmp(opt, "-err") == 0)
        {
            /* This file contains the output of stderr. We want to keep it in
             * cases where there have been problems. */
            fn           = opt2fn(opt, nfile, fnm);
            numstring[0] = '\0';
            if (nr > 0)
            {
                sprintf(numstring, "_%d", nr);
            }
            sprintf(newfilename, "%s_no%d_np%d_npme%d%s", fn, k, nnodes, nPMEnodes, numstring);
            if (gmx_fexist(fn))
            {
                if (bKeepStderr)
                {
                    fprintf(stdout, "Saving stderr output in %s\n", newfilename);
                    make_backup(newfilename);
                    rename(fn, newfilename);
                }
                else
                {
                    fprintf(stdout, "Deleting %s\n", fn);
                    remove(fn);
                }
            }
        }
        /* Delete the files which are created for each benchmark run: (options -b*) */
        else if ( (0 == std::strncmp(opt, "-b", 2)) && (opt2bSet(opt, nfile, fnm) || !is_optional(&fnm[i])) )
        {
            remove_if_exists(opt2fn(opt, nfile, fnm));
        }
    }
}


enum {
    eNpmeAuto, eNpmeAll, eNpmeReduced, eNpmeSubset, eNpmeNr
};

/* Create a list of numbers of PME nodes to test */
static void make_npme_list(
        const char *npmevalues_opt, /* Make a complete list with all
                                     * possibilities or a short list that keeps only
                                     * reasonable numbers of PME nodes                  */
        int        *nentries,       /* Number of entries we put in the nPMEnodes list   */
        int        *nPMEnodes[],    /* Each entry contains the value for -npme          */
        int         nnodes,         /* Total number of nodes to do the tests on         */
        int         minPMEnodes,    /* Minimum number of PME nodes                      */
        int         maxPMEnodes)    /* Maximum number of PME nodes                      */
{
    int i, npme, npp;
    int min_factor = 1;   /* We request that npp and npme have this minimal
                           * largest common factor (depends on npp)           */
    int nlistmax;         /* Max. list size                                   */
    int nlist;            /* Actual number of entries in list                 */
    int eNPME = 0;


    /* Do we need to check all possible values for -npme or is a reduced list enough? */
    if (!std::strcmp(npmevalues_opt, "all") )
    {
        eNPME = eNpmeAll;
    }
    else if (!std::strcmp(npmevalues_opt, "subset") )
    {
        eNPME = eNpmeSubset;
    }
    else /* "auto" or "range" */
    {
        if (nnodes <= 64)
        {
            eNPME = eNpmeAll;
        }
        else if (nnodes < 128)
        {
            eNPME = eNpmeReduced;
        }
        else
        {
            eNPME = eNpmeSubset;
        }
    }

    /* Calculate how many entries we could possibly have (in case of -npme all) */
    if (nnodes > 2)
    {
        nlistmax = maxPMEnodes - minPMEnodes + 3;
        if (0 == minPMEnodes)
        {
            nlistmax--;
        }
    }
    else
    {
        nlistmax = 1;
    }

    /* Now make the actual list which is at most of size nlist */
    snew(*nPMEnodes, nlistmax);
    nlist = 0; /* start counting again, now the real entries in the list */
    for (i = 0; i < nlistmax - 2; i++)
    {
        npme = maxPMEnodes - i;
        npp  = nnodes-npme;
        switch (eNPME)
        {
            case eNpmeAll:
                min_factor = 1;
                break;
            case eNpmeReduced:
                min_factor = 2;
                break;
            case eNpmeSubset:
                /* For 2d PME we want a common largest factor of at least the cube
                 * root of the number of PP nodes */
                min_factor = static_cast<int>(std::cbrt(npp));
                break;
            default:
                gmx_fatal(FARGS, "Unknown option for eNPME in make_npme_list");
                break;
        }
        if (gmx_greatest_common_divisor(npp, npme) >= min_factor)
        {
            (*nPMEnodes)[nlist] = npme;
            nlist++;
        }
    }
    /* We always test 0 PME nodes and the automatic number */
    *nentries             = nlist + 2;
    (*nPMEnodes)[nlist  ] =  0;
    (*nPMEnodes)[nlist+1] = -1;

    fprintf(stderr, "Will try the following %d different values for -npme:\n", *nentries);
    for (i = 0; i < *nentries-1; i++)
    {
        fprintf(stderr, "%d, ", (*nPMEnodes)[i]);
    }
    fprintf(stderr, "and %d (auto).\n", (*nPMEnodes)[*nentries-1]);
}


/* Allocate memory to store the performance data */
static void init_perfdata(t_perf *perfdata[], int ntprs, int datasets, int repeats)
{
    int i, j, k;


    for (k = 0; k < ntprs; k++)
    {
        snew(perfdata[k], datasets);
        for (i = 0; i < datasets; i++)
        {
            for (j = 0; j < repeats; j++)
            {
                snew(perfdata[k][i].Gcycles, repeats);
                snew(perfdata[k][i].ns_per_day, repeats);
                snew(perfdata[k][i].PME_f_load, repeats);
            }
        }
    }
}


/* Check for errors on mdrun -h */
static void make_sure_it_runs(char *mdrun_cmd_line, int length, FILE *fp,
                              const t_filenm *fnm, int nfile)
{
    char       *command, *msg;
    int         ret;

    snew(command, length +  15);
    snew(msg, length + 500);

    fprintf(stdout, "Making sure the benchmarks can be executed by running just 1 step...\n");
    sprintf(command, "%s -nsteps 1 -quiet", mdrun_cmd_line);
    fprintf(stdout, "Executing '%s' ...\n", command);
    ret = gmx_system_call(command);

    if (0 != ret)
    {
        /* To prevent confusion, do not again issue a gmx_fatal here since we already
         * get the error message from mdrun itself */
        sprintf(msg,
                "Cannot run the first benchmark simulation! Please check the error message of\n"
                "mdrun for the source of the problem. Did you provide a command line\n"
                "argument that neither gmx tune_pme nor mdrun understands? If you're\n"
                "sure your command line should work, you can bypass this check with \n"
                "gmx tune_pme -nocheck. The failing command was:\n"
                "\n%s\n\n", command);

        fprintf(stderr, "%s", msg);
        sep_line(fp);
        fprintf(fp, "%s", msg);

        exit(ret);
    }
    fprintf(stdout, "Benchmarks can be executed!\n");

    /* Clean up the benchmark output files we just created */
    fprintf(stdout, "Cleaning up ...\n");
    remove_if_exists(opt2fn("-bc", nfile, fnm));
    remove_if_exists(opt2fn("-be", nfile, fnm));
    remove_if_exists(opt2fn("-bcpo", nfile, fnm));
    remove_if_exists(opt2fn("-bg", nfile, fnm));
    remove_if_exists(opt2fn("-bo", nfile, fnm));
    remove_if_exists(opt2fn("-bx", nfile, fnm));

    sfree(command);
    sfree(msg    );
}

static void do_the_tests(
        FILE           *fp,               /* General tune_pme output file           */
        char          **tpr_names,        /* Filenames of the input files to test   */
        int             maxPMEnodes,      /* Max fraction of nodes to use for PME   */
        int             minPMEnodes,      /* Min fraction of nodes to use for PME   */
        int             npme_fixed,       /* If >= -1, test fixed number of PME
                                           * nodes only                             */
        const char     *npmevalues_opt,   /* Which -npme values should be tested    */
        t_perf        **perfdata,         /* Here the performace data is stored     */
        int            *pmeentries,       /* Entries in the nPMEnodes list          */
        int             repeats,          /* Repeat each test this often            */
        int             nnodes,           /* Total number of nodes = nPP + nPME     */
        int             nr_tprs,          /* Total number of tpr files to test      */
        gmx_bool        bThreads,         /* Threads or MPI?                        */
        char           *cmd_mpirun,       /* mpirun command string                  */
        char           *cmd_np,           /* "-np", "-n", whatever mpirun needs     */
        char           *cmd_mdrun,        /* mdrun command string                   */
        char           *cmd_args_bench,   /* arguments for mdrun in a string        */
        const t_filenm *fnm,              /* List of filenames from command line    */
        int             nfile,            /* Number of files specified on the cmdl. */
        int             presteps,         /* DLB equilibration steps, is checked    */
        gmx_int64_t     cpt_steps,        /* Time step counter in the checkpoint    */
        gmx_bool        bCheck,           /* Check whether benchmark mdrun works    */
        const char     *eligible_gpu_ids) /* GPU IDs for
                                           * constructing mdrun command lines */
{
    int      i, nr, k, ret, count = 0, totaltests;
    int     *nPMEnodes = nullptr;
    t_perf  *pd        = nullptr;
    int      cmdline_length;
    char    *command, *cmd_stub;
    char     buf[STRLEN];
    gmx_bool bResetProblem = FALSE;
    gmx_bool bFirst        = TRUE;

    /* This string array corresponds to the eParselog enum type at the start
     * of this file */
    const char* ParseLog[] = {
        "OK.",
        "Logfile not found!",
        "No timings, logfile truncated?",
        "Run was terminated.",
        "Counters were not reset properly.",
        "No DD grid found for these settings.",
        "TPX version conflict!",
        "mdrun was not started in parallel!",
        "Number of PP ranks has a prime factor that is too large.",
        "The number of PP ranks did not suit the number of GPUs.",
        "Some GPUs were not detected or are incompatible.",
        "An error occurred."
    };
    char        str_PME_f_load[13];


    /* Allocate space for the mdrun command line. 100 extra characters should
       be more than enough for the -npme etcetera arguments */
    cmdline_length =  std::strlen(cmd_mpirun)
        + std::strlen(cmd_np)
        + std::strlen(cmd_mdrun)
        + std::strlen(cmd_args_bench)
        + std::strlen(tpr_names[0]) + 100;
    snew(command, cmdline_length);
    snew(cmd_stub, cmdline_length);

    /* Construct the part of the command line that stays the same for all tests: */
    if (bThreads)
    {
        sprintf(cmd_stub, "%s%s", cmd_mdrun, cmd_np);
    }
    else
    {
        sprintf(cmd_stub, "%s%s%s ", cmd_mpirun, cmd_np, cmd_mdrun);
    }

    /* Create a list of numbers of PME nodes to test */
    if (npme_fixed < -1)
    {
        make_npme_list(npmevalues_opt, pmeentries, &nPMEnodes,
                       nnodes, minPMEnodes, maxPMEnodes);
    }
    else
    {
        *pmeentries  = 1;
        snew(nPMEnodes, 1);
        nPMEnodes[0] = npme_fixed;
        fprintf(stderr, "Will use a fixed number of %d PME-only ranks.\n", nPMEnodes[0]);
    }

    if (0 == repeats)
    {
        fprintf(fp, "\nNo benchmarks done since number of repeats (-r) is 0.\n");
        gmx_ffclose(fp);
        finalize(opt2fn("-p", nfile, fnm));
        exit(0);
    }

    /* Allocate one dataset for each tpr input file: */
    init_perfdata(perfdata, nr_tprs, *pmeentries, repeats);

    /*****************************************/
    /* Main loop over all tpr files to test: */
    /*****************************************/
    totaltests = nr_tprs*(*pmeentries)*repeats;
    for (k = 0; k < nr_tprs; k++)
    {
        fprintf(fp, "\nIndividual timings for input file %d (%s):\n", k, tpr_names[k]);
        fprintf(fp, "PME ranks      Gcycles       ns/day        PME/f    Remark\n");
        /* Loop over various numbers of PME nodes: */
        for (i = 0; i < *pmeentries; i++)
        {
            pd = &perfdata[k][i];

            auto cmd_gpu_ids = make_gpu_id_command_line(eligible_gpu_ids);

            /* Loop over the repeats for each scenario: */
            for (nr = 0; nr < repeats; nr++)
            {
                pd->nPMEnodes = nPMEnodes[i];

                /* Add -npme and -s to the command line and save it. Note that
                 * the -passall (if set) options requires cmd_args_bench to be
                 * at the end of the command line string */
                snew(pd->mdrun_cmd_line, cmdline_length);
                sprintf(pd->mdrun_cmd_line, "%s-npme %d -s %s %s %s",
                        cmd_stub, pd->nPMEnodes, tpr_names[k], cmd_args_bench, cmd_gpu_ids.c_str());

                /* To prevent that all benchmarks fail due to a show-stopper argument
                 * on the mdrun command line, we make a quick check first.
                 * This check can be turned off in cases where the automatically chosen
                 * number of PME-only ranks leads to a number of PP ranks for which no
                 * decomposition can be found (e.g. for large prime numbers) */
                if (bFirst && bCheck)
                {
                    /* TODO This check is really for a functional
                     * .tpr, and if we need it, it should take place
                     * for every .tpr, and the logic for it should be
                     * immediately inside the loop over k, not in
                     * this inner loop. */
                    char *temporary_cmd_line;

                    snew(temporary_cmd_line, cmdline_length);
                    /* TODO -npme 0 is more likely to succeed at low
                       parallelism than the default of -npme -1, but
                       is more likely to fail at the scaling limit
                       when the PP domains may be too small. "mpirun
                       -np 1 mdrun" is probably a reasonable thing to
                       do for this check, but it'll be easier to
                       implement that after some refactoring of how
                       the number of MPI ranks is managed. */
                    sprintf(temporary_cmd_line, "%s-npme 0 -nb cpu -s %s %s",
                            cmd_stub, tpr_names[k], cmd_args_bench);
                    make_sure_it_runs(temporary_cmd_line, cmdline_length, fp, fnm, nfile);
                }
                bFirst = FALSE;

                /* Do a benchmark simulation: */
                if (repeats > 1)
                {
                    sprintf(buf, ", pass %d/%d", nr+1, repeats);
                }
                else
                {
                    buf[0] = '\0';
                }
                fprintf(stdout, "\n=== Progress %2.0f%%, tpr %d/%d, run %d/%d%s:\n",
                        (100.0*count)/totaltests,
                        k+1, nr_tprs, i+1, *pmeentries, buf);
                make_backup(opt2fn("-err", nfile, fnm));
                sprintf(command, "%s 1> /dev/null 2>%s", pd->mdrun_cmd_line, opt2fn("-err", nfile, fnm));
                fprintf(stdout, "%s\n", pd->mdrun_cmd_line);
                gmx_system_call(command);

                /* Collect the performance data from the log file; also check stderr
                 * for fatal errors */
                ret = parse_logfile(opt2fn("-bg", nfile, fnm), opt2fn("-err", nfile, fnm),
                                    pd, nr, presteps, cpt_steps, nnodes);
                if ((presteps > 0) && (ret == eParselogResetProblem))
                {
                    bResetProblem = TRUE;
                }

                if (-1 == pd->nPMEnodes)
                {
                    sprintf(buf, "(%3d)", pd->guessPME);
                }
                else
                {
                    sprintf(buf, "     ");
                }

                /* Nicer output */
                if (pd->PME_f_load[nr] > 0.0)
                {
                    sprintf(str_PME_f_load, "%12.3f", pd->PME_f_load[nr]);
                }
                else
                {
                    sprintf(str_PME_f_load, "%s", "         -  ");
                }

                /* Write the data we got to disk */
                fprintf(fp, "%4d%s %12.3f %12.3f %s    %s", pd->nPMEnodes,
                        buf, pd->Gcycles[nr], pd->ns_per_day[nr], str_PME_f_load, ParseLog[ret]);
                if (!(ret == eParselogOK || ret == eParselogNoDDGrid || ret == eParselogNotFound) )
                {
                    fprintf(fp, " Check %s file for problems.", ret == eParselogFatal ? "err" : "log");
                }
                fprintf(fp, "\n");
                fflush(fp);
                count++;

                /* Do some cleaning up and delete the files we do not need any more */
                cleanup(fnm, nfile, k, nnodes, pd->nPMEnodes, nr, ret == eParselogFatal);

                /* If the first run with this number of processors already failed, do not try again: */
                if (pd->Gcycles[0] <= 0.0 && repeats > 1)
                {
                    fprintf(stdout, "Skipping remaining passes of unsuccessful setting, see log file for details.\n");
                    count += repeats-(nr+1);
                    break;
                }
            } /* end of repeats loop */
        }     /* end of -npme loop */
    }         /* end of tpr file loop */

    if (bResetProblem)
    {
        sep_line(fp);
        fprintf(fp, "WARNING: The cycle and time step counters could not be reset properly. ");
        sep_line(fp);
    }
    sfree(command);
    sfree(cmd_stub);
}


static void check_input(
        int             nnodes,
        int             repeats,
        int            *ntprs,
        real           *rmin,
        real            rcoulomb,
        real           *rmax,
        real            maxPMEfraction,
        real            minPMEfraction,
        int             npme_fixed,
        gmx_int64_t     bench_nsteps,
        const t_filenm *fnm,
        int             nfile,
        int             sim_part,
        int             presteps,
        int             npargs,
        t_pargs        *pa)
{
    int old;


    /* Make sure the input file exists */
    if (!gmx_fexist(opt2fn("-s", nfile, fnm)))
    {
        gmx_fatal(FARGS, "File %s not found.", opt2fn("-s", nfile, fnm));
    }

    /* Make sure that the checkpoint file is not overwritten during benchmarking */
    if ( (0 == std::strcmp(opt2fn("-cpi", nfile, fnm), opt2fn("-bcpo", nfile, fnm)) ) && (sim_part > 1) )
    {
        gmx_fatal(FARGS, "Checkpoint input (-cpi) and benchmark checkpoint output (-bcpo) files must not be identical.\n"
                  "The checkpoint input file must not be overwritten during the benchmarks.\n");
    }

    /* Make sure that repeats is >= 0 (if == 0, only write tpr files) */
    if (repeats < 0)
    {
        gmx_fatal(FARGS, "Number of repeats < 0!");
    }

    /* Check number of nodes */
    if (nnodes < 1)
    {
        gmx_fatal(FARGS, "Number of ranks/threads must be a positive integer.");
    }

    /* Automatically choose -ntpr if not set */
    if (*ntprs < 1)
    {
        if (nnodes < 16)
        {
            *ntprs = 1;
        }
        else
        {
            *ntprs = 3;
            /* Set a reasonable scaling factor for rcoulomb */
            if (*rmax <= 0)
            {
                *rmax = rcoulomb * 1.2;
            }
        }
        fprintf(stderr, "Will test %d tpr file%s.\n", *ntprs, *ntprs == 1 ? "" : "s");
    }
    else
    {
        if (1 == *ntprs)
        {
            fprintf(stderr, "Note: Choose ntpr>1 to shift PME load between real and reciprocal space.\n");
        }
    }

    /* Make shure that rmin <= rcoulomb <= rmax */
    if (*rmin <= 0)
    {
        *rmin = rcoulomb;
    }
    if (*rmax <= 0)
    {
        *rmax = rcoulomb;
    }
    if (!(*rmin <= *rmax) )
    {
        gmx_fatal(FARGS, "Please choose the Coulomb radii such that rmin <= rmax.\n"
                  "rmin = %g, rmax = %g, actual rcoul from .tpr file = %g\n", *rmin, *rmax, rcoulomb);
    }
    /* Add test scenarios if rmin or rmax were set */
    if (*ntprs <= 2)
    {
        if (!gmx_within_tol(*rmin, rcoulomb, GMX_REAL_EPS) && (*ntprs == 1) )
        {
            (*ntprs)++;
            fprintf(stderr, "NOTE: Setting -rmin to %g changed -ntpr to %d\n",
                    *rmin, *ntprs);
        }
        if (!gmx_within_tol(*rmax, rcoulomb, GMX_REAL_EPS) && (*ntprs == 1) )
        {
            (*ntprs)++;
            fprintf(stderr, "NOTE: Setting -rmax to %g changed -ntpr to %d\n",
                    *rmax, *ntprs);
        }
    }
    old = *ntprs;
    /* If one of rmin, rmax is set, we need 2 tpr files at minimum */
    if (!gmx_within_tol(*rmax, rcoulomb, GMX_REAL_EPS) || !gmx_within_tol(*rmin, rcoulomb, GMX_REAL_EPS) )
    {
        *ntprs = std::max(*ntprs, 2);
    }

    /* If both rmin, rmax are set, we need 3 tpr files at minimum */
    if (!gmx_within_tol(*rmax, rcoulomb, GMX_REAL_EPS) && !gmx_within_tol(*rmin, rcoulomb, GMX_REAL_EPS) )
    {
        *ntprs = std::max(*ntprs, 3);
    }

    if (old != *ntprs)
    {
        fprintf(stderr, "NOTE: Your rmin, rmax setting changed -ntpr to %d\n", *ntprs);
    }

    if (*ntprs > 1)
    {
        if (gmx_within_tol(*rmin, rcoulomb, GMX_REAL_EPS) && gmx_within_tol(rcoulomb, *rmax, GMX_REAL_EPS)) /* We have just a single rc */
        {
            fprintf(stderr, "WARNING: Resetting -ntpr to 1 since no Coulomb radius scaling is requested.\n"
                    "Please set rmin < rmax to test Coulomb radii in the [rmin, rmax] interval\n"
                    "with correspondingly adjusted PME grid settings\n");
            *ntprs = 1;
        }
    }

    /* Check whether max and min fraction are within required values */
    if (maxPMEfraction > 0.5 || maxPMEfraction < 0)
    {
        gmx_fatal(FARGS, "-max must be between 0 and 0.5");
    }
    if (minPMEfraction > 0.5 || minPMEfraction < 0)
    {
        gmx_fatal(FARGS, "-min must be between 0 and 0.5");
    }
    if (maxPMEfraction < minPMEfraction)
    {
        gmx_fatal(FARGS, "-max must be larger or equal to -min");
    }

    /* Check whether the number of steps - if it was set - has a reasonable value */
    if (bench_nsteps < 0)
    {
        gmx_fatal(FARGS, "Number of steps must be positive.");
    }

    if (bench_nsteps > 10000 || bench_nsteps < 100)
    {
        fprintf(stderr, "WARNING: steps=");
        fprintf(stderr, "%" GMX_PRId64, bench_nsteps);
        fprintf(stderr, ". Are you sure you want to perform so %s steps for each benchmark?\n", (bench_nsteps < 100) ? "few" : "many");
    }

    if (presteps < 0)
    {
        gmx_fatal(FARGS, "Cannot have a negative number of presteps.\n");
    }

    /* Check for rcoulomb scaling if more than one .tpr file is tested */
    if (*ntprs > 1)
    {
        if (*rmin/rcoulomb < 0.75 || *rmax/rcoulomb > 1.25)
        {
            fprintf(stderr, "WARNING: Applying extreme scaling factor. I hope you know what you are doing.\n");
        }
    }

    /* If a fixed number of PME nodes is set we do rcoulomb and PME gird tuning
     * only. We need to check whether the requested number of PME-only nodes
     * makes sense. */
    if (npme_fixed > -1)
    {
        /* No more than 50% of all nodes can be assigned as PME-only nodes. */
        if (2*npme_fixed > nnodes)
        {
            gmx_fatal(FARGS, "Cannot have more than %d PME-only ranks for a total of %d ranks (you chose %d).\n",
                      nnodes/2, nnodes, npme_fixed);
        }
        if ((npme_fixed > 0) && (5*npme_fixed < nnodes))
        {
            fprintf(stderr, "WARNING: Only %g percent of the ranks are assigned as PME-only ranks.\n",
                    (100.0*npme_fixed)/nnodes);
        }
        if (opt2parg_bSet("-min", npargs, pa) || opt2parg_bSet("-max", npargs, pa))
        {
            fprintf(stderr, "NOTE: The -min, -max, and -npme options have no effect when a\n"
                    "      fixed number of PME-only ranks is requested with -fix.\n");
        }
    }
}


/* Returns TRUE when "opt" is needed at launch time */
static gmx_bool is_launch_file(char *opt, gmx_bool bSet)
{
    if (0 == std::strncmp(opt, "-swap", 5))
    {
        return bSet;
    }

    /* Apart from the input .tpr and the output log files we need all options that
     * were set on the command line and that do not start with -b */
    if    (0 == std::strncmp(opt, "-b", 2) || 0 == std::strncmp(opt, "-s", 2)
           || 0 == std::strncmp(opt, "-err", 4) || 0 == std::strncmp(opt, "-p", 2) )
    {
        return FALSE;
    }

    return bSet;
}


/* Returns TRUE when "opt" defines a file which is needed for the benchmarks runs */
static gmx_bool is_bench_file(char *opt, gmx_bool bSet, gmx_bool bOptional, gmx_bool bIsOutput)
{
    /* Apart from the input .tpr, all files starting with "-b" are for
     * _b_enchmark files exclusively */
    if (0 == std::strncmp(opt, "-s", 2))
    {
        return FALSE;
    }

    if (0 == std::strncmp(opt, "-b", 2) || 0 == std::strncmp(opt, "-s", 2))
    {
        if (!bOptional || bSet)
        {
            return TRUE;
        }
        else
        {
            return FALSE;
        }
    }
    else
    {
        if (bIsOutput)
        {
            return FALSE;
        }
        else
        {
            if (bSet) /* These are additional input files like -cpi -ei */
            {
                return TRUE;
            }
            else
            {
                return FALSE;
            }
        }
    }
}


/* Adds 'buf' to 'str' */
static void add_to_string(char **str, const char *buf)
{
    int len;


    len = std::strlen(*str) + std::strlen(buf) + 1;
    srenew(*str, len);
    std::strcat(*str, buf);
}


/* Create the command line for the benchmark as well as for the real run */
static void create_command_line_snippets(
        gmx_bool  bAppendFiles,
        gmx_bool  bKeepAndNumCPT,
        gmx_bool  bResetHWay,
        int       presteps,
        int       nfile,
        t_filenm  fnm[],
        char     *cmd_args_bench[],  /* command line arguments for benchmark runs */
        char     *cmd_args_launch[], /* command line arguments for simulation run */
        char      extra_args[],      /* Add this to the end of the command line */
        char     *deffnm)            /* Default file names, or NULL if not set */
{
    int         i;
    char       *opt;
    const char *name;
    char        strbuf[STRLEN];


    /* strlen needs at least '\0' as a string: */
    snew(*cmd_args_bench, 1);
    snew(*cmd_args_launch, 1);
    *cmd_args_launch[0] = '\0';
    *cmd_args_bench[0]  = '\0';


    /*******************************************/
    /* 1. Process other command line arguments */
    /*******************************************/
    if (presteps > 0)
    {
        /* Add equilibration steps to benchmark options */
        sprintf(strbuf, "-resetstep %d ", presteps);
        add_to_string(cmd_args_bench, strbuf);
    }
    /* These switches take effect only at launch time */
    if (deffnm)
    {
        sprintf(strbuf, "-deffnm %s ", deffnm);
        add_to_string(cmd_args_launch, strbuf);
    }
    if (FALSE == bAppendFiles)
    {
        add_to_string(cmd_args_launch, "-noappend ");
    }
    if (bKeepAndNumCPT)
    {
        add_to_string(cmd_args_launch, "-cpnum ");
    }
    if (bResetHWay)
    {
        add_to_string(cmd_args_launch, "-resethway ");
    }

    /********************/
    /* 2. Process files */
    /********************/
    for (i = 0; i < nfile; i++)
    {
        opt  = (char *)fnm[i].opt;
        name = opt2fn(opt, nfile, fnm);

        /* Strbuf contains the options, now let's sort out where we need that */
        sprintf(strbuf, "%s %s ", opt, name);

        if (is_bench_file(opt, opt2bSet(opt, nfile, fnm), is_optional(&fnm[i]), is_output(&fnm[i])) )
        {
            /* All options starting with -b* need the 'b' removed,
             * therefore overwrite strbuf */
            if (0 == std::strncmp(opt, "-b", 2))
            {
                sprintf(strbuf, "-%s %s ", &opt[2], name);
            }

            add_to_string(cmd_args_bench, strbuf);
        }

        if (is_launch_file(opt, opt2bSet(opt, nfile, fnm)) )
        {
            add_to_string(cmd_args_launch, strbuf);
        }
    }

    add_to_string(cmd_args_bench, extra_args);
    add_to_string(cmd_args_launch, extra_args);
}


/* Set option opt */
static void setopt(const char *opt, int nfile, t_filenm fnm[])
{
    int i;

    for (i = 0; (i < nfile); i++)
    {
        if (std::strcmp(opt, fnm[i].opt) == 0)
        {
            fnm[i].flag |= ffSET;
        }
    }
}


/* This routine inspects the tpr file and ...
 * 1. checks for output files that get triggered by a tpr option. These output
 *    files are marked as 'set' to allow for a proper cleanup after each
 *    tuning run.
 * 2. returns the PME:PP load ratio
 * 3. returns rcoulomb from the tpr */
static float inspect_tpr(int nfile, t_filenm fnm[], real *rcoulomb)
{
    gmx_bool     bTpi;      /* Is test particle insertion requested?          */
    gmx_bool     bFree;     /* Is a free energy simulation requested?         */
    gmx_bool     bNM;       /* Is a normal mode analysis requested?           */
    gmx_bool     bSwap;     /* Is water/ion position swapping requested?      */
    t_state      state;
    gmx_mtop_t   mtop;


    /* Check tpr file for options that trigger extra output files */
    t_inputrec  irInstance;
    t_inputrec *ir = &irInstance;
    read_tpx_state(opt2fn("-s", nfile, fnm), ir, &state, &mtop);
    bFree = (efepNO  != ir->efep );
    bNM   = (eiNM    == ir->eI   );
    bSwap = (eswapNO != ir->eSwapCoords);
    bTpi  = EI_TPI(ir->eI);

    /* Set these output files on the tuning command-line */
    if (ir->bPull)
    {
        setopt("-pf", nfile, fnm);
        setopt("-px", nfile, fnm);
    }
    if (bFree)
    {
        setopt("-dhdl", nfile, fnm);
    }
    if (bTpi)
    {
        setopt("-tpi", nfile, fnm);
        setopt("-tpid", nfile, fnm);
    }
    if (bNM)
    {
        setopt("-mtx", nfile, fnm);
    }
    if (bSwap)
    {
        setopt("-swap", nfile, fnm);
    }

    *rcoulomb = ir->rcoulomb;

    /* Return the estimate for the number of PME nodes */
    float npme = pme_load_estimate(&mtop, ir, state.box);
    return npme;
}


static void couple_files_options(int nfile, t_filenm fnm[])
{
    int      i;
    gmx_bool bSet, bBench;
    char    *opt;
    char     buf[20];


    for (i = 0; i < nfile; i++)
    {
        opt    = (char *)fnm[i].opt;
        bSet   = ((fnm[i].flag & ffSET) != 0);
        bBench = (0 == std::strncmp(opt, "-b", 2));

        /* Check optional files */
        /* If e.g. -eo is set, then -beo also needs to be set */
        if (is_optional(&fnm[i]) && bSet && !bBench)
        {
            sprintf(buf, "-b%s", &opt[1]);
            setopt(buf, nfile, fnm);
        }
        /* If -beo is set, then -eo also needs to be! */
        if (is_optional(&fnm[i]) && bSet && bBench)
        {
            sprintf(buf, "-%s", &opt[2]);
            setopt(buf, nfile, fnm);
        }
    }
}


#define BENCHSTEPS (1000)

int gmx_tune_pme(int argc, char *argv[])
{
    const char     *desc[] = {
        "For a given number [TT]-np[tt] or [TT]-ntmpi[tt] of ranks, [THISMODULE] systematically",
        "times [gmx-mdrun] with various numbers of PME-only ranks and determines",
        "which setting is fastest. It will also test whether performance can",
        "be enhanced by shifting load from the reciprocal to the real space",
        "part of the Ewald sum. ",
        "Simply pass your [REF].tpr[ref] file to [THISMODULE] together with other options",
        "for [gmx-mdrun] as needed.[PAR]",
        "[THISMODULE] needs to call [gmx-mdrun] and so requires that you",
        "specify how to call mdrun with the argument to the [TT]-mdrun[tt]",
        "parameter. Depending how you have built GROMACS, values such as",
        "'gmx mdrun', 'gmx_d mdrun', or 'mdrun_mpi' might be needed.[PAR]",
        "The program that runs MPI programs can be set in the environment variable",
        "MPIRUN (defaults to 'mpirun'). Note that for certain MPI frameworks,",
        "you need to provide a machine- or hostfile. This can also be passed",
        "via the MPIRUN variable, e.g.[PAR]",
        "[TT]export MPIRUN=\"/usr/local/mpirun -machinefile hosts\"[tt]",
        "Note that in such cases it is normally necessary to compile",
        "and/or run [THISMODULE] without MPI support, so that it can call",
        "the MPIRUN program.[PAR]",
        "Before doing the actual benchmark runs, [THISMODULE] will do a quick",
        "check whether [gmx-mdrun] works as expected with the provided parallel settings",
        "if the [TT]-check[tt] option is activated (the default).",
        "Please call [THISMODULE] with the normal options you would pass to",
        "[gmx-mdrun] and add [TT]-np[tt] for the number of ranks to perform the",
        "tests on, or [TT]-ntmpi[tt] for the number of threads. You can also add [TT]-r[tt]",
        "to repeat each test several times to get better statistics. [PAR]",
        "[THISMODULE] can test various real space / reciprocal space workloads",
        "for you. With [TT]-ntpr[tt] you control how many extra [REF].tpr[ref] files will be",
        "written with enlarged cutoffs and smaller Fourier grids respectively.",
        "Typically, the first test (number 0) will be with the settings from the input",
        "[REF].tpr[ref] file; the last test (number [TT]ntpr[tt]) will have the Coulomb cutoff",
        "specified by [TT]-rmax[tt] with a somewhat smaller PME grid at the same time. ",
        "In this last test, the Fourier spacing is multiplied with [TT]rmax[tt]/rcoulomb. ",
        "The remaining [REF].tpr[ref] files will have equally-spaced Coulomb radii (and Fourier "
        "spacings) between these extremes. [BB]Note[bb] that you can set [TT]-ntpr[tt] to 1",
        "if you just seek the optimal number of PME-only ranks; in that case",
        "your input [REF].tpr[ref] file will remain unchanged.[PAR]",
        "For the benchmark runs, the default of 1000 time steps should suffice for most",
        "MD systems. The dynamic load balancing needs about 100 time steps",
        "to adapt to local load imbalances, therefore the time step counters",
        "are by default reset after 100 steps. For large systems (>1M atoms), as well as ",
        "for a higher accuracy of the measurements, you should set [TT]-resetstep[tt] to a higher value.",
        "From the 'DD' load imbalance entries in the md.log output file you",
        "can tell after how many steps the load is sufficiently balanced. Example call:[PAR]"
        "[TT]gmx tune_pme -np 64 -s protein.tpr -launch[tt][PAR]",
        "After calling [gmx-mdrun] several times, detailed performance information",
        "is available in the output file [TT]perf.out[tt].",
        "[BB]Note[bb] that during the benchmarks, a couple of temporary files are written",
        "(options [TT]-b*[tt]), these will be automatically deleted after each test.[PAR]",
        "If you want the simulation to be started automatically with the",
        "optimized parameters, use the command line option [TT]-launch[tt].[PAR]",
        "Basic support for GPU-enabled [TT]mdrun[tt] exists. Give a string containing the IDs",
        "of the GPUs that you wish to use in the optimization in the [TT]-gpu_id[tt]",
        "command-line argument. This works exactly like [TT]mdrun -gpu_id[tt], does not imply a mapping,",
        "and merely declares the eligible set of GPU devices. [TT]gmx-tune_pme[tt] will construct calls to",
        "mdrun that use this set appropriately. [TT]gmx-tune_pme[tt] does not support",
        "[TT]-gputasks[tt].[PAR]",
    };

    int             nnodes         = 1;
    int             repeats        = 2;
    int             pmeentries     = 0; /* How many values for -npme do we actually test for each tpr file */
    real            maxPMEfraction = 0.50;
    real            minPMEfraction = 0.25;
    int             maxPMEnodes, minPMEnodes;
    float           guessPMEratio;                    /* guessed PME:PP ratio based on the tpr file */
    float           guessPMEnodes;
    int             npme_fixed     = -2;              /* If >= -1, use only this number
                                                       * of PME-only nodes                */
    int             ntprs          = 0;
    real            rmin           = 0.0, rmax = 0.0; /* min and max value for rcoulomb if scaling is requested */
    real            rcoulomb       = -1.0;            /* Coulomb radius as set in .tpr file */
    gmx_bool        bScaleRvdw     = TRUE;
    gmx_int64_t     bench_nsteps   = BENCHSTEPS;
    gmx_int64_t     new_sim_nsteps = -1;   /* -1 indicates: not set by the user */
    gmx_int64_t     cpt_steps      = 0;    /* Step counter in .cpt input file   */
    int             presteps       = 1500; /* Do a full cycle reset after presteps steps */
    gmx_bool        bOverwrite     = FALSE, bKeepTPR;
    gmx_bool        bLaunch        = FALSE;
    char           *ExtraArgs      = nullptr;
    char          **tpr_names      = nullptr;
    const char     *simulation_tpr = nullptr;
    char           *deffnm         = nullptr;
    int             best_npme, best_tpr;
    int             sim_part = 1; /* For benchmarks with checkpoint files */
    char            bbuf[STRLEN];


    /* Default program names if nothing else is found */
    char               *cmd_mpirun = nullptr, *cmd_mdrun = nullptr;
    char               *cmd_args_bench, *cmd_args_launch;
    char               *cmd_np           = nullptr;

    /* IDs of GPUs that are eligible for computation */
    char               *eligible_gpu_ids = nullptr;

    t_perf            **perfdata = nullptr;
    t_inputinfo        *info;
    int                 i;
    FILE               *fp;

    /* Print out how long the tuning took */
    double          seconds;

    static t_filenm fnm[] = {
        /* tune_pme */
        { efOUT, "-p",      "perf",     ffWRITE },
        { efLOG, "-err",    "bencherr", ffWRITE },
        { efTPR, "-so",     "tuned",    ffWRITE },
        /* mdrun: */
        { efTPR, nullptr,      nullptr,       ffREAD },
        { efTRN, "-o",      nullptr,       ffWRITE },
        { efCOMPRESSED, "-x", nullptr,     ffOPTWR },
        { efCPT, "-cpi",    nullptr,       ffOPTRD },
        { efCPT, "-cpo",    nullptr,       ffOPTWR },
        { efSTO, "-c",      "confout",  ffWRITE },
        { efEDR, "-e",      "ener",     ffWRITE },
        { efLOG, "-g",      "md",       ffWRITE },
        { efXVG, "-dhdl",   "dhdl",     ffOPTWR },
        { efXVG, "-field",  "field",    ffOPTWR },
        { efXVG, "-table",  "table",    ffOPTRD },
        { efXVG, "-tablep", "tablep",   ffOPTRD },
        { efXVG, "-tableb", "table",    ffOPTRD },
        { efTRX, "-rerun",  "rerun",    ffOPTRD },
        { efXVG, "-tpi",    "tpi",      ffOPTWR },
        { efXVG, "-tpid",   "tpidist",  ffOPTWR },
        { efEDI, "-ei",     "sam",      ffOPTRD },
        { efXVG, "-eo",     "edsam",    ffOPTWR },
        { efXVG, "-devout", "deviatie", ffOPTWR },
        { efXVG, "-runav",  "runaver",  ffOPTWR },
        { efXVG, "-px",     "pullx",    ffOPTWR },
        { efXVG, "-pf",     "pullf",    ffOPTWR },
        { efXVG, "-ro",     "rotation", ffOPTWR },
        { efLOG, "-ra",     "rotangles", ffOPTWR },
        { efLOG, "-rs",     "rotslabs", ffOPTWR },
        { efLOG, "-rt",     "rottorque", ffOPTWR },
        { efMTX, "-mtx",    "nm",       ffOPTWR },
        { efXVG, "-swap",   "swapions", ffOPTWR },
        /* Output files that are deleted after each benchmark run */
        { efTRN, "-bo",     "bench",    ffWRITE },
        { efXTC, "-bx",     "bench",    ffWRITE },
        { efCPT, "-bcpo",   "bench",    ffWRITE },
        { efSTO, "-bc",     "bench",    ffWRITE },
        { efEDR, "-be",     "bench",    ffWRITE },
        { efLOG, "-bg",     "bench",    ffWRITE },
        { efXVG, "-beo",    "benchedo", ffOPTWR },
        { efXVG, "-bdhdl",  "benchdhdl", ffOPTWR },
        { efXVG, "-bfield", "benchfld", ffOPTWR },
        { efXVG, "-btpi",   "benchtpi", ffOPTWR },
        { efXVG, "-btpid",  "benchtpid", ffOPTWR },
        { efXVG, "-bdevout", "benchdev", ffOPTWR },
        { efXVG, "-brunav", "benchrnav", ffOPTWR },
        { efXVG, "-bpx",    "benchpx",  ffOPTWR },
        { efXVG, "-bpf",    "benchpf",  ffOPTWR },
        { efXVG, "-bro",    "benchrot", ffOPTWR },
        { efLOG, "-bra",    "benchrota", ffOPTWR },
        { efLOG, "-brs",    "benchrots", ffOPTWR },
        { efLOG, "-brt",    "benchrott", ffOPTWR },
        { efMTX, "-bmtx",   "benchn",   ffOPTWR },
        { efNDX, "-bdn",    "bench",    ffOPTWR },
        { efXVG, "-bswap",  "benchswp", ffOPTWR }
    };

    gmx_bool        bThreads     = FALSE;

    int             nthreads = 1;

    const char     *procstring[] =
    { nullptr, "np", "n", "none", nullptr };
    const char     *npmevalues_opt[] =
    { nullptr, "auto", "all", "subset", nullptr };

    gmx_bool          bAppendFiles          = TRUE;
    gmx_bool          bKeepAndNumCPT        = FALSE;
    gmx_bool          bResetCountersHalfWay = FALSE;
    gmx_bool          bBenchmark            = TRUE;
    gmx_bool          bCheck                = TRUE;

    gmx_output_env_t *oenv = nullptr;

    t_pargs           pa[] = {
        /***********************/
        /* tune_pme options: */
        /***********************/
        { "-mdrun",    FALSE, etSTR,  {&cmd_mdrun},
          "Command line to run a simulation, e.g. 'gmx mdrun' or 'mdrun_mpi'" },
        { "-np",       FALSE, etINT,  {&nnodes},
          "Number of ranks to run the tests on (must be > 2 for separate PME ranks)" },
        { "-npstring", FALSE, etENUM, {procstring},
          "Name of the [TT]$MPIRUN[tt] option that specifies the number of ranks to use ('np', or 'n'; use 'none' if there is no such option)" },
        { "-ntmpi",    FALSE, etINT,  {&nthreads},
          "Number of MPI-threads to run the tests on (turns MPI & mpirun off)"},
        { "-r",        FALSE, etINT,  {&repeats},
          "Repeat each test this often" },
        { "-max",      FALSE, etREAL, {&maxPMEfraction},
          "Max fraction of PME ranks to test with" },
        { "-min",      FALSE, etREAL, {&minPMEfraction},
          "Min fraction of PME ranks to test with" },
        { "-npme",     FALSE, etENUM, {npmevalues_opt},
          "Within -min and -max, benchmark all possible values for [TT]-npme[tt], or just a reasonable subset. "
          "Auto neglects -min and -max and chooses reasonable values around a guess for npme derived from the .tpr"},
        { "-fix",      FALSE, etINT,  {&npme_fixed},
          "If >= -1, do not vary the number of PME-only ranks, instead use this fixed value and only vary rcoulomb and the PME grid spacing."},
        { "-rmax",     FALSE, etREAL, {&rmax},
          "If >0, maximal rcoulomb for -ntpr>1 (rcoulomb upscaling results in fourier grid downscaling)" },
        { "-rmin",     FALSE, etREAL, {&rmin},
          "If >0, minimal rcoulomb for -ntpr>1" },
        { "-scalevdw",  FALSE, etBOOL, {&bScaleRvdw},
          "Scale rvdw along with rcoulomb"},
        { "-ntpr",     FALSE, etINT,  {&ntprs},
          "Number of [REF].tpr[ref] files to benchmark. Create this many files with different rcoulomb scaling factors depending on -rmin and -rmax. "
          "If < 1, automatically choose the number of [REF].tpr[ref] files to test" },
        { "-steps",    FALSE, etINT64, {&bench_nsteps},
          "Take timings for this many steps in the benchmark runs" },
        { "-resetstep", FALSE, etINT,  {&presteps},
          "Let dlb equilibrate this many steps before timings are taken (reset cycle counters after this many steps)" },
        { "-nsteps",   FALSE, etINT64, {&new_sim_nsteps},
          "If non-negative, perform this many steps in the real run (overwrites nsteps from [REF].tpr[ref], add [REF].cpt[ref] steps)" },
        { "-launch",   FALSE, etBOOL, {&bLaunch},
          "Launch the real simulation after optimization" },
        { "-bench",    FALSE, etBOOL, {&bBenchmark},
          "Run the benchmarks or just create the input [REF].tpr[ref] files?" },
        { "-check",    FALSE, etBOOL, {&bCheck},
          "Before the benchmark runs, check whether mdrun works in parallel" },
        { "-gpu_id",   FALSE, etSTR,  {&eligible_gpu_ids},
          "List of unique GPU device IDs that are eligible for use" },
        /******************/
        /* mdrun options: */
        /******************/
        /* We let tune_pme parse and understand these options, because we need to
         * prevent that they appear on the mdrun command line for the benchmarks */
        { "-append",   FALSE, etBOOL, {&bAppendFiles},
          "Append to previous output files when continuing from checkpoint instead of adding the simulation part number to all file names (for launch only)" },
        { "-cpnum",    FALSE, etBOOL, {&bKeepAndNumCPT},
          "Keep and number checkpoint files (launch only)" },
        { "-deffnm",   FALSE, etSTR,  {&deffnm},
          "Set the default filenames (launch only)" },
        { "-resethway", FALSE, etBOOL, {&bResetCountersHalfWay},
          "HIDDENReset the cycle counters after half the number of steps or halfway [TT]-maxh[tt] (launch only)" }
    };

#define NFILE asize(fnm)

    seconds = gmx_gettime();

    if (!parse_common_args(&argc, argv, PCA_NOEXIT_ON_ARGS,
                           NFILE, fnm, asize(pa), pa, asize(desc), desc,
                           0, nullptr, &oenv))
    {
        return 0;
    }

    // procstring[0] is used inside two different conditionals further down
    GMX_RELEASE_ASSERT(procstring[0] != nullptr, "Options inconsistency; procstring[0] is NULL");

    /* Store the remaining unparsed command line entries in a string which
     * is then attached to the mdrun command line */
    snew(ExtraArgs, 1);
    ExtraArgs[0] = '\0';
    for (i = 1; i < argc; i++) /* argc will now be 1 if everything was understood */
    {
        add_to_string(&ExtraArgs, argv[i]);
        add_to_string(&ExtraArgs, " ");
    }

    if (opt2parg_bSet("-ntmpi", asize(pa), pa))
    {
        bThreads = TRUE;
        if (opt2parg_bSet("-npstring", asize(pa), pa))
        {
            fprintf(stderr, "WARNING: -npstring has no effect when using threads.\n");
        }

        if (nnodes > 1)
        {
            gmx_fatal(FARGS, "Can't run multi-threaded MPI simulation yet!");
        }
        /* and now we just set this; a bit of an ugly hack*/
        nnodes = nthreads;
    }
    /* Check for PME:PP ratio and whether tpr triggers additional output files */
    guessPMEratio = inspect_tpr(NFILE, fnm, &rcoulomb);

    /* Automatically set -beo options if -eo is set etc. */
    couple_files_options(NFILE, fnm);

    /* Construct the command line arguments for benchmark runs
     * as well as for the simulation run */
    if (bThreads)
    {
        sprintf(bbuf, " -ntmpi %d ", nthreads);
    }
    else
    {
        /* This string will be used for MPI runs and will appear after the
         * mpirun command. */
        if (std::strcmp(procstring[0], "none") != 0)
        {
            sprintf(bbuf, " -%s %d ", procstring[0], nnodes);
        }
        else
        {
            sprintf(bbuf, " ");
        }
    }

    cmd_np = bbuf;

    create_command_line_snippets(bAppendFiles, bKeepAndNumCPT, bResetCountersHalfWay, presteps,
                                 NFILE, fnm, &cmd_args_bench, &cmd_args_launch, ExtraArgs, deffnm);

    /* Prepare to use checkpoint file if requested */
    sim_part = 1;
    if (opt2bSet("-cpi", NFILE, fnm))
    {
        const char *filename = opt2fn("-cpi", NFILE, fnm);
        int         cpt_sim_part;
        read_checkpoint_part_and_step(filename,
                                      &cpt_sim_part, &cpt_steps);
        if (cpt_sim_part == 0)
        {
            gmx_fatal(FARGS, "Checkpoint file %s could not be read!", filename);
        }
        /* tune_pme will run the next part of the simulation */
        sim_part = cpt_sim_part + 1;
    }

    /* Open performance output file and write header info */
    fp = gmx_ffopen(opt2fn("-p", NFILE, fnm), "w");

    /* Make a quick consistency check of command line parameters */
    check_input(nnodes, repeats, &ntprs, &rmin, rcoulomb, &rmax,
                maxPMEfraction, minPMEfraction, npme_fixed,
                bench_nsteps, fnm, NFILE, sim_part, presteps,
                asize(pa), pa);

    /* Determine the maximum and minimum number of PME nodes to test,
     * the actual list of settings is build in do_the_tests(). */
    if ((nnodes > 2) && (npme_fixed < -1))
    {
        if (0 == std::strcmp(npmevalues_opt[0], "auto"))
        {
            /* Determine the npme range automatically based on the PME:PP load guess */
            if (guessPMEratio > 1.0)
            {
                /* More PME than PP work, probably we do not need separate PME nodes at all! */
                maxPMEnodes = nnodes/2;
                minPMEnodes = nnodes/2;
            }
            else
            {
                /* PME : PP load is in the range 0..1, let's test around the guess */
                guessPMEnodes = static_cast<int>(nnodes/(1.0 + 1.0/guessPMEratio));
                minPMEnodes   = static_cast<int>(std::floor(0.7*guessPMEnodes));
                maxPMEnodes   = static_cast<int>(std::ceil(1.6*guessPMEnodes));
                maxPMEnodes   = std::min(maxPMEnodes, nnodes/2);
            }
        }
        else
        {
            /* Determine the npme range based on user input */
            maxPMEnodes = static_cast<int>(std::floor(maxPMEfraction*nnodes));
            minPMEnodes = std::max(static_cast<int>(std::floor(minPMEfraction*nnodes)), 0);
            fprintf(stdout, "Will try runs with %d ", minPMEnodes);
            if (maxPMEnodes != minPMEnodes)
            {
                fprintf(stdout, "- %d ", maxPMEnodes);
            }
            fprintf(stdout, "PME-only ranks.\n  Note that the automatic number of PME-only ranks and no separate PME ranks are always tested.\n");
        }
    }
    else
    {
        maxPMEnodes = 0;
        minPMEnodes = 0;
    }

    /* Get the commands we need to set up the runs from environment variables */
    get_program_paths(bThreads, &cmd_mpirun, &cmd_mdrun);
    if (bBenchmark && repeats > 0)
    {
        check_mdrun_works(bThreads, cmd_mpirun, cmd_np, cmd_mdrun, nullptr != eligible_gpu_ids);
    }

    /* Print some header info to file */
    sep_line(fp);
    fprintf(fp, "\n      P E R F O R M A N C E   R E S U L T S\n");
    sep_line(fp);
    fprintf(fp, "%s for GROMACS %s\n", output_env_get_program_display_name(oenv),
            gmx_version());
    if (!bThreads)
    {
        fprintf(fp, "Number of ranks         : %d\n", nnodes);
        fprintf(fp, "The mpirun command is   : %s\n", cmd_mpirun);
        if (std::strcmp(procstring[0], "none") != 0)
        {
            fprintf(fp, "Passing # of ranks via  : -%s\n", procstring[0]);
        }
        else
        {
            fprintf(fp, "Not setting number of ranks in system call\n");
        }
    }
    else
    {
        fprintf(fp, "Number of threads       : %d\n", nnodes);
    }

    fprintf(fp, "The mdrun  command is   : %s\n", cmd_mdrun);
    fprintf(fp, "mdrun args benchmarks   : %s\n", cmd_args_bench);
    fprintf(fp, "Benchmark steps         : ");
    fprintf(fp, "%" GMX_PRId64, bench_nsteps);
    fprintf(fp, "\n");
    fprintf(fp, "dlb equilibration steps : %d\n", presteps);
    if (sim_part > 1)
    {
        fprintf(fp, "Checkpoint time step    : ");
        fprintf(fp, "%" GMX_PRId64, cpt_steps);
        fprintf(fp, "\n");
    }
    fprintf(fp, "mdrun args at launchtime: %s\n", cmd_args_launch);

    if (new_sim_nsteps >= 0)
    {
        bOverwrite = TRUE;
        fprintf(stderr, "Note: Simulation input file %s will have ", opt2fn("-so", NFILE, fnm));
        fprintf(stderr, "%" GMX_PRId64, new_sim_nsteps+cpt_steps);
        fprintf(stderr, " steps.\n");
        fprintf(fp, "Simulation steps        : ");
        fprintf(fp, "%" GMX_PRId64, new_sim_nsteps);
        fprintf(fp, "\n");
    }
    if (repeats > 1)
    {
        fprintf(fp, "Repeats for each test   : %d\n", repeats);
    }

    if (npme_fixed >= -1)
    {
        fprintf(fp, "Fixing -npme at         : %d\n", npme_fixed);
    }

    fprintf(fp, "Input file              : %s\n", opt2fn("-s", NFILE, fnm));
    fprintf(fp, "   PME/PP load estimate : %g\n", guessPMEratio);

    /* Allocate memory for the inputinfo struct: */
    snew(info, 1);
    info->nr_inputfiles = ntprs;
    for (i = 0; i < ntprs; i++)
    {
        snew(info->rcoulomb, ntprs);
        snew(info->rvdw, ntprs);
        snew(info->rlist, ntprs);
        snew(info->nkx, ntprs);
        snew(info->nky, ntprs);
        snew(info->nkz, ntprs);
        snew(info->fsx, ntprs);
        snew(info->fsy, ntprs);
        snew(info->fsz, ntprs);
    }
    /* Make alternative tpr files to test: */
    snew(tpr_names, ntprs);
    for (i = 0; i < ntprs; i++)
    {
        snew(tpr_names[i], STRLEN);
    }

    /* It can be that ntprs is reduced by make_benchmark_tprs if not enough
     * different grids could be found. */
    make_benchmark_tprs(opt2fn("-s", NFILE, fnm), tpr_names, bench_nsteps+presteps,
                        cpt_steps, rmin, rmax, bScaleRvdw, &ntprs, info, fp);

    /********************************************************************************/
    /* Main loop over all scenarios we need to test: tpr files, PME nodes, repeats  */
    /********************************************************************************/
    snew(perfdata, ntprs);
    if (bBenchmark)
    {
        GMX_RELEASE_ASSERT(npmevalues_opt[0] != nullptr, "Options inconsistency; npmevalues_opt[0] is NULL");
        do_the_tests(fp, tpr_names, maxPMEnodes, minPMEnodes, npme_fixed, npmevalues_opt[0], perfdata, &pmeentries,
                     repeats, nnodes, ntprs, bThreads, cmd_mpirun, cmd_np, cmd_mdrun,
                     cmd_args_bench, fnm, NFILE, presteps, cpt_steps, bCheck, eligible_gpu_ids);

        fprintf(fp, "\nTuning took%8.1f minutes.\n", (gmx_gettime()-seconds)/60.0);

        /* Analyse the results and give a suggestion for optimal settings: */
        bKeepTPR = analyze_data(fp, opt2fn("-p", NFILE, fnm), perfdata, nnodes, ntprs, pmeentries,
                                repeats, info, &best_tpr, &best_npme);

        /* Take the best-performing tpr file and enlarge nsteps to original value */
        if (bKeepTPR && !bOverwrite)
        {
            simulation_tpr = opt2fn("-s", NFILE, fnm);
        }
        else
        {
            simulation_tpr = opt2fn("-so", NFILE, fnm);
            modify_PMEsettings(bOverwrite ? (new_sim_nsteps+cpt_steps) : info->orig_sim_steps,
                               info->orig_init_step, tpr_names[best_tpr], simulation_tpr);
        }

        /* Let's get rid of the temporary benchmark input files */
        for (i = 0; i < ntprs; i++)
        {
            fprintf(stdout, "Deleting temporary benchmark input file %s\n", tpr_names[i]);
            remove(tpr_names[i]);
        }

        /* Now start the real simulation if the user requested it ... */
        launch_simulation(bLaunch, fp, bThreads, cmd_mpirun, cmd_np, cmd_mdrun,
                          cmd_args_launch, simulation_tpr, best_npme, eligible_gpu_ids);
    }
    gmx_ffclose(fp);

    /* ... or simply print the performance results to screen: */
    if (!bLaunch)
    {
        finalize(opt2fn("-p", NFILE, fnm));
    }

    return 0;
}
