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
 * Copyright (c) 2001-2008, The GROMACS development team,
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
 * 
 * And Hey:
 * Gallium Rubidium Oxygen Manganese Argon Carbon Silicon
 */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif


#include <time.h>
#ifdef HAVE_SYS_TIME_H
#include <sys/time.h>
#endif



#include "statutil.h"
#include "typedefs.h"
#include "smalloc.h"
#include "vec.h"
#include "copyrite.h"
#include "statutil.h"
#include "tpxio.h"
#include "string2.h"
#include "readinp.h"
#include "calcgrid.h"
#include "checkpoint.h"
#include "gmx_ana.h"
#include "names.h"



enum {
  ddnoSEL, ddnoINTERLEAVE, ddnoPP_PME, ddnoCARTESIAN, ddnoNR
};

/* Enum for situations that can occur during log file parsing, the
 * corresponding string entries can be found in do_the_tests() in
 * const char* ParseLog[] */
enum {
    eParselogOK,
    eParselogNotFound,
    eParselogNoPerfData,
    eParselogTerm,
    eParselogResetProblem,
    eParselogNoDDGrid,
    eParselogTPXVersion,
    eParselogNotParallel,
    eParselogFatal,
    eParselogNr
};


typedef struct
{
    int  nPMEnodes;       /* number of PME-only nodes used in this test */
    int  nx, ny, nz;      /* DD grid */
    int  guessPME;        /* if nPMEnodes == -1, this is the guessed number of PME nodes */
    double *Gcycles;      /* This can contain more than one value if doing multiple tests */
    double Gcycles_Av;
    float *ns_per_day;
    float ns_per_day_Av;
    float *PME_f_load;    /* PME mesh/force load average*/
    float PME_f_load_Av;  /* Average average ;) ... */
    char *mdrun_cmd_line; /* Mdrun command line used for this test */
} t_perf;


typedef struct
{
    int  nr_inputfiles;         /* The number of tpr and mdp input files */
    gmx_large_int_t orig_sim_steps;  /* Number of steps to be done in the real simulation */
    gmx_large_int_t orig_init_step;  /* Init step for the real simulation */
    real *rcoulomb;             /* The coulomb radii [0...nr_inputfiles] */
    real *rvdw;                 /* The vdW radii */
    real *rlist;                /* Neighbourlist cutoff radius */
    real *rlistlong;
    int  *nkx, *nky, *nkz;
    real *fsx, *fsy, *fsz;      /* Fourierspacing in x,y,z dimension */
} t_inputinfo;


static void sep_line(FILE *fp)
{
    fprintf(fp, "\n------------------------------------------------------------\n");
}


/* Wrapper for system calls */
static int gmx_system_call(char *command)
{
#ifdef GMX_NO_SYSTEM
    gmx_fatal(FARGS,"No calls to system(3) supported on this platform. Attempted to call:\n'%s'\n",command);
#else
    return ( system(command) );
#endif
}
 

/* Check if string starts with substring */
static gmx_bool str_starts(const char *string, const char *substring)
{
    return ( strncmp(string, substring, strlen(substring)) == 0);
}


static void cleandata(t_perf *perfdata, int test_nr)
{
    perfdata->Gcycles[test_nr]    = 0.0;
    perfdata->ns_per_day[test_nr] = 0.0;
    perfdata->PME_f_load[test_nr] = 0.0;
    
    return;
}


static gmx_bool is_equal(real a, real b)
{
    real diff, eps=1.0e-7;


    diff = a - b;

    if (diff < 0.0) diff = -diff;

    if (diff < eps)
        return TRUE;
    else
        return FALSE;
}


static void finalize(const char *fn_out)
{
    char buf[STRLEN];
    FILE *fp;


    fp = fopen(fn_out,"r");
    fprintf(stdout,"\n\n");

    while( fgets(buf,STRLEN-1,fp) != NULL )
    {
        fprintf(stdout,"%s",buf);
    }
    fclose(fp);
    fprintf(stdout,"\n\n");
    thanx(stderr);
}


enum {eFoundNothing, eFoundDDStr, eFoundAccountingStr, eFoundCycleStr};

static int parse_logfile(const char *logfile, const char *errfile,
        t_perf *perfdata, int test_nr, int presteps, gmx_large_int_t cpt_steps,
        int nnodes)
{
    FILE  *fp;
    char  line[STRLEN], dumstring[STRLEN], dumstring2[STRLEN];
    const char matchstrdd[]="Domain decomposition grid";
    const char matchstrcr[]="resetting all time and cycle counters";
    const char matchstrbal[]="Average PME mesh/force load:";
    const char matchstring[]="R E A L   C Y C L E   A N D   T I M E   A C C O U N T I N G";
    const char errSIG[]="signal, stopping at the next";
    int   iFound;
    int   procs;
    float  dum1,dum2,dum3;
    int   npme;
    gmx_large_int_t resetsteps=-1;
    gmx_bool  bFoundResetStr = FALSE;
    gmx_bool  bResetChecked  = FALSE;


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
        iFound = eFoundDDStr; /* Skip some case statements */

    while (fgets(line, STRLEN, fp) != NULL)
    {
        /* Remove leading spaces */
        ltrim(line);

        /* Check for TERM and INT signals from user: */
        if ( strstr(line, errSIG) != NULL )
        {
            fclose(fp);
            cleandata(perfdata, test_nr);
            return eParselogTerm;
        }
        
        /* Check whether cycle resetting  worked */
        if (presteps > 0 && !bFoundResetStr)
        {
            if (strstr(line, matchstrcr) != NULL)
            {
                sprintf(dumstring, "Step %s", gmx_large_int_pfmt);
                sscanf(line, dumstring, &resetsteps);
                bFoundResetStr = TRUE;
                if (resetsteps == presteps+cpt_steps)
                {
                    bResetChecked = TRUE;
                }
                else
                {
                    sprintf(dumstring , gmx_large_int_pfmt, resetsteps);
                    sprintf(dumstring2, gmx_large_int_pfmt, presteps+cpt_steps);
                    fprintf(stderr, "WARNING: Time step counters were reset at step %s,\n"
                                    "         though they were supposed to be reset at step %s!\n", 
                            dumstring, dumstring2);
                }
            }
        }

        /* Look for strings that appear in a certain order in the log file: */
        switch(iFound)
        {
            case eFoundNothing:
                /* Look for domain decomp grid and separate PME nodes: */
                if (str_starts(line, matchstrdd))
                {
                    sscanf(line, "Domain decomposition grid %d x %d x %d, separate PME nodes %d",
                            &(perfdata->nx), &(perfdata->ny), &(perfdata->nz), &npme);
                    if (perfdata->nPMEnodes == -1)
                        perfdata->guessPME = npme;
                    else if (perfdata->nPMEnodes != npme)
                        gmx_fatal(FARGS, "PME nodes from command line and output file are not identical");
                    iFound = eFoundDDStr;
                }
                /* Catch a few errors that might have occured: */
                else if (str_starts(line, "There is no domain decomposition for"))
                {
                    fclose(fp);
                    return eParselogNoDDGrid;
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
                /* Look for PME mesh/force balance (not necessarily present, though) */
                if (str_starts(line, matchstrbal))
                    sscanf(&line[strlen(matchstrbal)], "%f", &(perfdata->PME_f_load[test_nr]));
                /* Look for matchstring */
                if (str_starts(line, matchstring))
                    iFound = eFoundAccountingStr;
                break;
            case eFoundAccountingStr:
                /* Already found matchstring - look for cycle data */
                if (str_starts(line, "Total  "))
                {
                    sscanf(line,"Total %d %lf",&procs,&(perfdata->Gcycles[test_nr]));
                    iFound = eFoundCycleStr;
                }
                break;
            case eFoundCycleStr:
                /* Already found cycle data - look for remaining performance info and return */
                if (str_starts(line, "Performance:"))
                {
                    sscanf(line,"%s %f %f %f %f", dumstring, &dum1, &dum2, &(perfdata->ns_per_day[test_nr]), &dum3);
                    fclose(fp);
                    if (bResetChecked || presteps == 0)
                        return eParselogOK;
                    else 
                        return eParselogResetProblem;
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
        while (fgets(line, STRLEN, fp) != NULL)
        {
            if ( str_starts(line, "Fatal error:") )
            {
                if (fgets(line, STRLEN, fp) != NULL)
                    fprintf(stderr, "\nWARNING: A fatal error has occured during this benchmark:\n"
                                    "%s\n", line);
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
        FILE        *fp,
        const char  *fn,
        t_perf      **perfdata,
        int         nnodes,
        int         ntprs,
        int         ntests,
        int         nrepeats,
        t_inputinfo *info,
        int         *index_tpr,    /* OUT: Nr of mdp file with best settings */
        int         *npme_optimal) /* OUT: Optimal number of PME nodes */
{
    int  i,j,k;
    int line=0, line_win=-1;
    int  k_win=-1, i_win=-1, winPME;
    double s=0.0;  /* standard deviation */
    t_perf *pd;
    char strbuf[STRLEN];
    char str_PME_f_load[13];
    gmx_bool bCanUseOrigTPR;
    gmx_bool bRefinedCoul, bRefinedVdW, bRefinedGrid;


    if (nrepeats > 1)
    {
        sep_line(fp);
        fprintf(fp, "Summary of successful runs:\n");
        fprintf(fp, "Line tpr PME nodes  Gcycles Av.     Std.dev.       ns/day        PME/f");
        if (nnodes > 1)
            fprintf(fp, "    DD grid");
        fprintf(fp, "\n");
    }


    for (k=0; k<ntprs; k++)
    {
        for (i=0; i<ntests; i++)
        {
            /* Select the right dataset: */
            pd = &(perfdata[k][i]);

            pd->Gcycles_Av    = 0.0;
            pd->PME_f_load_Av = 0.0;
            pd->ns_per_day_Av = 0.0;

            if (pd->nPMEnodes == -1)
                sprintf(strbuf, "(%3d)", pd->guessPME);
            else
                sprintf(strbuf, "     ");

            /* Get the average run time of a setting */
            for (j=0; j<nrepeats; j++)
            {
                pd->Gcycles_Av    += pd->Gcycles[j];
                pd->PME_f_load_Av += pd->PME_f_load[j];
            }
            pd->Gcycles_Av    /= nrepeats;
            pd->PME_f_load_Av /= nrepeats;

            for (j=0; j<nrepeats; j++)
            {
                if (pd->ns_per_day[j] > 0.0)
                    pd->ns_per_day_Av += pd->ns_per_day[j];
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
                sprintf(str_PME_f_load, "%12.3f", pd->PME_f_load_Av);
            else
                sprintf(str_PME_f_load, "%s", "         -  ");


            /* We assume we had a successful run if both averages are positive */
            if (pd->Gcycles_Av > 0.0 && pd->ns_per_day_Av > 0.0)
            {
                /* Output statistics if repeats were done */
                if (nrepeats > 1)
                {
                    /* Calculate the standard deviation */
                    s = 0.0;
                    for (j=0; j<nrepeats; j++)
                        s += pow( pd->Gcycles[j] - pd->Gcycles_Av, 2 );
                    s /= (nrepeats - 1);
                    s = sqrt(s);

                    fprintf(fp, "%4d %3d %4d%s %12.3f %12.3f %12.3f %s",
                            line, k, pd->nPMEnodes, strbuf, pd->Gcycles_Av, s,
                            pd->ns_per_day_Av, str_PME_f_load);
                    if (nnodes > 1)
                        fprintf(fp, "  %3d %3d %3d", pd->nx, pd->ny, pd->nz);
                    fprintf(fp, "\n");
                }
                /* Store the index of the best run found so far in 'winner': */
                if ( (k_win == -1) || (pd->Gcycles_Av < perfdata[k_win][i_win].Gcycles_Av) )
                {
                    k_win = k;
                    i_win = i;
                    line_win = line;
                }
                line++;
            }
        }
    }

    if (k_win == -1)
        gmx_fatal(FARGS, "None of the runs was successful! Check %s for problems.", fn);

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
            sprintf(strbuf, "%s", "the automatic number of PME nodes");
        else
            sprintf(strbuf, "%d PME nodes", winPME);
    }
    fprintf(fp, "Best performance was achieved with %s", strbuf);
    if ((nrepeats > 1) && (ntests > 1))
        fprintf(fp, " (see line %d)", line_win);
    fprintf(fp, "\n");

    /* Only mention settings if they were modified: */
    bRefinedCoul = !is_equal(info->rcoulomb[k_win], info->rcoulomb[0]);
    bRefinedVdW  = !is_equal(info->rvdw[k_win]    , info->rvdw[0]    );
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
        fprintf(fp, "   New Coulomb radius: %f nm (was %f nm)\n", info->rcoulomb[k_win], info->rcoulomb[0]);

    if (bRefinedVdW)
        fprintf(fp, "   New Van der Waals radius: %f nm (was %f nm)\n", info->rvdw[k_win], info->rvdw[0]);

    if (bRefinedGrid)
        fprintf(fp, "   New Fourier grid xyz: %d %d %d (was %d %d %d)\n", info->nkx[k_win], info->nky[k_win], info->nkz[k_win],
                                                                          info->nkx[0], info->nky[0], info->nkz[0]);

    if (bCanUseOrigTPR && ntprs > 1)
        fprintf(fp, "and original PME settings.\n");
    
    fflush(fp);
    
    /* Return the index of the mdp file that showed the highest performance
     * and the optimal number of PME nodes */
    *index_tpr    = k_win; 
    *npme_optimal = winPME;
    
    return bCanUseOrigTPR;
}


/* Get the commands we need to set up the runs from environment variables */
static void get_program_paths(gmx_bool bThreads, char *cmd_mpirun[], char cmd_np[],
                              char *cmd_mdrun[], int repeats)
{
    char *command=NULL;
    char *cp;
    char *cp2;
    char line[STRLEN];
    FILE *fp;
    const char def_mpirun[] = "mpirun";
    const char def_mdrun[]  = "mdrun";
    const char filename[]   = "benchtest.log";
    const char match_mpi[]  = "NNODES=";
    const char match_mdrun[]= "Program: ";
    const char empty_mpirun[] = "";
    gmx_bool  bMdrun = FALSE;
    gmx_bool  bMPI   = FALSE;
    

    /* Get the commands we need to set up the runs from environment variables */
    if (!bThreads)
    {
        if ( (cp = getenv("MPIRUN")) != NULL)
            *cmd_mpirun = strdup(cp);
        else
            *cmd_mpirun = strdup(def_mpirun);
    }
    else
    {
        *cmd_mpirun = strdup(empty_mpirun);
    }
     
    if ( (cp = getenv("MDRUN" )) != NULL )
        *cmd_mdrun  = strdup(cp);
    else
        *cmd_mdrun  = strdup(def_mdrun);


    /* If no simulations have to be performed, we are done here */
    if (repeats <= 0)
        return;

    /* Run a small test to see whether mpirun + mdrun work  */
    fprintf(stdout, "Making sure that mdrun can be executed. ");
    if (bThreads)
    {
        snew(command, strlen(*cmd_mdrun) + strlen(cmd_np) + strlen(filename) + 50);
        sprintf(command, "%s%s-version -maxh 0.001 1> %s 2>&1", *cmd_mdrun, cmd_np, filename);
    }	
    else
    {
        snew(command, strlen(*cmd_mpirun) + strlen(cmd_np) + strlen(*cmd_mdrun) + strlen(filename) + 50);
        sprintf(command, "%s%s%s -version -maxh 0.001 1> %s 2>&1", *cmd_mpirun, cmd_np, *cmd_mdrun, filename);
    }
    fprintf(stdout, "Trying '%s' ... ", command);
    make_backup(filename);
    gmx_system_call(command);

    /* Check if we find the characteristic string in the output: */
    if (!gmx_fexist(filename))
        gmx_fatal(FARGS, "Output from test run could not be found.");

    fp = fopen(filename, "r");
    /* We need to scan the whole output file, since sometimes the queuing system
     * also writes stuff to stdout/err */
    while ( !feof(fp) )
    {
        cp2=fgets(line, STRLEN, fp);
        if (cp2!=NULL)
        {
            if ( str_starts(line, match_mdrun) )
                bMdrun = TRUE;
            if ( str_starts(line, match_mpi) )
                bMPI = TRUE;
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
                    *cmd_mdrun);
        }
    }
    else
    {
        if (bMdrun && !bMPI)
        {
            gmx_fatal(FARGS, "Need an MPI-enabled version of mdrun. This one\n"
                    "(%s)\n"
                    "seems to have been compiled without MPI support.",
                    *cmd_mdrun);
        }
    }

    if (!bMdrun)
    {
        gmx_fatal(FARGS, "Cannot execute mdrun. Please check %s for problems!",
                filename);
    }

    fprintf(stdout, "passed.\n");

    /* Clean up ... */
    remove(filename);
    sfree(command);
}


static void launch_simulation(
        gmx_bool bLaunch,           /* Should the simulation be launched? */
        FILE *fp,               /* General log file */
        gmx_bool bThreads,          /* whether to use threads */
        char *cmd_mpirun,       /* Command for mpirun */
        char *cmd_np,           /* Switch for -np or -nt or empty */
        char *cmd_mdrun,        /* Command for mdrun */
        char *args_for_mdrun,   /* Arguments for mdrun */
        const char *simulation_tpr,   /* This tpr will be simulated */
        int  nnodes,            /* Number of nodes to run on */
        int  nPMEnodes)         /* Number of PME nodes to use */
{
    char  *command;
    
    
    /* Make enough space for the system call command, 
     * (100 extra chars for -npme ... etc. options should suffice): */
    snew(command, strlen(cmd_mpirun)+strlen(cmd_mdrun)+strlen(cmd_np)+strlen(args_for_mdrun)+strlen(simulation_tpr)+100);
   
    /* Note that the -passall options requires args_for_mdrun to be at the end
     * of the command line string */
    if (bThreads)
    {
        sprintf(command, "%s%s-npme %d -s %s %s",
                cmd_mdrun, cmd_np, nPMEnodes, simulation_tpr, args_for_mdrun);
    }
    else 
    {
        sprintf(command, "%s%s%s -npme %d -s %s %s",
                cmd_mpirun, cmd_np, cmd_mdrun, nPMEnodes, simulation_tpr, args_for_mdrun);
    }
        
    fprintf(fp, "%s this command line to launch the simulation:\n\n%s", bLaunch? "Using":"Please use", command);
    sep_line(fp);
    fflush(fp);

    /* Now the real thing! */
    if (bLaunch)
    {
        fprintf(stdout, "\nLaunching simulation with best parameters now.\nExecuting '%s'", command);
        sep_line(stdout);
        fflush(stdout);
        gmx_system_call(command);
        thanx(fp);
    }
}


static void modify_PMEsettings(
        gmx_large_int_t simsteps,  /* Set this value as number of time steps */
        gmx_large_int_t init_step, /* Set this value as init_step */
        const char *fn_best_tpr,   /* tpr file with the best performance */
        const char *fn_sim_tpr)    /* name of tpr file to be launched */
{
    t_inputrec   *ir;
    t_state      state;
    gmx_mtop_t   mtop;
    char         buf[200];
    
    snew(ir,1);
    read_tpx_state(fn_best_tpr,ir,&state,NULL,&mtop);
        
    /* Reset nsteps and init_step to the value of the input .tpr file */
    ir->nsteps = simsteps;
    ir->init_step = init_step;
    
    /* Write the tpr file which will be launched */
    sprintf(buf, "Writing optimized simulation file %s with nsteps=%s.\n", fn_sim_tpr, gmx_large_int_pfmt);
    fprintf(stdout,buf,ir->nsteps);
    fflush(stdout);
    write_tpx_state(fn_sim_tpr,ir,&state,&mtop);
        
    sfree(ir);
}


typedef struct
{
    int  nkx, nky, nkz;      /* Fourier grid                                  */
    real fac;                /* actual scaling factor                         */
    real fs;                 /* Fourierspacing                                */
} t_pmegrid;


static void copy_grid(t_pmegrid *ingrid, t_pmegrid *outgrid)
{
    outgrid->nkx = ingrid->nkx;
    outgrid->nky = ingrid->nky;
    outgrid->nkz = ingrid->nkz;
    outgrid->fac = ingrid->fac;
    outgrid->fs  = ingrid->fs;
}

/* Removes entry 'index' from the t_pmegrid list */
static void remove_from_list(
        t_pmegrid gridlist[],
        int *nlist,           /* Length of the list                           */
        int index)            /* Index to remove from the list                */
{
    int j;


    for (j = index; j < (*nlist - 1); j++)
    {
        copy_grid(&gridlist[j+1], &gridlist[j]);
    }
    *nlist = *nlist - 1;
}


/* Returns the index of the least necessary grid in the list.
 *
 * This is the one where the following conditions hold for the scaling factor:
 * 1. this grid has the smallest distance to its neighboring grid (distance
 *    measured by fac) -> this condition is true for two grids at the same time
 * 2. this grid (of the two) has the smaller distance to the other neighboring
 *    grid
 *
 * The extreme grids (the ones with the smallest and largest
 * scaling factor) are never thrown away.
 */
static int get_throwaway_index(
        t_pmegrid gridlist[],
        int nlist)           /* Length of the list                           */
{
    int i,index = -1;
    real dist,mindist,d_left,d_right;


    /* Find the two factors with the smallest mutual distance */
    mindist = GMX_FLOAT_MAX;
    for (i = 1; i < nlist; i++)
    {
        dist = fabs(gridlist[i].fac - gridlist[i-1].fac);
        if (dist < mindist)
        {
            index = i;
            mindist = dist;
        }
    }
    /* index and index-1 have the smallest mutual distance */
    if (index == 1)
    {
        /* Never return the first index, i.e. index == 0 */
        ;
    }
    else
    {
        d_left  = fabs(gridlist[index-1].fac - gridlist[index-2].fac);
        d_right = fabs(gridlist[index+1].fac - gridlist[index  ].fac);

        /* Return the left index if its distance to its left neighbor is shorter
         * than the distance of the right index to its right neighbor */
        if (d_left < d_right)
            index--;
    }
    /* Never return the last index */
    if (index == nlist-1)
        index--;

    return index;
}


static void make_grid_list(
        real fmin,            /* minimum scaling factor (downscaling fac)     */
        real fmax,            /* maximum scaling factor (upscaling fac)       */
        int *ntprs,           /* Number of tpr files to test                  */
        matrix box,           /* The box                                      */
        t_pmegrid *griduse[], /* List of grids that have to be tested         */
        real fs)              /* Requested fourierspacing at a scaling factor
                                 of unity                                     */
{
    real req_fac,act_fac=0,act_fs,eps;
    rvec box_size;
    int i,jmin=0,d,ngridall=0;
    int nkx=0,nky=0,nkz=0;
    int nkx_old=0,nky_old=0,nkz_old=0;
    t_pmegrid *gridall;
    int gridalloc,excess;


    /* Determine length of triclinic box vectors */
    for(d=0; d<DIM; d++)
    {
        box_size[d] = 0;
        for(i=0;i<DIM;i++)
            box_size[d] += box[d][i]*box[d][i];
        box_size[d] = sqrt(box_size[d]);
    }

    gridalloc = 25;
    snew(gridall, gridalloc);

    fprintf(stdout, "Possible PME grid settings (apart from input file settings):\n");
    fprintf(stdout, "  nkx  nky  nkz  max spacing  scaling factor\n");

    /* eps should provide a fine enough spacing not to miss any grid */
    if (fmax != fmin)
        eps = (fmax-fmin)/(100.0*(*ntprs - 1));
    else
        eps = 1.0/max( (*griduse)[0].nkz, max( (*griduse)[0].nkx, (*griduse)[0].nky ) );

    for (req_fac = fmin; act_fac < fmax; req_fac += eps)
    {
        nkx=0;
        nky=0;
        nkz=0;
        calc_grid(NULL,box,fs*req_fac,&nkx,&nky,&nkz);
        act_fs = max(box_size[XX]/nkx,max(box_size[YY]/nky,box_size[ZZ]/nkz));
        act_fac = act_fs/fs;
        if (    ! ( nkx==nkx_old           && nky==nky_old           && nkz==nkz_old           )    /* Exclude if grid is already in list */
             && ! ( nkx==(*griduse)[0].nkx && nky==(*griduse)[0].nky && nkz==(*griduse)[0].nkz ) )  /* Exclude input file grid */
        {
            /* We found a new grid that will do */
            nkx_old = nkx;
            nky_old = nky;
            nkz_old = nkz;
            gridall[ngridall].nkx = nkx;
            gridall[ngridall].nky = nky;
            gridall[ngridall].nkz = nkz;
            gridall[ngridall].fac = act_fac;
            gridall[ngridall].fs  = act_fs;
            fprintf(stdout, "%5d%5d%5d %12f %12f\n",nkx,nky,nkz,act_fs,act_fac);
            ngridall++;
            if (ngridall >= gridalloc)
            {
                gridalloc += 25;
                srenew(gridall, gridalloc);
            }
        }
    }

    /* excess is the number of grids we have to get rid of */
    excess = ngridall+1 - *ntprs;

    /* If we found less grids than tpr files were requested, simply test all
     * the grid there are ... */
    if (excess < 0)
    {
        fprintf(stdout, "NOTE: You requested %d tpr files, but apart from the input file grid only the\n"
                        "      above listed %d suitable PME grids were found. Will test all suitable settings.\n",
                        *ntprs, ngridall);
        *ntprs = ngridall+1;
    }
    else
    {
        if (2 == *ntprs)
        {
            /* We can only keep the input tpr file plus one extra tpr file.
             * We make that choice depending on the values of -upfac and -downfac */
            if (fmin < 1.0)
            {
                /* Keep the one with the -downfac as scaling factor. This is already
                 * stored in gridall[0] */
                ;
            }
            else
            {
                /* Keep the one with -upfac as scaling factor */
                copy_grid(&(gridall[ngridall-1]), &(gridall[0]));
            }

        }
        else
        {
            /* From the grid list throw away the unnecessary ones (keep the extremes) */
            for (i = 0; i < excess; i++)
            {
                /* Get the index of the least necessary grid from the list ... */
                jmin = get_throwaway_index(gridall, ngridall);
                /* ... and remove the grid from the list */
                remove_from_list(gridall, &ngridall, jmin);
            }
        }
    }

    /* The remaining list contains the grids we want to test */
    for (i=1; i < *ntprs; i++)
        copy_grid(&(gridall[i-1]), &((*griduse)[i]));

    sfree(gridall);
}


#define EPME_SWITCHED(e) ((e) == eelPMESWITCH || (e) == eelPMEUSERSWITCH)

/* Make additional TPR files with more computational load for the
 * direct space processors: */
static void make_benchmark_tprs(
        const char *fn_sim_tpr,     /* READ : User-provided tpr file                 */
        char *fn_bench_tprs[],      /* WRITE: Names of benchmark tpr files           */
        gmx_large_int_t benchsteps, /* Number of time steps for benchmark runs       */
        gmx_large_int_t statesteps, /* Step counter in checkpoint file               */
        real upfac,                 /* Scale rcoulomb inbetween downfac and upfac    */
        real downfac,
        int *ntprs,                 /* No. of TPRs to write, each with a different
                                       rcoulomb and fourierspacing. If not enough
                                       grids are found, ntprs is reduced accordingly */
        real fourierspacing,        /* Basic fourierspacing from tpr input file      */
        t_inputinfo *info,          /* Contains information about mdp file options   */
        FILE *fp)                   /* Write the output here                         */
{
    int          i,j,d;
    t_inputrec   *ir;
    t_state      state;
    gmx_mtop_t   mtop;
    real         nlist_buffer; /* Thickness of the buffer regions for PME-switch potentials: */
    char         buf[200];
    rvec         box_size;
    gmx_bool         bNote = FALSE;
    t_pmegrid    *pmegrid=NULL; /* Grid settings for the PME grids to test    */
    

    sprintf(buf, "Making benchmark tpr file%s with %s time step%s",
            *ntprs > 1? "s":"", gmx_large_int_pfmt, benchsteps>1?"s":"");
    fprintf(stdout, buf, benchsteps);
    if (statesteps > 0)
    {
        sprintf(buf, " (adding %s steps from checkpoint file)", gmx_large_int_pfmt);
        fprintf(stdout, buf, statesteps);
        benchsteps += statesteps;
    }
    fprintf(stdout, ".\n");
        
    
    snew(ir,1);
    read_tpx_state(fn_sim_tpr,ir,&state,NULL,&mtop);

    /* Check if some kind of PME was chosen */
    if (EEL_PME(ir->coulombtype) == FALSE)
        gmx_fatal(FARGS, "Can only do optimizations for simulations with %s electrostatics.",
                EELTYPE(eelPME));
    
    /* Check if rcoulomb == rlist, which is necessary for plain PME. */
    if (  (eelPME == ir->coulombtype) && !(ir->rcoulomb == ir->rlist) )
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

    /* Reduce the number of steps for the benchmarks */
    info->orig_sim_steps = ir->nsteps;
    ir->nsteps           = benchsteps;
    /* We must not use init_step from the input tpr file for the benchmarks */
    info->orig_init_step = ir->init_step;
    ir->init_step        = 0;

    /* For PME-switch potentials, keep the radial distance of the buffer region */
    nlist_buffer   = ir->rlist - ir->rcoulomb;

    /* Determine length of triclinic box vectors */
    for(d=0; d<DIM; d++)
    {
        box_size[d] = 0;
        for(i=0;i<DIM;i++)
            box_size[d] += state.box[d][i]*state.box[d][i];
        box_size[d] = sqrt(box_size[d]);
    }

    /* Reconstruct fourierspacing per dimension from the number of grid points and box size */
    info->fsx[0] = box_size[XX]/ir->nkx;
    info->fsy[0] = box_size[YY]/ir->nky;
    info->fsz[0] = box_size[ZZ]/ir->nkz;

    /* Put the input grid as first entry into the grid list */
    snew(pmegrid, *ntprs);
    pmegrid[0].fac = 1.00;
    pmegrid[0].nkx = ir->nkx;
    pmegrid[0].nky = ir->nky;
    pmegrid[0].nkz = ir->nkz;
    pmegrid[0].fs  = max(box_size[ZZ]/ir->nkz, max(box_size[XX]/ir->nkx, box_size[YY]/ir->nky));

    /* If no value for the fourierspacing was provided on the command line, we
     * use the reconstruction from the tpr file */
    if (fourierspacing <= 0)
    {
        fourierspacing = pmegrid[0].fs;
    }

    fprintf(stdout, "Calculating PME grid points on the basis of a fourierspacing of %f nm\n", fourierspacing);

    /* Print information about settings of which some are potentially modified: */
    fprintf(fp, "   Coulomb type         : %s\n", EELTYPE(ir->coulombtype));
    fprintf(fp, "   Grid spacing x y z   : %f %f %f\n",
            box_size[XX]/ir->nkx, box_size[YY]/ir->nky, box_size[ZZ]/ir->nkz);
    fprintf(fp, "   Van der Waals type   : %s\n", EVDWTYPE(ir->vdwtype));
    if (EVDW_SWITCHED(ir->vdwtype))
        fprintf(fp, "   rvdw_switch          : %f nm\n", ir->rvdw_switch);
    if (EPME_SWITCHED(ir->coulombtype))
        fprintf(fp, "   rlist                : %f nm\n", ir->rlist);
    if (ir->rlistlong != max_cutoff(ir->rvdw,ir->rcoulomb))
        fprintf(fp, "   rlistlong            : %f nm\n", ir->rlistlong);

    /* Print a descriptive line about the tpr settings tested */
    fprintf(fp, "\nWill try these real/reciprocal workload settings:\n");
    fprintf(fp, " No.   scaling  rcoulomb");
    fprintf(fp, "  nkx  nky  nkz");
    fprintf(fp, "   spacing");
    if (evdwCUT == ir->vdwtype)
        fprintf(fp, "      rvdw");
    if (EPME_SWITCHED(ir->coulombtype))
        fprintf(fp, "     rlist");
    if ( ir->rlistlong != max_cutoff(ir->rlist,max_cutoff(ir->rvdw,ir->rcoulomb)) )
        fprintf(fp, " rlistlong");
    fprintf(fp, "  tpr file\n");


    /* Get useful PME grids if requested, the actual number of entries is
     * returned in npmegrid */
    if (*ntprs > 1)
        make_grid_list(downfac, upfac, ntprs, state.box, &pmegrid, fourierspacing);
    
    /* Loop to create the requested number of tpr input files */
    for (j = 0; j < *ntprs; j++)
    {
        /* The first one is the provided tpr file, just need to modify the number
         * of steps, so skip the following block */
        if (j > 0)
        {
            /* Scale the Coulomb radius */
            ir->rcoulomb = info->rcoulomb[0]*pmegrid[j].fac;

            /* Adjust other radii since various conditions neet to be fulfilled */
            if (eelPME == ir->coulombtype)
            {
                /* plain PME, rcoulomb must be equal to rlist */
                ir->rlist = ir->rcoulomb;
            }
            else
            {
                /* rlist must be >= rcoulomb, we keep the size of the buffer region */
                ir->rlist = ir->rcoulomb + nlist_buffer;
            }

            if (evdwCUT == ir->vdwtype)
            {
                /* For vdw cutoff, rvdw >= rlist */
                ir->rvdw = max(info->rvdw[0], ir->rlist);
            }

            ir->rlistlong = max_cutoff(ir->rlist,max_cutoff(ir->rvdw,ir->rcoulomb));

            /* Copy the optimal grid dimensions to ir */
            ir->nkx = pmegrid[j].nkx;
            ir->nky = pmegrid[j].nky;
            ir->nkz = pmegrid[j].nkz;

        } /* end of "if (j != 0)" */

        /* for j==0: Save the original settings
         * for j >0: Save modified radii and fourier grids */
        info->rcoulomb[j]  = ir->rcoulomb;
        info->rvdw[j]      = ir->rvdw;
        info->nkx[j]       = ir->nkx;
        info->nky[j]       = ir->nky;
        info->nkz[j]       = ir->nkz;
        info->rlist[j]     = ir->rlist;
        info->rlistlong[j] = ir->rlistlong;
        info->fsx[j]       = pmegrid[j].fs;
        info->fsy[j]       = pmegrid[j].fs;
        info->fsz[j]       = pmegrid[j].fs;

        /* Write the benchmark tpr file */
        strncpy(fn_bench_tprs[j],fn_sim_tpr,strlen(fn_sim_tpr)-strlen(".tpr"));
        sprintf(buf, "_bench%.2d.tpr", j);
        strcat(fn_bench_tprs[j], buf);
        fprintf(stdout,"Writing benchmark tpr %s with nsteps=", fn_bench_tprs[j]);
        fprintf(stdout, gmx_large_int_pfmt, ir->nsteps);
        if (j > 0)
            fprintf(stdout,", scaling factor %f\n", pmegrid[j].fac);
        else
            fprintf(stdout,", unmodified settings\n");

        write_tpx_state(fn_bench_tprs[j],ir,&state,&mtop);
        
        /* Write information about modified tpr settings to log file */
        if (j==0)
            fprintf(fp, "%4d%10s%10f", j, "-input-", ir->rcoulomb);
        else
            fprintf(fp, "%4d%10f%10f", j, pmegrid[j].fac, ir->rcoulomb);
        fprintf(fp, "%5d%5d%5d", ir->nkx, ir->nky, ir->nkz);
        fprintf(fp, " %9f ", info->fsx[j]);
        if (evdwCUT == ir->vdwtype)
            fprintf(fp, "%10f", ir->rvdw);
        if (EPME_SWITCHED(ir->coulombtype))
            fprintf(fp, "%10f", ir->rlist);
        if ( info->rlistlong[0] != max_cutoff(info->rlist[0],max_cutoff(info->rvdw[0],info->rcoulomb[0])) )
            fprintf(fp, "%10f", ir->rlistlong);
        fprintf(fp, "  %-14s\n",fn_bench_tprs[j]);

        /* Make it clear to the user that some additional settings were modified */
        if (   !is_equal(ir->rvdw     , info->rvdw[0])
            || !is_equal(ir->rlistlong, info->rlistlong[0]) )
        {
            bNote = TRUE;
        }
    }
    if (bNote)
        fprintf(fp, "\nNote that in addition to rcoulomb and the fourier grid\n"
                    "also other input settings were changed (see table above).\n"
                    "Please check if the modified settings are appropriate.\n");
    fflush(stdout);
    fflush(fp);
    sfree(ir);
}


/* Rename the files we want to keep to some meaningful filename and
 * delete the rest */
static void cleanup(const t_filenm *fnm, int nfile, int k, int nnodes, 
                    int nPMEnodes, int nr, gmx_bool bKeepStderr)
{
    char numstring[STRLEN];
    char newfilename[STRLEN];
    const char *fn=NULL;
    int i;
    const char *opt;


    fprintf(stdout, "Cleaning up, deleting benchmark temp files ...\n");

    for (i=0; i<nfile; i++)
    {
        opt = (char *)fnm[i].opt;
        if ( strcmp(opt, "-p") == 0 )
        {
            /* do nothing; keep this file */
            ;
        }
        else if (strcmp(opt, "-bg") == 0)
        {
            /* Give the log file a nice name so one can later see which parameters were used */
            numstring[0] = '\0';
            if (nr > 0)
                sprintf(numstring, "_%d", nr);
            sprintf(newfilename, "%s_no%d_np%d_npme%d%s", opt2fn("-bg",nfile,fnm), k, nnodes, nPMEnodes, numstring);
            if (gmx_fexist(opt2fn("-bg",nfile,fnm)))
            {
                fprintf(stdout, "renaming log file to %s\n", newfilename);
                make_backup(newfilename);
                rename(opt2fn("-bg",nfile,fnm), newfilename);
            }
        }
        else if (strcmp(opt, "-err") == 0)
        {
            /* This file contains the output of stderr. We want to keep it in
             * cases where there have been problems. */
            fn = opt2fn(opt, nfile, fnm);
            numstring[0] = '\0';
            if (nr > 0)
                sprintf(numstring, "_%d", nr);
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
        else if ( (0 == strncmp(opt, "-b", 2)) && (opt2bSet(opt,nfile,fnm) || !is_optional(&fnm[i])) )
        {
            fn = opt2fn(opt, nfile, fnm);
            if (gmx_fexist(fn))
            {
                fprintf(stdout, "Deleting %s\n", fn);
                remove(fn);
            }
        }
    }
}


/* Returns the largest common factor of n1 and n2 */
static int largest_common_factor(int n1, int n2)
{
    int factor, nmax;

    nmax = min(n1, n2);
    for (factor=nmax; factor > 0; factor--)
    {
        if ( 0==(n1 % factor) && 0==(n2 % factor) )
        {
            return(factor);
        }
    }
    return 0; /* one for the compiler */
}

enum {eNpmeAuto, eNpmeAll, eNpmeReduced, eNpmeSubset, eNpmeNr};

/* Create a list of numbers of PME nodes to test */
static void make_npme_list(
        const char *npmevalues_opt,  /* Make a complete list with all
                           * possibilities or a short list that keeps only
                           * reasonable numbers of PME nodes                  */
        int *nentries,    /* Number of entries we put in the nPMEnodes list   */
        int *nPMEnodes[], /* Each entry contains the value for -npme          */
        int nnodes,       /* Total number of nodes to do the tests on         */
        int minPMEnodes,  /* Minimum number of PME nodes                      */
        int maxPMEnodes)  /* Maximum number of PME nodes                      */
{
    int i,npme,npp;
    int min_factor=1;     /* We request that npp and npme have this minimal
                           * largest common factor (depends on npp)           */
    int nlistmax;         /* Max. list size                                   */
    int nlist;            /* Actual number of entries in list                 */
    int eNPME=0;


    /* Do we need to check all possible values for -npme or is a reduced list enough? */
    if ( 0 == strcmp(npmevalues_opt, "all") )
    {
        eNPME = eNpmeAll;
    }
    else if ( 0 == strcmp(npmevalues_opt, "subset") )
    {
        eNPME = eNpmeSubset;
    }
    else if ( 0 == strcmp(npmevalues_opt, "auto") )
    {
        if (nnodes <= 64)
            eNPME = eNpmeAll;
        else if (nnodes < 128)
            eNPME = eNpmeReduced;
        else
            eNPME = eNpmeSubset;
    }
    else
    {
        gmx_fatal(FARGS, "Unknown option for -npme in make_npme_list");
    }

    /* Calculate how many entries we could possibly have (in case of -npme all) */
    if (nnodes > 2)
    {
        nlistmax = maxPMEnodes - minPMEnodes + 3;
        if (0 == minPMEnodes)
            nlistmax--;
    }
    else
        nlistmax = 1;

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
                min_factor = (int) pow(npp, 1.0/3.0);
                break;
            default:
                gmx_fatal(FARGS, "Unknown option for eNPME in make_npme_list");
                break;
        }
        if (largest_common_factor(npp, npme) >= min_factor)
        {
            (*nPMEnodes)[nlist] = npme;
            nlist++;
        }
    }
    /* We always test 0 PME nodes and the automatic number */
    *nentries = nlist + 2;
    (*nPMEnodes)[nlist  ] =  0;
    (*nPMEnodes)[nlist+1] = -1;

    fprintf(stderr, "Will try the following %d different values for -npme:\n", *nentries);
    for (i=0; i<*nentries-1; i++)
        fprintf(stderr, "%d, ", (*nPMEnodes)[i]);
    fprintf(stderr, "and %d (auto).\n", (*nPMEnodes)[*nentries-1]);
}


/* Allocate memory to store the performance data */
static void init_perfdata(t_perf *perfdata[], int ntprs, int datasets, int repeats)
{
    int i, j, k;


    for (k=0; k<ntprs; k++)
    {
        snew(perfdata[k], datasets);
        for (i=0; i<datasets; i++)
        {
            for (j=0; j<repeats; j++)
            {
                snew(perfdata[k][i].Gcycles   , repeats);
                snew(perfdata[k][i].ns_per_day, repeats);
                snew(perfdata[k][i].PME_f_load, repeats);
            }
        }
    }
}


static void do_the_tests(
        FILE *fp,                   /* General g_tune_pme output file         */
        char **tpr_names,           /* Filenames of the input files to test   */
        int maxPMEnodes,            /* Max fraction of nodes to use for PME   */
        int minPMEnodes,            /* Min fraction of nodes to use for PME   */
        int npme_fixed,             /* If >= -1, test fixed number of PME
                                     * nodes only                             */
        const char *npmevalues_opt, /* Which -npme values should be tested    */
        t_perf **perfdata,          /* Here the performace data is stored     */
        int *pmeentries,            /* Entries in the nPMEnodes list          */
        int repeats,                /* Repeat each test this often            */
        int nnodes,                 /* Total number of nodes = nPP + nPME     */
        int nr_tprs,                /* Total number of tpr files to test      */
        gmx_bool bThreads,          /* Threads or MPI?                        */
        char *cmd_mpirun,           /* mpirun command string                  */
        char *cmd_np,               /* "-np", "-n", whatever mpirun needs     */
        char *cmd_mdrun,            /* mdrun command string                   */
        char *cmd_args_bench,       /* arguments for mdrun in a string        */
        const t_filenm *fnm,        /* List of filenames from command line    */
        int nfile,                  /* Number of files specified on the cmdl. */
        int sim_part,               /* For checkpointing                      */
        int presteps,               /* DLB equilibration steps, is checked    */
        gmx_large_int_t cpt_steps)  /* Time step counter in the checkpoint    */
{
    int     i,nr,k,ret,count=0,totaltests;
    int     *nPMEnodes=NULL;
    t_perf  *pd=NULL;
    int     cmdline_length;
    char    *command, *cmd_stub;
    char    buf[STRLEN];
    gmx_bool    bResetProblem=FALSE;


    /* This string array corresponds to the eParselog enum type at the start
     * of this file */
    const char* ParseLog[] = {"OK.",
                              "Logfile not found!",
                              "No timings, logfile truncated?",
                              "Run was terminated.",
                              "Counters were not reset properly.",
                              "No DD grid found for these settings.",
                              "TPX version conflict!",
                              "mdrun was not started in parallel!",
                              "A fatal error occured!" };
    char    str_PME_f_load[13];


    /* Allocate space for the mdrun command line. 100 extra characters should 
       be more than enough for the -npme etcetera arguments */
    cmdline_length =  strlen(cmd_mpirun) 
                    + strlen(cmd_np)
                    + strlen(cmd_mdrun) 
                    + strlen(cmd_args_bench) 
                    + strlen(tpr_names[0]) + 100;
    snew(command , cmdline_length);
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
        fprintf(stderr, "Will use a fixed number of %d PME-only nodes.\n", nPMEnodes[0]);
    }

    if (0 == repeats)
    {
        fprintf(fp, "\nNo benchmarks done since number of repeats (-r) is 0.\n");
        fclose(fp);
        finalize(opt2fn("-p", nfile, fnm));
        exit(0);
    }

    /* Allocate one dataset for each tpr input file: */
    init_perfdata(perfdata, nr_tprs, *pmeentries, repeats);

    /*****************************************/
    /* Main loop over all tpr files to test: */
    /*****************************************/
    totaltests = nr_tprs*(*pmeentries)*repeats;
    for (k=0; k<nr_tprs;k++)
    {
        fprintf(fp, "\nIndividual timings for input file %d (%s):\n", k, tpr_names[k]);
        fprintf(fp, "PME nodes      Gcycles       ns/day        PME/f    Remark\n");
        /* Loop over various numbers of PME nodes: */
        for (i = 0; i < *pmeentries; i++)
        {
            pd = &perfdata[k][i];

            /* Loop over the repeats for each scenario: */
            for (nr = 0; nr < repeats; nr++)
            {
                pd->nPMEnodes = nPMEnodes[i];
                
                /* Add -npme and -s to the command line and save it. Note that
                 * the -passall (if set) options requires cmd_args_bench to be
                 * at the end of the command line string */
                snew(pd->mdrun_cmd_line, cmdline_length);
                sprintf(pd->mdrun_cmd_line, "%s-npme %d -s %s %s",
                        cmd_stub, pd->nPMEnodes, tpr_names[k], cmd_args_bench);

                /* Do a benchmark simulation: */
                if (repeats > 1)
                    sprintf(buf, ", pass %d/%d", nr+1, repeats);
                else
                    buf[0]='\0';
                fprintf(stdout, "\n=== Progress %2.0f%%, tpr %d/%d, run %d/%d%s:\n",
                        (100.0*count)/totaltests,
                        k+1, nr_tprs, i+1, *pmeentries, buf);
                make_backup(opt2fn("-err",nfile,fnm));
                sprintf(command, "%s 1> /dev/null 2>%s", pd->mdrun_cmd_line, opt2fn("-err",nfile,fnm));
                fprintf(stdout, "%s\n", pd->mdrun_cmd_line);
                gmx_system_call(command);

                /* Collect the performance data from the log file; also check stderr
                 * for fatal errors */
                ret = parse_logfile(opt2fn("-bg",nfile,fnm), opt2fn("-err",nfile,fnm),
                        pd, nr, presteps, cpt_steps, nnodes);
                if ((presteps > 0) && (ret == eParselogResetProblem))
                    bResetProblem = TRUE;

                if (-1 == pd->nPMEnodes)
                    sprintf(buf, "(%3d)", pd->guessPME);
                else
                    sprintf(buf, "     ");

                /* Nicer output */
                if (pd->PME_f_load[nr] > 0.0)
                    sprintf(str_PME_f_load, "%12.3f", pd->PME_f_load[nr]);
                else
                    sprintf(str_PME_f_load, "%s", "         -  ");
                
                /* Write the data we got to disk */
                fprintf(fp, "%4d%s %12.3f %12.3f %s    %s", pd->nPMEnodes,
                        buf, pd->Gcycles[nr], pd->ns_per_day[nr], str_PME_f_load, ParseLog[ret]);
                if (! (ret==eParselogOK || ret==eParselogNoDDGrid || ret==eParselogNotFound) )
                    fprintf(fp, " Check %s file for problems.", ret==eParselogFatal? "err":"log");
                fprintf(fp, "\n");
                fflush(fp);
                count++;

                /* Do some cleaning up and delete the files we do not need any more */
                cleanup(fnm, nfile, k, nnodes, pd->nPMEnodes, nr, ret==eParselogFatal);

                /* If the first run with this number of processors already failed, do not try again: */
                if (pd->Gcycles[0] <= 0.0 && repeats > 1)
                {
                    fprintf(stdout, "Skipping remaining passes of unsuccessful setting, see log file for details.\n");
                    count += repeats-(nr+1);
                    break;
                }
            } /* end of repeats loop */
        } /* end of -npme loop */
    } /* end of tpr file loop */
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
        int nnodes, 
        int repeats, 
        int *ntprs, 
        real *upfac,
        real *downfac,
        real maxPMEfraction,
        real minPMEfraction,
        int  npme_fixed,
        real fourierspacing,
        gmx_large_int_t bench_nsteps,
        const t_filenm *fnm,
        int nfile,
        int sim_part,
        int presteps,
        int npargs,
        t_pargs *pa)
{
    /* Make sure the input file exists */
    if (!gmx_fexist(opt2fn("-s",nfile,fnm)))
        gmx_fatal(FARGS, "File %s not found.", opt2fn("-s",nfile,fnm));
    
    /* Make sure that the checkpoint file is not overwritten during benchmarking */
    if ( (0 == strcmp(opt2fn("-cpi",nfile,fnm), opt2fn("-bcpo",nfile,fnm)) ) && (sim_part > 1) )
        gmx_fatal(FARGS, "Checkpoint input (-cpi) and benchmark checkpoint output (-bcpo) files must not be identical.\n"
                         "The checkpoint input file must not be overwritten during the benchmarks.\n");

    /* Make sure that repeats is >= 0 (if == 0, only write tpr files) */
    if (repeats < 0)
        gmx_fatal(FARGS, "Number of repeats < 0!");

    /* Check number of nodes */
    if (nnodes < 1)
        gmx_fatal(FARGS, "Number of nodes/threads must be a positive integer.");

    /* Automatically choose -ntpr if not set */
    if (*ntprs < 1)
    {
        if (nnodes < 16)
            *ntprs = 1;
        else 
            *ntprs = 3;
        fprintf(stderr, "Will test %d tpr file%s.\n", *ntprs, *ntprs==1?"":"s");
    }
    else
    {
        if (1 == *ntprs)
            fprintf(stderr, "Note: Choose ntpr>1 to shift PME load between real and reciprocal space.\n");
    }
    
    if ( is_equal(*downfac,*upfac) && (*ntprs > 2) && !(fourierspacing > 0.0))
    {
        fprintf(stderr, "WARNING: Resetting -ntpr to 2 since both scaling factors are the same.\n"
                        "Please choose upfac unequal to downfac to test various PME grid settings\n");
        *ntprs = 2;
    }

    /* Check whether max and min fraction are within required values */
    if (maxPMEfraction > 0.5 || maxPMEfraction < 0)
        gmx_fatal(FARGS, "-max must be between 0 and 0.5");
    if (minPMEfraction > 0.5 || minPMEfraction < 0)
        gmx_fatal(FARGS, "-min must be between 0 and 0.5");
    if (maxPMEfraction < minPMEfraction)
        gmx_fatal(FARGS, "-max must be larger or equal to -min");
    
    /* Check whether the number of steps - if it was set - has a reasonable value */
    if (bench_nsteps < 0)
        gmx_fatal(FARGS, "Number of steps must be positive.");

    if (bench_nsteps > 10000 || bench_nsteps < 100)
    {
        fprintf(stderr, "WARNING: steps=");
        fprintf(stderr, gmx_large_int_pfmt, bench_nsteps);
        fprintf(stderr, ". Are you sure you want to perform so %s steps for each benchmark?\n", (bench_nsteps < 100)? "few" : "many");
    }
    
    if (presteps < 0)
    {
        gmx_fatal(FARGS, "Cannot have a negative number of presteps.\n");
    }
    
    if (*upfac <= 0.0 || *downfac <= 0.0 || *downfac > *upfac)
        gmx_fatal(FARGS, "Both scaling factors must be larger than zero and upper\n"
                         "scaling limit (%f) must be larger than lower limit (%f).",
                         *upfac, *downfac);

    if (*downfac < 0.75 || *upfac > 1.25)
        fprintf(stderr, "WARNING: Applying extreme scaling factor. I hope you know what you are doing.\n");
    
    if (fourierspacing < 0)
        gmx_fatal(FARGS, "Please choose a positive value for fourierspacing.");

    /* Make shure that the scaling factor options are compatible with the number of tprs */
    if ( (*ntprs < 3) && ( opt2parg_bSet("-upfac",npargs,pa) || opt2parg_bSet("-downfac",npargs,pa) ) )
    {
        if (opt2parg_bSet("-upfac",npargs,pa) && opt2parg_bSet("-downfac",npargs,pa) && !is_equal(*upfac,*downfac))
        {
            fprintf(stderr, "NOTE: Specify -ntpr > 2 for both scaling factors to take effect.\n"
                             "(upfac=%f, downfac=%f)\n", *upfac, *downfac);
        }
        if (opt2parg_bSet("-upfac",npargs,pa))
            *downfac = *upfac;
        if (opt2parg_bSet("-downfac",npargs,pa))
            *upfac = *downfac;
    }
    if ( (2 == *ntprs) && (opt2parg_bSet("-downfac",npargs,pa)) && !is_equal(*downfac, 1.0))
    {
        *upfac = 1.0;
    }

    /* If a fixed number of PME nodes is set we do rcoulomb and PME gird tuning
     * only. We need to check whether the requested number of PME-only nodes
     * makes sense. */
    if (npme_fixed > -1)
    {
        /* No more than 50% of all nodes can be assigned as PME-only nodes. */
        if (2*npme_fixed > nnodes)
        {
            gmx_fatal(FARGS, "Cannot have more than %d PME-only nodes for a total of %d nodes (you chose %d).\n",
                             nnodes/2, nnodes, npme_fixed);
        }
        if ((npme_fixed > 0) && (5*npme_fixed < nnodes))
        {
            fprintf(stderr, "WARNING: Only %g percent of the nodes are assigned as PME-only nodes.\n",
                             100.0*((real)npme_fixed / (real)nnodes));
        }
        if (opt2parg_bSet("-min",npargs,pa) || opt2parg_bSet("-max",npargs,pa))
        {
            fprintf(stderr, "NOTE: The -min, -max, and -npme options have no effect when a\n"
                            "      fixed number of PME-only nodes is requested with -fix.\n");
        }
    }
}


/* Returns TRUE when "opt" is a switch for g_tune_pme itself */
static gmx_bool is_main_switch(char *opt)
{
    if ( (0 == strcmp(opt,"-s"        ))
      || (0 == strcmp(opt,"-p"        ))
      || (0 == strcmp(opt,"-launch"   ))
      || (0 == strcmp(opt,"-r"        ))
      || (0 == strcmp(opt,"-ntpr"     ))
      || (0 == strcmp(opt,"-max"      ))
      || (0 == strcmp(opt,"-min"      ))
      || (0 == strcmp(opt,"-upfac"    ))
      || (0 == strcmp(opt,"-downfac"  ))
      || (0 == strcmp(opt,"-fix"      ))
      || (0 == strcmp(opt,"-four"     ))
      || (0 == strcmp(opt,"-steps"    ))
      || (0 == strcmp(opt,"-simsteps" ))
      || (0 == strcmp(opt,"-resetstep"))
      || (0 == strcmp(opt,"-so"       ))
      || (0 == strcmp(opt,"-npstring" ))
      || (0 == strcmp(opt,"-npme"     ))
      || (0 == strcmp(opt,"-err"      ))
      || (0 == strcmp(opt,"-passall"  )) )
    return TRUE;
    
    return FALSE;
}


/* Returns TRUE when "opt" is needed at launch time */
static gmx_bool is_launch_option(char *opt, gmx_bool bSet)
{
    if (bSet)
        return TRUE;
    else    
        return FALSE;
}


/* Returns TRUE when "opt" is needed at launch time */
static gmx_bool is_launch_file(char *opt, gmx_bool bSet)
{
    /* We need all options that were set on the command line 
     * and that do not start with -b */
    if (0 == strncmp(opt,"-b", 2))
        return FALSE;

    if (bSet)
        return TRUE;
    else
        return FALSE;
}


/* Returns TRUE when "opt" gives an option needed for the benchmarks runs */
static gmx_bool is_bench_option(char *opt, gmx_bool bSet)
{
    /* If option is set, we might need it for the benchmarks.
     * This includes -cpi */
    if (bSet)
    {
        if ( (0 == strcmp(opt, "-append"   ))
          || (0 == strcmp(opt, "-maxh"     ))
          || (0 == strcmp(opt, "-deffnm"   ))
          || (0 == strcmp(opt, "-resethway"))
          || (0 == strcmp(opt, "-cpnum"    )) )
            return FALSE;
        else
            return TRUE;
    }
    else
        return FALSE;
}


/* Returns TRUE when "opt" defines a file which is needed for the benchmarks runs */
static gmx_bool is_bench_file(char *opt, gmx_bool bSet, gmx_bool bOptional, gmx_bool bIsOutput)
{
    /* All options starting with "-b" are for _b_enchmark files exclusively */
    if (0 == strncmp(opt,"-b", 2))
    { 
        if (!bOptional || bSet)
            return TRUE;
        else
            return FALSE;
    }
    else
    {
        if (bIsOutput)
            return FALSE;
        else
            if (bSet) /* These are additional input files like -cpi -ei */
                return TRUE;
            else 
                return FALSE;
    }
}


/* Adds 'buf' to 'str' */
static void add_to_string(char **str, char *buf)
{
    int len;
    
    
    len = strlen(*str) + strlen(buf) + 1;
    srenew(*str, len);
    strcat(*str, buf);
}


/* Create the command line for the benchmark as well as for the real run */
static void create_command_line_snippets(
        gmx_bool     bThreads,
        int      presteps,
        int      nfile,
        t_filenm fnm[],
        int      npargs,
        t_pargs  *pa,
        const char *procstring,      /* How to pass the number of processors to $MPIRUN */
        char     *cmd_np[],          /* Actual command line snippet, e.g. '-np <N>' */
        char     *cmd_args_bench[],  /* command line arguments for benchmark runs */
        char     *cmd_args_launch[], /* command line arguments for simulation run */
        char     extra_args[])       /* Add this to the end of the command line */
{
    int        i;
    char       *opt;
    const char *name;
    char       *np_or_nt;
#define BUFLENGTH 255
    char       buf[BUFLENGTH];
    char       strbuf[BUFLENGTH];
    char       strbuf2[BUFLENGTH];
    

    if (bThreads)
        np_or_nt=strdup("-nt");
    else
        np_or_nt=strdup("-np");
    
    /* strlen needs at least '\0' as a string: */
    snew(*cmd_args_bench ,1);
    snew(*cmd_args_launch,1);
    *cmd_args_launch[0]='\0';
    *cmd_args_bench[0] ='\0';
    
        
    /*******************************************/
    /* 1. Process other command line arguments */
    /*******************************************/
    for (i=0; i<npargs; i++)
    {
        /* What command line switch are we currently processing: */
        opt = (char *)pa[i].option;
        /* Skip options not meant for mdrun */        
        if (!is_main_switch(opt))
        {
            /* Boolean arguments need to be generated in the -[no]argname format */
            if (pa[i].type == etBOOL)
            {
                sprintf(strbuf,"-%s%s ",*pa[i].u.b ? "" : "no",opt+1);
            }
            else
            {
                /* Print it to a string buffer, strip away trailing whitespaces that pa_val also returns: */
                sprintf(strbuf2, "%s", pa_val(&pa[i],buf,BUFLENGTH));
                rtrim(strbuf2);
                sprintf(strbuf, "%s %s ", opt, strbuf2);
            }

            /* We need the -np (or -nt) switch in a separate buffer - whether or not it was set! */
            if (0 == strcmp(opt,np_or_nt))
            {
                if (strcmp(procstring, "none")==0 && !bThreads)
                {
                    /* Omit -np <N> entirely */
                    snew(*cmd_np, 2);
                    sprintf(*cmd_np, " ");
                }
                else
                {
                    /* This is the normal case with -np <N> */
                    snew(*cmd_np, strlen(procstring)+strlen(strbuf2)+4);
                    sprintf(*cmd_np, " %s %s ", bThreads? "-nt" : procstring, strbuf2);
                }
            }
            else 
            {
                if (is_bench_option(opt,pa[i].bSet))
                    add_to_string(cmd_args_bench, strbuf);

                if (is_launch_option(opt,pa[i].bSet))
                    add_to_string(cmd_args_launch, strbuf);
            }
        }
    }
    if (presteps > 0)
    {
        /* Add equilibration steps to benchmark options */
        sprintf(strbuf, "-resetstep %d ", presteps);
        add_to_string(cmd_args_bench, strbuf);
    }
    
    /********************/
    /* 2. Process files */
    /********************/
    for (i=0; i<nfile; i++)
    {
        opt  = (char *)fnm[i].opt;
        name = opt2fn(opt,nfile,fnm);
        
        /* Strbuf contains the options, now let's sort out where we need that */
        sprintf(strbuf, "%s %s ", opt, name);
        
        /* Skip options not meant for mdrun */        
        if (!is_main_switch(opt))
        {
            
            if ( is_bench_file(opt, opt2bSet(opt,nfile,fnm), is_optional(&fnm[i]), is_output(&fnm[i])) )
            {
                /* All options starting with -b* need the 'b' removed,
                 * therefore overwrite strbuf */
                if (0 == strncmp(opt, "-b", 2))     
                    sprintf(strbuf, "-%s %s ", &opt[2], name);
                
                add_to_string(cmd_args_bench, strbuf);
            }
            
            if ( is_launch_file(opt,opt2bSet(opt,nfile,fnm)) )
                add_to_string(cmd_args_launch, strbuf);
        }
    }

    add_to_string(cmd_args_bench , extra_args);
    add_to_string(cmd_args_launch, extra_args);
#undef BUFLENGTH
}


/* Set option opt */
static void setopt(const char *opt,int nfile,t_filenm fnm[])
{
  int i;
  
  for(i=0; (i<nfile); i++)
    if (strcmp(opt,fnm[i].opt)==0)
      fnm[i].flag |= ffSET;
}


/* This routine checks for output files that get triggered by a tpr option.
 * These output files are marked as set to allow for proper cleanup after 
 * each tuning run. */
static void get_tpr_outfiles(int nfile, t_filenm fnm[])
{
    gmx_bool     bPull;     /* Is pulling requested in .tpr file?             */
    gmx_bool     bTpi;      /* Is test particle insertion requested?          */
    gmx_bool     bFree;     /* Is a free energy simulation requested?         */
    gmx_bool     bNM;       /* Is a normal mode analysis requested?           */
    t_inputrec   ir;
    t_state      state;
    gmx_mtop_t   mtop;


    /* Check tpr file for options that trigger extra output files */
    read_tpx_state(opt2fn("-s",nfile,fnm),&ir,&state,NULL,&mtop);
    bPull = (epullNO != ir.ePull);
    bFree = (efepNO  != ir.efep );
    bNM   = (eiNM    == ir.eI   );
    bTpi  = EI_TPI(ir.eI);

    /* Set these output files on the tuning command-line */
    if (bPull)
    {
        setopt("-pf"  , nfile, fnm);
        setopt("-px"  , nfile, fnm);
    }
    if (bFree)
    {
        setopt("-dhdl", nfile, fnm);
    }
    if (bTpi)
    {
        setopt("-tpi" , nfile, fnm);
        setopt("-tpid", nfile, fnm);
    }
    if (bNM)
    {
        setopt("-mtx" , nfile, fnm);
    }
}


static void couple_files_options(int nfile, t_filenm fnm[])
{
    int i;
    gmx_bool bSet,bBench;
    char *opt;
    char buf[20];
    
    
    for (i=0; i<nfile; i++)
    {
        opt  = (char *)fnm[i].opt;
        bSet = ((fnm[i].flag & ffSET) != 0);
        bBench = (0 == strncmp(opt,"-b", 2));

        /* Check optional files */
        /* If e.g. -eo is set, then -beo also needs to be set */
        if (is_optional(&fnm[i]) && bSet && !bBench)
        {
            sprintf(buf, "-b%s", &opt[1]);
            setopt(buf,nfile,fnm);
        }
        /* If -beo is set, then -eo also needs to be! */
        if (is_optional(&fnm[i]) && bSet && bBench)
        {
            sprintf(buf, "-%s", &opt[2]);
            setopt(buf,nfile,fnm);
        }
    }
}


static double gettime()
{
#ifdef HAVE_GETTIMEOFDAY
    struct timeval t;
    double seconds;

    gettimeofday(&t,NULL);
    
    seconds = (double) t.tv_sec + 1e-6*(double)t.tv_usec;
    
    return seconds;
#else
    double  seconds;
    
    seconds = time(NULL);
    
    return seconds;
#endif
}


#define BENCHSTEPS (1000)

int gmx_tune_pme(int argc,char *argv[])
{
    const char *desc[] = {
            "For a given number [TT]-np[tt] or [TT]-nt[tt] of processors/threads, this program systematically",
            "times [TT]mdrun[tt] with various numbers of PME-only nodes and determines",
            "which setting is fastest. It will also test whether performance can",
            "be enhanced by shifting load from the reciprocal to the real space",
            "part of the Ewald sum. ",
            "Simply pass your [TT].tpr[tt] file to [TT]g_tune_pme[tt] together with other options",
            "for [TT]mdrun[tt] as needed.[PAR]",
            "Which executables are used can be set in the environment variables",
            "MPIRUN and MDRUN. If these are not present, 'mpirun' and 'mdrun'",
            "will be used as defaults. Note that for certain MPI frameworks you",
            "need to provide a machine- or hostfile. This can also be passed",
            "via the MPIRUN variable, e.g.[PAR]",
            "[TT]export MPIRUN=\"/usr/local/mpirun -machinefile hosts\"[tt][PAR]",
            "Please call [TT]g_tune_pme[tt] with the normal options you would pass to",
            "[TT]mdrun[tt] and add [TT]-np[tt] for the number of processors to perform the",
            "tests on, or [TT]-nt[tt] for the number of threads. You can also add [TT]-r[tt]",
            "to repeat each test several times to get better statistics. [PAR]",
            "[TT]g_tune_pme[tt] can test various real space / reciprocal space workloads",
            "for you. With [TT]-ntpr[tt] you control how many extra [TT].tpr[tt] files will be",
            "written with enlarged cutoffs and smaller fourier grids respectively.",
            "Typically, the first test (number 0) will be with the settings from the input",
            "[TT].tpr[tt] file; the last test (number [TT]ntpr[tt]) will have cutoffs multiplied",
            "by (and at the same time fourier grid dimensions divided by) the scaling",
            "factor [TT]-fac[tt] (default 1.2). The remaining [TT].tpr[tt] files will have about ",
            "equally-spaced values in between these extremes. [BB]Note[bb] that you can set [TT]-ntpr[tt] to 1",
            "if you just want to find the optimal number of PME-only nodes; in that case",
            "your input [TT].tpr[tt] file will remain unchanged.[PAR]",
            "For the benchmark runs, the default of 1000 time steps should suffice for most",
            "MD systems. The dynamic load balancing needs about 100 time steps",
            "to adapt to local load imbalances, therefore the time step counters",
            "are by default reset after 100 steps. For large systems",
            "(>1M atoms) you may have to set [TT]-resetstep[tt] to a higher value.",
            "From the 'DD' load imbalance entries in the md.log output file you",
            "can tell after how many steps the load is sufficiently balanced. Example call:[PAR]"
            "[TT]g_tune_pme -np 64 -s protein.tpr -launch[tt][PAR]",
            "After calling [TT]mdrun[tt] several times, detailed performance information",
            "is available in the output file [TT]perf.out.[tt] ",
            "[BB]Note[bb] that during the benchmarks, a couple of temporary files are written",
            "(options [TT]-b[tt]*), these will be automatically deleted after each test.[PAR]",
            "If you want the simulation to be started automatically with the",
            "optimized parameters, use the command line option [TT]-launch[tt].[PAR]",
    };

    int        nnodes =1;
    int        repeats=2;
    int        pmeentries=0; /* How many values for -npme do we actually test for each tpr file */
    real       maxPMEfraction=0.50;
    real       minPMEfraction=0.25;
    int        maxPMEnodes, minPMEnodes;
    int        npme_fixed=-2;             /* If >= -1, use only this number
                                           * of PME-only nodes                */
    real       downfac=1.0,upfac=1.2;
    int        ntprs=0;
    real       fs=0.0;                    /* 0 indicates: not set by the user */
    gmx_large_int_t bench_nsteps=BENCHSTEPS;
    gmx_large_int_t new_sim_nsteps=-1;   /* -1 indicates: not set by the user */
    gmx_large_int_t cpt_steps=0;         /* Step counter in .cpt input file   */
    int        presteps=100;    /* Do a full cycle reset after presteps steps */
    gmx_bool       bOverwrite=FALSE, bKeepTPR;
    gmx_bool       bLaunch=FALSE;
    gmx_bool       bPassAll=FALSE;
    char       *ExtraArgs=NULL;
    char       **tpr_names=NULL;
    const char *simulation_tpr=NULL;
    int        best_npme, best_tpr;
    int        sim_part = 1;     /* For benchmarks with checkpoint files */
    
    /* Default program names if nothing else is found */
    char        *cmd_mpirun=NULL, *cmd_mdrun=NULL;
    char        *cmd_args_bench, *cmd_args_launch;
    char        *cmd_np=NULL;

    t_perf      **perfdata=NULL;
    t_inputinfo *info;
    int         i;
    FILE        *fp;
    t_commrec   *cr;

    /* Print out how long the tuning took */
    double      seconds;
    
    static t_filenm fnm[] = {
      /* g_tune_pme */
      { efOUT, "-p",      "perf",     ffWRITE },
      { efLOG, "-err",    "errors",   ffWRITE },
      { efTPX, "-so",     "tuned",    ffWRITE },
      /* mdrun: */
      { efTPX, NULL,      NULL,       ffREAD },
      { efTRN, "-o",      NULL,       ffWRITE },
      { efXTC, "-x",      NULL,       ffOPTWR },
      { efCPT, "-cpi",    NULL,       ffOPTRD },
      { efCPT, "-cpo",    NULL,       ffOPTWR },
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
      { efEDO, "-eo",     "sam",      ffOPTWR },
      { efGCT, "-j",      "wham",     ffOPTRD },
      { efGCT, "-jo",     "bam",      ffOPTWR },
      { efXVG, "-ffout",  "gct",      ffOPTWR },
      { efXVG, "-devout", "deviatie", ffOPTWR },
      { efXVG, "-runav",  "runaver",  ffOPTWR },
      { efXVG, "-px",     "pullx",    ffOPTWR },
      { efXVG, "-pf",     "pullf",    ffOPTWR },
      { efMTX, "-mtx",    "nm",       ffOPTWR },
      { efNDX, "-dn",     "dipole",   ffOPTWR },
      /* Output files that are deleted after each benchmark run */
      { efTRN, "-bo",     "bench",    ffWRITE },
      { efXTC, "-bx",     "bench",    ffWRITE },
      { efCPT, "-bcpo",   "bench",    ffWRITE },
      { efSTO, "-bc",     "bench",    ffWRITE },
      { efEDR, "-be",     "bench",    ffWRITE },
      { efLOG, "-bg",     "bench",    ffWRITE },
      { efEDO, "-beo",    "bench",    ffOPTWR },
      { efXVG, "-bdhdl",  "benchdhdl",ffOPTWR },
      { efXVG, "-bfield", "benchfld" ,ffOPTWR },
      { efXVG, "-btpi",   "benchtpi", ffOPTWR },
      { efXVG, "-btpid",  "benchtpid",ffOPTWR },
      { efGCT, "-bjo",    "bench",    ffOPTWR },
      { efXVG, "-bffout", "benchgct", ffOPTWR },
      { efXVG, "-bdevout","benchdev", ffOPTWR },
      { efXVG, "-brunav", "benchrnav",ffOPTWR },
      { efXVG, "-bpx",    "benchpx",  ffOPTWR },
      { efXVG, "-bpf",    "benchpf",  ffOPTWR },
      { efMTX, "-bmtx",   "benchn",   ffOPTWR },
      { efNDX, "-bdn",    "bench",    ffOPTWR }
    };

    /* Command line options of mdrun */
    gmx_bool bDDBondCheck = TRUE;
    gmx_bool bDDBondComm  = TRUE;
    gmx_bool bVerbose     = FALSE;
    gmx_bool bCompact     = TRUE;
    gmx_bool bSepPot      = FALSE;
    gmx_bool bRerunVSite  = FALSE;
    gmx_bool bIonize      = FALSE;
    gmx_bool bConfout     = TRUE;
    gmx_bool bReproducible = FALSE;
    gmx_bool bThreads     = FALSE;

    int  nmultisim=0;
    int  nstglobalcomm=-1;
    int  repl_ex_nst=0;
    int  repl_ex_seed=-1;
    int  nstepout=100;
    int  nthreads=1;

    const char *ddno_opt[ddnoNR+1] =
      { NULL, "interleave", "pp_pme", "cartesian", NULL };
    const char *dddlb_opt[] =
      { NULL, "auto", "no", "yes", NULL };
    const char *procstring[] =
      { NULL, "-np", "-n", "none", NULL };
    const char *npmevalues_opt[] =
      { NULL, "auto", "all", "subset", NULL };
    real rdd=0.0,rconstr=0.0,dlb_scale=0.8,pforce=-1;
    char *ddcsx=NULL,*ddcsy=NULL,*ddcsz=NULL;
    char *deffnm=NULL;
#define STD_CPT_PERIOD (15.0)
    real cpt_period=STD_CPT_PERIOD,max_hours=-1;
    gmx_bool bAppendFiles=TRUE;
    gmx_bool bKeepAndNumCPT=FALSE;
    gmx_bool bResetCountersHalfWay=FALSE;
    output_env_t oenv=NULL;

    t_pargs pa[] = {
      /***********************/
      /* g_tune_pme options: */
      /***********************/
      { "-np",       FALSE, etINT,  {&nnodes},
        "Number of nodes to run the tests on (must be > 2 for separate PME nodes)" },
      { "-npstring", FALSE, etENUM, {procstring},
        "Specify the number of processors to [TT]$MPIRUN[tt] using this string"},
      { "-passall",  FALSE, etBOOL, {&bPassAll},
        "HIDDENPut arguments unknown to [TT]mdrun[tt] at the end of the command line. Can be used for debugging purposes. "},
      { "-nt",       FALSE, etINT,  {&nthreads},
        "Number of threads to run the tests on (turns MPI & mpirun off)"},
      { "-r",        FALSE, etINT,  {&repeats},
        "Repeat each test this often" },
      { "-max",      FALSE, etREAL, {&maxPMEfraction},
        "Max fraction of PME nodes to test with" },
      { "-min",      FALSE, etREAL, {&minPMEfraction},
        "Min fraction of PME nodes to test with" },
      { "-npme",     FALSE, etENUM, {npmevalues_opt},
        "Benchmark all possible values for [TT]-npme[tt] or just the subset that is expected to perform well"},
      { "-fix",      FALSE, etINT, {&npme_fixed},
        "If >= -1, do not vary the number of PME-only nodes, instead use this fixed value and only vary rcoulomb and the PME grid spacing."},
      { "-upfac",    FALSE, etREAL, {&upfac},
        "Upper limit for rcoulomb scaling factor (Note that rcoulomb upscaling results in fourier grid downscaling)" },
      { "-downfac",  FALSE, etREAL, {&downfac},
        "Lower limit for rcoulomb scaling factor" },
      { "-ntpr",     FALSE, etINT,  {&ntprs},
        "Number of [TT].tpr[tt] files to benchmark. Create this many files with scaling factors ranging from 1.0 to fac. If < 1, automatically choose the number of [TT].tpr[tt] files to test" },
      { "-four",     FALSE, etREAL, {&fs},
        "Use this fourierspacing value instead of the grid found in the [TT].tpr[tt] input file. (Spacing applies to a scaling factor of 1.0 if multiple [TT].tpr[tt] files are written)" },
      { "-steps",    FALSE, etGMX_LARGE_INT, {&bench_nsteps},
        "Take timings for this many steps in the benchmark runs" }, 
      { "-resetstep",FALSE, etINT,  {&presteps},
        "Let dlb equilibrate this many steps before timings are taken (reset cycle counters after this many steps)" },
      { "-simsteps", FALSE, etGMX_LARGE_INT, {&new_sim_nsteps},
        "If non-negative, perform this many steps in the real run (overwrites nsteps from [TT].tpr[tt], add [TT].cpt[tt] steps)" }, 
      { "-launch",   FALSE, etBOOL, {&bLaunch},
        "Launch the real simulation after optimization" },
      /******************/
      /* mdrun options: */
      /******************/
      { "-deffnm",    FALSE, etSTR, {&deffnm},
          "Set the default filename for all file options at launch time" },
      { "-ddorder",   FALSE, etENUM, {ddno_opt},
        "DD node order" },
      { "-ddcheck",   FALSE, etBOOL, {&bDDBondCheck},
        "Check for all bonded interactions with DD" },
      { "-ddbondcomm",FALSE, etBOOL, {&bDDBondComm},
        "HIDDENUse special bonded atom communication when [TT]-rdd[tt] > cut-off" },
      { "-rdd",       FALSE, etREAL, {&rdd},
        "The maximum distance for bonded interactions with DD (nm), 0 is determine from initial coordinates" },
      { "-rcon",      FALSE, etREAL, {&rconstr},
        "Maximum distance for P-LINCS (nm), 0 is estimate" },
      { "-dlb",       FALSE, etENUM, {dddlb_opt},
        "Dynamic load balancing (with DD)" },
      { "-dds",       FALSE, etREAL, {&dlb_scale},
        "Minimum allowed dlb scaling of the DD cell size" },
      { "-ddcsx",     FALSE, etSTR,  {&ddcsx},
        "HIDDENThe DD cell sizes in x" },
      { "-ddcsy",     FALSE, etSTR,  {&ddcsy},
        "HIDDENThe DD cell sizes in y" },
      { "-ddcsz",     FALSE, etSTR,  {&ddcsz},
        "HIDDENThe DD cell sizes in z" },
      { "-gcom",      FALSE, etINT,  {&nstglobalcomm},
        "Global communication frequency" },
      { "-v",         FALSE, etBOOL, {&bVerbose},
        "Be loud and noisy" },
      { "-compact",   FALSE, etBOOL, {&bCompact},
        "Write a compact log file" },
      { "-seppot",    FALSE, etBOOL, {&bSepPot},
        "Write separate V and dVdl terms for each interaction type and node to the log file(s)" },
      { "-pforce",    FALSE, etREAL, {&pforce},
        "Print all forces larger than this (kJ/mol nm)" },
      { "-reprod",    FALSE, etBOOL, {&bReproducible},
        "Try to avoid optimizations that affect binary reproducibility" },
      { "-cpt",       FALSE, etREAL, {&cpt_period},
        "Checkpoint interval (minutes)" },
      { "-cpnum",   FALSE, etBOOL, {&bKeepAndNumCPT},
        "Keep and number checkpoint files" },
      { "-append",    FALSE, etBOOL, {&bAppendFiles},
        "Append to previous output files when continuing from checkpoint instead of adding the simulation part number to all file names (for launch only)" },
      { "-maxh",      FALSE, etREAL, {&max_hours},
        "Terminate after 0.99 times this time (hours)" },
      { "-multi",     FALSE, etINT,  {&nmultisim},
        "Do multiple simulations in parallel" },
      { "-replex",    FALSE, etINT,  {&repl_ex_nst},
        "Attempt replica exchange periodically with this period (steps)" },
      { "-reseed",    FALSE, etINT,  {&repl_ex_seed},
        "Seed for replica exchange, -1 is generate a seed" },
      { "-rerunvsite", FALSE, etBOOL, {&bRerunVSite},
        "HIDDENRecalculate virtual site coordinates with [TT]-rerun[tt]" },
      { "-ionize",    FALSE, etBOOL, {&bIonize},
        "Do a simulation including the effect of an X-ray bombardment on your system" },
      { "-confout",   FALSE, etBOOL, {&bConfout},
        "HIDDENWrite the last configuration with [TT]-c[tt] and force checkpointing at the last step" },
      { "-stepout",   FALSE, etINT,  {&nstepout},
        "HIDDENFrequency of writing the remaining runtime" },
      { "-resethway", FALSE, etBOOL, {&bResetCountersHalfWay},
        "HIDDENReset the cycle counters after half the number of steps or halfway [TT]-maxh[tt] (launch only)" }
    };

    
#define NFILE asize(fnm)

    CopyRight(stderr,argv[0]);

    seconds = gettime();

    parse_common_args(&argc,argv,PCA_NOEXIT_ON_ARGS,
                      NFILE,fnm,asize(pa),pa,asize(desc),desc,
                      0,NULL,&oenv);        

    /* Store the remaining unparsed command line entries in a string */
    snew(ExtraArgs, 1);
    ExtraArgs[0] = '\0';
    for (i=1; i<argc; i++) /* argc will now be 1 if everything was understood */
    {
        add_to_string(&ExtraArgs, argv[i]);
        add_to_string(&ExtraArgs, " ");
    }
    if ( !bPassAll && (ExtraArgs[0] != '\0') )
    {
        fprintf(stderr, "\nWARNING: The following arguments you provided have no effect:\n"
                        "%s\n"
                        "Use the -passall option to force them to appear on the command lines\n"
                        "for the benchmark simulations%s.\n\n",
                        ExtraArgs, bLaunch? " and at launch time" : "");
    }

    if (opt2parg_bSet("-nt",asize(pa),pa))
    {
        bThreads=TRUE;
        if (opt2parg_bSet("-npstring",asize(pa),pa))
            fprintf(stderr, "WARNING: -npstring has no effect when using threads.\n");

        if (nnodes > 1)
            gmx_fatal(FARGS, "Can't run multi-threaded MPI simulation yet!");
        /* and now we just set this; a bit of an ugly hack*/
        nnodes=nthreads;
    }
    /* tpr-triggered output files */
    get_tpr_outfiles(NFILE,fnm);

    /* Automatically set -beo options if -eo is set etc. */
    couple_files_options(NFILE,fnm);
    
    /* Construct the command line arguments for benchmark runs 
     * as well as for the simulation run 
     */
    create_command_line_snippets(bThreads,presteps,NFILE,fnm,asize(pa),pa,procstring[0],
                                 &cmd_np, &cmd_args_bench, &cmd_args_launch,
                                 bPassAll? ExtraArgs : (char *)"");

    /* Read in checkpoint file if requested */
    sim_part = 1;
    if(opt2bSet("-cpi",NFILE,fnm))
    {
        snew(cr,1);
        cr->duty=DUTY_PP; /* makes the following routine happy */
        read_checkpoint_simulation_part(opt2fn("-cpi",NFILE,fnm),
					&sim_part,&cpt_steps,cr,
					FALSE,NFILE,fnm,NULL,NULL);
        sfree(cr);
        sim_part++;
        /* sim_part will now be 1 if no checkpoint file was found */
        if (sim_part<=1)
            gmx_fatal(FARGS, "Checkpoint file %s not found!", opt2fn("-cpi",
                                                                     NFILE,
                                                                     fnm));
    }
    
    /* Open performance output file and write header info */
    fp = ffopen(opt2fn("-p",NFILE,fnm),"w");
    
    /* Make a quick consistency check of command line parameters */
    check_input(nnodes, repeats, &ntprs, &upfac, &downfac,
                maxPMEfraction, minPMEfraction, npme_fixed,
                fs, bench_nsteps, fnm, NFILE, sim_part, presteps,
                asize(pa),pa);
    
    /* Determine the maximum and minimum number of PME nodes to test,
     * the actual list of settings is build in do_the_tests(). */
    if ((nnodes > 2) && (npme_fixed < -1))
    {
        maxPMEnodes = floor(maxPMEfraction*nnodes);
        minPMEnodes = max(floor(minPMEfraction*nnodes), 0);
        fprintf(stdout, "Will try runs with %d ", minPMEnodes);
        if (maxPMEnodes != minPMEnodes)
            fprintf(stdout, "- %d ", maxPMEnodes);
        fprintf(stdout, "PME-only nodes.\n  Note that the automatic number of PME-only nodes and no separate PME nodes are always tested.\n");
    }
    else
    {
        maxPMEnodes = 0;
        minPMEnodes = 0;
    }   
    
    /* Get the commands we need to set up the runs from environment variables */
    get_program_paths(bThreads, &cmd_mpirun, cmd_np, &cmd_mdrun, repeats);
    
    /* Print some header info to file */
    sep_line(fp);
    fprintf(fp, "\n      P E R F O R M A N C E   R E S U L T S\n");
    sep_line(fp);
    fprintf(fp, "%s for Gromacs %s\n", ShortProgram(),GromacsVersion());
    if (!bThreads)
    {
        fprintf(fp, "Number of nodes         : %d\n", nnodes);
        fprintf(fp, "The mpirun command is   : %s\n", cmd_mpirun);
        if ( strcmp(procstring[0], "none") != 0)
            fprintf(fp, "Passing # of nodes via  : %s\n", procstring[0]);
        else
            fprintf(fp, "Not setting number of nodes in system call\n");
    }
    else
        fprintf(fp, "Number of threads       : %d\n", nnodes);

    fprintf(fp, "The mdrun  command is   : %s\n", cmd_mdrun);
    fprintf(fp, "mdrun args benchmarks   : %s\n", cmd_args_bench);
    fprintf(fp, "Benchmark steps         : ");
    fprintf(fp, gmx_large_int_pfmt, bench_nsteps);
    fprintf(fp, "\n");
    fprintf(fp, "dlb equilibration steps : %d\n", presteps);
    if (sim_part > 1)
    {
        fprintf(fp, "Checkpoint time step    : ");
        fprintf(fp, gmx_large_int_pfmt, cpt_steps);
        fprintf(fp, "\n");
    }
    if (bLaunch)
        fprintf(fp, "mdrun args at launchtime: %s\n", cmd_args_launch);
    if (!bPassAll && ExtraArgs[0] != '\0')
        fprintf(fp, "Unused arguments        : %s\n", ExtraArgs);
    if (new_sim_nsteps >= 0)
    {
        bOverwrite = TRUE;
        fprintf(stderr, "Note: Simulation input file %s will have ", opt2fn("-so",NFILE,fnm));
        fprintf(stderr, gmx_large_int_pfmt, new_sim_nsteps+cpt_steps);
        fprintf(stderr, " steps.\n");
        fprintf(fp, "Simulation steps        : ");
        fprintf(fp, gmx_large_int_pfmt, new_sim_nsteps);
        fprintf(fp, "\n");
    }   
    if (repeats > 1)
        fprintf(fp, "Repeats for each test   : %d\n", repeats);
    
    if (fs > 0.0)
    {
        fprintf(fp, "Requested grid spacing  : %f (This will be the grid spacing at a scaling factor of 1.0)\n", fs);
    }
    
    if (npme_fixed >= -1)
    {
        fprintf(fp, "Fixing -npme at         : %d\n", npme_fixed);
    }

    fprintf(fp, "Input file              : %s\n", opt2fn("-s",NFILE,fnm));

    /* Allocate memory for the inputinfo struct: */
    snew(info, 1);
    info->nr_inputfiles = ntprs;
    for (i=0; i<ntprs; i++)
    {
        snew(info->rcoulomb , ntprs);
        snew(info->rvdw     , ntprs);
        snew(info->rlist    , ntprs);
        snew(info->rlistlong, ntprs);
        snew(info->nkx      , ntprs);
        snew(info->nky      , ntprs);
        snew(info->nkz      , ntprs);
        snew(info->fsx      , ntprs);
        snew(info->fsy      , ntprs);
        snew(info->fsz      , ntprs);
    }
    /* Make alternative tpr files to test: */
    snew(tpr_names, ntprs);
    for (i=0; i<ntprs; i++)
        snew(tpr_names[i], STRLEN);

    /* It can be that ntprs is reduced by make_benchmark_tprs if not enough
     * different grids could be found. */
    make_benchmark_tprs(opt2fn("-s",NFILE,fnm), tpr_names, bench_nsteps+presteps,
            cpt_steps, upfac, downfac, &ntprs, fs, info, fp);


    /********************************************************************************/
    /* Main loop over all scenarios we need to test: tpr files, PME nodes, repeats  */
    /********************************************************************************/
    snew(perfdata, ntprs);
    do_the_tests(fp, tpr_names, maxPMEnodes, minPMEnodes, npme_fixed, npmevalues_opt[0], perfdata, &pmeentries,
                 repeats, nnodes, ntprs, bThreads, cmd_mpirun, cmd_np, cmd_mdrun,
                 cmd_args_bench, fnm, NFILE, sim_part, presteps, cpt_steps);
    
    fprintf(fp, "\nTuning took%8.1f minutes.\n", (gettime()-seconds)/60.0);

    /* Analyse the results and give a suggestion for optimal settings: */
    bKeepTPR = analyze_data(fp, opt2fn("-p", NFILE, fnm), perfdata, nnodes, ntprs, pmeentries,
                            repeats, info, &best_tpr, &best_npme);
    
    /* Take the best-performing tpr file and enlarge nsteps to original value */
    if ( bKeepTPR && !bOverwrite && !(fs > 0.0) )
    {
        simulation_tpr = opt2fn("-s",NFILE,fnm);
    }
    else
    {
        simulation_tpr = opt2fn("-so",NFILE,fnm);
        modify_PMEsettings(bOverwrite? (new_sim_nsteps+cpt_steps) : info->orig_sim_steps,
                info->orig_init_step, tpr_names[best_tpr], simulation_tpr);
    }

    /* Now start the real simulation if the user requested it ... */
    launch_simulation(bLaunch, fp, bThreads, cmd_mpirun, cmd_np, cmd_mdrun,
                      cmd_args_launch, simulation_tpr, nnodes, best_npme);
    ffclose(fp);
        
    /* ... or simply print the performance results to screen: */
    if (!bLaunch)
        finalize(opt2fn("-p", NFILE, fnm));
    
    return 0;
}
