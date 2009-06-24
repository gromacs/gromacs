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
#include <gmx_ana.h>



enum {
  ddnoSEL, ddnoINTERLEAVE, ddnoPP_PME, ddnoCARTESIAN, ddnoNR
};

/* Enum for situations that can occur during log file parsing */
enum {
    eParselogOK, eParselogNotFound, eParselogNoPerfData, eParselogTerm, eParselogNr
};


typedef struct
{
    int  nPMEnodes;       /* number of PME only nodes used in this test */
    int  nx, ny, nz;      /* DD grid */
    int  guessPME;        /* if nPMEnodes == -1, this is the guessed number of PME nodes */
    real *Gcycles;        /* This can contain more than one value if doing multiple tests */
    real Gcycles_Av;
    real *ns_per_day;
    real ns_per_day_Av;
    real *PME_f_load;     /* PME mesh/force load average*/
    real PME_f_load_Av;   /* Average average ;) ... */
    char *mdrun_cmd_line; /* Mdrun command line used for this test */
} t_perf;


typedef struct
{
    int  nr_inputfiles;         /* The number of tpr and mdp input files */
    gmx_step_t orig_sim_steps;  /* Number of steps to be done in the real simulation */
    real *r_coulomb;            /* The coulomb radii [0...nr_inputfiles] */
    real *r_vdW;                /* The vdW radii */
    int  *fourier_nx, *fourier_ny, *fourier_nz;
    real *fourier_sp;           /* Fourierspacing */
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
static bool str_starts(const char *string, const char *substring)
{
    return ( strncmp(string, substring, strlen(substring)) == 0);
}


static int parse_logfile(char *logfile, t_perf *perfdata, int test_nr)
{
    FILE  *fp;
    char  line[STRLEN], dumstring[STRLEN];
    char  matchstrdd[]="Domain decomposition grid";
    char  matchstrbal[]="Average PME mesh/force load:";
    char  matchstring[]="R E A L   C Y C L E   A N D   T I M E   A C C O U N T I N G";
    char  errTERM[]="Received the TERM signal, stopping at the next step";
    char  errUSR1[]="Received the USR1 signal, stopping at the next NS step";
    int   iFound;
    int   procs;
    real  dum1,dum2,dum3;
    int   npme;


    if (!gmx_fexist(logfile))
    {
        fprintf(stderr, "WARNING: Could not find logfile %s.\n", logfile);
        return eParselogNotFound;
    }

    fp = fopen(logfile, "r");
    perfdata->PME_f_load[test_nr] = -1.0;
    perfdata->guessPME            = -1;
    iFound = 0;

    while (fgets(line, STRLEN, fp) != NULL)
    {
        /* Remove leading spaces */
        ltrim(line);

        /* Check for TERM and USR1 signals from user: */
        if ( str_starts(line, errTERM) || str_starts(line, errUSR1) )
        {
            fclose(fp);
            return eParselogTerm;
        }

        /* Look for strings that appear in a certain order in the log file: */
        switch(iFound)
        {
            case 0: /* Look for domain decomp grid and separate PME nodes: */
                if (str_starts(line, matchstrdd))
                {
                    sscanf(line, "Domain decomposition grid %d x %d x %d, separate PME nodes %d",
                            &(perfdata->nx), &(perfdata->ny), &(perfdata->nz), &npme);
                    if (perfdata->nPMEnodes == -1)
                        perfdata->guessPME = npme;
                    else if (perfdata->nPMEnodes != npme)
                        gmx_fatal(FARGS, "PME nodes from command line and output file are not identical");
                    iFound = 1;
                }
                break;
            case 1: /* Look for PME mesh/force balance (not necessarily present, though) */
                if (str_starts(line, matchstrbal))
                    sscanf(&line[strlen(matchstrbal)], "%f", &(perfdata->PME_f_load[test_nr]));
                /* Look for matchstring */
                if (str_starts(line, matchstring))
                    iFound = 2;
                break;
            case 2: /* Already found matchstring - look for cycle data */
                if (str_starts(line, "Total  "))
                {
                    sscanf(line,"Total %d %f",&procs,&(perfdata->Gcycles[test_nr]));
                    iFound = 3;
                }
                break;
            case 3: /* Already found cycle data - look for remaining performance info and return */
                if (str_starts(line, "Performance:"))
                {
                    sscanf(line,"%s %f %f %f %f", dumstring, &dum1, &dum2, &(perfdata->ns_per_day[test_nr]), &dum3);
                    fclose(fp);
                    return eParselogOK;
                }
                break;
        }
    } /* while */
    fprintf(stdout, "No performance data in log file.\n");
    fclose(fp);
    return eParselogNoPerfData;
}


/* This function returns the command line that we need to pass on to mdrun as
 * well as to grompp when we call those programs later. We have to strip away
 * the program name as well as the "-np" and "-r" options */
static void mdrun_args(int argc,char *argv[],char **arguments_mdrun)
{
    int  cmdlength=0;
    int  i;


    /* How long is the command line? */
    for (i=1; i<argc; i++)
        cmdlength += strlen(argv[i]);
    /* Allocate space for the command line, take into account separating spaces and terminal '\0' */
    /* More space than this is never needed: */
    snew(*arguments_mdrun ,cmdlength+argc+1);

    for (i=1; i<argc; i++)
    {
        /* Skip the options which are used by g_tune_pme itself! */
        if ( (0 == strcmp(argv[i],"-np"      ))
          || (0 == strcmp(argv[i],"-r"       ))
          || (0 == strcmp(argv[i],"-ntpr"    ))
          || (0 == strcmp(argv[i],"-max"     ))
          || (0 == strcmp(argv[i],"-min"     ))
          || (0 == strcmp(argv[i],"-fac"     ))
          || (0 == strcmp(argv[i],"-four"    ))
          || (0 == strcmp(argv[i],"-s"       ))
          || (0 == strcmp(argv[i],"-steps"   ))
          || (0 == strcmp(argv[i],"-simsteps"))
          || (0 == strcmp(argv[i],"-deffnm"  )) )
            i++;    /* Move one forward for the argument! */
        else if (0 == strcmp(argv[i],"-launch"))
            ;       /* Simply ignore this one here */
        else
        {
            /* All other options must be meant for mdrun: */
            strcat(*arguments_mdrun,argv[i]);
            strcat(*arguments_mdrun," "    );
        }
    }
}


static int analyze_data(
        FILE        *fp,
        t_perf      **perfdata,
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
    real s=0.0;  /* standard deviation */
    t_perf *pd;
    char strbuf[STRLEN];


    if (nrepeats > 1)
    {
        sep_line(fp);
        fprintf(fp, "Summary of successful runs:\n");
        fprintf(fp, "Line tpr PME nodes  Gcycles Av.     Std.dev.       ns/day        PME/f   DD grid\n");
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

                    fprintf(fp, "%4d %3d %4d%s %12.3f %12.3f %12.3f %12.3f %3d %3d %3d\n",
                            line, k, pd->nPMEnodes, strbuf, pd->Gcycles_Av, s,
                            pd->ns_per_day_Av, pd->PME_f_load_Av, pd->nx, pd->ny, pd->nz);
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
        gmx_fatal(FARGS, "None of the runs was successful! Exiting.");

    sep_line(fp);

    winPME = perfdata[k_win][i_win].nPMEnodes;
    if (winPME == -1)
        sprintf(strbuf, "%s", "the automatic number of");
    else
        sprintf(strbuf, "%d", winPME);
    fprintf(fp, "Best performance was achieved with %s PME nodes", strbuf);
    if (nrepeats > 1)
        fprintf(fp, " (see line %d) ", line_win);
    if (ntprs > 1)
        fprintf(fp, "\nand %s PME settings. ", (k_win ? "optimized" : "original"));
    fprintf(fp, "\n");

    /* Only mention settings if rcoulomb, rvdv, nkx, nky, or nkz was modified: */
    if (k_win)
    {
        fprintf(fp, "Optimized PME settings:\n");
        fprintf(fp, "r_coulomb = %f, r_vdW = %f, nx,ny,nz = %d %d %d\n",
                info->r_coulomb[k_win], info->r_vdW[k_win],
                info->fourier_nx[k_win], info->fourier_ny[k_win], info->fourier_nz[k_win]);
    }
    fflush(fp);
    
    /* Return the index of the mdp file that showed the highest performance
     * and the optimal number of PME nodes */
    *index_tpr    = k_win; 
    *npme_optimal = winPME;
    
    return 0;
}


/* Get the commands we need to set up the runs from environment variables */
static void get_program_paths(char *cmd_mpirun[], char *cmd_mdrun[], int repeats)
{
    char *command=NULL;
    char *cp;
    char *cp2;
    FILE *fp;
    const char filename[] = "tune.test";
    char line[STRLEN];
    const char match_mdrun[]="NNODES=";
    bool  bFound = FALSE;
    
    
    /* Get the commands we need to set up the runs from environment variables */
     if ( (cp = getenv("MPIRUN" )) != NULL)
         *cmd_mpirun  = cp;
     if ( (cp = getenv("MDRUN"  )) != NULL )
         *cmd_mdrun   = cp;

     /* Run a small test to see if mpirun and mdrun work if we intend to execute mdrun! */
     snew(command, strlen(*cmd_mpirun) + strlen(*cmd_mdrun) + strlen(filename) + 30);
     if (repeats > 0)
     {
         sprintf(command, "%s -np 1 %s -h -quiet >& %s", *cmd_mpirun, *cmd_mdrun, filename);
         fprintf(stdout, "Making shure that mdrun can be executed: '%s' ...", command);
         
         gmx_system_call(command);
             
         /* Check if we find the gromacs header in the log file: */
         fp = fopen(filename, "r");
         while ( (!feof(fp)) && (bFound==FALSE) )
         {
             cp2=fgets(line, STRLEN, fp);
             if (cp2!=NULL && str_starts(line, match_mdrun))
                 bFound = TRUE;
         }
         if (!bFound)
             gmx_fatal(FARGS, "Cannot execute mdrun. Please check %s for problems!", filename);
         fclose(fp);
         fprintf(stdout, "passed.\n");
     }

     /* Clean up ... */
	 remove(filename);
     gmx_system_call(command);
}


static void launch_simulation(
        bool bLaunch,           /* Should the simulation be launched? */
        FILE *fp,               /* General log file */
        char *cmd_mpirun,       /* Command for mpirun */
        char *cmd_mdrun,        /* Command for mdrun */
        char *args_for_mdrun,   /* Arguments for mdrun */
        const char *deffnm,     /* Extra argument just for launch */
        char *simulation_tpr,   /* This tpr will be simulated */
        int  nnodes,            /* Number of nodes to run on */
        int  nPMEnodes)         /* Number of PME nodes to use */
{
    char  *command;
    char  *command_extra;
    
    
    /* Make enough space for the system call command, 100 extra chars for the tpr filename: */
    snew(command_extra, strlen(deffnm) + strlen("-deffnm ") + 1);
    snew(command, strlen(cmd_mpirun)+strlen(cmd_mdrun)+strlen(args_for_mdrun)+strlen(simulation_tpr)+strlen(deffnm)+100);

    /* Construct the command line to call mdrun: */
    if (strlen(deffnm) > 0)
        sprintf(command_extra, "-deffnm %s", deffnm);

    sprintf(command, "%s -np %d %s %s-npme %d -s %s %s", 
            cmd_mpirun, nnodes, cmd_mdrun, args_for_mdrun, nPMEnodes, simulation_tpr, command_extra);
 
        
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
        gmx_step_t simsteps, /* Set this value as number of time steps */
        char *fn_best_tpr,   /* tpr file with the best performance */
        char *fn_sim_tpr)    /* name of tpr file to be launched */
{
    t_inputrec   *ir;
    t_state      state;
    gmx_mtop_t   mtop;
    char         buf[200];
    
    snew(ir,1);
    read_tpx_state(fn_best_tpr,ir,&state,NULL,&mtop);
        
    /* Set nsteps to the right value */
    ir->nsteps = simsteps;
    
    /* Write the tpr file which will be launched */
    sprintf(buf, "Writing optimized simulation file %s with %s steps.\n", fn_sim_tpr, gmx_step_pfmt);
    fprintf(stdout,buf,ir->nsteps);
    fflush(stdout);
    write_tpx_state(fn_sim_tpr,ir,&state,&mtop);
        
    sfree(ir);
}


/* Make additional TPR files with more computational load for the
 * direct space processors: */
static void make_benchmark_tprs(
        char *fn_sim_tpr,       /* READ : User-provided tpr file */
        char *fn_bench_tprs[],  /* WRITE: Names of benchmark tpr files */
        gmx_step_t benchsteps,  /* Number of time steps for benchmark runs */
        real maxfac,            /* Max scaling factor for rcoulomb and fourierspacing */
        int ntprs,              /* No. of TPRs to write, each with a different rcoulomb and fourierspacing */
        real fourierspacing,    /* Basic fourierspacing from tpr input file */
        t_inputinfo *info,      /* Contains information about mdp file options */
        FILE *fp)               /* Write the output here */
{
    int          i,j,d;
    t_inputrec   *ir;
    t_state      state;
    gmx_mtop_t   mtop;
    real         fac;
    real         orig_rcoulomb, orig_rvdw, orig_rlist;
    rvec         orig_fs;      /* original fourierspacing per dimension */
    ivec         orig_nk;      /* original number of grid points per dimension */
    char         buf[200];
    real         max_spacing;
    rvec         box_size;
    

    sprintf(buf, "Making benchmark tpr files with %s steps.\n", gmx_step_pfmt);
    fprintf(stdout, buf, benchsteps);
    
    snew(ir,1);
    read_tpx_state(fn_sim_tpr,ir,&state,NULL,&mtop);

    /* Check if PME was chosen */
    if (EEL_PME(ir->coulombtype) == FALSE)
        gmx_fatal(FARGS, "Can only do optimizations for simulations with PME");
    
    /* Check if rcoulomb == rlist, which is necessary for PME */
    if (!(ir->rcoulomb == ir->rlist))
        gmx_fatal(FARGS, "PME requires rcoulomb (%f) to be equal to rlist (%f).", ir->rcoulomb, ir->rlist);

    /* Reduce the number of steps for the benchmarks */
    info->orig_sim_steps = ir->nsteps;
    ir->nsteps           = benchsteps;
    
    /* Determine lenght of triclinic box vectors */
    for(d=0; d<DIM; d++)
    {
        box_size[d] = 0;
        for(i=0;i<DIM;i++)
            box_size[d] += state.box[d][i]*state.box[d][i];
        box_size[d] = sqrt(box_size[d]);
    }
    
    /* Remember the original values: */
    orig_rvdw     = ir->rvdw;
    orig_rcoulomb = ir->rcoulomb;
    orig_rlist    = ir->rlist;
    orig_nk[XX]   = ir->nkx;
    orig_nk[YY]   = ir->nky;
    orig_nk[ZZ]   = ir->nkz;
    orig_fs[XX]   = box_size[XX]/ir->nkx;  /* fourierspacing in x direction */
    orig_fs[YY]   = box_size[YY]/ir->nky;
    orig_fs[ZZ]   = box_size[ZZ]/ir->nkz;
     
    fprintf(fp, "\nWill try these real/reciprocal workload settings:\n");
    fprintf(fp, " No. scaling   r_coul   (r_vdW)     nkx  nky  nkz   (spacing)   tpr file\n");
    
    if (ntprs > 1)
    {
        fprintf(stdout, "Calculating PME grid points on the basis of ");
        if (fourierspacing > 0)
            fprintf(stdout, "a fourierspacing of %f nm\n", fourierspacing);
        else
            fprintf(stdout, "original nkx/nky/nkz settings from tpr file\n");
    }
    
    /* Loop to create the requested number of tpr input files */
    for (j = 0; j < ntprs; j++)
    {
        /* Rcoulomb scaling factor for this file: */
        if (ntprs == 1)
            fac = 1.0;
         else
            fac = (maxfac-1.0f)/(ntprs-1) * j +1;
        fprintf(stdout, "--- Scaling factor %f ---\n", fac);
        
        ir->rcoulomb = orig_rcoulomb*fac;
        ir->rlist    = orig_rlist   *fac;
        ir->rvdw     = orig_rvdw    *fac;
        
        /* Try to reduce the number of reciprocal grid points in a smart way */
        /* Did the user supply a value for fourierspacing on the command line? */
        if (fourierspacing > 0)
        {
            info->fourier_sp[j] = fourierspacing*fac;
            /* Calculate the optimal grid dimensions */
            ir->nkx = 0;
            ir->nky = 0;
            ir->nkz = 0;
            max_spacing = calc_grid(stdout,state.box,info->fourier_sp[j],&(ir->nkx),&(ir->nky),&(ir->nkz),1);
            /* Check consistency */
            if (0 == j)
                if ((ir->nkx != orig_nk[XX]) || (ir->nky != orig_nk[YY]) || (ir->nkz != orig_nk[ZZ]))
                    gmx_fatal(FARGS, "Wrong fourierspacing %f, actual grid = %dx%dx%d, original grid = %dx%dx%d", 
                            fourierspacing,ir->nkx,ir->nky,ir->nkz,orig_nk[XX],orig_nk[YY],orig_nk[ZZ]);
        }
        else
        {
            if (0 == j)
            {
                /* Print out fourierspacing from input tpr */
                fprintf(stdout, "Input file fourier grid is %dx%dx%d\n", orig_nk[XX], orig_nk[YY], orig_nk[ZZ]);
            }
            else
            {
                /* Reconstruct fourierspacing for each dimension from the input file */
                ir->nkx=0;
                max_spacing = calc_grid(stdout,state.box,orig_fs[XX]*fac,&(ir->nkx),&(ir->nky),&(ir->nkz),1);
                ir->nky=0;
                max_spacing = calc_grid(stdout,state.box,orig_fs[XX]*fac,&(ir->nkx),&(ir->nky),&(ir->nkz),1);
                ir->nkz=0;
                max_spacing = calc_grid(stdout,state.box,orig_fs[XX]*fac,&(ir->nkx),&(ir->nky),&(ir->nkz),1);
            }
        }
        /* r_vdw should only grow if necessary! */
        if (j > 0)
        {
            ir->rvdw = min(ir->rvdw, orig_rcoulomb*fac);
            ir->rvdw = max(ir->rvdw, orig_rvdw);
        }
        /* Save modified radii and fourier grid components for later output: */
        info->r_coulomb[j] = ir->rcoulomb;
        info->r_vdW[j]     = ir->rvdw;        
        info->fourier_nx[j]= ir->nkx;
        info->fourier_ny[j]= ir->nky;
        info->fourier_nz[j]= ir->nkz;

        /* Write the benchmark tpr file */
        strncpy(fn_bench_tprs[j],fn_sim_tpr,strlen(fn_sim_tpr)-strlen(".tpr"));
        sprintf(buf, "_tune%.2d.tpr", j);
        strcat(fn_bench_tprs[j], buf);
        fprintf(stdout,"Writing benchmark tpr %s with ", fn_bench_tprs[j]);
        fprintf(stdout, gmx_step_pfmt, ir->nsteps);
        fprintf(stdout," steps, scaling factor %f\n", fac);
        write_tpx_state(fn_bench_tprs[j],ir,&state,&mtop);
        
        /* Write some info to log file */
        fprintf(fp, "%3d %9f %9f (%7f) %4d %4d %4d   %9f   %-14s\n",
                j, fac, ir->rcoulomb, ir->rvdw, ir->nkx, ir->nky, ir->nkz, info->fourier_sp[j],fn_bench_tprs[j]);
    }
    fflush(stdout);
    fflush(fp);
    
    sfree(ir);
}


/* Rename the files we want to keep to some meaningful filename and
 * delete the rest */
static void cleanup(t_filenm *fnm, int nfile, int k, int nnodes, int nPMEnodes, int nr)
{
    char numstring[STRLEN];
    char newfilename[STRLEN];
    char *fn=NULL;
    int i;
    const char *opt;


    fprintf(stdout, "Cleaning up ...\n");

    for (i=0; i<nfile; i++)
    {
        opt = (char *)fnm[i].opt;
        if ( strcmp(opt, "-p") == 0 )
        {
            /* do nothing */
            ;
        }
        else if (strcmp(opt, "-g") == 0)
        {
            /* Give the log file a nice name so one can later see which parameters were used */
            numstring[0] = '\0';
            if (nr > 0)
                sprintf(numstring, "_%d", nr);
            sprintf(newfilename, "%s_no%d_np%d_npme%d%s", opt2fn("-g",nfile,fnm), k, nnodes, nPMEnodes, numstring);
            fprintf(stdout, "renaming log file to %s\n", newfilename);
            make_backup(newfilename);
            rename(opt2fn("-g",nfile,fnm), newfilename);
        }
        /* Delete the files which are created for any run: */
        else if (strcmp(opt, "-o"  ) == 0
              || strcmp(opt, "-x"  ) == 0
              || strcmp(opt, "-cpo") == 0
              || strcmp(opt, "-c"  ) == 0
              || strcmp(opt, "-e"  ) == 0 )
        {
            fn = opt2fn(opt, nfile, fnm);
            fprintf(stdout, "Deleting %s\n", fn);
            remove(fn);
        }
        else
        {
            /* Delete all other output files (except -PPo) that have been specified on 
             * the cmd line, since these were produced during the g_tune_pme run */
            if ( opt2bSet(opt, nfile, fnm) && is_output(&fnm[i]) && (0 != strcmp(opt, "-PPo")) )
            {
                fn = opt2fn(opt, nfile, fnm);
                fprintf(stdout, "Deleting %s\n", fn);
                remove(fn);
            }
        }
    }
}


static void do_the_tests(FILE *fp, char **tpr_names, int maxPMEnodes, int minPMEnodes,
        int datasets, t_perf **perfdata, int repeats, int nnodes, int nr_tprs,
        char *cmd_mpirun, char *cmd_mdrun, char *args_for_mdrun, t_filenm *fnm, int nfile)
{
    int     i,nr,k,ret;
    int     nPMEnodes;
    t_perf  *pd=NULL;
    int     cmdline_length;
    char    *command;
    char    buf[STRLEN];

    /* This string array corresponds to the eParselog enum type from above */
    const char* ParseLog[] = {"OK", "Logfile not found", "No timings in log file", "Run was terminated"};


    /* Allocate space for the mdrun command line. 100 extra characters should be more than enough
     * for the -npme etcetera arguments */
    cmdline_length = strlen(cmd_mpirun) + strlen(cmd_mdrun) + strlen(args_for_mdrun) +100;
    snew(command, cmdline_length);

    /* Loop over all tpr files to test: */
    for (k=0; k<nr_tprs;k++)
    {
        fprintf(fp, "\nIndividual timings for input file %d (%s):\n", k, tpr_names[k]);
        fprintf(fp, "PME nodes      Gcycles       ns/day        PME/f    Remark\n");
        i=0;
        /* Start with the maximum number of PME only nodes: */
        nPMEnodes = maxPMEnodes;

        /* Loop over various numbers of PME nodes: */
        for (i = 0; i<datasets; i++)
        {
            pd = &perfdata[k][i];

            /* Loop over the repeats for each scenario: */
            for (nr = 0; nr < repeats; nr++)
            {
                pd->nPMEnodes = nPMEnodes;

                /* Construct the command line to call mdrun (and save it): */
                sprintf(pd->mdrun_cmd_line, "%s -np %d %s %s-npme %d -s %s",
                        cmd_mpirun, nnodes, cmd_mdrun, args_for_mdrun, nPMEnodes, tpr_names[k]);

                /* Do a test simulation: */
                if (repeats > 1)
                    sprintf(buf, ", pass %d/%d", nr+1, repeats);
                else
                    buf[0]='\0';
                fprintf(stdout, "\n=== tpr %d/%d, run %d/%d%s:\n", k+1, nr_tprs, i+1, datasets, buf);
                sprintf(command, "%s >& tune_mdrun.out", pd->mdrun_cmd_line);
                fprintf(stdout, "%s\n", pd->mdrun_cmd_line);
                gmx_system_call(command);

                /* Collect the performance data from the log file */
                ret = parse_logfile(ftp2fn(efLOG,nfile,fnm), pd, nr);

                if (nPMEnodes == -1)
                    sprintf(buf, "(%3d)", pd->guessPME);
                else
                    sprintf(buf, "     ");

                /* Write the data we got to disk */
                fprintf(fp, "%4d%s %12.3f %12.3f %12.3f    %s\n", pd->nPMEnodes, buf, pd->Gcycles[nr], pd->ns_per_day[nr], pd->PME_f_load[nr], ParseLog[ret]);
                fflush(fp);

                /* Do some cleaning up and delete the files we do not need any more */
                cleanup(fnm, nfile, k, nnodes, nPMEnodes, nr);

                /* If the first run with this number of processors already failed, do not try again: */
                if (pd->Gcycles[0] <= 0.0 && repeats > 1)
                {
                    fprintf(stdout, "Skipping remaining passes of unsuccessful setting, see log file for details.\n");
                    break;
                }
            }
            /* Prepare for the next number of PME only nodes */
            /* The last but one check is always without MPMD PME ... */
            if ((nPMEnodes == minPMEnodes) && (0 != minPMEnodes)) 
                nPMEnodes = 0;
            /* ... and the last check with the guessed settings */
            else if (nPMEnodes == 0)
                nPMEnodes = -1;
            else
                nPMEnodes--;
        }
    }
}


static bool is_equal(real a, real b)
{
    real diff, eps=1.0e-6;
    
    
    diff = a - b;
    
    if (diff < 0.0) diff = -diff;
    
    if (diff < eps)
        return TRUE;
    else
        return FALSE;
}


static void check_input(
        int nnodes, 
        int repeats, 
        int *ntprs, 
        real maxfac,
        real maxPMEfraction,
        real minPMEfraction,
        real fourierspacing,
        gmx_step_t bench_nsteps,
        char *fn_top)
{
    /* Make shure the input file exists */
    if (!gmx_fexist(fn_top))
        gmx_fatal(FARGS, "File %s not found.", fn_top);
    
    /* Make shure that repeats is >= 0 (if == 0, only write tpr files) */
    if (repeats < 1)
        gmx_fatal(FARGS, "Number of repeats < 1. Need at least one test for each configuration.");

    /* We need at least one topology to test */
    if ((*ntprs) <= 0)
        gmx_fatal(FARGS, "-ntpr must be larger than zero.");

    if ( (1 == *ntprs) && !is_equal(maxfac,1.0) )
        fprintf(stderr, "Note: Choose ntpr>1 to shift PME load to real space.\n");
    
    if ( is_equal(1.0,maxfac) && (*ntprs > 1) )
    {
        fprintf(stderr, "WARNING: Resetting -ntpr to 1 since upscaling factor equals unity.\n  Please select -fac>1 if you want to test various PME grid settings\n");
        *ntprs = 1;
    }
    /* Check whether we have enough nodes */
    if (nnodes < 3)
        gmx_fatal(FARGS, "Can not have separate PME nodes with 2 or less nodes, so there is nothing to optimize here.");

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
        fprintf(stderr, gmx_step_pfmt, bench_nsteps);
        fprintf(stderr, ". Are you shure you want to perform so %s steps for each benchmark?\n", (bench_nsteps < 100)? "few" : "many");
    }
    
    if (maxfac <= 0.0)
        gmx_fatal(FARGS, "Scaling factor must be larger than zero.");
    
    if (maxfac < 1.0)
        fprintf(stderr, "WARNING: A scaling factor smaller than one means that load will be shifted to reciprocal space. Are you shure you want that?\n");

    if (maxfac < 0.75 || maxfac > 1.5)
        fprintf(stderr, "WARNING: Applying extreme scaling factor. I hope you know what you are doing.\n");
    
    if (fourierspacing < 0)
        gmx_fatal(FARGS, "Please choose a positive value for fourierspacing.");
}


#define BENCHSTEPS (2500)

int gmx_tune_pme(int argc,char *argv[])
{
    const char *desc[] = {
            "For a given number [TT]-np[tt] of processors this program systematically",
            "times mdrun with various numbers of PME-only nodes and determines",
            "which setting is fastest. It will also test whether performance can",
            "be enhanced by shifting load from the reciprocal to the real space",
            "part of the Ewald sum. "
            "Simply pass your [TT].tpr[tt] file to g_tune_pme together with other options",
            "for mdrun as needed.[PAR]",
            "Which executables are used can be set in the environment variables",
            "MPIRUN and MDRUN. If these are not present, 'mpirun' and 'mdrun'",
            "will be used as defaults. Note that for certain MPI frameworks you",
            "need to provide a machine- or hostfile. This can also be passed",
            "via the MPIRUN variable, e.g.",
            "'export MPIRUN=\"/usr/local/mpirun -machinefile hosts\"'[PAR]",
            "Please call g_tune_pme with the normal options you would pass to",
            "mdrun and add [TT]-np[tt] for the number of processors to perform the",
            "tests on. You can also add [TT]-r[tt] to repeat each test several times",
            "to get better statistics. [PAR]",
            "g_tune_pme can test various real space / reciprocal space workloads",
            "for you. With [TT]-ntpr[tt] you control how many extra [TT].tpr[tt] files will be",
            "written with enlarged cutoffs and smaller fourier grids respectively.",
            "The first test (no. 0) will be with the settings from the input",
            "[TT].tpr[tt] file; the last test (no. [TT]ntpr[tt]) will have cutoffs multiplied",
            "by (and at the same time fourier grid dimensions divided by) the scaling",
            "factor [TT]-fac[tt] (default 1.2). The remaining [TT].tpr[tt] files will have equally",
            "spaced values inbetween these extremes. Note that you can set [TT]-ntpr[tt] to 1",
            "if you just want to find the optimal number of PME-only nodes; in that case",
            "your input [TT].tpr[tt] file will remain unchanged[PAR]",
            "For the benchmark runs, 2500 time steps should suffice for most MD",
            "systems. Note that dynamic load balancing needs about 100 time steps",
            "to adapt to local load imbalances. To get clean benchmark numbers,",
            "[TT]-steps[tt] should therefore always be much larger than 100![PAR]",
            "Example call: [TT]g_tune_pme -np 64 -s protein.tpr -launch[tt][PAR]",
            "After calling mdrun several times, detailed performance information",
            "is available in the output file perf.out.[PAR]"
            "If you want the simulation to be started automatically with the",
            "optimized parameters, use the command line option [TT]-launch[tt]."
    };

    static int  nnodes =3;
    static int  repeats=1;
    static real maxPMEfraction=0.50;
    static real minPMEfraction=0.25;
    static int  maxPMEnodes, minPMEnodes;
    static real maxfac=1.2;
    static int  ntprs=3;
    static real fs=0.0;            /* 0 indicates: not set by the user */
    static gmx_step_t bench_nsteps=BENCHSTEPS;
    static gmx_step_t new_sim_nsteps=-1; /* -1 indicates: not set by the user */
    bool        bOverwrite=FALSE;
    static bool bLaunch=FALSE;
    char        **tpr_names=NULL;
    char        *simulation_tpr=NULL;
    int         best_npme, best_tpr;
    
    /* Default program names if nothing else is found */
    char        *cmd_mpirun, *cmd_mdrun;
    char        def_mpirun[]="mpirun", def_mdrun[]="mdrun";

    int         cmdline_length;
    char        *command;
    t_perf      **perfdata;
    t_inputinfo *info;
    int         datasets;
    int         i,j,k;
    char        *args_for_mdrun=NULL;
    FILE        *fp;

    static t_filenm fnm[] = {
      /* g_tune_pme */
      { efOUT, "-p",      "perf",     ffWRITE },
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
      { efXVG, "-dgdl",   "dgdl",     ffOPTWR },
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
      { efNDX, "-dn",     "dipole",   ffOPTWR }
    };

    /* Command line options of mdrun */
    static bool bDDBondCheck  = TRUE;
    static bool bDDBondComm   = TRUE;
    static bool bSumEner      = TRUE;
    static bool bVerbose      = FALSE;
    static bool bCompact      = TRUE;
    static bool bSepPot       = FALSE;
    static bool bGlas         = FALSE;
    static bool bIonize       = FALSE;
    static bool bConfout      = TRUE;
    static bool bReproducible = FALSE;

    static int  nmultisim=0;
    static int  repl_ex_nst=0;
    static int  repl_ex_seed=-1;
    static int  nstepout=100;

    static const char *ddno_opt[ddnoNR+1] =
      { NULL, "interleave", "pp_pme", "cartesian", NULL };
    static const char *dddlb_opt[] =
      { NULL, "auto", "no", "yes", NULL };
    static real rdd=0.0,rconstr=0.0,dlb_scale=0.8,pforce=-1;
    static char *ddcsx=NULL,*ddcsy=NULL,*ddcsz=NULL;
    static real cpt_period=15.0,max_hours=-1;
    static bool bAppendFiles=FALSE;
    static const char *deffnm="";

    static t_pargs pa[] = {
      /***********************/
      /* g_tune_pme options: */
      /***********************/
      { "-np",       FALSE, etINT,  {&nnodes},
        "Number of nodes to run the tests on (at least 3)" },
      { "-r",        FALSE, etINT,  {&repeats},
        "Repeat each test this often" },
      { "-max",      FALSE, etREAL, {&maxPMEfraction},
        "Max fraction of PME nodes to test with" },
      { "-min",      FALSE, etREAL, {&minPMEfraction},
        "Min fraction of PME nodes to test with" },
      { "-fac",      FALSE, etREAL, {&maxfac},
        "Max upscaling factor for rcoulomb (= downscaling factor for the fourier grid)" },
      { "-ntpr",     FALSE, etINT,  {&ntprs},
        "Number of tpr files to benchmark. If >1, create files with scaling factors ranging from 1.0 to fac" },
      { "-four",     FALSE, etREAL, {&fs},
          "Fourierspacing that was chosen to create the input tpr file" },        
      { "-steps",    FALSE, etGMX_STEP_T, {&bench_nsteps},
        "Use these many steps for the benchmarks" }, 
      { "-simsteps", FALSE, etGMX_STEP_T, {&new_sim_nsteps},
          "If non-negative, overwrite nsteps from tpr for the simulation" }, 
      { "-launch",   FALSE, etBOOL, {&bLaunch},
        "Lauch the real simulation after optimization" },
      /******************/
      /* mdrun options: */
      /******************/
      { "-deffnm",    FALSE, etSTR,  {&deffnm}, 
        "Set the default filename for mdrun file options" },
      { "-ddorder",   FALSE, etENUM, {ddno_opt},
        "DD node order" },
      { "-ddcheck",   FALSE, etBOOL, {&bDDBondCheck},
        "Check for all bonded interactions with DD" },
      { "-ddbondcomm",FALSE, etBOOL, {&bDDBondComm},
        "HIDDENUse special bonded atom communication when -rdd > cut-off" },
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
      { "-sum",       FALSE, etBOOL, {&bSumEner},
        "Sum the energies at every step" },
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
      { "-append",    FALSE, etBOOL, {&bAppendFiles},
        "Append to previous output files when restarting from checkpoint" },
      { "-maxh",      FALSE, etREAL, {&max_hours},
        "Terminate after 0.99 times this time (hours)" },
      { "-multi",     FALSE, etINT,  {&nmultisim},
        "Do multiple simulations in parallel" },
      { "-replex",    FALSE, etINT,  {&repl_ex_nst},
        "Attempt replica exchange every # steps" },
      { "-reseed",    FALSE, etINT,  {&repl_ex_seed},
        "Seed for replica exchange, -1 is generate a seed" },
      { "-glas",      FALSE, etBOOL, {&bGlas},
        "Do glass simulation with special long range corrections" },
      { "-ionize",    FALSE, etBOOL, {&bIonize},
        "Do a simulation including the effect of an X-Ray bombardment on your system" },
      { "-confout",   FALSE, etBOOL, {&bConfout},
        "HIDDENWrite the last configuration with -c" },
      { "-stepout",   FALSE, etINT,  {&nstepout},
        "HIDDENFrequency of writing the remaining runtime" },
    };

    cmd_mdrun  = def_mdrun;
    cmd_mpirun = def_mpirun;

#define NFILE asize(fnm)

    CopyRight(stderr,argv[0]);

    /* Extract the parts of the command line that are meant for mdrun: */
    mdrun_args(argc,argv,&args_for_mdrun);

    parse_common_args(&argc,argv,PCA_NOEXIT_ON_ARGS,
              NFILE,fnm,asize(pa),pa,asize(desc),desc,0,NULL);        

    /* Open performance output file and write header info */
    fp = ffopen(opt2fn("-p",NFILE,fnm),"w");
    
    /* Make a quick consistency check of command line parameters */
    check_input(nnodes, repeats, &ntprs, maxfac, maxPMEfraction, minPMEfraction, fs, bench_nsteps, opt2fn("-s",NFILE,fnm));
    
    /* Determine max and min number of PME nodes to test: */
    maxPMEnodes = floor(maxPMEfraction*nnodes);
    minPMEnodes = max(floor(minPMEfraction*nnodes), 0);
    fprintf(stdout, "Will try runs with %d ", minPMEnodes);
    if (maxPMEnodes != minPMEnodes)
        fprintf(stdout, "- %d ", maxPMEnodes);
    fprintf(stdout, "PME-only nodes.\n  Note that the automatic number of PME-only nodes and no separate PME nodes are always tested.\n");
    
    /* Get the commands we need to set up the runs from environment variables */
    get_program_paths(&cmd_mpirun, &cmd_mdrun, repeats);

    /* Print some header info to file */
    sep_line(fp);
    fprintf(fp, "\n      P E R F O R M A N C E   R E S U L T S\n");
    sep_line(fp);
    fprintf(fp, "%s for Gromacs %s\n", ShortProgram(),GromacsVersion());
    fprintf(fp, "Number of nodes         : %d\n", nnodes         );
    fprintf(fp, "The mpirun command is   : %s\n", cmd_mpirun     );
    fprintf(fp, "The mdrun  command is   : %s\n", cmd_mdrun      );
    fprintf(fp, "Arguments for mdrun are : %s\n", args_for_mdrun );
    if (strlen(deffnm) > 0)
        fprintf(fp, "     plus at launchtime : %s\n", deffnm     );
    if (fs > 0.0)
        fprintf(fp, "Basic fourierspacing    : %f\n", fs         );
    fprintf(fp, "Benchmark nsteps        : ");
    fprintf(fp, gmx_step_pfmt, bench_nsteps);
    fprintf(fp, "\n");
    if (new_sim_nsteps >= 0)
    {
        bOverwrite = TRUE;
        fprintf(stderr, "Note: Simulation input file %s will have ", opt2fn("-so",NFILE,fnm));
        fprintf(stderr, gmx_step_pfmt, new_sim_nsteps);
        fprintf(stderr, " steps.\n");
        fprintf(fp, "Simulation nsteps       : ");
        fprintf(fp, gmx_step_pfmt, new_sim_nsteps);
        fprintf(fp, "\n");
    }   
    if (repeats > 1)
        fprintf(fp, "Doing %d repeats for each test.\n", repeats);
        
    /* Allocate memory for the inputinfo struct: */
    snew(info, 1);
    info->nr_inputfiles = ntprs;
    for (i=0; i<ntprs; i++)
    {
        snew(info->r_coulomb , ntprs);
        snew(info->r_vdW     , ntprs);
        snew(info->fourier_nx, ntprs);
        snew(info->fourier_ny, ntprs);
        snew(info->fourier_nz, ntprs);
        snew(info->fourier_sp, ntprs);
    }
    /* Make alternative tpr files to test: */
    snew(tpr_names, ntprs);
    for (i=0; i<ntprs; i++)
        snew(tpr_names[i], STRLEN);

    make_benchmark_tprs(opt2fn("-s",NFILE,fnm), tpr_names, bench_nsteps, maxfac, ntprs, fs, info, fp);

    /* Allocate space for the mdrun command line. 100 extra characters should be more than enough
     * for the -npme etcetera arguments */
    cmdline_length = strlen(cmd_mpirun) + strlen(cmd_mdrun) + strlen(args_for_mdrun) +100;
    snew(command, cmdline_length);

    /* Memory allocation for performance data */
    datasets = maxPMEnodes - minPMEnodes + 3;
    if (0 == minPMEnodes)
        datasets--;

    /* Allocate one dataset for each tpr input file: */
    snew(perfdata, ntprs);

    /* Allocate a subset for each test with a given number of PME nodes */
    for (k=0; k<ntprs; k++)
    {
        snew(perfdata[k], datasets);
        for (i=0; i<datasets; i++)
        {
            snew(perfdata[k][i].mdrun_cmd_line, cmdline_length);
            for (j=0; j<repeats; j++)
            {
                snew(perfdata[k][i].Gcycles   , repeats);
                snew(perfdata[k][i].ns_per_day, repeats);
                snew(perfdata[k][i].PME_f_load, repeats);
            }
        }
    }

    /********************************************************************************/
    /* Main loop over all scenarios we need to test: tpr files, PME nodes, repeats  */
    /********************************************************************************/
    do_the_tests(fp, tpr_names, maxPMEnodes, minPMEnodes, datasets, perfdata, repeats,
            nnodes, ntprs, cmd_mpirun, cmd_mdrun, args_for_mdrun, fnm, NFILE);

    /* Analyse the results and give a suggestion for optimal settings: */
    analyze_data(fp, perfdata, ntprs, datasets, repeats, info, &best_tpr, &best_npme);
    
    /* Take the best-performing tpr file and enlarge nsteps to original value */
    if ((best_tpr > 0) || bOverwrite)
    {
        simulation_tpr = opt2fn("-so",NFILE,fnm);
        modify_PMEsettings(bOverwrite? new_sim_nsteps:info->orig_sim_steps, tpr_names[best_tpr], simulation_tpr);            
    }
    else
        simulation_tpr = opt2fn("-s",NFILE,fnm);
            
    /* Now start the real simulation if the user requested it ... */
    launch_simulation(bLaunch, fp, cmd_mpirun, cmd_mdrun, args_for_mdrun, deffnm, simulation_tpr, nnodes, best_npme);
    fclose(fp);

    /* ... or simply print the performance results to screen: */
    if (!bLaunch)
    {
		FILE *fp = fopen(opt2fn("-p", NFILE, fnm),"r");
		char buf[STRLEN];
        fprintf(stdout,"\n\n");
		
		while( fgets(buf,STRLEN-1,fp) != NULL )
		{
			fprintf(stdout,"%s",buf);
		}
		fclose(fp);
        fprintf(stdout,"\n\n");
        thanx(stderr);
    }
    
    return 0;
}
