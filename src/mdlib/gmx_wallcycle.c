/*  -*- mode: c; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; c-file-style: "stroustrup"; -*-
 *
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

#include <string.h>
#include "gmx_wallcycle.h"
#include "gmx_cyclecounter.h"
#include "smalloc.h"
#include "gmx_fatal.h"

#ifdef GMX_LIB_MPI
#include <mpi.h>
#endif
#ifdef GMX_THREADS
#include "tmpi.h"
#endif

typedef struct
{
    int          n;
    gmx_cycles_t c;
    gmx_cycles_t start;
    gmx_cycles_t last;
} wallcc_t;

typedef struct gmx_wallcycle
{
    wallcc_t     *wcc;
    /* variables for testing/debugging */
    gmx_bool         wc_barrier;
    wallcc_t     *wcc_all;
    int          wc_depth;
    int          ewc_prev;
    gmx_cycles_t cycle_prev;
    gmx_large_int_t   reset_counters;
#ifdef GMX_MPI
    MPI_Comm     mpi_comm_mygroup;
#endif
} gmx_wallcycle_t_t;

/* Each name should not exceed 19 characters */
static const char *wcn[ewcNR] =
{ "Run", "Step", "PP during PME", "Domain decomp.", "DD comm. load", "DD comm. bounds", "Vsite constr.", "Send X to PME", "Comm. coord.", "Neighbor search", "Born radii", "Force", "Wait + Comm. F", "PME mesh", "PME redist. X/F", "PME spread/gather", "PME 3D-FFT", "PME solve", "Wait + Comm. X/F", "Wait + Recv. PME F", "Vsite spread", "Write traj.", "Update", "Constraints", "Comm. energies", "Test" };

gmx_bool wallcycle_have_counter(void)
{
  return gmx_cycles_have_counter();
}

gmx_wallcycle_t wallcycle_init(FILE *fplog,int resetstep,t_commrec *cr)
{
    gmx_wallcycle_t wc;
    
    
    if (!wallcycle_have_counter())
    {
        return NULL;
    }

    snew(wc,1);

    wc->wc_barrier = FALSE;
    wc->wcc_all    = NULL;
    wc->wc_depth   = 0;
    wc->ewc_prev   = -1;
    wc->reset_counters = resetstep;

#ifdef GMX_MPI
    if (PAR(cr) && getenv("GMX_CYCLE_BARRIER") != NULL)
    {
        if (fplog) 
        {
            fprintf(fplog,"\nWill call MPI_Barrier before each cycle start/stop call\n\n");
        }
        wc->wc_barrier = TRUE;
        wc->mpi_comm_mygroup = cr->mpi_comm_mygroup;
    }
#endif

    snew(wc->wcc,ewcNR);
    if (getenv("GMX_CYCLE_ALL") != NULL)
    {
/*#ifndef GMX_THREADS*/
        if (fplog) 
        {
            fprintf(fplog,"\nWill time all the code during the run\n\n");
        }
        snew(wc->wcc_all,ewcNR*ewcNR);
/*#else*/
        gmx_fatal(FARGS, "GMX_CYCLE_ALL is incompatible with threaded code");
/*#endif*/
    }
    
    return wc;
}

void wallcycle_destroy(gmx_wallcycle_t wc)
{
    if (wc == NULL)
    {
        return;
    }
    
    if (wc->wcc != NULL)
    {
        sfree(wc->wcc);
    }
    if (wc->wcc_all != NULL)
    {
        sfree(wc->wcc_all);
    }
    sfree(wc);
}

static void wallcycle_all_start(gmx_wallcycle_t wc,int ewc,gmx_cycles_t cycle)
{
    wc->ewc_prev = ewc;
    wc->cycle_prev = cycle;
}

static void wallcycle_all_stop(gmx_wallcycle_t wc,int ewc,gmx_cycles_t cycle)
{
    wc->wcc_all[wc->ewc_prev*ewcNR+ewc].n += 1;
    wc->wcc_all[wc->ewc_prev*ewcNR+ewc].c += cycle - wc->cycle_prev;
}

void wallcycle_start(gmx_wallcycle_t wc, int ewc)
{
    gmx_cycles_t cycle;

    if (wc == NULL)
    {
        return;
    }

#ifdef GMX_MPI
    if (wc->wc_barrier)
    {
        MPI_Barrier(wc->mpi_comm_mygroup);
    }
#endif

    cycle = gmx_cycles_read();
    wc->wcc[ewc].start = cycle;
    if (wc->wcc_all != NULL)
    {
        wc->wc_depth++;
        if (ewc == ewcRUN)
        {
            wallcycle_all_start(wc,ewc,cycle);
        }
        else if (wc->wc_depth == 3)
        {
            wallcycle_all_stop(wc,ewc,cycle);
        }
    }
}

double wallcycle_stop(gmx_wallcycle_t wc, int ewc)
{
    gmx_cycles_t cycle,last;
    
    if (wc == NULL)
    {
        return 0;
    }
    
#ifdef GMX_MPI
    if (wc->wc_barrier)
    {
        MPI_Barrier(wc->mpi_comm_mygroup);
    }
#endif
    
    cycle = gmx_cycles_read();
    last = cycle - wc->wcc[ewc].start;
    wc->wcc[ewc].c += last;
    wc->wcc[ewc].n++;
    if (wc->wcc_all)
    {
        wc->wc_depth--;
        if (ewc == ewcRUN)
        {
            wallcycle_all_stop(wc,ewc,cycle);
        }
        else if (wc->wc_depth == 2)
        {
            wallcycle_all_start(wc,ewc,cycle);
        }
    }

    return last;
}

void wallcycle_reset_all(gmx_wallcycle_t wc)
{
    int i;

    if (wc == NULL)
    {
        return;
    }

    for(i=0; i<ewcNR; i++)
    {
        wc->wcc[i].n = 0;
        wc->wcc[i].c = 0;
        wc->wcc[i].start = 0;
        wc->wcc[i].last = 0;
    }
}

void wallcycle_sum(t_commrec *cr, gmx_wallcycle_t wc,double cycles[])
{
    wallcc_t *wcc;
    double cycles_n[ewcNR],buf[ewcNR],*cyc_all,*buf_all;
    int    i;

    if (wc == NULL)
    {
        return;
    }

    wcc = wc->wcc;

    if (wcc[ewcDDCOMMLOAD].n > 0)
    {
        wcc[ewcDOMDEC].c -= wcc[ewcDDCOMMLOAD].c;
    }
    if (wcc[ewcDDCOMMBOUND].n > 0)
    {
        wcc[ewcDOMDEC].c -= wcc[ewcDDCOMMBOUND].c;
    }
    if (cr->npmenodes == 0)
    {
        /* All nodes do PME (or no PME at all) */
        if (wcc[ewcPMEMESH].n > 0)
        {
            wcc[ewcFORCE].c -= wcc[ewcPMEMESH].c;
        }
    }
    else
    {
        /* The are PME-only nodes */
        if (wcc[ewcPMEMESH].n > 0)
        {
            /* This must be a PME only node, calculate the Wait + Comm. time */
            wcc[ewcPMEWAITCOMM].c = wcc[ewcRUN].c - wcc[ewcPMEMESH].c;
        }
    }
    
    /* Store the cycles in a double buffer for summing */
    for(i=0; i<ewcNR; i++)
    {
        cycles_n[i] = (double)wcc[i].n;
        cycles[i]   = (double)wcc[i].c;
    }
    
#ifdef GMX_MPI
    if (cr->nnodes > 1)
    {
        MPI_Allreduce(cycles_n,buf,ewcNR,MPI_DOUBLE,MPI_MAX,
                      cr->mpi_comm_mysim);
        for(i=0; i<ewcNR; i++)
        {
            wcc[i].n = (int)(buf[i] + 0.5);
        }
        MPI_Allreduce(cycles,buf,ewcNR,MPI_DOUBLE,MPI_SUM,
                      cr->mpi_comm_mysim);
        for(i=0; i<ewcNR; i++)
        {
            cycles[i] = buf[i];
        }

        if (wc->wcc_all != NULL)
        {
            snew(cyc_all,ewcNR*ewcNR);
            snew(buf_all,ewcNR*ewcNR);
            for(i=0; i<ewcNR*ewcNR; i++)
            {
                cyc_all[i] = wc->wcc_all[i].c;
            }
            MPI_Allreduce(cyc_all,buf_all,ewcNR*ewcNR,MPI_DOUBLE,MPI_SUM,
                          cr->mpi_comm_mysim);
            for(i=0; i<ewcNR*ewcNR; i++)
            {
                wc->wcc_all[i].c = buf_all[i];
            }
            sfree(buf_all);
            sfree(cyc_all);
        }
    }
#endif
}

static void print_cycles(FILE *fplog, double c2t, const char *name, int nnodes,
                         int n, double c, double tot)
{
    char num[11];
  
    if (c > 0)
    {
        if (n > 0)
        {
            sprintf(num,"%10d",n);
        }
        else
        {
            sprintf(num,"          ");
        }
        fprintf(fplog," %-19s %4d %10s %12.3f %10.1f   %5.1f\n",
                name,nnodes,num,c*1e-9,c*c2t,100*c/tot);
    }
}

static gmx_bool subdivision(int ewc)
{
    return (ewc >= ewcPME_REDISTXF && ewc <= ewcPME_SOLVE);
}

void wallcycle_print(FILE *fplog, int nnodes, int npme, double realtime,
                     gmx_wallcycle_t wc, double cycles[])
{
    double c2t,tot,sum;
    int    i,j,npp;
    char   buf[STRLEN];
    const char *myline = "-----------------------------------------------------------------------";
    
    if (wc == NULL)
    {
        return;
    }

    if (npme > 0)
    {
        npp = nnodes - npme;
    }
    else
    {
        npp  = nnodes;
        npme = nnodes;
    }
    tot = cycles[ewcRUN];
    /* Conversion factor from cycles to seconds */
    if (tot > 0)
    {
      c2t = nnodes*realtime/tot;
    }
    else
    {
      c2t = 0;
    }

    fprintf(fplog,"\n     R E A L   C Y C L E   A N D   T I M E   A C C O U N T I N G\n\n");

    fprintf(fplog," Computing:         Nodes     Number     G-Cycles    Seconds     %c\n",'%');
    fprintf(fplog,"%s\n",myline);
    sum = 0;
    for(i=ewcPPDURINGPME+1; i<ewcNR; i++)
    {
        if (!subdivision(i))
        {
            print_cycles(fplog,c2t,wcn[i],
                         (i==ewcPMEMESH || i==ewcPMEWAITCOMM) ? npme : npp,
                         wc->wcc[i].n,cycles[i],tot);
            sum += cycles[i];
        }
    }
    if (wc->wcc_all != NULL)
    {
        for(i=0; i<ewcNR; i++)
        {
            for(j=0; j<ewcNR; j++)
            {
                sprintf(buf,"%-9s",wcn[i]);
                buf[9] = ' ';
                sprintf(buf+10,"%-9s",wcn[j]);
                buf[19] = '\0';
                print_cycles(fplog,c2t,buf,
                             (i==ewcPMEMESH || i==ewcPMEWAITCOMM) ? npme : npp,
                             wc->wcc_all[i*ewcNR+j].n,
                             wc->wcc_all[i*ewcNR+j].c,
                             tot);
            }
        }
    }
    print_cycles(fplog,c2t,"Rest",npp,0,tot-sum,tot);
    fprintf(fplog,"%s\n",myline);
    print_cycles(fplog,c2t,"Total",nnodes,0,tot,tot);
    fprintf(fplog,"%s\n",myline);
    
    if (wc->wcc[ewcPMEMESH].n > 0)
    {
        fprintf(fplog,"%s\n",myline);
        for(i=ewcPPDURINGPME+1; i<ewcNR; i++)
        {
            if (subdivision(i))
            {
                print_cycles(fplog,c2t,wcn[i],
                             (i>=ewcPMEMESH && i<=ewcPME_SOLVE) ? npme : npp,
                             wc->wcc[i].n,cycles[i],tot);
            }
        }
        fprintf(fplog,"%s\n",myline);
    }

    if (cycles[ewcMoveE] > tot*0.05)
    {
        sprintf(buf,
                "NOTE: %d %% of the run time was spent communicating energies,\n"
                "      you might want to use the -gcom option of mdrun\n",
                (int)(100*cycles[ewcMoveE]/tot+0.5));
        if (fplog)
        {
            fprintf(fplog,"\n%s\n",buf);
        }
        /* Only the sim master calls this function, so always print to stderr */
        fprintf(stderr,"\n%s\n",buf);
    }
}

extern gmx_large_int_t wcycle_get_reset_counters(gmx_wallcycle_t wc)
{
    if (wc == NULL)
    {
        return -1;
    }
    
    return wc->reset_counters;
}

extern void wcycle_set_reset_counters(gmx_wallcycle_t wc, gmx_large_int_t reset_counters)
{
    if (wc == NULL)
        return;

    wc->reset_counters = reset_counters;
}
