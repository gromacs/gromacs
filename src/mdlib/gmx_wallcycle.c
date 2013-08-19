/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2008, The GROMACS development team,
 * check out http://www.gromacs.org for more information.
 * Copyright (c) 2012,2013, by the GROMACS development team, led by
 * David van der Spoel, Berk Hess, Erik Lindahl, and including many
 * others, as listed in the AUTHORS file in the top-level source
 * directory and at http://www.gromacs.org.
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


#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <string.h>
#include "gmx_wallcycle.h"
#include "gmx_cyclecounter.h"
#include "smalloc.h"
#include "gmx_fatal.h"
#include "md_logging.h"
#include "string2.h"

#ifdef GMX_LIB_MPI
#include <mpi.h>
#endif
#ifdef GMX_THREAD_MPI
#include "tmpi.h"
#endif

/* DEBUG_WCYCLE adds consistency checking for the counters.
 * It checks if you stop a counter different from the last
 * one that was opened and if you do nest too deep.
 */
/* #define DEBUG_WCYCLE */

typedef struct
{
    int          n;
    gmx_cycles_t c;
    gmx_cycles_t start;
    gmx_cycles_t last;
} wallcc_t;

typedef struct gmx_wallcycle
{
    wallcc_t        *wcc;
    /* variables for testing/debugging */
    gmx_bool         wc_barrier;
    wallcc_t        *wcc_all;
    int              wc_depth;
#ifdef DEBUG_WCYCLE
#define DEPTH_MAX 6
    int               counterlist[DEPTH_MAX];
    int               count_depth;
#endif
    int               ewc_prev;
    gmx_cycles_t      cycle_prev;
    gmx_large_int_t   reset_counters;
#ifdef GMX_MPI
    MPI_Comm          mpi_comm_mygroup;
#endif
    int               nthreads_pp;
    int               nthreads_pme;
#ifdef GMX_CYCLE_SUBCOUNTERS
    wallcc_t         *wcsc;
#endif
    double           *cycles_sum;
} gmx_wallcycle_t_t;

/* Each name should not exceed 19 printing characters
   (ie. terminating null can be twentieth) */
static const char *wcn[ewcNR] =
{
    "Run", "Step", "PP during PME", "Domain decomp.", "DD comm. load",
    "DD comm. bounds", "Vsite constr.", "Send X to PME", "Neighbor search", "Launch GPU ops.",
    "Comm. coord.", "Born radii", "Force", "Wait + Comm. F", "PME mesh",
    "PME redist. X/F", "PME spread/gather", "PME 3D-FFT", "PME 3D-FFT Comm.", "PME solve",
    "PME wait for PP", "Wait + Recv. PME F", "Wait GPU nonlocal", "Wait GPU local", "NB X/F buffer ops.",
    "Vsite spread", "Write traj.", "Update", "Constraints", "Comm. energies",
    "Enforced rotation", "Add rot. forces", "Test"
};

static const char *wcsn[ewcsNR] =
{
    "DD redist.", "DD NS grid + sort", "DD setup comm.",
    "DD make top.", "DD make constr.", "DD top. other",
    "NS grid local", "NS grid non-loc.", "NS search local", "NS search non-loc.",
    "Bonded F", "Nonbonded F", "Ewald F correction",
    "NB X buffer ops.", "NB F buffer ops."
};

gmx_bool wallcycle_have_counter(void)
{
    return gmx_cycles_have_counter();
}

gmx_wallcycle_t wallcycle_init(FILE *fplog, int resetstep, t_commrec *cr,
                               int nthreads_pp, int nthreads_pme)
{
    gmx_wallcycle_t wc;


    if (!wallcycle_have_counter())
    {
        return NULL;
    }

    snew(wc, 1);

    wc->wc_barrier          = FALSE;
    wc->wcc_all             = NULL;
    wc->wc_depth            = 0;
    wc->ewc_prev            = -1;
    wc->reset_counters      = resetstep;
    wc->nthreads_pp         = nthreads_pp;
    wc->nthreads_pme        = nthreads_pme;
    wc->cycles_sum          = NULL;

#ifdef GMX_MPI
    if (PAR(cr) && getenv("GMX_CYCLE_BARRIER") != NULL)
    {
        if (fplog)
        {
            fprintf(fplog, "\nWill call MPI_Barrier before each cycle start/stop call\n\n");
        }
        wc->wc_barrier       = TRUE;
        wc->mpi_comm_mygroup = cr->mpi_comm_mygroup;
    }
#endif

    snew(wc->wcc, ewcNR);
    if (getenv("GMX_CYCLE_ALL") != NULL)
    {
        if (fplog)
        {
            fprintf(fplog, "\nWill time all the code during the run\n\n");
        }
        snew(wc->wcc_all, ewcNR*ewcNR);
    }

#ifdef GMX_CYCLE_SUBCOUNTERS
    snew(wc->wcsc, ewcsNR);
#endif

#ifdef DEBUG_WCYCLE
    wc->count_depth = 0;
#endif

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
#ifdef GMX_CYCLE_SUBCOUNTERS
    if (wc->wcsc != NULL)
    {
        sfree(wc->wcsc);
    }
#endif
    sfree(wc);
}

static void wallcycle_all_start(gmx_wallcycle_t wc, int ewc, gmx_cycles_t cycle)
{
    wc->ewc_prev   = ewc;
    wc->cycle_prev = cycle;
}

static void wallcycle_all_stop(gmx_wallcycle_t wc, int ewc, gmx_cycles_t cycle)
{
    wc->wcc_all[wc->ewc_prev*ewcNR+ewc].n += 1;
    wc->wcc_all[wc->ewc_prev*ewcNR+ewc].c += cycle - wc->cycle_prev;
}


#ifdef DEBUG_WCYCLE
static void debug_start_check(gmx_wallcycle_t wc, int ewc)
{
    /* fprintf(stderr,"wcycle_start depth %d, %s\n",wc->count_depth,wcn[ewc]); */

    if (wc->count_depth < 0 || wc->count_depth >= DEPTH_MAX)
    {
        gmx_fatal(FARGS, "wallcycle counter depth out of range: %d",
                  wc->count_depth);
    }
    wc->counterlist[wc->count_depth] = ewc;
    wc->count_depth++;
}

static void debug_stop_check(gmx_wallcycle_t wc, int ewc)
{
    wc->count_depth--;

    /* fprintf(stderr,"wcycle_stop depth %d, %s\n",wc->count_depth,wcn[ewc]); */

    if (wc->count_depth < 0)
    {
        gmx_fatal(FARGS, "wallcycle counter depth out of range when stopping %s: %d", wcn[ewc], wc->count_depth);
    }
    if (wc->counterlist[wc->count_depth] != ewc)
    {
        gmx_fatal(FARGS, "wallcycle mismatch at stop, start %s, stop %s",
                  wcn[wc->counterlist[wc->count_depth]], wcn[ewc]);
    }
}
#endif

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

#ifdef DEBUG_WCYCLE
    debug_start_check(wc, ewc);
#endif

    cycle              = gmx_cycles_read();
    wc->wcc[ewc].start = cycle;
    if (wc->wcc_all != NULL)
    {
        wc->wc_depth++;
        if (ewc == ewcRUN)
        {
            wallcycle_all_start(wc, ewc, cycle);
        }
        else if (wc->wc_depth == 3)
        {
            wallcycle_all_stop(wc, ewc, cycle);
        }
    }
}

void wallcycle_start_nocount(gmx_wallcycle_t wc, int ewc)
{
    if (wc == NULL)
    {
        return;
    }

    wallcycle_start(wc, ewc);
    wc->wcc[ewc].n--;
}

double wallcycle_stop(gmx_wallcycle_t wc, int ewc)
{
    gmx_cycles_t cycle, last;

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

#ifdef DEBUG_WCYCLE
    debug_stop_check(wc, ewc);
#endif

    cycle           = gmx_cycles_read();
    last            = cycle - wc->wcc[ewc].start;
    wc->wcc[ewc].c += last;
    wc->wcc[ewc].n++;
    if (wc->wcc_all)
    {
        wc->wc_depth--;
        if (ewc == ewcRUN)
        {
            wallcycle_all_stop(wc, ewc, cycle);
        }
        else if (wc->wc_depth == 2)
        {
            wallcycle_all_start(wc, ewc, cycle);
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

    for (i = 0; i < ewcNR; i++)
    {
        wc->wcc[i].n = 0;
        wc->wcc[i].c = 0;
    }
    if (wc->wcc_all)
    {
        for (i = 0; i < ewcNR*ewcNR; i++)
        {
            wc->wcc_all[i].n = 0;
            wc->wcc_all[i].c = 0;
        }
    }
#ifdef GMX_CYCLE_SUBCOUNTERS
    for (i = 0; i < ewcsNR; i++)
    {
        wc->wcsc[i].n = 0;
        wc->wcsc[i].c = 0;
    }
#endif
}

static gmx_bool is_pme_counter(int ewc)
{
    return (ewc >= ewcPMEMESH && ewc <= ewcPMEWAITCOMM);
}

static gmx_bool is_pme_subcounter(int ewc)
{
    return (ewc >= ewcPME_REDISTXF && ewc < ewcPMEWAITCOMM);
}

void wallcycle_sum(t_commrec *cr, gmx_wallcycle_t wc)
{
    wallcc_t *wcc;
    double   *cycles;
    double    cycles_n[ewcNR+ewcsNR], buf[ewcNR+ewcsNR], *cyc_all, *buf_all;
    int       i, j;
    int       nsum;

    if (wc == NULL)
    {
        return;
    }

    snew(wc->cycles_sum, ewcNR+ewcsNR);
    cycles = wc->cycles_sum;

    wcc = wc->wcc;

    for (i = 0; i < ewcNR; i++)
    {
        if (is_pme_counter(i) || (i == ewcRUN && cr->duty == DUTY_PME))
        {
            wcc[i].c *= wc->nthreads_pme;

            if (wc->wcc_all)
            {
                for (j = 0; j < ewcNR; j++)
                {
                    wc->wcc_all[i*ewcNR+j].c *= wc->nthreads_pme;
                }
            }
        }
        else
        {
            wcc[i].c *= wc->nthreads_pp;

            if (wc->wcc_all)
            {
                for (j = 0; j < ewcNR; j++)
                {
                    wc->wcc_all[i*ewcNR+j].c *= wc->nthreads_pp;
                }
            }
        }
    }

    if (wcc[ewcDDCOMMLOAD].n > 0)
    {
        wcc[ewcDOMDEC].c -= wcc[ewcDDCOMMLOAD].c;
    }
    if (wcc[ewcDDCOMMBOUND].n > 0)
    {
        wcc[ewcDOMDEC].c -= wcc[ewcDDCOMMBOUND].c;
    }
    if (wcc[ewcPME_FFTCOMM].n > 0)
    {
        wcc[ewcPME_FFT].c -= wcc[ewcPME_FFTCOMM].c;
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
    for (i = 0; i < ewcNR; i++)
    {
        cycles_n[i] = (double)wcc[i].n;
        cycles[i]   = (double)wcc[i].c;
    }
    nsum = ewcNR;
#ifdef GMX_CYCLE_SUBCOUNTERS
    for (i = 0; i < ewcsNR; i++)
    {
        wc->wcsc[i].c    *= wc->nthreads_pp;
        cycles_n[ewcNR+i] = (double)wc->wcsc[i].n;
        cycles[ewcNR+i]   = (double)wc->wcsc[i].c;
    }
    nsum += ewcsNR;
#endif

#ifdef GMX_MPI
    if (cr->nnodes > 1)
    {
        MPI_Allreduce(cycles_n, buf, nsum, MPI_DOUBLE, MPI_MAX,
                      cr->mpi_comm_mysim);
        for (i = 0; i < ewcNR; i++)
        {
            wcc[i].n = (int)(buf[i] + 0.5);
        }
#ifdef GMX_CYCLE_SUBCOUNTERS
        for (i = 0; i < ewcsNR; i++)
        {
            wc->wcsc[i].n = (int)(buf[ewcNR+i] + 0.5);
        }
#endif

        MPI_Allreduce(cycles, buf, nsum, MPI_DOUBLE, MPI_SUM,
                      cr->mpi_comm_mysim);
        for (i = 0; i < nsum; i++)
        {
            cycles[i] = buf[i];
        }

        if (wc->wcc_all != NULL)
        {
            snew(cyc_all, ewcNR*ewcNR);
            snew(buf_all, ewcNR*ewcNR);
            for (i = 0; i < ewcNR*ewcNR; i++)
            {
                cyc_all[i] = wc->wcc_all[i].c;
            }
            MPI_Allreduce(cyc_all, buf_all, ewcNR*ewcNR, MPI_DOUBLE, MPI_SUM,
                          cr->mpi_comm_mysim);
            for (i = 0; i < ewcNR*ewcNR; i++)
            {
                wc->wcc_all[i].c = buf_all[i];
            }
            sfree(buf_all);
            sfree(cyc_all);
        }
    }
#endif
}

static void print_cycles(FILE *fplog, double c2t, const char *name,
                         int nthreads_tot,
                         int nnodes, int nthreads,
                         int n, double c, double tot)
{
    char   num[11];
    char   thstr[6];
    double wallt;

    if (c > 0)
    {
        if (n > 0)
        {
            snprintf(num, sizeof(num), "%10d", n);
            if (nthreads < 0)
            {
                snprintf(thstr, sizeof(thstr), "N/A");
            }
            else
            {
                snprintf(thstr, sizeof(thstr), "%4d", nthreads);
            }
        }
        else
        {
            sprintf(num, "          ");
            sprintf(thstr, "    ");
        }
        /* Convert the cycle count to wallclock time for this task */
        if (nthreads > 0)
        {
            /* Cycle count has been multiplied by the thread count,
             * correct for the number of threads used.
             */
            wallt = c*c2t*nthreads_tot/(double)(nnodes*nthreads);
        }
        else
        {
            /* nthreads=-1 signals total run time, no correction required */
            wallt = c*c2t;
        }
        fprintf(fplog, " %-19.19s %4d %4s %10s  %10.3f %12.3f   %5.1f\n",
                name, nnodes, thstr, num, wallt, c*1e-9, 100*c/tot);
    }
}

static void print_gputimes(FILE *fplog, const char *name,
                           int n, double t, double tot_t)
{
    char num[11];
    char avg_perf[11];

    if (n > 0)
    {
        snprintf(num, sizeof(num), "%10d", n);
        snprintf(avg_perf, sizeof(avg_perf), "%10.3f", t/n);
    }
    else
    {
        sprintf(num, "          ");
        sprintf(avg_perf, "          ");
    }
    if (t != tot_t)
    {
        fprintf(fplog, " %-29s %10s%12.3f   %s   %5.1f\n",
                name, num, t/1000, avg_perf, 100 * t/tot_t);
    }
    else
    {
        fprintf(fplog, " %-29s %10s%12.3f   %s   %5.1f\n",
                name, "", t/1000, avg_perf, 100.0);
    }
}

void wallcycle_print(FILE *fplog, int nnodes, int npme, double realtime,
                     gmx_wallcycle_t wc, wallclock_gpu_t *gpu_t)
{
    double     *cycles;
    double      c2t, tot, tot_gpu, tot_cpu_overlap, gpu_cpu_ratio, sum, tot_k;
    int         i, j, npp, nth_pp, nth_pme, nth_tot;
    char        buf[STRLEN];
    const char *hline = "-----------------------------------------------------------------------------";

    if (wc == NULL)
    {
        return;
    }

    nth_pp  = wc->nthreads_pp;
    nth_pme = wc->nthreads_pme;

    cycles = wc->cycles_sum;

    if (npme > 0)
    {
        npp = nnodes - npme;

        nth_tot = npp*nth_pp + npme*nth_pme;
    }
    else
    {
        npp  = nnodes;
        npme = nnodes;

        nth_tot = npp*nth_pp;
    }

    tot = cycles[ewcRUN];

    /* Conversion factor from cycles to seconds */
    if (tot > 0)
    {
        c2t = realtime/tot;
    }
    else
    {
        c2t = 0;
    }

    fprintf(fplog, "\n     R E A L   C Y C L E   A N D   T I M E   A C C O U N T I N G\n\n");

    fprintf(fplog, " Computing:         Nodes   Th.     Count  Wall t (s)     G-Cycles       %c\n", '%');
    fprintf(fplog, "%s\n", hline);
    sum = 0;
    for (i = ewcPPDURINGPME+1; i < ewcNR; i++)
    {
        if (!is_pme_subcounter(i))
        {
            print_cycles(fplog, c2t, wcn[i], nth_tot,
                         is_pme_counter(i) ? npme : npp,
                         is_pme_counter(i) ? nth_pme : nth_pp,
                         wc->wcc[i].n, cycles[i], tot);
            sum += cycles[i];
        }
    }
    if (wc->wcc_all != NULL)
    {
        for (i = 0; i < ewcNR; i++)
        {
            for (j = 0; j < ewcNR; j++)
            {
                snprintf(buf, 20, "%-9.9s %-9.9s", wcn[i], wcn[j]);
                print_cycles(fplog, c2t, buf, nth_tot,
                             is_pme_counter(i) ? npme : npp,
                             is_pme_counter(i) ? nth_pme : nth_pp,
                             wc->wcc_all[i*ewcNR+j].n,
                             wc->wcc_all[i*ewcNR+j].c,
                             tot);
            }
        }
    }
    print_cycles(fplog, c2t, "Rest", nth_tot, npp, -1, 0, tot-sum, tot);
    fprintf(fplog, "%s\n", hline);
    print_cycles(fplog, c2t, "Total", nth_tot, nnodes, -1, 0, tot, tot);
    fprintf(fplog, "%s\n", hline);

    if (wc->wcc[ewcPMEMESH].n > 0)
    {
        fprintf(fplog, "%s\n", hline);
        for (i = ewcPPDURINGPME+1; i < ewcNR; i++)
        {
            if (is_pme_subcounter(i))
            {
                print_cycles(fplog, c2t, wcn[i], nth_tot,
                             is_pme_counter(i) ? npme : npp,
                             is_pme_counter(i) ? nth_pme : nth_pp,
                             wc->wcc[i].n, cycles[i], tot);
            }
        }
        fprintf(fplog, "%s\n", hline);
    }

#ifdef GMX_CYCLE_SUBCOUNTERS
    fprintf(fplog, "%s\n", hline);
    for (i = 0; i < ewcsNR; i++)
    {
        print_cycles(fplog, c2t, wcsn[i], nth_tot, npp, nth_pp,
                     wc->wcsc[i].n, cycles[ewcNR+i], tot);
    }
    fprintf(fplog, "%s\n", hline);
#endif

    /* print GPU timing summary */
    if (gpu_t)
    {
        const char *k_log_str[2][2] = {
            {"Nonbonded F kernel", "Nonbonded F+ene k."},
            {"Nonbonded F+prune k.", "Nonbonded F+ene+prune k."}
        };

        tot_gpu = gpu_t->pl_h2d_t + gpu_t->nb_h2d_t + gpu_t->nb_d2h_t;

        /* add up the kernel timings */
        tot_k = 0.0;
        for (i = 0; i < 2; i++)
        {
            for (j = 0; j < 2; j++)
            {
                tot_k += gpu_t->ktime[i][j].t;
            }
        }
        tot_gpu += tot_k;

        tot_cpu_overlap = wc->wcc[ewcFORCE].c;
        if (wc->wcc[ewcPMEMESH].n > 0)
        {
            tot_cpu_overlap += wc->wcc[ewcPMEMESH].c;
        }
        tot_cpu_overlap *= c2t * 1000; /* convert s to ms */

        fprintf(fplog, "\n GPU timings\n%s\n", hline);
        fprintf(fplog, " Computing:                         Count  Wall t (s)      ms/step       %c\n", '%');
        fprintf(fplog, "%s\n", hline);
        print_gputimes(fplog, "Pair list H2D",
                       gpu_t->pl_h2d_c, gpu_t->pl_h2d_t, tot_gpu);
        print_gputimes(fplog, "X / q H2D",
                       gpu_t->nb_c, gpu_t->nb_h2d_t, tot_gpu);

        for (i = 0; i < 2; i++)
        {
            for (j = 0; j < 2; j++)
            {
                if (gpu_t->ktime[i][j].c)
                {
                    print_gputimes(fplog, k_log_str[i][j],
                                   gpu_t->ktime[i][j].c, gpu_t->ktime[i][j].t, tot_gpu);
                }
            }
        }

        print_gputimes(fplog, "F D2H",  gpu_t->nb_c, gpu_t->nb_d2h_t, tot_gpu);
        fprintf(fplog, "%s\n", hline);
        print_gputimes(fplog, "Total ", gpu_t->nb_c, tot_gpu, tot_gpu);
        fprintf(fplog, "%s\n", hline);

        gpu_cpu_ratio = tot_gpu/tot_cpu_overlap;
        fprintf(fplog, "\nForce evaluation time GPU/CPU: %.3f ms/%.3f ms = %.3f\n",
                tot_gpu/gpu_t->nb_c, tot_cpu_overlap/wc->wcc[ewcFORCE].n,
                gpu_cpu_ratio);

        /* only print notes related to CPU-GPU load balance with PME */
        if (wc->wcc[ewcPMEMESH].n > 0)
        {
            fprintf(fplog, "For optimal performance this ratio should be close to 1!\n");

            /* print note if the imbalance is high with PME case in which
             * CPU-GPU load balancing is possible */
            if (gpu_cpu_ratio < 0.75 || gpu_cpu_ratio > 1.2)
            {
                /* Only the sim master calls this function, so always print to stderr */
                if (gpu_cpu_ratio < 0.75)
                {
                    if (npp > 1)
                    {
                        /* The user could have used -notunepme,
                         * but we currently can't check that here.
                         */
                        md_print_warn(NULL, fplog,
                                      "\nNOTE: The GPU has >25%% less load than the CPU. This imbalance causes\n"
                                      "      performance loss. Maybe the domain decomposition limits the PME tuning.\n"
                                      "      In that case, try setting the DD grid manually (-dd) or lowering -dds.");
                    }
                    else
                    {
                        /* We should not end up here, unless the box is
                         * too small for increasing the cut-off for PME tuning.
                         */
                        md_print_warn(NULL, fplog,
                                      "\nNOTE: The GPU has >25%% less load than the CPU. This imbalance causes\n"
                                      "      performance loss.");
                    }
                }
                if (gpu_cpu_ratio > 1.2)
                {
                    md_print_warn(NULL, fplog,
                                  "\nNOTE: The GPU has >20%% more load than the CPU. This imbalance causes\n"
                                  "      performance loss, consider using a shorter cut-off and a finer PME grid.");
                }
            }
        }
    }

    if (wc->wcc[ewcNB_XF_BUF_OPS].n > 0 &&
        (cycles[ewcDOMDEC] > tot*0.1 ||
         cycles[ewcNS] > tot*0.1))
    {
        /* Only the sim master calls this function, so always print to stderr */
        if (wc->wcc[ewcDOMDEC].n == 0)
        {
            md_print_warn(NULL, fplog,
                          "NOTE: %d %% of the run time was spent in pair search,\n"
                          "      you might want to increase nstlist (this has no effect on accuracy)\n",
                          (int)(100*cycles[ewcNS]/tot+0.5));
        }
        else
        {
            md_print_warn(NULL, fplog,
                          "NOTE: %d %% of the run time was spent in domain decomposition,\n"
                          "      %d %% of the run time was spent in pair search,\n"
                          "      you might want to increase nstlist (this has no effect on accuracy)\n",
                          (int)(100*cycles[ewcDOMDEC]/tot+0.5),
                          (int)(100*cycles[ewcNS]/tot+0.5));
        }
    }

    if (cycles[ewcMoveE] > tot*0.05)
    {
        /* Only the sim master calls this function, so always print to stderr */
        md_print_warn(NULL, fplog,
                      "NOTE: %d %% of the run time was spent communicating energies,\n"
                      "      you might want to use the -gcom option of mdrun\n",
                      (int)(100*cycles[ewcMoveE]/tot+0.5));
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
    {
        return;
    }

    wc->reset_counters = reset_counters;
}

#ifdef GMX_CYCLE_SUBCOUNTERS

void wallcycle_sub_start(gmx_wallcycle_t wc, int ewcs)
{
    if (wc != NULL)
    {
        wc->wcsc[ewcs].start = gmx_cycles_read();
    }
}

void wallcycle_sub_stop(gmx_wallcycle_t wc, int ewcs)
{
    if (wc != NULL)
    {
        wc->wcsc[ewcs].c += gmx_cycles_read() - wc->wcsc[ewcs].start;
        wc->wcsc[ewcs].n++;
    }
}

#endif /* GMX_CYCLE_SUBCOUNTERS */
