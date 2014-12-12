/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2004, The GROMACS development team.
 * Copyright (c) 2013,2014,2015, by the GROMACS development team, led by
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

#include "gromacs/legacyheaders/nrnb.h"

#include <stdlib.h>
#include <string.h>

#include "gromacs/legacyheaders/macros.h"
#include "gromacs/legacyheaders/names.h"
#include "gromacs/legacyheaders/types/commrec.h"
#include "gromacs/utility/smalloc.h"

typedef struct {
    const char *name;
    int         flop;
} t_nrnb_data;


static const t_nrnb_data nbdata[eNRNB] = {
    /* These are re-used for different NB kernels, since there are so many.
     * The actual number of flops is set dynamically.
     */
    { "NB VdW [V&F]",                    1 },
    { "NB VdW [F]",                      1 },
    { "NB Elec. [V&F]",                  1 },
    { "NB Elec. [F]",                    1 },
    { "NB Elec. [W3,V&F]",               1 },
    { "NB Elec. [W3,F]",                 1 },
    { "NB Elec. [W3-W3,V&F]",            1 },
    { "NB Elec. [W3-W3,F]",              1 },
    { "NB Elec. [W4,V&F]",               1 },
    { "NB Elec. [W4,F]",                 1 },
    { "NB Elec. [W4-W4,V&F]",            1 },
    { "NB Elec. [W4-W4,F]",              1 },
    { "NB VdW & Elec. [V&F]",            1 },
    { "NB VdW & Elec. [F]",              1 },
    { "NB VdW & Elec. [W3,V&F]",         1 },
    { "NB VdW & Elec. [W3,F]",           1 },
    { "NB VdW & Elec. [W3-W3,V&F]",      1 },
    { "NB VdW & Elec. [W3-W3,F]",        1 },
    { "NB VdW & Elec. [W4,V&F]",         1 },
    { "NB VdW & Elec. [W4,F]",           1 },
    { "NB VdW & Elec. [W4-W4,V&F]",      1 },
    { "NB VdW & Elec. [W4-W4,F]",        1 },

    { "NB Generic kernel",               1 },
    { "NB Generic charge grp kernel",    1 },
    { "NB Generic AdResS kernel",        1 },
    { "NB Free energy kernel",           1 },
    { "NB All-vs-all",                   1 },
    { "NB All-vs-all, GB",               1 },

    { "Pair Search distance check",      9 }, /* nbnxn pair dist. check */
    /* nbnxn kernel flops are based on inner-loops without exclusion checks.
     * Plain Coulomb runs through the RF kernels, except with GPUs.
     * invsqrt is counted as 6 flops: 1 for _mm_rsqt_ps + 5 for iteration.
     * The flops are equal for plain-C, x86 SIMD and GPUs, except for:
     * - plain-C kernel uses one flop more for Coulomb-only (F) than listed
     * - x86 SIMD LJ geom-comb.rule kernels (fastest) use 2 more flops
     * - x86 SIMD LJ LB-comb.rule kernels (fast) use 3 (8 for F+E) more flops
     * - GPU always does exclusions, which requires 2-4 flops, but as invsqrt
     *   is always counted as 6 flops, this roughly compensates.
     */
    { "NxN RF Elec. + LJ [F]",          38 }, /* nbnxn kernel LJ+RF, no ener */
    { "NxN RF Elec. + LJ [V&F]",        54 },
    { "NxN QSTab Elec. + LJ [F]",       41 }, /* nbnxn kernel LJ+tab, no en */
    { "NxN QSTab Elec. + LJ [V&F]",     59 },
    { "NxN Ewald Elec. + LJ [F]",       66 }, /* nbnxn kernel LJ+Ewald, no en */
    { "NxN Ewald Elec. + LJ [V&F]",    107 },
    { "NxN LJ [F]",                     33 }, /* nbnxn kernel LJ, no ener */
    { "NxN LJ [V&F]",                   43 },
    { "NxN RF Electrostatics [F]",      31 }, /* nbnxn kernel RF, no ener */
    { "NxN RF Electrostatics [V&F]",    36 },
    { "NxN QSTab Elec. [F]",            34 }, /* nbnxn kernel tab, no ener */
    { "NxN QSTab Elec. [V&F]",          41 },
    { "NxN Ewald Elec. [F]",            61 }, /* nbnxn kernel Ewald, no ener */
    { "NxN Ewald Elec. [V&F]",          84 },
    /* The switch function flops should be added to the LJ kernels above */
    { "NxN LJ add F-switch [F]",        12 }, /* extra cost for LJ F-switch */
    { "NxN LJ add F-switch [V&F]",      22 },
    { "NxN LJ add P-switch [F]",        27 }, /* extra cost for LJ P-switch */
    { "NxN LJ add P-switch [V&F]",      20 },
    { "NxN LJ add LJ Ewald [F]",        36 }, /* extra cost for LJ Ewald */
    { "NxN LJ add LJ Ewald [V&F]",      33 },
    { "1,4 nonbonded interactions",     90 },
    { "Born radii (Still)",             47 },
    { "Born radii (HCT/OBC)",          183 },
    { "Born force chain rule",          15 },
    { "All-vs-All Still radii",          1 },
    { "All-vs-All HCT/OBC radii",        1 },
    { "All-vs-All Born chain rule",      1 },
    { "Calc Weights",                   36 },
    { "Spread Q",                        6 },
    { "Spread Q Bspline",                2 },
    { "Gather F",                      23  },
    { "Gather F Bspline",              6   },
    { "3D-FFT",                        8   },
    { "Convolution",                   4   },
    { "Solve PME",                     64  },
    { "NS-Pairs",                      21  },
    { "Reset In Box",                  3   },
    { "Shift-X",                       6   },
    { "CG-CoM",                        3   },
    { "Sum Forces",                    1   },
    { "Bonds",                         59  },
    { "G96Bonds",                      44  },
    { "FENE Bonds",                    58  },
    { "Tab. Bonds",                    62  },
    { "Restraint Potential",           86  },
    { "Linear Angles",                 57  },
    { "Angles",                        168 },
    { "G96Angles",                     150 },
    { "Quartic Angles",                160 },
    { "Tab. Angles",                   169 },
    { "Propers",                       229 },
    { "Impropers",                     208 },
    { "RB-Dihedrals",                  247 },
    { "Four. Dihedrals",               247 },
    { "Tab. Dihedrals",                227 },
    { "Dist. Restr.",                  200 },
    { "Orient. Restr.",                200 },
    { "Dihedral Restr.",               200 },
    { "Pos. Restr.",                   50  },
    { "Flat-bottom posres",            50  },
    { "Angle Restr.",                  191 },
    { "Angle Restr. Z",                164 },
    { "Morse Potent.",                 83  },
    { "Cubic Bonds",                   54  },
    { "Walls",                         31  },
    { "Polarization",                  59  },
    { "Anharmonic Polarization",       72  },
    { "Water Pol.",                    62  },
    { "Thole Pol.",                    296 },
    { "Virial",                        18  },
    { "Update",                        31  },
    { "Ext.ens. Update",               54  },
    { "Stop-CM",                       10  },
    { "P-Coupling",                    6   },
    { "Calc-Ekin",                     27  },
    { "Lincs",                         60  },
    { "Lincs-Mat",                     4   },
    { "Shake",                         30  },
    { "Constraint-V",                   8  },
    { "Shake-Init",                    10  },
    { "Constraint-Vir",                24  },
    { "Settle",                        323 },
    { "Virtual Site 2",                23  },
    { "Virtual Site 3",                37  },
    { "Virtual Site 3fd",              95  },
    { "Virtual Site 3fad",             176 },
    { "Virtual Site 3out",             87  },
    { "Virtual Site 4fd",              110 },
    { "Virtual Site 4fdn",             254 },
    { "Virtual Site N",                 15 },
    { "Mixed Generalized Born stuff",   10 }
};

static void pr_two(FILE *out, int c, int i)
{
    if (i < 10)
    {
        fprintf(out, "%c0%1d", c, i);
    }
    else
    {
        fprintf(out, "%c%2d", c, i);
    }
}

static void pr_difftime(FILE *out, double dt)
{
    int        ndays, nhours, nmins, nsecs;
    gmx_bool   bPrint, bPrinted;

    ndays    = dt/(24*3600);
    dt       = dt-24*3600*ndays;
    nhours   = dt/3600;
    dt       = dt-3600*nhours;
    nmins    = dt/60;
    dt       = dt-nmins*60;
    nsecs    = dt;
    bPrint   = (ndays > 0);
    bPrinted = bPrint;
    if (bPrint)
    {
        fprintf(out, "%d", ndays);
    }
    bPrint = bPrint || (nhours > 0);
    if (bPrint)
    {
        if (bPrinted)
        {
            pr_two(out, 'd', nhours);
        }
        else
        {
            fprintf(out, "%d", nhours);
        }
    }
    bPrinted = bPrinted || bPrint;
    bPrint   = bPrint || (nmins > 0);
    if (bPrint)
    {
        if (bPrinted)
        {
            pr_two(out, 'h', nmins);
        }
        else
        {
            fprintf(out, "%d", nmins);
        }
    }
    bPrinted = bPrinted || bPrint;
    if (bPrinted)
    {
        pr_two(out, ':', nsecs);
    }
    else
    {
        fprintf(out, "%ds", nsecs);
    }
    fprintf(out, "\n");
}

void init_nrnb(t_nrnb *nrnb)
{
    int i;

    for (i = 0; (i < eNRNB); i++)
    {
        nrnb->n[i] = 0.0;
    }
}

void cp_nrnb(t_nrnb *dest, t_nrnb *src)
{
    int i;

    for (i = 0; (i < eNRNB); i++)
    {
        dest->n[i] = src->n[i];
    }
}

void add_nrnb(t_nrnb *dest, t_nrnb *s1, t_nrnb *s2)
{
    int i;

    for (i = 0; (i < eNRNB); i++)
    {
        dest->n[i] = s1->n[i]+s2->n[i];
    }
}

void print_nrnb(FILE *out, t_nrnb *nrnb)
{
    int i;

    for (i = 0; (i < eNRNB); i++)
    {
        if (nrnb->n[i] > 0)
        {
            fprintf(out, " %-26s %10.0f.\n", nbdata[i].name, nrnb->n[i]);
        }
    }
}

void _inc_nrnb(t_nrnb *nrnb, int enr, int inc, char gmx_unused *file, int gmx_unused line)
{
    nrnb->n[enr] += inc;
#ifdef DEBUG_NRNB
    printf("nrnb %15s(%2d) incremented with %8d from file %s line %d\n",
           nbdata[enr].name, enr, inc, file, line);
#endif
}

/* Returns in enr is the index of a full nbnxn VdW kernel */
static gmx_bool nrnb_is_nbnxn_vdw_kernel(int enr)
{
    return (enr >= eNR_NBNXN_LJ_RF && enr <= eNR_NBNXN_LJ_E);
}

/* Returns in enr is the index of an nbnxn kernel addition (LJ modification) */
static gmx_bool nrnb_is_nbnxn_kernel_addition(int enr)
{
    return (enr >= eNR_NBNXN_ADD_LJ_FSW && enr <= eNR_NBNXN_ADD_LJ_EWALD_E);
}

void print_flop(FILE *out, t_nrnb *nrnb, double *nbfs, double *mflop)
{
    int           i, j;
    double        mni, frac, tfrac, tflop;
    const char   *myline = "-----------------------------------------------------------------------------";

    *nbfs = 0.0;
    for (i = 0; (i < eNR_NBKERNEL_ALLVSALLGB); i++)
    {
        if (strstr(nbdata[i].name, "W3-W3") != NULL)
        {
            *nbfs += 9e-6*nrnb->n[i];
        }
        else if (strstr(nbdata[i].name, "W3") != NULL)
        {
            *nbfs += 3e-6*nrnb->n[i];
        }
        else if (strstr(nbdata[i].name, "W4-W4") != NULL)
        {
            *nbfs += 10e-6*nrnb->n[i];
        }
        else if (strstr(nbdata[i].name, "W4") != NULL)
        {
            *nbfs += 4e-6*nrnb->n[i];
        }
        else
        {
            *nbfs += 1e-6*nrnb->n[i];
        }
    }
    tflop = 0;
    for (i = 0; (i < eNRNB); i++)
    {
        tflop += 1e-6*nrnb->n[i]*nbdata[i].flop;
    }

    if (tflop == 0)
    {
        fprintf(out, "No MEGA Flopsen this time\n");
        return;
    }
    if (out)
    {
        fprintf(out, "\n\tM E G A - F L O P S   A C C O U N T I N G\n\n");
    }

    if (out)
    {
        fprintf(out, " NB=Group-cutoff nonbonded kernels    NxN=N-by-N cluster Verlet kernels\n");
        fprintf(out, " RF=Reaction-Field  VdW=Van der Waals  QSTab=quadratic-spline table\n");
        fprintf(out, " W3=SPC/TIP3p  W4=TIP4p (single or pairs)\n");
        fprintf(out, " V&F=Potential and force  V=Potential only  F=Force only\n\n");

        fprintf(out, " %-32s %16s %15s  %7s\n",
                "Computing:", "M-Number", "M-Flops", "% Flops");
        fprintf(out, "%s\n", myline);
    }
    *mflop = 0.0;
    tfrac  = 0.0;
    for (i = 0; (i < eNRNB); i++)
    {
        mni     = 1e-6*nrnb->n[i];
        /* Skip empty entries and nbnxn additional flops,
         * which have been added to the kernel entry.
         */
        if (mni > 0 && !nrnb_is_nbnxn_kernel_addition(i))
        {
            int flop;

            flop    = nbdata[i].flop;
            if (nrnb_is_nbnxn_vdw_kernel(i))
            {
                /* Possibly add the cost of an LJ switch/Ewald function */
                for (j = eNR_NBNXN_ADD_LJ_FSW; j <= eNR_NBNXN_ADD_LJ_EWALD; j += 2)
                {
                    int e_kernel_add;

                    /* Select the force or energy flop count */
                    e_kernel_add = j + ((i - eNR_NBNXN_LJ_RF) % 2);

                    if (nrnb->n[e_kernel_add] > 0)
                    {
                        flop += nbdata[e_kernel_add].flop;
                    }
                }
            }
            *mflop += mni*flop;
            frac    = 100.0*mni*flop/tflop;
            tfrac  += frac;
            if (out != NULL)
            {
                fprintf(out, " %-32s %16.6f %15.3f  %6.1f\n",
                        nbdata[i].name, mni, mni*flop, frac);
            }
        }
    }
    if (out)
    {
        fprintf(out, "%s\n", myline);
        fprintf(out, " %-32s %16s %15.3f  %6.1f\n",
                "Total", "", *mflop, tfrac);
        fprintf(out, "%s\n\n", myline);

        if (nrnb->n[eNR_NBKERNEL_GENERIC] > 0)
        {
            fprintf(out,
                    "WARNING: Using the slow generic C kernel. This is fine if you are\n"
                    "comparing different implementations or MD software. Routine\n"
                    "simulations should use a different non-bonded setup for much better\n"
                    "performance.\n\n");
        }
    }
}

void print_perf(FILE *out, double time_per_thread, double time_per_node,
                gmx_int64_t nsteps, double delta_t,
                double nbfs, double mflop)
{
    double wallclocktime;

    fprintf(out, "\n");

    if (time_per_node > 0)
    {
        fprintf(out, "%12s %12s %12s %10s\n", "", "Core t (s)", "Wall t (s)", "(%)");
        fprintf(out, "%12s %12.3f %12.3f %10.1f\n", "Time:",
                time_per_thread, time_per_node, 100.0*time_per_thread/time_per_node);
        /* only print day-hour-sec format if time_per_node is more than 30 min */
        if (time_per_node > 30*60)
        {
            fprintf(out, "%12s %12s", "", "");
            pr_difftime(out, time_per_node);
        }
        if (delta_t > 0)
        {
            mflop         = mflop/time_per_node;
            wallclocktime = nsteps*delta_t;

            if (getenv("GMX_DETAILED_PERF_STATS") == NULL)
            {
                fprintf(out, "%12s %12s %12s\n",
                        "", "(ns/day)", "(hour/ns)");
                fprintf(out, "%12s %12.3f %12.3f\n", "Performance:",
                        wallclocktime*24*3.6/time_per_node, 1000*time_per_node/(3600*wallclocktime));
            }
            else
            {
                fprintf(out, "%12s %12s %12s %12s %12s\n",
                        "", "(Mnbf/s)", (mflop > 1000) ? "(GFlops)" : "(MFlops)",
                        "(ns/day)", "(hour/ns)");
                fprintf(out, "%12s %12.3f %12.3f %12.3f %12.3f\n", "Performance:",
                        nbfs/time_per_node, (mflop > 1000) ? (mflop/1000) : mflop,
                        wallclocktime*24*3.6/time_per_node, 1000*time_per_node/(3600*wallclocktime));
            }
        }
        else
        {
            if (getenv("GMX_DETAILED_PERF_STATS") == NULL)
            {
                fprintf(out, "%12s %14s\n",
                        "", "(steps/hour)");
                fprintf(out, "%12s %14.1f\n", "Performance:",
                        nsteps*3600.0/time_per_node);
            }
            else
            {
                fprintf(out, "%12s %12s %12s %14s\n",
                        "", "(Mnbf/s)", (mflop > 1000) ? "(GFlops)" : "(MFlops)",
                        "(steps/hour)");
                fprintf(out, "%12s %12.3f %12.3f %14.1f\n", "Performance:",
                        nbfs/time_per_node, (mflop > 1000) ? (mflop/1000) : mflop,
                        nsteps*3600.0/time_per_node);
            }
        }
    }
}

int cost_nrnb(int enr)
{
    return nbdata[enr].flop;
}

const char *nrnb_str(int enr)
{
    return nbdata[enr].name;
}

static const int    force_index[] = {
    eNR_BONDS,  eNR_ANGLES,  eNR_PROPER, eNR_IMPROPER,
    eNR_RB,     eNR_DISRES,  eNR_ORIRES, eNR_POSRES,
    eNR_FBPOSRES, eNR_NS,
};
#define NFORCE_INDEX asize(force_index)

static const int    constr_index[] = {
    eNR_SHAKE,     eNR_SHAKE_RIJ, eNR_SETTLE,       eNR_UPDATE,       eNR_PCOUPL,
    eNR_CONSTR_VIR, eNR_CONSTR_V
};
#define NCONSTR_INDEX asize(constr_index)

static double pr_av(FILE *log, t_commrec *cr,
                    double fav, double ftot[], const char *title)
{
    int    i, perc;
    double dperc, unb;

    unb = 0;
    if (fav > 0)
    {
        fav /= cr->nnodes - cr->npmenodes;
        fprintf(log, "\n %-26s", title);
        for (i = 0; (i < cr->nnodes); i++)
        {
            dperc = (100.0*ftot[i])/fav;
            unb   = max(unb, dperc);
            perc  = dperc;
            fprintf(log, "%3d ", perc);
        }
        if (unb > 0)
        {
            perc = 10000.0/unb;
            fprintf(log, "%6d%%\n\n", perc);
        }
        else
        {
            fprintf(log, "\n\n");
        }
    }
    return unb;
}

void pr_load(FILE *log, t_commrec *cr, t_nrnb nrnb[])
{
    int     i, j, perc;
    double  dperc, unb, uf, us;
    double *ftot, fav;
    double *stot, sav;
    t_nrnb *av;

    snew(av, 1);
    snew(ftot, cr->nnodes);
    snew(stot, cr->nnodes);
    init_nrnb(av);
    for (i = 0; (i < cr->nnodes); i++)
    {
        add_nrnb(av, av, &(nrnb[i]));
        /* Cost due to forces */
        for (j = 0; (j < eNR_NBKERNEL_ALLVSALLGB); j++)
        {
            ftot[i] += nrnb[i].n[j]*cost_nrnb(j);
        }
        for (j = 0; (j < NFORCE_INDEX); j++)
        {
            ftot[i] += nrnb[i].n[force_index[j]]*cost_nrnb(force_index[j]);
        }
        /* Due to shake */
        for (j = 0; (j < NCONSTR_INDEX); j++)
        {
            stot[i] += nrnb[i].n[constr_index[j]]*cost_nrnb(constr_index[j]);
        }
    }
    for (j = 0; (j < eNRNB); j++)
    {
        av->n[j] = av->n[j]/(double)(cr->nnodes - cr->npmenodes);
    }

    fprintf(log, "\nDetailed load balancing info in percentage of average\n");

    fprintf(log, " Type                 RANK:");
    for (i = 0; (i < cr->nnodes); i++)
    {
        fprintf(log, "%3d ", i);
    }
    fprintf(log, "Scaling\n");
    fprintf(log, "---------------------------");
    for (i = 0; (i < cr->nnodes); i++)
    {
        fprintf(log, "----");
    }
    fprintf(log, "-------\n");

    for (j = 0; (j < eNRNB); j++)
    {
        unb = 100.0;
        if (av->n[j] > 0)
        {
            fprintf(log, " %-26s", nrnb_str(j));
            for (i = 0; (i < cr->nnodes); i++)
            {
                dperc = (100.0*nrnb[i].n[j])/av->n[j];
                unb   = max(unb, dperc);
                perc  = dperc;
                fprintf(log, "%3d ", perc);
            }
            if (unb > 0)
            {
                perc = 10000.0/unb;
                fprintf(log, "%6d%%\n", perc);
            }
            else
            {
                fprintf(log, "\n");
            }
        }
    }
    fav = sav = 0;
    for (i = 0; (i < cr->nnodes); i++)
    {
        fav += ftot[i];
        sav += stot[i];
    }
    uf = pr_av(log, cr, fav, ftot, "Total Force");
    us = pr_av(log, cr, sav, stot, "Total Constr.");

    unb = (uf*fav+us*sav)/(fav+sav);
    if (unb > 0)
    {
        unb = 10000.0/unb;
        fprintf(log, "\nTotal Scaling: %.0f%% of max performance\n\n", unb);
    }
}
