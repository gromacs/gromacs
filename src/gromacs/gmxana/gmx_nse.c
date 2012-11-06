/*
 *
 *                This source code is part of
 *
 *                 G   R   O   M   A   C   S
 *
 *          GROningen MAchine for Chemical Simulations
 *
 *                        VERSION 3.2.0
 * Written by David van der Spoel, Erik Lindahl, Berk Hess, and others.
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2004, The GROMACS development team,
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
 * Green Red Orange Magenta Azure Cyan Skyblue
 */

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <ctype.h>
#include <string.h>
#include "smalloc.h"
#include "sysstuff.h"
#include "typedefs.h"
#include "macros.h"
#include "vec.h"
#include "pbc.h"
#include "xvgr.h"
#include "copyrite.h"
#include "futil.h"
#include "statutil.h"
#include "tpxio.h"
#include "index.h"
#include "gstat.h"
#include "matio.h"
#include "gmx_ana.h"
#include "nsfactor.h"
#include "gmx_omp.h"

int gmx_nse(int argc, char *argv[])
{
    const char          *desc[] = {
        "This is simple tool to compute Neutron Spin Echo (NSE) spectra.",
        "Besides the trajectory, the topology is required to assign elements to each atom.",
        "[PAR]",
        "You need to remove jumps in trajectory before you will use this tool, so",
        "preprocess trajectory with trjconv -pbc nojump",
        "[PAR]",
        "Note: This tools produces large number of sqt files (one file per needed q value)!"
    };
    static gmx_bool      bPBC       = TRUE;
    static gmx_bool      bNSE       = TRUE;
    static gmx_bool      bNORMALIZE = TRUE;
    static real          binwidth   = 0.2, grid = 0.05; /* bins shouldnt be smaller then bond (~0.1nm) length */
    static real          start_q    = 0.01, end_q = 2.0, q_step = 0.01;
    static real          mcover     = -1;
    static unsigned int  seed       = 0;
    static int           nthreads   = -1;

    static const char   *emode[]   = { NULL, "direct", "mc", NULL };
    static const char   *emethod[] = { NULL, "debye", "fft", NULL };

    gmx_neutron_atomic_structurefactors_t    *gnsf;
    gmx_nse_t                                *gnse;

#define NPA asize(pa)

    t_pargs      pa[] = {
        { "-bin", FALSE, etREAL, {&binwidth},
          "Binwidth (nm)" },
        { "-mode", FALSE, etENUM, {emode},
          "Mode for sans spectra calculation" },
        { "-mcover", FALSE, etREAL, {&mcover},
          "Number of iterations for Monte-Carlo run"},
        { "-method", FALSE, etENUM, {emethod},
          "[HIDDEN]Method for sans spectra calculation" },
        { "-pbc", FALSE, etBOOL, {&bPBC},
          "Use periodic boundary conditions for computing distances" },
        { "-normalize", FALSE, etBOOL, {&bNORMALIZE},
          "Normalize I(q,t) output to I(q,t=0)"},
        { "-grid", FALSE, etREAL, {&grid},
          "[HIDDEN]Grid spacing (in nm) for FFTs" },
        {"-startq", FALSE, etREAL, {&start_q},
         "Starting q (1/nm) "},
        {"-endq", FALSE, etREAL, {&end_q},
         "Ending q (1/nm)"},
        { "-qstep", FALSE, etREAL, {&q_step},
          "Stepping in q (1/nm)"},
        { "-seed",     FALSE, etINT,  {&seed},
          "Random seed for Monte-Carlo"},
#ifdef GMX_OPENMP
        { "-nt",    FALSE, etINT, {&nthreads},
          "Number of threads to start"},
#endif
    };
    FILE        *fp;
    const char  *fnTPX, *fnTRX, *fnNDX, *fnDAT = NULL;
    t_trxstatus *status;
    t_topology  *top   = NULL;
    t_atom      *atom  = NULL;
    t_atoms     *atoms = NULL;
    gmx_rmpbc_t  gpbc  = NULL;
    gmx_bool     bTPX;
    gmx_bool     bFFT = FALSE, bDEBYE = FALSE;
    gmx_bool     bMC  = FALSE;
    int          ePBC = -1;
    matrix       box;
    char         title[STRLEN];
    rvec        *x, *xf;
    int          natoms;
    int          nframes;
    int          nralloc = 1;
    int          pairs;
    real         t;
    char       **grpname = NULL;
    atom_id     *index   = NULL;
    int          isize;
    int          i, j, k;
    char        *hdr    = NULL;
    char        *suffix = NULL;
    t_filenm    *fnmdup = NULL;
    gmx_radial_distribution_histogram_t   *grc = NULL;
    output_env_t oenv;

#define NFILE asize(fnm)

    t_filenm   fnm[] = {
        { efTPX,  "-s",         NULL,   ffREAD },
        { efTRX,  "-f",         NULL,   ffREAD },
        { efNDX,  NULL,         NULL,   ffOPTRD },
        { efDAT,  "-d",   "nsfactor",   ffOPTRD },
        { efXVG, "-grt",       "grt",   ffOPTWR },
        { efXVG, "-sqt",       "sqt",   ffWRITE },
        { efXVG, "-stq",       "stq",   ffOPTWR }
    };


    nthreads = gmx_omp_get_max_threads();

    parse_common_args(&argc, argv, PCA_CAN_TIME | PCA_TIME_UNIT | PCA_BE_NICE,
                      NFILE, fnm, asize(pa), pa, asize(desc), desc, 0, NULL, &oenv);

    /* Check binwith and mcover */
    check_binwidth(binwidth);
    check_mcover(mcover);

    /* setting number of omp threads globaly */
    gmx_omp_set_num_threads(nthreads);

    /* Now try to parse opts for modes */
    switch (emethod[0][0])
    {
        case 'd':
            bDEBYE = TRUE;
            switch (emode[0][0])
            {
                case 'd':
                    bMC = FALSE;
                    break;
                case 'm':
                    bMC = TRUE;
                    break;
                default:
                    break;
            }
            break;
        case 'f':
            bFFT = TRUE;
            break;
        default:
            break;
    }

    if (bDEBYE)
    {
        if (bMC)
        {
            fprintf(stderr, "Using Monte Carlo Debye method to calculate spectrum\n");
        }
        else
        {
            fprintf(stderr, "Using direct Debye method to calculate spectrum\n");
        }
    }
    else if (bFFT)
    {
        gmx_fatal(FARGS, "FFT method not implemented!");
    }
    else
    {
        gmx_fatal(FARGS, "Unknown combination for mode and method!");
    }

    /* Try to read files */
    fnDAT = ftp2fn(efDAT, NFILE, fnm);
    fnTPX = ftp2fn(efTPX, NFILE, fnm);
    fnTRX = ftp2fn(efTRX, NFILE, fnm);

    gnsf = gmx_neutronstructurefactors_init(fnDAT);
    fprintf(stderr, "Read %d atom names from %s with neutron scattering parameters\n\n", gnsf->nratoms, fnDAT);

    snew(top, 1);
    snew(gnse, 1);
    snew(grpname, 1);
    snew(index, 1);

    bTPX = read_tps_conf(fnTPX, title, top, &ePBC, &x, NULL, box, TRUE);

    atoms = &(top->atoms);

    printf("\nPlease select group for SANS spectra calculation:\n");
    get_index(&(top->atoms), ftp2fn_null(efNDX, NFILE, fnm), 1, &isize, &index, grpname);

    gnse->sans = gmx_sans_init(top, gnsf);

    /* Prepare reference frame */
    if (bPBC)
    {
        gpbc = gmx_rmpbc_init(&top->idef, ePBC, top->atoms.nr, box);
    }

    natoms = read_first_x(oenv, &status, fnTRX, &t, &xf, box);
    if (natoms != atoms->nr)
    {
        fprintf(stderr, "\nWARNING: number of atoms in tpx (%d) and trajectory (%d) do not match\n", natoms, atoms->nr);
    }
    /* realy do calc */
    nframes = 0;
    /* allocate structures for gnse->x,t,box for first frame */
    snew(gnse->x, nframes+1);
    snew(gnse->t, nframes+1);
    snew(gnse->box, nframes+1);
    /* copy data for x, t , box */
    copy_mat(box, gnse->box[nframes]);
    gnse->t[nframes] = t;
    snew(gnse->x[nframes], isize);
    for (i = 0; i < isize; i++)
    {
        copy_rvec(xf[index[i]], gnse->x[nframes][i]);
    }
    /* Read whole trajectory into memroy and allocate gnse structure */
    do
    {
        if (bPBC)
        {
            gmx_rmpbc(gpbc, atoms->nr, box, xf);
        }
        /* resize arrays */
        srenew(gnse->x, nframes+1);
        srenew(gnse->t, nframes+1);
        srenew(gnse->box, nframes+1);

        copy_mat(box, gnse->box[nframes]);
        gnse->t[nframes] = t;
        snew(gnse->x[nframes], isize);
        for (i = 0; i < isize; i++)
        {
            copy_rvec(xf[index[i]], gnse->x[nframes][i]);
        }
        nframes++;
    }
    while (read_next_x(oenv, status, &t, natoms, xf, box));
    close_trj(status);

    gnse->nrframes = nframes;

    /* allocate gnse->sq,gr structures */
    snew(gnse->sq, gnse->nrframes);
    snew(gnse->gr, gnse->nrframes);
    for (i = 0; i < gnse->nrframes; i++)
    {
        snew(gnse->gr[i], 1);
    }

    /*
     * now try to populate avereged over time cross histograms t,t+dt
     * j = dt, frame = i
     */
    pairs = (int)floor(0.5*gnse->nrframes*(gnse->nrframes+1));
    snew(gnse->dt, gnse->nrframes);
    fprintf(stderr, "Total numer of pairs = %10d\n", pairs);

    for (i = 0; i < gnse->nrframes; i++)
    {
        for (j = 0; j+i < gnse->nrframes; j++)
        {
            if (grc == NULL)
            {
                snew(grc, 1);
            }
            grc = calc_radial_distribution_histogram(gnse->sans, gnse->x[j], gnse->x[j+i], gnse->box[j], gnse->box[j+i], index, isize, binwidth, bMC, bNSE, mcover, seed);
            /* Copy common things */
            gnse->gr[i]->binwidth = grc->binwidth;
            /* also we should be make sure that there will be no buffer overruns */
            if (gnse->gr[i]->gr == NULL)
            {
                snew(gnse->gr[i]->gr, grc->grn);
                snew(gnse->gr[i]->r, grc->grn);
                gnse->gr[i]->grn = grc->grn;
            }
            else
            {
                if (grc->grn > gnse->gr[i]->grn)
                {
                    gnse->gr[i]->grn = grc->grn;
                    srenew(gnse->gr[i]->gr, grc->grn);
                    srenew(gnse->gr[i]->r, grc->grn);
                }
            }

            /* now we need to summ up gr's */
            for (k = 0; k < grc->grn; k++)
            {
                gnse->gr[i]->gr[k] += grc->gr[k];
                gnse->gr[i]->r[k]   = grc->r[k];
            }
            /* we can free grc */
            sfree(grc->gr);
            sfree(grc->r);
            sfree(grc);
            pairs--;
            gnse->dt[i] = gnse->t[j+i]-gnse->t[j];
            fprintf(stderr, "\rdt = %10.2f pairs left to compute %10d", gnse->dt[i], pairs);
        }
        normalize_probability(gnse->gr[i]->grn, gnse->gr[i]->gr);
        gnse->sq[i] = convert_histogram_to_intensity_curve(gnse->gr[i], start_q, end_q, q_step);
        if (opt2fn_null("-stq", NFILE, fnm))
        {
            snew(hdr, 25);
            snew(suffix, GMX_PATH_MAX);
            /* prepare header */
            sprintf(hdr, "S(q,dt), dt = %10.2f", gnse->dt[i]);
            /* prepare output filename */
            fnmdup = dup_tfn(NFILE, fnm);
            sprintf(suffix, "-t%010.2f", gnse->dt[i]);
            add_suffix_to_output_names(fnmdup, NFILE, suffix);
            fp = xvgropen(opt2fn_null("-stq", NFILE, fnmdup), hdr, "q, nm^-1", "S(q,t) arb.u.", oenv);
            for (j = 0; j < gnse->sq[i]->qn; j++)
            {
                fprintf(fp, "%10.6lf%10.6lf\n", gnse->sq[i]->q[j], gnse->sq[i]->s[j]);
            }
            done_filenms(NFILE, fnmdup);
            fclose(fp);
            sfree(hdr);
            sfree(fnmdup);
        }
        if (opt2fn_null("-grt", NFILE, fnm))
        {
            snew(hdr, 25);
            snew(suffix, GMX_PATH_MAX);
            /* prepare header */
            sprintf(hdr, "G(r,dt), dt = %10.2f", gnse->dt[i]);
            /* prepare output filename */
            fnmdup = dup_tfn(NFILE, fnm);
            sprintf(suffix, "-t%010.2f", gnse->dt[i]);
            add_suffix_to_output_names(fnmdup, NFILE, suffix);
            fp = xvgropen(opt2fn_null("-grt", NFILE, fnmdup), hdr, "r, nm", "G(r,dt)", oenv);
            for (j = 0; j < gnse->gr[i]->grn; j++)
            {
                fprintf(fp, "%10.6lf%10.6lf\n", gnse->gr[i]->r[j], gnse->gr[i]->gr[j]);
            }
            done_filenms(NFILE, fnmdup);
            fclose(fp);
            sfree(hdr);
            sfree(fnmdup);
        }
        /* Now we can free gnse->gr[i] */
        sfree(gnse->gr[i]->gr);
        sfree(gnse->gr[i]->r);
        sfree(gnse->gr[i]);

    }
    /* now we can free gnse->x gnse->box */
    for (i = 0; i < gnse->nrframes; i++)
    {
        sfree(gnse->x[i]);
    }
    sfree(gnse->box);
    sfree(gnse->t);
    /* Also we can clean index and top */
    sfree(index);
    sfree(grpname);
    done_sans(gnse->sans);
    done_nsf(gnsf);
    sfree(top); // done_top already done via done_sans


    fprintf(stderr, "\n");
    gnse->sqtn = gnse->sq[0]->qn;
    snew(gnse->sqt, gnse->sqtn);

    /* now we will gather s(q(t)) from s(q) spectrums */
    for (i = 0; i < gnse->sq[0]->qn; i++)
    {
        snew(gnse->sqt[i], 1);
        gnse->sqt[i]->q = gnse->sq[0]->q[i];
        snew(gnse->sqt[i]->s, gnse->nrframes);
        for (j = 0; j < gnse->nrframes; j++)
        {
            gnse->sqt[i]->s[j] = gnse->sq[j]->s[i];
        }
    }
    /* Now we can free gnse->sq */
    for (i = 0; i < gnse->nrframes; i++)
    {
        sfree(gnse->sq[i]->q);
        sfree(gnse->sq[i]->s);
        sfree(gnse->sq[i]);
    }

    /* actualy print data */
    for (i = 0; i < gnse->sqtn; i++)
    {
        snew(hdr, 25);
        snew(suffix, GMX_PATH_MAX);
        /* prepare header */
        sprintf(hdr, "S(q,dt) , q = %2.2lf", gnse->sqt[i]->q);
        /* prepare output filename */
        fnmdup = dup_tfn(NFILE, fnm);
        sprintf(suffix, "-q%2.2lf", gnse->sqt[i]->q);
        add_suffix_to_output_names(fnmdup, NFILE, suffix);
        fp = xvgropen(opt2fn("-sqt", NFILE, fnmdup), hdr, "dt", "S(q,dt)", oenv);
        for (j = 0; j < gnse->nrframes; j++)
        {
            if (bNORMALIZE)
            {
                fprintf(fp, "%10.6lf%10.6lf\n", gnse->dt[j], gnse->sqt[i]->s[j]/gnse->sqt[i]->s[0]);
            }
            else
            {
                fprintf(fp, "%10.6lf%10.6lf\n", gnse->dt[j], gnse->sqt[i]->s[j]);
            }

        }
        done_filenms(NFILE, fnmdup);
        fclose(fp);
        sfree(hdr);
        sfree(fnmdup);
    }

    /* Now we need to clean up rest structures */
    sfree(gnse->dt);
    for (i = 0; i < gnse->nrframes; i++)
    {
        sfree(gnse->sqt[i]->s);
        sfree(gnse->sqt[i]);
    }
    sfree(gnse->sqt);
    sfree(gnse);

    // please_cite("Shvetsov2012"); /* to be published */

    thanx(stderr);

    return 0;
}
