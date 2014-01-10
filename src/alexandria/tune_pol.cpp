/*
 * $Id: tune_pol.c,v 1.32 2009/05/20 06:25:20 spoel Exp $
 *
 *                This source code is part of
 *
 *                 G   R   O   M   A   C   S
 *
 *          GROningen MAchine for Chemical Simulations
 *
 *                        VERSION 4.0.99
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
 * Groningen Machine for Chemical Simulation
 */

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "maths.h"
#include "gromacs/fileio/futil.h"
#include "smalloc.h"
#include "string2.h"
#include "vec.h"
#include "statutil.h"
#include "copyrite.h"
#include "gstat.h"
#include "gmx_fatal.h"
#include "names.h"
#include "gromacs/linearalgebra/matrix.h"
#include "gromacs/legacyheaders/gmx_statistics.h"
#include "gromacs/legacyheaders/gmx_random.h"
#include "gromacs/legacyheaders/xvgr.h"
#include "molselect.hpp"
#include "poldata.hpp"
#include "poldata_xml.hpp"
#include "molprop.hpp"
#include "molprop_xml.hpp"
#include "molprop_tables.hpp"
#include "molprop_util.hpp"

int check_matrix(double **a, double *x, int nrow, int ncol, char **atype)
{
    int nrownew = nrow;
    for (int i = 0; (i < ncol); i++)
    {
        for (int j = i+1; (j < ncol); j++)
        {
            bool bSame = true;
            for (int k = 0; bSame && (k < nrow); k++)
            {
                bSame = (a[k][i] == a[k][j]);
            }
            if (bSame)
            {
                gmx_fatal(FARGS, "Columns %d (%s) and %d (%s) are linearly dependent",
                          i, atype[i], j, atype[j]);
            }
        }
    }
    //! Check diagonal
    if (nrow == ncol)
    {
        for (int i = 0; (i < ncol); i++)
        {
            if (a[i][i] == 0)
            {
                gmx_fatal(FARGS, "a[%d][%d] = 0. Atom type = %s", i, i, atype[i]);
            }
        }
    }
    return nrow;
    for (int i = 0; (i < nrownew); i++)
    {
        for (int j = i+1; (j < nrownew); j++)
        {
            bool bSame = true;
            for (int k = 0; bSame && (k < ncol); k++)
            {
                bSame = (a[i][k] == a[j][k]);
            }
            if (bSame)
            {
                fprintf(stderr, "Rows %d and %d are linearly dependent. Removing %d.\n",
                        i, j, j);
                if (j < nrownew-1)
                {
                    for (int k = 0; (k < ncol); k++)
                    {
                        a[j][k] = a[nrownew-1][k];
                    }
                    x[j] = x[nrownew-1];
                }
                nrownew--;
            }
        }
    }
    return nrownew;
}

static void dump_csv(gmx_poldata_t pd,
                     std::vector<alexandria::MolProp> &mp,
                     char **atype,
                     int  nusemol,
                     bool bUseMol[],
                     double x[],
                     int  ntest,
                     double **a,
                     double **at)
{
    FILE *csv = ffopen("out.csv", "w");

    fprintf(csv, "\"molecule\",\"formula\",");
    for (int i = 0; (i < ntest); i++)
    {
        fprintf(csv, "\"%d %s\",", i, atype[i]);
    }
    fprintf(csv, "\"polarizability\"\n");
    int nn = 0;
    int j  = 0;
    for (alexandria::MolPropIterator mpi = mp.begin(); (mpi < mp.end()); mpi++, j++)
    {
        if (bUseMol[j])
        {
            alexandria::MolecularCompositionIterator mci =
                mpi->SearchMolecularComposition((char *)"spoel");
            fprintf(csv, "\"%d %s\",\"%s\",",
                    nn, mpi->GetMolname().c_str(), mpi->GetFormula().c_str());
            int *count;
            snew(count, ntest);
            for (alexandria::AtomNumIterator ani = mci->BeginAtomNum();
                 ani < mci->EndAtomNum(); ani++)
            {
                const char *atomname = ani->GetAtom().c_str();
                const char *ptype    = gmx_poldata_atype_to_ptype(pd, atomname);
                int         i;
                for (i = 0; (i < ntest); i++)
                {
                    if (strcmp(ptype, atype[i]) == 0)
                    {
                        break;
                    }
                }
                if (i < ntest)
                {
                    count[i] += ani->GetNumber();
                }
                else
                {
                    gmx_fatal(FARGS, "Supported molecule %s has unsupported atom %s (ptype %s)", mpi->GetMolname().c_str(), atomname, ptype);
                }
            }
            for (int i = 0; (i < ntest); i++)
            {
                a[nn][i] = at[i][nn] = count[i];
                fprintf(csv, "%d,", count[i]);
            }
            sfree(count);
            fprintf(csv, "%.3f\n", x[nn]);
            nn++;
        }
    }
    fclose(csv);
    if (nusemol != nn)
    {
        gmx_fatal(FARGS, "Consistency error: nusemol = %d, nn = %d", nusemol, nn);
    }

}
     
static int decompose_frag(FILE *fplog,
                          const char *hisfn,
                          gmx_poldata_t pd,
                          std::vector<alexandria::MolProp> &mp,
                          gmx_bool bQM, char *lot, 
                          int mindata, gmx_molselect_t gms,
                          gmx_bool bZero, gmx_bool bForceFit,
                          int nBootStrap, real fractionBootStrap,
                          int seed,
                          output_env_t oenv)
{
    double  *x, *atx;
    double **a, **at, **ata, *fpp;
    double   pol, poltot, a0, da0, ax, chi2;
    char   **atype = NULL, *ptype;
    int      j, niter = 0, nusemol;
    int     *test = NULL, ntest[2], ntmax, row, cur = 0;
    bool    *bUseMol, *bUseAtype;
#define prev (1-cur)

    snew(bUseMol, mp.size()+1);
    snew(x, mp.size()+1);
    ntmax = gmx_poldata_get_nptypes(pd);
    snew(bUseAtype, ntmax);
    snew(test, ntmax);
    snew(atype, ntmax);
    ntest[prev] = ntest[cur] = ntmax;
    j           = 0;
    // Copy all atom types into array. Set usage array.
    while (1 == gmx_poldata_get_ptype(pd, &ptype, NULL, NULL, NULL, NULL))
    {
        atype[j]     = strdup(ptype);
        bUseAtype[j] = true;
        j++;
    }
    int iter      = 1;
    int nmol_orig = mp.size();
    do
    {
        if (NULL != fplog)
        {
            fprintf(fplog, "iter %d ntest %d\n", iter++, ntest[cur]);
        }
        cur     = prev;
        nusemol = 0;
        poltot  = 0;
        j       = 0;
        for (alexandria::MolPropIterator mpi = mp.begin(); (mpi < mp.end()); mpi++, j++)
        {
            iMolSelect ims  = gmx_molselect_status(gms, mpi->GetIupac().c_str());

            bool       bPol = mpi->GetProp(MPO_POLARIZABILITY,
                                           bQM ? iqmBoth : iqmExp,
                                           lot, NULL, NULL, &pol);
            alexandria::MolecularCompositionIterator mci =
                mpi->SearchMolecularComposition((char *)"spoel");

            bUseMol[j] = ((imsTrain == ims) && bPol && (pol > 0) &&
                          (mci != mpi->EndMolecularComposition()));

            /* Now check for the right composition */
            for (alexandria::AtomNumIterator ani = mci->BeginAtomNum();
                 (bUseMol[j] && (ani < mci->EndAtomNum())); ani++)
            {
                double      apol, spol;
                const char *atomname = ani->GetAtom().c_str();
                const char *ptype    = gmx_poldata_atype_to_ptype(pd, atomname);
                if (NULL == ptype)
                {
                    if (NULL != fplog)
                    {
                        fprintf(fplog, "No ptype for atom %s in molecule %s\n",
                                atomname, mpi->GetMolname().c_str());
                    }
                    bUseMol[j] = false;
                }
                else if (0 == gmx_poldata_get_ptype_pol(pd, ptype, &apol, &spol))
                {
                    bUseMol[j] = false;
                }
                else if (apol > 0)
                {
                    bUseMol[j] = false;
                    pol       -= apol*ani->GetNumber();
                }
            }
            if (NULL != debug)
            {
                fprintf(debug, "Mol: %s, IUPAC: %s, ims: %s, bPol: %s, pol: %g - %s\n",
                        mpi->GetMolname().c_str(),
                        mpi->GetIupac().c_str(),
                        ims_names[ims], bool_names[bPol], pol,
                        bUseMol[j] ? "Used" : "Not used");
            }

            if (bUseMol[j])
            {
                x[nusemol] = pol;
                poltot    += pol;
                nusemol++;
            }
            else
            {
                mpi = mp.erase(mpi);
            }
        }
        ntest[cur] = 0;
        int ntp = 0;
        while (1 == gmx_poldata_get_ptype(pd, &ptype, NULL, NULL, &pol, NULL))
        {
            if ((pol == 0) || bForceFit)
            {
                int j = 0;
                for (alexandria::MolPropIterator mpi = mp.begin(); (mpi < mp.end()); mpi++, j++)
                {
                    if (bUseMol[j])
                    {
                        alexandria::MolecularCompositionIterator mci =
                            mpi->SearchMolecularComposition((char *)"spoel");
                        for (alexandria::AtomNumIterator ani = mci->BeginAtomNum();
                             ani < mci->EndAtomNum(); ani++)
                        {
                            if (strcmp(ptype,
                                       gmx_poldata_atype_to_ptype(pd, ani->GetAtom().c_str())) == 0)
                            {
                                test[ntp] += ani->GetNumber();
                            }
                        }
                    }
                }
                bUseAtype[ntp] = (test[ntp] >= mindata);
                if (bUseAtype[ntp])
                {
                    ntest[cur]++;
                }
                else
                {
                    if (NULL != fplog)
                    {
                        fprintf(fplog, "Not enough polarizability data points (%d out of %d required) to optimize %s\n", 
                                test[ntp], mindata, ptype);
                    }
                }
            }
            else
            {
                if (NULL != fplog)
                {
                    fprintf(fplog, "Polarizability for %s is %g. Not optimizing this one.\n",
                            ptype, pol);
                }
            }
            ntp++;
        }
    }
    while (ntest[cur] < ntest[prev]);
    if (ntest[cur] == 0)
    {
        gmx_fatal(FARGS, "Nothing to optimize. Check your input");
    }
    if ((int)mp.size() < nmol_orig)
    {
        printf("Reduced number of molecules from %d to %d\n", nmol_orig, (int)mp.size());
    }

    //! Now condense array of atom types
    for (int i = 0; (i < ntmax); i++)
    {
        double apol, spol;
        if (bUseAtype[i] &&
            (gmx_poldata_get_ptype_pol(pd, atype[i], &apol, &spol) == 1))
        {
            bUseAtype[i] = (apol == 0);
        }
    }

    int i, nt = 0;
    for (i = 0; (i < ntmax); i++)
    {
        if (bUseAtype[i])
        {
            nt++;
        }
    }
    ntest[cur] = nt;
    for (i = 0; (i < ntmax); i++)
    {
        for (int j = i+1; !bUseAtype[i] && (j < ntmax); j++)
        {
            if (bUseAtype[j])
            {
                bUseAtype[i] = true;
                sfree(atype[i]);
                atype[i]     = atype[j];
                atype[j]     = NULL;
                bUseAtype[j] = false;
            }
        }
    }
    /* Now we know how many atomtyopes there are (ntest[cur]) and 
     * how many molecules (nusemol), now we fill the matrix a and the
     * transposed matrix at 
     */
    a      = alloc_matrix(nusemol, ntest[cur]);
    at     = alloc_matrix(ntest[cur], nusemol);
    if (NULL != fplog)
    {
        fprintf(fplog, "There are %d different atomtypes to optimize the polarizabilities\n",
                ntest[cur]);
        fprintf(fplog, "There are %d (experimental) reference polarizabilities.\n", nusemol);
    }
    // As a side effect this function fills a and at and probably x.
    dump_csv(pd, mp, atype, nusemol, bUseMol, x, ntest[cur], a, at);
    
    if (fplog)
    {
        for (int i = 0; (i < ntest[cur]); i++)
        {
            fprintf(fplog, "Optimizing polarizability for %s with %d copies\n",
                    atype[i], test[i]);
        }
    }

    // Check for linear dependencies
    nusemol = check_matrix(a, x, nusemol, ntest[cur], atype);
    
    // Now loop over the number of bootstrap loops
    gmx_stats_t *polstats;
    snew(polstats, ntest[cur]);
    for(int ii=0; (ii<ntest[cur]); ii++)
    {
        polstats[ii] = gmx_stats_init();
    }
    int nUseBootStrap = 1+floor(fractionBootStrap*nusemol);
    if (seed <= 0)
    {
        seed = gmx_rng_make_seed();
    }
    gmx_rng_t rng = gmx_rng_init(seed);
    for(int kk=0; (kk<nBootStrap); kk++)
    {
        fprintf(stderr, "\rBootStrap %d", 1+kk);
        ata    = alloc_matrix(ntest[cur], ntest[cur]);

        double **a_copy  = alloc_matrix(nUseBootStrap, ntest[cur]);
        double **at_copy = alloc_matrix(ntest[cur], nUseBootStrap);
        double  *x_copy;
        snew(x_copy, nUseBootStrap);
        for(int ii=0; (ii<nUseBootStrap); ii++)
        {
            unsigned int uu = gmx_rng_uniform_uint32(rng) % nUseBootStrap;
            for(int jj=0; (jj<ntest[cur]); jj++)
            {
                a_copy[ii][jj] = at_copy[jj][ii] =a[uu][jj];
                x_copy[ii] = x[uu];
            }
        }
        matrix_multiply(fplog, nUseBootStrap, ntest[cur], a_copy, at_copy, ata);
        (void) check_matrix(ata, x_copy, ntest[cur], ntest[cur], atype);

        if ((row = matrix_invert(fplog, ntest[cur], ata)) == 0)
        {
            snew(atx, ntest[cur]);
            snew(fpp, ntest[cur]);
            a0 = 0;
            do
            {
                for (int i = 0; (i < ntest[cur]); i++)
                {
                    atx[i] = 0;
                    for (j = 0; (j < nUseBootStrap); j++)
                    {
                        atx[i] += at_copy[i][j]*(x_copy[j]-a0);
                    }
                }
                for (int i = 0; (i < ntest[cur]); i++)
                {
                    fpp[i] = 0;
                    for (j = 0; (j < ntest[cur]); j++)
                    {
                        fpp[i] += ata[i][j]*atx[j];
                    }
                }
                da0  = 0;
                chi2 = 0;
                if (bZero)
                {
                    for (j = 0; (j < nUseBootStrap); j++)
                    {
                        ax = a0;
                        for (int i = 0; (i < ntest[cur]); i++)
                        {
                            ax += fpp[i]*a_copy[j][i];
                        }
                        da0  += (x_copy[j]-ax);
                        chi2 += sqr(x_copy[j]-ax);
                    }
                    da0 = da0 / nusemol;
                    a0 += da0;
                    niter++;
                    printf("iter: %d <pol> = %g, a0 = %g, chi2 = %g\n",
                           niter, poltot/nusemol, a0, chi2/nusemol);
                }
            }
            while (bZero && (fabs(da0) > 1e-5) && (niter < 1000));
            for (int i = 0; (i < ntest[cur]); i++)
            {
                gmx_stats_add_point_ydy(polstats[i], fpp[i], 0);
            }
        }
    }
    fprintf(stderr, "\n");
    FILE *xp = xvgropen(hisfn, "Polarizability distribution", "alpha (A\\S3\\N)", "", oenv);
    for (int i = 0; (i < ntest[cur]); i++)
    {
        int result1, result2;
        real aver, sigma;
        // Extract data from statistics
        if ((estatsOK == (result1 = gmx_stats_get_average(polstats[i], &aver))) &&
            (estatsOK == (result2 = gmx_stats_get_sigma(polstats[i], &sigma))))
        {
            gmx_poldata_set_ptype_polarizability(pd, atype[i], aver, sigma);
            fprintf(fplog, "%-5s  %8.3f +/- %.3f\n", atype[i], aver, sigma);
            int nbins = 1+sqrt(nBootStrap);
            real *my_x, *my_y;
            gmx_stats_make_histogram(polstats[i], 0, &nbins,
                                     ehistoY, 1, &my_x, &my_y);
            fprintf(xp, "@type xy\n");
            for(int ll=0; (ll<nbins); ll++)
            {
                fprintf(xp, "%.3f  %.3f\n", my_x[i], my_y[i]);
            }
            fprintf(xp, "&\n");
            sfree(my_x);
            sfree(my_y);
        }
        else
        {
            fprintf(stderr, "Could not determine polarizability for %s\n",
                    atype[i]);
        }
        gmx_stats_done(polstats[i]);
    }
    fclose(xp);
    sfree(polstats);
    if (bZero)
    {
        const char *null = (const char *)"0";
        gmx_poldata_add_ptype(pd, null, NULL, null, a0, 0);
    }
    sfree(bUseMol);
    sfree(fpp);
    
    return ntest[cur];
}

int alex_tune_pol(int argc, char *argv[])
{
    static const char               *desc[] =
    {
        "tune_pol optimes atomic polarizabilities that together build",
        "an additive model for polarization. The set of atomtypes used",
        "is determined by the input force field file ([TT]-di[tt] option). All",
        "atomtypes for which the polarizability is zero, and for which",
        "there is sufficient data in the input data file ([TT]-f[tt] option)",
        "will be used in the least-squares fit (done by matrix inversion).[PAR]"
        "Bootstrapping can be used to estimate the error in the resulting parameters",
        "and this will be stored in the resulting file such that programs",
        "using the parameters can estimate an error in the polarizability.[PAR]",
        "A selection of molecules into a training set and a test set (or ignore set)",
        "can be made using option [TT]-sel[tt]. The format of this file is:[BR]",
        "iupac|Train[BR]",
        "iupac|Test[BR]",
        "iupac|Ignore[BR]",
        "and you should ideally have a line for each molecule in the molecule database",
        "([TT]-f[tt] option). Missing molecules will be ignored."
    };
    t_filenm                         fnm[] =
    {
        { efDAT, "-f",  "data",      ffRDMULT },
        { efDAT, "-o",  "allmols",   ffOPTWR  },
        { efDAT, "-di", "gentop",    ffOPTRD  },
        { efDAT, "-do", "tune_pol",  ffWRITE  },
        { efDAT, "-sel","molselect", ffREAD   },
        { efLOG, "-g",  "tune_pol",  ffWRITE  },
        { efXVG, "-his",  "polhisto",  ffWRITE  }
    };
    int                              NFILE    = (sizeof(fnm)/sizeof(fnm[0]));
    static char                     *sort[]   = { NULL, (char *)"molname", (char *)"formula", (char *)"composition", NULL };
    static gmx_bool                  bQM      = FALSE;
    static int                       mindata  = 1;
    static int                       nBootStrap = 1;
    static int                       seed;
    static real                      fractionBootStrap = 1;
    static char                     *lot      = (char *)"B3LYP/aug-cc-pVTZ";
    static real                      sigma    = 0;
    static gmx_bool                  bZero    = FALSE, bForceFit = FALSE, bCompress = TRUE;
    t_pargs                          pa[]     =
    {
        { "-sort",   FALSE, etENUM, {sort},
          "Key to sort the final data file on." },
        { "-qm",     FALSE, etBOOL, {&bQM},
          "Use QM data for optimizing the empirical polarizabilities as well." },
        { "-lot",    FALSE, etSTR,  {&lot},
          "Use this method and level of theory" },
        { "-mindata", FALSE, etINT, {&mindata},
          "Minimum number of data points to optimize a polarizability value" },
        { "-sigma",  FALSE, etREAL, {&sigma},
          "Assumed standard deviation (relative) in the experimental data points" },
        { "-zero",    FALSE, etBOOL, {&bZero},
          "Use a zero term (like in Bosque and Sales)" },
        { "-force",   FALSE, etBOOL, {&bForceFit},
          "Reset all polarizablities to zero in order to re-optimize based on a gentop.dat file with previous values" },
        { "-nBootStrap", FALSE, etINT, {&nBootStrap},
          "Number of trials for bootstrapping" },
        { "-fractionBootStrap", FALSE, etREAL, {&fractionBootStrap},
          "Fraction of data points to use in each trial for bootstrapping" },
        { "-seed", FALSE, etINT, {&seed},
          "Seed for random numbers in bootstrapping. If <= 0 a seed will be generated." },
        { "-compress", FALSE, etBOOL, {&bCompress},
          "Compress output XML files" }
    };
    int                              i, nalexandria_atypes;
    char                           **fns;
    int                              nfiles;
    std::vector<alexandria::MolProp> mp;
    alexandria::MolPropIterator      mpi;
    MolPropSortAlgorithm             mpsa;

    gmx_atomprop_t                   ap;
    gmx_poldata_t                    pd;
    output_env_t                     oenv;
    gmx_molselect_t                  gms;
    int                              npa = sizeof(pa)/sizeof(pa[0]);

    if (!parse_common_args(&argc, argv, PCA_NOEXIT_ON_ARGS, NFILE, fnm,
                           npa, pa, sizeof(desc)/sizeof(desc[0]), desc,
                           0, NULL, &oenv))
    {
        return 0;
    }
    if ((fractionBootStrap < 0) || (fractionBootStrap > 1))
    {
        gmx_fatal(FARGS, "fractionBootStrap should be in [0..1]");
    }
    if (nBootStrap <= 0) 
    {
        gmx_fatal(FARGS, "nBootStrap should be >= 1");
    }
    ap = gmx_atomprop_init();
    if ((pd = gmx_poldata_read(opt2fn_null("-di", NFILE, fnm), ap)) == NULL)
    {
        gmx_fatal(FARGS, "Can not read the force field information. File missing or incorrect.");
    }
    nfiles = opt2fns(&fns, "-f", NFILE, fnm);
    merge_xml(nfiles, fns, mp, NULL, NULL, (char *)"double_dip.dat", ap, pd, TRUE);
    for (mpi = mp.begin(); (mpi < mp.end()); mpi++)
    {
        mpi->CheckConsistency();
    }
    gms               = gmx_molselect_init(opt2fn("-sel", NFILE, fnm));
    FILE *fplog       = opt2FILE("-g", NFILE, fnm, "w");
    nalexandria_atypes = decompose_frag(fplog, opt2fn("-his", NFILE, fnm),
                                        pd, mp, bQM, lot, mindata,
                                        gms, bZero, bForceFit,
                                        nBootStrap, fractionBootStrap, seed,
                                        oenv);
    fprintf(fplog, "There are %d alexandria atom types\n", nalexandria_atypes);

    const char *pdout = opt2fn("-do", NFILE, fnm);
    fprintf(fplog, "Now writing force field file %s\n", pdout);
    gmx_poldata_write(pdout, pd, bCompress);
    
    ffclose(fplog);
    
    const char *mpfn = opt2fn_null("-o", NFILE, fnm);
    if (NULL != mpfn) 
    {
        mpsa = MPSA_NR;
        if (opt2parg_bSet("-sort", npa, pa))
        {
            for (i = 0; (i < MPSA_NR); i++)
            {
                if (strcasecmp(sort[0], sort[i+1]) == 0)
                {
                    mpsa = (MolPropSortAlgorithm) i;
                    break;
                }
            }
        }
        if (mpsa != MPSA_NR)
        {
            MolPropSort(mp, mpsa, ap, NULL);
        }
        
        MolPropWrite(mpfn, mp, bCompress);
    }
    
    return 0;
}
