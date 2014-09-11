/*
 * This source file is part of the Aleandria project.
 *
 * Copyright (C) 2014 David van der Spoel and Paul J. van Maaren
 *
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.
 */
/*! \internal \brief
 * Implements part of the alexandria program.
 * \author David van der Spoel <david.vanderspoel@icm.uu.se>
 */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "gromacs/utility/futil.h"
#include "gromacs/utility/smalloc.h"
#include "gromacs/utility/cstringutil.h"
#include "gromacs/math/vec.h"
#include "gromacs/legacyheaders/copyrite.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/legacyheaders/names.h"
#include "gromacs/commandline/pargs.h"
#include "gromacs/linearalgebra/matrix.h"
#include "gromacs/legacyheaders/gmx_statistics.h"
#include "gromacs/random/random.h"
#include "gromacs/fileio/xvgr.h"
#include "molselect.h"
#include "poldata.h"
#include "poldata_xml.h"
#include "molprop.h"
#include "molprop_xml.h"
#include "molprop_tables.h"
#include "molprop_util.h"
#include "composition.h"

bool check_matrix(double **a, double *x, int nrow, int ncol, char **atype)
{
    int nrownew = nrow;

    //printf("check_matrix called with nrow = %d ncol = %d\n", nrow, ncol);
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
                return false;
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
                return false;
                gmx_fatal(FARGS, "a[%d][%d] = 0. Atom type = %s", i, i, atype[i]);
            }
        }
    }
    return true;
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

static void dump_csv(gmx_poldata_t                     pd,
                     std::vector<alexandria::MolProp> &mp,
                     char                            **atype,
                     int                               nusemol,
                     double                            x[],
                     int                               ntest,
                     double                          **a,
                     double                          **at)
{
    alexandria::CompositionSpecs cs;
    FILE *csv = gmx_ffopen("out.csv", "w");

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
        alexandria::MolecularCompositionIterator mci =
            mpi->SearchMolecularComposition(cs.searchCS(alexandria::iCalexandria)->name());
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
                gmx_fatal(FARGS, "Supported molecule %s has unsupported atom %s (ptype %s)",
                          mpi->GetMolname().c_str(), atomname, ptype);
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
    double                      *x, *atx;
    double                     **a, **at, **ata, *fpp;
    double                       pol, poltot, a0, da0, ax, chi2;
    char                       **atype = NULL, *ptype;
    int                          j, niter = 0, nusemol;
    int                         *test = NULL, ntest[2], ntmax, row, cur = 0;
    bool                         bUseMol, *bUseAtype;
#define prev (1-cur)
    alexandria::CompositionSpecs cs;
    const char                  *alex = cs.searchCS(alexandria::iCalexandria)->name();

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
    // Check whether molecules have all supported atom types, remove the molecule otherwise
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
        for (alexandria::MolPropIterator mpi = mp.begin(); (mpi < mp.end()); j++)
        {
            iMolSelect ims  = gmx_molselect_status(gms, mpi->GetIupac().c_str());

            pol             = 0;
            bool       bPol = mpi->GetProp(MPO_POLARIZABILITY,
                                           bQM ? iqmBoth : iqmExp,
                                           lot, NULL, NULL, &pol);
            alexandria::MolecularCompositionIterator mci =
                mpi->SearchMolecularComposition(alex);

            bUseMol = ((imsTrain == ims) && bPol && (pol > 0) &&
                       (mci != mpi->EndMolecularComposition()));

            /* Now check whether all atoms have supported ptypes */
            for (alexandria::AtomNumIterator ani = mci->BeginAtomNum();
                 (bUseMol && (ani < mci->EndAtomNum())); ani++)
            {
                double      apol     = 0, spol = 0;
                const char *atomname = ani->GetAtom().c_str();
                const char *ptype    = gmx_poldata_atype_to_ptype(pd, atomname);
                if (NULL == ptype)
                {
                    if (NULL != fplog)
                    {
                        fprintf(fplog, "No ptype for atom %s in molecule %s\n",
                                atomname, mpi->GetMolname().c_str());
                    }
                    bUseMol = false;
                }
                else if (0 == gmx_poldata_get_ptype_pol(pd, ptype, &apol, &spol))
                {
                    /* No polarizability found for this one, seems unnecessary
                     * as we first lookup the polarizability ptype */
                    bUseMol = false;
                    fprintf(stderr, "Dit mag nooit gebeuren %s %d\n", __FILE__, __LINE__);
                }
                else if (apol > 0)
                {
                    bUseMol = false;
                }
            }
            if (NULL != debug)
            {
                fprintf(debug, "Mol: %s, IUPAC: %s, ims: %s, bPol: %s, pol: %g - %s\n",
                        mpi->GetMolname().c_str(),
                        mpi->GetIupac().c_str(),
                        ims_names[ims], bool_names[bPol], pol,
                        bUseMol ? "Used" : "Not used");
            }

            if (bUseMol)
            {
                x[nusemol] = pol;
                poltot    += pol;
                nusemol++;
                mpi++;
            }
            else
            {
                fprintf(fplog, "Removing %s. bPol = %s pol = %g composition found = %s\n",
                        mpi->GetMolname().c_str(), bool_names[bPol], pol,
                        (mci == mpi->EndMolecularComposition()) ? "true" : "false");
                mpi = mp.erase(mpi);
            }
        }
        // Now we have checked all molecules, now check whether the ptypes still have support
        // in experimental/QM polarizabilities
        ntest[cur] = 0;
        int ntp = 0;
        while (1 == gmx_poldata_get_ptype(pd, &ptype, NULL, NULL, &pol, NULL))
        {
            if ((pol == 0) || bForceFit)
            {
                int j = 0;
                for (alexandria::MolPropIterator mpi = mp.begin(); (mpi < mp.end()); mpi++, j++)
                {
                    alexandria::MolecularCompositionIterator mci =
                        mpi->SearchMolecularComposition(alex);
                    for (alexandria::AtomNumIterator ani = mci->BeginAtomNum();
                         ani < mci->EndAtomNum(); ani++)
                    {
                        if (strcmp(ptype,
                                   gmx_poldata_atype_to_ptype(pd, ani->GetAtom().c_str())) == 0)
                        {
                            test[ntp]++;
                            break;
                        }
                    }
                }
                printf("iter %2d ptype %s test[%d] = %d\n", iter, ptype, ntp, test[ntp]);
                bUseAtype[ntp] = (test[ntp] >= mindata);
                if (bUseAtype[ntp])
                {
                    ntest[cur]++;
                    ntp++;
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
        }
    }
    while (ntest[cur] < ntest[prev]);

    if ((int)mp.size() < nmol_orig)
    {
        printf("Reduced number of molecules from %d to %d\n", nmol_orig, (int)mp.size());
    }
    if (ntest[cur] == 0)
    {
        printf("No polarization types to optimize. Check your input.\n");
        return 0;
    }
    else if (ntest[cur] < ntmax)
    {
        printf("Reduced number of atomtypes from %d to %d\n", ntmax, ntest[cur]);
    }
    // Compact the arrays
    for (int i = 0; (i < ntmax); i++)
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
    dump_csv(pd, mp, atype, nusemol, x, ntest[cur], a, at);

    if (fplog)
    {
        for (int i = 0; (i < ntest[cur]); i++)
        {
            fprintf(fplog, "Optimizing polarizability for %s with %d copies\n",
                    atype[i], test[i]);
        }
    }

    // Check for linear dependencies
    if (!check_matrix(a, x, nusemol, ntest[cur], atype))
    {
        fprintf(stderr, "Matrix is linearly dependent. Sorry.\n");
    }

    // Now loop over the number of bootstrap loops
    gmx_stats_t *polstats;
    snew(polstats, ntest[cur]);
    for (int ii = 0; (ii < ntest[cur]); ii++)
    {
        polstats[ii] = gmx_stats_init();
    }
    int nUseBootStrap = std::min(nusemol, (int)(1+floor(fractionBootStrap*nusemol)));
    if (seed <= 0)
    {
        seed = gmx_rng_make_seed();
    }
    gmx_rng_t rng = gmx_rng_init(seed);
    for (int kk = 0; (kk < nBootStrap); kk++)
    {
        fprintf(stderr, "\rBootStrap %d", 1+kk);
        ata    = alloc_matrix(ntest[cur], ntest[cur]);

        double **a_copy  = alloc_matrix(nUseBootStrap, ntest[cur]);
        double **at_copy = alloc_matrix(ntest[cur], nUseBootStrap);
        double  *x_copy;
        snew(x_copy, nUseBootStrap);
        for (int ii = 0; (ii < nUseBootStrap); ii++)
        {
            // Pick random molecule uu out of stack
            unsigned int uu = gmx_rng_uniform_uint32(rng) % nUseBootStrap;

            for (int jj = 0; (jj < ntest[cur]); jj++)
            {
                // Make row ii equal to uu in the original matrix
                a_copy[ii][jj] = at_copy[jj][ii] = a[uu][jj];
                // Make experimenal (QM) value equal to the original
                x_copy[ii] = x[uu];
            }
        }
        matrix_multiply(debug, nUseBootStrap, ntest[cur], a_copy, at_copy, ata);
        if (check_matrix(ata, x_copy, ntest[cur], ntest[cur], atype) &&
            ((row = matrix_invert(debug, ntest[cur], ata)) == 0))
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
        free_matrix(a_copy);
        free_matrix(at_copy);
        free_matrix(ata);
        sfree(x_copy);
    }
    fprintf(stderr, "\n");
    FILE *xp = xvgropen(hisfn, "Polarizability distribution", "alpha (A\\S3\\N)", "", oenv);
    xvgr_legend(xp, ntest[cur], (const char **)atype, oenv);

    for (int i = 0; (i < ntest[cur]); i++)
    {
        int  result1, result2;
        real aver, sigma;
        // Extract data from statistics
        if ((estatsOK == (result1 = gmx_stats_get_average(polstats[i], &aver))) &&
            (estatsOK == (result2 = gmx_stats_get_sigma(polstats[i], &sigma))))
        {
            gmx_poldata_set_ptype_polarizability(pd, atype[i], aver, sigma);
            fprintf(fplog, "%-5s  %8.3f +/- %.3f\n", atype[i], aver, sigma);
            int   nbins = 1+sqrt(nBootStrap);
            real *my_x, *my_y;
            gmx_stats_make_histogram(polstats[i], 0, &nbins,
                                     ehistoY, 1, &my_x, &my_y);
            fprintf(xp, "@type xy\n");
            for (int ll = 0; (ll < nbins); ll++)
            {
                fprintf(xp, "%.3f  %.3f\n", my_x[ll], my_y[ll]);
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
    sfree(fpp);
    free_matrix(a);
    free_matrix(at);

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
        "([TT]-f[tt] option). Missing molecules will be ignored.[PAR]",
        "tune_pol produces a table and a figure to describe the data."
    };
    t_filenm                         fnm[] =
    {
        { efDAT, "-f",  "data",      ffRDMULT },
        { efDAT, "-o",  "allmols",   ffOPTWR  },
        { efDAT, "-di", "gentop",    ffOPTRD  },
        { efDAT, "-do", "tune_pol",  ffWRITE  },
        { efDAT, "-sel", "molselect", ffREAD   },
        { efLOG, "-g",  "tune_pol",  ffWRITE  },
        { efTEX, "-atype",  "atomtypes", ffWRITE  },
        { efXVG, "-his",  "polhisto",  ffWRITE  }
    };
    int                              NFILE       = (sizeof(fnm)/sizeof(fnm[0]));
    static char                     *sort[]      = { NULL, (char *)"molname", (char *)"formula", (char *)"composition", NULL };
    static char                     *exp_type    = (char *)"Polarizability";
    static gmx_bool                  bQM         = FALSE;
    static int                       mindata     = 1;
    static int                       nBootStrap  = 1;
    static int                       seed;
    static real                      fractionBootStrap = 1;
    static char                     *lot               = (char *)"B3LYP/aug-cc-pVTZ";
    static real                      sigma             = 0;
    static gmx_bool                  bZero             = FALSE, bForceFit = FALSE, bCompress = TRUE;
    t_pargs                          pa[]              =
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
        { "-exp_type", FALSE, etSTR, {&exp_type},
          "[HIDDEN]The experimental property type that all the stuff above has to be compared to." },
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

    const char *atype = opt2fn("-atype", NFILE, fnm);
    fprintf(fplog, "Now writing LaTeX description of force field to %s\n", atype);
    FILE       *tp = fopen(atype, "w");
    gmx_molprop_atomtype_table(tp, true, pd, mp, lot, exp_type);
    fclose(tp);

    gmx_ffclose(fplog);

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
