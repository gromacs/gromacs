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
#include "gmxpre.h"
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

class pType 
{
private:
    std::string name_;
    bool        bUse_;
    int         nCopies_;
    gmx_stats_t polstats;
public:
    pType(std::string name, const bool bUse, const int nCopies);
    ~pType();
    bool bUse() { return bUse_; }
    void setUse(bool bUse) { bUse_ = bUse; }
    void checkUse(int mindata) { bUse_ = (nCopies_ > mindata); }
    std::string name() { return name_; }
    int nCopies() { return nCopies_; }
    void resetCopies() { nCopies_ = 0; }
    void incCopies() { nCopies_++; }
    gmx_stats_t stats() { return polstats; }
};

pType::pType(std::string name, const bool bUse, const int nCopies) 
{
    name_    = name; 
    bUse_    = bUse; 
    nCopies_ = nCopies; 
    polstats = gmx_stats_init();
}

pType::~pType()
{
    gmx_stats_done(polstats);
}

bool check_matrix(double **a, double *x, unsigned int nrow,
                  std::vector<pType> &ptypes)
{
    int nrownew = nrow;

    //printf("check_matrix called with nrow = %d ncol = %d\n", nrow, ncol);
    for (unsigned int i = 0; (i < ptypes.size()); i++)
    {
        for (unsigned int j = i+1; (j < ptypes.size()); j++)
        {
            bool bSame = true;
            for (unsigned int k = 0; bSame && (k < nrow); k++)
            {
                bSame = (a[k][i] == a[k][j]);
            }
            if (bSame)
            {
                return false;
                gmx_fatal(FARGS, "Columns %d (%s) and %d (%s) are linearly dependent",
                          i, ptypes[i].name().c_str(), j, ptypes[j].name().c_str());
            }
        }
    }
    //! Check diagonal
    if (nrow == ptypes.size())
    {
        for (unsigned int i = 0; (i <ptypes.size()); i++)
        {
            if (a[i][i] == 0)
            {
                return false;
                gmx_fatal(FARGS, "a[%d][%d] = 0. Pol type = %s",
                          i, i, ptypes[i].name().c_str());
            }
        }
    }
    return true;
    for (int i = 0; (i < nrownew); i++)
    {
        for (int j = i+1; (j < nrownew); j++)
        {
            bool bSame = true;
            for (unsigned int k = 0; bSame && (k < ptypes.size()); k++)
            {
                bSame = (a[i][k] == a[j][k]);
            }
            if (bSame)
            {
                fprintf(stderr, "Rows %d and %d are linearly dependent. Removing %d.\n",
                        i, j, j);
                if (j < nrownew-1)
                {
                    for (unsigned int k = 0; (k < ptypes.size()); k++)
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

static bool bZeroPol(const char *ptype, const char *zeropol)
{
    if (NULL != zeropol)
    {
        if (NULL != strstr(zeropol, ptype))
        {
            return true;
        }
    }
    return false;
}

static void dump_csv(gmx_poldata_t                     pd,
                     std::vector<alexandria::MolProp> &mp,
                     gmx_molselect_t                   gms,
                     std::vector<pType>               &ptypes,
                     int                               nusemol,
                     double                            x[],
                     double                          **a,
                     double                          **at)
{
    alexandria::CompositionSpecs cs;
    FILE *csv = gmx_ffopen("out.csv", "w");

    fprintf(csv, "\"molecule\",\"formula\",");
    for (unsigned int i = 0; (i < ptypes.size()); i++)
    {
        fprintf(csv, "\"%u %s\",", i, ptypes[i].name().c_str());
    }
    fprintf(csv, "\"polarizability\"\n");
    int nn = 0;
    int j  = 0;
    for (alexandria::MolPropIterator mpi = mp.begin(); (mpi < mp.end()); mpi++, j++)
    {
        iMolSelect ims  = gmx_molselect_status(gms, mpi->GetIupac().c_str());
        
        if (imsTrain == ims)
        {
            alexandria::MolecularCompositionIterator mci =
                mpi->SearchMolecularComposition(cs.searchCS(alexandria::iCalexandria)->name());
            fprintf(csv, "\"%d %s\",\"%s\",",
                    nn, mpi->GetMolname().c_str(), mpi->GetFormula().c_str());
            int *count;
            snew(count, ptypes.size());
            for (alexandria::AtomNumIterator ani = mci->BeginAtomNum();
                 ani < mci->EndAtomNum(); ani++)
            {
                const char *atomname = ani->GetAtom().c_str();
                const char *ptype    = gmx_poldata_atype_to_ptype(pd, atomname);
                unsigned int i;
                for (i = 0; (i < ptypes.size()); i++)
                {
                    if (strcmp(ptype, ptypes[i].name().c_str()) == 0)
                    {
                        break;
                    }
                }
                if (i < ptypes.size())
                {
                    count[i] += ani->GetNumber();
                }
                else if (NULL != debug)
                {
                    fprintf(debug, "Supported molecule %s has unsupported or zeropol atom %s (ptype %s)",
                            mpi->GetMolname().c_str(), atomname, ptype);
                }
            }
            for (unsigned int i = 0; (i < ptypes.size()); i++)
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
                          const char *zeropol,
                          output_env_t oenv)
{
    double             *x, *atx, *fpp;
    double            **a, **at, **ata;
    double              pol, sig_pol, poltot, a0, da0, ax, chi2;
    std::vector<pType>  ptypes;
    int                 j, niter = 0;
    int                 row;
    unsigned int        ptsize, nusemol;
    alexandria::CompositionSpecs cs;
    const char                  *alex = cs.searchCS(alexandria::iCalexandria)->name();

    snew(x, mp.size()+1);
    // Copy all atom types into array. Set usage array.
    {
        char *ptype;
        while (1 == gmx_poldata_get_ptype(pd, &ptype, NULL, NULL, NULL, NULL))
        {
            if (!bZeroPol(ptype, zeropol))
            {
                ptypes.push_back(pType(ptype, false, 0));
            }
        }
    }
    int iter      = 1;
    int nmol_orig = mp.size();
    // Check whether molecules have all supported atom types, remove the molecule otherwise
    do
    {
        ptsize = ptypes.size();
        
        if (NULL != fplog)
        {
            fprintf(fplog, "iter %d %d ptypes left\n", iter++, ptsize);
        }
        nusemol = 0;
        poltot  = 0;
        for (alexandria::MolPropIterator mpi = mp.begin(); (mpi < mp.end()); )
        {
            iMolSelect ims  = gmx_molselect_status(gms, mpi->GetIupac().c_str());

            pol             = 0;
            bool       bPol = mpi->GetProp(MPO_POLARIZABILITY,
                                           bQM ? iqmBoth : iqmExp,
                                           lot, NULL, NULL, &pol, NULL);
            alexandria::MolecularCompositionIterator mci =
                mpi->SearchMolecularComposition(alex);
            
            bool bHaveComposition = mci != mpi->EndMolecularComposition();
            bool bUseMol = ((imsTrain == ims) && bPol && (pol > 0) &&
                            bHaveComposition);
            if (NULL != fplog)
            {
                fprintf(fplog, "%s pol %g Use: %s\n", mpi->GetMolname().c_str(), pol, bool_names[bUseMol]);
            }

            /* Now check whether all atoms have supported ptypes */
            unsigned int npolarizable = 0;
            for (alexandria::AtomNumIterator ani = mci->BeginAtomNum();
                 (bUseMol && (ani < mci->EndAtomNum())); ++ani)
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
                else if (!bZeroPol(ptype, zeropol))
                {
                    npolarizable++;
                }
#ifdef OLD
                else if (0 == gmx_poldata_get_ptype_pol(pd, ptype, &apol, &spol))
                {
                    /* No polarizability found for this one, seems unnecessary
                     * as we first lookup the polarizability ptype */
                    bUseMol = false;
                    fprintf(stderr, "Dit mag nooit gebeuren %s %d\n", __FILE__, __LINE__);
                }
                else if (apol > 0)
                {
                    if (NULL != fplog)
                    {
                        fprintf(fplog, "There already is a polarizbility %g for atom %s.\n", 
                                apol, atomname);
                    }
                    /* bUseMol = false; */
                }
#endif
            }
            if (NULL != fplog)
            {
                fprintf(fplog, "Mol: %s, IUPAC: %s, ims: %s, bPol: %s, pol: %g - %s\n",
                        mpi->GetMolname().c_str(),
                        mpi->GetIupac().c_str(),
                        ims_names[ims], bool_names[bPol], pol,
                        bUseMol ? "Used" : "Not used");
            }

            if (bUseMol && (npolarizable > 0))
            {
                x[nusemol] = pol;
                poltot    += pol;
                nusemol++;
                mpi++;
            }
            else
            {
                fprintf(fplog, "Removing %s. bPol = %s pol = %g composition found = %s npolarizable = %d\n",
                        mpi->GetMolname().c_str(), bool_names[bPol], pol,
                        (mci == mpi->EndMolecularComposition()) ? "true" : "false",
                        npolarizable);
                mpi = mp.erase(mpi);
            }
        }
        // Now we have checked all molecules, now check whether the
        // ptypes still have support
        // in experimental/QM polarizabilities
        for(std::vector<pType>::iterator pi = ptypes.begin();
            (pi < ptypes.end()); )
        {
            pi->resetCopies();
            
            if ((1 == gmx_poldata_get_ptype_pol(pd, pi->name().c_str(),
                                               &pol, &sig_pol)) &&
                ((pol == 0) || bForceFit))
            {
                for (alexandria::MolPropIterator mpi = mp.begin(); (mpi < mp.end()); mpi++)
                {
                    alexandria::MolecularCompositionIterator mci =
                        mpi->SearchMolecularComposition(alex);
                    for (alexandria::AtomNumIterator ani = mci->BeginAtomNum();
                         ani < mci->EndAtomNum(); ani++)
                    {
                        if (strcmp(pi->name().c_str(),
                                   gmx_poldata_atype_to_ptype(pd, ani->GetAtom().c_str())) == 0)
                        {
                            pi->incCopies();
                            break;
                        }
                    }
                }
            }
            else
            {
                if (NULL != fplog)
                {
                    fprintf(fplog, "Polarizability for %s is %g. Not optimizing this one.\n",
                            pi->name().c_str(), pol);
                }
            }
            pi->checkUse(mindata);
            if (!pi->bUse())
            {
                pi = ptypes.erase(pi);
            }
            else
            {
                if (NULL != debug)
                {
                    fprintf(debug, "Will optimize polarizability for ptype %s with %d copies\n", 
                            pi->name().c_str(),
                            pi->nCopies());
                }
                pi++;
            }
        }
    }
    while (ptypes.size() < ptsize);

    if ((int)mp.size() < nmol_orig)
    {
        printf("Reduced number of molecules from %d to %d\n", nmol_orig, (int)mp.size());
    }
    if (ptypes.size() == 0)
    {
        printf("No polarization types to optimize, but before compaction %d.\n",
               ptsize);
        return 0;
    }
    /* Now we know how many atomtypes there are (ptypes.size()) and
     * how many molecules (nusemol), now we fill the matrix a and the
     * transposed matrix at
     */
    a      = alloc_matrix(nusemol, ptypes.size());
    at     = alloc_matrix(ptypes.size(), nusemol);
    if (NULL != fplog)
    {
        fprintf(fplog, "There are %d different atomtypes to optimize the polarizabilities\n",
                (int)ptypes.size());
        fprintf(fplog, "There are %d (experimental) reference polarizabilities.\n", nusemol);
        fflush(fplog);
    }
    // As a side effect this function fills a and at and probably x.
    dump_csv(pd, mp, gms, ptypes, nusemol, x, a, at);

    if (fplog)
    {
        for (std::vector<pType>::iterator pi = ptypes.begin();
             pi < ptypes.end(); ++pi)
        {
            fprintf(fplog, "Optimizing polarizability for %s with %d copies\n",
                    pi->name().c_str(), pi->nCopies());
        }
    }

    // Check for linear dependencies
    if (!check_matrix(a, x, nusemol, ptypes))
    {
        fprintf(stderr, "Matrix is linearly dependent. Sorry.\n");
    }

    // Now loop over the number of bootstrap loops
    int nUseBootStrap = std::min(nusemol, (unsigned int)(1+floor(fractionBootStrap*nusemol)));
    if (seed <= 0)
    {
        seed = gmx_rng_make_seed();
    }
    gmx_rng_t rng = gmx_rng_init(seed);
    for (int kk = 0; (kk < nBootStrap); kk++)
    {
        fprintf(stderr, "\rBootStrap %d", 1+kk);
        ata    = alloc_matrix(ptypes.size(), ptypes.size());

        double **a_copy  = alloc_matrix(nUseBootStrap, ptypes.size());
        double **at_copy = alloc_matrix(ptypes.size(), nUseBootStrap);
        double  *x_copy;
        snew(x_copy, nUseBootStrap);
        for (int ii = 0; (ii < nUseBootStrap); ii++)
        {
            // Pick random molecule uu out of stack
            unsigned int uu = gmx_rng_uniform_uint32(rng) % nUseBootStrap;

            for (unsigned int jj = 0; (jj < ptypes.size()); jj++)
            {
                // Make row ii equal to uu in the original matrix
                a_copy[ii][jj] = at_copy[jj][ii] = a[uu][jj];
                // Make experimenal (QM) value equal to the original
                x_copy[ii] = x[uu];
            }
        }
        matrix_multiply(debug, nUseBootStrap, ptypes.size(),
                        a_copy, at_copy, ata);
        if (check_matrix(ata, x_copy, ptypes.size(), ptypes) &&
            ((row = matrix_invert(debug, ptypes.size(), ata)) == 0))
        {
            snew(atx, ptypes.size());
            snew(fpp, ptypes.size());
            a0 = 0;
            do
            {
                for (unsigned int i = 0; (i < ptypes.size()); i++)
                {
                    atx[i] = 0;
                    for (j = 0; (j < nUseBootStrap); j++)
                    {
                        atx[i] += at_copy[i][j]*(x_copy[j]-a0);
                    }
                }
                for (unsigned int i = 0; (i < ptypes.size()); i++)
                {
                    fpp[i] = 0;
                    for (unsigned int j = 0; (j < ptypes.size()); j++)
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
                        for (unsigned int i = 0; (i < ptypes.size()); i++)
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
            for (unsigned int i = 0; (i < ptypes.size()); i++)
            {
                gmx_stats_add_point_ydy(ptypes[i].stats(), fpp[i], 0);
            }
        }
        free_matrix(a_copy);
        free_matrix(at_copy);
        free_matrix(ata);
        sfree(x_copy);
    }
    fprintf(stderr, "\n");
    FILE *xp = xvgropen(hisfn, "Polarizability distribution", "alpha (A\\S3\\N)", "", oenv);
    //xvgr_legend(xp, ptypes.size(), (const char **)ptypes, oenv);

    for (unsigned int i = 0; (i < ptypes.size()); i++)
    {
        int  result1, result2;
        real aver, sigma;
        // Extract data from statistics
        if ((estatsOK == (result1 = gmx_stats_get_average(ptypes[i].stats(), &aver))) &&
            (estatsOK == (result2 = gmx_stats_get_sigma(ptypes[i].stats(), &sigma))))
        {
            gmx_poldata_set_ptype_polarizability(pd, ptypes[i].name().c_str(), aver, sigma);
            fprintf(fplog, "%-5s  %8.3f +/- %.3f\n", ptypes[i].name().c_str(), aver, sigma);
            int   nbins = 1+sqrt(nBootStrap);
            real *my_x, *my_y;
            gmx_stats_make_histogram(ptypes[i].stats(), 0, &nbins,
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
                    ptypes[i].name().c_str());
        }
    }
    fclose(xp);
    if (bZero)
    {
        const char *null = (const char *)"0";
        gmx_poldata_add_ptype(pd, null, NULL, null, a0, 0);
    }
    sfree(fpp);
    free_matrix(a);
    free_matrix(at);

    return ptypes.size();
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
    static char                     *zeropol     = (char *) NULL;
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
        { "-zeropol", FALSE, etSTR,  {&zeropol},
          "Polarizability for these poltypes (polarizability types) is fixed to zero. Note that all poltypes in gentop.dat that are non zero are fixed as well." },
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
        if (NULL != debug)
        {
            fprintf(debug, "%s", mpi->GetMolname().c_str());
            for(alexandria::MolecularCompositionIterator mci =
                    mpi->BeginMolecularComposition(); (mci < mpi->EndMolecularComposition()); ++mci)
            {
                fprintf(debug, "  %s", mci->GetCompName().c_str());
            }
            fprintf(debug, "\n");
        }
    }
    gms               = gmx_molselect_init(opt2fn("-sel", NFILE, fnm));
    FILE *fplog       = opt2FILE("-g", NFILE, fnm, "w");
    nalexandria_atypes = decompose_frag(fplog, opt2fn("-his", NFILE, fnm),
                                        pd, mp, bQM, lot, mindata,
                                        gms, bZero, bForceFit,
                                        nBootStrap, fractionBootStrap, seed,
                                        zeropol, oenv);
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
