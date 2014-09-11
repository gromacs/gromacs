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
#include <math.h>
#include "gromacs/math/utilities.h"
#include "gromacs/utility/futil.h"
#include "gromacs/utility/exceptions.h"
#include "gromacs/utility/smalloc.h"
#include "gromacs/utility/cstringutil.h"
#include "gromacs/math/vec.h"
#include "gromacs/fileio/xvgr.h"
#include "gromacs/legacyheaders/copyrite.h"
#include "gromacs/statistics/statistics.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/linearalgebra/matrix.h"
#include "gromacs/statistics/statistics.h"
#include "poldata.h"
#include "poldata_xml.h"
#include "molselect.h"
#include "molprop.h"
#include "molprop_util.h"
#include "molprop_tables.h"
#include "composition.h"
#include "stringutil.h"

typedef struct
{
    char *cat;
    bool  bPrint;
    int   count;
    real  rms;
} t_catli;

static int comp_catli(const void *a, const void *b)
{
    t_catli *ca = (t_catli *)a;
    t_catli *cb = (t_catli *)b;

    return strcmp(ca->cat, cb->cat);
}

typedef struct
{
    int      ncat;
    t_catli *catli;
} t_cats;

static void add_cat(t_cats *cats, const char *catname)
{
    int i;

    for (i = 0; (i < cats->ncat); i++)
    {
        if (strcasecmp(catname, cats->catli[i].cat) == 0)
        {
            cats->catli[i].count++;
            break;
        }
    }
    if (i == cats->ncat)
    {
        cats->ncat++;
        srenew(cats->catli, cats->ncat);
        cats->catli[cats->ncat-1].cat   = strdup(catname);
        cats->catli[cats->ncat-1].count = 1;
        cats->catli[cats->ncat-1].rms   = 0;
    }
}

namespace alexandria
{

class LongTable
{
    private:
        FILE                    *fp_;
        std::string              caption_;
        std::string              columns_;
        std::string              label_;
        std::vector<std::string> headLines_;
        bool                     bLandscape_;
    public:
        //! Constructor with a file pointer
        LongTable(FILE *fp, bool bLandscape);

        //! Constructor with a file name
        LongTable(const char *fn, bool bLandscape);

        //! Destructor
        ~LongTable() {};

        void setCaption(const char *caption) { caption_.assign(caption); }

        void setLabel(const char *label) { label_.assign(label); }

        //! Generate columns entry with first column left aligned and other center
        void setColumns(int nColumns);

        void setColumns(const char *columns) { columns_.assign(columns); }

        void addHeadLine(const char *headline) { headLines_.push_back(headline); }

        void printHeader();

        void printFooter();

        void printLine(const char * line);

        void printHLine();
};

LongTable::LongTable(FILE *fp, bool bLandscape)
{
    fp_         = fp;
    bLandscape_ = bLandscape;
    if (NULL == fp_)
    {
        GMX_THROW(gmx::FileIOError("File not open"));
    }
}

LongTable::LongTable(const char *fn, bool bLandscape)
{
    fp_         = fopen(fn, "w");
    bLandscape_ = bLandscape;
    if (NULL == fp_)
    {
        GMX_THROW(gmx::FileIOError("Could not open file"));
    }
}

void LongTable::setColumns(int nColumns)
{
    columns_.assign("l");
    for (int i = 1; (i < nColumns); i++)
    {
        columns_.append("c");
    }
}

void LongTable::printHeader()
{
    if (bLandscape_)
    {
        fprintf(fp_, "\\begin{landscape}\n");
    }
    fprintf(fp_, "\\begin{longtable}{%s}\n", columns_.c_str());
    printHLine();
    fprintf(fp_, "\\caption{%s}\n", caption_.c_str());
    fprintf(fp_, "\\label{%s}\\\\\n", label_.c_str());
    printHLine();
    for (unsigned int i = 0; (i < headLines_.size()); i++)
    {
        fprintf(fp_, "%s\\\\\n", headLines_[i].c_str());
    }
    printHLine();
    fprintf(fp_, "\\endfirsthead\n");
    printHLine();
    for (unsigned int i = 0; (i < headLines_.size()); i++)
    {
        fprintf(fp_, "%s\\\\\n", headLines_[i].c_str());
    }
    printHLine();
    fprintf(fp_, "\\endhead\n");
    printHLine();
    fprintf(fp_, "\\endfoot\n");
}

void LongTable::printFooter()
{
    fprintf(fp_, "\\end{longtable}\n");
    if (bLandscape_)
    {
        fprintf(fp_, "\\end{landscape}\n");
    }
    fflush(fp_);
}

void LongTable::printLine(const char *line)
{
    fprintf(fp_, "%s\\\\\n", line);
}

void LongTable::printHLine()
{
    fprintf(fp_, "\\hline\n");
}

}

static void stats_header(alexandria::LongTable &lt,
                         MolPropObservable      mpo,
                         t_qmcount             *qmc,
                         iMolSelect             ims)
{
    char caption[STRLEN];
    char label[STRLEN];
    char unit[32];

    lt.setColumns(1+qmc->n);

    for (int k = 0; (k < 2); k++)
    {
        if (0 == k)
        {
            switch (mpo)
            {
                case MPO_POLARIZABILITY:
                    sprintf(unit, "\\AA$^3$");
                    break;
                case MPO_DIPOLE:
                    sprintf(unit, "D");
                    break;
                case MPO_QUADRUPOLE:
                    sprintf(unit, "B");
                    break;
                case MPO_ENERGY:
                    sprintf(unit, "kJ/mol");
                    break;
                default:
                    gmx_fatal(FARGS, "Unknown property %s", mpo_name[mpo]);
            }
            snprintf(caption, STRLEN, "Performance of the different methods for predicting the molecular %s for molecules containing different chemical groups, given as the RMSD from experimental values (%s), and in brackets the number of molecules in this particular subset. {\\bf Data set: %s.} At the bottom the correlation coefficient R, the regression coefficient a and the intercept b are given as well as the normalized quality of the fit $\\chi^2$, the mean signed error (MSE) and the mean absolute error (MSA).",
                     mpo_name[mpo], unit, ims_names[ims]);
            snprintf(label, STRLEN, "%s_rmsd", mpo_name[mpo]);
            lt.setCaption(caption);
            lt.setLabel(label);
        }
        else
        {
            char hline[STRLEN];
            char koko[STRLEN];
            snprintf(hline, STRLEN, "Method ");
            for (int i = 0; (i < qmc->n); i++)
            {
                snprintf(koko, STRLEN, " & %s", qmc->method[i]);
                strncat(hline, koko, STRLEN-strlen(hline)-1);
            }
            lt.addHeadLine(hline);

            snprintf(hline, STRLEN, " ");
            for (int i = 0; (i < qmc->n); i++)
            {
                snprintf(koko, STRLEN, "& %s ", qmc->basis[i]);
                strncat(hline, koko, STRLEN-strlen(hline)-1);
            }
            lt.addHeadLine(hline);

            snprintf(hline, STRLEN, " ");
            for (int i = 0; (i < qmc->n); i++)
            {
                snprintf(koko, STRLEN, "& %s ", qmc->type[i]);
                strncat(hline, koko, STRLEN-strlen(hline)-1);
            }
            lt.addHeadLine(hline);
        }
    }
    lt.printHeader();
}

void gmx_molprop_stats_table(FILE                            *fp,
                             MolPropObservable                mpo,
                             std::vector<alexandria::MolProp> mp,
                             t_qmcount                       *qmc,
                             char                            *exp_type,
                             double                           outlier,
                             gmx_molselect_t                  gms,
                             iMolSelect                       ims)
{
    std::vector<alexandria::MolProp>::iterator mpi;
    std::vector<std::string>::iterator         si;
    int                                        i, line, k, m, N, nprint;
    double                                     exp_val, qm_val;
    real                                       rms, R, a, da, b, db, chi2;
    char                                       lbuf[256], tmp[32], catbuf[STRLEN];
    char                                       buf[32];
    t_cats                                    *cats;
    gmx_stats_t                                lsq, *lsqtot;
    alexandria::CompositionSpecs               cs;
    alexandria::LongTable                      lt(fp, true);
    const char                                *alex = cs.searchCS(alexandria::iCalexandria)->name();

    nprint = 0;
    snew(cats, 1);
    for (mpi = mp.begin(); (mpi < mp.end()); mpi++)
    {
        if ((ims == gmx_molselect_status(gms, mpi->GetIupac().c_str())) &&
            mpi->HasComposition(alex))
        {
            nprint++;
            for (si = mpi->BeginCategory(); (si < mpi->EndCategory()); si++)
            {
                add_cat(cats, si->c_str());
            }
        }
    }
    if (nprint <= 0)
    {
        sfree(cats);
        return;
    }

    stats_header(lt, mpo, qmc, ims);

    qsort(cats->catli, cats->ncat, sizeof(cats->catli[0]), comp_catli);
    line = 0;
    for (i = 0; (i < cats->ncat); i++)
    {
        snprintf(catbuf, STRLEN,
                 " %s (%d) ", cats->catli[i].cat, cats->catli[i].count);
        int nqmres  = 0;
        int nexpres = 0;
        for (k = 0; (k < qmc->n); k++)
        {
            snprintf(lbuf, 256, "%s/%s", qmc->method[k], qmc->basis[k]);
            lsq = gmx_stats_init();
            for (mpi = mp.begin(); (mpi < mp.end()); mpi++)
            {
                if ((gmx_molselect_status(gms, mpi->GetIupac().c_str()) == ims) &&
                    (mpi->HasComposition(alex)))
                {
                    if (mpi->SearchCategory(cats->catli[i].cat) == 1)
                    {
                        if (mpi->GetProp(mpo, iqmExp, NULL,
                                         NULL /*qmc->conf[m]*/, exp_type, &exp_val))
                        {
                            nexpres = 1;
                            if (mpi->GetProp(mpo, iqmQM, lbuf, NULL /*qmc->conf[m]*/,
                                             qmc->type[k], &qm_val))
                            {
                                if (debug)
                                {
                                    fprintf(debug, "%s %s k = %d - TAB4\n",
                                            mpi->GetMolname().c_str(), cats->catli[i].cat, k);
                                }
                                gmx_stats_add_point(lsq, exp_val, qm_val, 0, 0);
                            }
                        }
                    }
                }
            }
            if (outlier > 0)
            {
                gmx_stats_remove_outliers(lsq, outlier);
            }
            if (gmx_stats_get_rmsd(lsq, &rms) == estatsOK)
            {
                sprintf(tmp, "& %8.2f", rms);
                strncat(catbuf, tmp, STRLEN-strlen(catbuf)-1);
                nqmres++;
            }
            else
            {
                strncat(catbuf, "& -", STRLEN-strlen(catbuf)-1);
            }

            if (gmx_stats_done(lsq) == estatsOK)
            {
                sfree(lsq);
            }
        }
        cats->catli[i].bPrint = ((nqmres > 0) && (nexpres > 0));
        if (cats->catli[i].bPrint)
        {
            lt.printLine(catbuf);
            line++;
        }
    }
    snprintf(catbuf, STRLEN, "All");
    snew(lsqtot, qmc->n);
    for (k = 0; (k < qmc->n); k++)
    {
        snprintf(lbuf, 256, "%s/%s", qmc->method[k], qmc->basis[k]);
        lsqtot[k] = gmx_stats_init();
        for (mpi = mp.begin(); (mpi < mp.end()); mpi++)
        {
            if ((gmx_molselect_status(gms, mpi->GetIupac().c_str()) == ims) &&
                (mpi->HasComposition(alex)))
            {
                for (m = 0; (m < qmc->nconf); m++)
                {
                    if ((mpi->GetProp(mpo, iqmExp, NULL, qmc->conf[m], exp_type, &exp_val)) &&
                        (mpi->GetProp(mpo, iqmQM, lbuf, qmc->conf[m], qmc->type[k], &qm_val)))
                    {
                        gmx_stats_add_point(lsqtot[k], exp_val, qm_val, 0, 0);
                    }
                }
            }
        }
        if (gmx_stats_get_rmsd(lsqtot[k], &rms) == estatsOK)
        {
            snprintf(buf, 32,  "& %8.2f", rms);
        }
        else
        {
            snprintf(buf, 32, "& -");
        }
        strncat(catbuf, buf, STRLEN-strlen(catbuf)-1);
    }
    lt.printLine(catbuf);
    lt.printHLine();

    snprintf(catbuf, STRLEN, "N");
    for (k = 0; (k < qmc->n); k++)
    {
        /* gmx_stats_remove_outliers(lsqtot[k],3); */
        gmx_stats_get_npoints(lsqtot[k], &N);
        snprintf(buf, 32, "& %d", N);
        strncat(catbuf, buf, STRLEN-strlen(catbuf)-1);
    }
    lt.printLine(catbuf);

    snprintf(catbuf, STRLEN, "a");
    for (k = 0; (k < qmc->n); k++)
    {
        if (gmx_stats_get_ab(lsqtot[k], elsqWEIGHT_NONE, &a, &b, &da, &db, &chi2, &R) ==
            estatsOK)
        {
            snprintf(buf, 32, "& %8.2f(%4.2f)", a, da);
        }
        else
        {
            snprintf(buf, 32, "& -");
        }
        strncat(catbuf, buf, STRLEN-strlen(catbuf)-1);
    }
    lt.printLine(catbuf);

    snprintf(catbuf, STRLEN, "b");
    for (k = 0; (k < qmc->n); k++)
    {
        if (gmx_stats_get_ab(lsqtot[k], elsqWEIGHT_NONE, &a, &b, &da, &db, &chi2, &R) ==
            estatsOK)
        {
            snprintf(buf, 32, "& %8.2f(%4.2f)", b, db);
        }
        else
        {
            snprintf(buf, 32, "& -");
        }
        strncat(catbuf, buf, STRLEN-strlen(catbuf)-1);
    }
    lt.printLine(catbuf);

    snprintf(catbuf, STRLEN, "R (\\%%)");
    for (k = 0; (k < qmc->n); k++)
    {
        if (gmx_stats_get_corr_coeff(lsqtot[k], &R) == estatsOK)
        {
            snprintf(buf, 32, "& %8.2f", 100*R);
        }
        else
        {
            snprintf(buf, 32, "& -");
        }
        strncat(catbuf, buf, STRLEN-strlen(catbuf)-1);
    }
    lt.printLine(catbuf);

    snprintf(catbuf, STRLEN, "$\\chi^2$");
    for (k = 0; (k < qmc->n); k++)
    {
        if (gmx_stats_get_ab(lsqtot[k], elsqWEIGHT_NONE, &a, &b, &da, &db, &chi2, &R) ==
            estatsOK)
        {
            snprintf(buf, 32, "& %8.2f", chi2);
        }
        else
        {
            snprintf(buf, 32, "& -");
        }
        strncat(catbuf, buf, STRLEN-strlen(catbuf)-1);
    }
    lt.printLine(catbuf);

    snprintf(catbuf, STRLEN, "MSE");
    for (k = 0; (k < qmc->n); k++)
    {
        real mse;
        if (gmx_stats_get_mse_mae(lsqtot[k], &mse, NULL) ==
            estatsOK)
        {
            snprintf(buf, 32, "& %8.2f", mse);
        }
        else
        {
            snprintf(buf, 32, "& -");
        }
        strncat(catbuf, buf, STRLEN-strlen(catbuf)-1);
    }
    lt.printLine(catbuf);

    snprintf(catbuf, STRLEN, "MAE");
    for (k = 0; (k < qmc->n); k++)
    {
        real mae;
        if (gmx_stats_get_mse_mae(lsqtot[k], NULL, &mae) ==
            estatsOK)
        {
            snprintf(buf, 32, "& %8.2f", mae);
        }
        else
        {
            snprintf(buf, 32, "& -");
        }
        strncat(catbuf, buf, STRLEN-strlen(catbuf)-1);
    }
    lt.printLine(catbuf);
    lt.printFooter();

    for (k = 0; (k < qmc->n); k++)
    {
        if (gmx_stats_done(lsqtot[k]) == estatsOK)
        {
            sfree(lsqtot[k]);
        }
    }
    sfree(lsqtot);
}

static void composition_header(alexandria::LongTable &lt,
                               iMolSelect             ims)
{
    char caption[STRLEN];

    snprintf(caption, STRLEN, "Decomposition of molecules into Alexandria atom types. {\\bf Data set: %s.} Charge is given when not zero, multiplicity is given when not 1.",
             ims_names[ims]);
    lt.setCaption(caption);
    lt.setLabel("frag_defs");
    lt.setColumns("p{75mm}ll");
    lt.addHeadLine("Molecule & Formula  & Types");
    lt.printHeader();
}

void gmx_molprop_composition_table(FILE *fp, std::vector<alexandria::MolProp> mp,
                                   gmx_molselect_t gms, iMolSelect ims)
{
    std::vector<alexandria::MolProp>::iterator mpi;
    alexandria::MolecularCompositionIterator   mci;
    int                                        q, m, iline, nprint;
    char                                       qbuf[32], longbuf[STRLEN], buf[256];
    alexandria::LongTable                      lt(fp, true);
    alexandria::CompositionSpecs               cs;
    const char                                *alex = cs.searchCS(alexandria::iCalexandria)->name();

    nprint = 0;
    for (mpi = mp.begin(); (mpi < mp.end()); mpi++)
    {
        if ((ims == gmx_molselect_status(gms, mpi->GetIupac().c_str())) &&
            (mpi->HasComposition(alex)))
        {
            nprint++;
        }
    }
    if (nprint <= 0)
    {
        return;
    }

    composition_header(lt, ims);
    iline = 0;
    for (mpi = mp.begin(); (mpi < mp.end()); mpi++)
    {
        if ((ims == gmx_molselect_status(gms, mpi->GetIupac().c_str())) &&
            (mpi->HasComposition(alex)))
        {
            q = mpi->GetCharge();
            m = mpi->GetMultiplicity();
            if ((q != 0) || (m != 1))
            {
                if ((q != 0) && (m != 1))
                {
                    sprintf(qbuf, " (q=%c%d, mult=%d)", (q < 0) ? '-' : '+', abs(q), m);
                }
                else if (q != 0)
                {
                    sprintf(qbuf, " (q=%c%d)", (q < 0) ? '-' : '+', abs(q));
                }
                else
                {
                    sprintf(qbuf, " (mult=%d)", m);
                }
            }
            else
            {
                qbuf[0] = '\0';
            }
            snprintf(longbuf, STRLEN, "%3d. %s%s & %s & ",
                     ++iline,
                     mpi->GetIupac().c_str(),
                     qbuf,
                     mpi->GetTexFormula().c_str());
            mci = mpi->SearchMolecularComposition(cs.searchCS(alexandria::iCalexandria)->name());
            if (mci != mpi->EndMolecularComposition())
            {
                alexandria::AtomNumIterator ani;

                for (ani = mci->BeginAtomNum(); (ani < mci->EndAtomNum()); ani++)
                {
                    snprintf(buf, 256, " %d %s\t",
                             ani->GetNumber(), ani->GetAtom().c_str());
                    strncat(longbuf, buf, STRLEN-sizeof(longbuf)-1);
                }
            }
            lt.printLine(longbuf);
        }
    }
    lt.printFooter();
}

typedef struct {
    char  *cat;
    int    nmol;
    char **molec;
} t_cat_stat;

static int comp_cats(const void *a, const void *b)
{
    t_cat_stat *ca = (t_cat_stat *)a;
    t_cat_stat *cb = (t_cat_stat *)b;

    return strcmp(ca->cat, cb->cat);
}

static void add_cats(int *ncs, t_cat_stat **cs, const char *iupac, const char *mycat)
{
    int j, nmol;

    for (j = 0; (j < *ncs); j++)
    {
        if (strcmp((*cs)[j].cat, mycat) == 0)
        {
            break;
        }
    }
    if (j == *ncs)
    {
        srenew(*cs, ++(*ncs));
        (*cs)[j].cat   = strdup(mycat);
        (*cs)[j].nmol  = 0;
        (*cs)[j].molec = NULL;
    }
    nmol = ++(*cs)[j].nmol;
    srenew((*cs)[j].molec, nmol);
    (*cs)[j].molec[nmol-1] = strdup(iupac);
}

static void category_header(alexandria::LongTable &lt)
{
    char longbuf[STRLEN];

    lt.setColumns("lcp{120mm}");
    for (int k = 0; (k < 2); k++)
    {
        if (0 == k)
        {
            snprintf(longbuf, STRLEN, "Molecules that are part of each category used for statistics.");
            lt.setCaption(longbuf);
            lt.setLabel("stats");
        }
        lt.addHeadLine("Category & N & Molecule(s)");
    }
    lt.printHeader();
}

void gmx_molprop_category_table(FILE *fp,
                                std::vector<alexandria::MolProp> mp,
                                gmx_molselect_t gms, iMolSelect ims)
{
    std::vector<alexandria::MolProp>::iterator mpi;
    std::vector<std::string>::iterator         si;
    int                   i, j, ncs;
    t_cat_stat           *cs = NULL;
    const char           *iupac;
    alexandria::LongTable lt(fp, true);

    ncs = 0;
    for (mpi = mp.begin(); (mpi < mp.end()); mpi++)
    {
        iupac = mpi->GetIupac().c_str();
        if (ims == gmx_molselect_status(gms, iupac))
        {
            for (si = mpi->BeginCategory(); (si < mpi->EndCategory()); si++)
            {
                add_cats(&ncs, &cs, iupac, si->c_str());
            }
        }
    }

    if (ncs > 0)
    {
        char longbuf[STRLEN], buf[256];
        qsort(cs, ncs, sizeof(cs[0]), comp_cats);
        category_header(lt);
        for (i = 0; (i < ncs); i++)
        {
            snprintf(longbuf, STRLEN, "%s & %d &", cs[i].cat, cs[i].nmol);
            for (j = 0; (j < cs[i].nmol-1); j++)
            {
                snprintf(buf, STRLEN, "%s, ", cs[i].molec[j]);
                strncat(longbuf, buf, STRLEN-strlen(longbuf)-1);
            }
            snprintf(buf, STRLEN, "%s, ", cs[i].molec[j]);
            strncat(longbuf, buf, STRLEN-strlen(longbuf)-1);
            lt.printLine(longbuf);
        }
    }
    lt.printFooter();
}

static void atomtype_tab_header(alexandria::LongTable &lt)
{
    char longbuf[STRLEN];
    alexandria::CompositionSpecs cs;

    lt.setColumns("ccccccc");

    snprintf(longbuf, STRLEN, "Definition of atom types for polarization (AX column). The number of datapoints used for fitting the polarizability of atomtypes is given in N$_{Exp}$ and N$_{QM}$ for experimental data or quantum chemistry, respectively. The columns Ahc and Ahp contain group polarizabilites computed using Miller's equation~\\protect\\cite{Miller1979a} and parameters~\\protect\\cite{Miller1990a} and Kang and Jhon's method~\\cite{Kang1982a} with the parametrization of Miller~\\protect\\cite{Miller1990a} respectively. The column BS contains the equivalent using the polarizabilities of Bosque and Sales~\\protect\\cite{Bosque2002a}.");
    lt.setCaption(longbuf);
    lt.setLabel("fragments");
    snprintf(longbuf, STRLEN, "Name  & N$_{Exp}$ & N$_{QM}$ & \\multicolumn{4}{c}{Polarizability}");
    lt.addHeadLine(longbuf);
    snprintf(longbuf, STRLEN, "& & & %s ($\\sigma_{AX}$) & Ahc & Ahp & %s ",
             cs.searchCS(alexandria::iCalexandria)->abbreviation(),
             cs.searchCS(alexandria::iCbosque)->abbreviation());
    lt.addHeadLine(longbuf);
    lt.printHeader();
}

typedef struct {
    char       *ptype, *miller, *bosque;
    gmx_stats_t lsq;
    int         nexp, nqm;
} t_sm_lsq;

static void gmx_molprop_atomtype_polar_table(FILE                            *fp,
                                             gmx_poldata_t                    pd,
                                             std::vector<alexandria::MolProp> mp,
                                             char                            *lot,
                                             char                            *exp_type)
{
    std::vector<alexandria::MolProp>::iterator mpi;
    double                                     ahc, ahp, bos_pol, alexandria_pol, sig_pol;
    char                                      *ptype;
    char                                      *miller, *bosque;
    char                                       longbuf[STRLEN];
    MolPropObservable                          mpo = MPO_POLARIZABILITY;
    alexandria::LongTable                      lt(fp, false);
    alexandria::CompositionSpecs               cs;
    const char                                *alex = cs.searchCS(alexandria::iCalexandria)->name();

    /* Prepare printing it! */
    atomtype_tab_header(lt);

    /* First gather statistics from different input files.
     * Note that the input files
     * do not need to have the same sets of types,
     * as we check for the type name.
     */
    while (1 == gmx_poldata_get_ptype(pd,
                                      &ptype,
                                      &miller,
                                      &bosque,
                                      &alexandria_pol,
                                      &sig_pol))
    {
        if (alexandria_pol > 0)
        {
            int atomnumber;
            int nexp   = 0;
            int nqm    = 0;

            /* Now loop over the poltypes and count number of occurrences
             * in molecules
             */
            for (mpi = mp.begin(); (mpi < mp.end()); mpi++)
            {
                alexandria::MolecularCompositionIterator mci = mpi->SearchMolecularComposition(alex);

                if (mci != mpi->EndMolecularComposition())
                {
                    bool bFound = false;

                    for (alexandria::AtomNumIterator ani = mci->BeginAtomNum(); !bFound && (ani < mci->EndAtomNum()); ++ani)
                    {
                        const char *pt =
                            gmx_poldata_atype_to_ptype(pd,
                                                       ani->GetAtom().c_str());
                        if ((NULL != pt) && (strcasecmp(pt, ptype) == 0))
                        {
                            bFound = true;
                        }
                    }
                    if (bFound)
                    {
                        double val;
                        if (mpi->GetProp(mpo, iqmExp, lot, NULL, exp_type, &val))
                        {
                            nexp++;
                        }
                        else if (mpi->GetProp(mpo, iqmQM, lot, NULL, NULL, &val))
                        {
                            nqm++;
                        }
                    }
                }
            }

            /* Determine Miller and Bosque polarizabilities for this Alexandria element */
            ahc = ahp = bos_pol = 0;
            if (1 == gmx_poldata_get_miller_pol(pd, miller,
                                                &atomnumber, &ahc, &ahp))
            {
                ahc = (4.0/atomnumber)*sqr(ahc);
            }
            if (0 == gmx_poldata_get_bosque_pol(pd, bosque, &bos_pol))
            {
                bos_pol = 0;
            }
            /* Compute how many molecules contributed to the optimization
             * of this polarizability.
             */
            /* Construct group name from element composition */
            /* strncpy(group,smlsq[j].bosque,sizeof(group));*/

            snprintf(longbuf, STRLEN, "%s & %s & %s & %s (%s) & %s & %s & %s",
                     ptype,
                     (nexp > 0)     ? gmx_itoa(nexp).c_str()     : "",
                     (nqm > 0)      ? gmx_itoa(nqm).c_str()      : "",
                     (alexandria_pol > 0)  ? gmx_ftoa(alexandria_pol).c_str()  : "",
                     (sig_pol > 0)  ? gmx_ftoa(sig_pol).c_str() : "-",
                     (ahc > 0)         ? gmx_ftoa(ahc).c_str()         : "",
                     (ahp > 0)         ? gmx_ftoa(ahp).c_str()         : "",
                     (bos_pol > 0)     ? gmx_ftoa(bos_pol).c_str()     : "");
            lt.printLine(longbuf);
        }
    }
    lt.printFooter();
    fflush(fp);
}

static void gmx_molprop_atomtype_dip_table(FILE *fp, gmx_poldata_t pd)
{
    int     i, k, m, cur = 0;
    char   *elem, *gt_type[2] = { NULL, NULL };
    char   *spref, *desc, *ptype, *btype;
#define prev (1-cur)
#define NEQG 5
    ChargeGenerationModel eqgcol[NEQG] = { eqgAXp, eqgAXs, eqgAXg };
    char                  longbuf[STRLEN], buf[256];
    int                   npcol[NEQG]  = { 2, 3, 3 };
    const char           *clab[3]      = { "$J_0$", "$\\chi_0$", "$\\zeta$" };
    int                   ncol;
    alexandria::LongTable lt(fp, true);

    ncol = 1;
    for (i = 0; (i < NEQG); i++)
    {
        ncol += npcol[i];
    }
    lt.setCaption("Electronegativity equalization parameters for Alexandria models. $J_0$ and $\\chi_0$ in eV, $\\zeta$ in 1/nm.");
    lt.setLabel("eemparams");
    lt.setColumns(ncol);

    snprintf(longbuf, STRLEN, " ");
    for (i = 0; (i < NEQG); i++)
    {
        snprintf(buf, 256, " & \\multicolumn{%d}{c}{%s}", npcol[i],
                 get_eemtype_name(eqgcol[i]));
        strncat(longbuf, buf, STRLEN-strlen(longbuf)-1);
    }
    lt.addHeadLine(longbuf);

    snprintf(longbuf, STRLEN, " ");
    for (i = 0; (i < NEQG); i++)
    {
        for (m = 0; (m < npcol[i]); m++)
        {
            snprintf(buf, 256, " & %s", clab[m]);
            strncat(longbuf, buf, STRLEN-strlen(longbuf)-1);
        }
    }
    lt.addHeadLine(longbuf);
    lt.printHeader();

    while (1 == gmx_poldata_get_atype(pd,
                                      &elem,
                                      &desc,
                                      &(gt_type[cur]),
                                      &ptype,
                                      &btype,
                                      &spref))
    {
        if (((NULL == gt_type[prev]) || (strcmp(gt_type[cur], gt_type[prev]) != 0)))
        {
            snprintf(longbuf, STRLEN, "%s\n", gt_type[cur]);
            for (k = 0; (k < NEQG); k++)
            {
                if (gmx_poldata_have_eem_support(pd, eqgcol[k], gt_type[cur], false))
                {
                    snprintf(buf, 256, " & %.3f", gmx_poldata_get_j00(pd, eqgcol[k], gt_type[cur]));
                    strncat(longbuf, buf, STRLEN-strlen(longbuf)-1);
                    snprintf(buf, 256, " & %.3f", gmx_poldata_get_chi0(pd, eqgcol[k], gt_type[cur]));
                    strncat(longbuf, buf, STRLEN-strlen(longbuf)-1);
                    if (npcol[k] == 3)
                    {
                        snprintf(buf, 256, " & %.3f", gmx_poldata_get_zeta(pd, eqgcol[k], gt_type[cur], 1));
                        strncat(longbuf, buf, STRLEN-strlen(longbuf)-1);
                    }
                    if (npcol[k] == 4)
                    {
                        snprintf(buf, 256, " & %.3f", gmx_poldata_get_zeta(pd, eqgcol[k], gt_type[cur], 2));
                        strncat(longbuf, buf, STRLEN-strlen(longbuf)-1);
                    }
                }
                else
                {
                    for (m = 0; (m < npcol[k]); m++)
                    {
                        snprintf(buf, 256, " & ");
                        strncat(longbuf, buf, STRLEN-strlen(longbuf)-1);
                    }
                }
            }
            lt.printLine(longbuf);
        }
        cur = prev;
    }
    lt.printFooter();
}

void gmx_molprop_atomtype_table(FILE *fp, bool bPolar,
                                gmx_poldata_t pd,
                                std::vector<alexandria::MolProp> mp,
                                char *lot, char *exp_type)
{
    if (bPolar)
    {
        gmx_molprop_atomtype_polar_table(fp, pd, mp, lot, exp_type);
    }
    else
    {
        gmx_molprop_atomtype_dip_table(fp, pd);
    }
}

static void prop_header(alexandria::LongTable &lt,
                        const char            *property,
                        const char            *unit,
                        real                   rel_toler,
                        real                   abs_toler,
                        t_qmcount             *qmc,
                        iMolSelect             ims,
                        bool                   bPrintConf,
                        bool                   bPrintBasis,
                        bool                   bPrintMultQ)
{
    int  i, k, nc;
    char columns[256];
    char longbuf[STRLEN];
    char buf[256];

    snprintf(columns, 256, "p{75mm}");
    nc = 2 + qmc->n;
    if (bPrintMultQ)
    {
        nc += 2;
    }
    if (bPrintConf)
    {
        nc++;
    }
    for (i = 0; (i < nc); i++)
    {
        strncat(columns, "c", 256-strlen(columns)-1);
    }
    lt.setColumns(columns);

    for (k = 0; (k < 2); k++)
    {
        if (0 == k)
        {
            snprintf(longbuf, STRLEN, "Comparison of experimental %s to calculated values. {\\bf Data set: %s}. Calculated numbers that are more than %.0f%s off the experimental values are printed in bold, more than %.0f%s off in bold red.",
                     property, ims_names[ims],
                     (abs_toler > 0) ? abs_toler   : 100*rel_toler,
                     (abs_toler > 0) ? unit : "\\%",
                     (abs_toler > 0) ? 2*abs_toler : 200*rel_toler,
                     (abs_toler > 0) ? unit : "\\%");
            lt.setCaption(longbuf);
            snprintf(longbuf, STRLEN, "%s", ims_names[ims]);
            lt.setLabel(longbuf);
        }
        else
        {
            snprintf(longbuf, STRLEN, "Molecule & Form. %s %s & Exper. ",
                     bPrintMultQ ? "& q & mult" : "",
                     bPrintConf  ? "& Conf." : "");
            for (i = 0; (i < qmc->n); i++)
            {
                snprintf(buf, 256, "& %s", qmc->method[i]);
                strncat(longbuf, buf, STRLEN-strlen(longbuf)-1);
            }
            lt.addHeadLine(longbuf);

            if (bPrintBasis)
            {
                snprintf(longbuf, STRLEN, " & & %s %s",
                         bPrintMultQ ? "& &" : "",
                         bPrintConf ? "&" : "");
                for (i = 0; (i < qmc->n); i++)
                {
                    snprintf(buf, 256, "& %s", qmc->basis[i]);
                    strncat(longbuf, buf, STRLEN-strlen(longbuf)-1);
                }
                lt.addHeadLine(longbuf);
            }
            snprintf(longbuf, STRLEN, "Type & &%s %s",
                     bPrintMultQ ? "& &" : "",
                     bPrintConf ? "&" : "");
            for (i = 0; (i < qmc->n); i++)
            {
                snprintf(buf, 256, "& %s", qmc->type[i]);
                strncat(longbuf, buf, STRLEN-strlen(longbuf)-1);
            }
            lt.addHeadLine(longbuf);
        }
    }
    lt.printHeader();
}

static int outside(real vexp, real vcalc, real rel_toler, real abs_toler)
{
    real rdv, adv = fabs(vexp-vcalc);
    if (abs_toler > 0)
    {
        if (adv > 2*abs_toler)
        {
            return 2;
        }
        else if (adv > abs_toler)
        {
            return 1;
        }
        return 0;
    }
    else
    {
        if (vexp == 0)
        {
            return 0;
        }
        rdv = adv/vexp;

        if (rdv > 2*rel_toler)
        {
            return 2;
        }
        else if (rdv > rel_toler)
        {
            return 1;
        }
        return 0;
    }
}

class ExpData
{
    public:
        double      val_, err_;
        std::string ref_, conf_, type_;
        ExpData(double val, double err, std::string ref, std::string conf, std::string type)
        { val_ = val; err_ = err; ref_ = ref; conf_ = conf; type_ = type; }
};

class CalcData
{
    public:
        double val_, err_;
        int    found_;
        CalcData(double val, double err, int found)
        { val_ = val; err_ = err; found_ = found; }
};

void gmx_molprop_prop_table(FILE *fp, MolPropObservable mpo,
                            real rel_toler, real abs_toler,
                            std::vector<alexandria::MolProp> mp,
                            t_qmcount *qmc, bool bPrintAll,
                            bool bPrintBasis, bool bPrintMultQ,
                            gmx_molselect_t gms, iMolSelect ims)
{
    alexandria::MolecularQuadrupoleIterator qi;
    alexandria::MolecularEnergyIterator     mei;

    int                                     j, iprint = 0, nexp, max_nexp;
#define BLEN 1024
    char                                    lbuf[BLEN], myline[BLEN], mylbuf[BLEN], vbuf[BLEN];
    double                                  calc_val, calc_err, vc;
    int                                     nprint;
    double                                  dvec[DIM];
    tensor                                  quadrupole;
    bool                                    bPrintConf;
    alexandria::CompositionSpecs            cs;
    const char                             *alex = cs.searchCS(alexandria::iCalexandria)->name();

    alexandria::LongTable                   lt(fp, true);

    nprint = 0;
    for (alexandria::MolPropIterator mpi = mp.begin(); (mpi < mp.end()); mpi++)
    {
        if ((ims == gmx_molselect_status(gms, mpi->GetIupac().c_str())) &&
            (mpi->HasComposition(alex)))
        {
            nprint++;
        }
    }
    if (nprint <= 0)
    {
        return;
    }
    bPrintConf = false; //(mpo == MPO_DIPOLE);
    prop_header(lt, mpo_name[mpo], mpo_unit[mpo],
                rel_toler, abs_toler, qmc,
                ims, bPrintConf, bPrintBasis, bPrintMultQ);
    for (alexandria::MolPropIterator mpi = mp.begin(); (mpi < mp.end()); mpi++)
    {
        if ((ims == gmx_molselect_status(gms, mpi->GetIupac().c_str())) &&
            (mpi->HasComposition(alex)))
        {
            std::vector<ExpData>  ed;
            std::vector<CalcData> cd;
            for (alexandria::ExperimentIterator ei = mpi->BeginExperiment(); (ei < mpi->EndExperiment()); ei++)
            {
                switch (mpo)
                {
                    case MPO_DIPOLE:
                        for (alexandria::MolecularDipPolarIterator mdi = ei->BeginDipole(); (mdi < ei->EndDipole()); mdi++)
                        {
                            ed.push_back(ExpData(mdi->GetAver(),
                                                 mdi->GetError(),
                                                 ei->GetReference(),
                                                 ei->GetConformation(),
                                                 mdi->GetType()));
                        }
                        break;
                    case MPO_POLARIZABILITY:
                        for (alexandria::MolecularDipPolarIterator mdi = ei->BeginPolar(); (mdi < ei->EndPolar()); mdi++)
                        {
                            ed.push_back(ExpData(mdi->GetAver(),
                                                 mdi->GetError(),
                                                 ei->GetReference(),
                                                 ei->GetConformation(),
                                                 mdi->GetType()));
                        }
                        break;
                    case MPO_ENERGY:
                        for (mei = ei->BeginEnergy(); (mei < ei->EndEnergy()); mei++)
                        {
                            ed.push_back(ExpData(mei->GetValue(),
                                                 mei->GetError(),
                                                 ei->GetReference(),
                                                 ei->GetConformation(),
                                                 mei->GetType()));
                        }
                        break;
                    default:
                        gmx_fatal(FARGS, "No support for for mpo %d", mpo);
                        break;
                }
            }
            max_nexp = ed.size();
            if (max_nexp <= 0)
            {
                max_nexp = 1;
            }
            for (nexp = 0; (nexp < max_nexp); nexp++)
            {
                for (j = 0; (j < qmc->n); j++)
                {
                    sprintf(lbuf, "%s/%s", qmc->method[j], qmc->basis[j]);
                    if (mpi->GetPropRef(mpo, iqmQM, lbuf, NULL, //conf_exp[nexp].c_str(),
                                        qmc->type[j], &calc_val, &calc_err,
                                        NULL, NULL, dvec, quadrupole))
                    {
                        cd.push_back(CalcData(calc_val, calc_err, 1));
                    }
                    else
                    {
                        cd.push_back(CalcData(0, 0, 0));
                    }
                }
                if ((bPrintAll || (ed.size() > 0))  && (cd.size() > 0))
                {
                    if (0 == nexp)
                    {
                        iprint++;
                    }
                    myline[0] = '\0';
                    if (nexp == 0)
                    {
                        if (bPrintMultQ)
                        {
                            snprintf(myline, BLEN, "%d. %-15s & %s & %d & %d",
                                     iprint, mpi->GetIupac().c_str(),
                                     mpi->GetTexFormula().c_str(),
                                     mpi->GetCharge(),
                                     mpi->GetMultiplicity());
                        }
                        else
                        {
                            snprintf(myline, BLEN, "%d. %-15s & %s",
                                     iprint, mpi->GetIupac().c_str(),
                                     mpi->GetTexFormula().c_str());
                        }
                    }
                    else
                    {
                        snprintf(myline, BLEN, " & ");
                    }
                    if (bPrintConf)
                    {
                        snprintf(mylbuf, BLEN, "      & %s ", ((ed[nexp].conf_.size() > 0) ?
                                                               ed[nexp].conf_.c_str() : "-"));
                        strncat(myline, mylbuf, BLEN-strlen(myline)-1);
                    }
                    if (ed.size() > 0)
                    {
                        sprintf(mylbuf, "& %8.3f", ed[nexp].val_);
                        strncat(myline, mylbuf, BLEN-strlen(myline)-1);
                        if (ed[nexp].err_ > 0)
                        {
                            sprintf(mylbuf, "(%.3f)", ed[nexp].err_);
                            strncat(myline, mylbuf, BLEN-strlen(myline)-1);
                        }
                        if (strcmp(ed[nexp].ref_.c_str(), "Maaren2014a") == 0)
                        {
                            sprintf(mylbuf, " (*)");
                        }
                        else
                        {
                            sprintf(mylbuf, "~\\cite{%s} ", ed[nexp].ref_.c_str());
                        }
                        strncat(myline, mylbuf, BLEN-strlen(myline)-1);
                    }
                    else
                    {
                        sprintf(mylbuf, "& - ");
                        strncat(myline, mylbuf, BLEN-strlen(myline)-1);
                    }
                    for (j = 0; (j < qmc->n); j++)
                    {
                        if (cd[j].found_ > 0)
                        {
                            vc = cd[j].val_;
                            if (cd[j].err_ > 0)
                            {
                                sprintf(vbuf, "%8.2f(%.2f)", vc, cd[j].err_);
                            }
                            else
                            {
                                sprintf(vbuf, "%8.2f", vc);
                            }
                            if (ed.size() > 0)
                            {
                                int oo = outside(ed[nexp].val_, vc, rel_toler, abs_toler);
                                switch (oo)
                                {
                                    case 2:
                                        sprintf(mylbuf, "& \\textcolor{Red}{\\bf %s} ", vbuf);
                                        break;
                                    case 1:
                                        sprintf(mylbuf, "& {\\bf %s} ", vbuf);
                                        break;
                                    default:
                                        sprintf(mylbuf, "& %s ", vbuf);
                                }
                            }
                            else
                            {
                                sprintf(mylbuf, "& %s ", vbuf);
                            }
                            strncat(myline, mylbuf, BLEN-strlen(myline)-1);
                        }
                        else
                        {
                            sprintf(mylbuf, "& ");
                            strncat(myline, mylbuf, BLEN-strlen(myline));
                        }
                    }
                    lt.printLine(myline);
                }
            }
        }
    }
    lt.printFooter();
}
