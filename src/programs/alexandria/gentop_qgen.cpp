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
#include <ctype.h>
#include "gromacs/legacyheaders/macros.h"
#include "gromacs/legacyheaders/copyrite.h"
#include "gromacs/bonded/bonded.h"
#include "gromacs/utility/cstringutil.h"
#include "gromacs/utility/smalloc.h"
#include "gromacs/fileio/strdb.h"
#include "gromacs/fileio/confio.h"
#include "gromacs/math/units.h"
#include "gromacs/random/random.h"
#include "gromacs/legacyheaders/txtdump.h"
#include "gromacs/legacyheaders/readinp.h"
#include "gromacs/legacyheaders/names.h"
#include "gromacs/math/vec.h"
#include "gromacs/topology/atomprop.h"
#include "gromacs/gmxpreprocess/grompp.h"
#include "gromacs/linearalgebra/matrix.h"
#include "coulombintegrals/coulombintegrals.h"
#include "molprop.h"
#include "poldata.h"
#include "gentop_nm2type.h"
#include "gentop_qgen.h"
#include "gmx_resp.h"

typedef struct gentop_qgen
{
    gmx_bool              bWarned;
    ChargeGenerationModel iModel;
    int                   natom, eQGEN;
    real                  qtotal, chieq, hfac, epsr;
    /* For each atom i there is an elem, atomnr, chi0, rhs, j00 and x */
    char                **elem;
    int                  *atomnr;
    real                 *chi0, *rhs, *j00;
    rvec                 *x;
    /* Jab is a matrix over atom pairs */
    real                **Jab;
    /* For each atom i there are nZeta[i] row, q and zeta entries */
    int                  *nZeta;
    int                 **row;
    gmx_bool              bAllocSave;
    real                **q, **zeta, **qsave, **zetasave;
} gentop_qgen;

static void gentop_qgen_save_params(gentop_qgen_t qgen, gmx_resp_t gr)
{
    int i, j;

    if (!qgen->bAllocSave)
    {
        snew(qgen->qsave, qgen->natom);
        snew(qgen->zetasave, qgen->natom);
    }
    for (i = 0; (i < qgen->natom); i++)
    {
        if (!qgen->bAllocSave)
        {
            snew(qgen->qsave[i], qgen->nZeta[i]);
            snew(qgen->zetasave[i], qgen->nZeta[i]);
        }
        for (j = 0; (j < qgen->nZeta[i]); j++)
        {
            if (NULL != gr)
            {
                qgen->q[i][j]    = gmx_resp_get_q(gr, i, j);
                qgen->zeta[i][j] = gmx_resp_get_zeta(gr, i, j);
            }
            qgen->qsave[i][j]    = qgen->q[i][j];
            qgen->zetasave[i][j] = qgen->zeta[i][j];
        }
    }
    qgen->bAllocSave = TRUE;
}

static void gentop_qgen_get_params(gentop_qgen_t qgen, gmx_resp_t gr)
{
    int i, j;

    if (qgen->bAllocSave)
    {
        for (i = 0; (i < qgen->natom); i++)
        {
            for (j = 0; (j < qgen->nZeta[i]); j++)
            {
                qgen->q[i][j]    = qgen->qsave[i][j];
                qgen->zeta[i][j] = qgen->zetasave[i][j];
                if (NULL != gr)
                {
                    gmx_resp_set_q(gr, i, j, qgen->q[i][j]);
                    gmx_resp_set_zeta(gr, i, j, qgen->zeta[i][j]);
                }
            }
        }
    }
    else
    {
        fprintf(stderr, "WARNING: no ESP charges generated.\n");
    }
}

int gentop_qgen_get_nzeta(gentop_qgen_t qgen, int atom)
{
    if ((0 <= atom) && (atom < qgen->natom))
    {
        return qgen->nZeta[atom];
    }
    else
    {
        return NOTSET;
    }
}

int gentop_qgen_get_row(gentop_qgen_t qgen, int atom, int z)
{
    if ((0 <= atom) && (atom < qgen->natom) &&
        (0 <= z) && (z <= qgen->nZeta[atom]))
    {
        return qgen->row[atom][z];
    }
    else
    {
        return NOTSET;
    }
}

double gentop_qgen_get_q(gentop_qgen_t qgen, int atom, int z)
{
    if ((0 <= atom) && (atom < qgen->natom) &&
        (0 <= z) && (z <= qgen->nZeta[atom]))
    {
        return qgen->q[atom][z];
    }
    else
    {
        return NOTSET;
    }
}

double gentop_qgen_get_zeta(gentop_qgen_t qgen, int atom, int z)
{
    if ((0 <= atom) && (atom < qgen->natom) &&
        (0 <= z) && (z <= qgen->nZeta[atom]))
    {
        return qgen->zeta[atom][z];
    }
    else
    {
        return NOTSET;
    }
}

static real Coulomb_NN(real r)
{
    return 1/r;
}

static real calc_jab(ChargeGenerationModel iModel,
                     rvec xi, rvec xj,
                     int nZi, int nZj,
                     real *zeta_i, real *zeta_j,
                     int *rowi, int *rowj)
{
    int  i, j;
    rvec dx;
    real r;
    real eTot = 0;

    rvec_sub(xi, xj, dx);
    r = norm(dx);
    if (r == 0)
    {
        gmx_fatal(FARGS, "Zero distance between atoms!\n");
    }
    if ((*zeta_i <= 0) || (*zeta_j <= 0))
    {
        iModel = eqgAXp;
    }
    switch (iModel)
    {
        case eqgAXp:
            eTot = Coulomb_NN(r);
            break;
        case eqgAXs:
        case eqgRappe:
        case eqgYang:
            eTot = 0;
            for (i = nZi-1; (i < nZi); i++)
            {
                for (j = nZj-1; (j < nZj); j++)
                {
                    eTot += Coulomb_SS(r, rowi[i], rowj[j], zeta_i[i], zeta_j[j]);
                }
            }
            break;
        case eqgAXg:
            eTot = 0;
            for (i = nZi-1; (i < nZi); i++)
            {
                for (j = nZj-1; (j < nZj); j++)
                {
                    eTot += Coulomb_GG(r, zeta_i[i], zeta_j[j]);
                }
            }
            break;
        default:
            gmx_fatal(FARGS, "Unsupported model %d in calc_jab", iModel);
    }

    return ONE_4PI_EPS0*(eTot)/ELECTRONVOLT;
}

static void solve_q_eem(FILE *fp, gentop_qgen *qgen, real hardness_factor)
{
    double **a, qtot, q;
    int      i, j, n;

    n = qgen->natom+1;
    a = alloc_matrix(n, n);
    for (i = 0; (i < n-1); i++)
    {
        for (j = 0; (j < n-1); j++)
        {
            a[i][j] = qgen->Jab[i][j];
        }
        a[i][i] = hardness_factor*qgen->Jab[i][i];
    }
    for (j = 0; (j < n-1); j++)
    {
        a[n-1][j] = 1;
    }
    for (i = 0; (i < n-1); i++)
    {
        a[i][n-1] = -1;
    }
    a[n-1][n-1] = 0;

    if (matrix_invert(fp, n, a) == 0)
    {
        for (i = 0; (i < n); i++)
        {
            q = 0;
            for (j = 0; (j < n); j++)
            {
                q += a[i][j]*qgen->rhs[j];
            }
            qgen->q[i][qgen->nZeta[i]-1] = q;
            if (fp)
            {
                fprintf(fp, "%2d RHS = %10g Charge= %10g\n", i, qgen->rhs[i], q);
            }
        }
    }
    else
    {
        for (i = 0; (i < n); i++)
        {
            qgen->q[i][qgen->nZeta[i]] = 1;
        }
    }
    qgen->chieq = qgen->q[n-1][qgen->nZeta[n-1]-1];
    qtot        = 0;
    for (i = 0; (i < n-1); i++)
    {
        for (j = 0; (j < qgen->nZeta[i]); j++)
        {
            qtot += qgen->q[i][j];
        }
    }

    if (fp && (fabs(qtot - qgen->qtotal) > 1e-2))
    {
        fprintf(fp, "qtot = %g, it should be %g\n", qtot, qgen->qtotal);
    }
    free_matrix(a);
}

static void qgen_update_J00(gentop_qgen *qgen)
{
    int    i;
    double j0, qq;
    double zetaH = 1.0698;

    for (i = 0; (i < qgen->natom); i++)
    {
        j0 = qgen->j00[i]/qgen->epsr;
        if (((qgen->iModel == eqgYang) ||
             (qgen->iModel == eqgRappe)) &&
            (qgen->atomnr[i] == 1))
        {
            qq = qgen->q[i][qgen->nZeta[i]-1];
            j0 = (1+qq/zetaH)*j0;

            if (debug && (j0 < 0) && !qgen->bWarned)
            {
                fprintf(debug, "WARNING: J00 = %g for atom %d. The equations will be instable.\n", j0, i+1);
                qgen->bWarned = TRUE;
            }
        }
        qgen->Jab[i][i] = (j0 > 0) ? j0 : 0;
    }
}

static void qgen_debug(FILE *fp, gentop_qgen *qgen)
{
    int i, j;

    for (i = 0; (i < qgen->natom); i++)
    {
        fprintf(fp, "QGEN: i: %2d chi0: %8g J0: %8g q:",
                i+1, qgen->chi0[i], qgen->Jab[i][i]);
        for (j = 0; (j < qgen->nZeta[i]); j++)
        {
            fprintf(fp, " %8g", qgen->q[i][j]);
        }
        fprintf(fp, "\n");
    }
    fprintf(fp, "QGEN Jab matrix:\n");
    for (i = 0; (i < qgen->natom); i++)
    {
        for (j = 0; (j <= i); j++)
        {
            fprintf(fp, "  %6.2f", qgen->Jab[i][j]);
        }
        fprintf(fp, "\n");
    }
    fprintf(fp, "\n");
}

static real qgen_calc_Sij(gentop_qgen *qgen, int i, int j)
{
    real dist, dism, Sij = 1.0;
    rvec dx;
    int  l, m, tag;

    rvec_sub(qgen->x[i], qgen->x[j], dx);
    dist = norm(dx);
    if ((dist < 0.118) && (qgen->atomnr[i] != 1) && (qgen->atomnr[j] != 1))
    {
        Sij = Sij*1.64;
    }
    else if ((dist < 0.122) && (qgen->atomnr[i] != 1) && (qgen->atomnr[j] != 1))
    {
        if ((qgen->atomnr[i] != 8) && (qgen->atomnr[j] != 8))
        {
            Sij = Sij*2.23;
        }
        else
        {
            Sij = Sij*2.23;
        }
    }
    else if (dist < 0.125)
    {
        tag = 0;
        if ((qgen->atomnr[i] == 6) && (qgen->atomnr[j] == 8))
        {
            tag = i;
        }
        else if ((qgen->atomnr[i] == 8) && (qgen->atomnr[j] == 6))
        {
            tag = j;
        }
        if (tag != 0)
        {
            printf("found CO\n");
            for (l = 0; (l < qgen->natom); l++)
            {
                if (qgen->atomnr[l] == 1)
                {
                    printf("found H\n");
                    dism = 0.0;
                    for (m = 0; (m < DIM); m++)
                    {
                        dism = dism+sqr(qgen->x[tag][m]-qgen->x[l][m]);
                    }

                    printf("dist: %8.3f\n", sqrt(dism));
                    if (sqrt(dism) < 0.105)
                    {
                        printf("dist %5d %5d %5s  %5s %8.3f\n",
                               i, l, qgen->elem[tag], qgen->elem[l], sqrt(dism));
                        Sij = Sij*1.605;
                    }
                }
            }
        }
    }
    else if ((qgen->atomnr[i] == 6) && (qgen->atomnr[j] == 8))
    {
        Sij = Sij*1.03;
    }
    else if (((qgen->atomnr[j] == 6) && (qgen->atomnr[i] == 7) && (dist < 0.15)) ||
             ((qgen->atomnr[i] == 6) && (qgen->atomnr[j] == 7) && (dist < 0.15)))
    {
        if (qgen->atomnr[i] == 6)
        {
            tag = i;
        }
        else
        {
            tag = j;
        }
        for (l = 0; (l < qgen->natom); l++)
        {
            if (qgen->atomnr[l] == 8)
            {
                printf("found Oxy\n");
                dism = 0.0;
                for (m = 0; (m < DIM); m++)
                {
                    dism = dism+sqr(qgen->x[tag][m]-qgen->x[l][m]);
                }
                if (sqrt(dism) < 0.130)
                {
                    printf("found peptide bond\n");
                    Sij = Sij*0.66;
                }
                else
                {
                    Sij = Sij*1.1;
                }
            }
        }
    }
    return Sij;
}

static void qgen_calc_Jab(gentop_qgen *qgen)
{
    int    i, j;
    double Jab;

    for (i = 0; (i < qgen->natom); i++)
    {
        for (j = i+1; (j < qgen->natom); j++)
        {
            Jab = calc_jab(qgen->iModel,
                           qgen->x[i], qgen->x[j],
                           qgen->nZeta[i], qgen->nZeta[j],
                           qgen->zeta[i], qgen->zeta[j],
                           qgen->row[i], qgen->row[j]);
            if (qgen->iModel == eqgYang)
            {
                Jab = Jab*qgen_calc_Sij(qgen, i, j);
            }
            qgen->Jab[j][i] = qgen->Jab[i][j] = Jab/qgen->epsr;
        }
    }
}

static void qgen_calc_rhs(gentop_qgen *qgen)
{
    int    i, j, k, l;
    rvec   dx;
    real   r, j1, j1q, qcore;

    /* This right hand side is for all models */
    for (i = 0; (i < qgen->natom); i++)
    {
        qgen->rhs[i] = -qgen->chi0[i];
    }
    qgen->rhs[qgen->natom] = qgen->qtotal;

    /* In case the charge is split in nuclear charge and electronic charge
     * we need to add some more stuff. See paper for details.
     */
    for (i = 0; (i < qgen->natom); i++)
    {
        j1q   = 0;
        qcore = 0;
        for (k = 0; (k < qgen->nZeta[i]-1); k++)
        {
            j1q   += qgen->j00[i]*qgen->q[i][k];
            qcore += qgen->q[i][k];
        }
        j1 = 0;
        /* This assignment for k is superfluous because of the previous loop,
         * but if I take it out it will at some stage break the loop below where
         * exactly this value of k is needed.
         */
        k  = qgen->nZeta[i]-1;
        for (j = 0; (j < qgen->natom); j++)
        {
            if (i != j)
            {
                rvec_sub(qgen->x[i], qgen->x[j], dx);
                r = norm(dx);
                switch (qgen->iModel)
                {
                    case eqgAXs:
                    case eqgRappe:
                    case eqgYang:
                        for (l = 0; (l < qgen->nZeta[j]-1); l++)
                        {
                            j1 += qgen->q[j][l]*Coulomb_SS(r, k, l, qgen->zeta[i][k], qgen->zeta[j][l]);
                        }
                        break;
                    case eqgAXg:
                        for (l = 0; (l < qgen->nZeta[j]-1); l++)
                        {
                            j1 += qgen->q[j][l]*Coulomb_GG(r, qgen->zeta[i][k], qgen->zeta[j][l]);
                        }
                        break;
                    default:
                        break;
                }
            }
        }
        qgen->rhs[i]           -= j1q + ONE_4PI_EPS0*j1/ELECTRONVOLT;
        qgen->rhs[qgen->natom] -= qcore;
    }
}

int atomicnumber2rowXX(int elem)
{
    int row;

    /* Compute which row in the periodic table is this element */
    if (elem <= 2)
    {
        row = 1;
    }
    else if (elem <= 10)
    {
        row = 2;
    }
    else if (elem <= 18)
    {
        row = 3;
    }
    else if (elem <= 36)
    {
        row = 4;
    }
    else if (elem <= 54)
    {
        row = 5;
    }
    else if (elem <= 86)
    {
        row = 6;
    }
    else
    {
        row = 7;
    }

    return row;
}

gentop_qgen_t
gentop_qgen_init(gmx_poldata_t pd, t_atoms *atoms, gmx_atomprop_t aps,
                 rvec *x, ChargeGenerationModel iModel, real hfac, int qtotal, real epsr)
{
    gentop_qgen *qgen;
    char        *atp;
    gmx_bool     bSup = TRUE;
    int          i, j, k, atm, nz;

    snew(qgen, 1);
    qgen->iModel = iModel;
    qgen->hfac   = hfac;
    qgen->qtotal = qtotal;
    if (epsr <= 1)
    {
        epsr = 1;
    }
    qgen->epsr   = epsr;
    for (i = j = 0; (i < atoms->nr); i++)
    {
        if (atoms->atom[i].ptype == eptAtom)
        {
            qgen->natom++;
        }
    }
    snew(qgen->chi0, qgen->natom);
    snew(qgen->rhs, qgen->natom+1);
    snew(qgen->elem, qgen->natom);
    snew(qgen->atomnr, qgen->natom);
    snew(qgen->row, qgen->natom);
    snew(qgen->Jab, qgen->natom+1);
    snew(qgen->zeta, qgen->natom);
    snew(qgen->j00, qgen->natom);
    snew(qgen->q, qgen->natom+1);
    snew(qgen->x, qgen->natom);
    qgen->bAllocSave = FALSE;
    snew(qgen->nZeta, qgen->natom+1);
    /* Special case for chi_eq */
    qgen->nZeta[qgen->natom] = 1;
    snew(qgen->q[qgen->natom], qgen->nZeta[qgen->natom]);
    for (i = j = 0; (i < atoms->nr) && bSup; i++)
    {
        if (atoms->atom[i].ptype == eptAtom)
        {
            snew(qgen->Jab[j], qgen->natom+1);
            atm = atoms->atom[i].atomnumber;
            if (atm == NOTSET)
            {
                gmx_fatal(FARGS, "Don't know atomic number for %s %s",
                          *(atoms->resinfo[i].name),
                          *(atoms->atomname[j]));
            }
            atp = *atoms->atomtype[j];
            if (gmx_poldata_have_eem_support(pd, qgen->iModel, atp, TRUE) == 0)
            {
                atp = gmx_atomprop_element(aps, atm);
                if (gmx_poldata_have_eem_support(pd, qgen->iModel, atp, TRUE) == 0)
                {
                    fprintf(stderr, "No charge generation support for atom %s (element %s), model %s\n",
                            *atoms->atomtype[j], atp, get_eemtype_name(qgen->iModel));
                    bSup = FALSE;
                }
            }
            if (bSup)
            {
                qgen->elem[j]   = strdup(atp);
                qgen->atomnr[j] = atm;
                nz              = gmx_poldata_get_nzeta(pd, qgen->iModel, atp);
                qgen->nZeta[j]  = nz;
                snew(qgen->q[j], nz);
                snew(qgen->zeta[j], nz);
                snew(qgen->row[j], nz);
                for (k = 0; (k < nz); k++)
                {
                    qgen->q[j][k]    = gmx_poldata_get_q(pd, qgen->iModel, *atoms->atomtype[j], k);
                    qgen->zeta[j][k] = gmx_poldata_get_zeta(pd, qgen->iModel, *atoms->atomtype[j], k);
                    qgen->row[j][k]  = gmx_poldata_get_row(pd, qgen->iModel, *atoms->atomtype[j], k);
                    if (qgen->row[j][k] > SLATER_MAX)
                    {
                        if (debug)
                        {
                            fprintf(debug, "Can not handle higher slaters than %d for atom %s %s\n",
                                    SLATER_MAX,
                                    *(atoms->resinfo[i].name),
                                    *(atoms->atomname[j]));
                        }
                        qgen->row[j][k] = SLATER_MAX;
                    }
                }
                qgen->chi0[j]  = 0;
                qgen->j00[j]   = 0;
                copy_rvec(x[i], qgen->x[j]);
                j++;
            }
        }
    }
    if (bSup)
    {
        return (gentop_qgen_t) qgen;
    }
    else
    {
        gentop_qgen_done(qgen);
        sfree(qgen);

        return NULL;
    }
}

static void qgen_print(FILE *fp, t_atoms *atoms, gentop_qgen *qgen)
{
    int  i, j, k, m;
    rvec mu = { 0, 0, 0 };
    real qq;

    if (qgen->eQGEN == eQGEN_OK)
    {
        if (fp)
        {
            fprintf(fp, "Res  Atom   Nr       J0     chi0 row        q zeta (1/nm)\n");
        }
        for (i = j = 0; (i < atoms->nr); i++)
        {
            if (atoms->atom[i].ptype == eptAtom)
            {
                qq = 0;
                for (k = 0; (k < qgen->nZeta[j]); k++)
                {
                    qq += qgen->q[j][k];
                }

                atoms->atom[i].q = qq;
                for (m = 0; (m < DIM); m++)
                {
                    mu[m] += qq* qgen->x[i][m] * ENM2DEBYE;
                }
                if (fp)
                {
                    fprintf(fp, "%4s %4s%5d %8g %8g",
                            *(atoms->resinfo[atoms->atom[i].resind].name),
                            *(atoms->atomname[i]), i+1, qgen->j00[j], qgen->chi0[j]);
                    for (k = 0; (k < qgen->nZeta[j]); k++)
                    {
                        fprintf(fp, " %3d %8.5f %8.4f", qgen->row[j][k], qgen->q[j][k],
                                qgen->zeta[j][k]);
                    }
                    fprintf(fp, "\n");
                }
                j++;
            }
        }
        if (fp)
        {
            fprintf(fp, "<chieq> = %10g\n|mu| = %8.3f ( %8.3f  %8.3f  %8.3f )\n",
                    qgen->chieq, norm(mu), mu[XX], mu[YY], mu[ZZ]);
        }
    }
}

void
gentop_qgen_done(gentop_qgen *qgen)
{
    int  i;

    sfree(qgen->chi0);
    sfree(qgen->rhs);
    sfree(qgen->atomnr);
    sfree(qgen->j00);
    sfree(qgen->x);
    for (i = 0; (i < qgen->natom); i++)
    {
        sfree(qgen->row[i]);
        sfree(qgen->q[i]);
        sfree(qgen->zeta[i]);
        sfree(qgen->Jab[i]);
        sfree(qgen->elem[i]);
        if (qgen->bAllocSave)
        {
            sfree(qgen->qsave[i]);
            sfree(qgen->zetasave[i]);
        }
    }
    sfree(qgen->row);
    sfree(qgen->zeta);
    sfree(qgen->elem);
    sfree(qgen->q);
    sfree(qgen->Jab);
    sfree(qgen->nZeta);
    if (qgen->bAllocSave)
    {
        sfree(qgen->qsave);
        sfree(qgen->zetasave);
        qgen->bAllocSave = FALSE;
    }
}

void qgen_message(gentop_qgen_t qgen, int len, char buf[], gmx_resp_t gr)
{
    switch (qgen->eQGEN)
    {
        case eQGEN_OK:
            if (NULL != gr)
            {
                gmx_resp_calc_pot(gr);
                gmx_resp_calc_rms(gr);
                gmx_resp_statistics(gr, len, buf);
            }
            else
            {
                sprintf(buf, "Charge generation finished correctly.\n");
            }
            break;
        case eQGEN_NOTCONVERGED:
            sprintf(buf, "Charge generation did not converge.\n");
            break;
        case eQGEN_NOSUPPORT:
            sprintf(buf, "No charge generation support for (some of) the atomtypes.\n");
            break;
        case eQGEN_ERROR:
        default:
            sprintf(buf, "Unknown status %d in charge generation.\n", qgen->eQGEN);
    }
}

static void qgen_check_support(gentop_qgen *qgen, gmx_poldata_t pd, gmx_atomprop_t aps)
{
    int      i;
    gmx_bool bSup = TRUE;

    for (i = 0; (i < qgen->natom); i++)
    {
        if (gmx_poldata_have_eem_support(pd, qgen->iModel, qgen->elem[i], TRUE) == 0)
        {
            /*sfree(qgen->elem[i]);*/
            qgen->elem[i] = strdup(gmx_atomprop_element(aps, qgen->atomnr[i]));
            if (gmx_poldata_have_eem_support(pd, qgen->iModel, qgen->elem[i], TRUE) == 0)
            {
                fprintf(stderr, "No charge generation support for atom %s, model %s\n",
                        qgen->elem[i], get_eemtype_name(qgen->iModel));
                bSup = FALSE;
            }
        }
    }
    if (bSup)
    {
        qgen->eQGEN = eQGEN_OK;
    }
    else
    {
        qgen->eQGEN = eQGEN_NOSUPPORT;
    }
}

static void qgen_update_pd(t_atoms *atoms, gmx_poldata_t pd, gentop_qgen_t qgen)
{
    int i, j, n, nz;

    for (i = j = 0; (i < atoms->nr); i++)
    {
        if (atoms->atom[i].ptype == eptAtom)
        {
            qgen->chi0[j]  = gmx_poldata_get_chi0(pd, qgen->iModel, qgen->elem[j]);
            qgen->j00[j]   = gmx_poldata_get_j00(pd, qgen->iModel, qgen->elem[j]);
            nz             = gmx_poldata_get_nzeta(pd, qgen->iModel, qgen->elem[j]);
            for (n = 0; (n < nz); n++)
            {
                qgen->zeta[j][n] = gmx_poldata_get_zeta(pd, qgen->iModel, qgen->elem[j], n);
                qgen->q[j][n]    = gmx_poldata_get_q(pd, qgen->iModel, qgen->elem[j], n);
                qgen->row[j][n]  = gmx_poldata_get_row(pd, qgen->iModel, qgen->elem[j], n);
            }
            j++;
        }
    }
}

int generate_charges_sm(FILE *fp,
                        gentop_qgen *qgen, gmx_poldata_t pd,
                        t_atoms *atoms,
                        real tol, int maxiter, gmx_atomprop_t aps,
                        real *chieq)
{
    real       *qq = NULL;
    int         i, j, iter;
    real        rms;

    qgen_check_support(qgen, pd, aps);
    if (eQGEN_OK == qgen->eQGEN)
    {
        qgen_update_pd(atoms, pd, qgen);

        snew(qq, atoms->nr+1);
        for (i = j = 0; (i < atoms->nr); i++)
        {
            if (atoms->atom[i].ptype != eptShell)
            {
                qq[j] = qgen->q[j][qgen->nZeta[j]-1];
                j++;
            }
        }
        iter = 0;
        qgen_calc_Jab(qgen);
        qgen_calc_rhs(qgen);
        do
        {
            qgen_update_J00(qgen);
            if (debug)
            {
                qgen_debug(debug, qgen);
            }
            solve_q_eem(debug, qgen, 1.0);
            rms = 0;
            for (i = j = 0; (i < atoms->nr); i++)
            {
                if (atoms->atom[i].ptype != eptShell)
                {
                    rms  += sqr(qq[j] - qgen->q[j][qgen->nZeta[j]-1]);
                    qq[j] = qgen->q[j][qgen->nZeta[j]-1];
                    j++;
                }
            }
            rms = sqrt(rms/atoms->nr);
            iter++;
        }
        while ((rms > tol) && (iter < maxiter));

        if (iter < maxiter)
        {
            qgen->eQGEN = eQGEN_OK;
        }
        else
        {
            qgen->eQGEN = eQGEN_NOTCONVERGED;
        }

        if (fp)
        {
            if (qgen->eQGEN == eQGEN_OK)
            {
                fprintf(fp, "Converged to tolerance %g after %d iterations\n",
                        tol, iter);
            }
            else
            {
                fprintf(fp, "Did not converge within %d iterations. RMS = %g\n",
                        maxiter, rms);
            }
        }
        *chieq = qgen->chieq;
        sfree(qq);
    }

    if (eQGEN_OK == qgen->eQGEN)
    {
        qgen_print(fp, atoms, qgen);
    }

    return qgen->eQGEN;
}

static int generate_charges_bultinck(FILE *fp,
                                     gentop_qgen_t qgen,
                                     gmx_poldata_t pd, t_atoms *atoms,
                                     gmx_atomprop_t aps)
{
    qgen_check_support(qgen, pd, aps);
    if (eQGEN_OK == qgen->eQGEN)
    {
        qgen_update_pd(atoms, pd, qgen);

        qgen_calc_Jab(qgen);
        qgen_calc_rhs(qgen);
        qgen_update_J00(qgen);
        solve_q_eem(debug, qgen, 2.0);

        qgen_print(fp, atoms, qgen);
    }

    return qgen->eQGEN;
}

int generate_charges(FILE *fp,
                     gentop_qgen_t qgen, gmx_resp_t gr,
                     const char *molname, gmx_poldata_t pd,
                     t_atoms *atoms,
                     real tol, int maxiter, int maxcycle,
                     gmx_atomprop_t aps)
{
    int  cc, eQGEN_min = eQGEN_NOTCONVERGED;
    real chieq, chi2, chi2min = GMX_REAL_MAX;

    /* Generate charges */
    switch (qgen->iModel)
    {
        case eqgRESP:
            if (NULL == gr)
            {
                gmx_incons("No RESP data structure");
            }
            if (fp)
            {
                fprintf(fp, "Generating %s charges for %s using RESP algorithm\n",
                        get_eemtype_name(qgen->iModel), molname);
            }
            for (cc = 0; (cc < maxcycle); cc++)
            {
                if (fp)
                {
                    fprintf(fp, "Cycle %d/%d\n", cc+1, maxcycle);
                }
                /* Fit charges to electrostatic potential */
                qgen->eQGEN = gmx_resp_optimize_charges(fp, gr, maxiter, tol, &chi2);
                if (qgen->eQGEN == eQGEN_OK)
                {
                    eQGEN_min = qgen->eQGEN;
                    if (chi2 <= chi2min)
                    {
                        gentop_qgen_save_params(qgen, gr);
                        chi2min = chi2;
                    }

                    if (NULL != fp)
                    {
                        fprintf(fp, "chi2 = %g kJ/mol e\n", chi2);
                    }
                    qgen_print(fp, atoms, qgen);
                }
            }
            if (maxcycle > 1)
            {
                if (fp)
                {
                    fprintf(fp, "---------------------------------\nchi2 at minimum is %g\n", chi2min);
                }
                gentop_qgen_get_params(qgen, gr);
                qgen_print(fp, atoms, qgen);
            }
            qgen->eQGEN = eQGEN_min;
            break;
        default:
            /* Use empirical algorithms */
            if (fp)
            {
                fprintf(fp, "Generating charges for %s using %s algorithm\n",
                        molname, get_eemtype_name(qgen->iModel));
            }
            if (qgen->iModel == eqgBultinck)
            {
                (void) generate_charges_bultinck(fp, qgen, pd, atoms, aps);
            }
            else
            {
                (void) generate_charges_sm(fp, qgen, pd, atoms, tol, maxiter, aps, &chieq);
            }
            gentop_qgen_save_params(qgen, gr);
    }
    return qgen->eQGEN;
}
