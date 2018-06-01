/*
 * This source file is part of the Alexandria program.
 *
 * Copyright (C) 2014-2018 
 *
 * Developers:
 *             Mohammad Mehdi Ghahremanpour, 
 *             Paul J. van Maaren, 
 *             David van der Spoel (Project leader)
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
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, 
 * Boston, MA  02110-1301, USA.
 */
 
/*! \internal \brief
 * Implements part of the alexandria program.
 * \author Mohammad Mehdi Ghahremanpour <mohammad.ghahremanpour@icm.uu.se>
 * \author David van der Spoel <david.vanderspoel@icm.uu.se>
 */

#include <cctype>

#include "gromacs/coulombintegrals/coulombintegrals.h"
#include "gromacs/fileio/confio.h"
#include "gromacs/linearalgebra/matrix.h"
#include "gromacs/listed-forces/bonded.h"
#include "gromacs/math/units.h"
#include "gromacs/math/vec.h"
#include "gromacs/topology/atomprop.h"
#include "gromacs/topology/atoms.h"

#include "molprop.h"
#include "poldata.h"
#include "qgen_eem.h"
#include "regression.h"


namespace alexandria
{

void QgenEem::setInfo(const Poldata            &pd,
                      t_atoms                  *atoms,
                      ChargeDistributionModel   iChargeDistributionModel,
                      double                    hfac,
                      int                       qtotal,
                      bool                      haveShell)
{
    int          i, j, k, atm, nz;
    bool         bSupport = true;
    std::string  atp;
    
    bWarned_                   = false;
    bAllocSave_                = false;
    bHaveShell_                = haveShell;
    eQGEN_                     = eQGEN_OK;
    hardnessFactor_            = 1;
    chieq_                     = 0;
    Jcs_                       = 0;
    Jss_                       = 0;
    rms_                       = 0;
    natom_                     = 0;
    iChargeDistributionModel_  = iChargeDistributionModel;
    hfac_                      = hfac;
    qtotal_                    = qtotal;
   
    for (i = 0; i < atoms->nr; i++)
    {
        if (atoms->atom[i].ptype == eptAtom)
        {
            coreIndex_.push_back(i);
            natom_++;
        }
        else if (atoms->atom[i].ptype == eptShell)
        {
            shellIndex_.push_back(i);
        }
    }

    chi0_.resize(natom_, 0);
    rhs_.resize(natom_ + 1, 0);
    elem_.resize(natom_);
    atomnr_.resize(natom_, 0);
    row_.resize(natom_);
    Jcc_.resize(natom_ + 1);
    zeta_.resize(natom_);
    j00_.resize(natom_, 0);
    q_.resize(natom_ + 1);
    nZeta_.resize(natom_ + 1, 0);
    x_.resize(atoms->nr);

    /* Special case for chi_eq */
    nZeta_[natom_] = 1;
    q_[natom_].resize(nZeta_[natom_], 0);

    for (i = j = 0; (i < atoms->nr) && bSupport; i++)
    {
        if (atoms->atom[i].ptype == eptAtom)
        {
            Jcc_[j].resize(natom_ + 1, 0);
            atm = atoms->atom[i].atomnumber;
            if (atm == 0)
            {
                gmx_fatal(FARGS, "Don't know atomic number for %s %s",
                          *(atoms->resinfo[i].name),
                          *(atoms->atomname[j]));
            }
            atp.assign(*atoms->atomtype[i]);
            if (!pd.haveEemSupport(iChargeDistributionModel_, atp, true))
            {
                fprintf(stderr, "No charge distribution support for atom %s, model %s\n",
                        *atoms->atomtype[i],
                        getEemtypeName(iChargeDistributionModel_));
                bSupport = false;
            }
            if (bSupport)
            {
                auto eem   = pd.findEem(iChargeDistributionModel_, atp);                
                elem_[j]   = atp;
                atomnr_[j] = atm;                
                chi0_[j]   = eem->getChi0();
                j00_[j]    = eem->getJ0();
                nz         = eem->getNzeta();
                nZeta_[j]  = nz;
                q_[j].resize(nz, 0);
                zeta_[j].resize(nz, 0);
                row_[j].resize(nz, 0);
                for (k = 0; k < nz; k++)
                {
                    q_[j][k]    = eem->getQ(k);
                    zeta_[j][k] = eem->getZeta(k);
                    row_[j][k]  = eem->getRow(k);
                    char buf[256];
                    snprintf(buf, sizeof(buf), "Row (in the periodic table) should be at least 1. Here: atype = %s q = %g zeta = %g row = %d model = %s",
                             atp.c_str(), q_[j][k], zeta_[j][k], row_[j][k],
                             getEemtypeName(iChargeDistributionModel_));
                    GMX_RELEASE_ASSERT(iChargeDistributionModel == eqdAXp || row_[j][k] != 0, buf);
                    #if HAVE_LIBCLN
                    if (row_[j][k] > SLATER_MAX_CLN)
                    {
                        if (debug)
                        {
                            fprintf(debug, "Cannot handle higher slaters than %d for atom %s %s\n",
                                    SLATER_MAX_CLN,
                                    *(atoms->resinfo[i].name),
                                    *(atoms->atomname[j]));
                        }
                        row_[j][k] = SLATER_MAX_CLN;
                    }
                    #else
                    if (row_[j][k] > SLATER_MAX)
                    {
                        if (debug)
                        {
                            fprintf(debug, "Cannot handle higher slaters than %d for atom %s %s\n",
                                    SLATER_MAX,
                                    *(atoms->resinfo[i].name),
                                    *(atoms->atomname[j]));
                        }
                        row_[j][k] = SLATER_MAX;
                    }
                    #endif
                }
                j++;
            }
        }
    }
}

void QgenEem::updateInfo(const Poldata &pd)
{
    for (auto i = 0; i < natom_; i++)
    {
        auto ei  = pd.findEem(iChargeDistributionModel_, elem_[i]);
        chi0_[i] = ei->getChi0();
        j00_[i]  = ei->getJ0();
        for (auto k = 0; k < nZeta_[i]; k++)
        {
            zeta_[i][k] = ei->getZeta(k);
        }
    }
}

int QgenEem::getNzeta( int atom)
{
    if ((0 <= atom) && (atom < natom_))
    {
        return nZeta_[atom];
    }
    return 0;
}

double QgenEem::getQ(int atom, int z)
{
    if ((0 <= atom) && (atom < natom_) &&
        (0 <= z) && (z <= nZeta_[atom]))
    {
        return q_[atom][z];
    }
    return 0;
}

int QgenEem::getRow( int atom, int z)
{
    if ((0 <= atom) && (atom < natom_) &&
        (0 <= z) && (z <= nZeta_[atom]))
    {
        return row_[atom][z];
    }
    return 0;

}

double QgenEem::getZeta(int atom, int z)
{
    if ((0 <= atom) && (atom < natom_) &&
        (0 <= z) && (z <= nZeta_[atom]))
    {
        return zeta_[atom][z];
    }
    return 0;
}

double Coulomb_PP(double r)
{
    return 1/r;
}

void QgenEem::debugFun(FILE *fp)
{
    auto i = 0, j = 0;

    for (i = 0; i < natom_; i++)
    {
        fprintf(fp, "THIS: i: %2d _chi0: %8g J0: %8g q:",
                i+1, chi0_[i], Jcc_[i][i]);
        for (j = 0; j < nZeta_[i]; j++)
        {
            fprintf(fp, "Charge: %8g  Zeta: %8g", q_[i][j], zeta_[i][j]);
        }
        fprintf(fp, "\n");
    }
    fprintf(fp, "qgen _Jcc matrix:\n");
    for (i = 0; i < natom_; i++)
    {
        for (j = 0; j <= i; j++)
        {
            fprintf(fp, "  %6.2f", Jcc_[i][j]);
        }
        fprintf(fp, "\n");
    }
    fprintf(fp, "\n");
}

double QgenEem::calcSij(int i, int j)
{
    double dist, dism, Sij = 1.0;
    rvec   dx;
    int    l, m, tag;

    rvec_sub(x_[i], x_[j], dx);
    dist = norm(dx);
    if ((dist < 0.118) && (atomnr_[i] != 1) && (atomnr_[j] != 1))
    {
        Sij = Sij*1.64;
    }
    else if ((dist < 0.122) && (atomnr_[i] != 1) && (atomnr_[j] != 1))
    {
        if ((atomnr_[i] != 8) && (atomnr_[j] != 8))
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
        if ((atomnr_[i] == 6) && (atomnr_[j] == 8))
        {
            tag = i;
        }
        else if ((atomnr_[i] == 8) && (atomnr_[j] == 6))
        {
            tag = j;
        }
        if (tag != 0)
        {
            printf("found CO\n");
            for (l = 0; (l < natom_); l++)
            {
                if (atomnr_[l] == 1)
                {
                    printf("found H\n");
                    dism = 0.0;
                    for (m = 0; (m < DIM); m++)
                    {
                        dism = dism+gmx::square(x_[tag][m] - x_[l][m]);
                    }

                    printf("dist: %8.3f\n", sqrt(dism));
                    if (sqrt(dism) < 0.105)
                    {
                        printf("dist %5d %5d %5s  %5s %8.3f\n",
                               i, l, elem_[tag].c_str(), elem_[l].c_str(), sqrt(dism));
                        Sij = Sij*1.605;
                    }
                }
            }
        }
    }
    else if ((atomnr_[i] == 6) && (atomnr_[j] == 8))
    {
        Sij = Sij*1.03;
    }
    else if (((atomnr_[j] == 6) && (atomnr_[i] == 7) && (dist < 0.15)) ||
             ((atomnr_[i] == 6) && (atomnr_[j] == 7) && (dist < 0.15)))
    {
        if (atomnr_[i] == 6)
        {
            tag = i;
        }
        else
        {
            tag = j;
        }
        for (l = 0; (l < natom_); l++)
        {
            if (atomnr_[l] == 8)
            {
                printf("found Oxy\n");
                dism = 0.0;
                for (m = 0; (m < DIM); m++)
                {
                    dism = dism+gmx::square(x_[tag][m] - x_[l][m]);
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

double QgenEem::calcJ(ChargeDistributionModel iChargeDistributionModel,
                      rvec                    xI, 
                      rvec                    xJ,
                      double                  zetaI, 
                      double                  zetaJ,
                      int                     rowI, 
                      int                     rowJ)
{
    rvec   dx;
    double r    = 0;
    double eTot = 0;

    rvec_sub(xI, xJ, dx);
    r = norm(dx);
    if (r == 0)
    {
        gmx_fatal(FARGS, "Zero distance between atoms!\n");
    }
    if ((zetaI <= 0) || (zetaJ <= 0))
    {
        iChargeDistributionModel = eqdAXp;
    }
    switch (iChargeDistributionModel)
    {
        case eqdAXp:
        case eqdAXpp:
            eTot = Coulomb_PP(r);
            break;
        case eqdAXs:
        case eqdAXps:
        case eqdRappe:
        case eqdYang:        
            eTot = Coulomb_SS(r, rowI, rowJ, zetaI, zetaJ);
            break;
        case eqdAXg:
        case eqdAXpg:
            eTot = Coulomb_GG(r, zetaI, zetaJ);
            break;
        default:
            gmx_fatal(FARGS, "Unsupported model %d in calc_jcc", iChargeDistributionModel);
    }
    return (ONE_4PI_EPS0*eTot)/ELECTRONVOLT;
}

void QgenEem::calcJcc(t_atoms *atoms)
{
    auto Jcc = 0.0;
    auto i   = 0;
    for (auto n = 0; n < atoms->nr; n++)
    {
        if (atoms->atom[n].ptype == eptAtom)
        {
            auto j = i;
            for (auto m = n; m < atoms->nr; m++)
            {
                if (atoms->atom[m].ptype == eptAtom)
                {
                    if (n != m)
                    {
                        j++;
                        Jcc = calcJ(iChargeDistributionModel_,
                                    x_[n], 
                                    x_[m],
                                    zeta_[i][0], 
                                    zeta_[j][0],
                                    row_[i][0], 
                                    row_[j][0]);
                        if (iChargeDistributionModel_ == eqdYang)
                        {
                            Jcc *= calcSij(i, j);
                        }
                        Jcc_[i][j] = Jcc_[j][i] = Jcc;
                    }
                    else
                    {
                        auto j0 = j00_[i];
                        if (((iChargeDistributionModel_ == eqdYang) ||
                             (iChargeDistributionModel_ == eqdRappe)) &&
                            (atomnr_[i] == 1))
                        {
                            auto zetaH = 1.0698;
                            j0 = (1+q_[i][0]/zetaH)*j0;                            
                            if (j0 < 0 && !bWarned_)
                            {
                                bWarned_ = true;
                            }
                        }
                        Jcc_[i][i] = (j0 > 0) ? j0 : 0;
                    }
                }                
            }
            i++;
        }
    }
}

void QgenEem::calcJcs(t_atoms *atoms,
                      int      top_ndx,
                      int      eem_ndx)
{
    auto l   = 0;
    auto k   = 0;
    auto Jcs = 0.0;
    Jcs_     = 0;
    if (atoms->atom[top_ndx].ptype == eptAtom)
    {
        auto itsShell = top_ndx + 1;
        for (l = k = 0; l < atoms->nr; l++)
        {
            if (atoms->atom[l].ptype == eptShell && (l != itsShell))
            {
                k++;
                Jcs = calcJ(iChargeDistributionModel_,
                            x_[top_ndx], 
                            x_[l],
                            zeta_[eem_ndx][0], 
                            zeta_[k][1],
                            row_[eem_ndx][0], 
                            row_[k][1]);
                if (iChargeDistributionModel_ == eqdYang)
                {
                    Jcs *= calcSij(eem_ndx, k);
                }
                Jcs   *= q_[k][1];
                Jcs_  += Jcs;
            }
        }
        Jcs_ *= 2;
    }
    else
    {
        gmx_fatal(FARGS, "atom %d must be eptAtom, but it is not\n", top_ndx);
    }
}

void QgenEem::calcJss(t_atoms *atoms,
                      int      top_ndx,
                      int      eem_ndx)
{
    auto l     = 0;
    auto k     = 0;
    auto Jss   = 0.0;
    Jss_       = 0;
    if (atoms->atom[top_ndx].ptype == eptShell)
    {
        for (l = k = 0; l < atoms->nr; l++)
        {
            if (atoms->atom[l].ptype == eptShell && (l != top_ndx))
            {
                k++;
                Jss = calcJ(iChargeDistributionModel_,
                            x_[top_ndx], 
                            x_[l],
                            zeta_[eem_ndx][1], 
                            zeta_[k][1],
                            row_[eem_ndx][1], 
                            row_[k][1]);
                if (iChargeDistributionModel_ == eqdYang)
                {
                    Jss *= calcSij(eem_ndx, k);
                }
                Jss  *= q_[k][1];
                Jss_ += Jss;
            }
        }
        Jss_ *= 0.5;
    }
    else
    {
        gmx_fatal(FARGS, "atom %d must be eptShell, but it is not\n", top_ndx);
    }
}

void QgenEem::calcRhs(t_atoms *atoms)
{
    auto   qcore     = 0.0;
    auto   qshell    = 0.0;
      
    for (auto i = 0; i < natom_; i++)
    {   
        rhs_[i]   = 0;
        rhs_[i]  -= chi0_[i];
        if (bHaveShell_)
        {
            calcJcs(atoms, coreIndex_[i], i);
            calcJss(atoms, shellIndex_[i], i);
            rhs_[i]   -= j00_[i]*hardnessFactor_*q_[i][1];
            rhs_[i]   -= Jcs_;
            //rhs_[i]   -= Jss_;
            qshell    += q_[i][1];
        }       
    }    
    qcore = qtotal_ - qshell;
    rhs_[natom_] = qcore;
}

void QgenEem::copyChargesToAtoms(t_atoms *atoms)
{
    if (bHaveShell_)
    {
        auto j = 0;
        for (auto i = j = 0; i < atoms->nr; i++)
        {
            if (atoms->atom[i].ptype == eptAtom)
            {
                atoms->atom[i].q = atoms->atom[i].qB = q_[j][0];
            }
            else if (atoms->atom[i].ptype == eptShell)
            {
                atoms->atom[i].q = atoms->atom[i].qB = q_[j][1];
                j++;
            }
        }
    }
    else
    {
        for (auto i = 0; i < atoms->nr; i++)
        {
            atoms->atom[i].q = atoms->atom[i].qB = q_[i][0];
        }
    }
}

void QgenEem::updatePositions(PaddedRVecVector x,
                           t_atoms         *atoms)
{
    for (auto i = 0; i < atoms->nr; i++)
    {
        copy_rvec(x[i], x_[i]);
    }
}

void QgenEem::print(FILE *fp, t_atoms *atoms)
{
    auto  i = 0, j = 0;
    rvec mu = { 0, 0, 0 };

    if (eQGEN_ == eQGEN_OK)
    {
        if (fp)
        {
            fprintf(fp, "                                          Core                 Shell\n");
            fprintf(fp, "Res  Atom   Nr       J0    _chi0  row   q   zeta        row    q     zeta\n");
        }
        for (i = j = 0; i < atoms->nr; i++)
        {
            for (auto m = 0; m < DIM; m++)
            {
                mu[m] += atoms->atom[i].q * x_[i][m] * ENM2DEBYE;
            }
            if (atoms->atom[i].ptype == eptAtom)
            {
                if (fp)
                {
                    fprintf(fp, "%4s %4s%5d %8g %8g",
                            *(atoms->resinfo[atoms->atom[i].resind].name),
                            *(atoms->atomname[i]), i+1, j00_[j], chi0_[j]);
                    for (auto k = 0; k < nZeta_[j]; k++)
                    {
                        fprintf(fp, " %3d %8.5f %8.4f", row_[j][k], q_[j][k], zeta_[j][k]);
                    }
                    fprintf(fp, "\n");
                }
                j++;
            }
        }
        if (fp)
        {
            fprintf(fp, "chieq = %10g\n|mu| = %8.3f ( %8.3f  %8.3f  %8.3f )\n",
                    chieq_, norm(mu), mu[XX], mu[YY], mu[ZZ]);
        }
    }
}

const char *QgenEem::message() const
{
    switch (eQGEN_)
    {
        case eQGEN_OK:
            return "Charge generation finished correctly";
        case eQGEN_NOTCONVERGED:
            return "Charge generation did not converge.";
        case eQGEN_NOSUPPORT:
            return "No charge generation support for (some of) the atomtypes.";
        case eQGEN_ERROR:
        default:
            return "Unknown status %d in charge generation";
    }
    return nullptr;
}

void QgenEem::checkSupport(const Poldata &pd)
{
    bool bSupport = true;
    for (auto i = 0; i < natom_; i++)
    {
        if (!pd.haveEemSupport(iChargeDistributionModel_, elem_[i].c_str(), true))
        {
            fprintf(stderr, "No charge generation support for atom %s, model %s\n",
                    elem_[i].c_str(), getEemtypeName(iChargeDistributionModel_));
            bSupport = false;
        }
    }
    if (bSupport)
    {
        eQGEN_ = eQGEN_OK;
    }
    else
    {
        eQGEN_ = eQGEN_NOSUPPORT;
    }
}

void QgenEem::solveQEem(FILE *fp)
{
    double **lhs, qtot, q;
    int      i, j, n;

    n = natom_ + 1;
    lhs = alloc_matrix(n, n);
    for (i = 0; i < natom_; i++)
    {
        for (j = 0; j < natom_; j++)
        {
            lhs[i][j] = Jcc_[i][j];
        }
        lhs[i][i] = hardnessFactor_*Jcc_[i][i];
    }
    for (j = 0; j < natom_; j++)
    {
        lhs[natom_][j] = 1;
    }
    for (i = 0; i < natom_; i++)
    {
        lhs[i][natom_] = -1;
    }    
    lhs[natom_][natom_] = 0;
    if (matrix_invert(fp, n, lhs) == 0)
    {
        for (i = 0; i < n; i++)
        {
            q = 0;
            for (j = 0; j < n; j++)
            {
                q += lhs[i][j]*rhs_[j];
            }
            q_[i][0] = q;
            if (fp)
            {
                fprintf(fp, "%2d _RHS = %10g Charge= %10g\n", i, rhs_[i], q);
            }
        }
    }
    else
    {
        for (i = 0; i < n; i++)
        {
            q_[i][0] = 1;
        }
    }
    chieq_      = q_[natom_][0];
    
    qtot        = 0;
    if (bHaveShell_)
    {
        for (i = 0; i < natom_; i++)
        {
            for (j = 0; j < nZeta_[i]; j++)
            {
                qtot += q_[i][j];
            }
        }
    }
    else
    {
        for (i = 0; i < natom_; i++)
        {
            qtot += q_[i][0];
        }
    }

    if (fp && (fabs(qtot - qtotal_) > 1e-2))
    {
        fprintf(fp, "qtot = %g, it should be %g\n", qtot, qtotal_);
    }
    free_matrix(lhs);
}

int QgenEem::generateChargesSm(FILE              *fp,
                               const Poldata     &pd,
                               t_atoms           *atoms,
                               double            *chieq,
                               PaddedRVecVector   x)
{
    checkSupport(pd);
    if (eQGEN_OK == eQGEN_)
    {
        updateInfo(pd);
        updatePositions(x, atoms);
        calcJcc(atoms);
        calcRhs(atoms);     
        if (fp)
        {
            debugFun(fp);
        }
        solveQEem(fp);
        *chieq = chieq_;
        copyChargesToAtoms(atoms);
        print(fp, atoms);
    }
    return eQGEN_;   
}

int QgenEem::generateChargesBultinck(FILE              *fp,
                                     const Poldata     &pd,
                                     t_atoms           *atoms,
                                     PaddedRVecVector   x)
{
    checkSupport(pd);
    if (eQGEN_OK == eQGEN_)
    {
        updatePositions(x, atoms);
        calcJcc(atoms);
        calcRhs(atoms);
        if (fp)
        {
            debugFun(fp);
        }
        solveQEem(fp);
        copyChargesToAtoms(atoms);
        print(fp, atoms);
    }
    return eQGEN_;
}

int QgenEem::generateCharges(FILE              *fp,
                             const std::string  molname,
                             const Poldata     &pd,
                             t_atoms           *atoms,
                             PaddedRVecVector   x)
{
    auto chieq = 0.0;
    if (fp)
    {
        fprintf(fp, "Generating charges for %s using %s algorithm\n",
                molname.c_str(),
                getEemtypeName(iChargeDistributionModel_));
    }
    if (iChargeDistributionModel_ == eqdBultinck)
    {
        (void) generateChargesBultinck(fp, pd, atoms, x);
    }
    else
    {
        generateChargesSm(fp, pd, atoms, &chieq, x);
    }
    return eQGEN_;
}

}
