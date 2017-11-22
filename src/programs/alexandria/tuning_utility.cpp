/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2016, by the GROMACS development team, led by
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
/*! \internal \brief
 * Implements part of the alexandria program.
 * \author Mohammad Mehdi Ghahremanpour <m.ghahremanpour@hotmail.com>
 */

#include <cctype>
#include <cmath>
#include <cstdio>
#include <cstdlib>
 
#include "tuning_utility.h"
 
namespace alexandria
{

void print_stats(FILE        *fp,        
                 const char  *prop, 
                 gmx_stats_t  lsq, 
                 gmx_bool     bHeader,
                 char        *xaxis,     
                 char        *yaxis)
{
    real a    = 0, da  = 0, b    = 0, db   = 0;
    real mse  = 0, mae = 0, chi2 = 0, rmsd = 0;
    real Rfit = 0;
    int  n;

    if (bHeader)
    {
        fprintf(fp, "Fitting data to y = ax + b, where x = %s and y = %s\n", xaxis, yaxis);
        fprintf(fp, "%-12s %5s %13s %13s %8s %8s %8s %8s\n",
                "Property", "N", "a", "b", "R(%)", "RMSD", "MSE", "MAE");
        fprintf(fp, "---------------------------------------------------------------\n");
    }
    gmx_stats_get_ab(lsq, elsqWEIGHT_NONE, &a, &b, &da, &db, &chi2, &Rfit);
    gmx_stats_get_rmsd(lsq,    &rmsd);
    gmx_stats_get_mse_mae(lsq, &mse, &mae);
    gmx_stats_get_npoints(lsq, &n);
    fprintf(fp, "%-12s %5d %6.3f(%5.3f) %6.3f(%5.3f) %7.2f %8.4f %8.4f %8.4f\n",
            prop, n, a, da, b, db, Rfit*100, rmsd, mse, mae);
}

void print_lsq_set(FILE *fp, gmx_stats_t lsq)
{
    real   x, y;

    fprintf(fp, "@type xy\n");
    while (gmx_stats_get_point(lsq, &x, &y, nullptr, nullptr, 0) == estatsOK)
    {
        fprintf(fp, "%10g  %10g\n", x, y);
    }
    fprintf(fp, "&\n");
}

void xvgr_symbolize(FILE                   *xvgf, 
                    int                     nsym, 
                    const char             *leg[],
                    const gmx_output_env_t *oenv)
{
    int i;

    xvgr_legend(xvgf, nsym, leg, oenv);
    for (i = 0; (i < nsym); i++)
    {
        xvgr_line_props(xvgf, i, elNone, ecBlack+i, oenv);
        fprintf(xvgf, "@ s%d symbol %d\n", i, i+1);
    }
}

void print_polarizability(FILE              *fp, 
                          alexandria::MyMol *mol,
                          char              *calc_name,
                          real               q_toler)
{
    tensor dalpha;
    real   delta = 0;
    
    if (nullptr != calc_name)
    {
        if (strcmp(calc_name, (char *)"Calc") == 0)
        {
            m_sub(mol->alpha_elec_, mol->alpha_calc_, dalpha);
            delta = sqrt(gmx::square(dalpha[XX][XX])+gmx::square(dalpha[XX][YY])+gmx::square(dalpha[XX][ZZ])+
                         gmx::square(dalpha[YY][YY])+gmx::square(dalpha[YY][ZZ]));
            fprintf(fp,
                    "%-4s (%6.2f %6.2f %6.2f) Dev: (%6.2f %6.2f %6.2f) Delta: %6.2f %s\n"
                    "     (%6s %6.2f %6.2f)      (%6s %6.2f %6.2f)\n"
                    "     (%6s %6s %6.2f)      (%6s %6s %6.2f)\n",
                    calc_name,
                    mol->alpha_calc_[XX][XX], mol->alpha_calc_[XX][YY], mol->alpha_calc_[XX][ZZ],
                    dalpha[XX][XX], dalpha[XX][YY], dalpha[XX][ZZ], delta, (delta > q_toler) ? "YYY" : "",
                    "", mol->alpha_calc_[YY][YY], mol->alpha_calc_[YY][ZZ],
                    "", dalpha[YY][YY], dalpha[YY][ZZ],
                    "", "", mol->alpha_calc_[ZZ][ZZ],
                    "", "", dalpha[ZZ][ZZ]);
        }
        else
        {
            fprintf(fp, "Polarizability analysis\n");
            fprintf(fp,
                    "Electronic   (%6.2f %6.2f %6.2f)\n"
                    "             (%6s %6.2f %6.2f)\n"
                    "             (%6s %6s %6.2f)\n",
                    mol->alpha_elec_[XX][XX], mol->alpha_elec_[XX][YY], mol->alpha_elec_[XX][ZZ],
                    "", mol->alpha_elec_[YY][YY], mol->alpha_elec_[YY][ZZ],
                    "", "", mol->alpha_elec_[ZZ][ZZ]);
        }
    }
}

void print_quadrapole(FILE              *fp, 
                      alexandria::MyMol *mol,
                      char              *calc_name,
                      real               q_toler)
{
    tensor dQ;
    real   delta = 0;
    
    if (nullptr != calc_name)
    {
        if (strcmp(calc_name, (char *)"Calc") == 0)
        {
            m_sub(mol->Q_elec_, mol->Q_calc_, dQ);
            delta = sqrt(gmx::square(dQ[XX][XX])+gmx::square(dQ[XX][YY])+gmx::square(dQ[XX][ZZ])+
                         gmx::square(dQ[YY][YY])+gmx::square(dQ[YY][ZZ]));
            fprintf(fp,
                    "%-4s (%6.2f %6.2f %6.2f) Dev: (%6.2f %6.2f %6.2f) Delta: %6.2f %s\n"
                    "     (%6s %6.2f %6.2f)      (%6s %6.2f %6.2f)\n"
                    "     (%6s %6s %6.2f)      (%6s %6s %6.2f)\n",
                    calc_name,
                    mol->Q_calc_[XX][XX], mol->Q_calc_[XX][YY], mol->Q_calc_[XX][ZZ],
                    dQ[XX][XX], dQ[XX][YY], dQ[XX][ZZ], delta, (delta > q_toler) ? "YYY" : "",
                    "", mol->Q_calc_[YY][YY], mol->Q_calc_[YY][ZZ],
                    "", dQ[YY][YY], dQ[YY][ZZ],
                    "", "", mol->Q_calc_[ZZ][ZZ],
                    "", "", dQ[ZZ][ZZ]);
        }
        else if (strcmp(calc_name, (char *)"ESP") == 0)
        {
            m_sub(mol->Q_elec_, mol->Q_esp_, dQ);
            delta = sqrt(gmx::square(dQ[XX][XX])+gmx::square(dQ[XX][YY])+gmx::square(dQ[XX][ZZ])+
                         gmx::square(dQ[YY][YY])+gmx::square(dQ[YY][ZZ]));
            fprintf(fp,
                    "%-4s (%6.2f %6.2f %6.2f) Dev: (%6.2f %6.2f %6.2f) Delta: %6.2f %s\n"
                    "     (%6s %6.2f %6.2f)      (%6s %6.2f %6.2f)\n"
                    "     (%6s %6s %6.2f)      (%6s %6s %6.2f)\n",
                    calc_name,
                    mol->Q_esp_[XX][XX], mol->Q_esp_[XX][YY], mol->Q_esp_[XX][ZZ],
                    dQ[XX][XX], dQ[XX][YY], dQ[XX][ZZ], delta, (delta > q_toler) ? "YYY" : "",
                    "", mol->Q_esp_[YY][YY], mol->Q_esp_[YY][ZZ],
                    "", dQ[YY][YY], dQ[YY][ZZ],
                    "", "", mol->Q_esp_[ZZ][ZZ],
                    "", "", dQ[ZZ][ZZ]);
        }
        else if (strcmp(calc_name, (char *)"MPA") == 0)
        {
            m_sub(mol->Q_elec_, mol->Q_mulliken_, dQ);
            delta = sqrt(gmx::square(dQ[XX][XX])+gmx::square(dQ[XX][YY])+gmx::square(dQ[XX][ZZ])+
                         gmx::square(dQ[YY][YY])+gmx::square(dQ[YY][ZZ]));
            fprintf(fp,
                    "%-4s (%6.2f %6.2f %6.2f) Dev: (%6.2f %6.2f %6.2f) Delta: %6.2f %s\n"
                    "     (%6s %6.2f %6.2f)      (%6s %6.2f %6.2f)\n"
                    "     (%6s %6s %6.2f)      (%6s %6s %6.2f)\n",
                    calc_name,
                    mol->Q_mulliken_[XX][XX], mol->Q_mulliken_[XX][YY], mol->Q_mulliken_[XX][ZZ],
                    dQ[XX][XX], dQ[XX][YY], dQ[XX][ZZ], delta, (delta > q_toler) ? "YYY" : "",
                    "", mol->Q_mulliken_[YY][YY], mol->Q_mulliken_[YY][ZZ],
                    "", dQ[YY][YY], dQ[YY][ZZ],
                    "", "", mol->Q_mulliken_[ZZ][ZZ],
                    "", "", dQ[ZZ][ZZ]);
        }
        else if (strcmp(calc_name, (char *)"HPA") == 0)
        {
            m_sub(mol->Q_elec_, mol->Q_hirshfeld_, dQ);
            delta = sqrt(gmx::square(dQ[XX][XX])+gmx::square(dQ[XX][YY])+gmx::square(dQ[XX][ZZ])+
                         gmx::square(dQ[YY][YY])+gmx::square(dQ[YY][ZZ]));
            fprintf(fp,
                    "%-4s (%6.2f %6.2f %6.2f) Dev: (%6.2f %6.2f %6.2f) Delta: %6.2f %s\n"
                    "     (%6s %6.2f %6.2f)      (%6s %6.2f %6.2f)\n"
                    "     (%6s %6s %6.2f)      (%6s %6s %6.2f)\n",
                    calc_name,
                    mol->Q_hirshfeld_[XX][XX], mol->Q_hirshfeld_[XX][YY], mol->Q_hirshfeld_[XX][ZZ],
                    dQ[XX][XX], dQ[XX][YY], dQ[XX][ZZ], delta, (delta > q_toler) ? "YYY" : "",
                    "", mol->Q_hirshfeld_[YY][YY], mol->Q_hirshfeld_[YY][ZZ],
                    "", dQ[YY][YY], dQ[YY][ZZ],
                    "", "", mol->Q_hirshfeld_[ZZ][ZZ],
                    "", "", dQ[ZZ][ZZ]);
        }
        else if (strcmp(calc_name, (char *)"CM5") == 0)
        {
            m_sub(mol->Q_elec_, mol->Q_cm5_, dQ);
            delta = sqrt(gmx::square(dQ[XX][XX])+gmx::square(dQ[XX][YY])+gmx::square(dQ[XX][ZZ])+
                         gmx::square(dQ[YY][YY])+gmx::square(dQ[YY][ZZ]));
            fprintf(fp,
                    "%-4s (%6.2f %6.2f %6.2f) Dev: (%6.2f %6.2f %6.2f) Delta: %6.2f %s\n"
                    "     (%6s %6.2f %6.2f)      (%6s %6.2f %6.2f)\n"
                    "     (%6s %6s %6.2f)      (%6s %6s %6.2f)\n",
                    calc_name,
                    mol->Q_cm5_[XX][XX], mol->Q_cm5_[XX][YY], mol->Q_cm5_[XX][ZZ],
                    dQ[XX][XX], dQ[XX][YY], dQ[XX][ZZ], delta, (delta > q_toler) ? "YYY" : "",
                    "", mol->Q_cm5_[YY][YY], mol->Q_cm5_[YY][ZZ],
                    "", dQ[YY][YY], dQ[YY][ZZ],
                    "", "", mol->Q_cm5_[ZZ][ZZ],
                    "", "", dQ[ZZ][ZZ]);
        }
        else
        {
            fprintf(fp, "Quadrupole analysis (6 independent components only)\n");
            fprintf(fp,
                    "Electronic   (%6.2f %6.2f %6.2f)\n"
                    "             (%6s %6.2f %6.2f)\n"
                    "             (%6s %6s %6.2f)\n",
                    mol->Q_elec_[XX][XX], mol->Q_elec_[XX][YY], mol->Q_elec_[XX][ZZ],
                    "", mol->Q_elec_[YY][YY], mol->Q_elec_[YY][ZZ],
                    "", "", mol->Q_elec_[ZZ][ZZ]);
        }
    }
}

void print_dipole(FILE              *fp, 
                  alexandria::MyMol *mol, 
                  char              *calc_name, 
                  real               toler)
{
    rvec dmu;
    real ndmu, cosa;
    char ebuf[32];

    if (nullptr != calc_name)
    {
        if (strcmp(calc_name, (char *)"Calc") == 0)
        {
            rvec_sub(mol->mu_elec_, mol->mu_calc_, dmu);
            ndmu = norm(dmu);
            cosa = cos_angle(mol->mu_elec_, mol->mu_calc_);
            if (ndmu > toler)
            {
                sprintf(ebuf, "XXX");
            }
            else if (fabs(cosa) < 0.1)
            {
                sprintf(ebuf, "YYY");
            }
            else
            {
                ebuf[0] = '\0';
            }
            fprintf(fp, "%-4s (%6.2f,%6.2f,%6.2f) |Mu| = %5.2f Dev: (%6.2f,%6.2f,%6.2f) |%5.2f|%s\n",
                    calc_name, mol->mu_calc_[XX], mol->mu_calc_[YY], mol->mu_calc_[ZZ], 
                    norm(mol->mu_calc_), dmu[XX], dmu[YY], dmu[ZZ], ndmu, ebuf);
        }
        else if (strcmp(calc_name, (char *)"ESP") == 0)
        {
            rvec_sub(mol->mu_elec_, mol->mu_esp_, dmu);
            ndmu = norm(dmu);
            cosa = cos_angle(mol->mu_elec_, mol->mu_esp_);
            if (ndmu > toler)
            {
                sprintf(ebuf, "XXX");
            }
            else if (fabs(cosa) < 0.1)
            {
                sprintf(ebuf, "YYY");
            }
            else
            {
                ebuf[0] = '\0';
            }
            fprintf(fp, "%-4s (%6.2f,%6.2f,%6.2f) |Mu| = %5.2f Dev: (%6.2f,%6.2f,%6.2f) |%5.2f|%s\n",
                    calc_name, mol->mu_esp_[XX], mol->mu_esp_[YY], mol->mu_esp_[ZZ], 
                    norm(mol->mu_esp_), dmu[XX], dmu[YY], dmu[ZZ], ndmu, ebuf);
        }
        else if (strcmp(calc_name, (char *)"MPA") == 0)
        {
            rvec_sub(mol->mu_elec_, mol->mu_mulliken_, dmu);
            ndmu = norm(dmu);
            cosa = cos_angle(mol->mu_elec_, mol->mu_mulliken_);
            if (ndmu > toler)
            {
                sprintf(ebuf, "XXX");
            }
            else if (fabs(cosa) < 0.1)
            {
                sprintf(ebuf, "YYY");
            }
            else
            {
                ebuf[0] = '\0';
            }
            fprintf(fp, "%-4s (%6.2f,%6.2f,%6.2f) |Mu| = %5.2f Dev: (%6.2f,%6.2f,%6.2f) |%5.2f|%s\n",
                    calc_name, mol->mu_mulliken_[XX], mol->mu_mulliken_[YY], mol->mu_mulliken_[ZZ], 
                    norm(mol->mu_mulliken_), dmu[XX], dmu[YY], dmu[ZZ], ndmu, ebuf);
        }
        else if (strcmp(calc_name, (char *)"HPA") == 0)
        {
            rvec_sub(mol->mu_elec_, mol->mu_hirshfeld_, dmu);
            ndmu = norm(dmu);
            cosa = cos_angle(mol->mu_elec_, mol->mu_hirshfeld_);
            if (ndmu > toler)
            {
                sprintf(ebuf, "XXX");
            }
            else if (fabs(cosa) < 0.1)
            {
                sprintf(ebuf, "YYY");
            }
            else
            {
                ebuf[0] = '\0';
            }
            fprintf(fp, "%-4s (%6.2f,%6.2f,%6.2f) |Mu| = %5.2f Dev: (%6.2f,%6.2f,%6.2f) |%5.2f|%s\n",
                    calc_name, mol->mu_hirshfeld_[XX], mol->mu_hirshfeld_[YY], mol->mu_hirshfeld_[ZZ], 
                    norm(mol->mu_hirshfeld_), dmu[XX], dmu[YY], dmu[ZZ], ndmu, ebuf);
        }
        else if (strcmp(calc_name, (char *)"CM5") == 0)
        {
            rvec_sub(mol->mu_elec_, mol->mu_cm5_, dmu);
            ndmu = norm(dmu);
            cosa = cos_angle(mol->mu_elec_, mol->mu_cm5_);
            if (ndmu > toler)
            {
                sprintf(ebuf, "XXX");
            }
            else if (fabs(cosa) < 0.1)
            {
                sprintf(ebuf, "YYY");
            }
            else
            {
                ebuf[0] = '\0';
            }
            fprintf(fp, "%-4s (%6.2f,%6.2f,%6.2f) |Mu| = %5.2f Dev: (%6.2f,%6.2f,%6.2f) |%5.2f|%s\n",
                    calc_name, mol->mu_cm5_[XX], mol->mu_cm5_[YY], mol->mu_cm5_[ZZ], 
                    norm(mol->mu_cm5_), dmu[XX], dmu[YY], dmu[ZZ], ndmu, ebuf);
        }
        else
        {
            fprintf(fp, "Dipole analysis\n");
            fprintf(fp, "Electronic   (%6.2f,%6.2f,%6.2f) |Mu| = %5.2f\n",
                    mol->mu_elec_[XX], mol->mu_elec_[YY], mol->mu_elec_[ZZ], norm(mol->mu_elec_));
        }
    }
}

void print_electric_props(FILE                           *fp, 
                          std::vector<alexandria::MyMol>  mymol,
                          const char                     *qhisto,
                          const char                     *DipCorr,
                          const char                     *MuCorr, 
                          const char                     *Qcorr,
                          const char                     *EspCorr,
                          const char                     *alphaCorr,
                          const char                     *isopolCorr,
                          const char                     *anisopolCorr, 
                          real                            dip_toler, 
                          real                            quad_toler,
                          real                            alpha_toler, 
                          const gmx_output_env_t         *oenv,
                          bool                            bPolar,
                          bool                            bDipole,
                          bool                            bQuadrupole,
                          bool                            bfullTensor,
                          IndexCount                     *indexCount,
                          real                            hfac,
                          t_commrec                      *cr,
                          real                            efield)
{
    int           i    = 0, j     = 0, n     = 0;
    int           nout = 0, mm    = 0, nn    = 0;
    real          sse  = 0, rms   = 0, sigma = 0;
    real          aver = 0, error = 0, qCalc = 0;
    
    FILE          *dipc, *muc,  *Qc;
    FILE          *hh,   *espc, *alphac, *isopolc, *anisopolc;   
        
    struct AtomTypeLsq {
        std::string atomtype;
        gmx_stats_t lsq;
    };    
    enum {
        eprCalc, eprESP, eprMPA, eprHPA, eprCM5, eprNR
    };
    
    gmx_stats_t               lsq_mu[eprNR], lsq_dip[eprNR], lsq_quad[eprNR];
    gmx_stats_t               lsq_esp, lsq_alpha, lsq_isoPol, lsq_anisoPol;
    const char               *eprnm[eprNR] = {"Calc", "ESP", "MPA", "HPA", "CM5"};
    std::vector<AtomTypeLsq>  lsqt;

    for (int i = 0; i < eprNR; i++)
    {
        lsq_quad[i] = gmx_stats_init();
        lsq_dip[i]  = gmx_stats_init();
        lsq_mu[i]   = gmx_stats_init();
    }
    lsq_esp      = gmx_stats_init();
    lsq_alpha    = gmx_stats_init();
    lsq_isoPol   = gmx_stats_init();
    lsq_anisoPol = gmx_stats_init();
    n            = 0;
    
    for (auto ai = indexCount->beginIndex(); ai < indexCount->endIndex(); ++ai)
    {
        AtomTypeLsq k;
        k.atomtype.assign(ai->name());
        k.lsq = gmx_stats_init();
        lsqt.push_back(std::move(k));
    }       
    for (auto &mol: mymol)
    {
        if (mol.eSupp_ != eSupportNo)
        {
            fprintf(fp, "Molecule %d: %s Qtot: %d, Multiplicity %d\n", n+1,
                    mol.molProp()->getMolname().c_str(),
                    mol.molProp()->getCharge(),
                    mol.molProp()->getMultiplicity());
                        
            mol.CalcDipole();
            print_dipole(fp, &mol, (char *)"Electronic",   dip_toler);
            print_dipole(fp, &mol, (char *)"Calc", dip_toler);
            print_dipole(fp, &mol, (char *)"ESP",  dip_toler);
            print_dipole(fp, &mol, (char *)"MPA",  dip_toler);
            print_dipole(fp, &mol, (char *)"HPA",  dip_toler);
            print_dipole(fp, &mol, (char *)"CM5",  dip_toler);
            
            gmx_stats_add_point(lsq_dip[eprCalc], mol.dip_elec_, mol.dip_calc_,      0, 0);
            gmx_stats_add_point(lsq_dip[eprESP],  mol.dip_elec_, mol.dip_esp_,       0, 0);
            gmx_stats_add_point(lsq_dip[eprMPA],  mol.dip_elec_, mol.dip_mulliken_,  0, 0);
            gmx_stats_add_point(lsq_dip[eprHPA],  mol.dip_elec_, mol.dip_hirshfeld_, 0, 0);
            gmx_stats_add_point(lsq_dip[eprCM5],  mol.dip_elec_, mol.dip_cm5_,       0, 0);

            sse += gmx::square(mol.dip_elec_ - mol.dip_calc_);
            
            mol.CalcQuadrupole();
            print_quadrapole(fp, &mol, (char *)"Electronic",   quad_toler);
            print_quadrapole(fp, &mol, (char *)"Calc", quad_toler);
            print_quadrapole(fp, &mol, (char *)"ESP",  quad_toler);
            print_quadrapole(fp, &mol, (char *)"MPA",  quad_toler);
            print_quadrapole(fp, &mol, (char *)"HPA",  quad_toler);
            print_quadrapole(fp, &mol, (char *)"CM5",  quad_toler);
            
            for (mm = 0; mm < DIM; mm++)
            {
                gmx_stats_add_point(lsq_mu[eprCalc], mol.mu_elec_[mm], mol.mu_calc_[mm],      0, 0);
                gmx_stats_add_point(lsq_mu[eprESP],  mol.mu_elec_[mm], mol.mu_esp_[mm],       0, 0);
                gmx_stats_add_point(lsq_mu[eprMPA],  mol.mu_elec_[mm], mol.mu_mulliken_[mm],  0, 0);
                gmx_stats_add_point(lsq_mu[eprHPA],  mol.mu_elec_[mm], mol.mu_hirshfeld_[mm], 0, 0);
                gmx_stats_add_point(lsq_mu[eprCM5],  mol.mu_elec_[mm], mol.mu_cm5_[mm],       0, 0);
                
                if (bfullTensor)
                {
                    for (nn = 0; nn < DIM; nn++)
                    {
                        gmx_stats_add_point(lsq_quad[eprCalc], mol.Q_elec_[mm][nn], mol.Q_calc_[mm][nn],      0, 0);
                        gmx_stats_add_point(lsq_quad[eprESP],  mol.Q_elec_[mm][nn], mol.Q_esp_[mm][nn],       0, 0);
                        gmx_stats_add_point(lsq_quad[eprMPA],  mol.Q_elec_[mm][nn], mol.Q_mulliken_[mm][nn],  0, 0);
                        gmx_stats_add_point(lsq_quad[eprHPA],  mol.Q_elec_[mm][nn], mol.Q_hirshfeld_[mm][nn], 0, 0);
                        gmx_stats_add_point(lsq_quad[eprCM5],  mol.Q_elec_[mm][nn], mol.Q_cm5_[mm][nn],       0, 0);
                    }
                }
                else
                {
                    gmx_stats_add_point(lsq_quad[eprCalc], mol.Q_elec_[mm][mm], mol.Q_calc_[mm][mm],      0, 0);
                    gmx_stats_add_point(lsq_quad[eprESP],  mol.Q_elec_[mm][mm], mol.Q_esp_[mm][mm],       0, 0);
                    gmx_stats_add_point(lsq_quad[eprMPA],  mol.Q_elec_[mm][mm], mol.Q_mulliken_[mm][nn],  0, 0);
                    gmx_stats_add_point(lsq_quad[eprHPA],  mol.Q_elec_[mm][mm], mol.Q_hirshfeld_[mm][nn], 0, 0);
                    gmx_stats_add_point(lsq_quad[eprCM5],  mol.Q_elec_[mm][mm], mol.Q_cm5_[mm][nn],       0, 0);

                }
            }
            
            if(bPolar)
            {
                mol.CalcPolarizability(efield, cr, nullptr);
                print_polarizability(fp, &mol, (char *)"Electronic", alpha_toler);
                print_polarizability(fp, &mol, (char *)"Calc",       alpha_toler);
                gmx_stats_add_point(lsq_isoPol, mol.isoPol_elec_, mol.isoPol_calc_,       0, 0);
                gmx_stats_add_point(lsq_anisoPol, mol.anisoPol_elec_, mol.anisoPol_calc_, 0, 0);
                for (mm = 0; mm < DIM; mm++)
                {
                    gmx_stats_add_point(lsq_alpha, mol.alpha_elec_[mm][mm], mol.alpha_calc_[mm][mm], 0, 0);
                }
            }
            
            rms = mol.espRms();
            fprintf(fp,   "ESP rms: %g (Hartree/e)\n", rms);                                  
            auto nEsp     = mol.Qgresp_.nEsp();
            auto EspPoint = mol.Qgresp_.espPoint();
            for (size_t i = 0; i < nEsp; i++)
            {
                gmx_stats_add_point(lsq_esp, gmx2convert(EspPoint[i].v(),eg2cHartree_e), gmx2convert(EspPoint[i].vCalc(), eg2cHartree_e), 0, 0);
            }
            
            fprintf(fp, "Atom   Type      q_Calc     q_ESP     q_MPA     q_HPA     q_CM5       x       y       z\n");
            for (j = i = 0; j < mol.topology_->atoms.nr; j++)
            {
                if (mol.topology_->atoms.atom[j].ptype == eptAtom)
                {               
                    const char *at = *(mol.topology_->atoms.atomtype[j]);
                    if(indexCount->isOptimized(at))
                    {
                        auto        k  = std::find_if(lsqt.begin(), lsqt.end(),
                                                      [at](const AtomTypeLsq &atlsq)
                                                      {
                                                          return atlsq.atomtype.compare(at) == 0;
                                                      });                                                 
                        if (k != lsqt.end())
                        {
                            qCalc = mol.topology_->atoms.atom[j].q;
                            if(nullptr != mol.shellfc_)
                            {
                                qCalc += mol.topology_->atoms.atom[j+1].q;
                            }
                            gmx_stats_add_point(k->lsq, mol.qESP_[i], qCalc, 0, 0);
                        }                        
                        fprintf(fp, "%-2d%3d  %-5s  %8.4f  %8.4f  %8.4f  %8.4f  %8.4f%8.3f%8.3f%8.3f\n",
                                mol.topology_->atoms.atom[j].atomnumber,
                                j+1,
                                *(mol.topology_->atoms.atomtype[j]),
                                qCalc, 
                                mol.qESP_[i],
                                mol.qMulliken_[i],
                                mol.qHirshfeld_[i],
                                mol.qCM5_[i],
                                mol.state_->x[j][XX], 
                                mol.state_->x[j][YY], 
                                mol.state_->x[j][ZZ]); 
                    }
                    i++;
                }
            }
            fprintf(fp, "\n");
            n++;
        }
    }

    fprintf(fp, "Dipoles are %s in Calc Parametrization.\n",     (bDipole ?     "used" : "not used"));
    fprintf(fp, "Quadrupoles are %s in Calc Parametrization.\n", (bQuadrupole ? "used" : "not used"));
    fprintf(fp, "\n");

    print_stats(fp, (char *)"ESP           (Hartree/e)",  lsq_esp,           true,  (char *)"Electronic", (char *)"Calc");
    print_stats(fp, (char *)"Dipoles       (Debye)",      lsq_mu[eprCalc],   false, (char *)"Electronic", (char *)"Calc");
    print_stats(fp, (char *)"Dipole Moment (Debye)",      lsq_dip[eprCalc],  false, (char *)"Electronic", (char *)"Calc");
    print_stats(fp, (char *)"Quadrupoles   (Buckingham)", lsq_quad[eprCalc], false, (char *)"Electronic", (char *)"Calc");
    if (bPolar)
    {
        print_stats(fp, (char *)"Principal Components of Polarizability (A^3)",  lsq_alpha, false,  (char *)"Electronic", (char *)"Calc");
        print_stats(fp, (char *)"Isotropic Polarizability (A^3)",    lsq_isoPol,   false,  (char *)"Electronic", (char *)"Calc");
        print_stats(fp, (char *)"Anisotropic Polarizability (A^3)",  lsq_anisoPol, false,  (char *)"Electronic", (char *)"Calc");
    }
    fprintf(fp, "\n");

    print_stats(fp, (char *)"Dipoles",       lsq_mu[eprESP],   true,  (char *)"Electronic", (char *)"ESP");
    print_stats(fp, (char *)"Dipole Moment", lsq_dip[eprESP],  false, (char *)"Electronic", (char *)"ESP");
    print_stats(fp, (char *)"Quadrupoles",   lsq_quad[eprESP], false, (char *)"Electronic", (char *)"ESP");       
    fprintf(fp, "\n");

    print_stats(fp, (char *)"Dipoles",       lsq_mu[eprMPA],   true,  (char *)"Electronic", (char *)"MPA");
    print_stats(fp, (char *)"Dipole Moment", lsq_dip[eprMPA],  false, (char *)"Electronic", (char *)"MPA");
    print_stats(fp, (char *)"Quadrupoles",   lsq_quad[eprMPA], false, (char *)"Electronic", (char *)"MPA");   
    fprintf(fp, "\n");

    print_stats(fp, (char *)"Dipoles",       lsq_mu[eprHPA],   true,  (char *)"Electronic", (char *)"HPA");
    print_stats(fp, (char *)"Dipole Moment", lsq_dip[eprHPA],  false, (char *)"Electronic", (char *)"HPA");
    print_stats(fp, (char *)"Quadrupoles",   lsq_quad[eprHPA], false, (char *)"Electronic", (char *)"HPA");   
    fprintf(fp, "\n");

    print_stats(fp, (char *)"Dipoles",       lsq_mu[eprCM5],   true,  (char *)"Electronic", (char *)"CM5");
    print_stats(fp, (char *)"Dipole Moment", lsq_dip[eprCM5],  false, (char *)"Electronic", (char *)"CM5");
    print_stats(fp, (char *)"Quadrupoles",   lsq_quad[eprCM5], false, (char *)"Electronic", (char *)"CM5");   
    fprintf(fp, "\n");
    
    std::vector<const char*> atypes;
    for (const auto &k : lsqt)
    {
        atypes.push_back(k.atomtype.c_str());
    }
    
    hh = xvgropen(qhisto, "Histogram for charges", "q (e)", "a.u.", oenv);
    xvgr_legend(hh, atypes.size(), atypes.data(), oenv);
    
    fprintf(fp, "\nParameters are optimized for %zu atom types:\n", atypes.size());    
    for (auto k = lsqt.begin(); k < lsqt.end(); ++k)
    {
        int   nbins;
        if (gmx_stats_get_npoints(k->lsq, &nbins) == estatsOK)
        {
            real *x, *y;
            fprintf(fp, "%-4d copies for %4s\n", nbins, k->atomtype.c_str());
            if (gmx_stats_make_histogram(k->lsq, 0, &nbins, ehistoY, 1, &x, &y) == estatsOK)
            {
                fprintf(hh, "@type xy\n");
                for (int i = 0; i < nbins; i++)
                {
                    fprintf(hh, "%10g  %10g\n", x[i], y[i]);
                }
                fprintf(hh, "&\n");
                
                free(x);
                free(y);
            }
        }
        gmx_stats_free(k->lsq);
    }
    fclose(hh);
    fprintf(fp, "\n");
    
    dipc = xvgropen(DipCorr, "Dipole Moment (Debye)", "Electronic", "Empirical", oenv);
    xvgr_symbolize(dipc, 5, eprnm, oenv);
    print_lsq_set(dipc, lsq_dip[eprCalc]);
    print_lsq_set(dipc, lsq_dip[eprESP]);
    print_lsq_set(dipc, lsq_dip[eprMPA]);
    print_lsq_set(dipc, lsq_dip[eprHPA]);
    print_lsq_set(dipc, lsq_dip[eprCM5]);
    fclose(dipc);
    
    muc = xvgropen(MuCorr, "Dipoles (Debye)", "Electronic", "Empirical", oenv);
    xvgr_symbolize(muc, 5, eprnm, oenv);
    print_lsq_set(muc, lsq_mu[eprCalc]);
    print_lsq_set(muc, lsq_mu[eprESP]);
    print_lsq_set(muc, lsq_mu[eprMPA]);
    print_lsq_set(muc, lsq_mu[eprHPA]);
    print_lsq_set(muc, lsq_mu[eprCM5]);
    fclose(muc);

    Qc = xvgropen(Qcorr, "Quadrupoles (Buckingham)", "Electronic", "Empirical", oenv);
    xvgr_symbolize(Qc, 5, eprnm, oenv);
    print_lsq_set(Qc, lsq_quad[eprCalc]);
    print_lsq_set(Qc, lsq_quad[eprESP]);
    print_lsq_set(Qc, lsq_quad[eprMPA]);
    print_lsq_set(Qc, lsq_quad[eprHPA]);
    print_lsq_set(Qc, lsq_quad[eprCM5]);
    fclose(Qc);
    
    espc = xvgropen(EspCorr, "Electrostatic Potential (Hartree/e)", "Electronic", "Calc", oenv);
    xvgr_symbolize(espc, 1, eprnm, oenv);
    print_lsq_set(espc, lsq_esp);
    fclose(espc);
    
    if (bPolar)
    {
        alphac = xvgropen(alphaCorr, "Pricipal Components of Polarizability Tensor (A\\S3\\N)", "Electronic", "Calc", oenv);
        xvgr_symbolize(alphac, 1, eprnm, oenv);
        print_lsq_set(alphac, lsq_alpha);
        fclose(alphac);
        
        isopolc = xvgropen(isopolCorr, "Isotropic Polarizability (A\\S3\\N)", "Electronic", "Calc", oenv);
        xvgr_symbolize(isopolc, 1, eprnm, oenv);
        print_lsq_set(isopolc, lsq_isoPol);
        fclose(isopolc);
        
        anisopolc = xvgropen(anisopolCorr, "Anisotropic Polarizability (A\\S3\\N)", "Electronic", "Calc", oenv);
        xvgr_symbolize(anisopolc, 1, eprnm, oenv);
        print_lsq_set(anisopolc, lsq_anisoPol);
        fclose(anisopolc);
    }

    fprintf(fp, "hfac = %g\n", hfac);
    gmx_stats_get_ase(lsq_mu[eprCalc], &aver, &sigma, &error);
    sigma = sqrt(sse/n);
    nout  = 0;
    fprintf(fp, "Overview of dipole moment outliers (> %.3f off)\n", 2*sigma);
    fprintf(fp, "----------------------------------\n");
    fprintf(fp, "%-20s  %12s  %12s  %12s\n", "Name", "Calc", "Electronic", "Deviation (Debye)");
            
    for (auto &mol : mymol)
    {
        auto deviation = std::abs(mol.dip_calc_ - mol.dip_elec_);
        if ((mol.eSupp_ != eSupportNo) &&
            (mol.dip_elec_ > sigma) &&
            (deviation > 2*sigma))
        {
            fprintf(fp, "%-20s  %12.3f  %12.3f  %12.3f\n",
                    mol.molProp()->getMolname().c_str(),
                    mol.dip_calc_, mol.dip_elec_, deviation);
            nout++;
        }
    }
    if (nout)
    {
        printf("There were %d outliers. See at the very bottom of the log file\n", nout);
    }
    else
    {
        printf("No outliers! Well done.\n");
    }
    for (int i = 0; i < eprNR; i++)
    {
        gmx_stats_free(lsq_quad[i]);
        gmx_stats_free(lsq_mu[i]);
        gmx_stats_free(lsq_dip[i]);    
    }
    gmx_stats_free(lsq_esp);
    gmx_stats_free(lsq_alpha);
    gmx_stats_free(lsq_isoPol);
    gmx_stats_free(lsq_anisoPol);
}
}


