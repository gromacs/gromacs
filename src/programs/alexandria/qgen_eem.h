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
 
 
#ifndef QGEN_EEM_H
#define QGEN_EEM_H

#include <cstdio>

#include "gromacs/math/vectypes.h"
#include "gromacs/mdtypes/state.h"

#include "poldata.h"

struct t_atoms;

enum {
    eQGEN_OK, eQGEN_NOTCONVERGED, eQGEN_NOSUPPORT, eQGEN_ERROR, eQGEN_NR
};

enum ChargeGenerationAlgorithm {
    eqgNONE, eqgEEM, eqgESP
};

namespace alexandria
{

class QgenEem
{
    public:
    
        QgenEem(){};
    
        void setInfo(const Poldata            &pd,
                     t_atoms                 *atoms,
                     ChargeDistributionModel  iChargeDistributionModel,
                     double                   hfac, 
                     int                      qtotal, 
                     bool                     haveShell);
                     
        void updateInfo(const Poldata &pd);

        int generateChargesSm(FILE            *fp,
                              const Poldata    &pd,
                              t_atoms         *atoms,  
                              double          *chieq,
                              PaddedRVecVector x);

        int generateCharges(FILE              *fp,
                            const std::string  molname,
                            const Poldata     &pd,
                            t_atoms           *atoms,
                            PaddedRVecVector   x);     
                            
        double rms() { return rms_; }

        const char *message() const;
        
        gmx_bool SplitQ(ChargeDistributionModel iDistributionModel);

        int getNzeta(int atom);

        int getRow(int atom, int z);

        std::vector<std::vector<double> > q() { return q_;}
        
        int natom() {return natom_;}
        
        double getQ(int atom, int z);

        void checkSupport(const Poldata &pd);

        double getZeta(int atom, int z);

        void print(FILE *fp, t_atoms *atoms);

        void debugFun(FILE *fp);

    private:
        gmx_bool                                           bWarned_;
        ChargeDistributionModel                            iChargeDistributionModel_;
        int                                                natom_, eQGEN_;
        double                                             qtotal_, chieq_, hfac_;
        double                                             Jcs_, Jss_, rms_, hardnessFactor_;
        /* For each atom i there is an elem, atomnr, chi0, rhs, j00 and x */
        std::vector<std::string>                           elem_;
        std::vector<int>                                   atomnr_;
        std::vector<double>                                chi0_, rhs_, j00_;
        std::vector<gmx::RVec>                             x_;
        /* Jab is a matrix over atom pairs */
        std::vector<std::vector<double>>                   Jcc_;
        /* For each atom i there are nZeta[i] row, q and zeta entries */
        std::vector<int>                                   nZeta_, coreIndex_, shellIndex_;
        std::vector<std::vector<int> >                     row_;
        bool                                               bAllocSave_, bHaveShell_;
        std::vector<std::vector<double>>                   q_, zeta_, qsave_, zetasave_;


        double calcJ(ChargeDistributionModel iChargeDistributionModel,
                     rvec                    xI, 
                     rvec                    xJ,
                     double                  zetaI,
                     double                  zetaJ,
                     int                     rowI,
                     int                     rowJ);

        void copyChargesToAtoms(t_atoms *atoms);
        
        void calcJcc(t_atoms *atoms);
        
        void calcJcs(t_atoms *atoms,
                     int      top_ndx,
                     int      eem_ndx);
        
        void calcJss(t_atoms *atoms,
                     int      top_ndx,
                     int      eem_ndx);

        void solveQEem(FILE *fp);
        
        void updatePositions(PaddedRVecVector x, t_atoms *atoms);

        double calcSij(int i, int j);

        int generateChargesBultinck(FILE              *fp,
                                    const Poldata     &pd,
                                    t_atoms           *atoms,
                                    PaddedRVecVector   x);
        void calcRhs(t_atoms *atoms);
};
}
#endif
