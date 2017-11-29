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
 * \author  Mohammad Mehdi Ghahremanpour <mohammad.ghahremanpour@icm.uu.se>
 * \author David van der Spoel <david.vanderspoel@icm.uu.se>
 */
#ifndef MOLDIP_H
#define MOLDIP_H

#include "gromacs/mdrunutility/mdmodules.h"
#include "gromacs/utility/real.h"

#include "mymol.h"

typedef struct {
    int       n;
    int       nopt;
    int       nconst;
    int       nopt_c;
    int      *tot_count;
    int      *count;
    char    **name;
    gmx_bool *bConst;
} t_index_count;

extern char *opt_index_count(t_index_count *ic);

enum {
    ermsBOUNDS = 0, 
    ermsMU     = 1, 
    ermsQUAD   = 2, 
    ermsCHARGE = 3, 
    ermsESP    = 4,
    ermsEPOT   = 5, 
    ermsForce2 = 6, 
    ermsPolar  = 7,
    ermsTOT    = 8, 
    ermsNR     = 9
};

namespace alexandria
{

class AtomIndex
{
    private:
        std::string name_;
        int         count_;
        bool        const_;
    public:
        AtomIndex(const std::string name, bool bConst) : name_(name), count_(1), const_(bConst) {};

        const std::string name() const { return name_; }

        int count() const { return count_; }

        bool isConst() const { return const_; }

        void increment() { count_++; }

        void decrement()
        {
            if (count_ > 0)
            {
                count_--;
            }
            else
            {
                fprintf(stderr, "Trying to decrease number of atoms %s below zero\n",
                        name().c_str());
            }
        }
};

using AtomIndexIterator      = typename std::vector<AtomIndex>::iterator;
using AtomIndexConstIterator = typename std::vector<AtomIndex>::const_iterator;

class IndexCount
{
    private:
        std::vector<AtomIndex> atomIndex_;
        std::vector<int>       totCount_;
    public:
        IndexCount() {};

        int nOpt() const
        {
            return std::count_if(atomIndex_.begin(), atomIndex_.end(),
                                 [](AtomIndex ai) { return ai.isConst(); });
        }
        int nConst() const { return atomIndex_.size() - nOpt(); }

        int nName(const std::string &name) const;

        void sumCount(t_commrec *cr);

        int cleanIndex(int   minimum_data,
                       FILE *fp);

        void addName(const std::string &name,
                     bool               bConst);

        void incrementName(const std::string &name);

        void decrementName(const std::string &name);

        int count(const std::string &name);
        
        bool isOptimized(const std::string &name);

        std::vector<AtomIndex>::iterator beginIndex() { return atomIndex_.begin(); }

        std::vector<AtomIndex>::iterator endIndex() { return atomIndex_.end(); }

        std::vector<AtomIndex>::const_iterator beginIndex() const { return atomIndex_.begin(); }

        std::vector<AtomIndex>::const_iterator endIndex() const { return atomIndex_.end(); }

        std::vector<AtomIndex>::iterator findName(const std::string &name)
        {
            return std::find_if(atomIndex_.begin(), atomIndex_.end(),
                                [name](const AtomIndex a)
                                { return a.name().compare(name) == 0; });
        }
};

using IndexCountIterator      = typename std::vector<IndexCount>::iterator;
using IndexCountConstIterator = typename std::vector<IndexCount>::const_iterator;

class MolDip
{
    private:
    public:
        int                             nmol_support_;
        int                             mindata_;
        real                            J0_min_;
        real                            Chi0_min_;
        real                            zeta_min_;
        real                            J0_max_;
        real                            Chi0_max_;
        real                            zeta_max_;
        real                            hfac_;
        real                            hfac0_;
        real                            decrzeta_;
        real                            ener_[ermsNR];
        real                            fc_[ermsNR];
        char                           *fixchi_;
        gmx_bool                        bOptHfac_;
        gmx_bool                        bPol_;
        gmx_bool                        bQM_;
        gmx_bool                        bDone_;
        gmx_bool                        bFinal_;
        gmx_bool                        bGaussianBug_;
        gmx_bool                        bFitZeta_;
        gmx_bool                        bfullTensor_;
        gmx_bool                        bGenViste_;
        Poldata                         pd_;
        t_commrec                      *cr_;
        t_inputrec                     *inputrec_;
        IndexCount                      indexCount_;
        gmx_hw_info_t                  *hwinfo_;
        gmx_atomprop_t                  atomprop_;
        ChargeDistributionModel         iChargeDistributionModel_;
        ChargeGenerationAlgorithm       iChargeGenerationAlgorithm_;
        gmx::MDModules                  mdModules_;
        std::vector<alexandria::MyMol>  mymol_;

        MolDip();

        ~MolDip() {};

        IndexCount *indexCount() { return &indexCount_; }
        
        immStatus check_data_sufficiency(alexandria::MyMol  mymol, 
                                         IndexCount        *ic);
                                                 
        void Init(t_commrec                *cr,
                  gmx_bool                  bQM,
                  gmx_bool                  bGaussianBug,
                  ChargeDistributionModel   iChargeDistributionModel,
                  ChargeGenerationAlgorithm iChargeGenerationAlgorithm,
                  real                      rDecrZeta,
                  real                      J0_min,
                  real                      Chi0_min,
                  real                      zeta_min,
                  real                      J0_max,
                  real                      Chi0_max,
                  real                      zeta_max,
                  real                      fc_bound,
                  real                      fc_mu,
                  real                      fc_quad,
                  real                      fc_charge,
                  real                      fc_esp,
                  real                      fc_epot,
                  real                      fc_force,
                  real                      fc_polar,
                  char                     *fixchi,
                  gmx_bool                  bOptHfac,
                  real                      hfac,
                  gmx_bool                  bPol,
                  gmx_bool                  bFitZeta, 
                  gmx_hw_info_t            *hwinfo,
                  gmx_bool                  bfullTensor,
                  int                       mindata,
                  gmx_bool                  bGenViste);

        void Read(FILE                      *fp,
                  const char                *fn,
                  const char                *pd_fn,
                  gmx_bool                   bZero,
                  char                      *opt_elem,
                  char                      *const_elem,
                  char                      *lot,
                  const MolSelect           &gms,
                  real                       watoms,
                  gmx_bool                   bCheckSupport,
                  bool                       bPairs,
                  bool                       bDihedral,
                  bool                       bPolar,
                  bool                       bZPE,
                  const char                *tabfn,
                  int                        qcycle,
                  real                       qtol,
                  bool                       qsymm);

};

}

#endif
