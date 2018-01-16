/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2016,2017, by the GROMACS development team, led by
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
#ifndef MOLGEN_H
#define MOLGEN_H

#include "gromacs/commandline/pargs.h"
#include "gromacs/gmxlib/network.h"
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

enum eRMS {
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

//! \brief Return string corresponding to eRMS
const char *rmsName(int e);

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

class MolGen
{
    private:
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
        real                            watoms_;
        real                            qtol_;
        int                             qcycle_;
        real                            ener_[ermsNR];
        real                            fc_[ermsNR];
        char                           *fixchi_;
        gmx_bool                        bOptHfac_;
        gmx_bool                        bQM_;
        gmx_bool                        bDone_;
        gmx_bool                        bFinal_;
        gmx_bool                        bGenVsite_;
        gmx_bool                        qsymm_;
        Poldata                         pd_;
        t_commrec                      *cr_;
        gmx::MDLogger                   mdlog_;
        t_inputrec                     *inputrec_;
        IndexCount                      indexCount_;
        gmx_hw_info_t                  *hwinfo_;
        gmx_atomprop_t                  atomprop_;
        ChargeDistributionModel         iChargeDistributionModel_;
        ChargeGenerationAlgorithm       iChargeGenerationAlgorithm_;
        gmx::MDModules                  mdModules_;
        std::vector<alexandria::MyMol>  mymol_;
        const char                     *lot_;
    public:

        MolGen();

        ~MolGen();

        IndexCount *indexCount() { return &indexCount_; }

        //! \brief Add options to the command line
        void addOptions(std::vector<t_pargs> *pargs);

        //! \brief Process options after parsing
        void optionsFinished();

        immStatus check_data_sufficiency(alexandria::MyMol  mymol,
                                         IndexCount        *ic);
        //! \brief Return the poldata as const variable
        const Poldata &poldata() const { return pd_; }
        //! \brief Return the poldata
        Poldata &poldata() { return pd_; }

        //! \brief Return the atomprop structure
        gmx_atomprop_t atomprop() const { return atomprop_; }

        //! \brief Return the const vector of molecules
        const std::vector<MyMol> &mymols() const { return mymol_; }

        //! \brief Return the mutable vector of molecules
        std::vector<MyMol> &mymols() { return mymol_; }

        //! \brief Return the ChargeDistributionModel
        ChargeDistributionModel iChargeDistributionModel() const { return iChargeDistributionModel_; }

        //! \brief Return the ChargeGenerationAlgorithm
        ChargeGenerationAlgorithm iChargeGenerationAlgorithm() const { return iChargeGenerationAlgorithm_; }

        //! \brief Return whether to optimize the H factor
        bool optHfac() const { return bOptHfac_; }

        //! \brief Return the H factor
        real hfac() const { return hfac_; }

        //! \brief Set the H factor
        void setHfac(real hfac) { hfac_ = hfac; }

        real hfacDiff()
        {
            if (hfac_ > hfac0_)
            {
                return hfac_ - hfac0_;
            }
            else if (hfac_ < -hfac0_)
            {
                return (hfac() + hfac0_);
            }
            return 0;
        }
        //! \brief The atom to fix
        const char *fixchi() const { return fixchi_; }

        //! \brief Return ESP weighting factor for atoms
        real watoms() const { return watoms_; }

        //! \brief Return minimum amount of data needed
        int mindata() const { return mindata_; }

        //! \brief Return const communication record
        const t_commrec *commrec() const { return cr_; }

        //! \brief Return non-const communication record
        t_commrec *commrec() { return cr_; }

        //! \brief Is this the last calculation?
        bool final() const { return bFinal_; }

        //! \brief Are we using QM only?
        bool bQM() const { return bQM_; }

        //! \brief Return level of theory
        const char *lot() const { return lot_; }

        void setFinal() { bFinal_ = true; }

        double zetaMin() const { return zeta_min_; }
        double zetaMax() const { return zeta_max_; }
        double J0Min() const { return J0_min_; }
        double J0Max() const { return J0_max_; }
        double chi0Min() const { return Chi0_min_; }
        double chi0Max() const { return Chi0_max_; }
        
        void setEnergy(int qt, real ener)
        {
            ener_[qt] = ener;
        }
        
        double energy(int qt) const 
        { 
            return fc_[qt]*ener_[qt]; 
        }
        
        void resetEnergies()
        {
            for (int j = 0; j < ermsNR; j++)
            {
                ener_[j] = 0;
            }
        }

        void increaseEnergy(int qt, real delta)
        {
            ener_[qt] += delta;
        }

        
        void sumEnergies()
        {
            if (PAR(commrec()) && !final())
            {
                gmx_sum(ermsNR, ener_, commrec());
            }
        }

        void normalizeEnergies()
        {
            if (MASTER(commrec()))
            {
                for (int j = 0; j < ermsTOT; j++)
                {
                    ener_[ermsTOT] += ((fc_[j]*ener_[j])/nmol_support_);
                }
            }
        }

        void printEnergies(FILE *fp)
        {
            if (nullptr != fp && MASTER(commrec()))
            {
                fprintf(fp, "ENER:");
                for (int j = 0; j < ermsNR; j++)
                {
                    fprintf(fp, "  %8.3f", ener_[j]);
                }
                fprintf(fp, "\n");
            }
        }

        void Read(FILE                      *fp,
                  const char                *fn,
                  const char                *pd_fn,
                  gmx_bool                   bZero,
                  char                      *opt_elem,
                  char                      *const_elem,
                  const MolSelect           &gms,
                  gmx_bool                   bCheckSupport,
                  bool                       bPairs,
                  bool                       bDihedral,
                  bool                       bZPE,
                  bool                       bFitZeta,
                  const char                *tabfn);

};

}

#endif
