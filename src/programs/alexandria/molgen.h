/*
 * This source file is part of the Alexandria Chemistry Toolkit.
 *
 * Copyright (C) 2014-2019 
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
    ermsBOUNDS  = 0,
    ermsMU      = 1,
    ermsQUAD    = 2,
    ermsCHARGE  = 3,
    ermsESP     = 4,
    ermsEPOT    = 5,
    ermsForce2  = 6,
    ermsPolar   = 7,
    ermsPENALTY = 8,
    ermsTOT     = 9,
    ermsNR      = 10
};

//! \brief Return string corresponding to eRMS
const char *rmsName(int e);

namespace alexandria
{

enum eTune {
    etuneEEM = 0,
    etuneZETA = 1,
    etuneFC = 2,
    etuneNR = 3
};

class AtomIndex
{
    private:
        std::string      name_;
        int              count_;
        bool             const_;
        EempropsIterator eem_;
    public:
        AtomIndex(const std::string name, bool bConst) :
            name_(name), count_(0), const_(bConst) {}

        const std::string name() const { return name_; }

        int count() const { return count_; }

        bool isConst() const { return const_; }
        
        void setEemProps(EempropsIterator eem) { eem_ = eem; }
    
        EempropsIterator eemProps() const { return eem_; }
        
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
        int                             maxESP_;
        int                             nexcl_;
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
        real                            ener_[ermsNR] = { 0 };
        real                            fc_[ermsNR] = { 0 };
        char                           *fixchi_;
        gmx_bool                        bOptHfac_;
        gmx_bool                        bQM_;
        gmx_bool                        bDone_;
        gmx_bool                        bFinal_;
        gmx_bool                        bGenVsite_;
        gmx_bool                        qsymm_;
        gmx_bool                        constrain_; 
        Poldata                         pd_;
        t_commrec                      *cr_;
        gmx::MDLogger                   mdlog_;
        t_inputrec                     *inputrec_;
        IndexCount                      indexCount_;
        gmx_hw_info_t                  *hwinfo_;
        gmx_atomprop_t                  atomprop_;
        gmx::MDModules                  mdModules_;
        std::vector<alexandria::MyMol>  mymol_;
        const char                     *lot_;
    public:

        /*! \brief 
         * Constructor of MolGen class.
         */ 
        MolGen();

        /*! \brief 
         * Deconstructor of MolGen class.
         */
        ~MolGen();

        /*! \brief 
         * Return the indexCount.
         */
        IndexCount *indexCount() { return &indexCount_; }

        //! \brief Add options to the command line
        void addOptions(std::vector<t_pargs> *pargs, eTune etune);

        //! \brief Process options after parsing
        void optionsFinished();

        //! \brief Check that we have enough data for all parameters tp optimize. 
        immStatus check_data_sufficiency(alexandria::MyMol  mymol,
                                         IndexCount        *ic);
                                         
        //! \brief Return the poldata as const variable
        const Poldata *poldata() const { return &pd_; }
        
        //! \brief Return the poldata
        Poldata *poldata() { return &pd_; }

        //! \brief Return the atomprop structure
        gmx_atomprop_t atomprop() const { return atomprop_; }

        //! \brief Return the const vector of molecules
        const std::vector<MyMol> &mymols() const { return mymol_; }

        //! \brief Return the mutable vector of molecules
        std::vector<MyMol> &mymols() { return mymol_; }

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
        
        gmx::MDLogger  mdlog()  const {return mdlog_; }
        
        int qcycle() const { return qcycle_;}
        
        int maxPot() const { return maxESP_;}
        
        double qtol() const { return qtol_;}
        
        gmx_bool bQsym() const { return qsymm_;}
        
        gmx_bool bConstrain() const { return constrain_;}
        
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
        
        gmx_hw_info_t *hwinfo() {return hwinfo_;}
        
        //! \brief Is this the last calculation?
        bool final() const { return bFinal_; }

        //! \brief Are we using QM only?
        bool bQM() const { return bQM_; }

        //! \brief Return level of theory
        const char *lot() const { return lot_; }

        void setFinal() { bFinal_ = true; }
        
        //! \brief Return lower boundary for Gaussian and Slater exponent.
        double zetaMin() const { return zeta_min_; }
        
        //! \brief Return upper boundary for Gaussian and Slater exponent.
        double zetaMax() const { return zeta_max_; }
        
        //! \brief Return lower boundary for atomic hardness.
        double J0Min() const { return J0_min_; }
        
        //! \brief Return upper boundary for atomic hardness.
        double J0Max() const { return J0_max_; }
        
        //! \brief Return lower boundary for atomic electronegativity.
        double chi0Min() const { return Chi0_min_; }
        
        //! \brief Return upper boundary for atomic electronegativity.
        double chi0Max() const { return Chi0_max_; }
        
        //! \brief Return the number of compounds in the data set
        int nMolSupport() const { return nmol_support_; }

        //! \brief Set the energy value for the corresponding rms.
        void setEnergy(int rms, real ener)
        {
            ener_[rms] = ener;
        }
        
        //! \brief Return the energy of the corresponding rms. 
        double energy(int rms) const 
        { 
            return ener_[rms]; 
        }
        
        //! \brief Return the weighting factor of the energy.
        double weight(int rms) const 
        { 
            return fc_[rms]; 
        }
        
        //! \brief Set all the energies to zero.
        void resetEnergies()
        {
            for (int rms = 0; rms < ermsNR; rms++)
            {
                ener_[rms] = 0;
            }
        }

        //! \brief Increase the rms component of the energy vector by delta. 
        void increaseEnergy(int rms, real delta)
        {
            ener_[rms] += delta;
        }

        //! \brief Sum over the energies of the cores. 
        void sumEnergies()
        {
            // Now sum over processors
            if (PAR(commrec()) && !final())
            {
                gmx_sum(ermsNR, ener_, commrec());
            }
            ener_[ermsTOT] = 0;
            for (auto e = 0; e < ermsTOT; e++)
            {
                ener_[ermsTOT] += ener_[e];
            }
        }

        //! \brief Normalize the energy by the number of molecules.
        void normalizeEnergies()
        {
            if (MASTER(commrec()))
            {
                double normFactor = nmol_support_; // total number of molecules used in fitting
                if (weight(ermsMU) && bQM())
                {
                    normFactor *= DIM; // three dipole residuals added per molecule if compared to QM dipole vector
                }
                if (weight(ermsQUAD))
                {
                    normFactor *= DIM; // three quadrupole residuals added per molecule
                }
                for (int rms = 0; rms < ermsTOT; rms++)
                {
                    ener_[ermsTOT] += (ener_[rms]/normFactor);
                }
            }           
        }

        //! \brief Print the energy components.
        void printEnergies(FILE *fp)
        {
            if (nullptr != fp && MASTER(commrec()))
            {
                fprintf(fp, "Components of fitting function\n");
                for (int j = 0; j < ermsNR; j++)
                {
                    auto eee = energy(j);
                    if (eee > 0)
                    {
                        fprintf(fp, "%-8s  %10.3f  weight: %g\n",
                                rmsName(j), eee, fc_[j]);
                    }
                }
            }
        }

        //! \brief Read the molecular property data file to generate molecues.
        void Read(FILE            *fp,
                  const char      *fn,
                  const char      *pd_fn,
                  gmx_bool         bZero,
                  char            *opt_elem,
                  const MolSelect &gms,
                  gmx_bool         bCheckSupport,
                  bool             bPairs,
                  bool             bDihedral,
                  bool             bZPE,
                  bool             bFitZeta,
                  bool             bDHform,
                  const char      *tabfn);

};

}

#endif
