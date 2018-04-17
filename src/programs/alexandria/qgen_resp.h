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
 
 
#ifndef GMX_QGEN_RESP_H
#define GMX_QGEN_RESP_H

#include <cstdio>
#include <vector>

#include "gromacs/math/vectypes.h"
#include "gromacs/mdtypes/state.h"
#include "gromacs/statistics/statistics.h"
#include "gromacs/topology/atomprop.h"

#include "molprop.h"
#include "poldata.h"
#include "qgen_resp_atom.h"

struct gmx_output_env_t;
struct t_atoms;
struct t_symtab;

namespace alexandria
{
class EspPoint
{
    public:
        EspPoint(gmx::RVec esp, double v) : esp_(esp), v_(v)
        {
            vCalc_  = 0;
            rho_    = 0;
        }
        const gmx::RVec &esp() const { return esp_; }

        double v() const { return v_; }

        void setV(double v) { v_ = v; }

        double vCalc() const { return vCalc_; }
        
        void setVCalc(double vcalc) { vCalc_ = vcalc; }
        
        double rho() const { return rho_; }

        void setRho(double rho) { rho_ = rho; }

    private:
        //! The coordinates of a point
        gmx::RVec esp_;
        //! The measured potential
        double    v_;
        //! The calculated potential
        double    vCalc_;
        //! The electron density in the point
        double    rho_;
};

class QgenResp
{
    public:
        QgenResp();

        ChargeDistributionModel chargeDistributionModel()
        { return iDistributionModel_; }

        /*! \brief Set options for ESP charge generation
         *
         * \param[in] c          Charge distribution model eqdAXp, eqdAXg or eqdAXs
         * \param[in] watoms     Weighting factor for atoms in ESP fit
         */
        void setChargeDistributionModel(ChargeDistributionModel qd)
        { iDistributionModel_ = qd; }

        void setAtomWeight(real watoms) { watoms_ = watoms; }

        real getMolecularCharge() const { return qtot_; }

        void setMolecularCharge(int qtot) { qtot_ = qtot; }

        size_t nRespAtom() const { return ra_.size(); }

        size_t nRespAtomType() const { return ratype_.size(); }

        size_t nEsp() const { return ep_.size(); }
        
        std::vector<EspPoint> &espPoint() {return ep_;}

        void summary(FILE *gp);

        RespAtomTypeIterator beginRAT() { return ratype_.begin(); }

        RespAtomTypeIterator endRAT() { return ratype_.end(); }

        RespAtomTypeConstIterator beginRAT() const { return ratype_.begin(); }

        RespAtomTypeConstIterator endRAT() const { return ratype_.end(); }

        RespAtomTypeIterator findRAT(int atype)
        {
            return std::find_if(ratype_.begin(), ratype_.end(),
                                [atype](RespAtomType const &rat)
                                { return rat.getAtype() == atype; });
        }

        RespAtomTypeConstIterator findRAT(int atype) const
        {
            return std::find_if(ratype_.begin(), ratype_.end(),
                                [atype](RespAtomType const &rat)
                                { return rat.getAtype() == atype; });
        }

        void setAtomInfo(t_atoms                *atoms,
                         const Poldata          &pd,
                         const PaddedRVecVector  x,
                         const int               qtotal);

        void updateAtomCoords(const PaddedRVecVector x);
        
        void updateAtomCharges(t_atoms  *atoms);

        const std::string &getStoichiometry() const { return stoichiometry_; }

        void setAtomSymmetry(const std::vector<int> &symmetricAtoms);

        void addEspPoint(double x,
                         double y,
                         double z,
                         double V);

        void makeGrid(real   spacing,
                      matrix box,
                      rvec   x[]);

        void copyGrid(QgenResp &src);

        void calcRms();

        real getRms(real *wtot, real *rrms);

        void plotLsq(const gmx_output_env_t *oenv);
        
        void calcRho();

        void calcPot();
        
        void calcVShell();

        void readCube(const std::string &fn,
                      bool               bESPonly);

        void writeCube(const std::string &fn,
                       const std::string &title);

        void writeRho(const std::string &fn,
                      const std::string &title);

        void writeDiffCube(QgenResp                   &src,
                           const std::string          &cubeFn,
                           const std::string          &histFn,
                           const std::string          &title,
                           const gmx_output_env_t     *oenv,
                           int                         rho);

        void writeHisto(const std::string      &fn,
                        const std::string      &title,
                        const gmx_output_env_t *oenv);

        /*!  brief Do the ESP optimization
         *
         * Optimizes the charges using matrix inversion. No restraints are
         * taken into account, except total charge and charge symmetries.
         */
        void optimizeCharges();

        // Make sure the total charge is correct and that symmetry is obeyed
        void regularizeCharges();

        void potcomp(const std::string      &potcomp,
                     const std::string      &pdbdiff,
                     const gmx_output_env_t *oenv);

        //! Return the net charge for an atom
        double getAtomCharge(int atom) const;

        real myWeight(int iatom) const;

        void updateZeta(t_atoms *atoms, const Poldata &pd);

    private:
        //! Return the charge for one "shell" of an atom
        double getCharge(int atom, size_t zz) const;

        double getZeta(int atom, int zz) const;

        void setCharge(int atom, int zz, double q);

        void setZeta(int atom, int zz, double zeta);

        double calcPenalty();

        ChargeDistributionModel   iDistributionModel_;
        double                    watoms_;
        int                       qtot_;
        int                       qshell_;
        double                    rms_, rrms_, penalty_, pfac_, wtot_;
        dvec                      origin_, space_;
        bool                      bFitZeta_;
        bool                      bRandZeta_, bRandQ_;
        ivec                      nxyz_;
        real                      qfac_, zmin_, zmax_, deltaZ_, qmin_, qmax_, rDecrZeta_;
        int                       uniqueQ_;
        int                       fitQ_;
        int                       nAtom_;
        int                       nShell_;

        //! Total number of parameters
        std::vector<RespAtom>     ra_;
        std::vector<RespAtomType> ratype_;
        std::vector<std::string>  dzatoms_;
        std::string               stoichiometry_;
        std::vector<EspPoint>     ep_;
        std::vector<int>          symmetricAtoms_;
};

} // namespace

#endif
