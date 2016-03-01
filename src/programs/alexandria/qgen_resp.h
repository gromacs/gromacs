/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2015,2016, by the GROMACS development team, led by
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
 * \author David van der Spoel <david.vanderspoel@icm.uu.se>
 */
#ifndef GMX_RESP_H
#define GMX_RESP_H

#include <cstdio>

#include <vector>

#include "gromacs/math/vectypes.h"
#include "gromacs/random/random.h"
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
        vCalc_ = 0;
        rho_ = 0;
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
    double v_;
    //! The calculated potential
    double vCalc_;
    //! The electron density in the point
    double rho_;
};

class QgenResp
{
    public:
        QgenResp();

        ~QgenResp();

        ChargeDistributionModel chargeDistributionModel()
        { return _iDistributionModel; }

        /*! \brief Set options for ESP charge generation
         *
         * \param[in] c          Charge distribution model eqdAXp, eqdAXg or eqdAXs
         * \param[in] seed       Random number generator seed. If <= 0 a seed will be generated.
         * \param[in] fitZeta    Do we fit the zeta (AXg/AXs only)
         * \param[in] zetaMin    Minimum allowed zeta value
         * \param[in] zetaMax    Maximum allowed zeta value
         * \param[in] deltaZeta  Minimum difference between zeta values in one atom
         * \param[in] randomZeta Use random starting values for zeta optimization
         * \param[in] qmin       Minimum allowed atomic charge
         * \param[in] qmax       Maximum allowed atomic charge
         * \param[in] randomQ    Use random starting values for charge optimization
         * \param[in] watoms     Weighting factor for atoms in ESP fit
         */
        void setOptions(ChargeDistributionModel c,
                        unsigned int            seed,
                        bool                    fitZeta,
                        real                    zetaMin,
                        real                    zetaMax,
                        real                    deltaZeta,
                        bool                    randomZeta,
                        real                    qmin,
                        real                    qmax,
                        bool                    randomQ,
                        real                    watoms);

        void setBAXpRESP(bool bAXpRESP) { _bAXpRESP = bAXpRESP; }

        void setRDecrZeta(real rDecrZeta) { _rDecrZeta = rDecrZeta; }

        void setBEntropy(bool bEntropy) { _bEntropy  = bEntropy; }

        real getMolecularCharge() const { return _qtot; }

        void setMolecularCharge(int qtot) { _qtot = qtot; }

        // int atomicnumber2row(int elem);

        size_t nAtom() const { return ra_.size(); }

        size_t nAtomType() const { return ratype_.size(); }

        size_t nEsp() const { return ep_.size(); }

        size_t nParam() const { return raparam_.size(); }

        void statistics(int len, char buf[]);

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

        void setAtomInfo(t_atoms       *atoms,
                         const Poldata &pd,
                         const rvec     x[]);
                         
        void updateAtomCoords(const rvec x[]);

        const std::string &getStoichiometry() const { return _stoichiometry; }

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

        void potLsq(gmx_stats_t lsq);

        void calcRho();

        void calcPot();

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
         * taken into account.
         */
        void optimizeCharges();

        // Make sure the total charge is correct and that symmetry is obeyed
        void regularizeCharges();

        int optimizeZeta(int maxiter, real *rms);

        void potcomp(const std::string      &potcomp,
                     const std::string      &pdbdiff,
                     const gmx_output_env_t *oenv);

        //! Return the net charge for an atom
        double getAtomCharge(int atom) const;

        double calcPenalty();

        void getVector(double *params);

        real myWeight(int iatom) const;

    private:
        void setVector(double *params);
        //! Return the charge for one "shell" of an atom
        double getCharge(int atom, size_t zz) const;

        double getZeta(int atom, int zz) const;

        void setCharge(int atom, int zz, double q);

        void setZeta(int atom, int zz, double zeta);

        ChargeDistributionModel   _iDistributionModel;
        double                    _watoms;
        int                       _qtot;
        double                    _rms, _rrms, _penalty, _pfac, _entropy, _wtot;
        dvec                      _origin, _space;
        bool                      _bFitZeta, _bEntropy;
        bool                      _bRandZeta, _bRandQ;
        bool                      _bAXpRESP;
        ivec                      _nxyz;
        real                      _qfac, _bHyper, _zmin, _zmax, _deltaZ, _qmin, _qmax, _rDecrZeta;
        gmx_rng_t                 rnd_;
        int                       uniqueQ_;
        int                       fitQ_;
        int                       nAtom_;
        int                       nShell_;

        //! Total number of parameters
        std::vector<RespAtom>     ra_;
        std::vector<RespAtomType> ratype_;
        std::vector<RespParam>    raparam_;
        std::vector<std::string>  _dzatoms;
        std::string               _stoichiometry;
        std::vector<EspPoint>     ep_;
        std::vector<int>          symmetricAtoms_;

        void warning(const std::string fn, int line);

        /*! \brief Adds a parameter to the optimization list
         *
         * \param[in] atom  Either the atom index or the atom type
         * \param[in] eparm Either eparmQ, in which case atom is the
         *                  atom index, or eparmZ, in which case atom
         *                  is the atom type.
         * \param[in] zz    The zeta index
         * \return The index in the parameter list.
         */
        int addParam(size_t atom, eParm eparm, size_t zz);
};

} // namespace

#endif
