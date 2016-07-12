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
 * \author David van der Spoel <david.vanderspoel@icm.uu.se>
 */
#ifndef POLDATA_H
#define POLDATA_H

#include "gmxpre.h"

#include <algorithm>
#include <vector>

#include "gromacs/utility/smalloc.h"
#include "gromacs/utility/stringutil.h"

#include "poldata-low.h"
#include "stringutil.h"

/* This source code file is part of the Alexandria project */

struct t_commrec;

namespace alexandria
{

class Poldata
{
    public:

        Poldata() {};

        void  setFilename(const std::string &fn2);

        void  setVdwFunction(const std::string &func);

        void  addPtype(const std::string &ptype,
                       const std::string &miller,
                       const std::string &bosque,
                       double             polarizability,
                       double             sigPol);

        void  addAtype(const std::string &elem,
                       const std::string &desc,
                       const std::string &atype,
                       const std::string &ptype,
                       const std::string &btype,
                       const std::string &vdwparams,
                       const std::string &ref_enthalpy);

        bool setPtypePolarizability(const std::string &ptype,
                                    double             polarizability,
                                    double             sigPol);

        void setPolarUnit(const std::string &polarUnit)
        {
            alexandriaPolarUnit_ = polarUnit;
        }

        void setPolarRef(const std::string &polarRef)
        {
            alexandriaPolarRef_ = polarRef;
        }

        const std::string &getForceField() const { return alexandriaForcefield_; }

        void setForceField(const std::string &forcefield)
        {
            alexandriaForcefield_ = forcefield;
        }

        int getVdwFtype() const { return gtVdwFtype_; }

        void setNexcl(int nexcl) { nexcl_ = nexcl; }

        int getNexcl() const { return nexcl_; }

        void setFudgeQQ(double fudgeQQ) { fudgeQQ_ = fudgeQQ; }

        size_t getNatypes() const { return alexandria_.size(); }

        size_t getNptypes() const { return ptype_.size(); }

        size_t getNgtBonds() const { return gtBonds_.size(); }

        size_t getNgtAngles() const { return gtAngles_.size(); }

        size_t getNgtDihedrals() const { return gtDihedrals_.size(); }

        double getFudgeQQ() const { return fudgeQQ_; }

        void setFudgeLJ(double fudgeLJ) { fudgeLJ_ = fudgeLJ; }

        double getFudgeLJ() const { return fudgeLJ_; }

        bool getAtypeRefEnthalpy(const std::string &atype,
                                 double            *Href) const;

        void setCombinationRule(const std::string &func);

        const std::string &getCombinationRule() const { return gtCombinationRule_; }

        int  getCombRule() const { return gtCombRule_; }

        std::string  getGeometry(  std::string gtBrule);

        std::string  getDesc(  std::string atype);

        /* Get the charge from the gentop.dat file */
        std::string  getCharge(  std::string atype);

        FfatypeIterator getAtypeBegin() { return alexandria_.begin(); }

        FfatypeIterator getAtypeEnd() { return alexandria_.end(); }

        FfatypeConstIterator getAtypeBegin() const { return alexandria_.begin(); }

        FfatypeConstIterator getAtypeEnd() const { return alexandria_.end(); }

        FfatypeIterator findAtype(const std::string &atype)
        {
            return std::find_if(alexandria_.begin(), alexandria_.end(),
                                [atype](Ffatype const &f)
                                { return (atype.compare(f.getType()) == 0); });
        }

        FfatypeConstIterator findAtype(const std::string &atype) const
        {
            return std::find_if(alexandria_.begin(), alexandria_.end(),
                                [atype](Ffatype const &f)
                                { return (atype.compare(f.getType()) == 0); });
        }

        FfatypeIterator btypeToAtype(const std::string &btype)
        {
            return std::find_if(alexandria_.begin(), alexandria_.end(),
                                [btype](Ffatype const &f)
                                { return (f.getBtype().compare(btype) == 0); });
        }

        bool haveBtype(const std::string &btype)
        {
            return (btypeToAtype(btype) != alexandria_.end());
        }

        FfatypeIterator ptypeToAtype(const std::string &ptype)
        {
            return std::find_if(alexandria_.begin(), alexandria_.end(),
                                [ptype](Ffatype const &f)
                                { return (f.getPtype().compare(ptype) == 0); });
        }

        PtypeConstIterator getPtypeBegin() const { return ptype_.begin(); }

        PtypeConstIterator getPtypeEnd() const { return ptype_.end(); }

        PtypeIterator findPtype(const std::string &ptype)
        {
            return std::find_if(ptype_.begin(), ptype_.end(),
                                [ptype](Ptype const &p)
                                { return (ptype.compare(p.getType()) == 0); });
        }

        //! Return the poltype corresponding to atype and true if successful
        bool atypeToPtype(const std::string &atype,
                          std::string       &ptype) const;

        //! Return the bondtype corresponding to atype and true if successful
        bool atypeToBtype(const std::string &atype,
                          std::string       &btype) const;

        /* Return 1 if OK, 0 if not found */
        bool getPtypePol(const std::string &ptype,
                         double *polarizability, double *sigPol) const;
        bool getAtypePol(const std::string &atype,
                         double *polarizability, double *sigPol) const;

        void addMiller(const std::string &miller,
                       int                atomnumber,
                       double             tauAhc,
                       double             alphaAhp,
                       const std::string &alexandria_equiv);

        /* Return true if "miller" was found */
        bool getMillerPol(const std::string &miller,
                          int               *atomnumber,
                          double            *tauAhc,
                          double            *alphaAhp,
                          std::string       &alexandria_equiv) const;

        MillerIterator getMillerBegin() { return miller_.begin(); }

        MillerIterator getMillerEnd() { return miller_.end(); }

        MillerConstIterator getMillerBegin() const { return miller_.begin(); }

        MillerConstIterator getMillerEnd() const { return miller_.end(); }

        void setMillerFlags(const std::string &tauUnit,
                            const std::string &ahpUnit,
                            const std::string &ref)
        {
            millerTauUnit_ = tauUnit;
            millerAhpUnit_ = ahpUnit;
            millerRef_     = ref;
        }

        void getMillerFlags(std::string &tauUnit,
                            std::string &ahpUnit,
                            std::string &ref) const
        {
            tauUnit = millerTauUnit_;
            ahpUnit = millerAhpUnit_;
            ref     = millerRef_;
        }

        //! Convert poltype to miller name. Return true if found
        bool ptypeToMiller(const std::string &ptype,
                           std::string       &mil_type) const;

        void  addBosque(const std::string &bosque,
                        double             polarizability)
        {
            Bosque bos(bosque, polarizability);
            bosque_.push_back(bos);
        }

        BosqueIterator getBosqueBegin() { return bosque_.begin(); }

        BosqueIterator getBosqueEnd() { return bosque_.end(); }

        BosqueConstIterator getBosqueBegin() const { return bosque_.begin(); }

        BosqueConstIterator getBosqueEnd() const { return bosque_.end(); }

        void setBosqueFlags(const std::string &polarUnit,
                            const std::string &ref)
        {
            bosquePolarUnit_ = polarUnit;
            bosqueRef_       = ref;
        }

        void getBosqueFlags(std::string &polarUnit,
                            std::string &ref) const
        {
            polarUnit = bosquePolarUnit_;
            ref       = bosqueRef_;
        }

        //! Convert poltype to bosque name or nullptr if not found
        bool ptypeToBosque(const std::string &ptype,
                           std::string       &bosque) const;

        bool getBosquePol(const std::string &bosque,
                          double            *polarizability) const;

        /* Return true on success or false otherwise */
        void addBond(GtBonds           *gtbs,
                     const std::string &atom1,
                     const std::string &atom2,
                     double             length,
                     double             sigma,
                     int                ntrain,
                     double             bondorder,
                     const std::string &params);

        void addGtBonds(GtBonds gtbs) {gtBonds_.push_back(gtbs); }

        GtBonds &getLastGtBonds() { return gtBonds_.back(); }

        const GtBonds &getLastGtBonds() const { return gtBonds_.back(); }

        GtBondsIterator getBondsBegin() { return gtBonds_.begin(); }

        GtBondsIterator getBondsEnd() { return gtBonds_.end(); }

        GtBondsConstIterator getBondsBegin() const { return gtBonds_.begin(); }

        GtBondsConstIterator getBondsEnd() const { return gtBonds_.end(); }

        GtBondsIterator findGtBonds(const std::string bondFunction)
        {
            return std::find_if(gtBonds_.begin(), gtBonds_.end(),
                                [bondFunction](GtBonds const &gtbs)
                                { return (bondFunction.compare(gtbs.getBondFunction()) == 0); });
        }

        GtBondsConstIterator findGtBonds(const std::string bondFunction) const
        {
            return std::find_if(gtBonds_.begin(), gtBonds_.end(),
                                [bondFunction](GtBonds const &gtbs)
                                { return (bondFunction.compare(gtbs.getBondFunction()) == 0); });
        }

        bool findBond(const std::string &btype1,
                      const std::string &btype2,
                      double             bondorder,
                      GtBondIterator    *gtb);

        bool findBond(const std::string &btype1,
                      const std::string &btype2,
                      double             bondorder,
                      GtBondIterator    *gtb,
                      int               *index);

        bool findBond(const std::string   &btype1,
                      const std::string   &btype2,
                      double               bondorder,
                      GtBondConstIterator *gtb) const;

        bool findBond(const std::string   &btype1,
                      const std::string   &btype2,
                      double               bondorder,
                      GtBondConstIterator *gtb,
                      int                 *index) const;

        bool searchBond(const std::string &btype1,
                        const std::string &btype2,
                        double            *length,
                        double            *sigma,
                        int               *ntrain,
                        double            *bondorder,
                        std::string       &params) const;

        const std::string &getVdwFunction() const { return gtVdwFunction_; }

        const std::string &getPolarUnit() const { return alexandriaPolarUnit_; }

        const std::string &getPolarRef() const { return alexandriaPolarRef_; }

        void addAngle(GtAngles          *gtas,
                      const std::string &btype1,
                      const std::string &btype2,
                      const std::string &btype3,
                      double             angle,
                      double             sigma,
                      int                ntrain,
                      const std::string &params);

        void addGtAngles(GtAngles gtas) {gtAngles_.push_back(gtas); }

        GtAngles &getLastGtAngles() { return gtAngles_.back(); }

        GtAnglesIterator getAnglesBegin() { return gtAngles_.begin(); }

        GtAnglesIterator getAnglesEnd() { return gtAngles_.end(); }

        GtAnglesConstIterator getAnglesBegin() const { return gtAngles_.begin(); }

        GtAnglesConstIterator getAnglesEnd() const { return gtAngles_.end(); }

        GtAnglesIterator findGtAngles(const std::string angleFunction)
        {
            return std::find_if(gtAngles_.begin(), gtAngles_.end(),
                                [angleFunction](GtAngles const &gtas)
                                { return (angleFunction.compare(gtas.getAngleFunction()) == 0); });
        }

        GtAnglesConstIterator findGtAngles(const std::string angleFunction) const
        {
            return std::find_if(gtAngles_.begin(), gtAngles_.end(),
                                [angleFunction](GtAngles const &gtas)
                                { return (angleFunction.compare(gtas.getAngleFunction()) == 0); });
        }

        bool findAngle(const std::string &btype1,
                       const std::string &btype2,
                       const std::string &btype3,
                       GtAngleIterator   *gta);

        bool findAngle(const std::string &btype1,
                       const std::string &btype2,
                       const std::string &btype3,
                       GtAngleIterator   *gta,
                       int               *index);

        bool findAngle(const std::string    &btype1,
                       const std::string    &btype2,
                       const std::string    &btype3,
                       GtAngleConstIterator *gta) const;

        bool findAngle(const std::string    &btype1,
                       const std::string    &btype2,
                       const std::string    &btype3,
                       GtAngleConstIterator *gta,
                       int                  *index) const;

        bool searchAngle(const std::string &btype1,
                         const std::string &btype2,
                         const std::string &btype3,
                         double            *angle,
                         double            *sigma,
                         int               *ntrain,
                         std::string       &params) const;

        void addDihedral(GtDihedrals       *gtds,
                         const std::string &btype1,
                         const std::string &btype2,
                         const std::string &btype3,
                         const std::string &btype4,
                         double             dihedral,
                         double             sigma,
                         int                ntrain,
                         const std::string &params);

        void addGtDihedrals(GtDihedrals gtds) {gtDihedrals_.push_back(gtds); }

        GtDihedrals &getLastGtDihedrals() { return gtDihedrals_.back(); }

        GtDihedralsIterator getDihedralsBegin() { return gtDihedrals_.begin(); }

        GtDihedralsIterator getDihedralsEnd() { return gtDihedrals_.end(); }

        GtDihedralsConstIterator getDihedralsBegin() const { return gtDihedrals_.begin(); }

        GtDihedralsConstIterator getDihedralsEnd() const { return gtDihedrals_.end(); }

        GtDihedralsIterator findGtDihedrals(const std::string dihedralFunction)
        {
            return std::find_if(gtDihedrals_.begin(), gtDihedrals_.end(),
                                [dihedralFunction](GtDihedrals const &gtds)
                                { return (dihedralFunction.compare(gtds.getDihedralFunction()) == 0); });
        }

        GtDihedralsConstIterator findGtDihedrals(const std::string dihedralFunction) const
        {
            return std::find_if(gtDihedrals_.begin(), gtDihedrals_.end(),
                                [dihedralFunction](GtDihedrals const &gtds)
                                { return (dihedralFunction.compare(gtds.getDihedralFunction()) == 0); });
        }

        bool findDihedral(const std::string  &btype1,
                          const std::string  &btype2,
                          const std::string  &btype3,
                          const std::string  &btype4,
                          GtDihedralIterator *gtd);

        bool findDihedral(const std::string  &btype1,
                          const std::string  &btype2,
                          const std::string  &btype3,
                          const std::string  &btype4,
                          GtDihedralIterator *gtd,
                          int                *index);

        bool findDihedral(const std::string       &btype1,
                          const std::string       &btype2,
                          const std::string       &btype3,
                          const std::string       &btype4,
                          GtDihedralConstIterator *gtd) const;

        bool findDihedral(const std::string       &btype1,
                          const std::string       &btype2,
                          const std::string       &btype3,
                          const std::string       &btype4,
                          GtDihedralConstIterator *gtd,
                          int                     *index) const;

        bool searchDihedral(const std::string &btype1,
                            const std::string &btype2,
                            const std::string &btype3,
                            const std::string &btype4,
                            double            *dihedral,
                            double            *sigma,
                            int               *ntrain,
                            std::string       &params) const;

        void eraseBonded()
        {
            gtBonds_.clear();
            gtAngles_.clear();
            gtDihedrals_.clear();
        }

        void addSymcharges(const std::string &central,
                           const std::string &attached,
                           int                numattach);

        SymchargesIterator getSymchargesBegin() { return symcharges_.begin(); }

        SymchargesIterator getSymchargesEnd() { return symcharges_.end(); }

        SymchargesConstIterator getSymchargesBegin() const { return symcharges_.begin(); }

        SymchargesConstIterator getSymchargesEnd() const { return symcharges_.end(); }

        int getNumprops(ChargeDistributionModel eqdModel) const;

        int havePolSupport(const std::string &atype) const;

        bool haveEemSupport(ChargeDistributionModel  eqdModel,
                            const std::string       &name,
                            gmx_bool                 bAllowZeroParameters) const;

        double getJ00(ChargeDistributionModel  eqdModel,
                      const std::string       &name) const;

        int getNzeta(ChargeDistributionModel eqdModel,
                     const std::string      &name) const;

        double getZeta(ChargeDistributionModel eqdModel,
                       const std::string &name, int zz) const;

        const char *getQstr(ChargeDistributionModel  eqdModel,
                            const std::string       &name) const;

        const char *getRowstr(ChargeDistributionModel  eqdModel,
                              const std::string       &name) const;

        double getQ(ChargeDistributionModel eqdModel,
                    const std::string      &name,
                    int                     zz) const;

        int getRow(ChargeDistributionModel eqdModel,
                   const std::string      &name,
                   int                     zz) const;

        double getChi0(ChargeDistributionModel eqdModel,
                       const std::string      &name) const;

        const char *getOpts(ChargeDistributionModel eqdModel,
                            const std::string      &name) const;

        void  addEemprops(Eemprops eep) { eep_.push_back(eep); }

        EempropsConstIterator BeginEemprops() const { return eep_.begin(); }

        EempropsConstIterator EndEemprops() const { return eep_.end(); }

        EempropsIterator BeginEemprops() { return eep_.begin(); }

        EempropsIterator EndEemprops() { return eep_.end(); }

        EempropsConstIterator findEem(ChargeDistributionModel  eqdModel,
                                      const std::string       &name) const;

        EempropsIterator findEem(ChargeDistributionModel  eqdModel,
                                 const std::string       &name);

        void  setEpref(ChargeDistributionModel eqdModel,
                       const std::string      &epref);

        const char *getEpref(ChargeDistributionModel eqdModel) const;

        //! Spread from master to slave nodes
        void  broadcast(t_commrec *cr);

        EprefConstIterator epRefBegin() const { return epr_.begin(); }

        EprefConstIterator epRefEnd() const { return epr_.end(); }

    private:
        std::string                           filename_;
        std::vector<Ptype>                    ptype_;
        std::vector<Ffatype>                  alexandria_;
        std::vector<std::string>              btype_;
        std::string                           alexandriaPolarUnit_;
        std::string                           alexandriaPolarRef_;
        std::string                           alexandriaForcefield_;
        int                                   nexcl_;
        double                                fudgeQQ_, fudgeLJ_;
        std::string                           gtVdwFunction_, gtCombinationRule_;
        int                                   gtVdwFtype_, gtCombRule_;
        std::vector<GtBonds>                  gtBonds_;
        std::vector<GtAngles>                 gtAngles_;
        std::vector<GtDihedrals>              gtDihedrals_;
        std::vector<Miller>                   miller_;
        std::string                           millerTauUnit_, millerAhpUnit_;
        std::string                           millerRef_;
        std::vector<Bosque>                   bosque_;
        std::string                           bosquePolarUnit_;
        std::string                           bosqueRef_;
        std::vector<Symcharges>               symcharges_;
        std::vector<Eemprops>                 eep_;
        std::vector<Epref>                    epr_;

        void addBtype(const std::string &btype);

        gmx_bool strcasestrStart(std::string needle, std::string haystack);

        template<class Type>
        int indexOfPointInVector(Type * pointer, std::vector<Type> vector)
        {
            return (pointer - &(vector[0]));
        }
};

const char *getEemtypeName(ChargeDistributionModel eem);

ChargeDistributionModel name2eemtype(const std::string name);

}
#endif
