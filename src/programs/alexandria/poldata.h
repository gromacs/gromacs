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
#ifndef POLDATA_H
#define POLDATA_H

#include <algorithm>
#include <vector>

#include "gromacs/utility/stringutil.h"

#include "poldata_low.h"
#include "stringutil.h"

/* This source code file is part of the Alexandria project */

struct t_commrec;

namespace alexandria
{

class Poldata
{
    public:

        Poldata() {};

        /*! \brief
         * Set the file name gentop.dat
         *
         */
        void  setFilename(const std::string &fn2);

        /*! \brief
         * Set the force field
         */
        void setForceField(const std::string &forcefield)
        {
            alexandriaForcefield_ = forcefield;
        }

        /*! \brief
         * Set the potential energy function for VDW interaction.
         * The VDW potentials supported by Alexandria are LJ, Buckingham, and Wang_Buckingham.
         *
         * \param[in] func  The name of the VDW potential function
         */
        void  setVdwFunction(const std::string &func);

        /*! \brief
         * Set the combination rule
         *
         */
        void setCombinationRule(const std::string &func);

        /*! \brief
         * Set the number of exclusion
         */
        void setNexcl(int nexcl) { nexcl_ = nexcl; }

        /*! \brief
         * Set the scaling factor for 1-4 electrostatic interaction
         */
        void setFudgeQQ(double fudgeQQ) { fudgeQQ_ = fudgeQQ; }

        /*! \brief
         * Set the scaling factor for 1-4 LJ interaction
         */
        void setFudgeLJ(double fudgeLJ) { fudgeLJ_ = fudgeLJ; }

        /*! \brief
         * Add the atom types used in Alexandria FF
         *
         **\param[in] elem          The element of the atom type
         **\param[in] desc          The description of the atom type
         **\param[in] atype         The atom type defined for this element in Alexandria FF
         **\param[in] ptype         The polarizability type defined for this element in Alexandria FF
         **\param[in] btype         The bond type defined for elem in Alexandria FF
         **\param[in] ztype         The zeta type defined for elem in Alexandria FF
         **\param[in] vdwparams     The VDW parameters for elem in Alexandria FF. The number of VDW parameters is 2 for LJ and 3 for
         *                         Buckingham and Wang_Buckingham potentials
         **\param[in] ref_enthalpy  The reference enthalpy for elem
         */
        void  addAtype(const std::string &elem,
                       const std::string &desc,
                       const std::string &atype,
                       const std::string &ptype,
                       const std::string &btype,
                       const std::string &ztype,
                       std::string       &vdwparams,
                       const std::string &ref_enthalpy);

        /*! \brief
         *  Add the polarizability types
         *
         **\param[in] ptype           The name specifying the polarizability type in Alexandria FF
         **\param[in] miller          Miller polarizability type
         **\param[in] bosque          Bosque polarizability type
         **\param[in] polarizability  The calulated value of polarizability
         **\param[in] sigPol          The uncertainty of the calculated polarizability
         */
        void  addPtype(const std::string &ptype,
                       const std::string &miller,
                       const std::string &bosque,
                       double             polarizability,
                       double             sigPol);

        void  addVsite(const std::string &atype,
                       const std::string &type,
                       int                number,
                       double             distance,
                       double             angle,
                       int                ncontrolatoms);

        /*! \brief
         * Set the value and the associated error for the given poltype
         *
         **\param[in] ptype            Polarizabilty type
         **\param[in] polarizability   The value of polarizabilty
         **\param[in] sigPol           The error
         */
        bool setPtypePolarizability(const std::string &ptype,
                                    double             polarizability,
                                    double             sigPol);

        /*! \brief
         * Set the unit of polarizability.
         */
        void setPolarUnit(const std::string &polarUnit)
        {
            alexandriaPolarUnit_ = polarUnit;
        }

        /*! \brief
         * Set the reference polarizability value.
         */
        void setPolarRef(const std::string &polarRef)
        {
            alexandriaPolarRef_ = polarRef;
        }

        /*! \brief
         * Set the vsite angle unit.
         */
        void setVsite_angle_unit(const std::string &angle_unit)
        {
            vsite_angle_unit_ = angle_unit;
        }

        /*! \brief
         * Set the vsite angle unit.
         */
        void setVsite_length_unit(const std::string &length_unit)
        {
            vsite_length_unit_ = length_unit;
        }

        std::vector<Vsite> &getVsite() {return vsite_; }

        const std::string &getForceField() const { return alexandriaForcefield_; }

        int getVdwFtype() const { return gtVdwFtype_; }

        int getNexcl() const { return nexcl_; }

        size_t getNatypes() const { return alexandria_.size(); }

        size_t getNptypes() const { return ptype_.size(); }

        double getFudgeQQ() const { return fudgeQQ_; }

        double getFudgeLJ() const { return fudgeLJ_; }

        /*! \brief
         * Return the reference enthalpy for the given atom type
         *
         **\param[in] atype Atom type
         **\param[ou] Href  Reference enthalpy
         */
        bool getAtypeRefEnthalpy(const std::string &atype,
                                 double            *Href) const;

        /*! \brief
         * Return the combination rule
         *
         */
        const std::string &getCombinationRule() const { return gtCombinationRule_; }

        /*! \brief
         * Return the combination rule.
         */
        int  getCombRule() const { return gtCombRule_; }

        std::string  getGeometry(std::string gtBrule);


        /*! \brief
         * Return the discription corresponding to the atom type
         *
         * \param[in] atype  Atom Type
         */
        const std::string &getDesc(const std::string &atype) const;

        /*! \brief
         * Return the element name corresponding to the atom type
         *
         * \param[in] atype  Atom Type
         */
        const std::string &getElem(const std::string &atype) const;
        
        /*! \brief
         * Return the element name corresponding to the zeta type
         *
         * \param[in] ztype  zeta Type
         */
        const std::string &ztype2elem(const std::string &ztype) const;


        /*! \brief
         * Return the charge corresponding to the atyom type
         * from the gentop.dat file
         *
         * \param[in] atype  Atom type
         */
        std::string  getCharge(  std::string atype);
        
        std::vector<Ffatype> &getAtypes() {return alexandria_; }

        FfatypeIterator getAtypeBegin() { return alexandria_.begin(); }

        FfatypeIterator getAtypeEnd() { return alexandria_.end(); }

        FfatypeConstIterator getAtypeBegin() const { return alexandria_.begin(); }

        FfatypeConstIterator getAtypeEnd() const { return alexandria_.end(); }

        /*! \brief
         * Return the iterator corresponding to the atom type
         *
         * \param[in] atype  Atom Type
         */
        FfatypeIterator findAtype(const std::string &atype)
        {
            return std::find_if(alexandria_.begin(), alexandria_.end(),
                                [atype](Ffatype const &f)
                                { return (atype == f.getType()); });
        }

        /*! \brief
         * Return the const_iterator corresponding to the atom type
         *
         * \param[in] atype  Atom Type
         */
        FfatypeConstIterator findAtype(const std::string &atype) const
        {
            return std::find_if(alexandria_.begin(), alexandria_.end(),
                                [atype](Ffatype const &f)
                                { return (atype == f.getType()); });
        }

        /*! \brief
         * Return the atom type corresponding to the bond type
         *
         * \param[in] btype  Bond Type
         */
        FfatypeIterator btypeToAtype(const std::string &btype)
        {
            return std::find_if(alexandria_.begin(), alexandria_.end(),
                                [btype](Ffatype const &f)
                                { return (btype == f.getBtype()); });
        }

        /*! \brief
         * Return true if a given bond type exists in Alexandria
         *
         * \param[in] btype  Bond Type
         */
        bool haveBtype(const std::string &btype)
        {
            return (btypeToAtype(btype) != alexandria_.end());
        }

        /*! \brief
         * Return the atom type corresponding to the polarizability type
         *
         * \param[in] ptype  Polarizability Type
         */
        FfatypeConstIterator ptypeToAtype(const std::string &ptype) const
        {
            return std::find_if(alexandria_.begin(), alexandria_.end(),
                                [ptype](Ffatype const &f)
                                { return (ptype == f.getPtype()); });
        }

        PtypeConstIterator getPtypeBegin() const { return ptype_.begin(); }

        PtypeConstIterator getPtypeEnd() const { return ptype_.end(); }

        VsiteIterator getVsiteBegin()  { return vsite_.begin(); }

        VsiteConstIterator getVsiteBegin()  const { return vsite_.begin(); }

        VsiteIterator getVsiteEnd() { return vsite_.end(); }

        VsiteConstIterator getVsiteEnd() const { return vsite_.end(); }

        VsiteIterator findVsite(std::string  atype)
        {

            return std::find_if(vsite_.begin(), vsite_.end(),
                                [atype](const Vsite &vs)
                                {
                                    return (atype == vs.atype());
                                });
        }

        VsiteConstIterator findVsite(std::string atype) const
        {

            return std::find_if(vsite_.begin(), vsite_.end(),
                                [atype](const Vsite &vs)
                                {
                                    return (atype == vs.atype());
                                });
        }

        /*! \brief
         * Return the iterator corresponding to the polarizability type
         *
         * \param[in] ptype  Polarizability Type
         */
        PtypeIterator findPtype(const std::string &ptype)
        {
            return std::find_if(ptype_.begin(), ptype_.end(),
                                [ptype](Ptype const &p)
                                { return (ptype == p.getType()); });
        }


        /*! \brief
         * Return the poltype corresponding to atype and true if successful
         *
         * \param[in]  atype  Atom type
         * \param[out] ptype  Polarizability type.
         */
        bool atypeToPtype(const std::string &atype,
                          std::string       &ptype) const;


        /*! \brief
         * Return the bond type corresponding to atom type and true if successful
         *
         * \param[in]  atype  Atom type
         * \param[out] btype  Polarizability type.
         */
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

        /*! \brief
         * Convert poltype to miller name. Return true if found
         */
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

        /*! \brief
         * Convert poltype to bosque name.  Return true if found.
         */
        bool ptypeToBosque(const std::string &ptype,
                           std::string       &bosque) const;

        bool getBosquePol(const std::string &bosque,
                          double            *polarizability) const;

        /* Return true on success or false otherwise */


        void addForces(ListedForces forces) {forces_.push_back(forces); }

        size_t nforces() const { return forces_.size(); }

        ListedForces &lastForces() { return forces_.back(); }

        const ListedForces &lastForces() const { return forces_.back(); }

        ListedForcesIterator forcesBegin() { return forces_.begin(); }

        ListedForcesIterator forcesEnd() { return forces_.end(); }

        ListedForcesConstIterator forcesBegin() const { return forces_.begin(); }

        ListedForcesConstIterator forcesEnd() const { return forces_.end(); }

        ListedForcesIterator findForces(InteractionType iType)
        {
            return std::find_if(forces_.begin(), forces_.end(),
                                [iType](ListedForces const &forces)
                                { return (iType == forces.iType()); });
        }

        ListedForcesConstIterator findForces(InteractionType iType) const
        {
            return std::find_if(forces_.begin(), forces_.end(),
                                [iType](ListedForces const &forces)
                                {return (iType == forces.iType()); });
        }

        bool findForce(std::vector<std::string> &atoms,
                       ListedForceIterator      *force);


        bool findForce(const std::vector<std::string> &atoms,
                       ListedForceConstIterator       *force) const;


        bool searchForce(std::vector<std::string> &atoms,
                         std::string              &params,
                         double                   *refValue,
                         double                   *sigma,
                         size_t                   *ntrain) const;

        bool searchForce(std::vector<std::string> &atoms,
                         std::string              &params,
                         double                   *refValue,
                         double                   *sigma,
                         size_t                   *ntrain,
                         InteractionType           iType) const;

        const std::string &getVdwFunction() const { return gtVdwFunction_; }

        const std::string &getPolarUnit() const { return alexandriaPolarUnit_; }

        const std::string &getPolarRef() const { return alexandriaPolarRef_; }

        const std::string &getVsite_angle_unit() const { return vsite_angle_unit_; }

        const std::string &getVsite_length_unit() const { return vsite_length_unit_; }

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
                     const std::string      &atype) const;

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
                                      const std::string       &atype) const;

        EempropsIterator findEem(ChargeDistributionModel  eqdModel,
                                 const std::string       &atype);

        std::vector<Eemprops> &getEemprops() {return eep_; }

        void  setEpref(ChargeDistributionModel eqdModel,
                       const std::string      &epref);

        const char *getEpref(ChargeDistributionModel eqdModel) const;

        //! Spread from master to slave nodes
        void  broadcast(const t_commrec *cr);

        EprefConstIterator epRefBegin() const { return epr_.begin(); }

        EprefConstIterator epRefEnd() const { return epr_.end(); }

        CommunicationStatus Send(const t_commrec *cr, int dest);

        CommunicationStatus Receive(const t_commrec *cr, int src);

    private:
        std::string                           filename_;
        std::vector<Ptype>                    ptype_;
        std::vector<Ffatype>                  alexandria_;
        std::vector<Vsite>                    vsite_;
        std::vector<std::string>              btype_;
        std::string                           alexandriaPolarUnit_;
        std::string                           alexandriaPolarRef_;
        std::string                           alexandriaForcefield_;
        std::string                           vsite_angle_unit_;
        std::string                           vsite_length_unit_;
        int                                   nexcl_;
        double                                fudgeQQ_, fudgeLJ_;
        std::string                           gtVdwFunction_, gtCombinationRule_;
        int                                   gtVdwFtype_, gtCombRule_;
        std::vector<ListedForces>             forces_;
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

}
#endif
