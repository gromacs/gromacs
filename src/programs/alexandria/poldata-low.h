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
#ifndef POLDATA_LOW_H
#define POLDATA_LOW_H

#include "gmxpre.h"

#include <string>
#include <vector>

#include "gromacs/utility/exceptions.h"
#include "gromacs/utility/fatalerror.h"

#include "plistwrapper.h"
#include "coulombintegrals/coulombintegrals.h"

/*! \brief
 * Enumerated type holding the charge distribution models used in PolData
 *
 * \inpublicapi
 * \ingroup module_alexandria
 */
enum ChargeDistributionModel {
    eqdAXp, eqdAXg, eqdAXs, eqdYang, eqdBultinck, eqdRappe, eqdNR
};

namespace alexandria
{

const char *iType2string(InteractionType iType);

InteractionType string2iType(const char *string);


/*! \brief
 * Contains all the information realted to
 * polarizability.
 *
 * \inpublicapi
 * \ingroup module_alexandria
 */
class Ptype
{
    public:

        /*! \brief
         * Ptype constructor
         *
         * \param[in] ptype           Polarizability type
         * \param[in] miller          Miller equivalent of the polarizability type
         * \param[in] bosque          Bosque equivalent of the polarizability type
         * \param[in] polarizability  Polarizability value
         * \param[in] sigPol          Uncertainty in the calulated polarizability
         */
        Ptype(const std::string &ptype,
              const std::string &miller,
              const std::string &bosque,
              double             polarizability,
              double             sigPol);

        /*! \brief
         * Return the polarizability type
         */
        const std::string &getType() const { return type_; }

        /*! \brief
         * Return Miller equivalent of the polarizability type
         */
        const std::string &getMiller() const { return miller_; }

        /*! \brief
         * Return Bosque equivalent of the polarizability type
         */
        const std::string &getBosque() const { return bosque_; }

        /*! \brief
         * Set the polarizability value
         *
         * \param[in] polarizability  Polarizability value
         */
        void setPolarizability(double polarizability) { polarizability_ = polarizability; }

        /*! \brief
         * Return the polarizability value
         *
         */
        double getPolarizability() const { return polarizability_; }

        /*! \brief
         * Set the uncertainty in the calulated polarizability
         *
         * \param[in] sigPol  uncertainty
         */
        void setSigPol(double sigPol) { sigPol_ = sigPol; }

        /*! \brief
         * Return the uncertainty in the calulated polarizability
         *
         */
        double getSigPol() const { return sigPol_; }

    private:
        std::string type_;
        std::string miller_;
        std::string bosque_;
        double      polarizability_;
        double      sigPol_;


};

using PtypeIterator      = typename std::vector<Ptype>::iterator;
using PtypeConstIterator = typename std::vector<Ptype>::const_iterator;

/*! \brief
 * Contains all the information realted to
 * alexandria force field atom types.
 *
 * \inpublicapi
 * \ingroup module_alexandria
 */
class Ffatype
{
    public:

        /*! \brief
         * Fftype constructor
         *
         * \param[in] desc        Description
         * \param[in] type        Atom type
         * \param[in] ptype       Polarizability type
         * \param[in] btype       Bond type
         * \param[in] elem        Element name
         * \param[in] vdwparams   Van der Waals parameters
         * \param[in] refEnthalpy Reference Enthalpy of Formation
         */
        Ffatype(const std::string &desc,
                const std::string &type,
                const std::string &ptype,
                const std::string &btype,
                const std::string &elem,
                const std::string &vdwparams,
                const std::string &refEnthalpy);

        /*! \brief
         * Fftype default constructor
         */
        Ffatype () {}

        /*! \brief
         * Return the decription of atoms
         */
        const std::string &getDesc() const { return desc_; }

        /*! \brief
         * Return the type of atoms
         */
        const std::string &getType() const { return type_; }

        /*! \brief
         * Return the polarizability type of atoms
         */
        const std::string &getPtype() const { return ptype_; }

        /*! \brief
         * Return the bond type of atoms
         */
        const std::string &getBtype() const { return btype_; }

        /*! \brief
         * Return the element name of atoms
         */
        const std::string &getElem() const { return elem_; }

        /*! \brief
         * Return the Van der Waals parameters of atoms
         */
        const std::string &getVdwparams() const { return vdwparams_; }

        /*! \brief
         * Return the reference enthalpy of formation of atoms
         */
        const std::string &getRefEnthalpy() const { return refEnthalpy_; }

    private:
        std::string desc_;
        std::string type_;
        std::string ptype_;
        std::string btype_;
        std::string elem_;
        std::string vdwparams_;
        std::string refEnthalpy_;
};

using FfatypeIterator      = typename std::vector<Ffatype>::iterator;
using FfatypeConstIterator = typename std::vector<Ffatype>::const_iterator;


/*! \brief
 * Contains the listed force including
 * bonds or angles, or dihedrals
 *
 * \inpublicapi
 * \ingroup module_alexandria
 */
class ListedForce
{

    public:

        /*! \brief
         * ListedForce constructor
         *
         * \param[in] atoms       A vector of atom bond types invloved in the listed force
         * \param[in] params      Listed force paramters
         * \param[in] refvalue    The reference value such as the reference bond length
         * \param[in] sigma       Uncertainty in the calculated reference value
         * \param[in] ntrain      Number of molecules in the training set
         */
        ListedForce(const std::vector<std::string> atoms,
                    std::string                    params,
                    double                         refValue,
                    double                         sigma,
                    size_t                         ntrain);

        /*! \brief
         * Return the vector of atom bond types invloved in the listed force
         */
        const std::vector<std::string> &atoms() const { return atoms_; }

        /*! \brief
         * Set the listed force paramters
         */
        void setParams(const std::string &params) { params_ = params; }

        /*! \brief
         * Return the listed force paramters
         */
        const std::string &params() const { return params_; }

        /*! \brief
         * Set reference value
         */
        void setRefValue(double value) { refValue_ = value; }

        /*! \brief
         * Return reference value
         */
        double refValue() const { return refValue_; }

        /*! \brief
         * Set uncertainty in reference value
         */
        void setSigma(double sigma) { sigma_ = sigma; }

        /*! \brief
         * Return uncertainty in reference value
         */
        double sigma() const { return sigma_; }

        /*! \brief
         * Set the number of molecules in the training set
         */
        void setNtrain(int ntrain) { ntrain_ = ntrain; }

        /*! \brief
         * Return the number of molecules in the training set
         */
        int ntrain() const { return ntrain_; }

    private:

        const std::vector<std::string> atoms_;
        std::string                    params_;
        double                         refValue_;
        double                         sigma_;
        size_t                         ntrain_;
};

using ListedForceIterator      = typename std::vector<ListedForce>::iterator;
using ListedForceConstIterator = typename std::vector<ListedForce>::const_iterator;


/*! \brief
 * Contains the listed forces including
 * bonds and angles, and dihedrals
 *
 * \inpublicapi
 * \ingroup module_alexandria
 */
class ListedForces
{

    public:

        /*! \brief
         * ListedForces constructor
         *
         * \param[in] iType       Interaction type
         * \param[in] function    The function used to calculate the force corresponding to the iType
         * \param[in] Unit        The unit of the reference value
         */
        ListedForces(const std::string      iType,
                     const std::string     &function,
                     const std::string     &unit);

        /* \brief
         * Return the unit of the reference value
         */
        const std::string &unit() const {return unit_; }

        /* \brief
         * Return the function name
         */
        const std::string &function() const {return function_; }

        /* \brief
         * Return the function type
         */
        const unsigned int &fType() const {return fType_; }

        /* \brief
         * Return the interaction type
         */
        const InteractionType &iType() const {return iType_; }

        /* \brief
         * Find a listed force in the force_ vector
         *
         * \param[in] atoms   A vector of atom bond types invloved in the listed force
         */
        ListedForceIterator findForce(const std::vector<std::string> &atoms);

        /* \brief
         * Find a listed force in the force_ vector
         *
         * \param[in] atoms   A vector of atom bond types invloved in the listed force
         */
        ListedForceConstIterator findForce(const std::vector<std::string> &atoms) const;

        /* \brief
         * Set the listed force paramters
         *
         * \param[in] atoms       A vector of atom bond types invloved in the listed force
         * \param[in] params      Listed force paramters
         * \param[in] refvalue    The reference value such as the reference bond length
         * \param[in] sigma       Uncertainty in the calculated reference value
         * \param[in] ntrain      Number of molecules in the training set
         */
        bool setForceParams(const std::vector<std::string> &atoms,
                            const std::string              &params,
                            double                          refValue,
                            double                          sigma,
                            size_t                          ntrain);

        /* \brief
         * Generate a ListedForce object and
         * add it to the force_ vector, if it does not exist already
         *
         * \param[in] atoms       A vector of atom bond types invloved in the listed force
         * \param[in] params      Listed force paramters
         * \param[in] refvalue    The reference value such as the reference bond length
         * \param[in] sigma       Uncertainty in the calculated reference value
         * \param[in] ntrain      Number of molecules in the training set
         */
        void addForce(const std::vector<std::string> &atoms,
                      const std::string              &params,
                      double                          refValue,
                      double                          sigma,
                      size_t                          ntrain);

        /* \brief
         * Search a listed force in the force_ vector
         *
         * \param[in] atoms       A vector of atom bond types invloved in the listed force
         * \param[in] params      Listed force paramters
         * \param[in] refvalue    The reference value such as the reference bond length
         * \param[in] sigma       Uncertainty in the calculated reference value
         * \param[in] ntrain      Number of molecules in the training set
         */
        bool searchForce(std::vector<std::string> &atoms,
                         std::string              &params,
                         double                   *refValue,
                         double                   *sigma,
                         size_t                   *ntrain) const;


        /*! \brief
         * Erase the force_ vector
         */
        void eraseListedForce()
        {
            force_.clear();
        }

        ListedForceIterator forceBegin() {return force_.begin(); }

        ListedForceConstIterator forceBegin() const {return force_.begin(); }

        ListedForceIterator forceEnd() {return force_.end(); }

        ListedForceConstIterator forceEnd() const {return force_.end(); }

    private:

        InteractionType          iType_;
        const std::string        function_;
        const std::string        unit_;
        std::vector<ListedForce> force_;
        unsigned int             fType_;
};

using ListedForcesIterator      = typename std::vector<ListedForces>::iterator;
using ListedForcesConstIterator = typename std::vector<ListedForces>::const_iterator;


class Bosque
{
    public:

        Bosque(const std::string &bosque, double polarizability);

        const std::string &getBosque() const { return bosque_; }

        double getPolarizability() const { return polarizability_; }

    private:
        std::string bosque_;
        double      polarizability_;

};

using BosqueIterator      = typename std::vector<Bosque>::iterator;
using BosqueConstIterator = typename std::vector<Bosque>::const_iterator;

class Miller
{
    public:
        Miller(const std::string &miller,
               int                atomnumber,
               double             tauAhc,
               double             alphaAhp,
               const std::string &alexandria_equiv);

        const std::string &getMiller() const { return miller_; }

        int getAtomnumber() const { return atomnumber_; }

        double getTauAhc() const { return tauAhc_; }

        double getAlphaAhp() const { return alphaAhp_; }

        const std::string &getAlexandriaEquiv() const { return alexandria_equiv_; }

    private:
        //! Atom type name
        std::string miller_;
        //! Atomic number
        int         atomnumber_;
        //! Polarizability description tau
        double      tauAhc_;
        //! Polarizability description alpha
        double      alphaAhp_;
        //! Alexandria type
        std::string alexandria_equiv_;
};

using MillerIterator      = typename std::vector<Miller>::iterator;
using MillerConstIterator = typename std::vector<Miller>::const_iterator;

class Symcharges
{
    public:
        Symcharges(const std::string &central,
                   const std::string &attached,
                   int                numattach);

        const std::string &getCentral() const { return central_; }

        const std::string &getAttached() const { return attached_; }

        int getNumattach() const { return numattach_; }

    private:
        const std::string central_;
        const std::string attached_;
        int               numattach_;
};

using SymchargesIterator      = typename std::vector<Symcharges>::iterator;
using SymchargesConstIterator = typename std::vector<Symcharges>::const_iterator;

class Epref
{
    public:
        Epref(ChargeDistributionModel  eqdModel,
              const std::string       &epref);

        ChargeDistributionModel getEqdModel() const { return eqdModel_; }

        const char *getEpref() const { return epref_.c_str(); }

        void setEpref(std::string epref) { epref_ = epref; }

    private:
        ChargeDistributionModel eqdModel_;
        std::string             epref_;
};

using EprefIterator      = typename std::vector<Epref>::iterator;
using EprefConstIterator = typename std::vector<Epref>::const_iterator;

class RowZetaQ
{
    public:

        RowZetaQ(int row, double zeta, double q);

        int row() const { return row_; };

        void setRow(int row) { row_ = row; }

        double q() const { return q_; }

        void setQ(double q)
        {
            GMX_RELEASE_ASSERT(!fixedQ_, "Trying to modify a fixed charge");
            q_ = q;
        }

        double zeta() const { return zeta_; }

        void setZeta(double z) { zeta_ = z; }

        double zetaRef() const { return zetaRef_; }

        void setZetaRef(double z) { zetaRef_ = z; }

        int zIndex() const { return zindex_; }

        void setZindex(int zi) { zindex_ = zi; }

        bool fixedQ() const { return fixedQ_; }

    private:
        //! The row in the periodic table for each of the charge components
        int    row_;
        //! Inverse screening length of each of the components
        double zeta_;
        //! Charge of each of the components
        double q_;
        //! Reference (starting) value for zeta
        double zetaRef_;
        //! Parameter optimization index
        int    zindex_;
        //! Determines whether a charge is fixed (not to be changed)
        bool   fixedQ_;
};

using RowZetaQIterator      = typename std::vector<RowZetaQ>::iterator;
using RowZetaQConstIterator = typename std::vector<RowZetaQ>::const_iterator;

class Eemprops
{
    public:
        Eemprops(ChargeDistributionModel eqdModel,
                 const std::string      &name,
                 const std::string      &zetastr,
                 const std::string      &qstr,
                 const std::string      &rowstr,
                 double                  J0,
                 double                  chi0);

        ChargeDistributionModel getEqdModel() const { return eqdModel_; }

        int getNzeta() const { return rzq_.size(); }

        const char *getName() const { return name_.c_str(); }

        const char *getZetastr() const { return zetastr_.c_str(); }

        const char *getQstr() const { return qstr_.c_str(); }

        const char *getRowstr() const { return rowstr_.c_str(); }

        double getJ0() const { return J0_; }

        double getChi0() const { return chi0_; }

        void setEqdModel(ChargeDistributionModel eqdModel) { eqdModel_ = eqdModel; }

        void setName(const std::string &name) { name_ = name; }

        void setRowZetaQ(const std::string &rowstr,
                         const std::string &zetastr,
                         const std::string &qstr);

        void setJ0(double J0) { J0_ = J0; }

        void setChi0(double chi0) { chi0_ = chi0; }

        double getZeta(int index) const { return rzq_[index].zeta(); }

        double getQ(int index) const { return rzq_[index].q(); }

        int getRow(int index) const { return rzq_[index].row(); }

        void setZeta(int index, double zeta) { rzq_[index].setZeta(zeta); }

        void setQ(int index, double q) { rzq_[index].setQ(q); }

        void setRow(int index, int row) { rzq_[index].setRow(row); }

    private:
        ChargeDistributionModel eqdModel_;
        std::string             name_;
        std::string             rowstr_;
        std::string             zetastr_;
        std::string             qstr_;
        double                  J0_;
        double                  chi0_;
        std::vector<RowZetaQ>   rzq_;
};

using EempropsIterator      = typename std::vector<Eemprops>::iterator;
using EempropsConstIterator = typename std::vector<Eemprops>::const_iterator;

} // namespace aleaxndria
#endif
