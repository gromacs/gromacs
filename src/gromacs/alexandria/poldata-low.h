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

enum DihedrlalType {
    egdPDIHS, egdIDIHS, egdNR
};

namespace alexandria
{

class Ptype
{
    private:
        //! Polarizability type
        std::string type_;
        //! Miller equivalent
        std::string miller_;
        //! Bosque equivalent
        std::string bosque_;
        //! Polarizability value
        double      polarizability_;
        //! Standard deviation
        double      sigPol_;

    public:
        Ptype(const std::string &ptype,
              const std::string &miller,
              const std::string &bosque,
              double             polarizability,
              double             sigPol) :
            type_(ptype),
            miller_(miller),
            bosque_(bosque),
            polarizability_(polarizability),
            sigPol_(sigPol)
        {}

        const std::string &getType() const { return type_; }

        const std::string &getMiller() const { return miller_; }

        const std::string &getBosque() const { return bosque_; }

        double getPolarizability() const { return polarizability_; }

        void setPolarizability(double polarizability) { polarizability_ = polarizability; }

        double getSigPol() const { return sigPol_; }

        void setSigPol(double sigPol) { sigPol_ = sigPol; }

};

class Ffatype
{
    private:
        std::string desc_;
        std::string type_;
        std::string ptype_;
        std::string btype_;
        std::string elem_;
        std::string vdwparams_;
        std::string refEnthalpy_;
    public:
        Ffatype(const std::string &desc,
                const std::string &type,
                const std::string &ptype,
                const std::string &btype,
                const std::string &elem,
                const std::string &vdwparams,
                const std::string &refEnthalpy) :
            desc_(desc),
            type_(type),
            ptype_(ptype),
            btype_(btype),
            elem_(elem),
            vdwparams_(vdwparams),
            refEnthalpy_(refEnthalpy) {}

        Ffatype() {}

        const std::string &getDesc() const { return desc_; }

        const std::string &getType() const { return type_; }

        const std::string &getPtype() const { return ptype_; }

        const std::string &getBtype() const { return btype_; }

        const std::string &getElem() const { return elem_; }

        const std::string &getVdwparams() const { return vdwparams_; }

        const std::string &getRefEnthalpy() const { return refEnthalpy_; }
};

class GtBond
{
    private:
        std::string atom1_;
        std::string atom2_;
        std::string params_;
        std::string elem1_, elem2_;
        double      length_;
        double      sigma_;
        double      bondorder_;
        int         ntrain_;
    public:
        GtBond(const std::string atom1,
               const std::string atom2,
               const std::string params,
               const std::string elem1,
               const std::string elem2,
               double            length,
               double            sigma,
               double            bondorder,
               int               ntrain)
            :
              atom1_(atom1),
              atom2_(atom2),
              params_(params),
              elem1_(elem1),
              elem2_(elem2),
              length_(length),
              sigma_(sigma),
              bondorder_(bondorder),
              ntrain_(ntrain)
        {}

        //GtBond() {}

        const std::string &getAtom1() const { return atom1_; }

        const std::string &getAtom2() const { return atom2_; }

        const std::string &getParams() const { return params_; }

        const std::string &getElem1() const { return elem1_; }

        const std::string &getElem2() const { return elem2_; }

        void setElem1(const std::string elem) { elem1_ = elem; }

        void setElem2(const std::string elem) { elem2_ = elem; }

        double getLength() const { return length_; }

        double getSigma() const { return sigma_; }

        double getBondorder() const { return bondorder_; }

        int getNtrain() const { return ntrain_; }

        void setAtom1(const std::string &atom) { atom1_ = atom; }

        void setAtom2(const std::string &atom) { atom2_ = atom; }

        void setBondorder(double bondorder) { bondorder_ = bondorder; }

        void setLength(double length) { length_ = length; }

        void setParams(const std::string &params) { params_ = params; }

        void setSigma(double sigma) { sigma_ = sigma; }

        void setNtrain(int ntrain) { ntrain_ = ntrain; }
};

class GtAngle
{
    private:
        std::string atom1_;
        std::string atom2_;
        std::string atom3_;
        std::string params_;
        double      angle_;
        double      sigma_;
        int         ntrain_;
    public:
        GtAngle(const std::string &atom1,
                const std::string &atom2,
                const std::string &atom3,
                const std::string &params,
                double             angle,
                double             sigma,
                int                ntrain)
            :
              atom1_(atom1),
              atom2_(atom2),
              atom3_(atom3),
              params_(params),
              angle_(angle),
              sigma_(sigma),
              ntrain_(ntrain)
        {}

        //GtAngle() {}

        const std::string &getAtom1() const { return atom1_; }

        const std::string &getAtom2() const { return atom2_; }

        const std::string &getAtom3() const { return atom3_; }

        const std::string &getParams() const { return params_; }

        double getAngle() const { return angle_; }

        double getSigma() const { return sigma_; }

        int getNtrain() const { return ntrain_; }

        void setParams(std::string params) { params_ = params; }

        void setAngle(double angle) { angle_ = angle; }

        void setSigma(double sigma) { sigma_ = sigma; }

        void setNtrain(int ntrain) { ntrain_ = ntrain; }
};

class GtDihedral
{
    private:
        std::string atom1_;
        std::string atom2_;
        std::string atom3_;
        std::string atom4_;
        std::string params_;
        double      dihedral_;
        double      sigma_;
        int         ntrain_;
    public:
        GtDihedral(const std::string &atom1,
                   const std::string &atom2,
                   const std::string &atom3,
                   const std::string &atom4,
                   const std::string &params,
                   double             dihedral,
                   double             sigma,
                   int                ntrain)
            :
              atom1_(atom1),
              atom2_(atom2),
              atom3_(atom3),
              atom4_(atom4),
              params_(params),
              dihedral_(dihedral),
              sigma_(sigma),
              ntrain_(ntrain)
        {}

        //GtDihedral() {}

        const std::string &getAtom1() const { return atom1_; }

        const std::string &getAtom2() const { return atom2_; }

        const std::string &getAtom3() const { return atom3_; }

        const std::string &getAtom4() const { return atom4_; }

        void setAtom1(const std::string &atom) { atom1_ = atom; }

        void setAtom2(const std::string &atom) { atom2_ = atom; }

        void setAtom3(const std::string &atom) { atom3_ = atom; }

        void setAtom4(const std::string &atom) { atom4_ = atom; }

        const std::string &getParams() const { return params_; }

        double getDihedral() const { return dihedral_; }

        double getSigma() const { return sigma_; }

        int getNtrain() const { return ntrain_; }

        bool compare(const GtDihedral &gtB) const;

        void setParams(const std::string &params) { params_ = params; }

        void setDihedral(double dihedral) { dihedral_ = dihedral; }

        void setSigma(double sigma) { sigma_ = sigma; }

        void setNtrain(int ntrain) { ntrain_ = ntrain; }
};

class Bosque
{
    private:
        std::string bosque_;
        double      polarizability_;
    public:
        Bosque(const std::string &bosque,
               double             polarizability)
            :
              bosque_(bosque),
              polarizability_(polarizability)
        {}

        const std::string &getBosque() const { return bosque_; }

        double getPolarizability() const { return polarizability_; }
};

class Miller
{
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
    public:
        Miller(const std::string &miller,
               int                atomnumber,
               double             tauAhc,
               double             alphaAhp,
               const std::string &alexandria_equiv)
            :
              miller_(miller),
              atomnumber_(atomnumber),
              tauAhc_(tauAhc),
              alphaAhp_(alphaAhp),
              alexandria_equiv_(alexandria_equiv) {}

        const std::string &getMiller() const { return miller_; }

        int getAtomnumber() const { return atomnumber_; }

        double getTauAhc() const { return tauAhc_; }

        double getAlphaAhp() const { return alphaAhp_; }

        const std::string &getAlexandriaEquiv() const { return alexandria_equiv_; }
};

class Symcharges
{
    private:
        std::string central_;
        std::string attached_;
        int         numattach_;
    public:
        Symcharges(const std::string &central,
                   const std::string &attached,
                   int                numattach) :
            central_(central),
            attached_(attached),
            numattach_(numattach)
        {
        }

        const std::string &getCentral() const { return central_; }

        const std::string &getAttached() const { return attached_; }

        int getNumattach() const { return numattach_; }
};

class Epref
{
    private:
        ChargeDistributionModel eqdModel_;
        std::string             epref_;
    public:
        Epref(ChargeDistributionModel  eqdModel,
              const std::string       &epref)
            : eqdModel_(eqdModel),
              epref_(epref) {}

        ChargeDistributionModel getEqdModel() const { return eqdModel_; }

        const char *getEpref() const { return epref_.c_str(); }

        void setEpref(std::string epref) { epref_ = epref; }
};
//! Loop over EEMprop references
typedef std::vector<Epref>::iterator EprefIterator;
//! Loop over EEMprop references
typedef std::vector<Epref>::const_iterator EprefConstIterator;

class RowZetaQ
{
    public:
        RowZetaQ(int row, double zeta, double q) : row_(row), zeta_(zeta), q_(q),
                                                   zetaRef_(zeta)
        {
            zindex_ = -1;
            char buf[256];
            row_ = std::min(row_, SLATER_MAX);
            if (row_ < row && debug)
            {
                fprintf(debug, "Reducing row from %d to %d\n", row, row_);
            }
            snprintf(buf, sizeof(buf), "Row (%d) in the periodic table must be > 0 and <= %d",
                     row_, SLATER_MAX);
            GMX_RELEASE_ASSERT(row_ > 0 && row_ <= SLATER_MAX, buf);
            fixedQ_ = (q != 0);
        }

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
//! Loop over RowZetaQ
typedef std::vector<RowZetaQ>::iterator RowZetaQIterator;

//! Loop over RowZetaQ
typedef std::vector<RowZetaQ>::const_iterator RowZetaQConstIterator;

class Eemprops
{
    private:
        ChargeDistributionModel eqdModel_;
        std::string             name_;
        std::string             rowstr_;
        std::string             zetastr_;
        std::string             qstr_;
        double                  J0_;
        double                  chi0_;
        std::vector<RowZetaQ>   rzq_;

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
};
typedef std::vector<Eemprops>::iterator EempropsIterator;
typedef std::vector<Eemprops>::const_iterator EempropsConstIterator;

} // namespace aleaxndria

#endif
