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

namespace alexandria
{

class Ptype
{
    public:
        Ptype(const std::string &ptype,
              const std::string &miller,
              const std::string &bosque,
              double             polarizability,
              double             sigPol);

        const std::string &getType() const { return type_; }

        const std::string &getMiller() const { return miller_; }

        const std::string &getBosque() const { return bosque_; }

	void setPolarizability(double polarizability) { polarizability_ = polarizability; }

        double getPolarizability() const { return polarizability_; }    

	void setSigPol(double sigPol) { sigPol_ = sigPol; }

        double getSigPol() const { return sigPol_; }

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


};

using PtypeIterator      = typename std::vector<Ptype>::iterator;
using PtypeConstIterator = typename std::vector<Ptype>::const_iterator;

class Ffatype
{
    public:
        Ffatype(const std::string &desc,
                const std::string &type,
                const std::string &ptype,
                const std::string &btype,
                const std::string &elem,
                const std::string &vdwparams,
                const std::string &refEnthalpy);

        Ffatype () {}

        const std::string &getDesc() const { return desc_; }

        const std::string &getType() const { return type_; }

        const std::string &getPtype() const { return ptype_; }

        const std::string &getBtype() const { return btype_; }

        const std::string &getElem() const { return elem_; }

        const std::string &getVdwparams() const { return vdwparams_; }

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

class GtBond
{
    public:
        GtBond(const std::string atom1,
               const std::string atom2,
               const std::string params,
               double            length,
               double            sigma,
               double            bondorder,
               int               ntrain);

        const std::string &getAtom1() const { return atom1_; }

        const std::string &getAtom2() const { return atom2_; }

	void setParams(const std::string &params) { params_ = params; }

        const std::string &getParams() const { return params_; }   

        void setLength(double length) { length_ = length; }

        double getLength() const { return length_; }

	void setSigma(double sigma) { sigma_ = sigma; }

        double getSigma() const { return sigma_; }

	void setBondorder(double bondorder) { bondorder_ = bondorder; }

        double getBondorder() const { return bondorder_; }

	void setNtrain(int ntrain) { ntrain_ = ntrain; }

        int getNtrain() const { return ntrain_; }
	
    private:
        const std::string atom1_;
        const std::string atom2_;
        std::string       params_;
        std::string       elem1_;
	std::string       elem2_;
        double            length_;
        double            sigma_;
        double            bondorder_;
        int               ntrain_;
};

using GtBondIterator      = typename std::vector<GtBond>::iterator;
using GtBondConstIterator = typename std::vector<GtBond>::const_iterator;

class GtBonds 
{
    public:

        GtBonds(const std::string function, const std::string unit);

	const std::string &getLengthUnit() const { return gtLengthUnit_;}

	void setBondFtype();

	const std::string &getBondFunction() const{ return gtBondFunction_;}

	unsigned int getBondFtype() const { return gtBondFtype_; }

        size_t getNgtBond() const { return gtBond_.size(); }

	GtBondIterator getBondBegin() { return gtBond_.begin(); }  

	GtBondIterator getBondEnd() { return gtBond_.end(); }

	GtBondConstIterator getBondBegin() const { return gtBond_.begin(); }

	GtBondConstIterator getBondEnd() const { return gtBond_.end(); }

	GtBondIterator findBond(const std::string &btype1,
				const std::string &btype2,
				double             bondorder);
                                
        GtBondConstIterator findBond(const std::string &btype1,
				     const std::string &btype2,
				     double             bondorder) const;

	bool setBondParams(const std::string &btype1,
                           const std::string &btype2,
                           double             length, 
                           double             sigma,
                           int                ntrain,
                           double             bondorder, 
                           const std::string &params);

	void addBond(const std::string &btype1,
                     const std::string &btype2,
                     double             length, 
                     double             sigma, 
                     int                ntrain,
                     double             bondorder,
                     const std::string &params);   

	bool searchBond(const std::string &btype1,
                        const std::string &btype2,
                        double           *length, 
                        double           *sigma, 
                        int              *ntrain,
                        double           *bondorder, 
                        std::string      &params) const;
	
    private:
	std::vector<GtBond> gtBond_;
        const std::string   gtBondFunction_;
        const std::string   gtLengthUnit_; 
	unsigned int        gtBondFtype_;   
};

using GtBondsIterator      = typename std::vector<GtBonds>::iterator;
using GtBondsConstIterator = typename std::vector<GtBonds>::const_iterator;

class GtAngle
{
    public:
        GtAngle(const std::string &atom1,
                const std::string &atom2,
                const std::string &atom3,
                const std::string &params,
                double             angle,
                double             sigma,
                int                ntrain);

        const std::string &getAtom1() const { return atom1_; }

        const std::string &getAtom2() const { return atom2_; }

        const std::string &getAtom3() const { return atom3_; }
	
	void setParams(std::string params) { params_ = params; }

        const std::string &getParams() const { return params_; }

	void setAngle(double angle) { angle_ = angle; }

        double getAngle() const { return angle_; }

	void setSigma(double sigma) { sigma_ = sigma; }

        double getSigma() const { return sigma_; }

	void setNtrain(int ntrain) { ntrain_ = ntrain; }

        int getNtrain() const { return ntrain_; }

    private:
        const std::string atom1_;
        const std::string atom2_;
        const std::string atom3_;
        std::string       params_;
        double            angle_;
        double            sigma_;
        int               ntrain_;
};

using GtAngleIterator      = typename std::vector<GtAngle>::iterator;
using GtAngleConstIterator = typename std::vector<GtAngle>::const_iterator;

class GtAngles
{
    public:

        GtAngles(const std::string function,const std::string unit);

	const std::string &getAngleUnit() const { return gtAngleUnit_;}

	void setAngleFtype();

        const std::string &getAngleFunction() const { return gtAngleFunction_; }

	unsigned int getAngleFtype() const { return gtAngleFtype_; }

	size_t getNgtAngle() const { return gtAngle_.size(); }

	GtAngleIterator getAngleBegin() { return gtAngle_.begin(); }

	GtAngleIterator getAngleEnd() { return gtAngle_.end(); }

	GtAngleConstIterator getAngleBegin() const { return gtAngle_.begin(); }

        GtAngleConstIterator getAngleEnd() const { return gtAngle_.end(); }

	GtAngleIterator findAngle(const std::string &btype1,
                                  const std::string &btype2,
                                  const std::string &btype3);
                                  
        GtAngleConstIterator findAngle(const std::string &btype1,
                                       const std::string &btype2,
                                       const std::string &btype3) const;

	bool setAngleParams(const std::string &btype1,
                            const std::string &btype2,
                            const std::string &btype3, 
                            double             angle, 
                            double             sigma, 
                            int                ntrain, 
                            const std::string &params);

	void addAngle(const std::string &btype1, 
                      const std::string &btype2,
                      const std::string &btype3,
                      double             angle, 
                      double             sigma,
                      int                ntrain,
                      const std::string &params);

	bool searchAngle(const std::string &btype1,
                         const std::string &btype2,
                         const std::string &btype3,
                         double            *angle, 
                         double            *sigma,
                         int               *ntrain, 
                         std::string       &params) const;
         
    private:
	std::vector<GtAngle> gtAngle_;
        const std::string    gtAngleFunction_;
	const std::string    gtAngleUnit_;
        unsigned int         gtAngleFtype_;   
};

using GtAnglesIterator      = typename std::vector<GtAngles>::iterator;
using GtAnglesConstIterator = typename std::vector<GtAngles>::const_iterator;

class GtDihedral
{
    public:
        GtDihedral(const std::string &atom1,
                   const std::string &atom2,
                   const std::string &atom3,
                   const std::string &atom4,
                   const std::string &params,
                   double             dihedral,
                   double             sigma,
                   int                ntrain);

        const std::string &getAtom1() const { return atom1_; }

        const std::string &getAtom2() const { return atom2_; }

        const std::string &getAtom3() const { return atom3_; }

        const std::string &getAtom4() const { return atom4_; }

        void setParams(const std::string &params) { params_ = params; }

        const std::string &getParams() const { return params_; }

	void setDihedral(double dihedral) { dihedral_ = dihedral; }

        double getDihedral() const { return dihedral_; }

	void setSigma(double sigma) { sigma_ = sigma; }

        double getSigma() const { return sigma_; }
	
	void setNtrain(int ntrain) { ntrain_ = ntrain; }

        int getNtrain() const { return ntrain_; }

    private:
        const std::string atom1_;
        const std::string atom2_;
        const std::string atom3_;
        const std::string atom4_;
        std::string       params_;
        double            dihedral_;
        double            sigma_;
        int               ntrain_;
};

using GtDihedralIterator      = typename std::vector<GtDihedral>::iterator;
using GtDihedralConstIterator = typename std::vector<GtDihedral>::const_iterator;

class GtDihedrals
{
    public:

        GtDihedrals(const std::string function, const std::string unit);

	const std::string &getDihedralUnit() const {return gtDihedralUnit_; }

	void setDihedralFtype();

        const std::string &getDihedralFunction() const { return gtDihedralFunction_; }

	unsigned int getDihedralFtype() const { return gtDihedralFtype_; }

	size_t getNgtDihedral() const { return gtDihedral_.size(); }

	GtDihedralIterator getDihedralBegin() { return gtDihedral_.begin(); }

        GtDihedralIterator getDihedralEnd() { return gtDihedral_.end(); }

        GtDihedralConstIterator getDihedralBegin() const { return gtDihedral_.begin(); }

        GtDihedralConstIterator getDihedralEnd() const { return gtDihedral_.end(); }

	GtDihedralIterator findDihedral(const std::string &btype1, 
					const std::string &btype2,
					const std::string &btype3, 
					const std::string &btype4);

        GtDihedralConstIterator findDihedral(const std::string &btype1, 
					     const std::string &btype2,
					     const std::string &btype3, 
					     const std::string &btype4) const;

	bool setDihedralParams(const std::string &btype1,
                               const std::string &btype2,
                               const std::string &btype3, 
                               const std::string &btype4,
                               double             angle, 
                               double             sigma,
                               int                ntrain, 
                               const std::string &params);


	void addDihedral(const std::string &btype1,
                         const std::string &btype2,
                         const std::string &btype3, 
                         const std::string &btype4,
                         double             dihedral, 
                         double             sigma,
                         int                ntrain,
                         const std::string &params);

        
	bool searchDihedral(const std::string &btype1, 
                            const std::string &btype2,
                            const std::string &btype3,
                            const std::string &btype4,
                            double            *dihedral, 
                            double            *sigma,
                            int               *ntrain, 
                            std::string       &params) const;
         
   private:
	std::vector<GtDihedral> gtDihedral_;
        const std::string       gtDihedralFunction_;
	const std::string       gtDihedralUnit_;
        unsigned int            gtDihedralFtype_; 
};

using GtDihedralsIterator      = typename std::vector<GtDihedrals>::iterator;
using GtDihedralsConstIterator = typename std::vector<GtDihedrals>::const_iterator;

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
