/*! \internal \brief
 * Implements part of the alexandria program.
 * \author David van der Spoel <david.vanderspoel@icm.uu.se>
 */
#ifndef POLDATA_LOW_H
#define POLDATA_LOW_H

#include "gmxpre.h"

#include <string>
#include <vector>

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
        double      refEnthalpy_;
    public:
        Ffatype(const std::string &desc,
                const std::string &type,
                const std::string &ptype,
                const std::string &btype,
                const std::string &elem,
                const std::string &vdwparams,
                double             refEnthalpy) :
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

        double getRefEnthalpy() const { return refEnthalpy_; }
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
               double      length,
               double      sigma,
               double      bondorder,
               int         ntrain)
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
                double      angle,
                double      sigma,
                int         ntrain)
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
                   double      dihedral,
                   double      sigma,
                   int         ntrain)
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
        std::string miller_;
        int         atomnumber_;
        double      tauAhc_;
        double      alphaAhp_;
    public:
        Miller(const std::string &miller,
               int               atomnumber,
               double            tauAhc,
               double            alphaAhp)
            :
              miller_(miller),
              atomnumber_(atomnumber),
              tauAhc_(tauAhc),
              alphaAhp_(alphaAhp) {}

        const std::string &getMiller() const { return miller_; }

        int getAtomnumber() const { return atomnumber_; }

        double getTauAhc() const { return tauAhc_; }

        double getAlphaAhp() const { return alphaAhp_; }
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
 
class EempropsData
{
    private:
        std::vector<double> q_;
        std::vector<double> zeta_;
        std::vector<int>    row_;
    public:

        EempropsData(int size) : q_(size), zeta_(size), row_(size) {}

        //EempropsData() {}

        double getQ(int index) const { return q_[index]; }

        const std::vector<double> &getQ() const { return q_; }

        double getZeta(int index) const { return zeta_[index]; }

        const std::vector<double> &getZeta() const { return zeta_; }

        int getRow(int index) const { return row_[index]; }

        const std::vector<int> &getRow() const { return row_; }

        void setZeta(double zeta, int index) { zeta_[index] = zeta; }

        void setRow(double row, int index) { row_[index] = row; }

        void setQ(double q, int index) { q_[index] = q; }
};
typedef std::vector<EempropsData>::iterator EempropsDataIterator;

#define MAXZETA    12
class Eemprops
{
    private:
        ChargeDistributionModel eqdModel_;
        int                     nzeta_;
        std::string             name_;
        std::string             zetastr_;
        std::string             qstr_;
        std::string             rowstr_;
        double                  J0_;
        double                  chi0_;
        EempropsData            data_;

    public:
        Eemprops(ChargeDistributionModel eqdModel,
                 int                     nzeta,
                 const std::string      &name,
                 const std::string      &zetastr,
                 const std::string      &qstr,
                 const std::string      &rowstr,
                 double                  J0,
                 double                  chi0)
            : eqdModel_(eqdModel),
              nzeta_(nzeta),
              name_(name),
              zetastr_(zetastr),
              qstr_(qstr),
              rowstr_(rowstr),
              J0_(J0),
              chi0_(chi0),
              data_(MAXZETA) {}

        Eemprops() : data_(MAXZETA) {}

        ChargeDistributionModel getEqdModel() const { return eqdModel_; }

        int getNzeta() const { return nzeta_; }

        const char *getName() const { return name_.c_str(); }

        const char *getZetastr() const { return zetastr_.c_str(); }

        const char *getQstr() const { return qstr_.c_str(); }

        const char *getRowstr() const { return rowstr_.c_str(); }

        double getJ0() const { return J0_; }

        double getChi0() const { return chi0_; }

        void setEqdModel(ChargeDistributionModel eqdModel) { eqdModel_ = eqdModel; }

        void setNzeta(int nzeta) { nzeta_ = nzeta; }

        void setName(const std::string &name) { name_ = name; }

        void setZetastr(const std::string &zetastr) { zetastr_ = zetastr; }

        void setQstr(const std::string &qstr) { qstr_ = qstr; }

        void setRowstr(const std::string &rowstr) { rowstr_ = rowstr; }

        void setJ0(double J0) { J0_ = J0; }

        void setChi0(double chi0) { chi0_ = chi0; }

        double getZeta(int index) const { return data_.getZeta(index); }

        double getQ(int index) const { return data_.getQ(index); }

        int getRow(int index) const { return data_.getRow(index); }

        std::vector<double> getAllZeta() { return data_.getZeta(); }

        std::vector<double> getAllQ() { return data_.getQ(); }

        std::vector<int> getAllRow() { return data_.getRow(); }

        void setZeta(int index, double zeta) { data_.setZeta(zeta, index); }

        void setQ(int index, double q) { data_.setQ(q, index); }

        void setRow(int index, int row) { data_.setRow(row, index); }
};
typedef std::vector<Eemprops>::iterator EempropsIterator;
typedef std::vector<Eemprops>::const_iterator EempropsConstIterator;

} // namespace aleaxndria

#endif
