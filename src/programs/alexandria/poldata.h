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

typedef std::vector<Ptype>::iterator PtypeIterator;
typedef std::vector<Ptype>::const_iterator PtypeConstIterator;

typedef std::vector<Ffatype>::iterator FfatypeIterator;
typedef std::vector<Ffatype>::const_iterator FfatypeConstIterator;

typedef std::vector<GtBond>::iterator GtBondIterator;
typedef std::vector<GtBond>::const_iterator GtBondConstIterator;

typedef std::vector<GtAngle>::iterator GtAngleIterator;
typedef std::vector<GtAngle>::const_iterator GtAngleConstIterator;

typedef std::vector<GtDihedral>::iterator DihedralIterator;
typedef std::vector<GtDihedral>::const_iterator DihedralConstIterator;

typedef std::vector<Bosque>::iterator BosqueIterator;
typedef std::vector<Bosque>::const_iterator BosqueConstIterator;

typedef std::vector<Miller>::iterator MillerIterator;
typedef std::vector<Miller>::const_iterator MillerConstIterator;

typedef std::vector<Symcharges>::iterator SymchargesIterator;
typedef std::vector<Symcharges>::const_iterator SymchargesConstIterator;

class Poldata
{
    public:
        Poldata(); //Constructor
        ~Poldata(){}

        void  setFilename(const std::string &fn2);

        void  setVdwFunction(const std::string &func);

        void  addPtype(const std::string &ptype,
                       const std::string &miller,
                       const std::string &bosque,
                       double            polarizability,
                       double            sigPol);

        void  addAtype(const std::string &elem,
                       const std::string &desc,
                       const std::string &atype,
                       const std::string &ptype,
                       const std::string &btype,
                       const std::string &vdwparams,
                       const std::string &ref_enthalpy);

        bool setPtypePolarizability(const std::string &ptype,
                                    double polarizability, 
                                    double sigPol);

        void setPolarUnit(const std::string &polarUnit)
        {
            _alexandriaPolarUnit = polarUnit;
        }

        void setPolarRef(const std::string &polarRef)
        {
            _alexandriaPolarRef = polarRef;
        }

        const std::string &getForceField() const { return _alexandriaForcefield; }

        void setForceField(const std::string &forcefield)
        {
            _alexandriaForcefield = forcefield;
        }

        void setLengthUnit(const std::string length_unit)
        {
            _gtLengthUnit = length_unit;
        }

        int getVdwFtype() const { return _gtVdwFtype; }

        void setNexcl(int nexcl) { _nexcl = nexcl; }

        int getNexcl() const { return _nexcl; }

        void setFudgeQQ(double fudgeQQ) { _fudgeQQ = fudgeQQ; }

        size_t getNatypes() const { return _alexandria.size(); }

        size_t getNptypes() const { return _ptype.size(); }

        size_t getNgtBond() const { return _gtBond.size(); }

        size_t getNgtAngle() const { return _gtAngle.size(); }

        size_t getNgtDihedral(int egd) const { return _gtDihedral[egd].size(); }

        double getFudgeQQ() const { return _fudgeQQ; }

        void setFudgeLJ(double fudgeLJ) { _fudgeLJ = fudgeLJ; }

        double getFudgeLJ() const { return _fudgeLJ; }

        bool getAtypeRefEnthalpy(const std::string &atype,
                                 double            *Href) const;

        void setCombinationRule(const std::string &func);

        const std::string &getCombinationRule() const { return _gtCombinationRule; }

        int  getCombRule() const { return _gtCombRule; }

        const std::string &getLengthUnit() const { return _gtLengthUnit; }

        std::string  getGeometry(  std::string gtBrule);

        std::string  getDesc(  std::string atype);

        /* Get the charge from the gentop.dat file */
        std::string  getCharge(  std::string atype);

        FfatypeIterator getAtypeBegin() { return _alexandria.begin(); }

        FfatypeIterator getAtypeEnd() { return _alexandria.end(); }
        
        FfatypeConstIterator getAtypeBegin() const { return _alexandria.begin(); }

        FfatypeConstIterator getAtypeEnd() const { return _alexandria.end(); }
        
        FfatypeIterator findAtype(const std::string &atype)
        {
            return std::find_if(_alexandria.begin(), _alexandria.end(),
                                [atype](Ffatype const &f)
                                { return (atype.compare(f.getType()) == 0); });
        }

        FfatypeConstIterator findAtype(const std::string &atype) const
        {
            return std::find_if(_alexandria.begin(), _alexandria.end(),
                                [atype](Ffatype const &f)
                                { return (atype.compare(f.getType()) ==0); });
        }

        FfatypeIterator btypeToAtype(const std::string &btype)
        {
            return std::find_if(_alexandria.begin(), _alexandria.end(), 
                                [btype](Ffatype const &f) 
                                { return (f.getBtype().compare(btype) == 0); });
        }

        bool haveBtype(const std::string &btype)
        {
            return (btypeToAtype(btype) != _alexandria.end());
        }
        
        FfatypeIterator ptypeToAtype(const std::string &ptype)
        {
            return std::find_if(_alexandria.begin(), _alexandria.end(), 
                                [ptype](Ffatype const &f) 
                                { return (f.getPtype().compare(ptype) ==0); });
        }
        
        PtypeConstIterator getPtypeBegin() const { return _ptype.begin(); }

        PtypeConstIterator getPtypeEnd() const { return _ptype.end(); }
        
        PtypeIterator findPtype(const std::string &ptype)
        {
            return std::find_if(_ptype.begin(), _ptype.end(),
                                [ptype](Ptype const &p)
                                { return (ptype.compare(p.getType()) == 0); });
        }

        //! Return the poltype corresponding to atype and true if successful
        bool atypeToPtype(const std::string &atype,
                          std::string       &ptype) const;

        //! Return the poltype corresponding to atype and true if successful
        bool atypeToBtype(const std::string &atype,
                          std::string       &btype) const;
        
        /* Return 1 if OK, 0 if not found */
        bool getPtypePol(const std::string &ptype,
                         double *polarizability, double *sigPol) const;
        bool getAtypePol(const std::string &atype,
                         double *polarizability, double *sigPol) const ;

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

        MillerIterator getMillerBegin() { return _miller.begin(); }

        MillerIterator getMillerEnd() { return _miller.end(); }

        MillerConstIterator getMillerBegin() const { return _miller.begin(); }

        MillerConstIterator getMillerEnd() const { return _miller.end(); }

        void setMillerFlags(const std::string &tauUnit,
                            const std::string &ahpUnit,
                            const std::string &ref)
        {
            _millerTauUnit = tauUnit;
            _millerAhpUnit = ahpUnit;
            _millerRef     = ref;
        }

        void getMillerFlags(std::string &tauUnit,
                            std::string &ahpUnit,
                            std::string &ref) const
        {
            tauUnit = _millerTauUnit;
            ahpUnit = _millerAhpUnit;
            ref     = _millerRef;
        }

        //! Convert poltype to miller name. Return true if found
        bool ptypeToMiller(const std::string &ptype,
                           std::string       &mil_type) const;

        void  addBosque(const std::string &bosque,
                        double             polarizability)
        {
            Bosque bos(bosque, polarizability);
            _bosque.push_back(bos);
        }

        BosqueIterator getBosqueBegin() { return _bosque.begin(); }

        BosqueIterator getBosqueEnd() { return _bosque.end(); }

        BosqueConstIterator getBosqueBegin() const { return _bosque.begin(); }

        BosqueConstIterator getBosqueEnd() const { return _bosque.end(); }

        void setBosqueFlags(const std::string &polarUnit,
                            const std::string &ref)
        {
            _bosquePolarUnit = polarUnit;
            _bosqueRef      = ref;
        }

        void getBosqueFlags(std::string &polarUnit,
                            std::string &ref) const
        { 
            polarUnit = _bosquePolarUnit; 
            ref       = _bosqueRef;
        }

        //! Convert poltype to bosque name or nullptr if not found
        bool ptypeToBosque(const std::string &ptype,
                           std::string       &bosque) const;

        bool getBosquePol(const std::string &bosque, 
                          double            *polarizability) const;

        /* Return true on success or false otherwise */
        void addBond(const std::string &atom1,
                     const std::string &atom2,
                     double length, 
                     double sigma, 
                     int ntrain,
                     double bondorder,
                     const std::string &params);

        bool setBondParams(const std::string &atom1,
                           const std::string &atom2,
                           double length, 
                           double sigma,
                           int ntrain,
                           double bondorder, 
                           const std::string &params);

        GtBondIterator getBondBegin() { return _gtBond.begin(); }

        GtBondIterator getBondEnd() { return _gtBond.end(); }

        GtBondConstIterator getBondBegin() const { return _gtBond.begin(); }

        GtBondConstIterator getBondEnd() const { return _gtBond.end(); }

        void  setBondFunction(const std::string &fn);

        const std::string &getBondFunction() const
        {
            return _gtBondFunction;
        }

        GtBondIterator findBond(const std::string &atom1,
                                const std::string &atom2,
                                double bondorder);
                                
        GtBondConstIterator findBond(const std::string &atom1,
                                     const std::string &atom2,
                                     double bondorder) const;
            
        /* Return true or false */
        bool searchBond(const std::string &atom1,
                        const std::string &atom2,
                        double *length, 
                        double *sigma, 
                        int *ntrain,
                        double *bondorder, 
                        std::string &params) const;

        void setAngleFunction(const std::string &fn);

        const std::string &getAngleFunction() const { return _gtAngleFunction; }

        int getBondFtype() const { return _gtBondFtype; }

        int getAngleFtype() const { return _gtAngleFtype; }

        int getDihedralFtype(int egd) const { return _gtDihedralFtype[egd]; }

        const std::string &getVdwFunction() const { return _gtVdwFunction; }

        const std::string &getPolarUnit() const { return _alexandriaPolarUnit; }

        const std::string &getPolarRef() const { return _alexandriaPolarRef; }

        /* Return true on success, false otherwise */
        void addAngle(const std::string &atom1, 
                      const std::string &atom2,
                      const std::string &atom3,
                      double angle, 
                      double sigma,
                      int ntrain,
                      const std::string &params);

        bool setAngleParams(const std::string &atom1,
                            const std::string &atom2,
                            const std::string &atom3, 
                            double angle, 
                            double sigma, 
                            int ntrain, 
                            const std::string &params);

        GtAngleIterator getAngleBegin() { return _gtAngle.begin(); }

        GtAngleIterator getAngleEnd() { return _gtAngle.end(); }

        GtAngleConstIterator getAngleBegin() const { return _gtAngle.begin(); }

        GtAngleConstIterator getAngleEnd() const { return _gtAngle.end(); }

        GtAngleIterator findAngle(const std::string &atom1,
                                  const std::string &atom2,
                                  const std::string &atom3);
                                  
        GtAngleConstIterator findAngle(const std::string &atom1,
                                       const std::string &atom2,
                                       const std::string &atom3) const;
                                  
                                  
        /* Return true or false */
        bool searchAngle(const std::string &atom1,
                         const std::string &atom2,
                         const std::string &atom3,
                         double *angle, 
                         double *sigma,
                         int *ntrain, 
                         std::string &params) const;

        void setAngleUnit(const std::string &angleUnit) { _gtAngleUnit= angleUnit; }

        const std::string &getAngleUnit() const { return _gtAngleUnit; }

        void  setDihedralFunction(int egd, const std::string &fn);

        const std::string &getDihedralFunction(int egd) const
        {
            return _gtDihedralFunction[egd];
        }

        /* Return 1 on success or 0 otherwise */
        void addDihedral(int egd, 
                         const std::string &atom1,
                         const std::string &atom2,
                         const std::string &atom3, 
                         const std::string &atom4,
                         double dihedral, 
                         double sigma,
                         int ntrain,
                         const std::string &params);

        bool setDihedralParams(int egd,
                               const std::string &atom1,
                               const std::string &atom2,
                               const std::string &atom3, 
                               const std::string &atom4,
                               double angle, 
                               double sigma,
                               int ntrain, 
                               const std::string &params);

        DihedralIterator getDihedralBegin(int egd) { return _gtDihedral[egd].begin(); }

        DihedralIterator getDihedralEnd(int egd) { return _gtDihedral[egd].end(); }

        DihedralConstIterator getDihedralBegin(int egd) const { return _gtDihedral[egd].begin(); }

        DihedralConstIterator getDihedralEnd(int egd) const { return _gtDihedral[egd].end(); }

        /* Return iterator */
        DihedralIterator findDihedral(int egd,
                                      const std::string &atom1, 
                                      const std::string &atom2,
                                      const std::string &atom3, 
                                      const std::string &atom4);

        DihedralConstIterator findDihedral(int egd,
                                           const std::string &atom1, 
                                           const std::string &atom2,
                                           const std::string &atom3, 
                                           const std::string &atom4) const;

        bool searchDihedral(int egd, 
                            const std::string &atom1, 
                            const std::string &atom2,
                            const std::string &atom3,
                            const std::string &atom4,
                            double *dihedral, 
                            double *sigma,
                            int *ntrain, 
                            std::string &params) const;

        const std::string &getDihedralUnit(int egd) const;

        void addSymcharges(const std::string &central,
                           const std::string &attached, 
                           int numattach);

        SymchargesIterator getSymchargesBegin() { return _symcharges.begin(); }

        SymchargesIterator getSymchargesEnd() { return _symcharges.end(); }

        SymchargesConstIterator getSymchargesBegin() const { return _symcharges.begin(); }

        SymchargesConstIterator getSymchargesEnd() const { return _symcharges.end(); }

        int getNumprops(ChargeDistributionModel eqdModel) const;

        int havePolSupport(const std::string &atype) const;

        bool haveEemSupport(ChargeDistributionModel  eqdModel,
                            const std::string       &name,
                            gmx_bool                 bAllowZeroParameters) const;

        double getJ00(ChargeDistributionModel  eqdModel, 
                      const std::string       &name) const;

        int getNzeta(ChargeDistributionModel eqdModel, 
                     const std::string &name) const;

        double getZeta(ChargeDistributionModel eqdModel, 
                       const std::string &name, int zz) const;

        const char *getQstr(ChargeDistributionModel  eqdModel, 
                            const std::string       &name) const;

        const char *getRowstr(ChargeDistributionModel  eqdModel, 
                              const std::string       &name) const;

        double getQ(ChargeDistributionModel eqdModel,
                    const std::string &name, 
                    int zz) const;

        int getRow(ChargeDistributionModel eqdModel,
                   const std::string &name,
                   int zz) const;

        double getChi0(ChargeDistributionModel eqdModel,
                       const std::string &name) const;

        const char *getOpts(ChargeDistributionModel eqdModel, 
                            const std::string &name) const;

        void  addEemprops(Eemprops eep) { _eep.push_back(eep); }
        
        EempropsConstIterator BeginEemprops() const { return _eep.begin(); }

        EempropsConstIterator EndEemprops() const { return _eep.end(); }

        EempropsIterator BeginEemprops() { return _eep.begin(); }

        EempropsIterator EndEemprops() { return _eep.end(); }
        
        EempropsConstIterator findEem(ChargeDistributionModel  eqdModel,
                                      const std::string       &name) const;

        EempropsIterator findEem(ChargeDistributionModel  eqdModel,
                                 const std::string       &name);
        
        void  setEpref(ChargeDistributionModel eqdModel, 
                       const std::string &epref);

        const char *getEpref(ChargeDistributionModel eqdModel) const;

        //! Spread from master to slave nodes
        void  broadcast(t_commrec *cr);

        EprefConstIterator epRefBegin() const { return _epr.begin(); }
        
        EprefConstIterator epRefEnd() const { return _epr.end(); }

    private:
        std::string                           _filename;
        std::vector<Ptype>                    _ptype;
        std::vector<Ffatype>                  _alexandria;
        std::vector<std::string>              _btype;
        std::string                           _alexandriaPolarUnit;
        std::string                           _alexandriaPolarRef;
        std::string                           _alexandriaForcefield;
        int                                   _nexcl;
        double                                _fudgeQQ, _fudgeLJ;
        std::string                           _gtVdwFunction, _gtCombinationRule;
        int                                   _gtVdwFtype, _gtCombRule;
        std::vector<GtBonds>                  _gtBonds;
        std::string                           _gtAngleFunction;
        unsigned int                          _gtAngleFtype;
        std::string                           _gtAngleUnit;
        std::vector<GtAngle>                  _gtAngle;
        std::vector<std::string>              _gtDihedralFunction;
        std::vector<unsigned int>             _gtDihedralFtype;
        std::vector<std::vector<GtDihedral> > _gtDihedral;
        std::vector<Miller>                   _miller;
        std::string                           _millerTauUnit, _millerAhpUnit;
        std::string                           _millerRef;
        std::vector<Bosque>                   _bosque;
        std::string                           _bosquePolarUnit;
        std::string                           _bosqueRef;
        std::vector<Symcharges>               _symcharges;
        std::vector<Eemprops>                 _eep;
        std::vector<Epref>                    _epr;

        void addBtype(const std::string &btype);

        gmx_bool strcasestrStart(std::string needle, std::string haystack);

        GtBondIterator searchBond(const std::string &atom1, 
                                  const std::string &atom2,
                                  double bondorder);

        static int gtbComp(const void *a, const void *b);

        static int gtdComp(const void *a, const void *b);

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
