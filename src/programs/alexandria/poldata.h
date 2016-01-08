/*! \internal \brief
 * Implements part of the alexandria program.
 * \author David van der Spoel <david.vanderspoel@icm.uu.se>
 */
#ifndef POLDATA_H
#define POLDATA_H

#include "gmxpre.h"

#include "gromacs/utility/smalloc.h"

#include "stringutil.h"

/* This source code file is part of the Alexandria project */

struct t_commrec;

/*! \brief
 * Enumerated type holding the charge distribution models used in PolData
 *
 * \inpublicapi
 * \ingroup module_alexandria
 */
enum ChargeDistributionModel {
    eqdAXp,
    eqdAXg,
    eqdAXs,
    eqdYang,
    eqdBultinck,
    eqdRappe,
    eqdNR
};

enum DihedralType {
    egdPDIHS,
    egdIDIHS,
    egdNR
};


namespace alexandria
{
class Ptype
{
    private:
        std::string type;
        std::string miller;
        std::string bosque;
        double      polarizability;
        double      sigPol;

    public:
        Ptype(std::string ptype, std::string miller, std::string bosque,
              double polarizability, double sigPol)
            :
              type(ptype),
              miller(miller),
              bosque(bosque),
              polarizability(polarizability),
              sigPol(sigPol)
        {
        }

        std::string getType()
        {
            return type;
        }

        std::string getMiller()
        {
            return miller;
        }

        std::string getBosque()
        {
            return bosque;
        }

        double getPolarizability()
        {
            return polarizability;
        }

        void setPolarizability(double polarizability)
        {
            this->polarizability = polarizability;
        }

        double getSigPol()
        {
            return sigPol;
        }

        void setSigPol(double sigPol)
        {
            this->sigPol = sigPol;
        }

};
typedef std::vector<Ptype>::iterator PtypeIterator;

class Ffatype
{
    private:
        std::string desc;
        std::string type;
        std::string ptype;
        std::string btype;
        std::string elem;
        std::string vdwparams;
        double      refEnthalpy;
    public:
        Ffatype(std::string desc,
                std::string type,
                std::string ptype,
                std::string btype,
                std::string elem,
                std::string vdwparams,
                double      refEnthalpy) :
            desc(desc),
            type(type),
            ptype(ptype),
            btype(btype),
            elem(elem),
            vdwparams(vdwparams),
            refEnthalpy(refEnthalpy)
        {
        }

        Ffatype()
        {
        }

        std::string getDesc() const { return desc; }

        std::string getType() const { return type; }

        std::string getPtype() const { return ptype; }

        std::string getBtype() const { return btype; }

        std::string getElem() const { return elem; }

        std::string getVdwparams() const { return vdwparams; }

        double getRefEnthalpy() const { return refEnthalpy; }
};
typedef std::vector<Ffatype>::iterator FfatypeIterator;


class Brule
{
    private:
        std::string              elem;
        std::string              rule;
        std::string              type;
        std::string              neighbors;
        std::string              geometry;
        int                      numbonds;
        int iAromatic;
        double                   valence;
        std::vector<std::string> nb;

    public:
        Brule(std::string              elem,
              std::string              rule,
              std::string              type,
              std::string              neighbors,
              std::string              geometry,
              int                      numbonds,
              int                      iAromatic,
              double                   valence,
              std::vector<std::string> nb) :
            elem(elem),
            rule(rule),
            type(type),
            neighbors(neighbors),
            geometry(geometry),
            numbonds(numbonds),
            iAromatic(iAromatic),
            valence(valence),
            nb(nb)
        {
        }

        std::string getElem()
        {
            return elem;
        }

        std::string getRule()
        {
            return rule;
        }

        std::string getType()
        {
            return type;
        }

        std::string getNeighbors()
        {
            return neighbors;
        }

        std::string getGeometry()
        {
            return geometry;
        }

        int getNumbonds()
        {
            return numbonds;
        }

        int getIAromatic()
        {
            return iAromatic;
        }

        double getValence()
        {
            return valence;
        }

        std::vector<std::string> getNb()
        {
            return nb;
        }

};
typedef std::vector<Brule>::iterator BruleIterator;

class GtBond
{
    private:
        std::string atom1;
        std::string atom2;
        std::string params;
        std::string elem1, elem2;
        double      length;
        double      sigma;
        double      bondorder;
        int         ntrain;
    public:
        GtBond(std::string atom1,
               std::string atom2,
               std::string params,
               std::string elem1,
               std::string elem2,
               double      length,
               double      sigma,
               double      bondorder,
               int         ntrain)
            :
              atom1(atom1),
              atom2(atom2),
              params(params),
              elem1(elem1),
              elem2(elem2),
              length(length),
              sigma(sigma),
              bondorder(bondorder),
              ntrain(ntrain)
        {
        }

        GtBond() {}

        std::string getAtom1() const { return atom1; }

        std::string getAtom2() const { return atom2; }

        std::string getParams() const { return params; }

        std::string getElem1() const { return elem1; }

        std::string getElem2() const { return elem2; }

        void setElem1(std::string elem) { elem1 = elem; }

        void setElem2(std::string elem) { elem2 = elem; }

        double getLength() { return length; }

        double getSigma() { return sigma; }

        double getBondorder() { return bondorder; }

        int getNtrain() { return ntrain; }

        void setAtom1(std::string atom) { atom1 = atom; }

        void setAtom2(std::string atom) { atom2 = atom; }

        void setBondorder(double bondorder) { this->bondorder = bondorder; }

        void setLength(double length) { this->length = length; }

        void setParams(std::string params) { this->params = params; }

        void setSigma(double sigma) { this->sigma = sigma; }

        void setNtrain(int ntrain) { this->ntrain = ntrain; }
};

typedef std::vector<GtBond>::iterator GtBondIterator;

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
        GtAngle(std::string atom1,
                std::string atom2,
                std::string atom3,
                std::string params,
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

        GtAngle() {}

        std::string getAtom1() const { return atom1_; }

        std::string getAtom2() const { return atom2_; }

        std::string getAtom3() const { return atom3_; }
        
        std::string getParams() const { return params_; }

        double getAngle() const { return angle_; }

        double getSigma() const { return sigma_; }

        int getNtrain() const { return ntrain_; }

        void setParams(std::string params) { params_ = params; }

        void setAngle(double angle) { angle_ = angle; }

        void setSigma(double sigma) { sigma_ = sigma; }

        void setNtrain(int ntrain) { ntrain_ = ntrain; }
};
typedef std::vector<GtAngle>::iterator GtAngleIterator;

class GtDihedral
{
    private:
        std::string atom1;
        std::string atom2;
        std::string atom3;
        std::string atom4;
        std::string params;
        double      dihedral;
        double      sigma;
        int         ntrain;
    public:
        GtDihedral(std::string atom1,
                   std::string atom2,
                   std::string atom3,
                   std::string atom4,
                   std::string params,
                   double      dihedral,
                   double      sigma,
                   int         ntrain)
            :
              atom1(atom1),
              atom2(atom2),
              atom3(atom3),
              atom4(atom4),
              params(params),
              dihedral(dihedral),
              sigma(sigma),
              ntrain(ntrain)
        {
        }

        GtDihedral()
        {
        }

        std::string getAtom1() const { return atom1; }

        std::string getAtom2() const { return atom2; }

        std::string getAtom3() const { return atom3; }

        std::string getAtom4() const { return atom4; }

        void setAtom1(std::string atom) { atom1 = atom; }

        void setAtom2(std::string atom) { atom2 = atom; }

        void setAtom3(std::string atom) { atom3 = atom; }

        void setAtom4(std::string atom) { atom4 = atom; }

        std::string getParams() { return params; }

        double getDihedral() { return dihedral; }

        double getSigma() { return sigma; }

        int getNtrain() { return ntrain; }
        
        bool compare(const GtDihedral &gtB) const;

        void setParams(std::string params) { this->params = params; }

        void setDihedral(double dihedral) { this->dihedral = dihedral; }

        void setSigma(double sigma) { this->sigma = sigma; }

        void setNtrain(int ntrain) { this->ntrain = ntrain; }
};
typedef std::vector<GtDihedral>::iterator DihedralIterator;

class Bosque
{
    private:
        std::string bosque;
        double      polarizability;
    public:
        Bosque(std::string bosque,
               double      polarizability)
            :
              bosque(bosque),
              polarizability(polarizability)
        {
        }

        std::string getBosque()
        {
            return bosque;
        }

        double getPolarizability()
        {
            return polarizability;
        }

};
typedef std::vector<Bosque>::iterator BosqueIterator;

class Miller
{
    private:
        std::string miller;
        int         atomnumber;
        double      tauAhc;
        double      alphaAhp;
    public:
        Miller(std::string miller,
               int         atomnumber,
               double      tauAhc,
               double      alphaAhp)
            :
              miller(miller),
              atomnumber(atomnumber),
              tauAhc(tauAhc),
              alphaAhp(alphaAhp)
        {
        }

        std::string getMiller()
        {
            return miller;
        }

        int getAtomnumber()
        {
            return atomnumber;
        }

        double getTauAhc()
        {
            return tauAhc;
        }

        double getAlphaAhp()
        {
            return alphaAhp;
        }
};
typedef std::vector<Miller>::iterator MillerIterator;

class Symcharges
{
    private:
        std::string central;
        std::string attached;
        int         numattach;
    public:
        Symcharges(std::string central,
                   std::string attached,
                   int         numattach) :
            central(central),
            attached(attached),
            numattach(numattach)
        {
        }

        std::string getCentral()
        {
            return central;
        }

        std::string getAttached()
        {
            return attached;
        }

        int getNumattach()
        {
            return numattach;
        }
};

typedef std::vector<Symcharges>::iterator SymchargesIterator;

class Epref
{
    private:
        ChargeDistributionModel eqdModel;
        std::string             epref;
    public:
        Epref(ChargeDistributionModel eqdModel,
              std::string             epref)
            :
              eqdModel(eqdModel),
              epref(epref)
        {
        }


        ChargeDistributionModel getEqdModel()
        {
            return eqdModel;
        }

        std::string getEpref()
        {
            return epref;
        }

        void setEpref(std::string epref)
        {
            this->epref = epref;
        }

};
typedef std::vector<Epref>::iterator EprefIterator;


class EempropsData
{
    private:
        std::vector<double> q;
        std::vector<double> zeta;
        std::vector<int>    row;
    public:

        EempropsData(int size)
            :
              q(size),
              zeta(size),
              row(size)
        {
        }

        EempropsData()
        {
        }


        std::vector<double> getQ()
        {
            return q;
        }

        std::vector<double> getZeta()
        {
            return zeta;
        }

        std::vector<int> getRow()
        {
            return row;
        }


        void setZeta(double zeta, int index)
        {
            this->zeta[index] = zeta;
        }

        void setRow(double row, int index)
        {
            this->row[index] = row;
        }


        void setQ(double q, int index)
        {
            this->q[index] = q;
        }
};

typedef std::vector<EempropsData>::iterator EempropsDataIterator;

//#define EEMBUFSIZE 256
#define MAXZETA    12
class Eemprops
{
    private:
        ChargeDistributionModel eqdModel;
        int                     nzeta;
        std::string             name;
        std::string             zetastr;
        std::string             qstr;
        std::string             rowstr;
        double                  J0;
        double                  chi0;
        EempropsData            data;

    public:
        Eemprops(ChargeDistributionModel eqdModel,
                 int                     nzeta,
                 std::string             name,
                 std::string             zetastr,
                 std::string             qstr,
                 std::string             rowstr,
                 double                  J0,
                 double                  chi0)
            :
              eqdModel(eqdModel),
              nzeta(nzeta),
              name(name),
              zetastr(zetastr),
              qstr(qstr),
              rowstr(rowstr),
              J0(J0),
              chi0(chi0),
              data(MAXZETA)
        {
        }

        Eemprops() :
            data(MAXZETA)
        {
        }

        ChargeDistributionModel getEqdModel()
        {
            return eqdModel;
        }

        int getNzeta()
        {
            return nzeta;
        }

        std::string getName()
        {
            return name;
        }

        std::string getZetastr()
        {
            return zetastr;
        }

        std::string getQstr()
        {
            return qstr;
        }

        std::string getRowstr()
        {
            return rowstr;
        }

        double getJ0()
        {
            return J0;
        }

        double getChi0()
        {
            return chi0;
        }

        void setEqdModel(ChargeDistributionModel eqdModel)
        {
            this->eqdModel = eqdModel;
        }

        void setNzeta(int nzeta)
        {
            this->nzeta = nzeta;
        }

        void setName(std::string name)
        {
            this->name = name;
        }

        void setZetastr(std::string zetastr)
        {
            this->zetastr = zetastr;
        }

        void setQstr(std::string qstr)
        {
            this->qstr = qstr;
        }

        void setRowstr(std::string rowstr)
        {
            this->rowstr = rowstr;
        }

        void setJ0(double J0)
        {
            this->J0 = J0;
        }

        void setChi0(double chi0)
        {
            this->chi0 = chi0;
        }

        double getZeta(int index)
        {
            return data.getZeta()[index];
        }

        double getQ(int index)
        {
            return data.getQ()[index];
        }

        int getRow(int index)
        {
            return data.getRow()[index];
        }

        std::vector<double> getAllZeta()
        {
            return data.getZeta();
        }

        std::vector<double> getAllQ()
        {
            return data.getQ();
        }

        std::vector<int> getAllRow()
        {
            return data.getRow();
        }

        void setZeta(int index, double zeta)
        {
            data.setZeta(zeta, index);
        }

        void setQ(int index, double q)
        {
            data.setQ(q, index);
        }

        void setRow(int index, int row)
        {
            data.setRow(row, index);
        }
};
typedef std::vector<Eemprops>::iterator EempropsIterator;


class Poldata
{
    public:

        Poldata(); //Constructor
        ~Poldata(){}


        void  setFilename(std::string fn2);

        void  addBondingRule(std::string gtBrule, std::string atype,
                             std::string geometry, int numbonds,
                             double valence, int iAromatic,
                             std::string neighbors);

        BruleIterator getBruleBegin()
        {
            return _brule.begin();
        }

        BruleIterator getBruleEnd()
        {
            return _brule.end();
        }

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
                       double            ref_enthalpy);

        void  setPtypePolarizability(  const std::string ptype,
                                       double polarizability, double sigPol);



        void setPolarUnit( const std::string polarUnit)
        {
            _alexandriaPolarUnit.assign(polarUnit);

        }

        void setPolarRef( const std::string polarRef)
        {
            _alexandriaPolarRef.assign(polarRef);
        }

        std::string getForceField()
        {
            return _alexandriaForcefield;
        }

        void setForceField( const std::string forcefield)
        {
            _alexandriaForcefield.assign(forcefield);
        }

        void setLengthUnit( const std::string length_unit)
        {
            _gtLengthUnit.assign(length_unit);
        }

        void  setVdwFunction(  const std::string func);

        int  getVdwFtype( );

        void setNexcl( int nexcl)
        {
            _nexcl = nexcl;
        }

        int getNexcl( );


        void setFudgeQQ( double fudgeQQ)
        {
            _fudgeQQ = fudgeQQ;
        }

        int getNatypes()
        {
            return _alexandria.size();
        }

        int getNptypes()
        {
            return _ptype.size();
        }

        unsigned int getNgtBond()
        {
            return _gtBond.size();
        }

        int getNgtAngle()
        {
            return _gtAngle.size();
        }

        int getNgtDihedral( int egd)
        {
            return _gtDihedral[egd].size();
        }

        double  getFudgeQQ( );

        void setFudgeLJ(double fudgeLJ)
        {
            _fudgeLJ = fudgeLJ;
        }

        double  getFudgeLJ( );

        int getAtypeRefEnthalpy( const std::string atype,
                                 double           *Href);

        void  setCombinationRule(  std::string func);

        std::string  getCombinationRule( );

        int  getCombRule( );

        std::string  getLengthUnit( );

        /* Return array of atomtypes compatible with the bonded neighbors.
           The array should be freed, but not the contents of the elements.
         */
        std::string * getBondingRules(  std::string elem,
                                        int nbond, std::string neighbors[],
                                        const std::string geometry,
                                        int iAromatic);

        std::string  getGeometry(  std::string gtBrule);

        std::string  getDesc(  std::string atype);

        /* Get the charge from the gentop.dat file */
        std::string  getCharge(  std::string atype);

        FfatypeIterator getAtypeBegin()
        {
            return _alexandria.begin();
        }

        FfatypeIterator getAtypeEnd()
        {
            return _alexandria.end();
        }
        
        FfatypeIterator searchType(const std::string &type);
        
        FfatypeIterator searchBtype(const std::string &btype);
        
        FfatypeIterator searchPtype(const std::string &ptype);
        
        PtypeIterator getPtypeBegin()
        {
            return _ptype.begin();
        }

        PtypeIterator getPtypeEnd()
        {
            return _ptype.end();
        }

        std::string  atypeToPtype(  const std::string atype);

        std::string  atypeToBtype(  const std::string atype);

        /* Return 1 if OK, 0 if not found */
        int  searchAtype(std::string   key,
                         Ffatype     * type);

        /* Return 1 if OK, 0 if not found */
        int  getPtypePol(  const std::string ptype,
                           double *polarizability, double *sigPol);
        int  getAtypePol(  const std::string atype,
                           double *polarizability, double *sigPol);

        int getAtypeRefEnthalpy(Poldata * pd, const std::string atype,
                                double *Href);

        /* Return 1 if OK, 0 if not found */
        int  bondingRuleValence(  std::string gtBrule, double *valence);

        void  addMiller(std::string   miller,
                        int           atomnumber,
                        double        tauAhc,
                        double        alphaAhp);

        /* Return 1 if "miller" was found */
        int  getMillerPol(std::string   miller,
                          int          *atomnumber,
                          double       *tauAhc,
                          double       *alphaAhp);



        MillerIterator getMillerBegin()
        {
            return _miller.begin();
        }

        MillerIterator getMillerEnd()
        {
            return _miller.end();
        }



        void  setMillerUnits(  std::string tauUnit,
                               std::string ahpUnit);

        void  getMillerUnits(  std::string *tauUnit,
                               std::string *ahpUnit);

        /* Returns miller name or NULL if not found */
        std::string  ptypeToMiller(  const std::string ptype);

        void  addBosque(std::string   bosque,
                        double        polarizability);

        BosqueIterator getBosqueBegin()
        {
            return _bosque.begin();
        }

        BosqueIterator getBosqueEnd()
        {
            return _bosque.end();
        }

        void setBosqueUnit( std::string polarUnit)
        {
            _bosquePolarUnit.assign(polarUnit);
        }

        std::string getBosqueUnit()
        {
            return _bosquePolarUnit;
        }

        /* Returns bosque name or NULL if not found */
        std::string  ptypeToBosque(  const std::string ptype);

        int  getBosquePol(  std::string bosque, double *polarizability);

        /* Return 1 on success or 0 otherwise */
        int  addBond(  std::string atom1, std::string atom2,
                       double length, double sigma, int ntrain,
                       double bondorder, std::string params);

        int  setBondParams(  std::string atom1, std::string atom2,
                             double length, double sigma, int ntrain,
                             double bondorder, std::string params);

        GtBondIterator getBondBegin()
        {
            return _gtBond.begin();
        }

        GtBondIterator getBondEnd()
        {
            return _gtBond.end();
        }

        void  setBondFunction(  std::string fn);

        std::string getBondFunction()
        {
            return _gtBondFunction;
        }

        /* Return bond-index 1-N or 0 if not found */
        int  searchBond(  std::string atom1, std::string atom2,
                          double *length, double *sigma, int *ntrain,
                          double *bondorder, std::string *params);

        /* Returns 1 if there is a bond, 0 if not. Toler is absolute in length-units. */
        int  elemIsBond(  std::string elem1, std::string elem2,
                          double distance, double toler);

        /* Return maximal valence for a give element */
        double  elemGetMaxValence(  std::string elem);

        /* Return NULL-terminated array of potential bondorders */
        double * elemGetBondorders(  std::string elem1, std::string elem2,
                                     double distance, double toler);
        /* Returns the bondorder. Toler is absolute in length-units. */
        double  atypeBondorder(  std::string atype1, std::string atype2,
                                 double distance, double toler);

        void  setAngleFunction(  std::string fn);

        std::string getAngleFunction()
        {
            return _gtAngleFunction;
        }

        int getBondFtype()
        {
            return _gtBondFtype;
        }

        int getAngleFtype()
        {
            return _gtAngleFtype;
        }

        int getDihedralFtype( int egd)
        {
            return _gtDihedralFtype[egd];
        }

        std::string getVdwFunction()
        {
            return _gtVdwFunction;
        }

        std::string getPolarUnit()
        {
            return _alexandriaPolarUnit;
        }

        std::string getPolarRef()
        {
            return _alexandriaPolarRef;
        }

        /* Return 1 on success, 0 otherwise */
        int  addAngle(
            std::string atom1, std::string atom2,
            std::string atom3, double angle, double sigma,
            int ntrain, std::string params);

        int  setAngleParams(std::string atom1,
                            std::string atom2,
                            std::string atom3, 
                            double angle, 
                            double sigma, 
                            int ntrain, 
                            std::string params);


        GtAngleIterator getAngleBegin()
        {
            return _gtAngle.begin();
        }

        GtAngleIterator getAngleEnd()
        {
            return _gtAngle.end();
        }

        /* Return angle-index 1-N or 0 if not found */
        int  searchAngle(  std::string atom1, std::string atom2,
                           std::string atom3, double *angle, double *sigma,
                           int *ntrain, std::string *params);

        void setAngleUnit( std::string angleUnit)
        {
            _gtAngleUnit.assign(angleUnit);
        }

        std::string getAngleUnit()
        {
            return _gtAngleUnit;
        }

        void  setDihedralFunction(  int egd, std::string fn);

        std::string getDihedralFunction( int egd)
        {
            return _gtDihedralFunction[egd];
        }

        /* Return 1 on success or 0 otherwise */
        int  addDihedral(  int egd, std::string atom1, std::string atom2,
                           std::string atom3, std::string atom4,
                           double dihedral, double sigma,
                           int ntrain, std::string params);

        int  setDihedralParams(  int egd,
                                 std::string atom1, std::string atom2,
                                 std::string atom3, std::string atom4,
                                 double angle, double sigma,
                                 int ntrain, std::string params);


        DihedralIterator getDihedralBegin(int egd)
        {
            return _gtDihedral[egd].begin();
        }

        DihedralIterator getDihedralEnd(int egd)
        {
            return _gtDihedral[egd].end();
        }

        /* Return dihedral-index 1-N or 0 if not found */
        DihedralIterator searchDihedral(int egd,
                                        const std::string &atom1, 
                                        const std::string &atom2,
                                        const std::string &atom3, 
                                        const std::string &atom4);

        int searchDihedral(int egd, 
                           std::string atom1, std::string atom2,
                           std::string atom3, std::string atom4,
                           double *dihedral, double *sigma,
                           int *ntrain, std::string *params);

        void  setDihedralUnit(  int         egd,
                                std::string dihedralUnit);

        std::string  getDihedralUnit(  int egd);

        void  addSymcharges(  std::string central,
                              std::string attached, int numattach);

        int  getSymcharges(  std::string *central,
                             std::string *attached, int *numattach);

        SymchargesIterator getSymchargesBegin()
        {
            return _symcharges.begin();
        }

        SymchargesIterator getSymchargesEnd()
        {
            return _symcharges.end();
        }

        int  searchSymcharges(  std::string central,
                                std::string attached, int numattach);

        static ChargeDistributionModel name2eemtype(const std::string name);

        static  std::string getEemtypeName(ChargeDistributionModel eem);

        std::string  getEemref(  ChargeDistributionModel eqdModel);

        int  getNumprops(  ChargeDistributionModel eqdModel);

        int  havePolSupport(  const std::string atype);

        int  haveEemSupport(  ChargeDistributionModel eqdModel,
                              const std::string       name,
                              gmx_bool                bAllowZeroParameters);

        double  getJ00(  ChargeDistributionModel eqdModel, const std::string name);

        int  getNzeta(  ChargeDistributionModel eqdModel, const std::string name);

        double  getZeta(  ChargeDistributionModel eqdModel, const std::string name, int zz);

        std::string  getQstr(  ChargeDistributionModel eqdModel, std::string name);

        std::string  getRowstr(  ChargeDistributionModel eqdModel, std::string name);

        double  getQ(  ChargeDistributionModel eqdModel, const std::string name, int zz);

        int  getRow(  ChargeDistributionModel eqdModel, const std::string name, int zz);

        double  getChi0(  ChargeDistributionModel eqdModel, const std::string name);

        std::string  getOpts(  ChargeDistributionModel eqdModel, std::string name);

        void  setEemprops(ChargeDistributionModel eqdModel, const std::string name,
                          double J0, double chi0,
                          const std::string zeta, const  std::string q, const std::string row);

        EempropsIterator getEempropsBegin()
        {
            return _eep.begin();
        }

        EempropsIterator getEempropsEnd()
        {
            return _eep.end();
        }

        void  setEpref(  ChargeDistributionModel eqdModel, std::string epref);

        std::string getEpref(ChargeDistributionModel eqdModel);

        EprefIterator getEprefBegin()
        {
            return _epr.begin();
        }

        EprefIterator getEprefEnd()
        {
            return _epr.end();
        }

        void  commEemprops(  t_commrec *cr);

        void  commForceParameters( t_commrec *cr);

    private:
        std::string                           _filename;
        std::vector<Ptype>                    _ptype;
        std::vector<Ffatype>                  _alexandria;
        std::vector<std::string>              _btype;
        std::vector<Brule>                    _brule;
        std::string                           _alexandriaPolarUnit;
        std::string                           _alexandriaPolarRef;
        std::string                           _alexandriaForcefield;
        int                                   _nexcl;
        double                                _fudgeQQ, _fudgeLJ;
        std::string                           _gtVdwFunction, _gtCombinationRule;
        int                                   _gtVdwFtype, _gtCombRule;
        std::string                           _gtBondFunction;
        unsigned int                          _gtBondFtype;
        std::string                           _gtLengthUnit;
        std::vector<GtBond>                   _gtBond;
        std::string                           _gtAngleFunction;
        unsigned int                          _gtAngleFtype;
        std::string                           _gtAngleUnit;
        std::vector<GtAngle>                  _gtAngle;
        std::vector<std::string>              _gtDihedralFunction;
        std::vector<unsigned int>             _gtDihedralFtype;
        std::vector<std::vector<GtDihedral> > _gtDihedral;
        std::vector<Miller>                   _miller;
        std::string                           _millerTauUnit, _millerAhpUnit;
        std::vector<Bosque>                   _bosque;
        std::string                           _bosquePolarUnit;
        std::vector<Symcharges>               _symcharges;
        std::vector<Eemprops>                 _eep;
        std::vector<Epref>                    _epr;

        void addBtype(const std::string &btype);

        gmx_bool strcasestrStart(std::string needle, std::string haystack);

        int countNeighbors(Brule *brule, int nbond, std::string nbhybrid[], int *score);

        GtBond *searchBond(  std::string atom1, std::string atom2,
                             double bondorder);

        int searchBondtype(  std::string atom);

        Eemprops *getEep(ChargeDistributionModel       eqdModel,
                         const std::string             name);

        static int gtbComp(const void *a, const void *b);

        static int gtdComp(const void *a, const void *b);

        template<class Type>
        int indexOfPointInVector(Type * pointer, std::vector<Type> vector)
        {
            return (pointer - &(vector[0]));
        }
};
}
#endif
