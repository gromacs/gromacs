/*
 * This source file is part of the Alexandria project.
 *
 * Copyright (C) 2014 David van der Spoel and Paul J. van Maaren
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
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.
 */
/*! \internal \brief
 * Implements part of the alexandria program.
 * \author David van der Spoel <david.vanderspoel@icm.uu.se>
 */
#ifndef POLDATA_H
#define POLDATA_H

#include "gromacs/legacyheaders/types/commrec.h"
#include "gromacs/legacyheaders/typedefs.h"
#include "gromacs/legacyheaders/network.h"
#include "gmxpre.h"
#include "gromacs/legacyheaders/typedefs.h"
#include "gromacs/legacyheaders/network.h"
#include "gromacs/utility/smalloc.h"
#include "stringutil.h"
/* This source code file is part of the Alexandria project */

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

class Ptype {
public:
  std::string type;
  std::string miller; 
  std::string bosque;
  double polarizability;
  double sig_pol;
}; 

class Ffatype {
public:
  std::string desc;
  std::string type;
  std::string ptype;
  std::string btype;
  std::string elem; 
  std::string vdwparams;
  double  ref_enthalpy;
};

class Brule {
public:
  std::string elem;
  std::string rule;
  std::string type;
  std::string neighbors;
  std::string geometry;
  int                      numbonds;
  int iAromatic;
  double                   valence;
  std::vector<std::string> nb;
}; 

class  GtBond {
public:
  std::string atom1;
  std::string atom2;
  std::string params;
  std::string    elem1, elem2;
  double  length;
  double sigma;
  double bondorder;
  int     ntrain;
};

class GtAngle {
public:
  std::string   atom1; 
  std::string atom2;
  std::string atom3; 
  std::string params;
  double  angle;
  double  sigma;
  int     ntrain;
};

class GtDihedral {
public:
  std::string atom1; 
  std::string atom2;
  std::string atom3;
  std::string atom4;
  std::string params;
  double  dihedral;
  double  sigma;
  int     ntrain;
};

class Bosque {
public:
  std::string bosque;
  double  polarizability;
};

class Miller {
public:
  std::string miller;
  int     atomnumber;
  double  tau_ahc;
  double alpha_ahp;
};

class Symcharges {
public:
  std::string central;
  std::string attached;
  int   numattach;
};

class Epref {
public:
  ChargeDistributionModel eqd_model;
  std::string epref;
};

#define EEMBUFSIZE 256
#define MAXZETA    12
class Eemprops {
 public:
  ChargeDistributionModel eqd_model;
  int nzeta;
  int row[MAXZETA];
  std::string name;
  std::string zetastr;
  std::string qstr;
  std::string rowstr;
  double J0;
  double chi0;
  double q[MAXZETA];
  double zeta[MAXZETA];
};


namespace alexandria
{
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

    int  getBondingRule(std::string *gtBrule, std::string *atype,
			std::string *geometry, int *numbonds,
			double *valence, int *iAromatic,
			std::string *neighbors);

    void  addPtype(const std::string ptype,
		   const std::string miller,
		   const std::string bosque,
		   double        polarizability,
		   double        sigPol);

    void  addAtype(const std::string elem,
		   const std::string desc,
		   const std::string atype,
		   const std::string ptype,
		   const std::string btype,
		   const std::string vdwparams,
		   double      ref_enthalpy);

    void  setPtypePolarizability(  const std::string ptype,
				   double polarizability, double sigPol);

   

    void setPolarUnit( const std::string polarUnit){
      _alexandriaPolarUnit.assign(polarUnit);

    }

    void setPolarRef( const std::string polarRef){
      _alexandriaPolarRef.assign(polarRef);
    }

     std::string getForceField(){
      return _alexandriaForcefield;
    }

    void setForceField( const std::string forcefield){
      _alexandriaForcefield.assign(forcefield);
    }

    void setLengthUnit( const std::string length_unit){
      _gtLengthUnit.assign(length_unit);
    }

    void  setVdwFunction(  const std::string func);

    int  getVdwFtype( );

    void setNexcl( int nexcl){
      _nexcl = nexcl;
    }

    int  getNexcl( );


    void setFudgeQQ( double fudgeQQ){
      _fudgeQQ = fudgeQQ;
    }

    int getNatypes(){
      return _alexandria.size();
    }

    int getNptypes(){
      return _ptype.size();
    }

    unsigned int getNgtBond(){
      return _gtBond.size();
    }

    int getNgtAngle(){
      return _gtAngle.size();
    }

    int getNgtDihedral( int egd){
      return _gtDihedral[egd].size();
    }

    double  getFudgeQQ( );

    void setFudgeLJ(double fudgeLJ){
      _fudgeLJ = fudgeLJ;
    }

    double  getFudgeLJ( );

    int getAtypeRefEnthalpy( const std::string atype,
			     double     *Href);

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

    /* Returns 1 if OK, 0 if last */
    int  getAtype(
		  std::string  *elem,
		  std::string  *desc,
		  std::string  *atype,
		  std::string  *ptype,
		  std::string  *btype,
		  std::string  *vdwparams,
		  double       *refEnthalpy);

    int  getPtype(
		  std::string  *ptype,
		  std::string  *miller,
		  std::string  *bosque,
		  double       *polarizability,
		  double       *sigPol);

     std::string  atypeToPtype(  const std::string atype);

     std::string  atypeToBtype(  const std::string atype);

    /* Return 1 if OK, 0 if not found */
    int  searchAtype(
		     std::string   key,
		     std::string  *elem,
		     std::string  *desc,
		     std::string  *atype,
		     std::string  *ptype,
		     std::string  *btype,
		     std::string  *vdwparams);

    /* Return 1 if OK, 0 if not found */
    int  getPtypePol(  const std::string ptype,
		       double *polarizability, double *sigPol);
    int  getAtypePol(  const std::string atype,
		       double *polarizability, double *sigPol);

    int getAtypeRefEnthalpy(Poldata * pd, const std::string atype,
			    double *Href);

    /* Return 1 if OK, 0 if not found */
    int  bondingRuleValence(  std::string gtBrule, double *valence);

    void  addMiller(
		    std::string   miller,
		    int           atomnumber,
		    double        tauAhc,
		    double        alphaAhp);

    /* Return 1 if "miller" was found */
    int  getMillerPol(
		      std::string   miller,
		      int          *atomnumber,
		      double       *tauAhc,
		      double       *alphaAhp);

    int  getMiller(
		   std::string  *miller,
		   int          *atomnumber,
		   double       *tauAhc,
		   double       *alphaAhp);

    void  setMillerUnits(  std::string tauUnit,
			   std::string ahpUnit);

    void  getMillerUnits(  std::string *tauUnit,
			   std::string *ahpUnit);

    /* Returns miller name or NULL if not found */
     std::string  ptypeToMiller(  const std::string ptype);

    void  addBosque(
		    std::string   bosque,
		    double        polarizability);

    int  getBosque(
		   std::string  *bosque,
		   double       *polarizability);

    void setBosqueUnit( std::string polarUnit){
      _bosquePolarUnit.assign(polarUnit);
    }

     std::string getBosqueUnit(){
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

    /* Return bond-index 1-N or 0 if not found */
    int  getBond(  std::string *atom1, std::string *atom2,
		   double *length, double *sigma, int *ntrain,
		   double *bondorder, std::string *params);

    void  setBondFunction(  std::string fn);

     std::string getBondFunction(){
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

     std::string getAngleFunction() {
      return _gtAngleFunction;
    }

    int getBondFtype() {
      return _gtBondFtype;
    }

    int getAngleFtype() {
      return _gtAngleFtype;
    }

    int getDihedralFtype( int egd){
      return _gtDihedralFtype[egd];
    }

     std::string getVdwFunction(){
      return _gtVdwFunction;
    }

     std::string getPolarUnit(){
      return _alexandriaPolarUnit;
    }

     std::string getPolarRef(){
      return _alexandriaPolarRef;
    }

    /* Return 1 on success, 0 otherwise */
    int  addAngle(
		  std::string atom1, std::string atom2,
		  std::string atom3, double angle, double sigma,
		  int ntrain, std::string params);

    int  setAngleParams(  std::string atom1, std::string atom2,
			  std::string atom3,
			  double angle, double sigma, int ntrain, std::string params);

    /* Return angle-index 1-N or 0 if not found */
    int  getAngle(  std::string *atom1, std::string *atom2,
		    std::string *atom3, double *angle, double *sigma,
		    int *ntrain, std::string *params);

    /* Return angle-index 1-N or 0 if not found */
    int  searchAngle(  std::string atom1, std::string atom2,
		       std::string atom3, double *angle, double *sigma,
		       int *ntrain, std::string *params);

    void setAngleUnit( std::string angleUnit){
      _gtAngleUnit.assign(angleUnit);
    }

      std::string getAngleUnit(){
      return _gtAngleUnit;
    }


    void  setDihedralFunction(  int egd, std::string fn);

      std::string getDihedralFunction( int egd){
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

    /* Return dihedral-index 1-N or 0 if not found */
    int  getDihedral(  int egd,
		       std::string *atom1, std::string *atom2,
		       std::string *atom3, std::string *atom4,
		       double *dihedral, double *sigma,
		       int *ntrain, std::string *params);

    /* Return dihedral-index 1-N or 0 if not found */
    int  searchDihedral(  int egd,
			  std::string atom1, std::string atom2,
			  std::string atom3, std::string atom4,
			  double *dihedral, double *sigma,
			  int *ntrain, std::string *params);

    void  setDihedralUnit(  int   egd,
			    std::string dihedralUnit);

     std::string  getDihedralUnit(  int egd);

    void  addSymcharges(  std::string central,
			  std::string attached, int numattach);

    int  getSymcharges(  std::string *central,
			 std::string *attached, int *numattach);

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

    double  getJ00(  ChargeDistributionModel eqdModel,const std::string name);

    int  getNzeta(  ChargeDistributionModel eqdModel,const std::string name);

    double  getZeta(  ChargeDistributionModel eqdModel,const std::string name, int zz);

     std::string  getQstr(  ChargeDistributionModel eqdModel, std::string name);

     std::string  getRowstr(  ChargeDistributionModel eqdModel, std::string name);

    double  getQ(  ChargeDistributionModel eqdModel,const std::string name, int zz);

    int  getRow(  ChargeDistributionModel eqdModel,const std::string name, int zz);

    double  getChi0(  ChargeDistributionModel eqdModel,const std::string name);

     std::string  getOpts(  ChargeDistributionModel eqdModel, std::string name);

    void  setEemprops(ChargeDistributionModel eqdModel,const std::string name,
		      double J0, double chi0,
		      const std::string zeta,const  std::string q, const std::string row);

    int  getEemprops(
		     ChargeDistributionModel *eqdModel, std::string *name,
		     double *J0, double *chi0,
		     std::string *zeta, std::string *q, std::string *row);

    void  setEpref(  ChargeDistributionModel eqdModel, std::string epref);

     std::string  getEpref(ChargeDistributionModel eqdModel);

    int listEpref(  ChargeDistributionModel *eqdModel, std::string *epref);

    void  commEemprops(  t_commrec *cr);

    void  commForceParameters( t_commrec *cr);


  private:
    std::string          _filename;
    unsigned int            _nptypeC;
    std::vector<Ptype>       _ptype;
    unsigned int            _nalexandriaC;
    std::vector<Ffatype>  _alexandria;
    std::vector<std::string>         _btype;
    unsigned int            _nbruleC;
    std::vector<Brule> _brule;
    std::string          _alexandriaPolarUnit;
    std::string          _alexandriaPolarRef;
    std::string          _alexandriaForcefield;
    int            _nexcl;
    double         _fudgeQQ, _fudgeLJ;
    std::string          _gtVdwFunction, _gtCombinationRule;
    int            _gtVdwFtype, _gtCombRule;
    std::string          _gtBondFunction;
    unsigned int            _ngtBondC, _gtBondFtype;
    std::string          _gtLengthUnit;
    std::vector<GtBond>  _gtBond;
    std::string          _gtAngleFunction;
    unsigned int            _ngtAngleC, _gtAngleFtype;
    std::string          _gtAngleUnit;
    std::vector<GtAngle>     _gtAngle;
    std::vector<std::string>          _gtDihedralFunction;
    std::vector<unsigned int>             _ngtDihedralC, _gtDihedralFtype;
    std::vector<std::vector<GtDihedral> > _gtDihedral;
    unsigned int            _nmillerC;
    std::vector<Miller>   _miller;
    std::string          _millerTauUnit, _millerAhpUnit;
    unsigned int             _nbosqueC;
    std::vector<Bosque>     _bosque;
    std::string          _bosquePolarUnit;
    unsigned int              _nsymchargesC;
    std::vector<Symcharges>  _symcharges;
    unsigned int            _nepC;
    std::vector<Eemprops>  _eep;
    unsigned int             _nerC;
    std::vector<Epref>   _epr;

    void addBtype(const std::string btype);

    gmx_bool strcasestrStart(std::string needle, std::string haystack);

    int countNeighbors(Brule *brule, int nbond, std::string nbhybrid[], int *score);

    GtBond *searchBond(  std::string atom1, std::string atom2,
			 double bondorder);

    int searchBondtype(  std::string atom);

    Eemprops *getEep(ChargeDistributionModel eqdModel,
		     const std::string             name);

    GtDihedral *searchDihedral( int egd, std::string atom1, std::string atom2,
				std::string atom3, std::string atom4);

    static int gtbComp(const void *a, const void *b);

    static int gtdComp(const void *a, const void *b);

    template<class Type>
      int indexOfPointInVector(Type * pointer, std::vector<Type> vector){
      return (pointer - &(vector[0]));
    }
  };
}
#endif
