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

typedef struct {
    char   *type, *miller, *bosque;
    double  polarizability, sig_pol;
} t_ptype;

typedef struct {
    char   *desc, *type, *ptype, *btype, *elem, *vdwparams;
    double  ref_enthalpy;
} t_ffatype;

typedef struct {
    char                    *elem, *rule, *type, *neighbors, *geometry;
    int                      numbonds, iAromatic;
    double                   valence;
    std::vector<std::string> nb;
} t_brule;

typedef struct {
    char   *atom1, *atom2, *params;
    char    elem1[4], elem2[4];
    double  length, sigma, bondorder;
    int     ntrain;
} t_gt_bond;

typedef struct {
    char   *atom1, *atom2, *atom3, *params;
    double  angle, sigma;
    int     ntrain;
} t_gt_angle;

typedef struct {
    char   *atom1, *atom2, *atom3, *atom4, *params;
    double  dihedral, sigma;
    int     ntrain;
} t_gt_dihedral;

typedef struct {
    char   *bosque;
    double  polarizability;
} t_bosque;

typedef struct {
    char   *miller;
    int     atomnumber;
    double  tau_ahc, alpha_ahp;
} t_miller;

typedef struct {
    char *central;
    char *attached;
    int   numattach;
} t_symcharges;

typedef struct {
    ChargeDistributionModel eqd_model;
    char                   *epref;
} t_epref;

#define EEMBUFSIZE 256
#define MAXZETA    12
typedef struct {
    ChargeDistributionModel eqd_model;
    int                     nzeta, row[MAXZETA];
    char                    name[EEMBUFSIZE], zetastr[EEMBUFSIZE], qstr[EEMBUFSIZE], rowstr[EEMBUFSIZE];
    double                  J0, chi0, q[MAXZETA], zeta[MAXZETA];
} t_eemprops;


namespace alexandria
{
  class Poldata
  {
  public:

    Poldata(); //Constructor
    ~Poldata(){}


    void  setFilename(char *fn2);


    int  getNatypes( );
    int  getNptypes( );
    int  getNgtBond( );
    int  getNgtAngle( );
    int  getNgtDihedral(  int egd);

    void  addBondingRule(
			   char *gtBrule, char *atype,
			   char *geometry, int numbonds,
			   double valence, int iAromatic,
			   char *neighbors);

    int  getBondingRule(
			  char **gtBrule, char **atype,
			  char **geometry, int *numbonds,
			  double *valence, int *iAromatic,
			  char **neighbors);

    void  addPtype(
		    const char   *ptype,
		    const char   *miller,
		    const char   *bosque,
		    double        polarizability,
		    double        sigPol);

    void  addAtype(  const char *elem,
		      const char *desc,
		      const char *atype,
		      const char *ptype,
		      const char *btype,
		      const char *vdwparams,
		      double      ref_enthalpy);

    void  setPtypePolarizability(  const char *ptype,
				     double polarizability, double sigPol);

    void  setForceField(  const char *forcefield);

    void  setPolarUnit(  const char *polarUnit);

    void  setPolarRef(  const char *polarRef);

    void  setLengthUnit(  const char *lengthUnit);

    void  setVdwFunction(  const char *func);

    char * getVdwFunction( );

    int  getVdwFtype( );

    void  setNexcl(  int nexcl);

    int  getNexcl( );

    void  setFudgeQQ(  double fudgeQQ);

    double  getFudgeQQ( );

    void  setFudgeLJ(  double fudgeLJ);

    double  getFudgeLJ( );

    int getAtypeRefEnthalpy( const char *atype,
				double     *Href);

    void  setCombinationRule(  char *func);

    char * getCombinationRule( );

    int  getCombRule( );

    char * getPolarUnit( );

    char * getPolarRef( );

    char * getForceField( );

    char * getLengthUnit( );


    /* Return array of atomtypes compatible with the bonded neighbors.
       The array should be freed, but not the contents of the elements.
    */
    char ** getBondingRules(  char *elem,
				int nbond, char *neighbors[],
				const char *geometry,
				int iAromatic);

    char * getGeometry(  char *gtBrule);

    char * getDesc(  char *atype);

    /* Get the charge from the gentop.dat file */
    char * getCharge(  char *atype);

    /* Returns 1 if OK, 0 if last */
    int  getAtype(
		   char        **elem,
		   char        **desc,
		   char        **atype,
		   char        **ptype,
		   char        **btype,
		   char        **vdwparams,
		   double       *refEnthalpy);

    int  getPtype(
		   char        **ptype,
		   char        **miller,
		   char        **bosque,
		   double       *polarizability,
		   double       *sigPol);

    const char * atypeToPtype(  const char *atype);

    const char * atypeToBtype(  const char *atype);

    /* Return 1 if OK, 0 if not found */
    int  searchAtype(
		      char         *key,
		      char        **elem,
		      char        **desc,
		      char        **atype,
		      char        **ptype,
		      char        **btype,
		      char        **vdwparams);

    /* Return 1 if OK, 0 if not found */
    int  getPtypePol(  const char *ptype,
			 double *polarizability, double *sigPol);
    int  getAtypePol(  const char *atype,
			 double *polarizability, double *sigPol);

    int getAtypeRefEnthalpy(Poldata * pd, const char *atype,
			       double *Href);

    /* Return 1 if OK, 0 if not found */
    int  bondingRuleValence(  char *gtBrule, double *valence);

    void  addMiller(
		     char         *miller,
		     int           atomnumber,
		     double        tauAhc,
		     double        alphaAhp);

    /* Return 1 if "miller" was found */
    int  getMillerPol(
			char         *miller,
			int          *atomnumber,
			double       *tauAhc,
			double       *alphaAhp);

    int  getMiller(
		    char        **miller,
		    int          *atomnumber,
		    double       *tauAhc,
		    double       *alphaAhp);

    void  setMillerUnits(  char *tauUnit,
			     char *ahpUnit);

    void  getMillerUnits(  char **tauUnit,
			     char **ahpUnit);

    /* Returns miller name or NULL if not found */
    char * ptypeToMiller(  const char *ptype);

    void  addBosque(
		     char         *bosque,
		     double        polarizability);

    int  getBosque(
		    char        **bosque,
		    double       *polarizability);

    void  setBosqueUnit(  char *polarUnit);

    char * getBosqueUnit( );

    /* Returns bosque name or NULL if not found */
    char * ptypeToBosque(  const char *ptype);

    int  getBosquePol(  char *bosque, double *polarizability);

    /* Return 1 on success or 0 otherwise */
    int  addBond(  char *atom1, char *atom2,
		    double length, double sigma, int ntrain,
		    double bondorder, char *params);

    int  setBondParams(  char *atom1, char *atom2,
			   double length, double sigma, int ntrain,
			   double bondorder, char *params);

    /* Return bond-index 1-N or 0 if not found */
    int  getBond(  char **atom1, char **atom2,
		    double *length, double *sigma, int *ntrain,
		    double *bondorder, char **params);

    void  setBondFunction(  char *fn);
    char * getBondFunction( );
    int  getBondFtype( );

    /* Return bond-index 1-N or 0 if not found */
    int  searchBond(  char *atom1, char *atom2,
		       double *length, double *sigma, int *ntrain,
		       double *bondorder, char **params);

    /* Returns 1 if there is a bond, 0 if not. Toler is absolute in length-units. */
    int  elemIsBond(  char *elem1, char *elem2,
			double distance, double toler);

    /* Return maximal valence for a give element */
    double  elemGetMaxValence(  char *elem);

    /* Return NULL-terminated array of potential bondorders */
    double * elemGetBondorders(  char *elem1, char *elem2,
				   double distance, double toler);
    /* Returns the bondorder. Toler is absolute in length-units. */
    double  atypeBondorder(  char *atype1, char *atype2,
			      double distance, double toler);

    void  setAngleFunction(  char *fn);
    char * getAngleFunction( );
    int  getAngleFtype( );

    /* Return 1 on success, 0 otherwise */
    int  addAngle(
		   char *atom1, char *atom2,
		   char *atom3, double angle, double sigma,
		   int ntrain, char *params);

    int  setAngleParams(  char *atom1, char *atom2,
			    char *atom3,
			    double angle, double sigma, int ntrain, char *params);

    /* Return angle-index 1-N or 0 if not found */
    int  getAngle(  char **atom1, char **atom2,
		     char **atom3, double *angle, double *sigma,
		     int *ntrain, char **params);

    /* Return angle-index 1-N or 0 if not found */
    int  searchAngle(  char *atom1, char *atom2,
			char *atom3, double *angle, double *sigma,
			int *ntrain, char **params);

    void  setAngleUnit(  char *angleUnit);

    char * getAngleUnit();

    void  setDihedralFunction(  int egd, char *fn);
    char * getDihedralFunction(  int egd);
    int  getDihedralFtype(  int egd);

    /* Return 1 on success or 0 otherwise */
    int  addDihedral(  int egd, char *atom1, char *atom2,
			char *atom3, char *atom4,
			double dihedral, double sigma,
			int ntrain, char *params);

    int  setDihedralParams(  int egd,
			       char *atom1, char *atom2,
			       char *atom3, char *atom4,
			       double angle, double sigma,
			       int ntrain, char *params);

    /* Return dihedral-index 1-N or 0 if not found */
    int  getDihedral(  int egd,
			char **atom1, char **atom2,
			char **atom3, char **atom4,
			double *dihedral, double *sigma,
			int *ntrain, char **params);

    /* Return dihedral-index 1-N or 0 if not found */
    int  searchDihedral(  int egd,
			   char *atom1, char *atom2,
			   char *atom3, char *atom4,
			   double *dihedral, double *sigma,
			   int *ntrain, char **params);

    void  setDihedralUnit(  int   egd,
			      char *dihedralUnit);

    char * getDihedralUnit(  int egd);

    void  addSymcharges(  char *central,
			   char *attached, int numattach);

    int  getSymcharges(  char **central,
			  char **attached, int *numattach);

    int  searchSymcharges(  char *central,
			     char *attached, int numattach);

    static ChargeDistributionModel name2eemtype(const char *name);

    static const char *getEemtypeName(ChargeDistributionModel eem);

    char * getEemref(  ChargeDistributionModel eqdModel);

    int  getNumprops(  ChargeDistributionModel eqdModel);

    int  havePolSupport(  const char *atype);

    int  haveEemSupport(  ChargeDistributionModel eqdModel,
			    const char             *name,
			    gmx_bool                bAllowZeroParameters);

    double  getJ00(  ChargeDistributionModel eqdModel, char *name);

    int  getNzeta(  ChargeDistributionModel eqdModel, char *name);

    double  getZeta(  ChargeDistributionModel eqdModel, char *name, int zz);

    char * getQstr(  ChargeDistributionModel eqdModel, char *name);

    char * getRowstr(  ChargeDistributionModel eqdModel, char *name);

    double  getQ(  ChargeDistributionModel eqdModel, char *name, int zz);

    int  getRow(  ChargeDistributionModel eqdModel, char *name, int zz);

    double  getChi0(  ChargeDistributionModel eqdModel, char *name);

    char * getOpts(  ChargeDistributionModel eqdModel, char *name);

    void  setEemprops(
		       ChargeDistributionModel eqdModel, char *name,
		       double J0, double chi0,
		       char *zeta, char *q, char *row);

    int  getEemprops(
		      ChargeDistributionModel *eqdModel, char **name,
		      double *J0, double *chi0,
		      char **zeta, char **q, char **row);

    void  setEpref(  ChargeDistributionModel eqdModel, char *epref);

    char * getEpref(ChargeDistributionModel eqdModel);

    int listEpref(  ChargeDistributionModel *eqdModel, char **epref);

    void  commEemprops(  t_commrec *cr);

    void  commForceParameters( t_commrec *cr);


  private:

    char          *_filename;
    int            _nptype, _nptypeC;
    std::vector<t_ptype>       _ptype;
    int            _nalexandria, _nalexandriaC;
    std::vector<t_ffatype>  _alexandria;
    int            _nbtype;
    std::vector<char *>         _btype;
    int            _nbrule, _nbruleC;
    std::vector<t_brule> _brule;
    char          *_alexandriaPolarUnit;
    char          *_alexandriaPolarRef;
    char          *_alexandriaForcefield;
    int            _nexcl;
    double         _fudgeQQ, _fudgeLJ;
    char          *_gtVdwFunction, *_gtCombinationRule;
    int            _gtVdwFtype, _gtCombRule;
    char          *_gtBondFunction;
    int            _ngtBond, _ngtBondC, _gtBondFtype;
    char          *_gtLengthUnit;
    std::vector<t_gt_bond>  _gtBond;
    char          *_gtAngleFunction;
    int            _ngtAngle, _ngtAngleC, _gtAngleFtype;
    char          *_gtAngleUnit;
    std::vector<t_gt_angle>     _gtAngle;
    std::vector<char*>          _gtDihedralFunction;
    std::vector<int>            _ngtDihedral, _ngtDihedralC, _gtDihedralFtype;
    std::vector<std::vector<t_gt_dihedral> > _gtDihedral;
    int            _nmiller, _nmillerC;
    std::vector<t_miller>   _miller;
    char          *_millerTauUnit, *_millerAhpUnit;
    int            _nbosque, _nbosqueC;
    std::vector<t_bosque>     _bosque;
    char          *_bosquePolarUnit;
    int            _nsymcharges, _nsymchargesC;
    std::vector<t_symcharges>  _symcharges;
    int            _nep, _nepC;
    std::vector<t_eemprops>  _eep;
    int            _ner, _nerC;
    std::vector<t_epref>   _epr;

    void addBtype(const char   *btype);

    gmx_bool strcasestrStart(char *needle, char *haystack);

    int countNeighbors(t_brule *brule, int nbond, char *nbhybrid[], int *score);

    t_gt_bond *searchBond(  char *atom1, char *atom2,
			     double bondorder);

    int searchBondtype(  char *atom);

    t_eemprops *getEep(ChargeDistributionModel eqdModel,
			const char             *name);

    t_gt_dihedral *searchDihedral( int egd, char *atom1, char *atom2,
				    char *atom3, char *atom4);

    static int gtbComp(const void *a, const void *b);

    static int gtdComp(const void *a, const void *b);

    template<class Type>
    int indexOfPointInVector(Type * pointer, std::vector<Type> vector){
    return (pointer - &(vector[0]));
    }
  };
}
#endif
