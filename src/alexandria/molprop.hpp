/*
 * $Id: molprop.h,v 1.16 2009/05/29 15:01:18 spoel Exp $
 * 
 *                This source code is part of
 * 
 *                 G   R   O   M   A   C   S
 * 
 *          GROningen MAchine for Chemical Simulations
 * 
 *                        VERSION 4.0.99
 * Written by David van der Spoel, Erik Lindahl, Berk Hess, and others.
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2008, The GROMACS development team,
 * check out http://www.gromacs.org for more information.

 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 * 
 * If you want to redistribute modifications, please consider that
 * scientific software is very special. Version control is crucial -
 * bugs must be traceable. We will be happy to consider code for
 * inclusion in the official distribution, but derived work must not
 * be called official GROMACS. Details are found in the README & COPYING
 * files - if they are missing, get the official version at www.gromacs.org.
 * 
 * To help us fund GROMACS development, we humbly ask that you cite
 * the papers on the package - you can find them in the top README file.
 * 
 * For more info, check our website at http://www.gromacs.org
 * 
 * And Hey:
 * Groningen Machine for Chemical Simulation
 */

#ifndef _molprop_h
#define _molprop_h

enum { eMOLPROP_Exp, eMOLPROP_Calc, eMOLPROP_Any, eMOLPROP_NR };

enum { empPOTENTIAL, empDIPOLE, empQUADRUPOLE, empPOLARIZABILITY, 
       empENERGY, empNR };

extern const char *emp_name[empNR];

#define assign_str(char *dst,char *src)  if (NULL != src) { if (NULL != dst) *dst = strdup(src); } else { *dst = NULL; }
#define assign_scal(dst,src) if (NULL != dst) *dst = src

namespace gmx 
{
  class MolecularComposition 
  {
  private:
    std::string _compname;
    int _catom_index;
    std::vector<std::string> _catom;
    std::vector<int> _cnumber;
  public:
    MolecularComposition(char *compname) { _compname.assign(compname); }
    MolecularComposition(std::string compname) { _compname.assign(compname); }
    ~MolecularComposition() {}
    void AddAtom(const char *catom,const int cnumber);
    void AddAtom(std::string catom,const int cnumber);
    void DeleteAtom(const char *catom);
    void DeleteAtom(std::string catom);
    void ReplaceAtom(const char *oldatom,const char *newatom);
    void ReplaceAtom(std::string oldatom,std::string newatom);
    void Reset() { _catom_index = 0; }
    const char *GetAtom(int *cnumber);
    const std::string GetAtom(int *cnumber);
    int CountAtoms(const char *atom);
    int CountAtoms(std::string atom);
  };

  class GenericProperty 
  {
  private:
    std::string _type,_unit;
  public:
    GenericProperty(char *type,char *unit) { _type.assign(type); _unit.assign(unit); };
    GenericProperty(std::string type,std::string unit) { _type = type; _unit = unit; };
    ~GenericProperty() {};
    const char *GetType() { return _type.c_str(); }
    const char *GetUnit() { return _unit.c_str(); }
    std::string GetType() { return _type; }
    std::string GetUnit() { return _unit; }
  };
  
  class MolecularQuadrupole : GenericProperty
  {
  private:
    double _xx,_yy,_zz,_xy,_xz,_yz;
  public:
    MolecularQuadrupole(char *type,char *unit,double xx,double yy,double zz,double xy,double xz,double yz) : GenericProperty(type,unit) { Set(xx,yy,zz,xy,xz,yz); };
    MolecularQuadrupole(std::string type,std::string unit,double xx,double yy,double zz,double xy,double xz,double yz) : GenericProperty(type,unit) { Set(xx,yy,zz,xy,xz,yz); };
    ~MolecularQuadrupole() {};
    Set(double xx,double yy,double zz,double xy,double xz,double yz) { _xx=xx; _yy=yy; _zz=zz; _xy=xy; _xz=xz; _yz=yz; };
    Get(double *xx,double *yy,double *zz,double *xy,double *xz,double *yz) { *xx=_xx; *yy=_yy; *zz=_zz; *xy=_xy; *xz=_xz; *yz=_yz; };
  };

  class MolecularEnergy : GenericProperty
  {
  private:
    double _value,_error;
  public:
    MolecularEnergy(char *type,char *unit,double value,double error) : GenericProperty(type,unit) { Set(value,error); };
    MolecularEnergy(std::string type,std::string unit,double value,double error) : GenericProperty(type,unit) { Set(value,error); };
    ~MolecularEnergy() {};
    Set(double value,double error) { _value = value; _error = error; };
    Get(double *value,double *error) { *vaue = _value; *error = _error };
  };

  class MolecularDipPolar : GenericProperty
  {
  private:
    double _x,_y,_z;
    double _aver,_error;
  public:
    MolecularDipPolar(char *type,char *unit,double x,double y,double z,double aver,double error) : GenericProperty(type,unit) { Set(x,y,z,aver,error); }
    MolecularDipPolar(std::string type,std::string unit,double x,double y,double z,double aver,double error) : GenericProperty(type,unit) { Set(x,y,z,aver,error); }
    ~MolecularDipPolar() {};
    Set(double x,double y,double z,double aver,double error) { _x = x; _y = y; _z = z; _aver = aver; _error = error; };
    Get(double *x,double *y,double *z,double *aver,double *error) { *x = _x; *y = _y; *z = _z; *aver = _aver; *error = _error; };
  };

  class ElectrostaticPotential
  {  
  private:
    std::string _xyz_unit,_V_unit;
    int _espid;
    double _x,_y,_z,_V;
  public:
    ElectrostaticPotential(const char *xyz_unit,const char *V_unit,int espid,double x,double y,double z,double V) { Set(xyz_unit,V_unit,espid,x,y,z,V); };
    ~ElectrostaticPotential() {};
    Set(const char *xyz_unit,const char *V_unit,int espid,double x,double y,double z,double V) { _xyz_unit.assign(xyz_unit); _V_unit.assign(V_unit); _espid = espid; _x = x; _y = y; _z = z; _V = V; };
    Get(char **xyz_unit,char **V_unit,int *espid,double *x,double *y,double *z,double *V) { _xyz_unit.assign(xyz_unit); _V_unit.assign(V_unit); _espid = espid; _x = x; _y = y; _z = z; _V = V; };
  };

  class Bond
  {
  private:
    int _ai,_aj,_bondorder;
  public:
    Bond(int ai,int aj,int bondorder) { Set(ai,aj,bondorder); }
    ~Bond() {};
    Set(int ai,int aj,int bondorder) {_ai = ai; _aj = aj; _bondorder = bondorder };
    Get(int *ai,int *aj,int *bondorder) { *ai = _ai; *aj = _aj; *bondorder = _bondorder };
  };

  class AtomicCharge : GenericProperty
  {
  private:
    double _q;
  public:
    AtomicCharge(const char *type,const char *unit,double q) : GenericProperty(type,unit) { Set(q); };
    AtomicCharge(std::string type,std::string unit,double q) : GenericProperty(type,unit) { Set(q); };
    ~AtomicCharge() {};
    void Set(double q) { _q = q };
    void Get(double *q) { *q = _q; }
  };
  typedef std::vector<AtomicCharge>::iterator AtomicChargeIterator; 

  class CalcAtom 
  {
  private:
    std::string _name,_obtype,_unit;
    double _x,_y,_z;
    int _atomid,_charge_index;
    std::vector<AtomicCharge> _q;
  public:
    CalcAtom(const char *name,const char *obtype,const char *unit) {};
    ~CalcAtom() {};
    void AddCharge(const char *type,const char *unit,double q) { _q.push_back(type,unit,q); }
    /* Get through iterator */
    AtomicChargeIterator BeginQ() { return _q.begin(); }
    AtomicChargeIterator EndQ() { return _q.end(); }
  };

  class Experiment
  {
  private:
    int _eMP;
    std::string _reference,_conformation;
    int _nprop_c[empNR];
    std::vector<MolecularDipPolar> _polar,_dipole;
    std::vector<MolecularEnergy> _energy;
    std::vector<MolecularQuadrupole> _quadrupole;
  public:
    Experiment(std::string reference,std::string conformation) 
    { 
      _reference = reference;
      _conformation = conformation;
    };
    Experiment(const char *reference,const char *conformation) 
    { 
      _reference.assign(reference);
      _conformation.assign(conformation);
    };
    ~Experiment() {};
    void AddPolar(const char *type,const char *unit,
                  double xx,double yy,double zz,
                  double aver,double error) { _polar.push_back(type,unit,xx,yy,zz,aver,error); }
    void AddEnergy(const char *type,const char *unit,
                   double value,double error) { _energy.push_back(type,unit,value,error); };
    /* Iterator ? */
    int GetPolar(char **type,char **unit,
                 double *xx,double *yy,double *zz,double *aver,double *error);
    int GetEnergy(char **type,char **unit,
                  double *value,double *error);
  };
  
  class Calculation : Experiment
  {
  private:
    std::string _program,_method,_basisset,_datafile;
    int _catom_index;
    std::vector<CalcAtom> _catom;
    std::vector<ElectrostaticPotential> _potential;
  public:
    Calculation() {};
    ~Calculation() {};
    CalcAtom GetAtom();

    void AddPotential(const char *xyzUnit,const char *VUnit,int espid,
                      double x,double y,double z,
                      double V) { _potential.push_back(xyzUnit,VUnit,espid,x,y,z,V); };
    int GetNPotential() { return _potential.size(); };

    /* Iterator! */
    int getPotential(char **xyzUnit,char **VUnit,int *espid,
                     double *x,double *y,double *z,
                     double *V);
  };

  class MolProp {
  private:
    double _mass;
    int _charge,_multiplicity;
    std::string _formula,_molname,_iupac,_cas,_cid,_inchi;
    int _category_index;
    std::vector<std::string> _category;
    int _mol_comp_index;
    std::vector<MolecularComposition> _mol_comp;
    int _calc_index;
    std::vector<Calculation> _calc;
    int _exper_index;
    std::vector<Experiment> _exper;
    int _bond_index;
    std::vector<Bond> _bond;
  public:
    MolProp() {}
    ~MolProp() {}
    
    void CheckConsistency();
    
    void SetMass(double mass) { _mass = mass; }
    double GetMass() { return _mass; }
    
    void SetCharge(double charge) { _charge = charge; }
    int GetCharge() { return _charge; }
    
    void SetMultiplicity(int multiplicity) { _multiplicity = multiplicity; }
    int GetMultiplicity() { return _multiplicity; }
    
    void SetFormula(const char *formula) { _formula.assign(formula); }
    void SetFormula(std::string formula) { _formula.assign(formula); }
    char *GetFormula() { return _formula.c_str(); }
    std::string GetFormula() { return _formula; }
    
    void SetMolname(const char *molname) { _molname.assign(molname); }
    void SetMolname(std::string molname) { _molname.assign(molname); }
    char *GetMolname() { return _molname.c_str(); }
    std::string GetMolname() { return _molname; }
    
    void SetIupac(const char *iupac) { _iupac.assign(iupac); }
    void SetIupac(std::string iupac) { _iupac.assign(iupac); }
    
    /* Return IUPAC name or, if not found, the molname */
    char *GetIupac() { if (_iupac.size() > 0) return _iupac.c_str() else return _molname.c_str(); }
    std::string GetIupac() { if (_iupac.size() > 0) return _iupac else return _molname; }
    
    void SetCas(const char *cas) { _cas.assign(cas); }
    void SetCas(std::string cas) { _cas.assign(cas); }
    char *GetCas() { return _cas.c_str(); }
    std::string GetCas() { return _cas; }
    
    void SetCid(const char *cid) { _cid.assign(cid); }
    void SetCid(std::string cid) { _cid.assign(cid); }
    char *GetCid() { return _cid.c_str(); }
    std::string GetCid() { return _cid; }
    
    void SetInchi(const char *inchi) { _inchi.assign(inchi); }
    void SetInchi(std::string inchi) { _inchi.assign(inchi); }
    char *GetInchi() { return _inchi.c_str(); }
    std::string GetInchi() { return _inchi; }
    
    void AddCategory(const char *category) { _category.push_back(category); }
    void AddCategory(std::string category) { _category.push_back(category); }
    int NCategory() { return _category.size(); }
    
    /* Returns one category at a time. If NULL, you got them all previously.
     * Implement with iterators.
     */
    const char *GetCategory();
    std::string GetCategory();
    void ResetCategory() { _category_index = 0; }
    int SearchCategory(const char *catname);
    int SearchCategory(std::string catname);
    
    void DeleteComposition(const char *compname);
    MolecularComposition AddComposition(const char *compname) { _mol_comp.push_back(compname); return _mol_comp.back(); }
    void ResetComposition() { _mol_comp_index = 0; }
    MolecularComposition GetComposition() { if (_mol_comp_index < _mol_comp.size()) return _mol_comp[_mol_comp_index++]; else { ResetComposition(); return NULL; } }
      
    int GetNAtom() { if (_mol_comp.size() > 0) return _mol_comp[0].CountAtoms(); else return 0; };
  
    int GetBond(int *ai,int *aj,int *bondorder) { if(_bond_index < _bonds.size()) { _bonds[_bond_index].Get(ai,aj,*bondorder); _bond_index++; return 1 } else { _bond_index=0; return 0; } };
    void AddBond(int ai,int aj,int bondorder) { _bonds.push_back(ai,aj,bondorder); }

    Experiment AddExperiment(const char *reference,const char *conformation) { _exper.push_back(reference,conformation); return _exper.back(); };

    /* Replace by iterator ? */
    int GetExperiment(char **reference,char **conformation);
    void ResetExperiment() {};


/* Potential. Ref is a calc ref. from addCalc, or getCalc. */	    

/* Dipoles. Ref is either exp or calc ref. from addExper,
   addCalc, getExper or getCalc. */	    
extern void addDipole(int ref,
				   const char *type,const char *unit,
				   double x,double y,double z,
				   double aver,double error);
				       
extern int getDipole(int ref,
				  char **type,char **unit,
				  double *x,double *y,double *z,
				  double *aver,double *error);

/* Quadrupoles. Ref is either exp or calc ref. from addExper,
   addCalc, getExper or getCalc. */	    
extern void addQuadrupole(int ref,
				       const char *type,const char *unit,
				       double xx,double yy,double zz,
				       double xy,double xz,double yz);
				       
extern int getQuadrupole(int ref,
				      char **type,char **unit,
				      double *xx,double *yy,double *zz,
				      double *xy,double *xz,double *yz);

extern void resetCalculation(t mpt);

extern void reset(t mpt);

extern int getNbond(t mpt);

  
/* Returns calcref that can be used to add properties later on */
extern void addCalculation(
					const char *program,const char *method,
					const char *basisset,const char *reference,
					const char *conformation,const char *datafile,
                                        int *calcref);

extern int getCalculation(char **program,char **method,
				       char **basisset,char **reference,
				       char **conformation,char **datafile,
                                       int *calcref);

/* Returns atomref that can be used to set coordinates and charges */
extern void calcAddAtom(int calcref,
				      const char *atomname,const char *obtype,int atomid,int *atomref);

extern int calcGetNatom(int calcref);

extern int calcGetAtom(int calcref,
				     char **atomname,char **obtype,int *atomid,int *atomref);

extern void calcSetAtomcoords(int calcref,int atomref,
					    const char *unit,double x,double y,double z);
					   
extern int calcGetAtomcoords(int calcref,int atomref,
					   char **unit,double *x,double *y,double *z);
					   
extern void calcSetAtomcharge(int calcref,int atomref,
					   const char *type,const char *unit,double q);
					   
extern int calcGetAtomcharge(int calcref,int atomref,
					   char **type,char **unit,double *q);

/* Generic routines */					   
extern t copy(t mpt);

extern void merge(t dst,t src);

#ifdef __cplusplus
}
#endif

#endif
