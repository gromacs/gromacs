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

#include <string.h>
#include <string>
#include <vector>
#include "typedefs.h"

enum { eMOLPROP_Exp, eMOLPROP_Calc, eMOLPROP_Any, eMOLPROP_NR };

enum { empPOTENTIAL, empDIPOLE, empQUADRUPOLE, empPOLARIZABILITY, 
       empENERGY, empNR };

extern const char *emp_name[empNR];

#define assign_str(dst,src)  if (NULL != src) { if (NULL != dst) *dst = strdup(src); } else { *dst = NULL; }
#define assign_scal(dst,src) if (NULL != dst) *dst = src

namespace alexandria 
{
  class AtomNum 
  {
  private:
    std::string _catom;
    int _cnumber;
  public:
    AtomNum(const char *catom,int cnumber) { SetAtom(catom); SetNumber(cnumber); }
    AtomNum(std::string catom,int cnumber) { SetAtom(catom); SetNumber(cnumber); }
    ~AtomNum() {};
    std::string GetAtom() { return _catom; }
    void SetAtom(std::string catom) { _catom = catom; }
    void SetAtom(const char *catom) { _catom.assign(catom); }
    //const char *GetAtom() { return _catom.c_str(); }
    int GetNumber() { return _cnumber; }
    void SetNumber(int cnumber) { _cnumber = cnumber; }
  };
  typedef std::vector<AtomNum>::iterator AtomNumIterator; 
  
  class MolecularComposition 
  {
  private:
    std::string _compname;
    std::vector<AtomNum> _atomnum;
  public:
    MolecularComposition(const char *compname) { _compname.assign(compname); }
    MolecularComposition(std::string compname) { _compname = compname; }
    ~MolecularComposition() {}
    std::string CompName() { return _compname; }
    void AddAtom(const char *catom,int cnumber) { _atomnum.push_back(AtomNum(catom,cnumber); }
    void AddAtom(std::string catom,int cnumber) { _atomnum.push_back(AtomNum(catom,cnumber)); }
    void DeleteAtom(const char *catom);
    void DeleteAtom(std::string catom);
    void ReplaceAtom(const char *oldatom,const char *newatom) { std::string so(oldatom),sn(newatom); ReplaceAtom(so,sn); }
    void ReplaceAtom(std::string oldatom,std::string newatom);
    AtomNumIterator BeginAtomNum() { return _atomnum.begin(); }
    AtomNumIterator EndAtomNum() { return _atomnum.end(); }
    int CountAtoms(const char *atom);
    int CountAtoms(std::string atom);
    int CountAtoms();
  };
  typedef std::vector<MolecularComposition>::iterator MolecularCompositionIterator;
  
  class GenericProperty 
  {
  private:
    std::string _type,_unit;
  public:
    GenericProperty(char *type,char *unit) { SetType(type); SetUnit(unit); };
    GenericProperty(std::string type,std::string unit) { SetType(type); SetUnit(unit); };
    ~GenericProperty() {};
    //const char *GetType() { return _type.c_str(); }
    //const char *GetUnit() { return _unit.c_str(); }
    std::string GetType() { return _type; }
    std::string GetUnit() { return _unit; }
    void SetType(std::string type) { _type = type; }
    void SetUnit(std::string unit) { _unit = unit; }
    void SetType(const char *type) { _type.assign(type); }
    void SetUnit(const char *unit) { _unit.assign(unit); }
  };
  
  class MolecularQuadrupole : public GenericProperty
  {
  private:
    double _xx,_yy,_zz,_xy,_xz,_yz;
  public:
    MolecularQuadrupole(char *type,char *unit,double xx,double yy,double zz,double xy,double xz,double yz) : GenericProperty(type,unit) { Set(xx,yy,zz,xy,xz,yz); };
    MolecularQuadrupole(std::string type,std::string unit,double xx,double yy,double zz,double xy,double xz,double yz) : GenericProperty(type,unit) { Set(xx,yy,zz,xy,xz,yz); };
    ~MolecularQuadrupole() {};
    void Set(double xx,double yy,double zz,double xy,double xz,double yz) { _xx=xx; _yy=yy; _zz=zz; _xy=xy; _xz=xz; _yz=yz; };
    void Get(double *xx,double *yy,double *zz,double *xy,double *xz,double *yz) { *xx=_xx; *yy=_yy; *zz=_zz; *xy=_xy; *xz=_xz; *yz=_yz; };
  };
  typedef std::vector<MolecularQuadrupole>::iterator MolecularQuadrupoleIterator;
  
  class MolecularEnergy : public GenericProperty
  {
  private:
    double _value,_error;
  public:
    MolecularEnergy(char *type,char *unit,double value,double error) : GenericProperty(type,unit) { Set(value,error); };
    MolecularEnergy(std::string type,std::string unit,double value,double error) : GenericProperty(type,unit) { Set(value,error); };
    ~MolecularEnergy() {};
    void Set(double value,double error) { _value = value; _error = error; };
    void Get(double *value,double *error) { *value = _value; *error = _error; };
  };
  typedef std::vector<MolecularEnergy>::iterator MolecularEnergyIterator;
  
  class MolecularDipPolar : public GenericProperty
  {
  private:
    double _x,_y,_z;
    double _aver,_error;
  public:
    MolecularDipPolar(char *type,char *unit,double x,double y,double z,double aver,double error) : GenericProperty(type,unit) { Set(x,y,z,aver,error); }
    MolecularDipPolar(std::string type,std::string unit,double x,double y,double z,double aver,double error) : GenericProperty(type,unit) { Set(x,y,z,aver,error); }
    ~MolecularDipPolar() {};
    void Set(double x,double y,double z,double aver,double error) { _x = x; _y = y; _z = z; _aver = aver; _error = error; };
    void Get(double *x,double *y,double *z,double *aver,double *error) { *x = _x; *y = _y; *z = _z; *aver = _aver; *error = _error; };
  };
  typedef std::vector<MolecularDipPolar>::iterator MolecularDipPolarIterator;
  
  class ElectrostaticPotential
  {  
  private:
    std::string _xyz_unit,_V_unit;
    int _espid;
    double _x,_y,_z,_V;
  public:
    ElectrostaticPotential(const char *xyz_unit,const char *V_unit,int espid,double x,double y,double z,double V) { Set(xyz_unit,V_unit,espid,x,y,z,V); };
    ~ElectrostaticPotential() {};
    void Set(const char *xyz_unit,const char *V_unit,int espid,double x,double y,double z,double V) { _xyz_unit.assign(xyz_unit); _V_unit.assign(V_unit); _espid = espid; _x = x; _y = y; _z = z; _V = V; };
    void Get(char **xyz_unit,char **V_unit,int *espid,double *x,double *y,double *z,double *V) { 
      *xyz_unit = strdup(_xyz_unit.c_str()); *V_unit = strdup(_V_unit.c_str()); *espid = _espid; *x = _x; *y = _y; *z = _z; *V = _V; 
    };
  };
  typedef std::vector<ElectrostaticPotential>::iterator ElectrostaticPotentialIterator;

  class Bond
  {
  private:
    int _ai,_aj,_bondorder;
  public:
    Bond(int ai,int aj,int bondorder) { Set(ai,aj,bondorder); }
    ~Bond() {};
    void Set(int ai,int aj,int bondorder) {_ai = ai; _aj = aj; _bondorder = bondorder; };
    void Get(int *ai,int *aj,int *bondorder) { *ai = _ai; *aj = _aj; *bondorder = _bondorder; };
    int GetAi() { return _ai; }
    int GetAj() { return _aj; }
    int GetBondOrder() { return _bondorder; }
  };
  typedef std::vector<Bond>::iterator BondIterator; 

  class AtomicCharge : public GenericProperty
  {
  private:
    double _q;
  public:
    AtomicCharge(const char *type,const char *unit,double q) : GenericProperty(type,unit) { SetQ(q); };
    AtomicCharge(std::string type,std::string unit,double q) : GenericProperty(type,unit) { SetQ(q); };
    ~AtomicCharge() {};
    void SetQ(double q) { _q = q; };
    double GetQ() { return _q; }
  };
  typedef std::vector<AtomicCharge>::iterator AtomicChargeIterator; 

  class CalcAtom 
  {
  private:
    std::string _name,_obtype,_unit;
    double _x,_y,_z;
    int _atomid;
    std::vector<AtomicCharge> _q;
  public:
    CalcAtom(const char *name,const char *obtype,int atomid) { 
      _name.assign(name); _obtype.assign(obtype); _atomid = atomid;
    };
    ~CalcAtom() {};
    
    void AddCharge(const char *type,const char *unit,double q) { 
      _q.push_back(AtomicCharge(type,unit,q)); 
    }
    /* Get through iterator */
    AtomicChargeIterator BeginQ() { return _q.begin(); }
    AtomicChargeIterator EndQ() { return _q.end(); }
    
    //char *GetUnit() { return _unit.c_str(); }
    int GetAtomid() { return _atomid; }
    std::string GetName() { return _name; }
    std::string GetObtype() { return _obtype; }
    std::string GetUnit() { return _unit; }
    void SetUnit(std::string unit) { _unit = unit; }
    void SetUnit(const char *unit) { _unit.assign(unit); }
        
    void SetCoords(double x,double y,double z) { _x = x; _y = y; _z = z; }
    void GetCoords(double *x,double *y,double *z) { *x = _x; *y = _y; *z = _z; }
  };
  typedef std::vector<CalcAtom>::iterator CalcAtomIterator;
   
  class Experiment
  {
  private:
    gmx_bool _bExperiment;
    std::string _reference,_conformation;
    std::vector<MolecularDipPolar> _polar,_dipole;
    std::vector<MolecularEnergy> _energy;
    std::vector<MolecularQuadrupole> _quadrupole;
  public:
    Experiment(std::string reference,std::string conformation) { 
      _reference = reference; _conformation = conformation; _bExperiment = TRUE;
    };
    Experiment(const char *reference,const char *conformation) { 
      _reference.assign(reference); _conformation.assign(conformation);
    };
    ~Experiment() {};
    gmx_bool GetExperiment() { return _bExperiment; }
    void SetExperiment(gmx_bool bExperiment) { _bExperiment = bExperiment; }

    void AddPolar(MolecularDipPolar mdp) { _polar.push_back(mdp); }
    MolecularDipPolarIterator BeginPolar() { return _polar.begin(); }
    MolecularDipPolarIterator EndPolar()   { return _polar.end(); }
                  
    void AddDipole(MolecularDipPolar mdp) { _dipole.push_back(mdp); }
    MolecularDipPolarIterator BeginDipole() { return _dipole.begin(); }
    MolecularDipPolarIterator EndDipole()   { return _dipole.end(); }
                  
    void AddQuadrupole(MolecularQuadrupole mq) { _quadrupole.push_back(mq); }
    MolecularQuadrupoleIterator BeginQuadrupole() { return _quadrupole.begin(); }
    MolecularQuadrupoleIterator EndQuadrupole() { return _quadrupole.end(); }

    void AddEnergy(MolecularEnergy me) { _energy.push_back(me); }
    MolecularEnergyIterator BeginEnergy() { return _energy.begin(); }
    MolecularEnergyIterator EndEnergy()   { return _energy.end(); }
    
    std::string GetConformation() { return _conformation; }
    //const char *GetConformation() { return _conformation.c_str(); }
    
    std::string GetReference() { return _reference; }
    //const char *GetReference() { return _reference.c_str(); }
  };
  typedef std::vector<Experiment>::iterator ExperimentIterator; 
  
  class Calculation : public Experiment
  {
  private:
    std::string _program,_method,_basisset,_datafile;
    std::vector<CalcAtom> _catom;
    std::vector<ElectrostaticPotential> _potential;
  public:
    Calculation(const char *program,const char *method,
                const char *basisset,const char *reference,
                const char *conformation,const char *datafile) : Experiment(reference,conformation) { 
      _program.assign(program); _method.assign(method); 
      _basisset.assign(basisset); _datafile.assign(datafile); 
      SetExperiment(FALSE);
    };
    ~Calculation() {};
    
    void AddAtom(CalcAtom ca) { _catom.push_back(ca); }
    int NAtom() { return _catom.size(); }
    CalcAtomIterator BeginAtom() { return _catom.begin(); }
    CalcAtomIterator EndAtom() { return _catom.end(); }

    void AddPotential(ElectrostaticPotential ep) { _potential.push_back(ep); }
    int NPotential() { return _potential.size(); };
    ElectrostaticPotentialIterator BeginPotential() { return _potential.begin(); }
    ElectrostaticPotentialIterator EndPotential() { return _potential.end(); }
    
    std::string GetProgram() { return _program; }
    //const char *GetProgram() { return _program.c_str(); }
    
    std::string GetBasisset() { return _basisset; }
    //const char *GetBasisset() { return _basisset.c_str(); }
    
    std::string GetMethod() { return _method; }
    //const char *GetMethod() { return _method.c_str(); }
    
    std::string GetDatafile() { return _datafile; }
    //const char *GetDatafile() { return _datafile.c_str(); }
  };
  typedef std::vector<Calculation>::iterator CalculationIterator;
  
  class MolProp {
  private:
    double _mass;
    int _charge,_multiplicity;
    std::string _formula,_molname,_iupac,_cas,_cid,_inchi;
    std::vector<std::string> _category;
    std::vector<MolecularComposition> _mol_comp;
    std::vector<Calculation> _calc;
    std::vector<Experiment> _exper;
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
    //char *GetFormula() { return _formula.c_str(); }
    std::string GetFormula() { return _formula; }
    
    void SetMolname(const char *molname) { _molname.assign(molname); }
    void SetMolname(std::string molname) { _molname.assign(molname); }
    //char *GetMolname() { return _molname.c_str(); }
    std::string GetMolname() { return _molname; }
    
    void SetIupac(const char *iupac) { _iupac.assign(iupac); }
    void SetIupac(std::string iupac) { _iupac.assign(iupac); }
    
    /* Return IUPAC name or, if not found, the molname */
    //char *GetIupac() { if (_iupac.size() > 0) return _iupac.c_str() else return _molname.c_str(); }
    std::string GetIupac() { if (_iupac.size() > 0) return _iupac; else return _molname; }
    
    void SetCas(const char *cas) { _cas.assign(cas); }
    void SetCas(std::string cas) { _cas.assign(cas); }
    //char *GetCas() { return _cas.c_str(); }
    std::string GetCas() { return _cas; }
    
    void SetCid(const char *cid) { _cid.assign(cid); }
    void SetCid(std::string cid) { _cid.assign(cid); }
    //char *GetCid() { return _cid.c_str(); }
    std::string GetCid() { return _cid; }
    
    void SetInchi(const char *inchi) { _inchi.assign(inchi); }
    void SetInchi(std::string inchi) { _inchi.assign(inchi); }
    //char *GetInchi() { return _inchi.c_str(); }
    std::string GetInchi() { return _inchi; }
    
    void AddCategory(const char *category) { _category.push_back(category); }
    void AddCategory(std::string category) { _category.push_back(category); }
    int NCategory() { return _category.size(); }
    
    std::vector<std::string>::iterator BeginCategory() { return _category.begin(); }
    std::vector<std::string>::iterator EndCategory() { return _category.end(); }

    int SearchCategory(const char *catname) { std::string str(catname); return SearchCategory(str); }
    int SearchCategory(std::string catname);
    
    void DeleteComposition(std::string compname);
    void DeleteComposition(const char *compname) { 
      std::string str(compname); DeleteComposition(compname); 
    }
    MolecularCompositionIterator AddComposition(const char *compname) { 
      _mol_comp.push_back(MolecularComposition(compname)); return _mol_comp.end(); 
    }
    MolecularCompositionIterator BeginMolecularComposition() { return _mol_comp.begin(); }
    MolecularCompositionIterator EndMolecularComposition()   { return _mol_comp.end(); }
    
    int NAtom() { if (_mol_comp.size() > 0) return _mol_comp[0].CountAtoms(); else return 0; };

    void AddBond(Bond b) { _bond.push_back(b); }
    int Nbond() { return _bond.size(); }
    BondIterator BeginBond() { return _bond.begin(); }
    BondIterator EndBond() { return _bond.end(); }
    
    void AddExperiment(Experiment myexp) { _exper.push_back(myexp); };
    int NExperiment() { return _exper.size(); }

    ExperimentIterator BeginExperiment() { return _exper.begin(); }
    ExperimentIterator EndExperiment() { return _exper.end(); }
    Experiment *LastExperiment() { 
      if (NExperiment() > 0) return &(_exper.back()); else return NULL; 
    }

    void AddCalculation(Calculation calc) { _calc.push_back(calc); }
    int NCalculation() { return _calc.size(); }
    CalculationIterator BeginCalculation() { return _calc.begin(); }
    CalculationIterator EndCalculation() { return _calc.end(); }
    Calculation *LastCalculation() { 
      if (NCalculation() > 0) return &(_calc.back()); else return NULL; 
    }
  };
  typedef std::vector<MolProp>::iterator MolPropIterator;
}
  
#endif
