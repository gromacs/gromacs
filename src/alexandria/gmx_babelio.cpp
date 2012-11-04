#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include "typedefs.h"
#include "gmx_fatal.h"
#include "smalloc.h"
#include "futil.h"
#include "poldata.h"
#include "atomprop.h"
#include "molprop.h"
#include "string2.h"
#include "gmx_babelio.h"
#include "vec.h"

// Include Open Babel classes for OBMol and OBConversion
#ifdef HAVE_LIBOPENBABEL2
#include <iostream>
#include <fstream>
#include <openbabel/mol.h>
#include <openbabel/atom.h>
#include <openbabel/residue.h>
#include <openbabel/obiter.h>
#include <openbabel/obconversion.h>

using namespace std;

int read_babel(const char *g98,OpenBabel::OBMol *mol)
{
  ifstream         g98f;
  char             *g98z;
  
  if (FALSE == gmx_fexist(g98)) 
    {
      snew(g98z,strlen(g98)+4);
      strcpy(g98z,g98);
      strcat(g98z,".gz");
      g98f.open(g98z,ios::in);
      sfree(g98z);
    }
  else 
    {
      g98f.open(g98,ios::in);
    }
  if (!g98f.is_open()) 
    {
      gmx_fatal(FARGS,"Can not open file %s for reading",g98);
    }
  
  // Read from g98f
  OpenBabel::OBConversion conv(&g98f,&cout);
  
  // Try to set input format to G98 file
  if (conv.SetInFormat("G98"))
    {
      if (conv.Read(mol))
        {
          g98f.close();
          return 0; // exit with success
        }
      else 
        {
          cerr << "Could not read input file " << g98 << " with OpenBabel2." << endl;
        }
    }
  else 
    {
      cerr << "Input file " << g98 << " has incorrect Gaussian98 format." << endl;
    }
  g98f.close();
  return -1;
}

gmx_molprop_t gmx_molprop_read_gauss(const char *g98,
                                     gmx_atomprop_t aps,gmx_poldata_t pd,
                                     char *molnm,char *iupac,char *conformation,
                                     gau_atomprop_t gaps,
                                     real th_toler,real ph_toler,
                                     int maxpot,gmx_bool bVerbose)
{
  /* Read a gaussian log file */
  OpenBabel::OBMol mol;
  OpenBabel::OBAtomIterator OBai;
  OpenBabel::OBAtom *OBa;
  OpenBabel::OBGenericData *OBdata;
  std::string formula;
  
  const char *program="Gaussian",*method="GX",*basisset="GX",
    *reference="Spoel2013a",*atomname;
  int calcref,atomref,atomid;
  gmx_molprop_t mpt;  
  int i;
  
  if (0 != read_babel(g98,&mol))
    gmx_fatal(FARGS,"Failed reading %s",g98);

  //  ...manipulate molecule
  cout << " Molecule has: " << mol.NumAtoms()
       << " atoms." << endl;
  mol.PerceiveBondOrders();
  /* Create new calculation */ 
  mpt = gmx_molprop_init();

  gmx_molprop_add_calculation(mpt,program,method,basisset,reference,
                              conformation,&calcref);
  gmx_molprop_set_charge(mpt,mol.GetTotalCharge());
  gmx_molprop_set_mass(mpt,mol.GetMolWt());
  gmx_molprop_set_multiplicity(mpt,mol.GetTotalSpinMultiplicity());
  formula = mol.GetFormula();
  gmx_molprop_set_formula(mpt,formula.c_str());
  gmx_molprop_add_energy(mpt,calcref,
                         "Heat of Formation","kcal/mol",
                         mol.GetEnergy(),0);
  OBdata = mol.GetData("G4 Enthalpy");
  if (NULL != OBdata) 
    {
      gmx_molprop_add_energy(mpt,calcref,
                             "G4 Enthalpy","Hartree",
                             atof(OBdata->GetValue().c_str()),0);
      
    }
  /* Now add properties by extracting them from the OpenBabel structure */
  OBai = mol.BeginAtoms();
  atomid = 1;
  for (OBa = mol.BeginAtom(OBai); (NULL != OBa); OBa = mol.NextAtom(OBai)) {
    gmx_molprop_calc_add_atom(mpt,calcref,OBa->GetType(),atomid,&atomref);
    gmx_molprop_calc_set_atomcoords(mpt,calcref,atomref,"Angstrom",
                                    OBa->x(),OBa->y(),OBa->z());
    
    atomid++;
  }

  return mpt;
}

#endif
