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
#include <openbabel/math/vector3.h>

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
  OpenBabel::OBConversion conv;
  OpenBabel::OBAtom *OBa;
  OpenBabel::OBPairData *OBpd;
  OpenBabel::OBGenericData *OBdhf;
  OpenBabel::OBVectorData *dipole;
  OpenBabel::OBMatrixData *quadrupole,*pol_tensor;
  OpenBabel::OBFreeGrid *esp;
  std::string formula,attr,value,inchi;
  
  const char *reference="Spoel2013a",*unknown="unknown";
  char *program,*method,*basis,*charge_model;
  int calcref,atomref,atomid;
  gmx_molprop_t mpt;  
  int k;
  const char *etypes[] = { "DHf(0K)", "DHf(298.15K)" };
  
  if (0 != read_babel(g98,&mol))
    gmx_fatal(FARGS,"Failed reading %s",g98);

  //  ...manipulate molecule
  //cout << " Molecule has: " << mol.NumAtoms()
  //     << " atoms." << endl;
  mol.PerceiveBondOrders();
  /* Create new calculation */ 
  mpt = gmx_molprop_init();

  OBpd = (OpenBabel::OBPairData *)mol.GetData("basis");
  if (NULL != OBpd)
    basis = strdup(OBpd->GetValue().c_str());
  else
    basis = strdup(unknown);
  
  OBpd = (OpenBabel::OBPairData *)mol.GetData("program");
  if (NULL != OBpd)
    program = strdup(OBpd->GetValue().c_str());
  else
    program = strdup(unknown);
  
  OBpd = (OpenBabel::OBPairData *)mol.GetData("method");
  if (NULL != OBpd)
    method = strdup(OBpd->GetValue().c_str());
  else 
    method = strdup(unknown);
  
  gmx_molprop_add_calculation(mpt,program,method,basis,reference,
                              conformation,&calcref);
  gmx_molprop_set_charge(mpt,mol.GetTotalCharge());
  gmx_molprop_set_mass(mpt,mol.GetMolWt());
  gmx_molprop_set_multiplicity(mpt,mol.GetTotalSpinMultiplicity());
  
  formula = mol.GetFormula();
  if (formula.size() > 0)
    gmx_molprop_set_formula(mpt,formula.c_str());
  else
    gmx_molprop_set_formula(mpt,unknown);
    
  
  conv.SetOutFormat("inchi");
  inchi = conv.WriteString(&mol);

  if (inchi.size() > 0)
    gmx_molprop_set_inchi(mpt,inchi.c_str());
  else
    gmx_molprop_set_inchi(mpt,unknown);
  
  if (NULL != molnm)
    gmx_molprop_set_molname(mpt,molnm);
  else
    gmx_molprop_set_molname(mpt,unknown);
    
  if (NULL != iupac)
    gmx_molprop_set_iupac(mpt,iupac);
  else
    gmx_molprop_set_iupac(mpt,unknown);
    
  for(k=0; (k<2); k++) 
    {
      OBdhf = mol.GetData(etypes[k]);
      if (NULL != OBdhf)
        {
          value = OBdhf->GetValue();
          gmx_molprop_add_energy(mpt,calcref,etypes[k],"kJ/mol",
                                 convert2gmx(atof(value.c_str()),
                                             eg2cKcal_Mole),0);
        }
    }
  
  /* Now add properties by extracting them from the OpenBabel structure */
  OBpd = (OpenBabel::OBPairData *) mol.GetData("PartialCharges");
  if (NULL != OBpd)
    charge_model = strdup(OBpd->GetValue().c_str());
  else
    charge_model = strdup(unknown);
    
  OBai = mol.BeginAtoms();
  atomid = 1;
  for (OBa = mol.BeginAtom(OBai); (NULL != OBa); OBa = mol.NextAtom(OBai)) {
    gmx_molprop_calc_add_atom(mpt,calcref,OBa->GetType(),atomid,&atomref);
    gmx_molprop_calc_set_atomcoords(mpt,calcref,atomref,"Angstrom",
                                    OBa->x(),OBa->y(),OBa->z());
    gmx_molprop_calc_set_atomcharge(mpt,calcref,atomref,charge_model,
                                    "e",OBa->GetPartialCharge());
    atomid++;
  }
  
  // Dipole
  dipole = (OpenBabel::OBVectorData *) mol.GetData("Dipole Moment");
  if (NULL != dipole)
    {
      OpenBabel::vector3 v3 = dipole->GetData();
      gmx_molprop_add_dipole(mpt,calcref,"electronic","Debye",
                             v3.GetX(),v3.GetY(),v3.GetZ(),
                             v3.length(),0.0);
    }
  
  // Quadrupole
  quadrupole = (OpenBabel::OBMatrixData *) mol.GetData("Traceless Quadrupole Moment");
  if (NULL != quadrupole)
    {
      OpenBabel::matrix3x3 m3 = quadrupole->GetData();
      double mm[9];
      m3.GetArray(mm);
      gmx_molprop_add_quadrupole(mpt,calcref,"electronic","Debye-Ang",
                                 mm[0],mm[4],mm[8],mm[1],mm[2],mm[5]);
    }
  
  // Polarizability
  pol_tensor = (OpenBabel::OBMatrixData *) mol.GetData("Exact polarizability");
  if (NULL != pol_tensor)
    {
      OpenBabel::matrix3x3 m3 = pol_tensor->GetData();
      double mm[9],alpha,fac;
      int i;
      m3.GetArray(mm);
      fac = 1000*pow(convert2gmx(1,eg2cBohr),3);
      for(i=0; (i<9); i++)
        mm[i] *= fac; 
      //cout << "fac = " << fac << "\n";
      alpha = (mm[0]+mm[4]+mm[8])/3.0;
      gmx_molprop_add_polar(mpt,calcref,"electronic","Angstrom^3",
                            mm[0],mm[4],mm[8],alpha,0);
    }
  
  // Electrostatic potential
  esp = (OpenBabel::OBFreeGrid *) mol.GetData("Electrostatic Potential");
  if (NULL != esp)
    {
      OpenBabel::OBFreeGridPoint *fgp;
      OpenBabel::OBFreeGridPointIterator fgpi;
      int espid=1;
      
      fgpi = esp->BeginPoints();
      for(fgp = esp->BeginPoint(fgpi); (NULL != fgp); fgp = esp->NextPoint(fgpi))
        {
          gmx_molprop_add_potential(mpt,calcref,
                                    unit2string(eg2cAngstrom),
                                    unit2string(eg2cHartree_e),espid,
                                    fgp->GetX(),fgp->GetY(),fgp->GetZ(),fgp->GetV());
          espid++;
        }
    }
  
  
  return mpt;
}

#endif
