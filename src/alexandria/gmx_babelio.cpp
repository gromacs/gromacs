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

typedef struct {
  int  id;
  real x,y,z;
  real V;
} t_espv;

static int esp_comp(const void *a,const void *b)
{
  t_espv *ea = (t_espv *)a;
  t_espv *eb = (t_espv *)b;
  
  if (ea->V < eb->V)
    return -1;
  else if (ea->V > ea->V)
    return 1;
  return 0;
}
 
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
                                     char *basisset,gau_atomprop_t gaps,
                                     real th_toler,real ph_toler,
                                     int maxpot,gmx_bool bVerbose)
{
  /* Read a gaussian log file */
  OpenBabel::OBMol mol;
  OpenBabel::OBAtomIterator OBai;
  OpenBabel::OBBondIterator OBbi;
  OpenBabel::OBConversion conv;
  OpenBabel::OBAtom *OBa,*OBa1,*OBa2;
  OpenBabel::OBBond *OBb;
  OpenBabel::OBPairData *OBpd;
  OpenBabel::OBGenericData *OBdhf;
  OpenBabel::OBVectorData *dipole;
  OpenBabel::OBMatrixData *quadrupole,*pol_tensor;
  OpenBabel::OBFreeGrid *esp;
  OpenBabel::OBElementTable *OBet;
  std::string formula,attr,value,inchi;
  
  t_espv *espv = NULL;
  const char *reference="Spoel2013a",*unknown="unknown";
  char *program,*method,*basis,*charge_model,*ptr,*g98ptr;
  int i,ii0,calcref,atomref,atomid,bondid,natom;
  gmx_molprop_t mpt;  
  double ii,deltai,dval;
  int k;
  const char *etypes[] = { "DHf(0K)", "DHf(298.15K)" };
  
  if (0 != read_babel(g98,&mol))
    {
      fprintf(stderr,"Failed reading %s\n",g98);
      return NULL;
    }

  mol.PerceiveBondOrders();
  /* Create new calculation */ 
  mpt = gmx_molprop_init();

  OBpd = (OpenBabel::OBPairData *)mol.GetData("basis");
  if (NULL != basisset)
    basis = basisset;
  else if (NULL != OBpd)
    {
      basis = strdup(OBpd->GetValue().c_str());
      if (NULL != (ptr = strstr(basis," (5D, 7F)")))
        *ptr = '\0';
    }
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
  
  g98ptr = (char *) strrchr(g98,'/');
  if (NULL == g98ptr) 
    g98ptr = (char *)g98;
  else {
    g98ptr++;
    if (strlen(g98ptr) == 0)
      g98ptr = (char *)g98;
  }
  gmx_molprop_add_calculation(mpt,program,method,basis,reference,
                              conformation,g98ptr,&calcref);
  
  gmx_molprop_set_charge(mpt,mol.GetTotalCharge());
  gmx_molprop_set_mass(mpt,mol.GetMolWt());
  gmx_molprop_set_multiplicity(mpt,mol.GetTotalSpinMultiplicity());
  
  formula = mol.GetFormula();
  if (formula.size() > 0)
    gmx_molprop_set_formula(mpt,formula.c_str());
  else
    gmx_molprop_set_formula(mpt,unknown);
  
  conv.SetOutFormat("inchi");
  inchi = conv.WriteString(&mol,true);

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
          dval = convert2gmx(atof(value.c_str()),eg2c_kcal_mole);
          gmx_molprop_add_energy(mpt,calcref,etypes[k],
                                 unit2string(eg2c_kJ_mole),dval,0);
        }
    }
  
  /* Now add properties by extracting them from the OpenBabel structure */
  OBpd = (OpenBabel::OBPairData *) mol.GetData("PartialCharges");
  if (NULL != OBpd)
    charge_model = strdup(OBpd->GetValue().c_str());
  else
    charge_model = strdup(unknown);
  
  OBet = new OpenBabel::OBElementTable();
  
  OBai = mol.BeginAtoms();
  atomid = 1;
  for (OBa = mol.BeginAtom(OBai); (NULL != OBa); OBa = mol.NextAtom(OBai)) {
    gmx_molprop_calc_add_atom(mpt,calcref,
                              OBet->GetSymbol(OBa->GetAtomicNum()),
                              OBa->GetType(),atomid,&atomref);
    gmx_molprop_calc_set_atomcoords(mpt,calcref,atomref,unit2string(eg2c_pm),
                                    100*OBa->x(),100*OBa->y(),100*OBa->z());
    gmx_molprop_calc_set_atomcharge(mpt,calcref,atomref,charge_model,
                                    "e",OBa->GetPartialCharge());
    atomid++;
  }
  delete OBet;
  
  OBbi = mol.BeginBonds();
  bondid = 1;
  for (OBb = mol.BeginBond(OBbi); (NULL != OBb); OBb = mol.NextBond(OBbi)) {
    OBa1 = OBb->GetBeginAtom();
    OBa2 = OBb->GetEndAtom(); 
    gmx_molprop_add_bond(mpt,1+OBa1->GetIndex(),1+OBa2->GetIndex(),OBb->GetBondOrder());
    bondid++;
  }
  
  // Dipole
  dipole = (OpenBabel::OBVectorData *) mol.GetData("Dipole Moment");
  if (NULL != dipole)
    {
      OpenBabel::vector3 v3 = dipole->GetData();
      gmx_molprop_add_dipole(mpt,calcref,"electronic",
                             unit2string(eg2c_Angstrom),
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
      gmx_molprop_add_quadrupole(mpt,calcref,"electronic",
                                 unit2string(eg2c_Buckingham),
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
      fac = 1000*pow(convert2gmx(1,eg2c_Bohr),3);
      for(i=0; (i<9); i++)
        mm[i] *= fac; 
      //cout << "fac = " << fac << "\n";
      alpha = (mm[0]+mm[4]+mm[8])/3.0;
      gmx_molprop_add_polar(mpt,calcref,"electronic",
                            unit2string(eg2c_Angstrom3),
                            mm[0],mm[4],mm[8],alpha,0);
    }
  
  // Electrostatic potential
  esp = (OpenBabel::OBFreeGrid *) mol.GetData("Electrostatic Potential");
  if (NULL != esp)
    {
      OpenBabel::OBFreeGridPoint *fgp;
      OpenBabel::OBFreeGridPointIterator fgpi;
      int espid=0;
      
      fgpi = esp->BeginPoints();
      snew(espv,esp->NumPoints());
      for(fgp = esp->BeginPoint(fgpi); (NULL != fgp); fgp = esp->NextPoint(fgpi))
        {
          espv[espid].id = espid+1;
          espv[espid].x = 100*fgp->GetX();
          espv[espid].y = 100*fgp->GetY();
          espv[espid].z = 100*fgp->GetZ();
          espv[espid].V = fgp->GetV();
          espid++;
        }
      natom = mol.NumAtoms();
      
      if ((maxpot > 0) && (maxpot < esp->NumPoints())) 
        {
          qsort(espv+natom,espid-natom,sizeof(espv[0]),esp_comp);
          deltai = (esp->NumPoints()-natom)/maxpot;
        }
      else 
        {
          maxpot = esp->NumPoints();
          deltai = 1;
        }
      for(i=0; (i<natom); i++) 
        {
          gmx_molprop_add_potential(mpt,calcref,
                                    unit2string(eg2c_pm),
                                    unit2string(eg2c_Hartree_e),
                                    espv[i].id,
                                    espv[i].x,espv[i].y,espv[i].z,
                                    espv[i].V);
        }
      ii = natom;
      for(; (i<maxpot); i++) 
        {
          /* Convert to integer */
          ii0 = ii;
          gmx_molprop_add_potential(mpt,calcref,
                                    unit2string(eg2c_pm),
                                    unit2string(eg2c_Hartree_e),
                                    i,
                                    espv[ii0].x,espv[ii0].y,espv[ii0].z,
                                    espv[ii0].V);
          ii+=deltai;
      }
      sfree(espv);
    }
  
  return mpt;
}

#endif
