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
#include "gmxpre.h"
#include <iostream>
#include <fstream>
#include <algorithm>
#include "gromacs/utility/real.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/smalloc.h"
#include "gromacs/fileio/strdb.h"
#include "gromacs/utility/futil.h"
#include "gromacs/topology/symtab.h"
#include "gromacs/utility/cstringutil.h"
#include "gromacs/math/vec.h"
#include "gromacs/topology/atomprop.h"
#include "poldata.h"
#include "molprop.h"
#include "molprop_util.h"
#include "gauss_io.h"
#include "stringutil.h"

using namespace std;

static bool comp_esp(alexandria::ElectrostaticPotential ea,
                     alexandria::ElectrostaticPotential eb)
{
    if (ea.getV() < eb.getV())
    {
        return true;
    }
    else
    {
        return false;
    }
}

static void merge_electrostatic_potential(alexandria::MolProp &mpt,
                                          std::vector<alexandria::ElectrostaticPotential> &espv,
                                          int natom, int maxpot)
{
    alexandria::ElectrostaticPotentialIterator esi;
    int i;

    if ((maxpot > 0) && (maxpot < (int)espv.size()))
    {
        std::sort(espv.begin()+natom, espv.end(), comp_esp);
    }
    else
    {
        maxpot = 1;
    }
    i  = 0;
    for (esi = espv.begin(); (esi < espv.end()); esi++, i++)
    {
        if ((i < natom) || (((i-natom) % maxpot) == 0))
        {
            mpt.LastCalculation()->AddPotential(*esi);
        }
    }
}

// Include Open Babel classes for OBMol and OBConversion
#ifdef HAVE_LIBOPENBABEL2
// Hack to make this compile!
#undef ANGSTROM
#include <openbabel/babelconfig.h>
#include <openbabel/obmolecformat.h>
#include <openbabel/mol.h>
#include <openbabel/atom.h>
#include <openbabel/data.h>
#include <openbabel/residue.h>
#include <openbabel/obiter.h>
#include <openbabel/obconversion.h>
#include <openbabel/forcefield.h>
#include <openbabel/math/vector3.h>

static OpenBabel::OBConversion *read_babel(const char *g98, OpenBabel::OBMol *mol)
{
    ifstream          g98f;
    char             *g98z, *ptr;

    if (FALSE == gmx_fexist(g98))
    {
        snew(g98z, strlen(g98)+4);
        strcpy(g98z, g98);
        strcat(g98z, ".gz");
        g98f.open(g98z, ios::in);
        sfree(g98z);
    }
    else
    {
        g98f.open(g98, ios::in);
    }
    if (!g98f.is_open())
    {
        gmx_fatal(FARGS, "Can not open file %s for reading", g98);
    }

    // Read from g98f
    OpenBabel::OBConversion *conv = new OpenBabel::OBConversion(&g98f, &cout);

    // Try to set input format to G98 file if it is not clear from the extension,
    // that means, this code will equally well read sdf, pdb etc.
    ptr = (char *)strrchr(g98, '.');
    if ((NULL != ptr) && (strlen(ptr) >= 2))
    {
        ptr++;
    }
    else
    {
        ptr = (char *)"g98";
    }
    if (conv->SetInFormat(ptr))
    {
        if (conv->Read(mol))
        {
            g98f.close();

            return conv; // exit with success
        }
        else
        {
            cerr << "Could not read input file " << g98 << " with OpenBabel2." << endl;
        }
    }
    else
    {
        cerr << "Input file " << g98 << " has incomprehensible format." << endl;
    }
    g98f.close();

    return NULL;
}

static void gmx_molprop_read_babel(const char *g98,
                                   alexandria::MolProp &mpt,
                                   char *molnm, char *iupac, char *conformation,
                                   char *basisset, int maxpot,
                                   int nsymm,
                                   const char *forcefield)
{
    /* Read a gaussian log file */
    OpenBabel::OBMol           mol, mol2;
    OpenBabel::OBAtomIterator  OBai;
    OpenBabel::OBBondIterator  OBbi;
    OpenBabel::OBConversion   *conv;
    //OpenBabel::OBAtom *OBa;
    OpenBabel::OBBond         *OBb;
    OpenBabel::OBPairData     *OBpd;
    OpenBabel::OBVectorData   *dipole;
    OpenBabel::OBMatrixData   *quadrupole, *pol_tensor;
    OpenBabel::OBFreeGrid     *esp;
    OpenBabel::OBElementTable *OBet;
    std::string                formula, attr, value, inchi;

    std::vector<alexandria::ElectrostaticPotential> espv;

    const char *reference = "Ghahremanpour2015a", *unknown = "unknown";
    char       *program, *method, *basis, *charge_model, *ptr, *g98ptr;
    int         bondid;

    conv = read_babel(g98, &mol);
    if (NULL == conv)
    {
        fprintf(stderr, "Failed reading %s\n", g98);
        return;
    }
    delete conv;

    conv = new OpenBabel::OBConversion(&cin, &cout);
    // Now extract classification info.
    if (conv->SetOutFormat("fpt"))
    {
        vector<string> vs;
        string         ss;
        const char    *exclude[] = { ">", "C_ONS_bond", "Rotatable_bond", "Conjugated_double_bond", "Conjugated_triple_bond", "Chiral_center_specified", "Cis_double_bond", "Bridged_rings", "Conjugated_tripple_bond", "Trans_double_bond" };
#define nexclude (sizeof(exclude)/sizeof(exclude[0]))
        char          *dup, *ptr;
        unsigned int   i, j;

        conv->AddOption("f", OpenBabel::OBConversion::OUTOPTIONS, "FP4");
        conv->AddOption("s");
        conv->Convert();
        mol2 = mol;
        ss   = conv->WriteString(&mol2, false);
        if (OpenBabel::tokenize(vs, ss))
        {
            for (i = 0; (i < vs.size()); i++)
            {
                for (j = 0; (j < nexclude); j++)
                {
                    if (strcasecmp(exclude[j], vs[i].c_str()) == 0)
                    {
                        break;
                    }
                }
                if (j == nexclude)
                {
                    dup = strdup(vs[i].c_str());
                    while (NULL != (ptr = strchr(dup, '_')))
                    {
                        *ptr = ' ';
                    }
                    mpt.AddCategory(dup);
                }
            }
        }
    }

    // get bondorders.
    mol.PerceiveBondOrders();

    OBpd = (OpenBabel::OBPairData *)mol.GetData("basis");
    if ((NULL != basisset) && (strlen(basisset) > 0))
    {
        basis = basisset;
    }
    else if (NULL != OBpd)
    {
        basis = strdup(OBpd->GetValue().c_str());
        if (NULL != (ptr = strstr(basis, " (5D, 7F)")))
        {
            *ptr = '\0';
        }
    }
    else
    {
        basis = strdup(unknown);
    }

    OBpd = (OpenBabel::OBPairData *)mol.GetData("program");
    if (NULL != OBpd)
    {
        program = strdup(OBpd->GetValue().c_str());
    }
    else
    {
        program = strdup(unknown);
    }

    OBpd = (OpenBabel::OBPairData *)mol.GetData("method");
    if (NULL != OBpd)
    {
        method = strdup(OBpd->GetValue().c_str());
    }
    else
    {
        method = strdup(unknown);
    }
    g98ptr = (char *) strrchr(g98, '/');
    if (NULL == g98ptr)
    {
        g98ptr = (char *)g98;
    }
    else
    {
        g98ptr++;
        if (strlen(g98ptr) == 0)
        {
            g98ptr = (char *)g98;
        }
    }
    alexandria::Calculation ca(program, method, basis, reference,
                               conformation, g98ptr);
    mpt.AddCalculation(ca);

    mpt.SetCharge(mol.GetTotalCharge());
    mpt.SetMass(mol.GetMolWt());
    mpt.SetMultiplicity(mol.GetTotalSpinMultiplicity());
    mpt.SetFormula(mol.GetFormula());

    conv->SetOutFormat("inchi");
    inchi = conv->WriteString(&mol, true);

    delete conv;

    mpt.SetInchi(inchi);

    if (NULL != molnm)
    {
        mpt.SetMolname(molnm);
    }
    else
    {
        mpt.SetMolname(unknown);
    }

    if (NULL != iupac)
    {
        mpt.SetIupac(iupac);
    }
    else
    {
        mpt.SetIupac(unknown);
    }

    {
        double temperature, DeltaHf0, DeltaHfT, DeltaGfT, DeltaSfT, S0T;
        if (extract_thermochemistry(mol, false, &nsymm,
                                    &temperature,
                                    &DeltaHf0,
                                    &DeltaHfT,
                                    &DeltaGfT,
                                    &DeltaSfT,
                                    &S0T))
        {
            alexandria::MolecularEnergy me1("DeltaHform",
                                            unit2string(eg2cKj_Mole),
                                            0,
                                            epGAS,
                                            convert2gmx(DeltaHf0, eg2cKcal_Mole),
                                            0);
            mpt.LastCalculation()->AddEnergy(me1);
            alexandria::MolecularEnergy me2("DeltaHform",
                                            unit2string(eg2cKj_Mole),
                                            temperature,
                                            epGAS,
                                            convert2gmx(DeltaHfT, eg2cKcal_Mole),
                                            0);
            mpt.LastCalculation()->AddEnergy(me2);
            alexandria::MolecularEnergy me3("DeltaGform",
                                            unit2string(eg2cKj_Mole),
                                            temperature,
                                            epGAS,
                                            convert2gmx(DeltaGfT, eg2cKcal_Mole),
                                            0);
            mpt.LastCalculation()->AddEnergy(me3);
            alexandria::MolecularEnergy me4("DeltaSform",
                                            unit2string(eg2cJ_MolK),
                                            temperature,
                                            epGAS,
                                            convert2gmx(DeltaSfT, eg2cCal_MolK),
                                            0);
            mpt.LastCalculation()->AddEnergy(me4);
            alexandria::MolecularEnergy me5("S0",
                                            unit2string(eg2cJ_MolK),
                                            temperature,
                                            epGAS,
                                            convert2gmx(S0T, eg2cCal_MolK),
                                            0);
            mpt.LastCalculation()->AddEnergy(me5);
        }
    }

    /* Now add properties by extracting them from the OpenBabel structure */
    OBpd = (OpenBabel::OBPairData *) mol.GetData("PartialCharges");
    if (NULL != OBpd)
    {
        charge_model = strdup(OBpd->GetValue().c_str());
    }
    else
    {
        charge_model = strdup(unknown);
    }

    OBet = new OpenBabel::OBElementTable();

    OpenBabel::OBForceField *ff = OpenBabel::OBForceField::FindForceField(forcefield);
    if ((NULL != ff) && (ff->Setup(mol)))
    {
        ff->GetAtomTypes(mol);
        FOR_ATOMS_OF_MOL (atom, mol) {
            OpenBabel::OBPairData *type = (OpenBabel::OBPairData*) atom->GetData("FFAtomType");
            if (NULL == type)
            {
                gmx_fatal(FARGS, "Could not find %s atom type for atom %s",
                          forcefield, atom->GetIdx());
            }
            if (NULL != debug)
            {
                fprintf(debug, "XXX atom %d gafftype %s OBtype %s\n",
                        atom->GetIdx(), type->GetValue().c_str(), atom->GetType());
            }
            alexandria::CalcAtom     ca(OBet->GetSymbol(atom->GetAtomicNum()),
                                        type->GetValue(), atom->GetIdx());
            alexandria::AtomicCharge aq(charge_model, "e", 0.0,
                                        atom->GetPartialCharge());

            ca.SetUnit(unit2string(eg2cPm));
            ca.SetCoords(100*atom->x(), 100*atom->y(), 100*atom->z());
            ca.AddCharge(aq);
            mpt.LastCalculation()->AddAtom(ca);
        }
        // Not necessary to delete?
        //delete ff;
    }
    else
    {
        gmx_fatal(FARGS, "Can not read %s force field", forcefield);
    }
    delete OBet;

    OBbi   = mol.BeginBonds();
    bondid = 1;
    for (OBb = mol.BeginBond(OBbi); (NULL != OBb); OBb = mol.NextBond(OBbi))
    {
        alexandria::Bond ab(1+OBb->GetBeginAtom()->GetIndex(),
                            1+OBb->GetEndAtom()->GetIndex(),
                            OBb->GetBondOrder());
        mpt.AddBond(ab);
        bondid++;
    }

    // Dipole
    dipole = (OpenBabel::OBVectorData *) mol.GetData("Dipole Moment");
    if (NULL != dipole)
    {
        OpenBabel::vector3            v3 = dipole->GetData();
        alexandria::MolecularDipole   dp("electronic",
                                         unit2string(eg2cDebye),
                                         0.0,
                                         v3.GetX(), v3.GetY(), v3.GetZ(),
                                         v3.length(), 0.0);
        mpt.LastCalculation()->AddDipole(dp);
    }

    // Quadrupole
    quadrupole = (OpenBabel::OBMatrixData *) mol.GetData("Traceless Quadrupole Moment");
    if (NULL != quadrupole)
    {
        OpenBabel::matrix3x3            m3 = quadrupole->GetData();
        double                          mm[9];
        m3.GetArray(mm);
        alexandria::MolecularQuadrupole mq("electronic",
                                           unit2string(eg2cBuckingham),
                                           0.0,
                                           mm[0], mm[4], mm[8],
                                           mm[1], mm[2], mm[5]);
        mpt.LastCalculation()->AddQuadrupole(mq);
    }

    // Polarizability
    pol_tensor = (OpenBabel::OBMatrixData *) mol.GetData("Exact polarizability");
    if (NULL != pol_tensor)
    {
        OpenBabel::matrix3x3 m3 = pol_tensor->GetData();
        double               mm[9], alpha, fac;
        int                  i;
        m3.GetArray(mm);
        fac = 1000*pow(convert2gmx(1, eg2cBohr), 3);
        for (i = 0; (i < 9); i++)
        {
            mm[i] *= fac;
        }
        //cout << "fac = " << fac << "\n";
        alpha = (mm[0]+mm[4]+mm[8])/3.0;

        alexandria::MolecularPolarizability mdp("electronic",
                                                unit2string(eg2cAngstrom3),
                                                0.0,
                                                mm[0], mm[4], mm[8], mm[1], mm[2], mm[5], alpha, 0);
        mpt.LastCalculation()->AddPolar(mdp);
    }

    // Electrostatic potential
    esp = (OpenBabel::OBFreeGrid *) mol.GetData("Electrostatic Potential");
    if (NULL != esp)
    {
        OpenBabel::OBFreeGridPoint        *fgp;
        OpenBabel::OBFreeGridPointIterator fgpi;
        std::string xyz_unit(unit2string(eg2cPm));
        std::string V_unit(unit2string(eg2cHartree_e));
        int         espid = 0;

        fgpi = esp->BeginPoints();
        for (fgp = esp->BeginPoint(fgpi); (NULL != fgp); fgp = esp->NextPoint(fgpi))
        {
            alexandria::ElectrostaticPotential ep(xyz_unit, V_unit, ++espid,
                                                  100*fgp->GetX(),
                                                  100*fgp->GetY(),
                                                  100*fgp->GetZ(),
                                                  fgp->GetV());
            espv.push_back(ep);
        }
        merge_electrostatic_potential(mpt, espv, mol.NumAtoms(), maxpot);
    }
}

#endif

/*******************************************************************
 *******************************************************************
 * Code for in case we do not have openbabel
 *******************************************************************
 *******************************************************************/

static int get_lib_file(const char *db, char ***strings)
{
    FILE  *in;
    char **ptr = NULL;
    char   buf[STRLEN];
    int    i, nstr, maxi;

    in = libopen(db);

    i = maxi = 0;
    while (fgets2(buf, STRLEN-1, in))
    {
        if (i >= maxi)
        {
            maxi += 50;
            srenew(ptr, maxi);
        }
        ptr[i] = strdup(buf);
        i++;
    }
    nstr = i;
    gmx_ffclose(in);
    srenew(ptr, nstr);
    *strings = ptr;

    return nstr;
}

/* read composite method atom data */
alexandria::GaussAtomProp::GaussAtomProp()
{
    char **strings = NULL;
    int    nstrings, i;

    nstrings = get_lib_file("atomization_energies.dat", &strings);

    /* First loop over strings to deduct basic stuff */
    for (i = 0; (i < nstrings); i++)
    {
        if (strings[i][0] == '#')
        {
            continue;
        }
        std::vector<std::string> ptr = split(strings[i], '|');
        if ((ptr.size() >= 4) &&
            (ptr[0].length() > 0) && (ptr[1].length() > 0) &&
            (ptr[2].length() > 0) && (ptr[3].length() > 0) &&
            (ptr[4].length() > 0))
        {
            alexandria::GaussAtomPropVal gapv(ptr[0], ptr[1], ptr[2],
                                              atof(ptr[3].c_str()),
                                              atof(ptr[4].c_str()));

            _gapv.push_back(gapv);
        }
        sfree(strings[i]);
    }
    sfree(strings);
}

int alexandria::GaussAtomProp::getValue(const char *element,
                                        const char *method,
                                        const char *desc,
                                        double      temp,
                                        double     *value)
{
    std::vector<alexandria::GaussAtomPropVal>::iterator g;
    int    found;
    double ttol = 0.01; /* K */

    found = 0;

    for (g = _gapv.begin(); (g < _gapv.end()); g++)
    {
        if ((0 == strcasecmp(g->getElement().c_str(), element)) &&
            (0 == strcasecmp(g->getMethod().c_str(), method)) &&
            (0 == strcasecmp(g->getDesc().c_str(), desc)) &&
            (fabs(temp - g->getTemp()) < ttol))
        {
            *value = g->getValue();
            found  = 1;
            break;
        }
    }
    return found;
}

void ReadGauss(const char *g98,
               alexandria::MolProp &mp,
               alexandria::GaussAtomProp &gap,
               gmx_atomprop_t aps, gmx_poldata_t pd,
               char *molnm, char *iupac, char *conf,
               char *basis,
               int maxpot, int nsymm, gmx_bool bVerbose,
               const char *forcefield)
{
#ifdef HAVE_LIBOPENBABEL2
    gmx_molprop_read_babel(g98, mp, molnm, iupac, conf, basis,
                           maxpot, nsymm, forcefield);
#else
    fprintf(stderr, "For reading Gaussian input you need to link to OpenBabel\n");
#endif
}
