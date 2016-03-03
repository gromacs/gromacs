/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2016, by the GROMACS development team, led by
 * Mark Abraham, David van der Spoel, Berk Hess, and Erik Lindahl,
 * and including many others, as listed in the AUTHORS file in the
 * top-level source directory and at http://www.gromacs.org.
 *
 * GROMACS is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public License
 * as published by the Free Software Foundation; either version 2.1
 * of the License, or (at your option) any later version.
 *
 * GROMACS is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with GROMACS; if not, see
 * http://www.gnu.org/licenses, or write to the Free Software Foundation,
 * Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA.
 *
 * If you want to redistribute modifications to GROMACS, please
 * consider that scientific software is very special. Version
 * control is crucial - bugs must be traceable. We will be happy to
 * consider code for inclusion in the official distribution, but
 * derived work must not be called official GROMACS. Details are found
 * in the README & COPYING files - if they are missing, get the
 * official version at http://www.gromacs.org.
 *
 * To help us fund GROMACS development, we humbly ask that you cite
 * the research papers on the package. Check out http://www.gromacs.org.
 */
/*! \internal \brief
 * Implements part of the alexandria program.
 * \author David van der Spoel <david.vanderspoel@icm.uu.se>
 */
#include "gmxpre.h"

#include "gauss_io.h"

#include "config.h"

#include <cstdio>

#include <algorithm>
#include <fstream>
#include <iostream>

#include "gromacs/math/units.h"
#include "gromacs/math/vec.h"
#include "gromacs/topology/atomprop.h"
#include "gromacs/topology/symtab.h"
#include "gromacs/topology/topology.h"
#include "gromacs/utility/cstringutil.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/futil.h"
#include "gromacs/utility/real.h"
#include "gromacs/utility/stringutil.h"

#include "molprop.h"
#include "molprop_util.h"
#include "poldata.h"
#include "stringutil.h"

static void merge_electrostatic_potential(alexandria::MolProp                             &mpt,
                                          std::vector<alexandria::ElectrostaticPotential> &espv,
                                          int                                              natom, 
                                          int                                              maxpot)
{
    int mymod = 1;
    if ((maxpot > 0) && (maxpot < (int)espv.size()))
    {
        std::sort(espv.begin()+natom, espv.end(), 
                  [](const alexandria::ElectrostaticPotential &a,
                     const alexandria::ElectrostaticPotential &b)
                  {  return (a.getV() < b.getV()); });
                  
        int npot = espv.size() - natom;
        mymod = npot / maxpot;
    }
    
    int i  = 0;
    for (auto esi = espv.begin(); (esi < espv.end()); esi++, i++)
    {
        if ((i < natom) || (((i-natom) % mymod) == 0))
        {
            mpt.LastExperiment()->AddPotential(*esi);
        }
    }
}

// Include Open Babel classes for OBMol and OBConversion
#if HAVE_LIBOPENBABEL2
// Hack to make this compile!
#undef ANGSTROM
#ifdef HAVE_SYS_TIME_H
#define KOKO HAVE_SYS_TIME_H
#undef HAVE_SYS_TIME_H
#endif
#include <openbabel/atom.h>
#include <openbabel/babelconfig.h>
#include <openbabel/data_utilities.h>
#include <openbabel/forcefield.h>
#include <openbabel/mol.h>
#include <openbabel/obconversion.h>
#include <openbabel/obiter.h>
#include <openbabel/obmolecformat.h>
#include <openbabel/residue.h>
#include <openbabel/math/vector3.h>
#ifdef KOKO
#ifndef HAVE_SYS_TIME_H
#define HAVE_SYS_TIME_H KOKO
#endif
#undef KOKO
#endif

static OpenBabel::OBConversion *read_babel(const char *g98, OpenBabel::OBMol *mol)
{
    std::ifstream g98f;
    bool          isGzip = false;

    if (!gmx_fexist(g98))
    {
        std::string g98z(g98);
        g98z += ".gz";
        g98f.open(g98z.c_str(), std::ios::in);
        isGzip = g98f.is_open();
    }
    else
    {
        g98f.open(g98, std::ios::in);
    }
    if (!g98f.is_open())
    {
        gmx_fatal(FARGS, "Can not open file %s for reading", g98);
    }

    // Read from g98f
    OpenBabel::OBConversion *conv = new OpenBabel::OBConversion(&g98f, &std::cout);

    if (conv->SetInFormat("g09", isGzip))
    {
        if (conv->Read(mol, &g98f))
        {
            g98f.close();

            return conv; // exit with success
        }
        else
        {
            fprintf(stderr, "Could not read input file %s with OpenBabel2.\n",
                    g98);
        }
    }
    else
    {
        fprintf(stderr, "Input file %s has incomprehensible format.\n", g98);
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
    OpenBabel::OBMol           mol;
    OpenBabel::OBAtomIterator  OBai;
    OpenBabel::OBBondIterator  OBbi;
    //OpenBabel::OBAtom *OBa;
    OpenBabel::OBBond         *OBb;
    OpenBabel::OBPairData     *OBpd;
    OpenBabel::OBVectorData   *dipole;
    OpenBabel::OBMatrixData   *quadrupole, *pol_tensor;
    OpenBabel::OBFreeGrid     *esp;
    OpenBabel::OBElementTable *OBet;
    std::string                formula, attr, value;

    std::vector<alexandria::ElectrostaticPotential> espv;

    const char *reference = "Ghahremanpour2016a", *unknown = "unknown";
    char       *program, *method, *basis, *charge_model, *ptr, *g98ptr;
    int         bondid;

    OpenBabel::OBConversion *conv = read_babel(g98, &mol);
    if (NULL == conv)
    {
        fprintf(stderr, "Failed reading %s\n", g98);
        return;
    }
    int qtot = mol.GetTotalCharge();
    int mult = mol.GetTotalSpinMultiplicity();
    printf("Total charge %d\n", qtot);
    delete conv;

    conv = new OpenBabel::OBConversion(&std::cin, &std::cout);
    // Now extract classification info.
    if (conv->SetOutFormat("fpt"))
    {
        const char    *exclude[] = { ">", "C_ONS_bond", "Rotatable_bond", "Conjugated_double_bond", "Conjugated_triple_bond", "Chiral_center_specified", "Cis_double_bond", "Bridged_rings", "Conjugated_tripple_bond", "Trans_double_bond" };
#define nexclude (sizeof(exclude)/sizeof(exclude[0]))
        
        conv->AddOption("f", OpenBabel::OBConversion::OUTOPTIONS, "FP4");
        conv->AddOption("s");
        conv->Convert();
        // We need a copy here because WriteString removes the H.
        OpenBabel::OBMol         mol2 = mol;
        std::string              ss = conv->WriteString(&mol2, false);
        std::vector<std::string> vs = gmx::splitString(ss);
        for (const auto &i : vs)
        {
            size_t j;
            for (j = 0; (j < nexclude); j++)
            {
                if (strcasecmp(exclude[j], i.c_str()) == 0)
                {
                    break;
                }
            }
            if (j == nexclude)
            {
                char *ptr;
                char *dup = strdup(i.c_str());
                while (NULL != (ptr = strchr(dup, '_')))
                {
                    *ptr = ' ';
                }
                mpt.AddCategory(dup);
                sfree(dup);
            }
        }
    }

    // get bondorders.
    //mol.PerceiveBondOrders();
    //mol.ConnectTheDots();

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
    alexandria::Experiment ca(program, method, basis, reference,
                              conformation, g98ptr);
    mpt.AddExperiment(ca);
    printf("Charge finally is %d\n", mol.GetTotalCharge());
    mpt.SetCharge(mol.GetTotalCharge());
    mpt.SetMass(mol.GetMolWt());
    mpt.SetMultiplicity(mol.GetTotalSpinMultiplicity());
    mpt.SetFormula(mol.GetFormula());

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
        double              temperature, DeltaHf0, DeltaHfT, DeltaGfT, DeltaSfT, S0T, CVT, CPT;
        std::vector<double> Scomponents;
        if (extract_thermochemistry(mol, false, &nsymm,
                                    0, 0.0,
                                    &temperature,
                                    &DeltaHf0,
                                    &DeltaHfT,
                                    &DeltaGfT,
                                    &DeltaSfT,
                                    &S0T,
                                    &CVT,
                                    &CPT,
                                    Scomponents))
        {
            alexandria::MolecularEnergy me1("DeltaHform",
                                            mpo_unit[MPO_ENERGY],
                                            0,
                                            epGAS,
                                            convert2gmx(DeltaHf0, eg2cKcal_Mole),
                                            0);
            mpt.LastExperiment()->AddEnergy(me1);
            alexandria::MolecularEnergy me2("DeltaHform",
                                            mpo_unit[MPO_ENERGY],
                                            temperature,
                                            epGAS,
                                            convert2gmx(DeltaHfT, eg2cKcal_Mole),
                                            0);
            mpt.LastExperiment()->AddEnergy(me2);
            alexandria::MolecularEnergy me3("DeltaGform",
                                            mpo_unit[MPO_ENERGY],
                                            temperature,
                                            epGAS,
                                            convert2gmx(DeltaGfT, eg2cKcal_Mole),
                                            0);
            mpt.LastExperiment()->AddEnergy(me3);
            alexandria::MolecularEnergy me4("DeltaSform",
                                            mpo_unit[MPO_ENTROPY],
                                            temperature,
                                            epGAS,
                                            convert2gmx(DeltaSfT, eg2cCal_MolK),
                                            0);
            mpt.LastExperiment()->AddEnergy(me4);
            alexandria::MolecularEnergy me5("S0",
                                            mpo_unit[MPO_ENTROPY],
                                            temperature,
                                            epGAS,
                                            convert2gmx(S0T, eg2cCal_MolK),
                                            0);
            mpt.LastExperiment()->AddEnergy(me5);
            alexandria::MolecularEnergy me6("cp",
                                            mpo_unit[MPO_ENTROPY],
                                            temperature,
                                            epGAS,
                                            convert2gmx(CPT, eg2cCal_MolK),
                                            0);
            mpt.LastExperiment()->AddEnergy(me6);
            const char *scomp[3] = { "Strans", "Srot", "Svib" };
            for (int i = 0; (i < 3); i++)
            {
                alexandria::MolecularEnergy mes(scomp[i],
                                                mpo_unit[MPO_ENTROPY],
                                                temperature,
                                                epGAS,
                                                convert2gmx(Scomponents[i], eg2cCal_MolK),
                                                0);
                mpt.LastExperiment()->AddEnergy(mes);
            }
        }
    }

    // Get the energy as well.
    alexandria::MolecularEnergy mes("HF", mpo_unit[MPO_ENERGY], 0, epGAS,
                                    convert2gmx( mol.GetEnergy(), eg2cKcal_Mole), 0);
    mpt.LastExperiment()->AddEnergy(mes);
    
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
    if (ff && (ff->Setup(mol)))
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
            mpt.LastExperiment()->AddAtom(ca);
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
        mpt.LastExperiment()->AddDipole(dp);
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
        mpt.LastExperiment()->AddQuadrupole(mq);
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
        mpt.LastExperiment()->AddPolar(mdp);
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

    found = 0;

    for (g = _gapv.begin(); (g < _gapv.end()); g++)
    {
        if ((0 == strcasecmp(g->getElement().c_str(), element)) &&
            (0 == strcasecmp(g->getMethod().c_str(), method)) &&
            (0 == strcasecmp(g->getDesc().c_str(), desc)) &&
            bCheckTemperature(temp, g->getTemp()))
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
               char *molnm, char *iupac, char *conf,
               char *basis, int maxpot, int nsymm,
               const char *forcefield)
{
#if HAVE_LIBOPENBABEL2
    gmx_molprop_read_babel(g98, mp, molnm, iupac, conf, basis,
                           maxpot, nsymm, forcefield);
#else
    gmx_fatal(FARGS, "For reading Gaussian input you need to link to OpenBabel");
#endif
}
