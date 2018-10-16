/*
 * This source file is part of the Alexandria program.
 *
 * Copyright (C) 2014-2018
 *
 * Developers:
 *             Mohammad Mehdi Ghahremanpour,
 *             Paul J. van Maaren,
 *             David van der Spoel (Project leader)
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
 * Foundation, Inc., 51 Franklin Street, Fifth Floor,
 * Boston, MA  02110-1301, USA.
 */

/*! \internal \brief
 * Implements part of the alexandria program.
 * \author Mohammad Mehdi Ghahremanpour <mohammad.ghahremanpour@icm.uu.se>
 * \author David van der Spoel <david.vanderspoel@icm.uu.se>
 */

#include "gmxpre.h"
#include "babel_io.h"
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
#include "gromacs/utility/exceptions.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/futil.h"
#include "gromacs/utility/real.h"
#include "gromacs/utility/stringutil.h"

#include "molprop.h"
#include "molprop_util.h"
#include "poldata.h"
#include "stringutil.h"

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
#include <openbabel/elements.h>
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


static inline double A2PM(double a) {return a*1.0e+2;} /* Angstrom to pm */

static inline double NM_cubed_to_A_cubed(double a) {return a*1.0e+3;} /* nm^3 to A^3 */

enum einformat{
    einfGaussian    = 0,
    einfNotGaussian = 1,
    einfNR
};

BabelFile::BabelFile(BabelFileType      ftype,
                     const std::string &ext,
                     const std::string &InFormat)
    :
      ftype_(ftype),
      ext_(ext),
      InFormat_(InFormat)
{}

BabelFiles::BabelFiles ()
{
    bfiles_.push_back(BabelFile(ebftPDB,  ".pdb",  "pdb"));
    bfiles_.push_back(BabelFile(ebftXYZ,  ".xyz",  "xyz"));
    bfiles_.push_back(BabelFile(ebftSDF,  ".sdf",  "sdf"));
    bfiles_.push_back(BabelFile(ebftMOL,  ".mol",  "mol"));
    bfiles_.push_back(BabelFile(ebftMOL2, ".mol2", "mol2"));
    bfiles_.push_back(BabelFile(ebftG09,  ".log",  "g03"));
}

BabelFileIterator BabelFiles::findBabelFile(const std::string &fn)
{
    auto              extension = fn.substr(fn.find_last_of("."));
    BabelFileIterator fb        = bfiles_.begin(), fe = bfiles_.end();
    return std::find_if(fb, fe, [extension](const BabelFile &bf)
                        {
                            return (extension == bf.ext());
                        });
}

static void merge_electrostatic_potential(alexandria::MolProp                             &mpt,
                                          std::vector<alexandria::ElectrostaticPotential> &espv,
                                          int                                              natom,
                                          int                                              maxPotential)
{
    maxPotential = std::max(0, std::min(maxPotential, 100));
    int npot   = espv.size() - natom;
    int maxpot = (npot * maxPotential)/100;
    int mod    = npot / maxpot;
    int i      = 0;
    for (auto esi = espv.begin(); esi < espv.end(); esi++, i++)
    {
        if ((i < natom) || (((i-natom) % mod) == 0))
        {
            mpt.LastExperiment()->AddPotential(*esi);
        }
    }
}

static OpenBabel::OBConversion *readBabel(const       char *g09, 
                                          OpenBabel::OBMol *mol,
                                          einformat        *inputformat)
{
    std::ifstream g09f;
    bool          isGzip = false;
    if (!gmx_fexist(g09))
    {
        std::string g09z(g09);
        g09z += ".gz";
        g09f.open(g09z.c_str(), std::ios::in);
        isGzip = g09f.is_open();
    }
    else
    {
        g09f.open(g09, std::ios::in);
    }
    if (!g09f.is_open())
    {
        gmx_fatal(FARGS, "Cannot open file %s for reading", g09);
    }
    OpenBabel::OBConversion *conv       = new OpenBabel::OBConversion(&g09f, &std::cout); // Read from g09f
    auto                     babelfiles = BabelFiles();
    auto                     informat   = babelfiles.findBabelFile(g09)->informat().c_str();
    
    if (strcmp (informat, "g03") == 0 || strcmp (informat, "g09") == 0)
    {
        *inputformat = einfGaussian;
    }
    
    if (conv->SetInFormat(informat, isGzip))
    {
        bool read_ok = false;
        try
        {
            read_ok = conv->Read(mol, &g09f);
        }
        catch (const std::exception& ex)
        {        
            gmx::printFatalErrorMessage(stderr, ex);
        }
        
        if (read_ok)
        {
            g09f.close();
            return conv; // exit with success
        }
        else
        {
            fprintf(stderr, "Could not read input file %s with OpenBabel2.\n", g09);
        }
    }
    else
    {
        fprintf(stderr, "Input file %s has incomprehensible format.\n", g09);
    }
    g09f.close();
    return nullptr;
}

void readBabel(const char          *g09,
               alexandria::MolProp &mpt,
               const char          *molnm,
               const char          *iupac,
               const char          *conformation,
               const char          *basisset,
               int                  maxPotential,
               int                  nsymm,
               std::string          forcefield,
               const char          *jobType,
               double               qtot)
{
    std::string                formula;
    std::string                attr;
    std::string                value;
    einformat                  inputformat = einfNotGaussian;
    const char                *reference   = "Ghahremanpour2016a";
    const char                *mymol       = "AMM";
    const char                *myprogram   = "Alexandria-2018";
    const char                *mymethod    = "AFF";
    const char                *mybasis     = "ACM";
   
    
    /* Variables to read a Gaussian log file */
    char                      *QM_charge_model  = nullptr;
    char                      *g09ptr;
    OpenBabel::OBMol           mol;
    OpenBabel::OBAtomIterator  OBai;
    OpenBabel::OBBondIterator  OBbi;
    OpenBabel::OBBond         *OBb;
    OpenBabel::OBPairData     *OBpd;
    OpenBabel::OBPcharge      *OBpc;
    OpenBabel::OBVectorData   *dipole;
    OpenBabel::OBMatrixData   *quadrupole;
    OpenBabel::OBMatrixData   *pol_tensor;
    OpenBabel::OBFreeGrid     *esp;
    std::vector<alexandria::ElectrostaticPotential> espv;

    alexandria::jobType        jobtype = alexandria::string2jobType(jobType);    
    OpenBabel::OBConversion    *conv   = readBabel(g09, &mol, &inputformat);
    if (nullptr == conv)
    {
        fprintf(stderr, "Failed reading %s\n", g09);
        return;
    }
    delete conv;
    conv = new OpenBabel::OBConversion(&std::cin, &std::cout);

    // Chemical Categories
    if (conv->SetOutFormat("fpt"))
    {
        std::vector<std::string> excluded_categories = {
            ">", "C_ONS_bond", "Rotatable_bond", "Conjugated_double_bond", "Conjugated_triple_bond",
            "Chiral_center_specified", "Cis_double_bond", "Bridged_rings", "Conjugated_tripple_bond",
            "Trans_double_bond"
        };

        conv->AddOption("f", OpenBabel::OBConversion::OUTOPTIONS, "FP4");
        conv->AddOption("s");
        conv->Convert();
        OpenBabel::OBMol         mol_copy   = mol;  // We need a copy here because WriteString removes the H.
        std::vector<std::string> categories = gmx::splitString(conv->WriteString(&mol_copy, false));
        for (const auto &category : categories)
        {
            size_t j;
            for (j = 0; j < excluded_categories.size(); j++)
            {
                if (excluded_categories[j] == category)
                {
                    break;
                }
            }
            if (j == excluded_categories.size())
            {
                std::string dup = category;
                std::replace_if(dup.begin(), dup.end(), [](const char c) {return c == '_';}, ' ');
                mpt.AddCategory(dup);
            }
        }
    }

    // Basis Set
    std::string basis;
    OBpd = (OpenBabel::OBPairData *)mol.GetData("basis");
    if ((nullptr != basisset) && (strlen(basisset) > 0))
    {
        basis.assign(basisset);
    }
    else if (nullptr != OBpd)
    {
        basis = OBpd->GetValue();
        size_t p = basis.find(" (5D, 7F)");
        if (p != basis.npos)
        {
            basis.erase(p, basis.npos);
        }
    }
    else
    {
        basis.assign(mybasis);
    }

    // QM Program
    std::string program(myprogram);
    OBpd = (OpenBabel::OBPairData *)mol.GetData("program");
    if (nullptr != OBpd)
    {
        program.assign(OBpd->GetValue());
    }

    // Method
    std::string method(mymethod);
    OBpd = (OpenBabel::OBPairData *)mol.GetData("method");
    if (nullptr != OBpd)
    {
        method.assign(OBpd->GetValue());
    }
    g09ptr = (char *) strrchr(g09, '/');
    if (nullptr == g09ptr)
    {
        g09ptr = (char *)g09;
    }
    else
    {
        g09ptr++;
        if (strlen(g09ptr) == 0)
        {
            g09ptr = (char *)g09;
        }
    }

    alexandria::Experiment exp(program, method, basis, reference, conformation, g09ptr, jobtype);
    mpt.AddExperiment(exp);
    mpt.SetMass(mol.GetMolWt());
    mpt.SetMultiplicity(mol.GetTotalSpinMultiplicity());
    mpt.SetFormula(mol.GetFormula());
    
    if (inputformat == einfGaussian)
    {
        mpt.SetCharge(mol.GetTotalCharge());
        if (debug)
        {
            fprintf(debug, "The total charge of the molecule (%0.2f) is taken from %s\n", qtot, g09);
        }
    }
    else
    {
        mpt.SetCharge(qtot);
        if (debug)
        {
            fprintf(debug, "The total charge of the molecule (%0.2f) is assigned from the gentop command line\n", qtot);
        }
    }

    if (nullptr != molnm)
    {
        mpt.SetMolname(molnm);
    }
    else
    {
        mpt.SetMolname(mymol);
    }

    if (nullptr != iupac)
    {
        mpt.SetIupac(iupac);
    }
    else
    {
        mpt.SetIupac(mymol);
    }

    // Thermochemistry
    if (inputformat == einfGaussian)
    {
        double              temperature = 0, DeltaHf0 = 0, DeltaHfT = 0;
        double              DeltaGfT    = 0, DeltaSfT = 0, S0T      = 0;
        double              CVT         = 0, CPT      = 0, ZPE      = 0;
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
                                    Scomponents,
                                    &ZPE))
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
            alexandria::MolecularEnergy me7("ZPE",
                                            mpo_unit[MPO_ENERGY],
                                            0,
                                            epGAS,
                                            convert2gmx(ZPE, eg2cKcal_Mole),
                                            0);
            mpt.LastExperiment()->AddEnergy(me7);
        }
    }

    // HF Eenergy
    alexandria::MolecularEnergy mes("HF", mpo_unit[MPO_ENERGY], 0, epGAS, convert2gmx( mol.GetEnergy(), eg2cKcal_Mole), 0);
    mpt.LastExperiment()->AddEnergy(mes);

    // Atoms
    auto *ff = OpenBabel::OBForceField::FindForceField(forcefield);
    if (ff && (ff->Setup(mol)))
    {
        ff->GetAtomTypes(mol);
        FOR_ATOMS_OF_MOL (atom, mol)
        {
            OpenBabel::OBPairData *type = (OpenBabel::OBPairData*) atom->GetData("FFAtomType");
            if (nullptr == type)
            {
                gmx_fatal(FARGS, "Cannot find %s atom type for atom %s",
                          forcefield.c_str(), static_cast<int>(atom->GetIdx()));
            }
            if (nullptr != debug)
            {
                fprintf(debug, "atom %d gafftype %s OBtype %s\n", atom->GetIdx(), type->GetValue().c_str(), atom->GetType());
            }
            
            alexandria::CalcAtom ca(OpenBabel::OBElements::GetSymbol(atom->GetAtomicNum()), type->GetValue(), atom->GetIdx());
            
            ca.SetUnit(unit2string(eg2cPm));
            ca.SetCoords(A2PM(atom->x()), A2PM(atom->y()), A2PM(atom->z()));

            if (inputformat == einfGaussian)
            {   
                std::vector<std::string> QM_charge_models = {"Mulliken charges", "ESP charges", "Hirshfeld charges", "CM5 charges"};
                for (const auto &cs : QM_charge_models)
                {
                    OBpd = (OpenBabel::OBPairData *) mol.GetData(cs);
                    if (nullptr != OBpd)
                    {
                        QM_charge_model    = strdup(OBpd->GetValue().c_str());
                        OBpc               = (OpenBabel::OBPcharge *) mol.GetData(QM_charge_model);
                        auto PartialCharge = OBpc->GetPartialCharge();
                        alexandria::AtomicCharge aq(QM_charge_model, "e", 0.0, PartialCharge[atom->GetIdx()-1]);
                        ca.AddCharge(aq);
                    }
                }
                if (nullptr == QM_charge_model)
                {
                    printf("\n"
                           "WARNING: None of the QM charge models known to Alexandria found in %s\n"
                           "         Partial charge is assigned to zero for atom %d\n\n",
                           g09, atom->GetIdx());
                    alexandria::AtomicCharge aq("UnknownQMCharge", "e", 0.0, 0.0);
                    ca.AddCharge(aq);
                }
            }
            mpt.LastExperiment()->AddAtom(ca);
        }
    }
    else
    {
        gmx_fatal(FARGS, "Cannot read %s force field", forcefield.c_str());
    }

    // Bonds
    OBbi = mol.BeginBonds();    
    if (OBbi != mol.EndBonds())
    {
        for (OBb = mol.BeginBond(OBbi); (nullptr != OBb); OBb = mol.NextBond(OBbi))
        {
            alexandria::Bond ab(1+OBb->GetBeginAtom()->GetIndex(),
                                1+OBb->GetEndAtom()->GetIndex(),
                                OBb->GetBondOrder());
            mpt.AddBond(ab);
        }
    }
    else
    {
        gmx_fatal(FARGS, "No bond is found for %s\n", g09);
    }

    // Dipole
    dipole = (OpenBabel::OBVectorData *) mol.GetData("Dipole Moment");
    if (nullptr != dipole)
    {
        OpenBabel::vector3            v3 = dipole->GetData();
        alexandria::MolecularDipole   dp("electronic",
                                         unit2string(eg2cDebye),
                                         0.0,
                                         v3.GetX(), 
                                         v3.GetY(), 
                                         v3.GetZ(),
                                         v3.length(), 
                                         0.0);
        mpt.LastExperiment()->AddDipole(dp);
    }

    // Quadrupole
    quadrupole = (OpenBabel::OBMatrixData *) mol.GetData("Traceless Quadrupole Moment");
    if (nullptr != quadrupole)
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
    if (nullptr != pol_tensor)
    {
        OpenBabel::matrix3x3 m3 = pol_tensor->GetData();
        double               mm[9], alpha, fac;
        int                  i;
        m3.GetArray(mm);
        fac = NM_cubed_to_A_cubed(pow(convert2gmx(1, eg2cBohr), 3));
        for (i = 0; i < 9; i++)
        {
            mm[i] *= fac;
        }
        alpha = (mm[0]+mm[4]+mm[8])/3.0;

        alexandria::MolecularPolarizability mdp("electronic",
                                                unit2string(eg2cAngstrom3),
                                                0.0,
                                                mm[0], mm[4], mm[8],
                                                mm[1], mm[2], mm[5],
                                                alpha, 0);
        mpt.LastExperiment()->AddPolar(mdp);
    }

    // Electrostatic potential
    esp = (OpenBabel::OBFreeGrid *) mol.GetData("Electrostatic Potential");
    if (nullptr != esp && maxPotential > 0)
    {
        OpenBabel::OBFreeGridPoint        *fgp;
        OpenBabel::OBFreeGridPointIterator fgpi;
        std::string                        xyz_unit(unit2string(eg2cPm));
        std::string                        V_unit(unit2string(eg2cHartree_e));
        int                                espid = 0;

        fgpi = esp->BeginPoints();
        for (fgp = esp->BeginPoint(fgpi); (nullptr != fgp); fgp = esp->NextPoint(fgpi))
        {
            alexandria::ElectrostaticPotential ep(xyz_unit, V_unit, ++espid,
                                                  A2PM(fgp->GetX()),
                                                  A2PM(fgp->GetY()),
                                                  A2PM(fgp->GetZ()),
                                                  fgp->GetV());
            espv.push_back(ep);
        }
        merge_electrostatic_potential(mpt, espv, mol.NumAtoms(), maxPotential);
    }
}
#endif
