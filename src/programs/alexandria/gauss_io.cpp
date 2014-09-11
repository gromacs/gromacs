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
    if (ea.GetV() < eb.GetV())
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

void translate_atomtypes(t_atoms *atoms, t_symtab *tab, const char *forcefield)
{
    //#define XHACK
#ifdef XHACK
    typedef struct {
        const char *ob, *ax;
    } t_xhack;
    static const t_xhack xh[] = {
        { "H", "HA" },
        { "HC", "HA" },
        { "H", "HA" },
        { "HO", "HP" },
        { "HS", "HP" },
        { "HB", "HP" },
        { "HN", "HP" },
        { "C3", "CTA" },
        { "C", "CA" },
        { "C2", "CTR" },
        { "Car", "CRA" },
        { "C1", "CDI" },
        { "Nr", "NPI2" },
        { "N3", "NTE" },
        { "Nox", "NTN" },
        { "N2", "NTR" },
        { "Npl", "NTE" },
        { "O3", "OTE" },
        { "O2", "OPI2" },
        { "F", "F" }
    };
    int                  i, j, nxh = sizeof(xh)/sizeof(xh[0]);
    for (i = 0; (i < atoms->nr); i++)
    {
        for (j = 0; (j < nxh); j++)
        {
            if (strcasecmp(xh[j].ob, *(atoms->atomtype[i])) == 0)
            {
                atoms->atomtype[i]  = put_symtab(tab, xh[j].ax);
                atoms->atomtypeB[i] = atoms->atomtype[i];
                break;
            }
        }
        if (j == nxh)
        {
            gmx_fatal(FARGS, "No atom type translation found for %s",
                      *(atoms->atomtype[i]));
        }
    }
#else
    OpenBabel::OBTypeTable obt;
    std::string            src, dst;
    int                    i;

    if (NULL == forcefield)
    {
        return;
    }
    obt.SetFromType("INT");
    if (obt.SetToType(forcefield))
    {
        for (i = 0; (i < atoms->nr); i++)
        {
            src.assign(*(atoms->atomtype[i]));
            obt.Translate(src, dst);
            atoms->atomtype[i]  = put_symtab(tab, dst.c_str());
            atoms->atomtypeB[i] = atoms->atomtype[i];
        }
    }
#endif
}

static void gmx_molprop_read_babel(const char *g98,
                                   alexandria::MolProp &mpt,
                                   char *molnm, char *iupac, char *conformation,
                                   char *basisset, int maxpot,
                                   const char *forcefield)
{
    /* Read a gaussian log file */
    OpenBabel::OBMol           mol, mol2;
    OpenBabel::OBAtomIterator  OBai;
    OpenBabel::OBBondIterator  OBbi;
    OpenBabel::OBConversion   *conv;
    //OpenBabel::OBAtom *OBa;
    OpenBabel::OBBond         *OBb;
    OpenBabel::OBPairData     *OBpd, *OBdhf;
    OpenBabel::OBVectorData   *dipole;
    OpenBabel::OBMatrixData   *quadrupole, *pol_tensor;
    OpenBabel::OBFreeGrid     *esp;
    OpenBabel::OBElementTable *OBet;
    std::string                formula, attr, value, inchi;

    std::vector<alexandria::ElectrostaticPotential> espv;

    const char *reference = "Spoel2013a", *unknown = "unknown";
    char       *program, *method, *basis, *charge_model, *ptr, *g98ptr;
    int         bondid;
    double      dval;
    int         k;
    const char *etypes[] = { "DHf(0K)", "DHf(298.15K)" };

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

    // Get bondorders.
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

    for (k = 0; (k < 2); k++)
    {
        OBdhf = (OpenBabel::OBPairData *)mol.GetData(etypes[k]);
        if (NULL != OBdhf)
        {
            value = OBdhf->GetValue();
            dval  = convert2gmx(atof(value.c_str()), eg2cKcal_Mole);

            alexandria::MolecularEnergy me(etypes[k], unit2string(eg2cKj_Mole), dval, 0);
            mpt.LastCalculation()->AddEnergy(me);
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
            alexandria::AtomicCharge aq(charge_model, "e", atom->GetPartialCharge());

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
        alexandria::MolecularDipPolar dp("electronic",
                                         unit2string(eg2cDebye),
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

        alexandria::MolecularDipPolar mdp("electronic",
                                          unit2string(eg2cAngstrom3),
                                          mm[0], mm[4], mm[8], alpha, 0);
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

int alexandria::GaussAtomProp::GetValue(const char *element,
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
        if ((0 == strcasecmp(g->GetElement().c_str(), element)) &&
            (0 == strcasecmp(g->GetMethod().c_str(), method)) &&
            (0 == strcasecmp(g->GetDesc().c_str(), desc)) &&
            (fabs(temp - g->GetTemp()) < ttol))
        {
            *value = g->GetValue();
            found  = 1;
            break;
        }
    }
    return found;
}

static int gmx_molprop_add_dhform(alexandria::MolProp &mpt,
                                  alexandria::GaussAtomProp &gap,
                                  const char *method,
                                  double eg34, double ezpe, double etherm,
                                  gmx_bool bVerbose)
{
    alexandria::MolecularCompositionIterator mci;
    alexandria::AtomNumIterator              ani;
    char   desc[128];
    int    natom, ntot;
    double vm, ve, vdh, vvm, vve, vvdh, ee, e0Hartree, eTHartree;

    if (mpt.NAtom() == 0)
    {
        return 0;
    }

    ntot = 0;
    vm   = ve = vdh = 0;
    sprintf(desc, "%s(0K)", method);
    mci = mpt.SearchMolecularComposition("bosque");
    for (ani = mci->BeginAtomNum(); (ani < mci->EndAtomNum()); ani++)
    {
        /* There are natom atoms with name atomname */
        if ((1 == gap.GetValue(ani->GetAtom().c_str(),
                               method, desc, 0, &vvm)) &&
            (1 == gap.GetValue(ani->GetAtom().c_str(),
                               (char *)"exp",
                               (char *)"DHf(0K)", 0, &vve)) &&
            (1 == gap.GetValue(ani->GetAtom().c_str(),
                               (char *)"exp",
                               (char *)"H(0K)-H(298.15K)", 298.15, &vvdh)))
        {
            natom = ani->GetNumber();
            vm   += natom*vvm;
            ve   += natom*vve;
            vdh  += natom*vvdh;
            ntot += natom;
        }
        else
        {
            return 0;
        }
    }
    /* Make sure units match! */
    e0Hartree = eg34+ve-vm;
    eTHartree = eg34-ezpe+etherm+ve-vm-vdh;

    ee = convert2gmx(e0Hartree, eg2cHartree);

    alexandria::MolecularEnergy me("DHf(0K)", unit2string(eg2cKj_Mole), ee, 0);
    mpt.LastCalculation()->AddEnergy(me);

    if (bVerbose)
    {
        printf("natom %d e0Hartree %f eTHartree %f eg34 %f ve %f vm %f vdh %f ezpe %f etherm %f \n",
               ntot, e0Hartree, eTHartree, eg34, ve, vm, vdh, ezpe, etherm);
    }

    ee = convert2gmx(eTHartree, eg2cHartree);

    alexandria::MolecularEnergy me2("DHf(298.15K)", unit2string(eg2cKj_Mole), ee, 0);
    mpt.LastCalculation()->AddEnergy(me2);

    return 1;
}

/* Read a line from a G03/G09 composite method (G3, G4, etc) record */
static gmx_bool gau_comp_meth_read_line(char *line, real *temp, real *pres)
{
    std::vector<std::string> ptr1 = split(line, '=');
    if ((ptr1.size() >= 2) && (ptr1[1].length() > 0))
    {
        std::vector<std::string> ptr2 = split(ptr1[1], ' ');
        if ((ptr2.size() >= 1) &&
            (ptr2[0].length() > 0) &&
            (ptr1[2].length() > 0))
        {
            *temp = atof(ptr2[0].c_str());
            *pres = atof(ptr1[2].c_str());
        }
        return TRUE;
    }
    return FALSE;
}

static gmx_bool read_polar(char *str, tensor T)
{
    std::vector<std::string> ptr = split(str, ' ');
    if (ptr.size() >= 7)
    {
        T[XX][XX] = atof(ptr[2].c_str());
        T[YY][XX] = T[XX][YY] = atof(ptr[3].c_str());
        T[YY][YY] = atof(ptr[4].c_str());
        T[XX][ZZ] = T[ZZ][XX] = atof(ptr[5].c_str());
        T[YY][ZZ] = T[ZZ][YY] = atof(ptr[6].c_str());
        T[ZZ][ZZ] = atof(ptr[7].c_str());
        return TRUE;
    }
    return FALSE;
}

static gmx_bool read_dipole(char *str, rvec mu)
{
    std::vector<std::string> ptr = split(str, ' ');
    if (ptr.size() >= 5)
    {
        mu[XX] = atof(ptr[1].c_str());
        mu[YY] = atof(ptr[3].c_str());
        mu[ZZ] = atof(ptr[5].c_str());
        return TRUE;
    }
    return FALSE;
}

static gmx_bool read_quad(char *str1, char *str2, tensor Q)
{
    gmx_bool bRes;
    rvec     q1, q2;

    clear_mat(Q);
    bRes = read_dipole(str1, q1);
    if (bRes)
    {
        bRes = read_dipole(str2, q2);
        if (bRes)
        {
            Q[XX][XX] = q1[XX];
            Q[YY][YY] = q1[YY];
            Q[ZZ][ZZ] = q1[ZZ];
            Q[XX][YY] = Q[YY][XX] = q2[XX];
            Q[XX][ZZ] = Q[ZZ][XX] = q2[YY];
            Q[YY][ZZ] = Q[ZZ][YY] = q2[ZZ];
        }
    }
    return bRes;
}

static void gmx_molprop_read_log(const char *fn,
                                 alexandria::MolProp &mpt,
                                 alexandria::GaussAtomProp &gap,
                                 gmx_atomprop_t aps, gmx_poldata_t pd,
                                 char *molnm, char *iupac, char *conformation,
                                 char *basisset,
                                 int maxpot, gmx_bool bVerbose)
{
    /* Read a gaussian log file */
    char   **strings = NULL;
    char     sbuf[STRLEN];
    int      nstrings, atomicnumber = -1;
    double   ee, x, y, z, V, myener, mass = 0;
    int      i, j, k, kk, anumber, natom, nesp, nelprop, charge, nfitpoints = -1;
    int      len;
    gmx_bool bWarnESP = FALSE;
    gmx_bool bAtomicCenter;
    real     mm;
    char    *atomname, *ginc, *hfener, *mp2ener, *g2ener, *g3ener, *g4ener, *cbsener;
    char    *reference = (char *)"This Work";
    char    *program   = NULL, *method = NULL, *basis = NULL;
    char    *mymeth;
    real     temp, pres, ezpe, ezpe2, etherm, etherm2, comp_0K, comp_energy, comp_enthalpy, comp_free_energy;
    gmx_bool bEtherm = FALSE, bEzpe = FALSE, bTemp = FALSE;
    gmx_bool bPolar, bQuad, bDipole;
    tensor   polar, quad;
    rvec     dipole;
    int      status;

    std::vector<alexandria::ElectrostaticPotential> espv;
    std::string xyz_unit(unit2string(eg2cPm));
    std::string V_unit(unit2string(eg2cHartree_e));

    nstrings = get_file(fn, &strings);

    natom   = 0;
    nesp    = 0;
    nelprop = NOTSET;
    charge  = NOTSET;
    hfener  = NULL;
    mp2ener = NULL;
    g2ener  = NULL;
    g3ener  = NULL;
    g4ener  = NULL;
    cbsener = NULL;
    bPolar  = bQuad = bDipole = FALSE;

    /* First loop over strings to deduct basic stuff */
    for (i = 0; (i < nstrings); i++)
    {
        if (NULL != strstr(strings[i], "Revision"))
        {
            program = strdup(strings[i]);
            trim(program);
            len = strlen(program);
            if ((len > 0) && (program[len-1] == ','))
            {
                program[len-1] = '\0';
            }
        }
        else if ((NULL != strstr(strings[i], "Standard basis:")) &&
                 (NULL == basisset))
        {
            std::vector<std::string> ptr = split(strings[i], ' ');
            if (ptr.size() >= 2)
            {
                basis = strdup(ptr[2].c_str());
            }
        }
        else if (NULL != strstr(strings[i], "Temperature="))
        {
            status = gau_comp_meth_read_line(strings[i], &temp, &pres);
            if (bVerbose)
            {
                printf("na gau_(): temp %f pres %f status = %d\n", temp, pres, status);
            }
            bTemp = TRUE;
        }
        else if (NULL != strstr(strings[i], "Exact polarizability"))
        {
            bPolar = read_polar(strings[i], polar);
        }
        else if (NULL != strstr(strings[i], "Dipole moment"))
        {
            bDipole = read_dipole(strings[i+1], dipole);
            i      += 1;
        }
        else if (NULL != strstr(strings[i], "Traceless Quadrupole moment"))
        {
            bQuad = read_quad(strings[i+1], strings[i+2], quad);
            i    += 2;
        }
        else if (NULL != strstr(strings[i], "E(ZPE)="))
        {
            status = gau_comp_meth_read_line(strings[i], &ezpe2, &etherm2);
            if (bVerbose)
            {
                printf("na gau_(): ezpe2 %f etherm2 %f status=%d\n",
                       ezpe2, etherm2, status);
            }
        }
        else if (NULL != strstr(strings[i], "Zero-point correction="))
        {
            std::vector<std::string> ptr = split(strings[i], ' ');
            if (ptr.size() >= 2)
            {
                ezpe  = atof(ptr[2].c_str());
                bEzpe = TRUE;
            }
            if (bVerbose)
            {
                printf("na gau_(): ezpe %f \n", ezpe);
            }
        }
        else if (NULL != strstr(strings[i], "Thermal correction to Enthalpy="))
        {
            std::vector<std::string> ptr = split(strings[i], ' ');
            if (ptr.size() >= 4)
            {
                etherm  = atof(ptr[4].c_str());
                bEtherm = TRUE;
            }
            if (bVerbose)
            {
                printf("na gau_(): etherm %f \n", etherm);
            }
        }
        else if ((NULL != strstr(strings[i], "G2(0 K)=")) ||
                 (NULL != strstr(strings[i], "G3(0 K)=")) ||
                 (NULL != strstr(strings[i], "G4(0 K)=")) ||
                 (NULL != strstr(strings[i], "CBS-QB3 (0 K)=")))
        {
            status = gau_comp_meth_read_line(strings[i], &comp_0K, &comp_energy);
            if (bVerbose)
            {
                printf("na gau_(): comp_0K %f comp_energy %f status=%d\n",
                       comp_0K, comp_energy, status);
            }
        }
        else if ((NULL != strstr(strings[i], "G2 Enthalpy=")) ||
                 (NULL != strstr(strings[i], "G3 Enthalpy=")) ||
                 (NULL != strstr(strings[i], "G4 Enthalpy=")) ||
                 (NULL != strstr(strings[i], "CBS-QB3 Enthalpy=")))
        {
            status = gau_comp_meth_read_line(strings[i], &comp_enthalpy, &comp_free_energy);
            if (bVerbose)
            {
                printf("na gau_(): comp_enthalpy %f comp_free_energy %f status=%d\n",
                       comp_enthalpy, comp_free_energy, status);
            }
        }
        else if (NULL != strstr(strings[i], "GINC"))
        {
            j = 0;
            while ((i+j < nstrings) && (strlen(strings[i+j]) > 0))
            {
                j++;
            }
            snew(ginc, 80*(j+1));
            for (k = i; (k < i+j); k++)
            {
                trim(strings[k]);
                strncat(ginc, strings[k], 80);
            }
            i = k-1;
            std::vector<std::string> ptr = split(ginc, '\\');
            for (std::vector<std::string>::iterator j = ptr.begin(); (j < ptr.end()); ++j)
            {
                const char *pj = j->c_str();
                if (NULL != strstr(pj, "HF="))
                {
                    hfener = strdup(pj+3);
                }
                if (NULL != strstr(pj, "MP2="))
                {
                    mp2ener = strdup(pj+4);
                }
                if (NULL != strstr(pj, "G2="))
                {
                    g2ener = strdup(pj+3);
                }
                if (NULL != strstr(pj, "G3="))
                {
                    g3ener = strdup(pj+3);
                }
                if (NULL != strstr(pj, "G4="))
                {
                    g4ener = strdup(pj+3);
                }
                if (NULL != strstr(pj, "CBSQB3="))
                {
                    cbsener = strdup(pj+7);
                }
            }
        }
        else if (NULL != strstr(strings[i], "#P"))
        {
            if (NULL == method)
            {
                /* The first method is the true method */
                std::vector<std::string> ptr = split(strings[i], ' ');
                if (ptr.size() >= 1)
                {
                    std::vector<std::string> qtr = split(ptr[1], '/');
                    if (qtr.size() > 0)
                    {
                        method = strdup(qtr[0].c_str());
                    }
                }
            }
        }
    }

    /* Add the calculation */
    if ((NULL == basis) && (NULL != basisset))
    {
        basis = basisset;
    }
    if ((NULL != program) && (NULL != method) && (NULL != basis))
    {
        alexandria::Calculation calc(program, method, basis,
                                     reference, conformation, fn);
        if (bPolar)
        {
            alexandria::MolecularDipPolar mdi(unit2string(eg2cElectron),
                                              unit2string(eg2cAngstrom3),
                                              polar[XX][XX], polar[YY][YY], polar[ZZ][ZZ],
                                              trace(polar)/3.0, 0);
            calc.AddPolar(mdi);
        }
        if (bDipole)
        {
            alexandria::MolecularDipPolar mdi(unit2string(eg2cElectron),
                                              unit2string(eg2cDebye),
                                              dipole[XX], dipole[YY], dipole[ZZ],
                                              norm(dipole), 0);
            calc.AddDipole(mdi);
        }
        if (bQuad)
        {
            alexandria::MolecularQuadrupole mq(unit2string(eg2cElectron),
                                               unit2string(eg2cBuckingham),
                                               quad[XX][XX], quad[YY][YY], quad[ZZ][ZZ],
                                               quad[XX][YY], quad[XX][ZZ], quad[YY][ZZ]);
            calc.AddQuadrupole(mq);
        }

        for (i = 0; (i < nstrings); i++)
        {
            bAtomicCenter = (NULL != strstr(strings[i], "Atomic Center"));

            if (NULL != strstr(strings[i], "fitting atomic charges"))
            {
                if (1 == sscanf(strings[i], "%d", &kk))
                {
                    if (nfitpoints == -1)
                    {
                        nfitpoints = kk;
                    }
                    else
                    {
                        if (!bWarnESP)
                        {
                            /*   gmx_resp_warning(fn,i+1);*/
                            bWarnESP = TRUE;
                        }
                        if (kk != nfitpoints)
                        {
                            fprintf(stderr, "nfitpoints was %d, now is %d\n",
                                    nfitpoints, kk);
                        }
                        nfitpoints = kk;
                    }
                }
            }
            else if (NULL != strstr(strings[i], "Stoichiometry"))
            {
                if (1 == sscanf(strings[i], "%*s%s", sbuf))
                {
                    mpt.SetFormula(sbuf);
                    if (NULL == molnm)
                    {
                        mpt.SetMolname(sbuf);
                    }
                    else
                    {
                        mpt.SetMolname(molnm);
                    }
                    if (NULL == iupac)
                    {
                        mpt.SetIupac(sbuf);
                    }
                    else
                    {
                        mpt.SetIupac(iupac);
                    }
                }
            }
            else if (NULL != strstr(strings[i], "Coordinates (Angstroms)"))
            {
                atomicnumber = i;
            }
            else if ((atomicnumber >= 0) && (i >= atomicnumber+3))
            {
                if (6 == sscanf(strings[i], "%d%d%d%lf%lf%lf",
                                &k, &anumber, &kk, &x, &y, &z))
                {
                    if (natom+1 != k)
                    {
                        if (!bWarnESP)
                        {
                            /* gmx_resp_warning(fn,i+1);*/
                            bWarnESP = TRUE;
                        }
                    }
                    else
                    {
                        natom++;
                        alexandria::CalcAtom caa(atomname, atomname,
                                                 natom);
                        atomname = gmx_atomprop_element(aps, anumber);

                        if (TRUE == gmx_atomprop_query(aps, epropMass, "", atomname, &mm))
                        {
                            mass += mm;
                        }
                        caa.SetUnit(unit2string(eg2cPm));
                        caa.SetCoords(100*x, 100*y, 100*z);
                        calc.AddAtom(caa);
                    }
                }
                else
                {
                    atomicnumber = -1;
                }
            }
            else if (NULL != strstr(strings[i], "Charge ="))
            {
                if (1 != sscanf(strings[i], "%*s%*s%d", &charge))
                {
                    fprintf(stderr, "Can not read the charge on line %d of file %s\n", i+1, fn);
                }
            }
            else if (bAtomicCenter || (NULL != strstr(strings[i], "ESP Fit Center")))
            {
                if ((bAtomicCenter  && (4 != sscanf(strings[i], "%*s%*s%d%*s%*s%lf%lf%lf", &k, &x, &y, &z))) ||
                    (!bAtomicCenter && (4 != sscanf(strings[i], "%*s%*s%*s%d%*s%*s%lf%lf%lf", &k, &x, &y, &z))))
                {
                    fprintf(stderr, "Warning something fishy on line %d of file %s\n", i+1, fn);
                }
                else
                {
                    if (bAtomicCenter)
                    {
                        if (NULL != debug)
                        {
                            fprintf(debug, "Coordinates for atom %d found on line %d %8.3f  %8.3f  %8.3f\n", k, i, x, y, z);
                        }
                    }
                    alexandria::ElectrostaticPotential ep(xyz_unit, V_unit, ++nesp, 100*x, 100*y, 100*z, 0);
                    espv.push_back(ep);

                    if (NULL != debug)
                    {
                        fprintf(debug, "Coordinates for fit %d found on line %d\n", k, i);
                    }
                }
            }
            else if (NULL != strstr(strings[i], "Electrostatic Properties (Atomic Units)"))
            {
                if ((NOTSET != nelprop) && (!bWarnESP))
                {
                    /*gmx_resp_warning(fn,i+1);*/
                    bWarnESP = TRUE;
                }
                nelprop = i;
            }
            else if ((NOTSET != nelprop) && (i >= nelprop+6) && (i < nelprop+6+natom+nfitpoints))
            {
                if (2 != sscanf(strings[i], "%d%*s%lf", &k, &V))
                {
                    fprintf(stderr, "Warning something fishy on line %d of file %s\n", i+1, fn);
                }
                if ((k-1) != (i-nelprop-6))
                {
                    fprintf(stderr, "Warning something fishy with fit center number on line %d of file %s\n", i+1, fn);
                }
                if (k-1 >= nesp)
                {
                    fprintf(stderr, "More potential points (%d) than fit centers (%d). Check line %d of file %s\n",
                            k, nesp, i+1, fn);
                }
                espv[k-1].SetV(V);
                if (NULL != debug)
                {
                    fprintf(debug, "Potential %d found on line %d\n", k, i);
                }
            }
            sfree(strings[i]);
        }
        if (debug)
        {
            fprintf(debug, "Found %d atoms, %d esp points, and %d fitpoints\n", natom, nesp, nfitpoints);
        }
        sfree(strings);

        if ((charge == NOTSET) || (natom == 0))
        {
            gmx_fatal(FARGS, "Error reading Gaussian log file. Charge = %d, natom = %d",
                      charge, natom);
        }
        mpt.SetCharge(charge);
        mpt.SetMass(mass);

        /* Generate atomic composition, needed for energetics */
        (void) mpt.GenerateComposition(pd);

        /* Add energies */
        if ((NULL != g4ener) || (NULL != g2ener) || (NULL != g3ener) || (NULL != cbsener))
        {
            if (NULL != g4ener)
            {
                mymeth = (char *)"G4";
                myener = atof(g4ener);
            }
            else if (NULL != g3ener)
            {
                mymeth = (char *)"G3";
                myener = atof(g3ener);
            }
            else if (NULL != g2ener)
            {
                mymeth = (char *)"G2";
                myener = atof(g2ener);
            }
            else
            {
                mymeth = (char *)"CBS-QB3";
                myener = atof(cbsener);
            }
            if (bVerbose)
            {
                printf("myener %f\n", myener);
            }
            if (bEtherm && bEzpe && bTemp)
            {
                if (0 == gmx_molprop_add_dhform(mpt, gap, mymeth, myener, ezpe, etherm, bVerbose))
                {
                    fprintf(stderr, "No support for atomic energies in %s, method %s\n",
                            mpt.GetMolname().c_str(), mymeth);
                }
            }
        }
        else if (NULL != mp2ener)
        {
            ee = convert2gmx(atof(mp2ener), eg2cHartree);
            alexandria::MolecularEnergy me("MP2", unit2string(eg2cKj_Mole), ee, 0);
            calc.AddEnergy(me);
        }
        else if (NULL != hfener)
        {
            ee = convert2gmx(atof(hfener), eg2cHartree);
            alexandria::MolecularEnergy me("HF", unit2string(eg2cKj_Mole), ee, 0);
            calc.AddEnergy(me);
        }
        mpt.AddCalculation(calc);

        // Now add the electrostatic potential points
        merge_electrostatic_potential(mpt, espv, natom, maxpot);
    }
    else
    {
        fprintf(stderr, "Error reading %s, program = '%s' basis = '%s' method = '%s'\n",
                fn, program, basis, method);
        for (i = 0; (i < nstrings); i++)
        {
            sfree(strings[i]);
        }
        sfree(strings);
    }
}

void ReadGauss(const char *g98,
               alexandria::MolProp &mp,
               alexandria::GaussAtomProp &gap,
               gmx_bool bBabel,
               gmx_atomprop_t aps, gmx_poldata_t pd,
               char *molnm, char *iupac, char *conf,
               char *basis,
               int maxpot, gmx_bool bVerbose,
               const char *forcefield)
{
#ifdef HAVE_LIBOPENBABEL2
    if (bBabel)
    {
        gmx_molprop_read_babel(g98, mp, molnm, iupac, conf, basis,
                               maxpot, forcefield);
    }
    else
    {
        gmx_molprop_read_log(g98, mp, gap, aps, pd, molnm, iupac, conf, basis,
                             maxpot, bVerbose);
    }
#else
    gmx_molprop_read_log(g98, mp, gap, aps, pd, molnm, iupac, conf, basis,
                         maxpot, bVerbose);
#endif
}

#ifndef HAVE_LIBOPENBABEL2
void translate_atomtypes(t_atoms *atoms, t_symtab *tab, const char *forcefield)
{
    fprintf(stderr, "The library function translate_atomtypes works only when OpenBabel is installed.\n");
}
#endif
