/*
 * This source file is part of the Aleandria project.
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
#include <string>
#include <vector>
#include "gromacs/math/utilities.h"
#include "gromacs/utility/smalloc.h"
#include "poldata.h"
#include "gmx_simple_comm.h"
#include "molprop.h"
#include "composition.h"
#include "stringutil.h"

const char *mpo_name[MPO_NR] =
{
    "potential", "dipole", "quadrupole", "polarizability", "energy"
};

const char *mpo_unit[MPO_NR] =
{
    "e/nm", "D", "B", "\\AA$^3$", "kJ/mol"
};

const char *cs_name(CommunicationStatus cs)
{
    switch (cs)
    {
        case CS_OK:
            static const char *ok = "Communication OK";
            return ok;
        case CS_ERROR:
            static const char *err = "Communication Error";
            return err;
        case CS_SEND_DATA:
            static const char *sd = "Communication sent data";
            return sd;
        case CS_RECV_DATA:
            static const char *rd = "Communication OK";
            return rd;
        default:
            gmx_fatal(FARGS, "Unknown communication status %d", (int) cs);
    }
    return NULL;
};

#define GMX_SEND_DATA 19823
#define GMX_SEND_DONE -666
static CommunicationStatus gmx_send_data(t_commrec *cr, int dest)
{
    gmx_send_int(cr, dest, GMX_SEND_DATA);

    return CS_OK;
}

static CommunicationStatus gmx_send_done(t_commrec *cr, int dest)
{
    gmx_send_int(cr, dest, GMX_SEND_DONE);

    return CS_OK;
}

static CommunicationStatus gmx_recv_data_(t_commrec *cr, int src, int line)
{
    int kk = gmx_recv_int(cr, src);

    if ((kk != GMX_SEND_DATA) && (kk != GMX_SEND_DONE))
    {
        gmx_fatal(FARGS, "Received %d in gmx_recv_data (line %d). Was expecting either %d or %d\n.", kk, line,
                  (int)GMX_SEND_DATA, (int)GMX_SEND_DONE);
    }
    return CS_OK;
}
#define gmx_recv_data(cr, src) gmx_recv_data_(cr, src, __LINE__)
#undef GMX_SEND_DATA
#undef GMX_SEND_DONE

namespace alexandria
{

void GenericProperty::SetType(std::string type)
{
    if ((_type.size() == 0) && (type.size() > 0))
    {
        _type = type;
    }
    else
    {
        if (_type.size() == 0)
        {
            fprintf(stderr, "Replacing GenericProperty type '%s' by '%s'\n", _type.c_str(), type.c_str());
        }
    }
}

void GenericProperty::SetUnit(std::string unit)
{
    if ((_unit.size() == 0) && (unit.size() > 0))
    {
        _unit = unit;
    }
    else
    {
        if (_unit.size() == 0)
        {
            fprintf(stderr, "Replacing GenericProperty unit '%s' by '%s'\n", _unit.c_str(), unit.c_str());
        }
    }
}

void CalcAtom::SetUnit(std::string unit)
{
    if ((_unit.size() == 0) && (unit.size() > 0))
    {
        _unit = unit;
    }
    else
    {
        if (_unit.size() == 0)
        {
            fprintf(stderr, "Replacing CalcAtom unit '%s' by '%s'\n", _unit.c_str(), unit.c_str());
        }
    }
}

void MolecularComposition::AddAtom(AtomNum an)
{
    AtomNumIterator mci = SearchAtom(an.GetAtom());
    if (mci == EndAtomNum())
    {
        _atomnum.push_back(an);
    }
    else
    {
        mci->SetNumber(mci->GetNumber()+an.GetNumber());
    }
}

int MolProp::NAtom()
{
    if (_mol_comp.size() > 0)
    {
        int nat = BeginMolecularComposition()->CountAtoms();
        return nat;
    }
    return 0;
}

void MolProp::AddBond(Bond b)
{
    BondIterator bi;
    bool         bFound = false;

    for (bi = BeginBond(); !bFound && (bi < EndBond()); bi++)
    {
        bFound = (((bi->GetAi() == b.GetAi()) && (bi->GetAj() == b.GetAj())) ||
                  ((bi->GetAi() == b.GetAj()) && (bi->GetAj() == b.GetAi())));
        if (bFound)
        {
            break;
        }
    }
    if (!bFound)
    {
        _bond.push_back(b);
    }
    else if ((NULL != debug) && (bi->GetBondOrder() != b.GetBondOrder()))
    {
        fprintf(debug, "Different bond orders in molecule %s\n", GetMolname().c_str());
        fflush(debug);
    }
}

void MolecularComposition::DeleteAtom(std::string catom)
{
    AtomNumIterator ani;

    if ((ani = SearchAtom(catom)) != EndAtomNum())
    {
        _atomnum.erase(ani);
    }
}

AtomNumIterator MolecularComposition::SearchAtom(std::string an)
{
    AtomNumIterator ani;

    for (ani = BeginAtomNum(); (ani < EndAtomNum()); ani++)
    {
        if (an.compare(ani->GetAtom()) == 0)
        {
            return ani;
        }
    }
    return EndAtomNum();
}

void MolecularComposition::ReplaceAtom(std::string oldatom, std::string newatom)
{
    AtomNumIterator i;

    for (i = BeginAtomNum(); (i < EndAtomNum()); i++)
    {
        if (oldatom.compare(i->GetAtom()) == 0)
        {
            i->SetAtom(newatom);
            break;
        }
    }
}

int MolecularComposition::CountAtoms(std::string atom)
{
    for (AtomNumIterator i = BeginAtomNum(); (i < EndAtomNum()); i++)
    {
        if (atom.compare(i->GetAtom()) == 0)
        {
            return i->GetNumber();
        }
    }
    return 0;
}

int MolecularComposition::CountAtoms(const char *atom)
{
    std::string str(atom);

    return CountAtoms(str);
}

int MolecularComposition::CountAtoms()
{
    int             nat = 0;
    AtomNumIterator i;

    for (i = BeginAtomNum(); (i < EndAtomNum()); i++)
    {
        nat += i->GetNumber();
    }
    return nat;
}

void MolProp::CheckConsistency()
{
}

bool MolProp::SearchCategory(std::string catname)
{
    std::vector<std::string>::iterator i;

    for (i = BeginCategory(); (i < EndCategory()); i++)
    {
        if (catname.compare(*i) == 0)
        {
            return true;
        }
    }
    return false;
}

void MolProp::DeleteComposition(std::string compname)
{
    MolecularCompositionIterator i;

    for (i = BeginMolecularComposition(); (i < EndMolecularComposition()); i++)
    {
        if (compname.compare(i->GetCompName()) == 0)
        {
            break;
        }
    }
    if (i < EndMolecularComposition())
    {
        _mol_comp.erase(i);
    }
}

void Experiment::Dump(FILE *fp)
{
    if (NULL != fp)
    {
        fprintf(fp, "Experiment %d polar %d dipole\n",
                NPolar(), NDipole());
        fprintf(fp, "reference    = %s\n", _reference.c_str());
        fprintf(fp, "conformation = %s\n", _conformation.c_str());
    }
}

void Calculation::Dump(FILE *fp)
{
    Experiment::Dump(fp);
    if (NULL != fp)
    {
        fprintf(fp, "Calculation %d polar %d dipole\n",
                NPolar(), NDipole());
        fprintf(fp, "program    = %s\n", _program.c_str());
        fprintf(fp, "method     = %s\n", _method.c_str());
        fprintf(fp, "basisset   = %s\n", _basisset.c_str());
        fprintf(fp, "datafile   = %s\n", _datafile.c_str());
    }
}

void Experiment::Merge(Experiment &src)
{
    MergeLow(&src);
}

void Experiment::MergeLow(Experiment *src)
{
    alexandria::MolecularEnergyIterator     mei;
    alexandria::MolecularDipPolarIterator   dpi;
    alexandria::MolecularQuadrupoleIterator mqi;

    for (mei = src->BeginEnergy(); (mei < src->EndEnergy()); mei++)
    {
        alexandria::MolecularEnergy me(mei->GetType(), mei->GetUnit(), mei->GetValue(), mei->GetError());
        AddEnergy(me);
    }

    for (dpi = src->BeginDipole(); (dpi < src->EndDipole()); dpi++)
    {
        alexandria::MolecularDipPolar dp(dpi->GetType(), dpi->GetUnit(),
                                         dpi->GetX(), dpi->GetY(), dpi->GetZ(),
                                         dpi->GetAver(), dpi->GetError());
        AddDipole(dp);
    }

    for (dpi = src->BeginPolar(); (dpi < src->EndPolar()); dpi++)
    {
        alexandria::MolecularDipPolar dp(dpi->GetType(), dpi->GetUnit(),
                                         dpi->GetX(), dpi->GetY(), dpi->GetZ(),
                                         dpi->GetAver(), dpi->GetError());
        AddPolar(dp);
    }

    for (mqi = src->BeginQuadrupole(); (mqi < src->EndQuadrupole()); mqi++)
    {
        alexandria::MolecularQuadrupole mq(mqi->GetType(), mqi->GetUnit(),
                                           mqi->GetXX(), mqi->GetYY(), mqi->GetZZ(),
                                           mqi->GetXY(), mqi->GetXZ(), mqi->GetYZ());
        AddQuadrupole(mq);
    }
}

void Calculation::Merge(Calculation &src)
{
    Experiment::MergeLow(&src);

    for (CalcAtomIterator cai = src.BeginAtom(); (cai < src.EndAtom()); cai++)
    {
        double   x, y, z;
        CalcAtom caa(cai->GetName(), cai->GetObtype(), cai->GetAtomid());

        cai->GetCoords(&x, &y, &z);
        caa.SetCoords(x, y, z);
        caa.SetUnit(cai->GetUnit());

        for (AtomicChargeIterator aci = cai->BeginQ(); (aci < cai->EndQ()); aci++)
        {
            AtomicCharge aq(aci->GetType(), aci->GetUnit(), aci->GetQ());
            caa.AddCharge(aq);
        }
        AddAtom(caa);
    }

    for (ElectrostaticPotentialIterator mep = src.BeginPotential(); (mep < src.EndPotential()); mep++)
    {
        alexandria::ElectrostaticPotential ep(mep->GetXYZunit(), mep->GetVunit(), mep->GetEspid(),
                                              mep->GetX(), mep->GetY(), mep->GetZ(), mep->GetV());
        AddPotential(ep);
    }
}

void CalcAtom::AddCharge(AtomicCharge q)
{
    AtomicChargeIterator aci;

    for (aci = BeginQ(); (aci < EndQ()); aci++)
    {
        if ((aci->GetType().compare(q.GetType()) == 0) &&
            (aci->GetUnit().compare(q.GetUnit()) == 0) &&
            (aci->GetQ() == q.GetQ()))
        {
            break;
        }
    }
    if (aci == EndQ())
    {
        _q.push_back(q);
    }
}

bool CalcAtom::Equal(CalcAtom ca)
{
    return !((_name.compare(ca.GetName()) != 0) ||
             (_obtype.compare(ca.GetObtype()) != 0) ||
             (_x != ca.GetX()) ||
             (_y != ca.GetY()) ||
             (_z != ca.GetZ()) ||
             (_atomid != ca.GetAtomid()));
}

CalcAtomIterator Calculation::SearchAtom(CalcAtom ca)
{
    CalcAtomIterator cai;
    for (cai = BeginAtom(); (cai < EndAtom()); cai++)
    {
        if (cai->Equal(ca))
        {
            break;
        }
    }
    return cai;
}

void Calculation::AddAtom(CalcAtom ca)
{
    CalcAtomIterator cai = SearchAtom(ca);

    if (cai == EndAtom())
    {
        _catom.push_back(ca);
    }
    else
    {
        printf("Trying to add identical atom %s (%s) twice. N = %d\n",
               ca.GetName().c_str(), ca.GetObtype().c_str(), (int)_catom.size());
    }
}

void MolProp::AddComposition(MolecularComposition mc)
{
    MolecularCompositionIterator mci = SearchMolecularComposition(mc.GetCompName());
    if (mci == EndMolecularComposition())
    {
        _mol_comp.push_back(mc);
    }
}

bool MolProp::BondExists(Bond b)
{
    for (alexandria::BondIterator bi = BeginBond(); (bi < EndBond()); bi++)
    {
        if (((bi->GetAi() == b.GetAi()) && (bi->GetAj() == b.GetAj())) ||
            ((bi->GetAi() == b.GetAj()) && (bi->GetAj() == b.GetAi())))
        {
            return true;
        }
    }
    return false;
}

void MolProp::Merge(MolProp &src)
{
    double      q, sq;
    std::string stmp, dtmp;

    for (std::vector<std::string>::iterator si = src.BeginCategory(); (si < src.EndCategory()); si++)
    {
        AddCategory(*si);
    }
    SetFormula(src.GetFormula());
    SetMass(src.GetMass());
    if (GetMultiplicity() <= 1)
    {
        SetMultiplicity(src.GetMultiplicity());
    }
    else
    {
        int smult = src.GetMultiplicity();
        if ((NULL != debug) && (smult != GetMultiplicity()))
        {
            fprintf(debug, "Not overriding multiplicity to %d when merging since it is %d (%s)\n",
                    smult, GetMultiplicity(), src.GetMolname().c_str());
            fflush(debug);
        }
    }
    q = GetCharge();
    if (q == 0)
    {
        SetCharge(src.GetCharge());
    }
    else
    {
        sq = src.GetCharge();
        if ((NULL != debug) && (sq != q))
        {
            fprintf(debug, "Not overriding charge to %g when merging since it is %g (%s)\n",
                    sq, q, GetMolname().c_str());
            fflush(debug);
        }
    }

    stmp = src.GetMolname();
    if ((GetMolname().size() == 0) && (stmp.size() != 0))
    {
        SetMolname(stmp);
    }
    stmp = src.GetIupac();
    if ((GetIupac().size() == 0) && (stmp.size() != 0))
    {
        SetIupac(stmp);
    }
    stmp = src.GetCas();
    if ((GetCas().size() == 0) && (stmp.size() != 0))
    {
        SetCas(stmp);
    }
    stmp = src.GetCid();
    if ((GetCid().size() == 0) && (stmp.size() != 0))
    {
        SetCid(stmp);
    }
    stmp = src.GetInchi();
    if ((GetInchi().size() == 0) && (stmp.size() != 0))
    {
        SetInchi(stmp);
    }
    if (NBond() == 0)
    {
        for (alexandria::BondIterator bi = src.BeginBond(); (bi < src.EndBond()); bi++)
        {
            alexandria::Bond bb(bi->GetAi(), bi->GetAj(), bi->GetBondOrder());
            AddBond(bb);
        }
    }
    else
    {
        for (alexandria::BondIterator bi = src.BeginBond(); (bi < src.EndBond()); bi++)
        {
            alexandria::Bond bb(bi->GetAi(), bi->GetAj(), bi->GetBondOrder());
            if (!BondExists(bb))
            {
                fprintf(stderr, "WARNING bond %d-%d not present in %s\n",
                        bi->GetAi(), bi->GetAj(), GetMolname().c_str());
            }
        }
    }

    for (alexandria::ExperimentIterator ei = src.BeginExperiment(); (ei < src.EndExperiment()); ei++)
    {
        Experiment ex(ei->GetReference(), ei->GetConformation());

        ex.Merge(*ei);
        AddExperiment(ex);
    }

    for (alexandria::CalculationIterator ci = src.BeginCalculation(); (ci < src.EndCalculation()); ci++)
    {
        Calculation ca(ci->GetProgram(), ci->GetMethod(),
                       ci->GetBasisset(), ci->GetReference(),
                       ci->GetConformation(), ci->GetDatafile());
        ca.Merge(*ci);

        AddCalculation(ca);
    }

    for (alexandria::MolecularCompositionIterator mci = src.BeginMolecularComposition();
         (mci < src.EndMolecularComposition()); mci++)
    {
        alexandria::MolecularComposition mc(mci->GetCompName());

        for (alexandria::AtomNumIterator ani = mci->BeginAtomNum(); (ani < mci->EndAtomNum()); ani++)
        {
            AtomNum an(ani->GetAtom(), ani->GetNumber());
            mc.AddAtom(an);
        }
        AddComposition(mc);
    }
}

MolecularCompositionIterator MolProp::SearchMolecularComposition(std::string str)
{
    MolecularCompositionIterator i;

    for (i = BeginMolecularComposition(); (i < EndMolecularComposition()); i++)
    {
        std::string s = i->GetCompName();
        if (s.compare(str) == 0)
        {
            break;
        }
    }
    return i;
}

void MolProp::Dump(FILE *fp)
{
    std::vector<std::string>::iterator si;
    ExperimentIterator                 ei;
    CalculationIterator                ci;

    if (fp)
    {
        fprintf(fp, "formula:      %s\n", GetFormula().c_str());
        fprintf(fp, "molname:      %s\n", GetMolname().c_str());
        fprintf(fp, "iupac:        %s\n", GetIupac().c_str());
        fprintf(fp, "CAS:          %s\n", GetCas().c_str());
        fprintf(fp, "cis:          %s\n", GetCid().c_str());
        fprintf(fp, "InChi:        %s\n", GetInchi().c_str());
        fprintf(fp, "mass:         %g\n", GetMass());
        fprintf(fp, "charge:       %d\n", GetCharge());
        fprintf(fp, "multiplicity: %d\n", GetMultiplicity());
        fprintf(fp, "category:    ");
        for (si = BeginCategory(); (si < EndCategory()); si++)
        {
            fprintf(fp, " '%s'", si->c_str());
        }
        fprintf(fp, "\n");
        for (ei = BeginExperiment(); (ei < EndExperiment()); ei++)
        {
            ei->Dump(fp);
        }
        for (ci = BeginCalculation(); (ci < EndCalculation()); ci++)
        {
            ci->Dump(fp);
        }
    }
}

bool MolProp::GenerateComposition(gmx_poldata_t pd)
{
    CalculationIterator  ci;
    CalcAtomIterator     cai;
    CompositionSpecs     cs;
    MolecularComposition mci_bosque(cs.searchCS(iCbosque)->name());
    MolecularComposition mci_alexandria(cs.searchCS(iCalexandria)->name());
    MolecularComposition mci_miller(cs.searchCS(iCmiller)->name());

    // Why was this again?
    //DeleteComposition(alexandria);
    //DeleteComposition(bosque);
    //DeleteComposition(miller);

    int natoms = 0;
    for (ci = BeginCalculation(); (mci_alexandria.CountAtoms() <= 0) && (ci < EndCalculation()); ci++)
    {
        /* This assumes we have either all atoms or none.
         * A consistency check could be
         * to compare the number of atoms to the formula */
        int nat = 0;
        for (cai = ci->BeginAtom(); (cai < ci->EndAtom()); cai++)
        {
            nat++;
            AtomNum ans(cai->GetObtype(), 1);
            mci_alexandria.AddAtom(ans);

            const char *ptype = gmx_poldata_atype_to_ptype(pd, cai->GetObtype().c_str());
            if (NULL != ptype)
            {
                const char *bos_type = gmx_poldata_ptype_to_bosque(pd, ptype);
                if (NULL != bos_type)
                {
                    AtomNum anb(bos_type, 1);
                    mci_bosque.AddAtom(anb);
                }
                const char *mil_type = gmx_poldata_ptype_to_miller(pd, ptype);

                if (NULL != mil_type)
                {
                    AtomNum anm(mil_type, 1);
                    mci_miller.AddAtom(anm);
                }
            }
        }
        natoms = std::max(natoms, nat);
    }

    if (natoms == mci_bosque.CountAtoms())
    {
        AddComposition(mci_bosque);
    }
    if (natoms == mci_miller.CountAtoms())
    {
        AddComposition(mci_miller);
    }
    if (natoms == mci_alexandria.CountAtoms())
    {
        AddComposition(mci_alexandria);
        if (NULL != debug)
        {
            fprintf(debug, "LO_COMP: ");
            for (AtomNumIterator ani = mci_alexandria.BeginAtomNum(); (ani < mci_alexandria.EndAtomNum()); ani++)
            {
                fprintf(debug, " %s:%d", ani->GetAtom().c_str(), ani->GetNumber());
            }
            fprintf(debug, "\n");
            fflush(debug);
        }
        return true;
    }
    return false;
}

static void add_element_to_formula(const char *elem, int number, char *formula, char *texform)
{
    if (number > 0)
    {
        strcat(formula, elem);
        strcat(texform, elem);
        if (number > 1)
        {
            char cnumber[32];

            sprintf(cnumber, "%d", number);
            strcat(formula, cnumber);
            sprintf(cnumber, "$_{%d}$", number);
            strcat(texform, cnumber);
        }
    }
}

bool MolProp::GenerateFormula(gmx_atomprop_t ap)
{
    char  formula[1280], texform[2560];
    int  *ncomp;
    alexandria::MolecularCompositionIterator mci;

    snew(ncomp, 110);
    formula[0] = '\0';
    texform[0] = '\0';
    mci        = SearchMolecularComposition("bosque");
    if (mci != EndMolecularComposition())
    {
        for (alexandria::AtomNumIterator ani = mci->BeginAtomNum(); (ani < mci->EndAtomNum()); ani++)
        {
            int         cnumber, an;
            real        value;
            std::string catom = ani->GetAtom();
            cnumber = ani->GetNumber();
            if (gmx_atomprop_query(ap, epropElement, "???", catom.c_str(), &value))
            {
                an = gmx_nint(value);
                range_check(an, 0, 110);
                if (an > 0)
                {
                    ncomp[an] += cnumber;
                }
            }
        }
    }
    add_element_to_formula("C", ncomp[6], formula, texform);
    add_element_to_formula("H", ncomp[1], formula, texform);
    ncomp[6] = ncomp[1] = 0;

    for (int j = 109; (j >= 1); j--)
    {
        add_element_to_formula(gmx_atomprop_element(ap, j), ncomp[j], formula, texform);
    }
    std::string mform = GetFormula();
    if (strlen(formula) > 0)
    {
        if (debug)
        {
            if ((mform.size() > 0) && (strcasecmp(formula, mform.c_str()) != 0))
            {
                fprintf(debug, "Formula '%s' does match '%s' based on composition for %s.\n",
                        mform.c_str(), formula, GetMolname().c_str());
                fflush(debug);
            }
        }
        SetFormula(formula);
        SetTexFormula(texform);
    }
    else if ((mform.size() == 0) && debug)
    {
        fprintf(debug, "Empty composition and formula for %s\n",
                GetMolname().c_str());
        fflush(debug);
    }
    sfree(ncomp);

    return (strlen(formula) > 0);
}

bool MolProp::HasComposition(std::string composition)
{
    bool bComp  = false;
    MolecularCompositionIterator mci;

    if (composition.size() > 0)
    {
        for (mci = BeginMolecularComposition(); !bComp && (mci < EndMolecularComposition()); mci++)
        {
            if (mci->GetCompName().compare(composition) == 0)
            {
                bComp = true;
            }
        }
    }
    if (debug && !bComp)
    {
        fprintf(debug, "No composition %s for molecule %s\n", composition.c_str(),
                GetMolname().c_str());
        fflush(debug);
    }

    return bComp;
}

bool Experiment::GetVal(const char *type, MolPropObservable mpo,
                        double *value, double *error, double vec[3],
                        tensor quadrupole)
{
    bool   done = false;
    double x, y, z;

    switch (mpo)
    {
        case MPO_ENERGY:
            for (MolecularEnergyIterator mei = BeginEnergy(); !done && (mei < EndEnergy()); mei++)
            {
                if ((NULL == type) || (strcasecmp(mei->GetType().c_str(), type) == 0))
                {
                    mei->Get(value, error);
                    done = true;
                }
            }
            break;
        case MPO_DIPOLE:
            for (MolecularDipPolarIterator mdp = BeginDipole(); !done && (mdp < EndDipole()); mdp++)
            {
                if ((NULL == type) || (strcasecmp(mdp->GetType().c_str(), type) == 0))
                {
                    mdp->Get(&x, &y, &z, value, error);
                    vec[XX] = x;
                    vec[YY] = y;
                    vec[ZZ] = z;
                    done    = true;
                }
            }
            break;
        case MPO_POLARIZABILITY:
            for (MolecularDipPolarIterator mdp = BeginPolar(); !done && (mdp < EndPolar()); mdp++)
            {
                if ((NULL == type) || (strcasecmp(mdp->GetType().c_str(), type) == 0))
                {
                    mdp->Get(&x, &y, &z, value, error);
                    vec[XX] = x;
                    vec[YY] = y;
                    vec[ZZ] = z;
                    done    = true;
                }
            }
            break;
        case MPO_QUADRUPOLE:
            for (MolecularQuadrupoleIterator mqi = BeginQuadrupole(); !done && (mqi < EndQuadrupole()); mqi++)
            {
                if ((NULL == type) || (strcasecmp(mqi->GetType().c_str(), type) == 0))
                {
                    double xx, yy, zz, xy, xz, yz;
                    mqi->Get(&xx, &yy, &zz, &xy, &xz, &yz);
                    quadrupole[XX][XX] = xx;
                    quadrupole[XX][YY] = xy;
                    quadrupole[XX][ZZ] = xz;
                    quadrupole[YY][XX] = 0;
                    quadrupole[YY][YY] = yy;
                    quadrupole[YY][ZZ] = yz;
                    quadrupole[ZZ][XX] = 0;
                    quadrupole[ZZ][YY] = 0;
                    quadrupole[ZZ][ZZ] = zz;
                    done               = true;
                }
            }
            break;
        default:
            break;
    }
    return done;
}

bool MolProp::GetPropRef(MolPropObservable mpo, iqmType iQM, char *lot,
                         const char *conf, const char *type, double *value, double *error,
                         char **ref, char **mylot,
                         double vec[3], tensor quadrupole)
{
    alexandria::ExperimentIterator  ei;
    alexandria::CalculationIterator ci;
    std::string                     reference, method, name, basisset, program, conformation, expconf;
    bool done = false;

    if ((iQM == iqmExp) || (iQM == iqmBoth))
    {
        for (ei = BeginExperiment(); (ei < EndExperiment()); ei++)
        {
            reference = ei->GetReference();
            expconf   = ei->GetConformation();

            if ((NULL == conf) || (strcasecmp(conf, expconf.c_str()) == 0))
            {
                done = ei->GetVal(type, mpo, value, error, vec, quadrupole);
            }
            if (done)
            {
                break;
            }
        }
        if (!done)
        {
            if (NULL != ref)
            {
                *ref = strdup(reference.c_str());
            }
            if (NULL != mylot)
            {
                *mylot = strdup("Experiment");
            }
        }
    }

    if ((!done) && ((iQM == iqmBoth) || (iQM == iqmQM)))
    {
        if (NULL != lot)
        {
            ci = GetLotPropType(lot, mpo, type);
            if (ci != EndCalculation())
            {
                basisset     = ci->GetBasisset();
                method       = ci->GetMethod();
                reference    = ci->GetReference();
                conformation = ci->GetConformation();

                if ((NULL == conf) || (strcasecmp(conf, conformation.c_str()) == 0))
                {
                    done = ci->GetVal(type, mpo, value, error, vec, quadrupole);
                }
                if (done && (NULL != ref))
                {
                    *ref = strdup(reference.c_str());
                }
                if (done && (NULL != mylot))
                {
                    *mylot = strdup(lot);
                }
            }
        }
    }
    return done;
}

bool MolProp::GetProp(MolPropObservable mpo, iqmType iQM, char *lot,
                      char *conf, char *type, double *value, double *error)
{
    double myerror, vec[3];
    tensor quad;
    bool   bReturn;
    
    bReturn = GetPropRef(mpo, iQM, lot, conf, type, value, &myerror,
                         NULL, NULL, vec, quad);
    if (NULL != error)
    {
        *error = myerror;
    }      
    return bReturn;
}


CalculationIterator MolProp::GetLotPropType(const char       *lot,
                                            MolPropObservable mpo,
                                            const char       *type)
{
    CalculationIterator      ci;

    std::vector<std::string> ll = split(lot, '/');
    if ((ll[0].length() > 0) && (ll[1].length() > 0))
    {
        for (ci = BeginCalculation(); (ci < EndCalculation()); ci++)
        {
            if ((strcasecmp(ci->GetMethod().c_str(), ll[0].c_str()) == 0) &&
                (strcasecmp(ci->GetBasisset().c_str(), ll[1].c_str()) == 0))
            {
                bool done = false;
                switch (mpo)
                {
                    case MPO_POTENTIAL:
                        done = ci->NPotential() > 0;
                        break;
                    case MPO_DIPOLE:
                        for (MolecularDipPolarIterator mdp = ci->BeginDipole(); !done && (mdp < ci->EndDipole()); mdp++)
                        {
                            done =  ((NULL == type) ||
                                     (strcasecmp(type, mdp->GetType().c_str()) == 0));
                        }
                        break;
                    case MPO_QUADRUPOLE:
                        for (MolecularQuadrupoleIterator mdp = ci->BeginQuadrupole(); !done && (mdp < ci->EndQuadrupole()); mdp++)
                        {
                            done =  ((NULL == type) ||
                                     (strcasecmp(type, mdp->GetType().c_str()) == 0));
                        }
                        break;
                    case MPO_POLARIZABILITY:
                        for (MolecularDipPolarIterator mdp = ci->BeginPolar(); !done && (mdp < ci->EndPolar()); mdp++)
                        {
                            done =  ((NULL == type) ||
                                     (strcasecmp(type, mdp->GetType().c_str()) == 0));
                        }
                        break;
                    case MPO_ENERGY:
                        for (MolecularEnergyIterator mdp = ci->BeginEnergy(); !done && (mdp < ci->EndEnergy()); mdp++)
                        {
                            done =  ((NULL == type) ||
                                     (strcasecmp(type, mdp->GetType().c_str()) == 0));
                        }
                        break;
                    default:
                        break;
                }
                if (done)
                {
                    break;
                }
            }
        }
        return ci;
    }
    else
    {
        return EndCalculation();
    }
}

CalculationIterator MolProp::GetLot(const char *lot)
{
    CalculationIterator      ci;

    std::vector<std::string> ll = split(lot, '/');
    if ((ll[0].length() > 0) && (ll[1].length() > 0))
    {
        bool done = false;
        for (ci = BeginCalculation(); (!done) && (ci < EndCalculation()); ci++)
        {
            done = ((strcasecmp(ci->GetMethod().c_str(), ll[0].c_str()) == 0) &&
                    (strcasecmp(ci->GetBasisset().c_str(), ll[1].c_str()) == 0));
            if (done)
            {
                break;
            }
        }
        return ci;
    }
    else
    {
        return EndCalculation();
    }
}

CommunicationStatus GenericProperty::Send(t_commrec *cr, int dest)
{
    CommunicationStatus cs;

    cs = gmx_send_data(cr, dest);
    if (CS_OK == cs)
    {
        gmx_send_str(cr, dest, _type.c_str());
        gmx_send_str(cr, dest, _unit.c_str());
    }
    else if (NULL != debug)
    {
        fprintf(debug, "Trying to send GenericProperty, status %s\n", cs_name(cs));
        fflush(debug);
    }
    return cs;
}

CommunicationStatus GenericProperty::Receive(t_commrec *cr, int src)
{
    CommunicationStatus cs;

    cs = gmx_recv_data(cr, src);
    if (CS_OK == cs)
    {
        _type.assign(gmx_recv_str(cr, src));
        _unit.assign(gmx_recv_str(cr, src));
    }
    else if (NULL != debug)
    {
        fprintf(debug, "Trying to receive GenericProperty, status %s\n", cs_name(cs));
        fflush(debug);
    }
    return cs;
}

CommunicationStatus MolecularQuadrupole::Send(t_commrec *cr, int dest)
{
    CommunicationStatus cs;

    cs = GenericProperty::Send(cr, dest);
    if (CS_OK == cs)
    {
        cs = gmx_send_data(cr, dest);
    }
    if (CS_OK == cs)
    {
        gmx_send_double(cr, dest, _xx);
        gmx_send_double(cr, dest, _yy);
        gmx_send_double(cr, dest, _zz);
        gmx_send_double(cr, dest, _xy);
        gmx_send_double(cr, dest, _xz);
        gmx_send_double(cr, dest, _yz);
    }
    else if (NULL != debug)
    {
        fprintf(debug, "Trying to send MolecularQuadrupole, status %s\n", cs_name(cs));
        fflush(debug);
    }
    return cs;
}

CommunicationStatus MolecularQuadrupole::Receive(t_commrec *cr, int src)
{
    CommunicationStatus cs;

    cs = GenericProperty::Receive(cr, src);
    if (CS_OK == cs)
    {
        cs = gmx_recv_data(cr, src);
    }
    if (CS_OK == cs)
    {
        _xx    = gmx_recv_double(cr, src);
        _yy    = gmx_recv_double(cr, src);
        _zz    = gmx_recv_double(cr, src);
        _xy    = gmx_recv_double(cr, src);
        _xz    = gmx_recv_double(cr, src);
        _yz    = gmx_recv_double(cr, src);
    }
    else if (NULL != debug)
    {
        fprintf(debug, "Trying to received MolecularQuadrupole, status %s\n", cs_name(cs));
        fflush(debug);
    }
    return cs;
}

CommunicationStatus MolecularEnergy::Receive(t_commrec *cr, int src)
{
    CommunicationStatus cs;

    cs = GenericProperty::Receive(cr, src);
    if (CS_OK == cs)
    {
        cs = gmx_recv_data(cr, src);
    }
    if (CS_OK == cs)
    {
        _value = gmx_recv_double(cr, src);
        _error = gmx_recv_double(cr, src);
    }
    else if (NULL != debug)
    {
        fprintf(debug, "Trying to receive MolecularEnergy, status %s\n", cs_name(cs));
        fflush(debug);
    }
    return cs;
}

CommunicationStatus MolecularEnergy::Send(t_commrec *cr, int dest)
{
    CommunicationStatus cs;

    cs = GenericProperty::Send(cr, dest);
    if (CS_OK == cs)
    {
        cs = gmx_send_data(cr, dest);
    }
    if (CS_OK == cs)
    {
        gmx_send_double(cr, dest, _value);
        gmx_send_double(cr, dest, _error);
    }
    else if (NULL != debug)
    {
        fprintf(debug, "Trying to send MolecularEnergy, status %s\n", cs_name(cs));
        fflush(debug);
    }
    return cs;
}

CommunicationStatus Bond::Send(t_commrec *cr, int dest)
{
    CommunicationStatus cs;

    cs = gmx_send_data(cr, dest);
    if (CS_OK == cs)
    {
        gmx_send_int(cr, dest, _ai);
        gmx_send_int(cr, dest, _aj);
        gmx_send_int(cr, dest, _bondorder);
    }
    else if (NULL != debug)
    {
        fprintf(debug, "Trying to send Bond, status %s\n", cs_name(cs));
        fflush(debug);
    }
    return cs;
}

CommunicationStatus Bond::Receive(t_commrec *cr, int src)
{
    CommunicationStatus cs;

    cs = gmx_recv_data(cr, src);
    if (CS_OK == cs)
    {
        _ai        = gmx_recv_int(cr, src);
        _aj        = gmx_recv_int(cr, src);
        _bondorder = gmx_recv_int(cr, src);
    }
    else if (NULL != debug)
    {
        fprintf(debug, "Trying to receive Bond, status %s\n", cs_name(cs));
        fflush(debug);
    }
    return cs;
}

CommunicationStatus MolecularDipPolar::Send(t_commrec *cr, int dest)
{
    CommunicationStatus cs;

    cs = GenericProperty::Send(cr, dest);
    if (CS_OK == cs)
    {
        cs = gmx_send_data(cr, dest);
    }
    if (CS_OK == cs)
    {
        gmx_send_double(cr, dest, _x);
        gmx_send_double(cr, dest, _y);
        gmx_send_double(cr, dest, _z);
        gmx_send_double(cr, dest, _aver);
        gmx_send_double(cr, dest, _error);
    }
    else if (NULL != debug)
    {
        fprintf(debug, "Trying to send MolecularDipPolar, status %s\n", cs_name(cs));
        fflush(debug);
    }
    return cs;
}

CommunicationStatus MolecularDipPolar::Receive(t_commrec *cr, int src)
{
    CommunicationStatus cs;

    cs = GenericProperty::Receive(cr, src);
    if (CS_OK == cs)
    {
        cs = gmx_recv_data(cr, src);
    }
    if (CS_OK == cs)
    {
        _x     = gmx_recv_double(cr, src);
        _y     = gmx_recv_double(cr, src);
        _z     = gmx_recv_double(cr, src);
        _aver  = gmx_recv_double(cr, src);
        _error = gmx_recv_double(cr, src);
    }
    else if (NULL != debug)
    {
        fprintf(debug, "Trying to receive MolecularDipPolar, status %s\n", cs_name(cs));
        fflush(debug);
    }
    return cs;
}

CommunicationStatus Experiment::Receive(t_commrec *cr, int src)
{
    CommunicationStatus cs;

    cs = gmx_recv_data(cr, src);
    if (CS_OK == cs)
    {
        //! Receive literature reference
        _reference    = gmx_recv_str(cr, src);

        //! Receive conformation
        _conformation = gmx_recv_str(cr, src);

        //! Receive polarizabilities
        do
        {
            MolecularDipPolar mdp;

            cs = mdp.Receive(cr, src);
            if (CS_OK == cs)
            {
                AddPolar(mdp);
            }
        }
        while (CS_OK == cs);

        //! Receive dipoles
        do
        {
            MolecularDipPolar mdp;

            cs = mdp.Receive(cr, src);
            if (CS_OK == cs)
            {
                AddDipole(mdp);
            }
        }
        while (CS_OK == cs);

        //! Receive energies
        do
        {
            MolecularEnergy me;

            cs = me.Receive(cr, src);
            if  (CS_OK == cs)
            {
                AddEnergy(me);
            }
        }
        while (CS_OK == cs);
    }
    if ((CS_OK != cs) && (NULL != debug))
    {
        fprintf(debug, "Trying to receive Experiment, status %s\n", cs_name(cs));
        fflush(debug);
    }
    return cs;
}

CommunicationStatus Experiment::Send(t_commrec *cr, int dest)
{
    CommunicationStatus cs;

    cs = gmx_send_data(cr, dest);
    if (CS_OK == cs)
    {
        //! Send literature reference
        gmx_send_str(cr, dest, _reference.c_str());

        //! Send conformation
        gmx_send_str(cr, dest, _conformation.c_str());

        //! Send polarizabilities
        for (MolecularDipPolarIterator dpi = BeginPolar(); (CS_OK == cs) && (dpi < EndPolar()); dpi++)
        {
            cs = dpi->Send(cr, dest);
        }
        cs = gmx_send_done(cr, dest);

        //! Send dipoles
        for (MolecularDipPolarIterator dpi = BeginDipole(); (CS_OK == cs) && (dpi < EndDipole()); dpi++)
        {
            cs = dpi->Send(cr, dest);
        }
        cs = gmx_send_done(cr, dest);

        //! Send energies
        for (MolecularEnergyIterator mei = BeginEnergy(); (CS_OK == cs) && (mei < EndEnergy()); mei++)
        {
            cs = mei->Send(cr, dest);
        }
        cs = gmx_send_done(cr, dest);
    }
    if ((CS_OK != cs) && (NULL != debug))
    {
        fprintf(debug, "Trying to send Experiment, status %s\n", cs_name(cs));
        fflush(debug);
    }
    return cs;
}

CommunicationStatus ElectrostaticPotential::Receive(t_commrec *cr, int src)
{
    CommunicationStatus cs;

    cs = gmx_recv_data(cr, src);
    if (CS_OK == cs)
    {
        _xyz_unit.assign(gmx_recv_str(cr, src));
        _V_unit.assign(gmx_recv_str(cr, src));
        _espid = gmx_recv_int(cr, src);
        _x     = gmx_recv_double(cr, src);
        _y     = gmx_recv_double(cr, src);
        _z     = gmx_recv_double(cr, src);
        _V     = gmx_recv_double(cr, src);
    }
    else if (NULL != debug)
    {
        fprintf(debug, "Trying to receive ElectrostaticPotential, status %s\n", cs_name(cs));
        fflush(debug);
    }
    return cs;
}

CommunicationStatus ElectrostaticPotential::Send(t_commrec *cr, int dest)
{
    CommunicationStatus cs;

    cs = gmx_send_data(cr, dest);
    if (CS_OK == cs)
    {
        gmx_send_str(cr, dest, _xyz_unit.c_str());
        gmx_send_str(cr, dest, _V_unit.c_str());
        gmx_send_int(cr, dest, _espid);
        gmx_send_double(cr, dest, _x);
        gmx_send_double(cr, dest, _y);
        gmx_send_double(cr, dest, _z);
        gmx_send_double(cr, dest, _V);
    }
    else if (NULL != debug)
    {
        fprintf(debug, "Trying to send ElectrostaticPotential, status %s\n", cs_name(cs));
        fflush(debug);
    }
    return cs;
}

CommunicationStatus AtomicCharge::Receive(t_commrec *cr, int src)
{
    CommunicationStatus cs;

    cs = GenericProperty::Receive(cr, src);

    if (CS_OK == cs)
    {
        cs = gmx_recv_data(cr, src);
    }
    if (CS_OK == cs)
    {
        _q = gmx_recv_double(cr, src);
    }
    else if (NULL != debug)
    {
        fprintf(debug, "Trying to receive AtomicCharge, status %s\n", cs_name(cs));
        fflush(debug);
    }
    return cs;
}

CommunicationStatus AtomicCharge::Send(t_commrec *cr, int dest)
{
    CommunicationStatus cs;

    cs = GenericProperty::Send(cr, dest);

    if (CS_OK == cs)
    {
        cs = gmx_send_data(cr, dest);
    }
    if (CS_OK == cs)
    {
        gmx_send_double(cr, dest, _q);
    }
    else if (NULL != debug)
    {
        fprintf(debug, "Trying to send AtomicCharge, status %s\n", cs_name(cs));
        fflush(debug);
    }
    return cs;
}

CommunicationStatus CalcAtom::Receive(t_commrec *cr, int src)
{
    CommunicationStatus cs;

    cs = gmx_recv_data(cr, src);
    if (CS_OK == cs)
    {
        _name.assign(gmx_recv_str(cr, src));
        _obtype.assign(gmx_recv_str(cr, src));
        _atomid = gmx_recv_int(cr, src);
        _unit.assign(gmx_recv_str(cr, src));
        _x = gmx_recv_double(cr, src);
        _y = gmx_recv_double(cr, src);
        _z = gmx_recv_double(cr, src);

        do
        {
            AtomicCharge aq;

            cs = aq.Receive(cr, src);
            if (CS_OK == cs)
            {
                AddCharge(aq);
            }
        }
        while (CS_OK == cs);

    }
    if (NULL != debug)
    {
        fprintf(debug, "Received CalcAtom, status %s\n", cs_name(cs));
        fflush(debug);
    }
    return cs;
}

CommunicationStatus CalcAtom::Send(t_commrec *cr, int dest)
{
    CommunicationStatus  cs;
    AtomicChargeIterator qi;

    cs = gmx_send_data(cr, dest);
    if (CS_OK == cs)
    {
        gmx_send_str(cr, dest, _name.c_str());
        gmx_send_str(cr, dest, _obtype.c_str());
        gmx_send_int(cr, dest, _atomid);
        gmx_send_str(cr, dest, _unit.c_str());
        gmx_send_double(cr, dest, _x);
        gmx_send_double(cr, dest, _y);
        gmx_send_double(cr, dest, _z);

        for (qi = BeginQ(); (CS_OK == cs) && (qi < EndQ()); qi++)
        {
            cs = qi->Send(cr, dest);
        }
        cs = gmx_send_done(cr, dest);
    }
    if (NULL != debug)
    {
        fprintf(debug, "Sent CalcAtom, status %s\n", cs_name(cs));
        fflush(debug);
    }
    return cs;
}

CommunicationStatus Calculation::Receive(t_commrec *cr, int src)
{
    CommunicationStatus            cs;
    ElectrostaticPotentialIterator epi;
    CalcAtomIterator               cai;

    cs = Experiment::Receive(cr, src);

    if (CS_OK == cs)
    {
        cs = gmx_recv_data(cr, src);
    }
    if (CS_OK == cs)
    {
        _program.assign(gmx_recv_str(cr, src));
        _method.assign(gmx_recv_str(cr, src));
        _basisset.assign(gmx_recv_str(cr, src));
        _datafile.assign(gmx_recv_str(cr, src));

        do
        {
            ElectrostaticPotential ep;

            cs = ep.Receive(cr, src);
            if (CS_OK == cs)
            {
                AddPotential(ep);
            }
        }
        while (CS_OK == cs);

        do
        {
            CalcAtom ca;

            cs = ca.Receive(cr, src);
            if (CS_OK == cs)
            {
                AddAtom(ca);
            }
        }
        while (CS_OK == cs);
    }
    if (NULL != debug)
    {
        fprintf(debug, "Received Calculation, status %s\n", cs_name(cs));
        fflush(debug);
    }
    return cs;
}

CommunicationStatus Calculation::Send(t_commrec *cr, int dest)
{
    CommunicationStatus            cs;
    ElectrostaticPotentialIterator epi;
    CalcAtomIterator               cai;

    cs = Experiment::Send(cr, dest);

    if (CS_OK == cs)
    {
        cs = gmx_send_data(cr, dest);
    }
    if (CS_OK == cs)
    {
        gmx_send_str(cr, dest, _program.c_str());
        gmx_send_str(cr, dest, _method.c_str());
        gmx_send_str(cr, dest, _basisset.c_str());
        gmx_send_str(cr, dest, _datafile.c_str());

        for (epi = BeginPotential(); (CS_OK == cs) && (epi < EndPotential()); epi++)
        {
            cs = epi->Send(cr, dest);
        }
        cs = gmx_send_done(cr, dest);

        for (cai = BeginAtom(); (CS_OK == cs) && (cai < EndAtom()); cai++)
        {
            cs = cai->Send(cr, dest);
        }
        cs = gmx_send_done(cr, dest);
    }
    if (NULL != debug)
    {
        fprintf(debug, "Sent Calculation, status %s\n", cs_name(cs));
        fflush(debug);
    }
    return cs;
}

CommunicationStatus AtomNum::Send(t_commrec *cr, int dest)
{
    CommunicationStatus cs;

    cs = gmx_send_data(cr, dest);
    if (CS_OK == cs)
    {
        gmx_send_str(cr, dest, _catom.c_str());
        gmx_send_int(cr, dest, _cnumber);
        if (NULL != debug)
        {
            fprintf(debug, "Sent AtomNum %s %d, status %s\n",
                    _catom.c_str(), _cnumber, cs_name(cs));
            fflush(debug);
        }
    }
    return cs;
}

CommunicationStatus AtomNum::Receive(t_commrec *cr, int src)
{
    CommunicationStatus cs;

    cs = gmx_recv_data(cr, src);
    if (CS_OK == cs)
    {
        _catom.assign(gmx_recv_str(cr, src));
        _cnumber = gmx_recv_int(cr, src);
        if (NULL != debug)
        {
            fprintf(debug, "Received AtomNum %s %d, status %s\n",
                    _catom.c_str(), _cnumber, cs_name(cs));
            fflush(debug);
        }
    }
    return cs;
}

CommunicationStatus MolecularComposition::Send(t_commrec *cr, int dest)
{
    CommunicationStatus cs = gmx_send_data(cr, dest);

    if (CS_OK == cs)
    {
        gmx_send_str(cr, dest, _compname.c_str());
        for (AtomNumIterator ani = BeginAtomNum(); (CS_OK == cs) && (ani < EndAtomNum()); ani++)
        {
            cs = ani->Send(cr, dest);
        }
        cs = gmx_send_done(cr, dest);
        if (NULL != debug)
        {
            fprintf(debug, "Sent MolecularComposition %s, status %s\n",
                    _compname.c_str(), cs_name(cs));
            fflush(debug);
        }
    }
    return cs;
}

CommunicationStatus MolecularComposition::Receive(t_commrec *cr, int src)
{
    CommunicationStatus cs;

    cs = gmx_recv_data(cr, src);
    if (CS_OK == cs)
    {
        _compname.assign(gmx_recv_str(cr, src));
        CommunicationStatus cs2;
        do
        {
            AtomNum an;

            cs2 = an.Receive(cr, src);
            if (CS_OK == cs2)
            {
                AddAtom(an);
            }
        }
        while (CS_OK == cs2);
        if (NULL != debug)
        {
            fprintf(debug, "Received MolecularComposition %s, status %s\n",
                    _compname.c_str(), cs_name(cs));
            fflush(debug);
        }
    }
    return cs;
}

CommunicationStatus MolProp::Send(t_commrec *cr, int dest)
{
    CommunicationStatus                cs;
    BondIterator                       bi;
    MolecularCompositionIterator       mci;
    std::vector<std::string>::iterator si;
    ExperimentIterator                 ei;
    CalculationIterator                ci;

    /* Generic stuff */
    cs = gmx_send_data(cr, dest);
    if (CS_OK == cs)
    {
        gmx_send_double(cr, dest, _mass);
        gmx_send_int(cr, dest, _charge);
        gmx_send_int(cr, dest, _multiplicity);
        gmx_send_str(cr, dest, _formula.c_str());
        gmx_send_str(cr, dest, _molname.c_str());
        gmx_send_str(cr, dest, _iupac.c_str());
        gmx_send_str(cr, dest, _cas.c_str());
        gmx_send_str(cr, dest, _cid.c_str());
        gmx_send_str(cr, dest, _inchi.c_str());

        /* Bonds */
        for (bi = BeginBond(); (CS_OK == cs) && (bi < EndBond()); bi++)
        {
            cs = bi->Send(cr, dest);
        }
        cs = gmx_send_done(cr, dest);

        /* Composition */
        for (mci = BeginMolecularComposition(); (CS_OK == cs) && (mci < EndMolecularComposition()); mci++)
        {
            cs = mci->Send(cr, dest);
        }
        cs = gmx_send_done(cr, dest);

        /* Categories */
        for (si = BeginCategory(); (CS_OK == cs) && (si < EndCategory()); si++)
        {
            cs = gmx_send_data(cr, dest);
            if (CS_OK == cs)
            {
                gmx_send_str(cr, dest, si->c_str());
                if (NULL != debug)
                {
                    fprintf(debug, "Sent category %s\n", si->c_str());
                    fflush(debug);
                }
            }
        }
        cs = gmx_send_done(cr, dest);

        /* Experiments */
        for (ei = BeginExperiment(); (CS_OK == cs) && (ei < EndExperiment()); ei++)
        {
            cs = ei->Send(cr, dest);
        }
        cs = gmx_send_done(cr, dest);

        /* Calculations */
        for (ci = BeginCalculation(); (CS_OK == cs) && (ci < EndCalculation()); ci++)
        {
            ci->Send(cr, dest);
        }
        cs = gmx_send_done(cr, dest);
        if (NULL != debug)
        {
            fprintf(debug, "Sent MolProp %s, status %s\n",
                    GetMolname().c_str(), cs_name(cs));
            fflush(debug);
        }
    }
    return cs;
}

CommunicationStatus MolProp::Receive(t_commrec *cr, int src)
{
    CommunicationStatus cs;

    /* Generic stuff */
    cs = gmx_recv_data(cr, src);
    if (CS_OK == cs)
    {
        //! Receive mass and more
        _mass         = gmx_recv_double(cr, src);
        _charge       = gmx_recv_int(cr, src);
        _multiplicity = gmx_recv_int(cr, src);
        _formula.assign(gmx_recv_str(cr, src));
        _molname.assign(gmx_recv_str(cr, src));
        _iupac.assign(gmx_recv_str(cr, src));
        _cas.assign(gmx_recv_str(cr, src));
        _cid.assign(gmx_recv_str(cr, src));
        _inchi.assign(gmx_recv_str(cr, src));
        if (NULL != debug)
        {
            fprintf(debug, "Got molname %s\n", GetMolname().c_str());
        }
        //! Receive Bonds
        do
        {
            Bond b;

            cs = b.Receive(cr, src);
            if (CS_OK == cs)
            {
                AddBond(b);
            }
        }
        while (CS_OK == cs);

        //! Receive Compositions
        do
        {
            MolecularComposition mc;

            cs = mc.Receive(cr, src);
            if (CS_OK == cs)
            {
                AddComposition(mc);
            }
        }
        while (CS_OK == cs);

        //! Receive Categories
        do
        {
            cs = gmx_recv_data(cr, src);
            if (CS_OK == cs)
            {
                char *ptr = gmx_recv_str(cr, src);
                if (NULL != ptr)
                {
                    AddCategory(ptr);
                    if (NULL != debug)
                    {
                        fprintf(debug, "Received a category %s\n", ptr);
                        fflush(debug);
                    }
                    sfree(ptr);
                }
                else
                {
                    gmx_fatal(FARGS, "A category was promised but I got a NULL pointer");
                }
            }
        }
        while (CS_OK == cs);

        //! Receive Experiments
        do
        {
            Experiment ex;

            cs = ex.Receive(cr, src);
            if (CS_OK == cs)
            {
                AddExperiment(ex);
            }
        }
        while (CS_OK == cs);

        //! Receive Calculations
        do
        {
            Calculation cc;

            cs = cc.Receive(cr, src);
            if (CS_OK == cs)
            {
                AddCalculation(cc);
            }
        }
        while (CS_OK == cs);
        if (NULL != debug)
        {
            fprintf(debug, "Reveived %d calculations from %d for mol %s\n",
                    NCalculation(), src, GetMolname().c_str());
            fprintf(debug, "Received MolProp %s, status %s\n",
                    GetMolname().c_str(), cs_name(cs));
            fflush(debug);
        }
    }
    return cs;
}

std::string MolProp::GetTexFormula()
{
    if (_texform.size() > 0)
    {
        return _texform;
    }
    else
    {
        return _formula;
    }
}

}
