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
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <vector>
#include <algorithm>
#include "gromacs/utility/smalloc.h"
#include "gromacs/math/utilities.h"
#include "gromacs/utility/real.h"
#include "gromacs/utility/stringutil.h"
#include "gromacs/utility/futil.h"
#include "gromacs/math/units.h"
#include "gromacs/math/vec.h"
#include "gromacs/pbcutil/pbc.h"
#include "gromacs/fileio/xvgr.h"
#include "gromacs/topology/atomprop.h"
#include "gromacs/gmxpreprocess/gpp_atomtype.h"
#include "gromacs/gmxpreprocess/gpp_nextnb.h"
#include "gromacs/gmxpreprocess/toputil.h"
#include "gromacs/gmxpreprocess/gen_ad.h"
#include "gmx_simple_comm.h"
#include "poldata.h"
#include "poldata_xml.h"
#include "molselect.h"
#include "molprop.h"
#include "molprop_util.h"
#include "molprop_xml.h"
#include "gentop_vsite.h"
#include "gentop_core.h"
#include "composition.h"
#include "stringutil.h"

void generate_composition(std::vector<alexandria::MolProp> &mp, gmx_poldata_t pd)
{
    int nOK = 0;
    alexandria::MolPropIterator  mpi;
    alexandria::CompositionSpecs cs;

    for (mpi = mp.begin(); (mpi < mp.end()); mpi++)
    {
        for (alexandria::CompositionSpecIterator csi = cs.beginCS(); (csi < cs.endCS()); ++csi)
        {
            mpi->DeleteComposition(csi->name());
        }
        if (true == mpi->GenerateComposition(pd))
        {
            nOK++;
        }
        else if (debug)
        {
            fprintf(debug, "Failed to make composition for %s\n",
                    mpi->GetMolname().c_str());
        }
    }
    if (mp.size() > 1)
    {
        printf("Generated composition for %d out of %d molecules.\n",
               nOK, (int)mp.size());
    }
}

void generate_formula(std::vector<alexandria::MolProp> &mp, gmx_atomprop_t ap)
{
    for (alexandria::MolPropIterator mpi = mp.begin(); (mpi < mp.end()); mpi++)
    {
        mpi->GenerateFormula(ap);
    }
}

alexandria::MolProp atoms_2_molprop(char *molname, int natoms, char **smnames,
                                    gmx_atomprop_t ap)
{
    alexandria::MolProp              mp;
    alexandria::CompositionSpecs     cs;
    alexandria::MolecularComposition mci(cs.searchCS(alexandria::iCalexandria)->name());
    int i;

    mp.SetMolname(molname);

    for (i = 0; (i < natoms); i++)
    {
        alexandria::AtomNum an(smnames[i], 1);
        mci.AddAtom(an);
    }
    mp.AddComposition(mci);
    (void) mp.GenerateFormula(ap);

    return mp;
}

static bool molprop_2_atoms(alexandria::MolProp mp, gmx_atomprop_t ap,
                            t_symtab *tab, const char *lot,
                            t_atoms *atoms, const char *q_algorithm,
                            rvec **x)
{
    alexandria::CalculationIterator  ci;
    alexandria::CalcAtomIterator     cai;
    alexandria::AtomicChargeIterator qi;
    std::string molnm;

    int         myunit;
    double      q, xx, yy, zz;
    int         natom;
    bool        bDone = false;

    (*x)  = NULL;
    molnm = mp.GetMolname();

    ci = mp.GetLot(lot);
    if (ci < mp.EndCalculation())
    {
        natom = 0;
        init_t_atoms(atoms, mp.NAtom(), false);
        snew(*x, mp.NAtom());

        for (cai = ci->BeginAtom(); (cai < ci->EndAtom()); cai++)
        {
            myunit = string2unit((char *)cai->GetUnit().c_str());
            cai->GetCoords(&xx, &yy, &zz);
            (*x)[natom][XX] = convert2gmx(xx, myunit);
            (*x)[natom][YY] = convert2gmx(yy, myunit);
            (*x)[natom][ZZ] = convert2gmx(zz, myunit);

            q = 0;
            for (qi = cai->BeginQ(); (qi < cai->EndQ()); qi++)
            {
                if (strcmp(qi->GetType().c_str(), q_algorithm) == 0)
                {
                    myunit = string2unit((char *)qi->GetUnit().c_str());
                    q      = convert2gmx(qi->GetQ(), myunit);
                    break;
                }
            }

            t_atoms_set_resinfo(atoms, natom, tab, molnm.c_str(), 1, ' ', 1, ' ');
            atoms->atomname[natom]        = put_symtab(tab, cai->GetName().c_str());
            atoms->atom[natom].atomnumber = gmx_atomprop_atomnumber(ap, cai->GetName().c_str());
            strcpy(atoms->atom[natom].elem, gmx_atomprop_element(ap, atoms->atom[natom].atomnumber));
            atoms->atom[natom].q      = q;
            atoms->atom[natom].resind = 0;
            natom++;
        }
        atoms->nr   = natom;
        atoms->nres = 1;
        bDone       = ((natom == mp.NAtom()) && (natom > 0));
    }
    if (NULL != debug)
    {
        fprintf(debug, "Tried to convert %s to gromacs. LOT is %s. Natoms is %d\n",
                molnm.c_str(), lot, natom);
    }
    return bDone;
}

bool molprop_2_topology2(alexandria::MolProp mp, gmx_atomprop_t ap,
                         t_symtab *tab, const char *lot,
                         t_topology **mytop, const char *q_algorithm,
                         rvec **x, t_params plist[F_NRE],
                         int nexcl, t_excls **excls)
{
    alexandria::BondIterator bi;

    int                      natom, ftb;
    t_param                  b;
    t_nextnb                 nnb;
    t_restp                 *rtp;
    t_topology              *top;

    snew(rtp, 1);
    snew(*mytop, 1);
    top = *mytop;
    init_top(top);

    /* Get atoms */
    if (false == molprop_2_atoms(mp, ap, tab, lot, &(top->atoms), q_algorithm, x))
    {
        return false;
    }

    /* Store bonds in harmonic potential list first, update type later */
    ftb = F_BONDS;
    memset(&b, 0, sizeof(b));
    for (bi = mp.BeginBond(); (bi < mp.EndBond()); bi++)
    {
        b.a[0] = bi->GetAi() - 1;
        b.a[1] = bi->GetAj() - 1;
        add_param_to_list(&(plist[ftb]), &b);
    }
    natom = mp.NAtom();

    /* Make Angles and Dihedrals */
    snew((*excls), natom);
    init_nnb(&nnb, natom, nexcl+2);
    gen_nnb(&nnb, plist);
    //detect_rings(&plist[F_BONDS],atoms->nr,bRing);
    //nbonds = plist[F_BONDS].nr;
    print_nnb(&nnb, "NNB");
    rtp->bKeepAllGeneratedDihedrals    = true;
    rtp->bRemoveDihedralIfWithImproper = true;
    rtp->bGenerateHH14Interactions     = true;
    rtp->nrexcl = nexcl;
    generate_excls(&nnb, nexcl, *excls);
    gen_pad(&nnb, &(top->atoms), rtp, plist, *excls, NULL, false);
    done_nnb(&nnb);
    sfree(rtp);

    return true;
}

int my_strcmp(const char *a, const char *b)
{
    if ((NULL != a) && (NULL != b))
    {
        return strcasecmp(a, b);
    }
    else
    {
        return 1;
    }
}

void merge_doubles(std::vector<alexandria::MolProp> &mp, char *doubles,
                   bool bForceMerge)
{
    alexandria::MolPropIterator mpi, mmm[2];
    std::string                 molname[2];
    std::string                 form[2];

    FILE *fp;
    int   i, ndouble = 0;
    bool  bForm, bName, bDouble;
    int   cur = 0;
#define prev (1-cur)

    if (NULL != doubles)
    {
        fp = fopen(doubles, "w");
    }
    else
    {
        fp = NULL;
    }
    i = 0;
    for (mpi = mp.begin(); (mpi < mp.end()); )
    {
        bDouble = false;
        mpi->Dump(debug);
        mmm[cur]     = mpi;
        molname[cur] = mpi->GetMolname();
        form[cur]    = mpi->GetFormula();
        if (i > 0)
        {
            bForm = (form[prev] == form[cur]);
            bName = (molname[prev] == molname[cur]);
            if (bName)
            {
                if (!bForm && debug)
                {
                    fprintf(debug, "%s %s with formulae %s - %s\n",
                            bForceMerge ? "Merging molecules" : "Found molecule",
                            molname[prev].c_str(), form[prev].c_str(), form[cur].c_str());
                }
                if (bForceMerge || bForm)
                {
                    if (fp)
                    {
                        fprintf(fp, "%5d  %s\n", ndouble+1, molname[prev].c_str());
                    }
                    mmm[prev]->Merge(*(mmm[cur]));
                    mpi = mp.erase(mmm[cur]);

                    bDouble = true;
                    ndouble++;
                    if (mpi != mp.end())
                    {
                        mmm[cur]      = mpi;
                        molname[cur]  = mmm[cur]->GetMolname();
                        form[cur]     = mmm[cur]->GetFormula();
                    }
                }
            }
        }
        if (!bDouble)
        {
            cur = prev;
            mpi++;
            i++;
        }
    }
    for (mpi = mp.begin(); (mpi < mp.end()); mpi++)
    {
        mpi->Dump(debug);
    }
    if (fp)
    {
        fclose(fp);
    }
    printf("There were %d double entries, leaving %d after merging.\n",
           ndouble, (int)mp.size());
}

static void dump_mp(std::vector<alexandria::MolProp> mp)
{
    alexandria::MolPropIterator mpi;
    FILE *fp;

    fp = fopen("dump_mp.dat", "w");

    for (mpi = mp.begin(); (mpi < mp.end()); mpi++)
    {
        fprintf(fp, "%-20s  %s\n", mpi->GetFormula().c_str(),
                mpi->GetMolname().c_str());
    }

    fclose(fp);
}

void merge_xml(int nfile, char **filens,
               std::vector<alexandria::MolProp> &mpout,
               char *outf, char *sorted, char *doubles,
               gmx_atomprop_t ap, gmx_poldata_t pd,
               bool bForceMerge)
{
    std::vector<alexandria::MolProp> mp;
    alexandria::MolPropIterator      mpi;
    int       i, npout = 0, tmp;

    for (i = 0; (i < nfile); i++)
    {
        if (!gmx_fexist(filens[i]))
        {
            continue;
        }
        MolPropRead(filens[i], mp);
        generate_composition(mp, pd);
        generate_formula(mp, ap);
        for (mpi = mp.begin(); (mpi < mp.end()); )
        {
            mpout.push_back(*mpi);
            mpi = mp.erase(mpi);
        }
    }
    tmp = mpout.size();
    if (debug)
    {
        fprintf(debug, "mpout.size() = %u mpout.max_size() = %u\n",
                (unsigned int)mpout.size(), (unsigned int)mpout.max_size());
    }
    for (mpi = mpout.begin(); (mpi < mpout.end()); mpi++)
    {
        mpi->Dump(debug);
    }
    MolPropSort(mpout, MPSA_MOLNAME, NULL, NULL);
    merge_doubles(mpout, doubles, bForceMerge);
    printf("There are %d total molecules before merging, %d after.\n",
           tmp, (int)mpout.size());

    if (outf)
    {
        printf("There are %d entries to store in output file %s\n", npout, outf);
        MolPropWrite(outf, mpout, false);
    }
    if (sorted)
    {
        MolPropSort(mpout, MPSA_FORMULA, NULL, NULL);
        MolPropWrite(sorted, mpout, false);
        dump_mp(mpout);
    }
}

static bool comp_mp_molname(alexandria::MolProp ma,
                            alexandria::MolProp mb)
{
    std::string mma = ma.GetMolname();
    std::string mmb = mb.GetMolname();

    return (mma.compare(mmb) < 0);
}

static bool comp_mp_formula(alexandria::MolProp ma,
                            alexandria::MolProp mb)
{
    std::string fma = ma.GetFormula();
    std::string fmb = mb.GetFormula();

    if (fma.compare(fmb) < 0)
    {
        return true;
    }
    else
    {
        return comp_mp_molname(ma, mb);
    }
}

gmx_atomprop_t my_aps;

static bool comp_mp_elem(alexandria::MolProp ma,
                         alexandria::MolProp mb)
{
    int         i;
    alexandria::MolecularCompositionIterator mcia, mcib;
    std::string bosque("bosque"), C("C");

    mcia = ma.SearchMolecularComposition(bosque);
    mcib = mb.SearchMolecularComposition(bosque);

    if ((mcia != ma.EndMolecularComposition()) &&
        (mcib != mb.EndMolecularComposition()))
    {
        int d = mcia->CountAtoms(C) - mcib->CountAtoms(C);

        if (d < 0)
        {
            return true;
        }
        else if (d > 0)
        {
            return false;
        }
        else
        {
            for (i = 1; (i <= 109); i++)
            {
                if (i != 6)
                {
                    const char *ee = gmx_atomprop_element(my_aps, i);
                    if (NULL != ee)
                    {
                        std::string elem(ee);
                        d = mcia->CountAtoms(elem) - mcib->CountAtoms(elem);
                        if (d == 0)
                        {
                            continue;
                        }
                        else
                        {
                            return (d < 0);
                        }
                    }
                }
            }
        }
    }
    return comp_mp_molname(ma, mb);
}

static bool comp_mp_index(alexandria::MolProp ma,
                          alexandria::MolProp mb)
{
    return (ma.GetIndex() < mb.GetIndex());
}

alexandria::MolPropIterator SearchMolProp(std::vector<alexandria::MolProp> &mp,
                                          const char                       *name)
{
    alexandria::MolPropIterator mpi;
    for (mpi = mp.begin(); (mpi < mp.end()); mpi++)
    {
        if (strcmp(mpi->GetMolname().c_str(), name) == 0)
        {
            break;
        }
    }
    return mpi;
}

void MolPropSort(std::vector<alexandria::MolProp> &mp,
                 MolPropSortAlgorithm mpsa, gmx_atomprop_t apt,
                 gmx_molselect_t gms)
{
    printf("There are %d molprops. Will now sort them.\n", (int)mp.size());
    for (alexandria::MolPropIterator mpi = mp.begin(); (mpi < mp.end()); mpi++)
    {
        if (NULL != debug)
        {
            mpi->Dump(debug);
        }
    }
    switch (mpsa)
    {
        case MPSA_MOLNAME:
            std::sort(mp.begin(), mp.end(), comp_mp_molname);
            break;
        case MPSA_FORMULA:
            std::sort(mp.begin(), mp.end(), comp_mp_formula);
            break;
        case MPSA_COMPOSITION:
            if (NULL != apt)
            {
                my_aps = apt;
                std::sort(mp.begin(), mp.end(), comp_mp_elem);
                my_aps = NULL;
            }
            else
            {
                gmx_fatal(FARGS, "Requesting a composition sort but atomprop is NULL");
            }
            break;
        case MPSA_SELECTION:
            if (NULL != gms)
            {
                for (alexandria::MolPropIterator mpi = mp.begin(); (mpi < mp.end()); mpi++)
                {
                    int index = gmx_molselect_index(gms, mpi->GetIupac().c_str());
                    mpi->SetIndex(index);
                }
                std::sort(mp.begin(), mp.end(), comp_mp_index);
            }
            else
            {
                gmx_fatal(FARGS, "Need molecule selection to sort on");
            }
            break;
        default:
            gmx_incons("Invalid algorithm for MolPropSort");
    }
}

static void add_qmc_conf(t_qmcount *qmc, const char *conformation)
{
    int j;

    for (j = 0; (j < qmc->nconf); j++)
    {
        if (strcasecmp(qmc->conf[j], conformation) == 0)
        {
            break;
        }
    }
    if (j == qmc->nconf)
    {
        qmc->nconf++;
        srenew(qmc->conf, qmc->nconf);
        qmc->conf[j] = strdup(conformation);
    }
}

static int get_qmc_count(t_qmcount *qmc, const char *method, const char *basis, const char *type)
{
    int j;

    for (j = 0; (j < qmc->n); j++)
    {
        if ((strcasecmp(qmc->method[j], method) == 0) &&
            (strcasecmp(qmc->basis[j], basis) == 0) &&
            (strcasecmp(qmc->type[j], type) == 0))
        {
            return qmc->count[j];
        }
    }
    return 0;
}

static void add_qmc_calc(t_qmcount *qmc, const char *method, const char *basis, const char *type)
{
    int j;

    for (j = 0; (j < qmc->n); j++)
    {
        if ((strcasecmp(qmc->method[j], method) == 0) &&
            (strcasecmp(qmc->basis[j], basis) == 0) &&
            (strcasecmp(qmc->type[j], type) == 0))
        {
            break;
        }
    }
    if (j == qmc->n)
    {
        srenew(qmc->method, qmc->n+1);
        srenew(qmc->basis, qmc->n+1);
        srenew(qmc->type, qmc->n+1);
        srenew(qmc->count, qmc->n+1);
        srenew(qmc->lot, qmc->n+1);
        qmc->method[qmc->n] = strdup(method);
        qmc->basis[qmc->n]  = strdup(basis);
        qmc->type[qmc->n]   = strdup(type);
        snew(qmc->lot[qmc->n], strlen(method)+strlen(basis)+strlen(type)+16);
        sprintf(qmc->lot[qmc->n], "%s/%s/%s", method, basis, type);
        qmc->count[qmc->n]  = 1;
        qmc->n++;
    }
    else
    {
        qmc->count[j]++;
    }
}

t_qmcount *find_calculations(std::vector<alexandria::MolProp> mp,
                             MolPropObservable                mpo,
                             const char                      *fc_str)
{
    alexandria::MolPropIterator        mpi;
    alexandria::ExperimentIterator     ei;
    alexandria::CalculationIterator    ci;

    const char                        *method, *basis;
    int                                i, n;
    double                             value, error, vec[3];
    tensor                             quadrupole;
    t_qmcount                         *qmc;
    std::vector<std::string>           types;
    std::vector<std::string>::iterator ti;

    snew(qmc, 1);
    for (mpi = mp.begin(); (mpi < mp.end()); mpi++)
    {
        for (ei = mpi->BeginExperiment(); (ei < mpi->EndExperiment()); ei++)
        {
            add_qmc_conf(qmc, ei->GetConformation().c_str());
        }
        for (ci = mpi->BeginCalculation(); (ci < mpi->EndCalculation()); ci++)
        {
            add_qmc_conf(qmc, ci->GetConformation().c_str());
        }
    }
    if (NULL != fc_str)
    {
        std::vector<std::string> qm = split(fc_str, ':');
        n = 0;
        for (std::vector<std::string>::iterator pqm = qm.begin(); (pqm < qm.end()); ++pqm)
        {
            if (pqm->length() > 0)
            {
                std::vector<std::string> ll = split(pqm->c_str(), '/');
                if (ll.size() == 3)
                {
                    add_qmc_calc(qmc, ll[0].c_str(), ll[1].c_str(), ll[2].c_str());
                    for (ti = types.begin(); (ti < types.end()); ti++)
                    {
                        if (0 == strcasecmp(ti->c_str(), ll[2].c_str()))
                        {
                            break;
                        }
                    }
                    if (ti == types.end())
                    {
                        types.push_back(ll[2].c_str());
                    }
                }
            }
            n++;
        }
    }

    for (mpi = mp.begin(); (mpi < mp.end()); mpi++)
    {
        for (ci = mpi->BeginCalculation(); (ci < mpi->EndCalculation()); ci++)
        {
            method = ci->GetMethod().c_str();
            basis  = ci->GetBasisset().c_str();
            for (ti = types.begin(); (ti < types.end()); ti++)
            {
                if ((NULL == fc_str) || (get_qmc_count(qmc, method, basis, ti->c_str()) > 0))
                {
                    if (ci->GetVal(ti->c_str(), mpo, &value, &error, vec, quadrupole))
                    {
                        add_qmc_calc(qmc, method, basis, ti->c_str());
                    }
                }
            }
        }
    }
    for (i = 0; (i < qmc->n); i++)
    {
        /* Since we initialized these we have counted one extra */
        if (NULL != fc_str)
        {
            qmc->count[i]--;
        }
    }

    return qmc;
}
