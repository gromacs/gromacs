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

#include "molprop_util.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <algorithm>
#include <vector>

#include "gromacs/fileio/xvgr.h"
#include "gromacs/gmxpreprocess/gen_ad.h"
#include "gromacs/gmxpreprocess/gpp_atomtype.h"
#include "gromacs/gmxpreprocess/gpp_nextnb.h"
#include "gromacs/gmxpreprocess/toputil.h"
#include "gromacs/math/units.h"
#include "gromacs/math/utilities.h"
#include "gromacs/math/vec.h"
#include "gromacs/pbcutil/pbc.h"
#include "gromacs/topology/atomprop.h"
#include "gromacs/topology/topology.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/futil.h"
#include "gromacs/utility/real.h"
#include "gromacs/utility/smalloc.h"
#include "gromacs/utility/stringutil.h"

#include "composition.h"
#include "gentop_core.h"
#include "gentop_vsite.h"
#include "gmx_simple_comm.h"
#include "molprop.h"
#include "molprop_xml.h"
#include "molselect.h"
#include "poldata.h"
#include "poldata_xml.h"
#include "stringutil.h"

namespace alexandria
{


void generate_composition(std::vector<MolProp> &mp,
                          const Poldata        &pd)
{
    int              nOK = 0;
    CompositionSpecs cs;

    for (auto &mpi : mp)
    {
        for (alexandria::CompositionSpecIterator csi = cs.beginCS(); (csi < cs.endCS()); ++csi)
        {
            mpi.DeleteComposition(csi->name());
        }
        if (true == mpi.GenerateComposition(pd))
        {
            nOK++;
        }
        else if (debug)
        {
            fprintf(debug, "Failed to make composition for %s\n",
                    mpi.getMolname().c_str());
        }
    }
    if (mp.size() > 1)
    {
        printf("Generated composition for %d out of %d molecules.\n",
               nOK, (int)mp.size());
    }
}

void generate_formula(std::vector<MolProp> &mp,
                      gmx_atomprop_t        ap)
{
    for (auto &mpi : mp)
    {
        mpi.GenerateFormula(ap);
    }
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

int merge_doubles(std::vector<alexandria::MolProp> &mp, char *doubles,
                  bool bForceMerge)
{
    alexandria::MolPropIterator mpi, mmm[2];
    std::string                 molname[2];
    std::string                 form[2];

    FILE *fp;
    int   i, ndouble = 0;
    bool  bForm, bName, bDouble;
    int   nwarn = 0;
    int   cur   = 0;
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
        molname[cur] = mpi->getMolname();
        form[cur]    = mpi->formula();
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
                    nwarn += mmm[prev]->Merge(mmm[cur]);
                    mpi    = mp.erase(mmm[cur]);

                    bDouble = true;
                    ndouble++;
                    if (mpi != mp.end())
                    {
                        mmm[cur]      = mpi;
                        molname[cur]  = mmm[cur]->getMolname();
                        form[cur]     = mmm[cur]->formula();
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
    return nwarn;
}

static void dump_mp(std::vector<alexandria::MolProp> mp)
{
    alexandria::MolPropIterator mpi;
    FILE *fp;

    fp = fopen("dump_mp.dat", "w");

    for (mpi = mp.begin(); (mpi < mp.end()); mpi++)
    {
        fprintf(fp, "%-20s  %s\n", mpi->formula().c_str(),
                mpi->getMolname().c_str());
    }

    fclose(fp);
}

int merge_xml(int nfile, char **filens,
              std::vector<alexandria::MolProp> &mpout,
              char *outf, char *sorted, char *doubles,
              gmx_atomprop_t ap,
              const Poldata &pd,
              bool bForceMerge)
{
    std::vector<alexandria::MolProp> mp;
    alexandria::MolPropIterator      mpi;
    int       i, npout = 0, tmp;
    int       nwarn = 0;

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
    if (NULL != debug)
    {
        fprintf(debug, "mpout.size() = %u mpout.max_size() = %u\n",
                (unsigned int)mpout.size(), (unsigned int)mpout.max_size());
        for (mpi = mpout.begin(); (mpi < mpout.end()); mpi++)
        {
            mpi->Dump(debug);
        }
    }
    MolSelect gms;
    MolPropSort(mpout, MPSA_MOLNAME, NULL, gms);
    nwarn += merge_doubles(mpout, doubles, bForceMerge);
    printf("There were %d total molecules before merging, %d after.\n",
           tmp, (int)mpout.size());

    if (outf)
    {
        printf("There are %d entries to store in output file %s\n", npout, outf);
        MolPropWrite(outf, mpout, false);
    }
    if (sorted)
    {
        MolPropSort(mpout, MPSA_FORMULA, NULL, gms);
        MolPropWrite(sorted, mpout, false);
        dump_mp(mpout);
    }
    return nwarn;
}

static bool comp_mp_molname(alexandria::MolProp ma,
                            alexandria::MolProp mb)
{
    std::string mma = ma.getMolname();
    std::string mmb = mb.getMolname();

    return (mma.compare(mmb) < 0);
}

static bool comp_mp_formula(alexandria::MolProp ma,
                            alexandria::MolProp mb)
{
    std::string fma = ma.formula();
    std::string fmb = mb.formula();

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
    return (ma.getIndex() < mb.getIndex());
}

alexandria::MolPropIterator SearchMolProp(std::vector<alexandria::MolProp> &mp,
                                          const char                       *name)
{
    alexandria::MolPropIterator mpi;
    for (mpi = mp.begin(); (mpi < mp.end()); mpi++)
    {
        if (strcmp(mpi->getMolname().c_str(), name) == 0)
        {
            break;
        }
    }
    return mpi;
}

void MolPropSort(std::vector<alexandria::MolProp> &mp,
                 MolPropSortAlgorithm mpsa, gmx_atomprop_t apt,
                 const MolSelect &gms)
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
            if (gms.nMol() > 0)
            {
                for (alexandria::MolPropIterator mpi = mp.begin(); (mpi < mp.end()); mpi++)
                {
                    int index = gms.index(mpi->getIupac());
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

void QmCount::addConf(const std::string &conformation)
{
    auto s = std::find_if(conf_.begin(), conf_.end(),
                          [conformation](const std::string &s)
                          { return (s.compare(conformation) == 0); });
    if (s == conf_.end())
    {
        conf_.push_back(conformation);
    }
}

int QmCount::qmCalcCount(const std::string &method,
                         const std::string &basis,
                         const std::string &type) const
{
    auto a = findCalc(method, basis, type);
    if (qmc_.end() == a)
    {
        return 0;
    }
    else
    {
        return a->count();
    }
}

void QmCount::addCalc(const std::string &method,
                      const std::string &basis,
                      const std::string &type)
{
    auto a = findCalc(method, basis, type);

    if (qmc_.end() == a)
    {
        QmCalc qmc(method, basis, type);
        qmc_.push_back(qmc);
    }
    else
    {
        a->increment();
    }
}

void find_calculations(std::vector<alexandria::MolProp> &mp,
                       MolPropObservable                 mpo,
                       const char                       *fc_str,
                       QmCount                          *qmc)
{
    std::vector<std::string> types;

    for (auto mpi = mp.begin(); (mpi < mp.end()); mpi++)
    {
        for (alexandria::ExperimentIterator ei = mpi->BeginExperiment(); (ei < mpi->EndExperiment()); ei++)
        {
            qmc->addConf(ei->getConformation());
        }
    }
    if (NULL != fc_str)
    {
        std::vector<std::string> qm = split(fc_str, ':');
        int n = 0;
        for (auto pqm = qm.begin(); (pqm < qm.end()); ++pqm)
        {
            if (pqm->length() > 0)
            {
                std::vector<std::string> ll = split(pqm->c_str(), '/');
                if (ll.size() == 3)
                {
                    std::vector<std::string>::iterator ti;
                    qmc->addCalc(ll[0], ll[1], ll[2]);
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

    for (auto mpi = mp.begin(); (mpi < mp.end()); mpi++)
    {
        for (auto ci = mpi->BeginExperiment(); (ci < mpi->EndExperiment()); ++ci)
        {
            if (dsExperiment ==  ci->dataSource())
            {
                continue;
            }
            for (auto ti = types.begin(); (ti < types.end()); ti++)
            {
                if ((NULL == fc_str) || (qmc->qmCalcCount(ci->getMethod(),
                                                          ci->getBasisset(),
                                                          *ti) > 0))
                {
                    double T, value, error, vec[3];
                    tensor quadrupole;
                    if (ci->getVal(ti->c_str(), mpo, &value, &error, &T,
                                   vec, quadrupole))
                    {
                        qmc->addCalc(ci->getMethod(),
                                     ci->getBasisset(),
                                     *ti);
                    }
                }
            }
        }
    }
    for (auto q = qmc->beginCalc(); q < qmc->endCalc(); ++q)
    {
        /* Since we initialized these we have counted one extra */
        if (NULL != fc_str)
        {
            q->decrement();
        }
    }
}

} // namespace alexandria
