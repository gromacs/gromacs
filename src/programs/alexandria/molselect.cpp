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

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <random>
#include <vector>
#include <strings.h>

#include "gromacs/commandline/filenm.h"
#include "gromacs/commandline/pargs.h"
#include "gromacs/gmxlib/network.h"
#include "gromacs/mdrunutility/mdmodules.h"
#include "gromacs/utility/arraysize.h"
#include "gromacs/utility/coolstuff.h"
#include "gromacs/utility/cstringutil.h"
#include "gromacs/utility/futil.h"
#include "gromacs/utility/smalloc.h"
#include "gromacs/utility/textreader.h"

#include "composition.h"
#include "getmdlogger.h"
#include "molgen.h"
#include "molprop.h"
#include "molprop_xml.h"
#include "molselect.h"
#include "mymol.h"
#include "poldata.h"
#include "poldata_xml.h"
#include "stringutil.h"

static const char *ims_names[imsNR] = { "Train", "Test", "Ignore", "Unknown" };

const char *iMolSelectName(iMolSelect ims)
{
    int i = static_cast<int>(ims);

    return ims_names[i];
}

namespace alexandria
{

static void sample_molecules(FILE                            *fp,
                             std::vector<alexandria::MyMol>   mols,
                             alexandria::Poldata              pd,
                             int                              minmol,
                             int                              maxatempt)
{

    int nmol      = 0;
    int atempt    = 0;

    CompositionSpecs                      cs;
    std::random_device                    rd;
    std::mt19937                          gen(rd());
    std::uniform_int_distribution<>       dis(0, mols.size()-1);
    std::vector<alexandria::MyMol>        sample;

    const char  *alexandria = cs.searchCS(alexandria::iCalexandria)->name();
    for (auto atp = pd.getAtypeBegin(); atp < pd.getAtypeEnd(); atp++)
    {
        if (atp->getElem() != "H")
        {
            nmol   = 0;
            atempt = 0;
            do
            {
                auto found = false;
                auto mol   = mols[dis(gen)];
                auto mci   = mol.molProp()->SearchMolecularComposition(alexandria);
                for (auto ani = mci->BeginAtomNum(); (!found) && (ani < mci->EndAtomNum()); ++ani)
                {
                    if (atp->getType() == ani->getAtom())
                    {
                        if (std::find(sample.begin(), sample.end(), mol) == sample.end())
                        {
                            sample.push_back(mol);
                            found = true;
                            nmol++;
                        }
                    }
                }
                atempt++;
            }
            while (nmol < minmol && atempt < maxatempt);
            if (debug && atempt >= maxatempt)
            {
                fprintf(debug, "Randomly picked only %d out of required %d molecules for %s after %d attempts\n", nmol, minmol, atp->getType().c_str(), atempt);
            }
        }
    }
    for (const auto &mol : sample)
    {
        fprintf(fp, "%s|Train\n", mol.molProp()->getMolname().c_str());
    }
}

void MolSelect::read(const char *fn)
{
    gmx::TextReader tr(fn);
    std::string     tmp;
    int             index = 0;

    while (tr.readLine(&tmp))
    {
        while (!tmp.empty() && tmp[tmp.length()-1] == '\n')
        {
            tmp.erase(tmp.length()-1);
        }
        auto ptr = split(tmp, '|');
        if ((ptr.size() == 2) && (ptr[0].length() > 1)
            && (ptr[1].length() > 1))
        {
            iMolSelect status;
            int        j;

            for (j = 0; (j < (int)imsNR); j++)
            {
                if (strcasecmp(ims_names[j], ptr[1].c_str()) == 0)
                {
                    break;
                }
            }
            if (j < imsNR)
            {
                status = static_cast<iMolSelect>(j);
            }
            else
            {
                status = imsUnknown;
                fprintf(stderr, "Unknown status '%s' for molecule %s on line %d in file %s\n",
                        ptr[1].c_str(), ptr[0].c_str(), index, fn);
            }
            ims_.push_back(IMolSelect(ptr[0], status, index++));
        }
        else
        {
            fprintf(stderr, "Invalid selection file\n");
        }
    }
}

iMolSelect MolSelect::status(const std::string &iupac) const
{
    auto imi = std::find_if(ims_.begin(), ims_.end(),
                            [iupac](IMolSelect const &i)
                            {
                                return i.iupac().compare(iupac) == 0;
                            });

    if (imi != ims_.end())
    {
        return imi->status();
    }

    return imsUnknown;
}

int MolSelect::index(const std::string &iupac) const
{
    auto imi = std::find_if(ims_.begin(), ims_.end(),
                            [iupac](IMolSelect const &i)
                            {
                                return i.iupac().compare(iupac) == 0;
                            });

    if (imi != ims_.end())
    {
        return imi->index();
    }

    return -1;
}

void printAtomtypeStatistics(FILE                                 *fp,
                             const alexandria::Poldata            &pd,
                             const std::vector<alexandria::MyMol> &mymol)
{
    struct NN
    {
        std::string name;
        int         count;
    };
    std::vector<NN> nn;
    for (auto atype = pd.getAtypeBegin(); atype < pd.getAtypeEnd(); ++atype)
    {
        struct NN n;
        n.name   = atype->getType();
        n.count  = 0;
        nn.push_back(n);
    }
    for (auto mol : mymol)
    {
        int ntypes = get_atomtype_ntypes(mol.atype_);
        for (int i = 0; i < ntypes; i++)
        {
            char *tp = get_atomtype_name(i, mol.atype_);
            for (auto &n : nn)
            {
                if (n.name.compare(tp) == 0)
                {
                    n.count += 1;
                    break;
                }
            }
        }
    }
    fprintf(fp, "Atomtype     Count\n");
    for (const auto &n : nn)
    {
        fprintf(fp, "%-8s  %8d\n", n.name.c_str(), n.count);
    }
}

}

int alex_molselect(int argc, char *argv[])
{
    static const char               *desc[] = {
        "molselect generates random samples from molprop database"
    };

    t_filenm                         fnm[] =
    {
        { efDAT, "-f",    "allmols",   ffOPTRD },
        { efDAT, "-d",    "gentop",    ffOPTRD },
        { efDAT, "-o",    "selection", ffWRITE },
        { efLOG, "-g",    "molselect", ffWRITE },
        { efDAT, "-sel",  "molselect", ffREAD  },
    };

    const  int                  NFILE     = asize(fnm);

    static int                  nsample   = 1;
    static int                  maxatempt = 5000;
    static char                *opt_elem  = nullptr;
    static gmx_bool             bZero     = TRUE;
    t_pargs                     pa[]      =
    {
        { "-nsample",   FALSE, etINT, {&nsample},
          "Number of replicas." },
        { "-zero_dipole",    FALSE, etBOOL, {&bZero},
          "Take into account molecules with zero dipoles." },
        { "-maxatempt", FALSE, etINT, {&maxatempt},
          "Maximum number of atempts to sample mindata molecules per atom types." },
        { "-opt_elem",  FALSE, etSTR, {&opt_elem},
          "Space-separated list of atom types to select molecules. If this variable is not set, all elements will be used." }
    };

    gmx_output_env_t       *oenv;
    alexandria::MolGen      mdp;
    alexandria::MolSelect   gms;
    time_t                  my_t;
    FILE                   *fp;

    std::vector<t_pargs>    pargs;
    for (size_t i = 0; i < sizeof(pa)/sizeof(pa[0]); i++)
    {
        pargs.push_back(pa[i]);
    }
    mdp.addOptions(&pargs);
    if (!parse_common_args(&argc, argv, PCA_CAN_VIEW, NFILE, fnm, 
                           pargs.size(), pargs.data(),
                           asize(desc), desc, 0, nullptr, &oenv))
    {
        return 0;
    }
    mdp.optionsFinished();

    fp = gmx_ffopen(opt2fn("-g", NFILE, fnm), "w");

    time(&my_t);
    fprintf(fp, "# This file was created %s", ctime(&my_t));
    fprintf(fp, "# alexandria is part of GROMACS:\n#\n");
    fprintf(fp, "# %s\n#\n", gmx::bromacs().c_str());

    gms.read(opt2fn_null("-sel", NFILE, fnm));

    mdp.Read(fp ? fp : (debug ? debug : nullptr),
             opt2fn("-f", NFILE, fnm),
             opt2fn_null("-d", NFILE, fnm),
             bZero, opt_elem, nullptr,
             gms, true, false, false,
             false, true, nullptr);

    printAtomtypeStatistics(fp, mdp.poldata(), mdp.mymols());
    
    for (int i = 0; i < nsample; i++)
    {
        char  buf[STRLEN];
        sprintf(buf, "%s_%d.dat", fnm[2].fn, i);
        FILE *dat = gmx_ffopen(buf, "w");
        sample_molecules(dat, mdp.mymols(), mdp.poldata(),
                         mdp.mindata(), maxatempt);
        gmx_ffclose(dat);
    }
    gmx_ffclose(fp);
    done_filenms(NFILE, fnm);

    return 0;
}
