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
 * \author  Mohammad Mehdi Ghahremanpour <mohammad.ghahremanpour@icm.uu.se>
 * \author David van der Spoel <david.vanderspoel@icm.uu.se>
 */

#include <vector>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <random>
#include <strings.h>

#include "gromacs/commandline/filenm.h"
#include "gromacs/commandline/pargs.h"
#include "gromacs/gmxlib/network.h"
#include "gromacs/hardware/detecthardware.h"
#include "gromacs/mdrunutility/mdmodules.h"
#include "gromacs/utility/arraysize.h"
#include "gromacs/utility/coolstuff.h"
#include "gromacs/utility/futil.h"
#include "gromacs/utility/cstringutil.h"
#include "gromacs/utility/smalloc.h"
#include "gromacs/utility/textreader.h"

#include "composition.h" 
#include "getmdlogger.h"
#include "moldip.h"
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

    while(tr.readLine(&tmp)) 
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

            for(j = 0; (j < (int)imsNR); j++)
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
    static int                  minmol    = 3;
    static int                  maxatempt = 5000;
    static char                *opt_elem  = nullptr;
    static char                *lot       = (char *)"B3LYP/aug-cc-pVTZ";
        
    t_pargs                     pa[]      = 
    {
        { "-nsample",   FALSE, etINT, {&nsample},
          "Number of replicas." },
        { "-minmol",    FALSE, etINT, {&minmol},
          "Minimum number of molecules per atom types." },
        { "-maxatempt", FALSE, etINT, {&maxatempt},
          "Maximum number of atempts to sample minmol molecules per atom types." },
        { "-opt_elem",  FALSE, etSTR, {&opt_elem},
          "Space-separated list of atom types to select molecules. If this variable is not set, all elements will be used." },
        { "-lot",       FALSE, etSTR,  {&lot},
          "Use this method and level of theory when selecting molecules." }        
    };
    
    gmx_output_env_t       *oenv;
    alexandria::MolDip      mdp;
    alexandria::MolSelect   gms;
    time_t                  my_t;
    FILE                   *fp;
    
    t_commrec              *cr     = init_commrec(); 
    gmx::MDLogger           mdlog  = getMdLogger(cr, stdout);
    gmx_hw_info_t          *hwinfo = gmx_detect_hardware(mdlog, cr, false);
    
    if (!parse_common_args(&argc, argv, PCA_CAN_VIEW, NFILE, fnm, asize(pa), pa,
                           asize(desc), desc, 0, nullptr, &oenv))
    {
        sfree(cr);
        return 0;
    }
    if (MASTER(cr))
    {
        printf("There are %d threads/processes.\n", cr->nnodes);
    }    
    if (MASTER(cr))
    {
        fp = gmx_ffopen(opt2fn("-g", NFILE, fnm), "w");

        time(&my_t);
        fprintf(fp, "# This file was created %s", ctime(&my_t));
        fprintf(fp, "# alexandria is part of GROMACS:\n#\n");
        fprintf(fp, "# %s\n#\n", gmx::bromacs().c_str());
    }
    else
    {
        fp = nullptr;
    }    
    if (MASTER(cr))
    {
        gms.read(opt2fn_null("-sel", NFILE, fnm));
                
        mdp.Init(cr, false, false, eqdAXp, eqgESP,
                 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
                 0, 0, 0, 0, 0, 0, false, 0, 
                 false, false, hwinfo,
                 false, minmol, false);
        
        mdp.Read(fp ? fp : (debug ? debug : nullptr),
                 opt2fn("-f", NFILE, fnm),
                 opt2fn_null("-d", NFILE, fnm),
                 false, opt_elem, nullptr, lot,
                 gms, 0, true, false, false, false,
                 false, nullptr, 0, 0);
        
        char  buf[STRLEN];      
        for (int i = 0; i < nsample; i++)
        {
            sprintf(buf, "%s_%d.dat", fnm[2].fn, i); 
            fp         = gmx_ffopen(buf, "w");
            sample_molecules(fp, mdp.mymol_, mdp.pd_, mdp.mindata_, maxatempt);
            gmx_ffclose(fp);
        }
        done_filenms(NFILE, fnm);
    }
    
    return 0;
}
