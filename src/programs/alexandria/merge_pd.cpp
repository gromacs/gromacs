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
#include "gmxpre.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "gromacs/commandline/pargs.h"
#include "gromacs/math/units.h"
#include "gromacs/math/vec.h"
#include "gromacs/statistics/statistics.h"
#include "gromacs/topology/topology.h"
#include "gromacs/utility/arraysize.h"
#include "gromacs/utility/cstringutil.h"
#include "gromacs/utility/exceptions.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/futil.h"
#include "gromacs/utility/smalloc.h"

#include "poldata.h"
#include "poldata_xml.h"




enum EemAtomProps {
    eEMJ00   = 0, 
    eEMChi   = 1, 
    eEMZeta  = 2,  
    eEMNR    = 3
};

typedef struct {
    EemAtomProps  eEM;
    const char   *name;
} t_eemAtom_props;

t_eemAtom_props eemAtom_props[eEMNR] = {
    {eEMJ00,   "j0"},
    {eEMChi,   "chi"},
    {eEMZeta,  "Zeta"}
};

static EemAtomProps name2eemprop(const std::string name)
{
    for (auto i = 0; i < eEMNR; i++)
    {
        if (strcasecmp(name.c_str(), eemAtom_props[i].name) == 0)
        {
            return eemAtom_props[i].eEM;
        }
    }
    return eEMNR;
}

static const char *getEemAtomName(EemAtomProps eem)
{
    for (auto i = 0; i < eEMNR; i++)
    {
        if (eem == eemAtom_props[i].eEM)
        {
            return eemAtom_props[i].name;
        }
    }
    return nullptr;
}

static void merge_J00Chi(std::vector<alexandria::Poldata>     pds,
                         alexandria::ChargeDistributionModel  ieqd,
                         alexandria::Poldata                 &pdout,
                         EemAtomProps                         eematp)
{
    double average = 0;
    double sigma   = 0;
    
    auto nAtypes = pdout.getNatypes();
    gmx_stats_t lsq[nAtypes];
    
    for (size_t i = 0; i < nAtypes; i++)
    {
        lsq[i] =  gmx_stats_init();
    }    
    int j = 0;
    for (auto atp = pdout.getAtypeBegin(); 
         atp < pdout.getAtypeEnd(); atp++, j++)
    {
        for (const auto& pd : pds)
        {
            auto ei = pd.findEem(ieqd, atp->getType());
            if (ei != pd.EndEemprops())
            {
                switch (eematp)
                {
                    case eEMJ00:
                        gmx_stats_add_point(lsq[j], 0, ei->getJ0(), 0, 0);
                    break;
                    case eEMChi:
                        gmx_stats_add_point(lsq[j], 0, ei->getChi0(), 0, 0);
                    break;
                    case eEMZeta:
                    case eEMNR:
                    break; 
                }
            }
            else
            {
                gmx_fatal(FARGS, "%s atomtype does not exists in Eeemprops.\n", 
                          atp->getType().c_str());
            }
        }
    }
    j = 0;
    for (auto atp = pdout.getAtypeBegin(); 
         atp < pdout.getAtypeEnd(); atp++, j++)
    {
        auto ei = pdout.findEem(ieqd, atp->getType());
        if (ei != pdout.EndEemprops())
        {
            if ((estatsOK == gmx_stats_get_average(lsq[j], &average))&&
                (estatsOK == gmx_stats_get_sigma(lsq[j], &sigma)))
            {
                switch (eematp)
                {
                    case eEMJ00:
                    {
                        ei->setJ0(average);
                        ei->setJ0_sigma(sigma);
                    }
                    break;
                    case eEMChi:
                    {
                        ei->setChi0(average);
                        ei->setChi0_sigma(sigma);
                    }
                    break;
                    case eEMZeta:
                    case eEMNR:
                    break; 
                }
            }
            else
            {
                gmx_fatal(FARGS, "estats is not OK for %s.\n", 
                          atp->getType().c_str());
            }
        }
        else
        {
            gmx_fatal(FARGS, "%s atomtype does not exists in Eeemprops.\n", 
                      atp->getType().c_str());
        }
    }    
    for (size_t i = 0; i < nAtypes; i++)
    {
         gmx_stats_free(lsq[i]);
    }
}

static void merge_zeta(std::vector<alexandria::Poldata>     pds,
                       alexandria::ChargeDistributionModel  ieqd,
                       alexandria::Poldata                 &pdout)
{
    real core_ave  = 0;
    real core_sig  = 0;
    real shell_ave = 0;
    real shell_sig = 0;
    
    char  zstr[STRLEN];
    char  z_sig[STRLEN];
    char  buf[STRLEN];
    char  buf_sig[STRLEN];
    
    auto nAtypes = pdout.getNatypes();
    
    gmx_stats_t core[nAtypes];
    gmx_stats_t shell[nAtypes];
 
       
    for (size_t i = 0; i < nAtypes; i++)
    {
        core[i]  =  gmx_stats_init();
        shell[i] =  gmx_stats_init();
    }     
    int j = 0;
    for (auto atp = pdout.getAtypeBegin(); 
         atp < pdout.getAtypeEnd(); atp++, j++)
    {
        for (const auto& pd : pds)
        {
            auto ei = pd.findEem(ieqd, atp->getType());
            if (ei != pd.EndEemprops())
            {
                auto nzeta = pd.getNzeta(ieqd, ei->getName());
                if (nzeta == 1)
                {
                    gmx_stats_add_point(core[j], 0, pd.getZeta(ieqd, ei->getName(), 0), 0, 0);
                }
                else if (nzeta == 2)
                {
                    gmx_stats_add_point(core[j],  0, pd.getZeta(ieqd, ei->getName(), 0), 0, 0);
                    gmx_stats_add_point(shell[j], 0, pd.getZeta(ieqd, ei->getName(), 1), 0, 0);
                }
                else
                {
                    gmx_fatal(FARGS, "The number of Zeta is wrong for %s atomtype.\n", 
                              ei->getName());
                }
                
            }   
        }
    }
    j = 0;
    for (auto atp = pdout.getAtypeBegin(); 
         atp < pdout.getAtypeEnd(); atp++, j++)
    {
        auto ei = pdout.findEem(ieqd, atp->getType());
        if (ei != pdout.EndEemprops())
        {
            zstr[0]  = '\0';
            z_sig[0] = '\0';
            auto nzeta  = pdout.getNzeta(ieqd, ei->getName());            
            if (nzeta == 1)
            {
                if ((estatsOK == gmx_stats_get_average(core[j], &core_ave))&&
                    (estatsOK == gmx_stats_get_sigma(core[j],   &core_sig)))
                {
                    sprintf(buf, "%g ", core_ave);
                    sprintf(buf_sig, "%g ", core_sig);
                    strcat(zstr, buf);
                    strcat(z_sig, buf_sig);
                    ei->setZetastr(zstr);
                    ei->setZeta_sigma(z_sig);
                }
                else
                {
                    gmx_fatal(FARGS, "estats is not OK for %s.\n", 
                              ei->getName());
                }
            }
            else if (nzeta == 2)
            {
                if ((estatsOK == gmx_stats_get_average(core[j],  &core_ave))  &&
                    (estatsOK == gmx_stats_get_sigma(core[j],    &core_sig))  &&
                    (estatsOK == gmx_stats_get_average(shell[j], &shell_ave)) &&
                    (estatsOK == gmx_stats_get_sigma(shell[j],   &shell_sig)))
                {
                    sprintf(buf, "%g %g ", core_ave, shell_ave);
                    sprintf(buf_sig, "%g %g ", core_sig, shell_sig);
                    strcat(zstr, buf);
                    strcat(z_sig, buf_sig);
                    ei->setZetastr(zstr);
                    ei->setZeta_sigma(z_sig);
                }
                else
                {
                    gmx_fatal(FARGS, "estats is not OK for %s.\n", 
                              ei->getName());
                }
            }
            else
            {
                gmx_fatal(FARGS, "The number of Zeta is wrong for %s atomtype.\n", 
                          ei->getName());
            }
        }        
    }    
    for (size_t i = 0; i < nAtypes; i++)
    {
         gmx_stats_free(core[i]);
         gmx_stats_free(shell[i]);
    }
}



int alex_merge_pd(int argc, char *argv[])
{
    static const char               *desc[] =
    {
        "merge_pd reads multiple gentop files and merges them",
        "into a single new gentop file.",
    };    
    t_filenm                         fnm[] =
    {
        { efDAT, "-di", "pdin",  ffRDMULT },
        { efDAT, "-do", "pdout", ffWRITE }
    };
    int                              NFILE       = asize(fnm);;
    
    static int                       compress    = 1;    
    static const char               *cqdist[]    = {nullptr, "AXp", "AXg", "AXs", "AXpp", "AXpg", "AXps", nullptr};
    static const char               *eemprop[]   = {nullptr, "j0", "chi", "zeta", nullptr};
    
    t_pargs                          pa[]        =
    {
        { "-compress", FALSE, etBOOL, {&compress},
          "Compress output XML files" },
        { "-qdist",   FALSE, etENUM, {cqdist},
          "Model used for charge distribution" },
        { "-eemprop", FALSE, etENUM, {eemprop},
          "Atomic property used to describe molecular eletric properties." }
    };
    char                           **fns;
    int                              nfiles;
    std::vector<alexandria::Poldata> pds;
    alexandria::Poldata              pdout;
    gmx_atomprop_t                   aps;
    gmx_output_env_t                *oenv;

    if (!parse_common_args(&argc, argv, PCA_NOEXIT_ON_ARGS, NFILE, fnm,
                           asize(pa), pa, asize(desc), desc, 0, nullptr, &oenv))
    {
        return 0;
    }

    aps = gmx_atomprop_init();
    
    /*
      Read all the gentop files.
     */
    nfiles = opt2fns(&fns, "-di", NFILE, fnm);
    for (int i = 0; i < nfiles; i++)
    {
        alexandria::Poldata pd;
        readPoldata(fns[i], pd, aps);
        pds.push_back(std::move(pd));
    }
    
    /* 
       Copy one of the gentop files into pdout.
       Later, we update different parts
       of the pdout. 
     */
    readPoldata(fns[0], pdout, aps);
    
    alexandria::ChargeDistributionModel   ieqd   = alexandria::name2eemtype(cqdist[0]);
    EemAtomProps                          eem    = name2eemprop(eemprop[0]);  
      
    if (eem == eEMJ00 || eem == eEMChi)
    {
        merge_J00Chi(pds, ieqd, pdout, eem);
    }
    else if (eem == eEMZeta)
    {
        merge_zeta(pds, ieqd, pdout);
    }
    else
    {
        gmx_fatal(FARGS, "There is no electric atomic property called %s in alexandria.\n", eemprop[0]);
    }
    
    alexandria::writePoldata(opt2fn("-do", NFILE, fnm), pdout, 0);
    
    return 0;
}

