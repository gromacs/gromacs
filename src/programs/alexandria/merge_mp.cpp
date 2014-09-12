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
#include <stdlib.h>
#include <string.h>
#include "gromacs/utility/futil.h"
#include "gromacs/utility/smalloc.h"
#include "gromacs/utility/cstringutil.h"
#include "gromacs/math/vec.h"
#include "gromacs/commandline/pargs.h"
#include "gromacs/legacyheaders/copyrite.h"
#include "gromacs/utility/fatalerror.h"
#include "poldata.h"
#include "poldata_xml.h"
#include "molprop.h"
#include "molprop_util.h"
#include "molprop_xml.h"
#include "molprop_sqlite3.h"
#include "stringutil.h"

typedef struct {
    const char *iupac;
    char       *prop;
    char       *value;
    char       *ref;
} t_prop;

static int tp_comp(const void *a, const void *b)
{
    t_prop *ta = (t_prop *)a;
    t_prop *tb = (t_prop *)b;

    return strcasecmp(ta->iupac, tb->iupac);
}

static void add_properties(const char *fn, std::vector<alexandria::MolProp> &mp)
{
    alexandria::MolPropIterator mpi;
    FILE                       *fp;
    int    nprop                   = 0;
    t_prop                     *tp = NULL, key, *tpp;
    char   buf[STRLEN];
    int    nadd = 0;

    if (NULL != fn)
    {
        fp = gmx_ffopen(fn, "r");
        while (!feof(fp))
        {
            fgets2(buf, STRLEN-1, fp);
            std::vector<std::string> ptr = split(buf, '|');
            if ((ptr.size() >= 3) &&
                (ptr[0].length() > 0) && (ptr[1].length() > 0) &&
                (ptr[2].length() > 0) && (ptr[3].length() > 0))
            {
                srenew(tp, ++nprop);
                tp[nprop-1].iupac = strdup(ptr[0].c_str());
                tp[nprop-1].prop  = strdup(ptr[1].c_str());
                tp[nprop-1].value = strdup(ptr[2].c_str());
                tp[nprop-1].ref   = strdup(ptr[3].c_str());
            }
        }
        printf("Read in %d properties from %s.\n", nprop, fn);
        qsort(tp, nprop, sizeof(tp[0]), tp_comp);
        fclose(fp);
        for (mpi = mp.begin(); (mpi < mp.end()); mpi++)
        {
            key.iupac = mpi->GetIupac().c_str();
            if (NULL != key.iupac)
            {
                tpp = (t_prop *) bsearch(&key, tp, nprop, sizeof(tp[0]), tp_comp);
                if (NULL != tpp)
                {
                    alexandria::Experiment      ex(tpp->ref, (char *)"minimum");
                    alexandria::MolecularEnergy me(tpp->prop,
                                                   unit2string(eg2cKj_Mole),
                                                   atof(tpp->value), 0);
                    ex.AddEnergy(me);
                    mpi->AddExperiment(ex);
                    nadd++;
                }
            }
        }
        printf("Added properties for %d out of %d molecules.\n",
               nadd, (int)mp.size());
    }
}

static void add_charges(const char *fn, std::vector<alexandria::MolProp> &mp)
{
    alexandria::MolPropIterator mpi;
    FILE                       *fp;
    int    nprop                   = 0;
    t_prop                     *tp = NULL, key, *tpp;
    char   buf[STRLEN];
    int    nadd = 0;

    if (NULL != fn)
    {
        fp = gmx_ffopen(fn, "r");
        while (!feof(fp))
        {
            fgets2(buf, STRLEN-1, fp);
            std::vector<std::string> ptr = split(buf, '|');
            if ((ptr.size() >= 3) &&
                (ptr[0].length() > 0) && (ptr[1].length() > 0) &&
                (ptr[2].length() > 0) && (ptr[3].length() > 0))
            {
                srenew(tp, ++nprop);
                tp[nprop-1].iupac = strdup(ptr[0].c_str());
                tp[nprop-1].prop  = strdup(ptr[1].c_str());
                tp[nprop-1].value = strdup(ptr[2].c_str());
                tp[nprop-1].ref   = strdup(ptr[3].c_str());
            }
        }
        printf("Read in %d charge sets from %s.\n", nprop, fn);
        qsort(tp, nprop, sizeof(tp[0]), tp_comp);
        fclose(fp);
        for (mpi = mp.begin(); (mpi < mp.end()); mpi++)
        {
            key.iupac = mpi->GetIupac().c_str();
            if (NULL != key.iupac)
            {
                tpp = (t_prop *) bsearch(&key, tp, nprop, sizeof(tp[0]), tp_comp);
                if (NULL != tpp)
                {
                    alexandria::Experiment      ex(tpp->ref, (char *)"minimum");
                    alexandria::MolecularEnergy me(tpp->prop,
                                                   unit2string(eg2cKj_Mole),
                                                   atof(tpp->value), 0);
                    ex.AddEnergy(me);
                    mpi->AddExperiment(ex);
                    nadd++;
                }
            }
        }
        printf("Added properties for %d out of %d molecules.\n",
               nadd, (int)mp.size());
    }
}

int alex_merge_mp(int argc, char *argv[])
{
    static const char               *desc[] =
    {
        "merge_mp reads multiple molprop files and merges the molecule descriptions",
        "into a single new file. By specifying the [TT]-db[TT] option additional experimental",
        "information will be read from a SQLite3 database.[PAR]",
    };
    t_filenm                         fnm[] =
    {
        { efDAT, "-f",  "data",      ffRDMULT },
        { efDAT, "-o",  "allmols",   ffWRITE },
        { efDAT, "-di", "gentop",    ffOPTRD },
        { efDAT, "-db", "sqlite",    ffOPTRD },
        { efDAT, "-x",  "extra",     ffOPTRD },
        { efDAT, "-c",  "charges",   ffOPTRD }
    };
    int                              NFILE    = (sizeof(fnm)/sizeof(fnm[0]));
    static const char               *sort[]   = { NULL, "molname", "formula", "composition", NULL };
    static int                       compress = 1;
    t_pargs                          pa[]     =
    {
        { "-sort",   FALSE, etENUM, {sort},
          "Key to sort the final data file on." },
        { "-compress", FALSE, etBOOL, {&compress},
          "Compress output XML files" }
    };
    char                           **fns;
    int                              nfiles;
    std::vector<alexandria::MolProp> mp;
    gmx_atomprop_t                   ap;
    gmx_poldata_t                    pd;
    output_env_t                     oenv;

    if (!parse_common_args(&argc, argv, PCA_NOEXIT_ON_ARGS, NFILE, fnm,
                           sizeof(pa)/sizeof(pa[0]), pa,
                           sizeof(desc)/sizeof(desc[0]), desc,
                           0, NULL, &oenv))
    {
        return 0;
    }

    ap = gmx_atomprop_init();
    if ((pd = gmx_poldata_read(opt2fn_null("-di", NFILE, fnm), ap)) == NULL)
    {
        gmx_fatal(FARGS, "Can not read the force field information. File missing or incorrect.");
    }
    nfiles = opt2fns(&fns, "-f", NFILE, fnm);
    merge_xml(nfiles, fns, mp, NULL, NULL, NULL, ap, pd, TRUE);

    ReadSqlite3(opt2fn_null("-db", NFILE, fnm), mp);

    add_properties(opt2fn_null("-x", NFILE, fnm), mp);

    add_charges(opt2fn_null("-c", NFILE, fnm), mp);

    MolPropWrite(opt2fn("-o", NFILE, fnm), mp, compress);

    return 0;
}
