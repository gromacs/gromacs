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

#include <stdlib.h>
#include <string.h>
#ifdef HAVE_CONFIG_H
#include "config.h"
#endif
#ifdef HAVE_LIBSQLITE3
#include <sqlite3.h>
#endif
#include "molprop_sqlite3.h"

#include "gromacs/utility/cstringutil.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/futil.h"
#include "gromacs/utility/smalloc.h"

#include "stringutil.h"

typedef struct {
    const char *molname, *iupac;
} t_synonym;

int syn_comp(const void *a, const void *b)
{
    t_synonym *sa, *sb;
    sa = (t_synonym *)a;
    sb = (t_synonym *)b;

    return strcmp(sa->molname, sb->molname);
}

#ifdef HAVE_LIBSQLITE3
static void check_sqlite3(sqlite3 *db, const char *extra, int rc)
{
    const char *msg;

    if (NULL != db)
    {
        if (SQLITE_OK != sqlite3_errcode(db))
        {
            msg = sqlite3_errmsg(db);
            sqlite3_close(db);
            sqlite3_shutdown();
            gmx_fatal(FARGS, "%s: %s", extra, msg);
        }
    }
    else if (SQLITE_OK != rc)
    {
        gmx_fatal(FARGS, "%s", extra);
    }
}
#endif

void ReadSqlite3(const char                       *sqlite_file,
                 std::vector<alexandria::MolProp> &mp)
{
#ifdef HAVE_LIBSQLITE3
    alexandria::MolPropIterator mpi;
    std::string                 cas2, csid2, iupac2;

    sqlite3                    *db   = NULL;
    sqlite3_stmt               *stmt = NULL, *stmt2 = NULL;
    char sql_str[1024];
    const char                 *iupac, *cas, *csid, *prop, *unit, *source;
    double                      value, error, temperature;
    int                         cidx, rc, nbind, nexp_prop, theory;
    t_synonym                  *syn  = NULL, key, *keyptr;
    int                         nsyn = 0, maxsyn = 0;

    if (NULL == sqlite_file)
    {
        return;
    }

    check_sqlite3(NULL, "Initializing sqlite",
                  sqlite3_initialize());

    check_sqlite3(NULL, "Opening sqlite database in read-only mode",
                  sqlite3_open_v2(sqlite_file, &db, SQLITE_OPEN_READONLY, NULL));

    /* Now database is open and everything is Hunky Dory */
    printf("Opened SQLite3 database %s\n", sqlite_file);

    /* Make renaming table */
    sprintf(sql_str, "SELECT syn.name,mol.iupac FROM molecules as mol,synonyms as syn WHERE syn.molid=mol.molid ORDER by syn.name");
    if (NULL != debug)
    {
        fprintf(debug, "sql_str = '%s'\n", sql_str);
    }

    check_sqlite3(db, "Preparing statement",
                  sqlite3_prepare_v2(db, sql_str, 1+strlen(sql_str), &stmt2, NULL));
    do
    {
        rc = sqlite3_step(stmt2);
        if (SQLITE_ROW == rc)
        {
            cidx   = 0;
            if (nsyn >= maxsyn)
            {
                maxsyn += 1000;
                srenew(syn, maxsyn);
            }
            syn[nsyn].molname = strdup((char *)sqlite3_column_text(stmt2, cidx));
            cidx++;
            syn[nsyn].iupac   = strdup((char *)sqlite3_column_text(stmt2, cidx));
            cidx++;
            nsyn++;
        }
        else if (SQLITE_DONE != rc)
        {
            check_sqlite3(db, "Stepping", rc);
        }
        else
        {
            printf("There are %d synonyms for %d molecules.\n", nsyn, (int)mp.size());
        }
    }
    while (SQLITE_ROW == rc);
    check_sqlite3(db, "Resetting sqlite3 statement",
                  sqlite3_reset(stmt2));
    check_sqlite3(db, "Finalizing sqlite3 statement",
                  sqlite3_finalize(stmt2));

    /* Now present a query statement */
    nexp_prop = 0;
    sprintf(sql_str, "SELECT mol.iupac,mol.cas,mol.csid,pt.prop,pt.unit_text,gp.temperature,gp.value,gp.error,ds.theory,ds.source FROM molecules as mol,molproperty as gp,proptypes as pt, datasource as ds,phasetype as ph WHERE ((gp.phaseid=ph.phaseid) AND (ph.phase='gas') AND (mol.molid = gp.molid) AND (gp.propid = pt.propid) AND (gp.srcid = ds.srcid) AND (upper(?) = upper(mol.iupac)));");
    check_sqlite3(db, "Preparing sqlite3 statement",
                  sqlite3_prepare_v2(db, sql_str, 1+strlen(sql_str), &stmt, NULL));

    if (NULL != debug)
    {
        fprintf(debug, "sql_str = '%s'\nvariable = '%s'\n", sql_str,
                sqlite3_bind_parameter_name(stmt, 1));
        nbind = sqlite3_bind_parameter_count(stmt);
        fprintf(debug, "%d binding parameter(s) in the statement\n%s\n", nbind, sql_str);
    }
    for (mpi = mp.begin(); (mpi < mp.end()); mpi++)
    {
        key.molname = mpi->getMolname().c_str();
        key.iupac   = mpi->getIupac().c_str();
        keyptr      = (t_synonym *) bsearch((const void *)&key,
                                            (const void *)syn,
                                            nsyn, sizeof(syn[0]), syn_comp);
        if (NULL == keyptr)
        {
            fprintf(stderr, "Warning: missing iupac for %s. Will be ignored.\n",
                    key.molname);
        }
        else
        {
            if ((NULL == key.iupac) || (strcmp(key.iupac, keyptr->iupac) != 0))
            {
                if (NULL != debug)
                {
                    fprintf(debug, "Warning: incorrect iupac %s for %s - changing to %s\n",
                            key.iupac, key.molname, keyptr->iupac);
                }
                mpi->SetIupac(keyptr->iupac);
            }
            iupac = keyptr->iupac;
            if (NULL != iupac)
            {
                if (NULL != debug)
                {
                    fprintf(debug, "Going to query for '%s'\n", iupac);
                }
                check_sqlite3(db, "Binding text",
                              sqlite3_bind_text(stmt, 1, iupac, -1, SQLITE_STATIC));
                do
                {
                    rc = sqlite3_step(stmt);
                    if (SQLITE_ROW == rc)
                    {
                        /* printf("Found a row\n"); */
                        cidx   = 0;
                        iupac2 = (char *)sqlite3_column_text(stmt, cidx++);
                        if (strcasecmp(iupac, iupac2.c_str()) != 0)
                        {
                            gmx_fatal(FARGS, "Selected '%s' from database but got '%s'. WTF?!",
                                      iupac, iupac2.c_str());
                        }
                        cas            = (char *)sqlite3_column_text(stmt, cidx++);
                        csid           = (char *)sqlite3_column_text(stmt, cidx++);
                        prop           = (char *)sqlite3_column_text(stmt, cidx++);
                        unit           = (char *)sqlite3_column_text(stmt, cidx++);
                        temperature    = sqlite3_column_double(stmt, cidx++);
                        value          = sqlite3_column_double(stmt, cidx++);
                        error          = sqlite3_column_double(stmt, cidx++);
                        theory         = sqlite3_column_int(stmt, cidx++);
                        source         = (char *)sqlite3_column_text(stmt, cidx++);

                        bool bExp = (0 == theory);
                        if (bExp)
                        {
                            nexp_prop++;
                        }
                        //printf("source = %s prop = %s value = %10g bExp = %d\n",
                        //      source, prop, value, (int) bExp);
                        if (bExp)
                        {
                            alexandria::Experiment exper("unknown", "minimum");
                            if (strcasecmp(prop, "Polarizability") == 0)
                            {
                                exper.AddPolar(alexandria::MolecularPolarizability(prop, unit, temperature, 0, 0, 0, 0, 0, 0, value, 0));

                            }
                            else if (strcasecmp(prop, "dipole") == 0)
                            {
                                exper.AddDipole(alexandria::MolecularDipole(prop, unit, temperature, 0, 0, 0, value, error));
                            }
                            else if ((strcasecmp(prop, "DeltaHform") == 0) ||
                                     (strcasecmp(prop, "DeltaGform") == 0) ||
                                     (strcasecmp(prop, "DeltaSform") == 0) ||
                                     (strcasecmp(prop, "S0") == 0) ||
                                     (strcasecmp(prop, "cp") == 0) ||
                                     (strcasecmp(prop, "cv") == 0))
                            {
                                exper.AddEnergy(alexandria::MolecularEnergy(prop, unit, temperature, epGAS, value, error));
                            }
                            mpi->AddExperiment(exper);
                        }
                        else
                        {
                            alexandria::Experiment calc("gentop", source,
                                                        "-", "unknown", "minimum",
                                                        "unknown" );
                            if (strcasecmp(prop, "Polarizability") == 0)
                            {
                                alexandria::MolecularPolarizability mp(prop, unit, temperature, 0, 0, 0, 0, 0, 0, value, 0);
                                calc.AddPolar(mp);
                            }
                            else if (strcasecmp(prop, "dipole") == 0)
                            {
                                calc.AddDipole(alexandria::MolecularDipole(prop, unit, temperature, 0, 0, 0, value, error));
                            }
                            else if ((strcasecmp(prop, "DeltaHform") == 0) ||
                                     (strcasecmp(prop, "DeltaGform") == 0) ||
                                     (strcasecmp(prop, "DeltaSform") == 0) ||
                                     (strcasecmp(prop, "S0") == 0) ||
                                     (strcasecmp(prop, "cp") == 0) ||
                                     (strcasecmp(prop, "cv") == 0))

                            {
                                calc.AddEnergy(alexandria::MolecularEnergy(prop, unit, temperature, epGAS, value, error));
                            }
                            mpi->AddExperiment(calc);
                        }
                        //mpi->Stats();
                        /*
                           if (0 && (strlen(classification) > 0))
                           {
                            std::vector<std::string> class_ptr = split(classification, ';');
                            for (std::vector<std::string>::iterator cp = class_ptr.begin();
                                 cp < class_ptr.end(); ++cp)
                            {
                                mpi->AddCategory(cp->c_str());
                            }
                            }*/
                        if (strlen(cas) > 0)
                        {
                            cas2 = mpi->getCas();
                            if ((cas2.length() > 0) &&
                                (strcmp(cas, cas2.c_str()) != 0))
                            {
                                fprintf(stderr, "cas in molprop %s not the same as database %s for %s\n", cas2.c_str(), cas, iupac);
                            }
                            mpi->SetCas(cas);
                        }
                        if (strlen(csid) > 0)
                        {
                            csid2 = mpi->getCid();
                            if ((csid2.length() > 0) &&
                                (strcmp(csid, csid2.c_str()) != 0))
                            {
                                fprintf(stderr, "csid in molprop %s not the same as database %s for %s\n", csid2.c_str(), csid, iupac);
                            }
                            mpi->SetCid(csid);
                        }
                    }
                    else if (SQLITE_DONE != rc)
                    {
                        check_sqlite3(db, "Stepping", rc);
                    }
                    else if (NULL != debug)
                    {
                        fprintf(debug, "Done finding rows for %s\n", iupac);
                    }
                }
                while (SQLITE_ROW == rc);
                sqlite3_clear_bindings(stmt);
                check_sqlite3(db, "Resetting sqlite3 statement",
                              sqlite3_reset(stmt));
            }
        }
    }
    check_sqlite3(db, "Finalizing sqlite3 statement",
                  sqlite3_finalize(stmt));

    /* Seems like we're done, close down and say goodbye */
    check_sqlite3(NULL, "Closing sqlite database",
                  sqlite3_close(db));

    check_sqlite3(NULL, "Shutting down sqlite. Sqlite3 code %d.",
                  sqlite3_shutdown());
    printf("Extracted %d data points from sql database\n", nexp_prop);
#else
    fprintf(stderr, "No support for sqlite3 database in this executable.\n");
    fprintf(stderr, "Please rebuild gromacs with cmake flag -DGMX_SQLITE3=ON set.\n");
#endif
}
