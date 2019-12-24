/*
 * This source file is part of the Alexandria Chemistry Toolkit.
 *
 * Copyright (C) 2014-2020 
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

#include <math.h> 
#include <stdlib.h>
#include <string.h>
#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <algorithm>
#include <vector>

#ifdef HAVE_LIBSQLITE3
#include <sqlite3.h>
#endif
#include "molprop_sqlite3.h"


class Synonym
{
    private:
        std::string molname_;
        std::string iupac_;
    public:
        Synonym(const std::string molname,
                const std::string iupac) : molname_(molname), iupac_(iupac) {}

        const std::string &iupac() const { return iupac_; }

        const std::string &molname() const { return molname_; }
};

class Classes
{
    private:
        std::string              iupac_;
        std::vector<std::string> classes_;
    public:
        Classes(const std::string iupac) : iupac_(iupac) {}

        void addClass(const std::string &klas) { classes_.push_back(klas); }

        const std::string &iupac() const { return iupac_; }

        const std::vector<std::string>::iterator classBegin() { return classes_.begin(); }

        const std::vector<std::string>::iterator classEnd() { return classes_.end(); }
};

#ifdef HAVE_LIBSQLITE3
static void check_sqlite3(sqlite3 *db, const char *extra, int rc)
{
    const char *msg;

    if (nullptr != db)
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

static void getSynonyms(sqlite3              *db,
                        std::vector<Synonym> &syn,
                        int                   nMol)
{
    sqlite3_stmt *stmt2 = nullptr;
    char          sql_str[1024];
    int           rc;

    /* Make renaming table */
    snprintf(sql_str, sizeof(sql_str),
             "SELECT syn.name,mol.iupac FROM molecules as mol,synonyms as syn WHERE syn.molid=mol.molid ORDER by syn.name");

    if (nullptr != debug)
    {
        fprintf(debug, "sql_str = '%s'\n", sql_str);
    }

    check_sqlite3(db, "Preparing statement",
                  sqlite3_prepare_v2(db, sql_str, 1+strlen(sql_str), &stmt2, nullptr));
    do
    {
        rc = sqlite3_step(stmt2);
        if (SQLITE_ROW == rc)
        {
            syn.push_back(Synonym((char *)sqlite3_column_text(stmt2, 0),
                                  (char *)sqlite3_column_text(stmt2, 1)));
        }
        else if (SQLITE_DONE != rc)
        {
            check_sqlite3(db, "Stepping", rc);
        }
        else
        {
            printf("There are %d synonyms for %d molecules.\n",
                   static_cast<int>(syn.size()), nMol);
        }
    }
    while (SQLITE_ROW == rc);
    check_sqlite3(db, "Resetting sqlite3 statement",
                  sqlite3_reset(stmt2));
    check_sqlite3(db, "Finalizing sqlite3 statement",
                  sqlite3_finalize(stmt2));
}

static void getClasses(sqlite3              *db,
                       std::vector<Classes> &classes,
                       int                   nMol)
{
    sqlite3_stmt *stmt2 = nullptr;
    char          sql_str[1024];
    int           rc;

    /* Make renaming table */
    snprintf(sql_str, sizeof(sql_str),
             "SELECT mol.iupac,class.class FROM molecules as mol,classification as class,link_mol_class as lmc WHERE (lmc.molid=mol.molid) and (lmc.classid=class.classid) ORDER by mol.iupac");

    if (nullptr != debug)
    {
        fprintf(debug, "sql_str = '%s'\n", sql_str);
    }

    check_sqlite3(db, "Preparing statement",
                  sqlite3_prepare_v2(db, sql_str, 1+strlen(sql_str), &stmt2, nullptr));
    do
    {
        rc = sqlite3_step(stmt2);
        if (SQLITE_ROW == rc)
        {
            const char *iupac = (char *)sqlite3_column_text(stmt2, 0);
            const char *klass = (char *)sqlite3_column_text(stmt2, 1);
            auto        s     = std::find_if(classes.begin(), classes.end(),
                                             [iupac](Classes const &c)
                                             { return c.iupac().compare(iupac) == 0; });
            if (s == classes.end())
            {
                std::string i(iupac);
                classes.push_back(i);
                classes.back().addClass(klass);
            }
            else
            {
                s->addClass(klass);
            }
        }
        else if (SQLITE_DONE != rc)
        {
            check_sqlite3(db, "Stepping", rc);
        }
        else
        {
            printf("There are %d classes for %d molecules.\n",
                   static_cast<int>(classes.size()), nMol);
        }
    }
    while (SQLITE_ROW == rc);
    check_sqlite3(db, "Resetting sqlite3 statement",
                  sqlite3_reset(stmt2));
    check_sqlite3(db, "Finalizing sqlite3 statement",
                  sqlite3_finalize(stmt2));
}

void ReadSqlite3(const char                       *sqlite_file,
                 std::vector<alexandria::MolProp> &mp,
                 double                            ref_temperature)
{
#ifdef HAVE_LIBSQLITE3
    std::string                 cas2, csid2;

    sqlite3                    *db   = nullptr;
    sqlite3_stmt               *stmt = nullptr;
    char sql_str[1024];
    const char                 *cas, *csid, *prop, *unit, *source;
    double                      value, error, temperature;
    int                         cidx, rc, nbind, nexp_prop, theory, preferred;
    std::vector<Synonym>        synonyms;
    std::vector<Classes>        classes;
    
    if (nullptr == sqlite_file)
    {
        return;
    }   
    check_sqlite3(nullptr, "Initializing sqlite", sqlite3_initialize());
    check_sqlite3(nullptr, "Opening sqlite database in read-only mode",
                  sqlite3_open_v2(sqlite_file, &db, SQLITE_OPEN_READONLY, nullptr));

    /* Now database is open and everything is Hunky Dory */
    printf("Opened SQLite3 database %s\n", sqlite_file);

    // First get the synonyms out.
    getSynonyms(db, synonyms, mp.size());

    // Now get the classes out.
    getClasses(db, classes, mp.size());

    /* Now present a query statement */
    nexp_prop = 0;
    sprintf(sql_str, "SELECT mol.iupac,mol.cas,mol.csid,pt.prop,pt.unit_text,gp.temperature,gp.value,gp.error, gp.preferred,ds.theory,ds.source FROM molecules as mol,molproperty as gp,proptypes as pt, datasource as ds,phasetype as ph WHERE ((gp.phaseid=ph.phaseid) AND (ph.phase='gas') AND (mol.molid = gp.molid) AND (gp.propid = pt.propid) AND (gp.srcid = ds.srcid) AND (upper(?) = upper(mol.iupac)));");
    check_sqlite3(db, "Preparing sqlite3 statement",
                  sqlite3_prepare_v2(db, sql_str, 1+strlen(sql_str), &stmt, nullptr));

    if (nullptr != debug)
    {
        fprintf(debug, "sql_str = '%s'\nvariable = '%s'\n", sql_str, sqlite3_bind_parameter_name(stmt, 1));
        nbind = sqlite3_bind_parameter_count(stmt);
        fprintf(debug, "%d binding parameter(s) in the statement\n%s\n", nbind, sql_str);
    }
    for (auto mpi = mp.begin(); (mpi < mp.end()); mpi++)
    {
        const std::string molname = mpi->getMolname();
        auto              keyptr  = std::find_if(synonyms.begin(), synonyms.end(),
                                                 [molname](Synonym const &s)
                                                 { return molname.compare(s.molname()) == 0; });
        if (synonyms.end() == keyptr)
        {
            fprintf(stderr, "Warning: missing iupac for %s (%s). Will be ignored.\n",
                    molname.c_str(), mpi->formula().c_str());
        }
        else
        {
            if (nullptr != debug)
            {
                fprintf(debug, "Going to query for '%s'\n", keyptr->iupac().c_str());
            }
            check_sqlite3(db, "Binding text",
                          sqlite3_bind_text(stmt, 1, keyptr->iupac().c_str(), -1, SQLITE_STATIC));
            do
            {
                rc = sqlite3_step(stmt);
                if (SQLITE_ROW == rc)
                {
                    cidx   = 0;
                    const char *iupac2 = (char *)sqlite3_column_text(stmt, cidx++);
                    if (strcasecmp(keyptr->iupac().c_str(), iupac2) != 0)
                    {
                        gmx_fatal(FARGS, "Selected '%s' from database but got '%s'. WTF?!",
                                  keyptr->iupac().c_str(), iupac2);
                    }
                    cas            = (char *)sqlite3_column_text(stmt, cidx++);
                    csid           = (char *)sqlite3_column_text(stmt, cidx++);
                    prop           = (char *)sqlite3_column_text(stmt, cidx++);
                    unit           = (char *)sqlite3_column_text(stmt, cidx++);
                    temperature    = sqlite3_column_double(stmt, cidx++);
                    value          = sqlite3_column_double(stmt, cidx++);
                    error          = sqlite3_column_double(stmt, cidx++);
                    preferred      = sqlite3_column_int(stmt, cidx++);
                    theory         = sqlite3_column_int(stmt, cidx++);
                    source         = (char *)sqlite3_column_text(stmt, cidx++);
                    
                    if (fabs(ref_temperature-temperature) < 0.1)
                    {
                        bool bExp = (0 == theory);
                        if (bExp)
                        {
                            if (preferred)
                            {
                                nexp_prop++;
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
                        }
                        else
                        {
                            alexandria::Experiment calc("gentop", 
                                                        source,
                                                        "-", 
                                                        "unknown", 
                                                        "minimum",
                                                        "unknown", 
                                                        alexandria::JOB_UNKNOWN);
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
                    }
                    const char *iupac = keyptr->iupac().c_str();
                    auto        cptr  = std::find_if(classes.begin(), classes.end(),
                                                     [iupac](Classes const &c)
                                                     { return c.iupac().compare(iupac) == 0; });
                    if (cptr != classes.end())
                    {
                        for (auto c = cptr->classBegin(); c < cptr->classEnd(); ++c)
                        {
                            mpi->AddCategory(*c);
                        }
                    }
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
                else if (nullptr != debug)
                {
                    fprintf(debug, "Done finding rows for %s\n", keyptr->iupac().c_str());
                }
            }
            while (SQLITE_ROW == rc);
            sqlite3_clear_bindings(stmt);
            check_sqlite3(db, "Resetting sqlite3 statement", sqlite3_reset(stmt));
        }
    }
    check_sqlite3(db, "Finalizing sqlite3 statement", sqlite3_finalize(stmt));
    check_sqlite3(nullptr, "Closing sqlite database", sqlite3_close(db));
    check_sqlite3(nullptr, "Shutting down sqlite. Sqlite3 code %d.", sqlite3_shutdown());
    printf("Extracted %d data points at %0.2f (K) from sql database\n", nexp_prop, ref_temperature);
    
#else
    fprintf(stderr, "No support for sqlite3 database in this executable.\n");
    fprintf(stderr, "Please rebuild gromacs with cmake flag -DGMX_SQLITE3=ON set.\n");
#endif
}
