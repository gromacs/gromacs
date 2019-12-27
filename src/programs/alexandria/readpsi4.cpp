/*
 * This source file is part of the Alexandria program.
 *
 * Copyright (C) 2019
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

#include "gmxpre.h"

#include "readpsi4.h"

#include "gromacs/math/units.h"
#include "gromacs/math/vec.h"
#include "gromacs/mdtypes/md_enums.h"
#include "gromacs/pbcutil/pbc.h"
#include "gromacs/topology/atomprop.h"
#include "gromacs/utility/arraysize.h"
#include "gromacs/utility/exceptions.h"
#include "gromacs/utility/real.h"
#include "gromacs/utility/strconvert.h"
#include "gromacs/utility/stringutil.h"
#include "gromacs/utility/textreader.h"

#include "molprop.h"
#include "molprop_util.h"
#include "poldata.h"
#include "stringutil.h"

namespace alexandria
{

bool readPsi4(const std::string &datafile, MolProp *mp)
{
    try
    {
        gmx::TextReader       tr(datafile);
        std::string           program("psi4");
        std::string           method, basisset, reference, conformation;
        jobType               jobtype = JOB_UNKNOWN;
        std::string           line;
        bool                  inputSection  = false;
        bool                  finalGeometry = false;
        int                   nLongLines    = 0;
        int                   charge        = -6;
        int                   multiplicity  = 0;
        double                energy        = 0;
        rvec                  mu            = { 0, 0, 0 };
        std::vector<CalcAtom> ca;
        while (tr.readLine(&line))
        {
            if (line.find("Input File") != std::string::npos)
            {
                inputSection = true;
            }
            else if (inputSection && line.find("--------") != std::string::npos)
            {
                nLongLines += 1;
                if (nLongLines == 2)
                {
                    inputSection = false;
                }
            }
            else if (inputSection && line.find("molecule") != std::string::npos)
            {
                auto words = gmx::splitString(line);
                if (words.size() == 3)
                {
                    auto wf = words[1].find("mol_");
                    if (wf != std::string::npos)
                    {
                        std::string molnm = words[1].substr(wf+4);
                        std::replace(molnm.begin(), molnm.end(), '_', '-');
                        mp->SetMolname(molnm);
                    }
                    else
                    {
                        mp->SetMolname(words[1]);
                    }
                }
                if (tr.readLine(&line))
                {
                    words = gmx::splitString(line);
                    if (words.size() == 2)
                    {
                        charge       = gmx::intFromString(words[0].c_str());
                        multiplicity = gmx::intFromString(words[1].c_str());
                    }
                }
            }
            else if (inputSection && line.find("basis") != std::string::npos)
            {
                auto words = gmx::splitString(line);
                if (words.size() == 2)
                {
                    basisset = words[1];
                }
            }
            else if (inputSection && line.find("optimize") != std::string::npos)
            {
                auto words = gmx::splitDelimitedString(line, '\'');
                if (words.size() == 3)
                {
                    method = words[1];
                }
                jobtype = JOB_OPT;
                conformation.assign("minimum");
            }
            else if (inputSection && line.find("energy") != std::string::npos)
            {
                auto words = gmx::splitDelimitedString(line, '\'');
                if (words.size() == 3)
                {
                    method = words[1];
                }
                jobtype = JOB_SP;
                conformation.assign("excited");
            }
            else if (line.find("Final optimized geometry") != std::string::npos)
            {
                finalGeometry = true;
            }
            else if (finalGeometry &&
                     line.find("Geometry (in Angstrom)") != std::string::npos &&
                     line.find("charge") != std::string::npos &&
                     line.find("multiplicity") != std::string::npos)
            {
                // Skip empty line
                if (tr.readLine(&line))
                {
                    bool cont = true;
                    do
                    {
                        if (tr.readLine(&line))
                        {
                            auto words = gmx::splitString(line);
                            if (words.size() == 4)
                            {
                                CalcAtom myatom(words[0], words[0], ca.size()+1);
                                myatom.SetUnit("Angstrom");
                                myatom.SetCoords(gmx::doubleFromString(words[1].c_str()),
                                                 gmx::doubleFromString(words[2].c_str()),
                                                 gmx::doubleFromString(words[3].c_str()));
                                ca.push_back(std::move(myatom));
                            }
                            else
                            {
                                cont = false;
                            }
                        }
                    }
                    while (cont);
                }
            }
            else if (jobtype == JOB_SP &&
                     line.find("Geometry (in Angstrom)") != std::string::npos &&
                     line.find("charge") != std::string::npos &&
                     line.find("multiplicity") != std::string::npos)
            {
                // Skip empty lines
                bool cont = true;
                for (int k = 0; cont && k < 3; k++)
                {
                    cont = tr.readLine(&line);
                }
                if (cont)
                {
                    do
                    {
                        cont = tr.readLine(&line);
                        if (cont)
                        {
                            auto words = gmx::splitString(line);
                            if (words.size() == 5)
                            {
                                CalcAtom myatom(words[0], words[0], ca.size()+1);
                                myatom.SetUnit("Angstrom");
                                myatom.SetCoords(gmx::doubleFromString(words[1].c_str()),
                                                 gmx::doubleFromString(words[2].c_str()),
                                                 gmx::doubleFromString(words[3].c_str()));
                                ca.push_back(std::move(myatom));
                            }
                            else
                            {
                                cont = false;
                            }
                        }
                    }
                    while (cont);
                }
            }
            else if (jobtype == JOB_SP &&
                     line.find("Final double-hybrid DFT total energy") != std::string::npos)
            {
                auto words = gmx::splitString(line);
                if (words.size() == 7)
                {
                    energy = convert2gmx(my_atof(words[6].c_str(), "energy"),
                                         eg2cHartree);
                }
            }
            else if (jobtype == JOB_OPT &&
                     line.find("Final energy is") != std::string::npos)
            {
                auto words = gmx::splitString(line);
                if (words.size() == 4)
                {
                    energy = convert2gmx(my_atof(words[3].c_str(), "energy"),
                                         eg2cHartree);
                }
            }
            else if (jobtype == JOB_OPT &&
                     line.find("Dipole Moment: [D]") != std::string::npos)
            {
                bool cont = tr.readLine(&line);
                if (cont)
                {
                    auto words = gmx::splitString(line);
                    if (words.size() == 8)
                    {
                        mu[XX] = my_atof(words[1].c_str(), "mu[XX]");
                        mu[YY] = my_atof(words[3].c_str(), "mu[YY]");
                        mu[ZZ] = my_atof(words[5].c_str(), "mu[ZZ]");
                    }
                }
            }
        }
        Experiment e(program, method, basisset, reference, conformation,
                     datafile, jobtype);
        for (auto &myatom : ca)
        {
            e.AddAtom(myatom);
        }
        e.AddEnergy(MolecularEnergy("HF", unit2string(eg2cKj_Mole),
                                    0, epGAS, energy, 0.0));
        double dip = norm(mu);
        if (dip > 0)
        {
            e.AddDipole(MolecularDipole("electronic", "Debye", 0,
                                        mu[XX], mu[YY], mu[ZZ], dip, 0));
        }
        mp->AddExperiment(std::move(e));
        mp->SetCharge(charge);
        mp->SetMultiplicity(multiplicity);
        return true;
    }
    GMX_CATCH_ALL_AND_EXIT_WITH_FATAL_ERROR;

    return false;
}

} // namespace alexandria
