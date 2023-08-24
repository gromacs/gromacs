/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2010- The GROMACS Authors
 * and the project initiators Erik Lindahl, Berk Hess and David van der Spoel.
 * Consult the AUTHORS/COPYING files and https://www.gromacs.org for details.
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
 * https://www.gnu.org/licenses, or write to the Free Software Foundation,
 * Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA.
 *
 * If you want to redistribute modifications to GROMACS, please
 * consider that scientific software is very special. Version
 * control is crucial - bugs must be traceable. We will be happy to
 * consider code for inclusion in the official distribution, but
 * derived work must not be called official GROMACS. Details are found
 * in the README & COPYING files - if they are missing, get the
 * official version at https://www.gromacs.org.
 *
 * To help us fund GROMACS development, we humbly ask that you cite
 * the research papers on the package. Check out https://www.gromacs.org.
 */
/*! \internal \file
 * \brief
 * Declares helper functions for assigning isotops from atom names
 *
 * \author Alexey Shvetsov <alexxyum@gmail.com>
 * \ingroup module_trajectoryanalysis
 */

#include "isotope.h"

#include "gromacs/utility/exceptions.h"
#include "gromacs/utility/stringutil.h"

namespace gmx
{

//! Convert isotope name from string to enum
Isotope getIsotopeFromString(const std::string& isotope)
{
    if (isotope == "H")
    {
        return Isotope::H;
    }
    else if (isotope == "D")
    {
        return Isotope::D;
    }
    else if (isotope == "He")
    {
        return Isotope::He;
    }
    else if (isotope == "Li")
    {
        return Isotope::Li;
    }
    else if (isotope == "Be")
    {
        return Isotope::Be;
    }
    else if (isotope == "B")
    {
        return Isotope::B;
    }
    else if (isotope == "C")
    {
        return Isotope::C;
    }
    else if (isotope == "N")
    {
        return Isotope::N;
    }
    else if (isotope == "O")
    {
        return Isotope::O;
    }
    else if (isotope == "F")
    {
        return Isotope::F;
    }
    else if (isotope == "Ne")
    {
        return Isotope::Ne;
    }
    else if (isotope == "Na")
    {
        return Isotope::Na;
    }
    else if (isotope == "Mg")
    {
        return Isotope::Mg;
    }
    else if (isotope == "Al")
    {
        return Isotope::Al;
    }
    else if (isotope == "Si")
    {
        return Isotope::Si;
    }
    else if (isotope == "P")
    {
        return Isotope::P;
    }
    else if (isotope == "S")
    {
        return Isotope::S;
    }
    else if (isotope == "Cl")
    {
        return Isotope::Cl;
    }
    else if (isotope == "Ar")
    {
        return Isotope::Ar;
    }
    else if (isotope == "K")
    {
        return Isotope::K;
    }
    else if (isotope == "Ca")
    {
        return Isotope::Ca;
    }
    else if (isotope == "Sc")
    {
        return Isotope::Sc;
    }
    else if (isotope == "Ti")
    {
        return Isotope::Ti;
    }
    else if (isotope == "V")
    {
        return Isotope::V;
    }
    else if (isotope == "Cr")
    {
        return Isotope::Cr;
    }
    else if (isotope == "Mn")
    {
        return Isotope::Mn;
    }
    else if (isotope == "Fe")
    {
        return Isotope::Fe;
    }
    else if (isotope == "Co")
    {
        return Isotope::Co;
    }
    else if (isotope == "Ni")
    {
        return Isotope::Ni;
    }
    else if (isotope == "Cu")
    {
        return Isotope::Cu;
    }
    else if (isotope == "Zn")
    {
        return Isotope::Zn;
    }
    else if (isotope == "Ga")
    {
        return Isotope::Ga;
    }
    else if (isotope == "Ge")
    {
        return Isotope::Ge;
    }
    else if (isotope == "As")
    {
        return Isotope::As;
    }
    else if (isotope == "Se")
    {
        return Isotope::Se;
    }
    else if (isotope == "Br")
    {
        return Isotope::Br;
    }
    else if (isotope == "Kr")
    {
        return Isotope::Kr;
    }
    else if (isotope == "Rb")
    {
        return Isotope::Rb;
    }
    else if (isotope == "Sr")
    {
        return Isotope::Sr;
    }
    else if (isotope == "Y")
    {
        return Isotope::Y;
    }
    else if (isotope == "Zr")
    {
        return Isotope::Zr;
    }
    else if (isotope == "Nb")
    {
        return Isotope::Nb;
    }
    else if (isotope == "Mo")
    {
        return Isotope::Mo;
    }
    else if (isotope == "Tc")
    {
        return Isotope::Tc;
    }
    else if (isotope == "Ru")
    {
        return Isotope::Ru;
    }
    else if (isotope == "Rh")
    {
        return Isotope::Rh;
    }
    else if (isotope == "Pd")
    {
        return Isotope::Pd;
    }
    else if (isotope == "Ag")
    {
        return Isotope::Ag;
    }
    else if (isotope == "Cd")
    {
        return Isotope::Cd;
    }
    else if (isotope == "In")
    {
        return Isotope::In;
    }
    else if (isotope == "Sn")
    {
        return Isotope::Sn;
    }
    else if (isotope == "Sb")
    {
        return Isotope::Sb;
    }
    else if (isotope == "Te")
    {
        return Isotope::Te;
    }
    else if (isotope == "I")
    {
        return Isotope::I;
    }
    else if (isotope == "Xe")
    {
        return Isotope::Xe;
    }
    else if (isotope == "Cs")
    {
        return Isotope::Cs;
    }
    else if (isotope == "Ba")
    {
        return Isotope::Ba;
    }
    else if (isotope == "La")
    {
        return Isotope::La;
    }
    else if (isotope == "Ce")
    {
        return Isotope::Ce;
    }
    else if (isotope == "Pr")
    {
        return Isotope::Pr;
    }
    else if (isotope == "Nd")
    {
        return Isotope::Nd;
    }
    else if (isotope == "Pm")
    {
        return Isotope::Pm;
    }
    else if (isotope == "Sm")
    {
        return Isotope::Sm;
    }
    else if (isotope == "Eu")
    {
        return Isotope::Eu;
    }
    else if (isotope == "Gd")
    {
        return Isotope::Gd;
    }
    else if (isotope == "Tb")
    {
        return Isotope::Tb;
    }
    else if (isotope == "Dy")
    {
        return Isotope::Dy;
    }
    else if (isotope == "Ho")
    {
        return Isotope::Ho;
    }
    else if (isotope == "Er")
    {
        return Isotope::Er;
    }
    else if (isotope == "Tm")
    {
        return Isotope::Tm;
    }
    else if (isotope == "Yb")
    {
        return Isotope::Yb;
    }
    else if (isotope == "Lu")
    {
        return Isotope::Lu;
    }
    else if (isotope == "Hf")
    {
        return Isotope::Hf;
    }
    else if (isotope == "Ta")
    {
        return Isotope::Ta;
    }
    else if (isotope == "W")
    {
        return Isotope::W;
    }
    else if (isotope == "Re")
    {
        return Isotope::Re;
    }
    else if (isotope == "Os")
    {
        return Isotope::Os;
    }
    else if (isotope == "Ir")
    {
        return Isotope::Ir;
    }
    else if (isotope == "Pt")
    {
        return Isotope::Pt;
    }
    else if (isotope == "Au")
    {
        return Isotope::Au;
    }
    else if (isotope == "Hg")
    {
        return Isotope::Hg;
    }
    else if (isotope == "Tl")
    {
        return Isotope::Tl;
    }
    else if (isotope == "Pb")
    {
        return Isotope::Pb;
    }
    else if (isotope == "Bi")
    {
        return Isotope::Bi;
    }
    else if (isotope == "Po")
    {
        return Isotope::Po;
    }
    else if (isotope == "At")
    {
        return Isotope::At;
    }
    else if (isotope == "Rn")
    {
        return Isotope::Rn;
    }
    else if (isotope == "Fr")
    {
        return Isotope::Fr;
    }
    else if (isotope == "Ra")
    {
        return Isotope::Ra;
    }
    else if (isotope == "Ac")
    {
        return Isotope::Ac;
    }
    else if (isotope == "Th")
    {
        return Isotope::Th;
    }
    else if (isotope == "Pa")
    {
        return Isotope::Pa;
    }
    else if (isotope == "U")
    {
        return Isotope::U;
    }
    else if (isotope == "Np")
    {
        return Isotope::Np;
    }
    else if (isotope == "Pu")
    {
        return Isotope::Pu;
    }
    else if (isotope == "Am")
    {
        return Isotope::Am;
    }
    else if (isotope == "Cm")
    {
        return Isotope::Cm;
    }
    else if (isotope == "Bk")
    {
        return Isotope::Bk;
    }
    else if (isotope == "Cf")
    {
        return Isotope::Cf;
    }
    else if (isotope == "Es")
    {
        return Isotope::Es;
    }
    else if (isotope == "Fm")
    {
        return Isotope::Fm;
    }
    else if (isotope == "Md")
    {
        return Isotope::Md;
    }
    else if (isotope == "No")
    {
        return Isotope::No;
    }
    else if (isotope == "Lr")
    {
        return Isotope::Lr;
    }
    else if (isotope == "Rf")
    {
        return Isotope::Rf;
    }
    else if (isotope == "Db")
    {
        return Isotope::Db;
    }
    else if (isotope == "Sg")
    {
        return Isotope::Sg;
    }
    else if (isotope == "Bh")
    {
        return Isotope::Bh;
    }
    else if (isotope == "Hs")
    {
        return Isotope::Hs;
    }
    else if (isotope == "Mt")
    {
        return Isotope::Mt;
    }
    else
    {
        std::string message =
                formatString("Isotope '%s' not found in scattering database\n", isotope.c_str());
        GMX_THROW(InconsistentInputError(message));
    }
}

//! Reads element list and returns vector  of isotopes
std::vector<Isotope> getIsotopes(const t_atoms* atoms)
{
    std::vector<Isotope> isotopes(atoms->nr);
    for (int i = 0; i < atoms->nr; ++i)
    {
        std::string elem = atoms->atom[i].elem;
        if (elem.empty())
        {
            continue;
        }
        isotopes.push_back(getIsotopeFromString(elem));
    }
    return isotopes;
}

} // namespace gmx
