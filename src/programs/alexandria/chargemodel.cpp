/*
 * This source file is part of the Alexandria Chemistry Toolkit.
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

#include "chargemodel.h"

#include <map>
#include <string>

namespace alexandria
{

std::map<ChargeModel, const std::string> cmNames =
    {
     { eqdESP_p,    "ESP-p"    },
     { eqdESP_pp,   "ESP-pp"   },
     { eqdESP_pg,   "ESP-pg"   },
     { eqdESP_ps,   "ESP-ps"   },
     { eqdACM_g,    "ACM-g"    },
     { eqdACM_s,    "ACM-s"    },
     { eqdACM_pg,   "ACM-pg"   },
     { eqdACM_ps,   "ACM-ps"   },
     { eqdYang,     "Yang"     },
     { eqdBultinck, "Bultinck" },
     { eqdRappe,    "Rappe"    }
    };

std::map<const std::string, ChargeModel> cmEEM;

ChargeModel name2eemtype(const std::string &name)
{
    if (cmEEM.size() == 0)
    {
        for(auto &k : cmNames)
        {
            cmEEM.emplace(k.second, k.first);
        }
    }
    auto cc = cmEEM.find(name);
    if (cc != cmEEM.end())
    {
        return cc->second;
    }
    return eqdNR;
}

const char *getEemtypeName(ChargeModel eem)
{
    auto cm = cmNames.find(eem);
    if (cm != cmNames.end())
    {
        return cm->second.c_str();
    }
    return nullptr;
}

} // namespace alexandria
