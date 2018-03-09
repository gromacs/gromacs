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
 
 
#include "latex_util.h"

namespace alexandria
{

LongTable::LongTable(FILE *fp, bool bLandscape, const char *font)
{
    fp_         = fp;
    bLandscape_ = bLandscape;
    font_       = font;
    if (nullptr == fp_)
    {
        GMX_THROW(gmx::FileIOError("File not open"));
    }
}

LongTable::LongTable(const char *fn, bool bLandscape)
{
    fp_         = fopen(fn, "w");
    bLandscape_ = bLandscape;
    if (nullptr == fp_)
    {
        GMX_THROW(gmx::FileIOError("Could not open file"));
    }
}

void LongTable::setColumns(int nColumns)
{
    columns_.assign("l");
    for (int i = 1; (i < nColumns); i++)
    {
        columns_.append("c");
    }
}

void LongTable::printHeader()
{
    if (bLandscape_)
    {
        fprintf(fp_, "\\begin{landscape}\n");
    }
    if (nullptr != font_)
    {
        fprintf(fp_, "\\begin{%s}\n", font_);
    }
    fprintf(fp_, "\\begin{spacing}{1}\n");
    fprintf(fp_, "\\begin{longtable}{%s}\n", columns_.c_str());
    fprintf(fp_, "\\caption{%s}\n",          caption_.c_str());
    fprintf(fp_, "\\label{%s}\\\\\n",        label_.c_str());
    printHLine();
    for (unsigned int i = 0; (i < headLines_.size()); i++)
    {
        fprintf(fp_, "%s\\\\\n", headLines_[i].c_str());
    }
    printHLine();
    fprintf(fp_, "\\endfirsthead\n");
    printHLine();
    for (unsigned int i = 0; (i < headLines_.size()); i++)
    {
        fprintf(fp_, "%s\\\\\n", headLines_[i].c_str());
    }
    printHLine();
    fprintf(fp_, "\\endhead\n");
    printHLine();
    fprintf(fp_, "\\endfoot\n");
}

void LongTable::printFooter()
{
    fprintf(fp_, "\\end{longtable}\n");
    fprintf(fp_, "\\end{spacing}\n");
    if (nullptr != font_)
    {
        fprintf(fp_, "\\end{%s}\n", font_);
    }
    if (bLandscape_)
    {
        fprintf(fp_, "\\end{landscape}\n");
    }
    fflush(fp_);
}

void LongTable::printLine(std::string line)
{
    fprintf(fp_, "%s\\\\\n", line.c_str());
}

void LongTable::printHLine()
{
    fprintf(fp_, "\\hline\n");
}

ExpData::ExpData(double val, double err, 
                 double temp, std::string ref, 
                 std::string conf, std::string type, 
                 std::string unit)
    :
        val_(val),
        err_(err), 
        temp_(temp), 
        ref_(ref), 
        conf_(conf), 
        type_(type), 
        unit_(unit)
{};

CalcData::CalcData(double val, double err, double temp, int found)
    :
        val_(val), 
        err_(err), 
        temp_(temp), 
        found_(found)
{};

} //namespace

