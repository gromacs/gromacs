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
 
#ifndef LATEX_UTIL_H
#define LATEX_UTIL_H

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "categories.h"
#include "molprop.h"
#include "molprop_util.h"


namespace alexandria
{

class LongTable
{
    private:
        FILE                    *fp_;
        const char              *font_;
        std::string              caption_;
        std::string              columns_;
        std::string              label_;
        std::vector<std::string> headLines_;
        bool                     bLandscape_;
    public:
        //! Constructor with a file pointer
        LongTable(FILE *fp, bool bLandscape, const char *font);

        //! Constructor with a file name
        LongTable(const char *fn, bool bLandscape);

        //! Destructor
        ~LongTable() {};

        void setCaption(const char *caption) { caption_.assign(caption); }

        void setLabel(const char *label) { label_.assign(label); }

        //! Generate columns entry with first column left aligned and other center
        void setColumns(int nColumns);

        void setColumns(const char *columns) { columns_.assign(columns); }

        void addHeadLine(const std::string &headline) { headLines_.push_back(headline); }

        void printHeader();

        void printFooter();

        void printLine(const std::string &line);

        void printHLine();
};

class ExpData
{
    public:
        double      val_, err_, temp_;
        std::string ref_, conf_, type_, unit_;

        ExpData(double val, double err, double temp, 
                std::string ref, std::string conf, 
                std::string type, std::string unit);
};

class CalcData
{
    public:
        double val_, err_, temp_;
        int    found_;
        CalcData(double val, double err, double temp, int found);
};




}// namespace

#endif
