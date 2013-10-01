/* -*- mode: c; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; c-file-style: "stroustrup"; -*-
 * $Id: gentop.c,v 1.26 2009/05/20 10:48:03 spoel Exp $
 *
 *                This source code is part of
 *
 *                 G   R   O   M   A   C   S
 *
 *          GROningen MAchine for Chemical Simulations
 *
 *                        VERSION 4.0.99
 * Written by David van der Spoel, Erik Lindahl, Berk Hess, and others.
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2008, The GROMACS development team,
 * check out http://www.gromacs.org for more information.

 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 *
 * If you want to redistribute modifications, please consider that
 * scientific software is very special. Version control is crucial -
 * bugs must be traceable. We will be happy to consider code for
 * inclusion in the official distribution, but derived work must not
 * be called official GROMACS. Details are found in the README & COPYING
 * files - if they are missing, get the official version at www.gromacs.org.
 *
 * To help us fund GROMACS development, we humbly ask that you cite
 * the papers on the package - you can find them in the top README file.
 *
 * For more info, check our website at http://www.gromacs.org
 *
 * And Hey:
 * Groningen Machine for Chemical Simulation
 */
#include <ios>
#include <iostream>
#include <fstream>
#include <algorithm>
#include <math.h>
#include "scattering_factors.h"
#include "gmx_fatal.h"


using namespace std;

//! trim a std::string from start
//! Needs fixing or removing
static inline std::string &ltrim(std::string &s)
{
    // s.erase(s.begin(), std::find_if(s.begin(), s.end(), std::not1(std::ptr_fun<int, int>(std::isspace))));
    return s;
}

namespace gmx
{

CromerMannSfactor::CromerMannSfactor(std::string element, int atomic_number,
                                     double a[4], double b[4], double c)
{
    int i;

    _element       = element;
    _atomic_number = atomic_number;
    for (i = 0; (i < 4); i++)
    {
        _a[i] = a[i];
        _b[i] = b[i];
    }
    _c = c;
}


double CromerMannSfactor::calc(double theta, double lambda)
{
    double q;

    if (lambda <= 0)
    {
        cerr << "CromerMannSfactor::calc called with lambda " << lambda << endl;
        return 0;
    }
    q = 4*M_PI*sin(theta)/lambda; // added 4pi here for q -AB
    return calc(q);
}


double CromerMannSfactor::calc(double q)
{
    int    i;
    double cm = _c;
    double q4 = q/(4*M_PI);

    for (i = 0; (i < 4); i++)
    {
        cm += _a[i]*exp(-_b[i]*q4*q4);
    }
    return cm;
}


FourierSfactor::FourierSfactor(std::string element, int atomic_number,
                               double q0, double qrange, double p[3], double a0,
                               std::vector<double> a, std::vector<double> b)
{
    int i;

    _element       = element;
    _atomic_number = atomic_number;
    _q0            = q0;
    _qrange        = qrange;
    _a0            = a0;
    for (i = 0; i < 3; i++)
    {
        _p[i] = p[i];
    }
    _a = a;
    _b = b;
}


double FourierSfactor::calc(double theta, double lambda)
{
    double q;

    if (lambda <= 0)
    {
        cerr << "FourierSfactor::calc called with lambda " << lambda << endl;
        return 0;
    }
    q = 4*M_PI*sin(theta)/lambda; // added 4pi here for q -AB
    return calc(q);
}

double FourierSfactor::calc(double q)
{
    unsigned int i;
    double       s = _a0;

    // add the polynomial
    s += _p[0]*q*q*q + _p[1]*q*q + _p[2]*q;

    // add the fourier series
    for (i = 0; i < _a.size(); i++)
    {
        s += _a[i]*cos(2*M_PI*(i+1)*(q-_q0)/_qrange);
        s += _b[i]*sin(2*M_PI*(i+1)*(q-_q0)/_qrange);
    }
    return s;
}


static bool sfactor_comp(ScatteringFactor* a, ScatteringFactor* b)
{
    if (a->atomic_number() < b->atomic_number() )
    {
        return 1;
    }
    else
    {
        return 0;
    }
}


ScatteringFactorTable::ScatteringFactorTable(const char *datafile)
{
    std::string         element;
    char                elem[32];
    int                 atomic_number;
    // Cromer-Mann parameters:
    double              a[4], b[4], c;
    // Fourier parameters:
    double              p[3], a0, q0, qrange;
    std::vector<double> fourier_a, fourier_b;
    //
    double              x[32]; // x[] is a buffer array
    ifstream            ifs;
    std::string         line;
    int                 i;

    /* Read the experimental data */
    ifs.open(datafile, ios_base::in);
    if (ifs.is_open())
    {
        while (!ifs.eof())
        {
            // read a line:
            std::getline(ifs, line);
            if (sscanf(line.c_str(), "%s%d%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf",
                       elem, &atomic_number, &q0, &qrange, &p[0], &p[1], &p[2], &a0, &x[0], &x[1], &x[2], &x[3], &x[4],
                       &x[5], &x[6], &x[7], &x[8], &x[9]) == 18)
            {
                // Fourier scattering factor with N=5
                element.assign(elem);
                ltrim(element);
                for (i = 0; i < 5; i++)
                {
                    fourier_a.push_back(x[i]);
                    fourier_b.push_back(x[5+i]);
                }
                FourierSfactor* temp_pointer = new FourierSfactor(element, atomic_number, q0, qrange, p, a0, fourier_a, fourier_b);
                fourier_a.clear();
                fourier_b.clear();
                if (element.substr(0, 1) != ";")
                {
                    _cm.push_back(temp_pointer);
                }
            }
            else if (sscanf(line.c_str(), "%s%d%lf%lf%lf%lf%lf%lf%lf%lf%lf",
                            elem, &atomic_number, &a[0], &a[1], &a[2], &a[3],
                            &b[0], &b[1], &b[2], &b[3], &c) == 11)
            {
                // Cromer-Mann scattering factor
                element.assign(elem);
                ltrim(element);
                CromerMannSfactor* temp_pointer = new CromerMannSfactor(element, atomic_number, a, b, c);
                if (element.substr(0, 1) != ";")
                {
                    _cm.push_back(temp_pointer);
                }
            }
        }
    }
    else
    {
        gmx_fatal(FARGS, "error opening file %s.", datafile);
    }
    ifs.close();
    sort(_cm.begin(), _cm.end(), sfactor_comp);
}


double ScatteringFactorTable::calc(int atomic_number, double q)
{
    unsigned int i = 0;
    // _cm is a vector of scattering_factor:s.
    for (i = 0; (i < _cm.size()); i++)
    {
        // class scattering_factor has a function atomic_number() that returns it's private variable _atomic_number:
        if (_cm[i]->atomic_number() == atomic_number)
        {
            // class scattering_factor has a method calc(), so call that when the right entry is found.
            return _cm[i]->calc(q);
        }
    }
    return -1;
}


double ScatteringFactorTable::calc(string element, double q)
{
    unsigned int i = 0;
    for (i = 0; (i < _cm.size()); i++)
    {
        if (_cm[i]->element() == element) // is it really this easy to compare strings? AB
        {
            return _cm[i]->calc(q);
        }
    }
    return -1;
}


int ScatteringFactorTable::max_atomic_number()
{
    unsigned int i         = 0;
    int          maxnumber = _cm[0]->atomic_number();
    for (i = 0; (i < _cm.size()); i++)
    {
        if (_cm[i]->atomic_number() > maxnumber)
        {
            maxnumber = _cm[i]->atomic_number();
        }
    }
    return maxnumber;
}

}
