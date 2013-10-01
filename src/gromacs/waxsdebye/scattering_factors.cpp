/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2013, by the GROMACS development team, led by
 * David van der Spoel, Berk Hess, Erik Lindahl, and including many
 * others, as listed in the AUTHORS file in the top-level source
 * directory and at http://www.gromacs.org.
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
#include <ios>
#include <iostream>
#include <fstream>
#include <algorithm>
#include <math.h>
#include "maths.h"
#include "scattering_factors.h"
#include "gromacs/utility/gmxassert.h"
#include "gromacs/utility/file.h"

using namespace std;

namespace gmx
{

/*! \brief
 * Heir of ScatteringFactor, for storing Cromer-Mann parameters.
 */
class CromerMannSfactor : public ScatteringFactor
{
    private:
        //! Parameters that describe the polynomial in the Cromer Mann description
        double a_[4], b_[4], c_;
    public:
        /*! \brief
         * Constructor
         *
         * \param[in] element        The element symbol
         * \param[in] atomic_number  The atomic number corresponding to this
         * \param[in] a              Cromer Man parameters a
         * \param[in] b              Cromer Man parameters b
         * \param[in] c              Cromer Man parameters c
         */
        CromerMannSfactor(const std::string &element, int atomic_number,
                          double a[4], double b[4], double c);

        //! Destructor
        virtual ~CromerMannSfactor() {}

        /*! \brief
         * Function to compute the scattering for a given q
         * \param[in] q scattering vector
         * \return scattering factor
         */
        virtual double computeScatteringFactor(double q);

        /*! \brief
         * Function to compute the scattering for a given q
         * \param[in] theta scattering angle
         * \param[in] lambda wave length
         * \return scattering factor
         */
        virtual double computeScatteringFactor(double theta, double lambda);
};

/*! \brief
 * Heir of ScatteringFactor, for storing parameters as polynomial+fourier series.
 */
class FourierSfactor : public ScatteringFactor
{
    private:
        //! Parameters describing the Fourier descrption of the scattering
        double p_[3], a0_, q0_, qrange_;

        //! More parameters
        std::vector<double> a_, b_;
    public:
        /*! \brief
         * Constructor
         *
         * \param[in] element        The element symbol
         * \param[in] atomic_number  The atomic number corresponding to this
         * \param[in] q0             Fourier parameters q0
         * \param[in] qrange         Fourier parameters qrange
         * \param[in] p              Fourier parameters p
         * \param[in] a0             Fourier parameters a0
         * \param[in] a              Fourier parameters a
         * \param[in] b              Fourier parameters b
         */
        FourierSfactor(const std::string &element, int atomic_number,
                       double q0, double qrange,
                       double p[3], double a0,
                       std::vector<double> &a,
                       std::vector<double> &b);

        //! Destructor
        virtual ~FourierSfactor() {}

        /*! \brief
         * Function to compute the scattering for a given q
         * \param[in] q scattering vector
         * \return scattering factor
         */
        virtual double computeScatteringFactor(double q);

        /*! \brief
         * Function to compute the scattering for a given q
         * \param[in] theta scattering angle
         * \param[in] lambda wave length
         * \return scattering factor
         */
        virtual double computeScatteringFactor(double theta, double lambda);
};


CromerMannSfactor::CromerMannSfactor(const std::string &element, int atomic_number,
                                     double a[4], double b[4], double c)
{
    int i;

    element_       = element;
    atomic_number_ = atomic_number;
    for (i = 0; (i < 4); i++)
    {
        a_[i] = a[i];
        b_[i] = b[i];
    }
    c_ = c;
}


double CromerMannSfactor::computeScatteringFactor(double theta, double lambda)
{
    double q;

    GMX_RELEASE_ASSERT((lambda <= 0),
                       "CromerMannSfactor::computeScatteringFactor called with lambda <= 0");

    q = 4*M_PI*sin(theta)/lambda; // added 4pi here for q -AB
    return computeScatteringFactor(q);
}


double CromerMannSfactor::computeScatteringFactor(double q)
{
    int    i;
    double cm = c_;
    double q4 = q/(4*M_PI);

    for (i = 0; (i < 4); i++)
    {
        cm += a_[i]*exp(-b_[i]*q4*q4);
    }
    return cm;
}


FourierSfactor::FourierSfactor(const std::string &element, int atomic_number,
                               double q0, double qrange, double p[3], double a0,
                               std::vector<double> &a, std::vector<double> &b)
{
    int i;

    element_       = element;
    atomic_number_ = atomic_number;
    q0_            = q0;
    qrange_        = qrange;
    a0_            = a0;
    for (i = 0; i < 3; i++)
    {
        p_[i] = p[i];
    }
    a_ = a;
    b_ = b;
}


double FourierSfactor::computeScatteringFactor(double theta, double lambda)
{
    double q;

    GMX_RELEASE_ASSERT((lambda <= 0),
                       "FourierSfactor::computeScatteringFactor called with lambda <= 0");

    q = 4*M_PI*sin(theta)/lambda; // added 4pi here for q -AB
    return computeScatteringFactor(q);
}

double FourierSfactor::computeScatteringFactor(double q)
{
    unsigned int i;
    double       s = a0_;

    // add the polynomial
    s += p_[0]*q*q*q + p_[1]*q*q + p_[2]*q;

    // add the fourier series
    for (i = 0; i < a_.size(); i++)
    {
        s += a_[i]*cos(2*M_PI*(i+1)*(q-q0_)/qrange_);
        s += b_[i]*sin(2*M_PI*(i+1)*(q-q0_)/qrange_);
    }
    return s;
}

//! Compare the atomic numbers in two ScatteringFactors
static bool sfactor_comp(ScatteringFactor* a, ScatteringFactor* b)
{
    if (a->atomic_number() < b->atomic_number() )
    {
        return true;
    }
    else
    {
        return false;
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
    std::string         line;
    int                 i;

    //! Open a GROMACS File object.
    gmx::File fp(datafile, "r");

    /* Read the experimental data */
    while (fp.readLine(&line))
    {
        if (sscanf(line.c_str(), "%12s%12d%12lf%12lf%12lf%12lf%12lf%12lf%12lf%12lf%12lf%12lf%12lf%12lf%12lf%12lf%12lf%12lf",
                   elem, &atomic_number, &q0, &qrange, &p[0], &p[1], &p[2], &a0, &x[0], &x[1], &x[2], &x[3], &x[4],
                   &x[5], &x[6], &x[7], &x[8], &x[9]) == 18)
        {
            // Fourier scattering factor with N=5
            element.assign(elem);
            for (i = 0; i < 5; i++)
            {
                fourier_a.push_back(x[i]);
                fourier_b.push_back(x[5+i]);
            }
            if (element.substr(0, 1) != ";")
            {
                ScatteringFactorPointer temp = new FourierSfactor(element, atomic_number, q0, qrange, p, a0, fourier_a, fourier_b);

                fourier_a.clear();
                fourier_b.clear();
                cm_.push_back(temp);
            }
        }
        else if (sscanf(line.c_str(), "%12s%12d%12lf%12lf%12lf%12lf%12lf%12lf%12lf%12lf%12lf",
                        elem, &atomic_number, &a[0], &a[1], &a[2], &a[3],
                        &b[0], &b[1], &b[2], &b[3], &c) == 11)
        {
            // Cromer-Mann scattering factor
            element.assign(elem);
            if (element.substr(0, 1) != ";")
            {
                ScatteringFactorPointer temp = new CromerMannSfactor(element, atomic_number, a, b, c);

                cm_.push_back(temp);
            }
        }
    }
    fp.close();
    sort(cm_.begin(), cm_.end(), sfactor_comp);
}


double ScatteringFactorTable::computeScatteringFactor(int atomic_number, double q)
{
    unsigned int i = 0;
    // cm_ is a vector of scattering_factor:s.
    for (i = 0; (i < cm_.size()); i++)
    {
        // class scattering_factor has a function atomic_number() that returns it's private variable atomic_number_:
        if (cm_[i]->atomic_number() == atomic_number)
        {
            // class scattering_factor has a method calc(), so call that when the right entry is found.
            return cm_[i]->computeScatteringFactor(q);
        }
    }
    return -1;
}


double ScatteringFactorTable::computeScatteringFactor(const std::string &element, double q)
{
    unsigned int i = 0;
    for (i = 0; (i < cm_.size()); i++)
    {
        if (cm_[i]->element() == element) // is it really this easy to compare strings? AB
        {
            return cm_[i]->computeScatteringFactor(q);
        }
    }
    return -1;
}


int ScatteringFactorTable::max_atomic_number()
{
    unsigned int i         = 0;
    int          maxnumber = cm_[0]->atomic_number();
    for (i = 0; (i < cm_.size()); i++)
    {
        if (cm_[i]->atomic_number() > maxnumber)
        {
            maxnumber = cm_[i]->atomic_number();
        }
    }
    return maxnumber;
}

}
