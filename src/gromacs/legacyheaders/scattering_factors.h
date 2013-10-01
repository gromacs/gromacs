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
#ifndef _scattering_factors_h
#define _scattering_factors_h

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include <string>
#include <vector>

namespace gmx {
  /*! \brief
   * Base class for storing scattering factor parameters for an atom type.
   */
  class ScatteringFactor {
  protected:
    std::string _element;
    int _atomic_number;
  public:
    ScatteringFactor() {}
    virtual ~ScatteringFactor() {};
    int atomic_number() { return _atomic_number; }
    std::string element() { return _element; } 
    virtual double calc(double q) = 0;
    virtual double calc(double theta, double lambda) = 0;
  };

  /*! \brief
   * Heir of ScatteringFactor, for storing Cromer-Mann parameters.
   */
  class CromerMannSfactor: public ScatteringFactor {
  private:
    double _a[4], _b[4], _c;
  public:
    CromerMannSfactor(std::string element,int atomic_number,
                        double a[4],double b[4],double c);
    ~CromerMannSfactor() {}
    double calc(double q);
    double calc(double theta,double lambda);
  };

  /*! \brief
   * Heir of ScatteringFactor, for storing parameters as polynomial+fourier series.
   */
  class FourierSfactor: public ScatteringFactor {
  private:
    double _p[3], _a0, _q0, _qrange;
    std::vector<double> _a, _b;
  public:
    FourierSfactor(std::string element,int atomic_number,double q0, double qrange,
                    double p[3], double a0, std::vector<double> a, std::vector<double> b);
    ~FourierSfactor() {}
    double calc(double q);
    double calc(double theta,double lambda);
  };
    
  /*! \brief
   * Class for loading and storing scattering factor parameters for several atom types, 
   * as well as calculating f(q) for an atom.
   */
  class ScatteringFactorTable {
  private:
    // Vectors of length number of elements in table
    std::vector<ScatteringFactor*> _cm;
  public:
    ScatteringFactorTable(const char *datafile);
    ~ScatteringFactorTable() {}
    unsigned int size() { return _cm.size(); }
    int atomic_number(int i) { return _cm[i]->atomic_number(); }
    double calc(int atomic_number,double q);
    double calc(std::string element,double q);
    double calc(char *element,double q) { std::string s(element); 
        return calc(s,q); }
    int max_atomic_number();
  };
}

#endif /* _scattering_factors check */
