/* 
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
      //! Element name
    std::string _element;
    //! Atomic number corresponding to element
    int _atomic_number;
  public:
    //! Constructor
    ScatteringFactor() {}
    
    //! Destructor
    virtual ~ScatteringFactor() {};
    
    //! Returns the atomic number of this element
    int atomic_number() { return _atomic_number; }
    
    //! Returns the element symbol as a string
    std::string element() { return _element; } 
    
    /*! \brief
     * Virtual function to compute the scattering for a given q
     * \param[in] q scattering vector
     * \return scattering factor
     */
    virtual double calc(double q) = 0;
    
    /*! \brief
     * Virtual function to compute the scattering for a given q
     * \param[in] theta scattering angle
     * \param[in] lambda wave length
     * \return scattering factor
     */
    virtual double calc(double theta, double lambda) = 0;
  };

  /*! \brief
   * Heir of ScatteringFactor, for storing Cromer-Mann parameters.
   */
  class CromerMannSfactor: public ScatteringFactor {
  private:
      //! Parameters that describe the polynomial in the Cromer Mann description
    double _a[4], _b[4], _c;
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
    CromerMannSfactor(std::string element,int atomic_number,
                        double a[4],double b[4],double c);
                        
    //! Destructor
    ~CromerMannSfactor() {}
    
    /*! \brief
     * Function to compute the scattering for a given q
     * \param[in] q scattering vector
     * \return scattering factor
     */
    double calc(double q);
    
    /*! \brief
     * Function to compute the scattering for a given q
     * \param[in] theta scattering angle
     * \param[in] lambda wave length
     * \return scattering factor
     */
    double calc(double theta,double lambda);
  };

  /*! \brief
   * Heir of ScatteringFactor, for storing parameters as polynomial+fourier series.
   */
  class FourierSfactor: public ScatteringFactor {
  private:
      //! Parameters describing the Fourier descrption of the scattering
    double _p[3], _a0, _q0, _qrange;
    
    //! More parameters
    std::vector<double> _a, _b;
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
    FourierSfactor(std::string element,int atomic_number,double q0, double qrange,
                    double p[3], double a0, std::vector<double> a, std::vector<double> b);
                    
    //! Destructor
    ~FourierSfactor() {}
    
    /*! \brief
     * Function to compute the scattering for a given q
     * \param[in] q scattering vector
     * \return scattering factor
     */
    double calc(double q);
    
    /*! \brief
     * Function to compute the scattering for a given q
     * \param[in] theta scattering angle
     * \param[in] lambda wave length
     * \return scattering factor
     */
    double calc(double theta,double lambda);
  };
    
  /*! \brief
   * Class for loading and storing scattering factor parameters for several atom types, 
   * as well as calculating f(q) for an atom. Note that the elements do not have to
   * chemical elements, they may be e.g. coarse grained particles or entire amino acids.
   */
  class ScatteringFactorTable {
  private:
    // Vectors of length number of elements in table
    std::vector<ScatteringFactor*> _cm;
  public:
    /*! \brief
     * Constructor
     * \param[in] datafile  File containing the scattering factors
     */
    ScatteringFactorTable(const char *datafile);
    
    //! Destructor
    ~ScatteringFactorTable() {}
    
    //! Return length of the table
    unsigned int size() { return _cm.size(); }
    
    /*! \brief
     * \param i index in table
     * \return the atomic number corresponding to index
     */
    int atomic_number(int i) { return _cm[i]->atomic_number(); }
    
    /*! \brief
     * Compute scattering factor
     * \param[in] atomic_number  Atomic Number of the element looked for
     * \param[in] q              Scattering vector
     * \return scattering factor
     */
    double calc(int atomic_number,double q);
    
    /*! \brief
     * Compute scattering factor
     * \param[in] element  Element name of the element looked for
     * \param[in] q        Scattering vector
     * \return scattering factor
     */
    double calc(std::string element,double q);
    
    /*! \brief
     * Compute scattering factor
     * \param[in] element  Element name of the element looked for
     * \param[in] q        Scattering vector
     * \return scattering factor
     */
    double calc(char *element,double q) { std::string s(element); 
        return calc(s,q); }
        
    //! Return largest atomic number in the table
    int max_atomic_number();
  };
}

#endif /* _scattering_factors check */
