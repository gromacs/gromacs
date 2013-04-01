/* 
 * $Id: gauss_io.h,v 1.8 2009/02/02 21:11:11 spoel Exp $
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

#ifndef _gauss_io_hpp
#define _gauss_io_hpp

#include "molprop.hpp"

namespace alexandria {

/*! \brief
 * Helper class for extracting atomization energies from Gaussian calcs.
 *
 * This contains data in order to extract enthalpy of
 * formation and Gibbs energy of formation from Thermochemistry methods
 * in Gaussian.
 * 
 * \inlibraryapi
 * \ingroup module_alexandria
 */
class GaussAtomPropVal {
private:
    std::string _element,_method,_desc;
    real _temp;
    real _value;
public:
    //! Constructor initiating all the values stored
    GaussAtomPropVal(std::string element,std::string method, std::string desc,
                     real temp,real value) {
        _element = element; _method == method; _desc = desc;
        _temp = temp; _value = value;
    }

    //! Default destructor
    ~GaussAtomPropVal() {}

    //! Return element name
    std::string GetElement() { return _element; }

    //! Return the method used
    std::string GetMethod() { return _method; }

    //! Return a description of the type of value stored
    std::string GetDesc() { return _desc; }

    //! Return the temperature
    real GetTemp() { return _temp; }

    //! Return the actual value
    real GetValue() { return _value; }
}; 

/*! \brief
 * Class for extracting atomization energies from Gaussian calcs.
 *
 * This contains data in order to extract enthalpy of
 * formation and Gibbs energy of formation from Thermochemistry methods
 * in Gaussian. The method and description must match that in a library file.
 * 
 * \inpublicapi
 * \ingroup module_alexandria
 */
class GaussAtomProp {
private:
    std::vector<GaussAtomPropVal> _gapv;
public:
    //! Default constructor
    GaussAtomProp();

    //! Default destructor
    ~GaussAtomProp() {}
    
    /*! \brief
     * Look up the value corresponding to input variables
     *
     * \param[in] element  From the periodic table
     * \param[in] method   Method used, e.g. G3, exp
     * \param[in] desc     Which type of name is 
     * \param[in] temp     Temperature
     * \param[out] value   The energy value
     * \return 1 on success, 0 otherwise
     * \ingroup module_alexandria
     */
    int GetValue(const char *element,
                 const char *method,
                 const char *desc,
                 double temp,
                 double *value);
};
}

/*! \brief
 * Read a Gaussian log file either using home grown methods or using OpenBabel 
 *
 *
 * \param[in] g98        The gaussian log file, or in case OpenBabel is used anything
 *                       that can be read by OpenBabel
 * \param[out] mpt       The MolProp
 * \param[in] gap        Helper data for reading atomization energies
 * \param[in] bBabel     Whether or not to use the OpenBabel library
 * \param[in] aps        Structure containing atomic data
 * \param[in] pd         The force field information
 * \param[in] molnm      Molecule name to override the one from the filename [ maybe NULL ]
 * \param[in] iupac      IUPAC name to override the one from the filename [ maybe NULL ]
 * \param[in] conformation  Conformation the molecule is in [ maybe NULL ]
 * \param[in] basisset   Basis set used for the calculation [ maybe NULL ]
 * \param[in] maxpot     Maximum number of electrostatic potential data points to store
 * \param[in] bVerbose   Whether or not to write to the terminal during processing
 * \ingroup module_alexandria
 */
void ReadGauss(const char *g98,
               alexandria::MolProp& mpt,
               alexandria::GaussAtomProp &gap,
               gmx_bool bBabel,
               gmx_atomprop_t aps,gmx_poldata_t pd,
               char *molnm,char *iupac,char *conformation,
               char *basisset,
               int maxpot,gmx_bool bVerbose);

/*! \brief
 * Convert the OpenBabel atomtypes to atomtypes corresponding to a force field
 *
 * \param[out] atoms        Atoms structure containing the input and output types 
 * \param[out] symtab       String handling structure
 * \param[in]    forcefield   Name of the desired force field
 * \todo Improve error handling, e.g. in case a non-existing force field is selected.
 * \ingroup module_alexandria
 */
void translate_atomtypes(t_atoms *atoms,t_symtab *tab,const char *forcefield);

#endif
