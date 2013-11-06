/*
 *
 *                This source code is part of
 *
 *                 G   R   O   M   A   C   S
 *
 *          GROningen MAchine for Chemical Simulations
 *
 * Written by David van der Spoel, Erik Lindahl, Berk Hess, and others.
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2012, The GROMACS development team,
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
 */
/*! \file
 * \brief
 * Declares code for calling GROMACS routines from an external program
 *
 * \author David van der Spoel <spoel@xray.bmc.uu.se>
 * \inpublicapi
 * \ingroup module_qmmm
 */
#ifndef GMX_MM_SLAVE_H
#define GMX_MM_SLAVE_H

#ifdef __cplusplus
extern "C" {
#endif

    //! Abstract type for the mm_slave code
    typedef struct gmx_mm_slave *gmx_mm_slave_t;

    //! Create the data structure
    extern gmx_mm_slave_t mm_slave_init(void);
    
    /*! \brief
     * Function to clean up
     * \param[in] gms  the data structure to be freed
     */
    extern gmx_mm_slave_t mm_slave_done(gmx_mm_slave_t gms);

    /*! \brief
     * Function to read a tpr file and store the information in the
     * data structur
     * \param[in] tpr  the file name
     * \param[in] gms  the data structure to be read
     * \return TRUE if successful
     */
    extern gmx_bool mm_slave_read_tpr(const char *tpr, 
                                      gmx_mm_slave_t gms);

    /*! \brief
     * Function to modify the charge of an atom
     * data structure
     * \param[in] gms  the data structure to be modified
     * \param[in] id   the atom id
     * \param[in] q    the new charge
     * \return TRUE if successful
     */
    extern gmx_bool mm_slave_set_q(gmx_mm_slave_t gms,
                                   atom_id id,
                                   double q);

    /*! \brief
     * Function to compute the energy and forces
     * \param[in]  gms the data structure to be modified
     * \param[in]  x   the atomic coordinates for the whole system (MM+QM)
     * \param[out] f   the forces on all atoms
     * \param[out] energy the total MM energy
     * \return TRUE if successful
     */
    extern gmx_bool mm_slave_calc_energy(gmx_mm_slave_t gms, 
                                         const rvec *x,
                                         rvec *f,
                                         double *energy);
    
    
#ifdef __cplusplus
}
#endif

#endif
