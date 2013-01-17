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
 * Declares gmx::Qmmm and related classes.
 *
 * \author David van der Spoel <spoel@xray.bmc.uu.se>
 * \inpublicapi
 * \ingroup module_qmmm
 */
#ifndef GMX_QMMM_GAUSSIAN_H
#define GMX_QMMM_GAUSSIAN_H

#include "../qmmm.h"

namespace gmx
{

/*! \brief
 * Base class for interfacing with Gaussian.
 *
 * \inpublicapi
 * \ingroup module_qmmm
 */
  class QmmmGaussian : public QMSystem
  {
    public:
        //! Creates a qmmm object.
        QmmmGaussian(const t_commrec *cr,
                     const matrix box,
                     const gmx_mtop_t *mtop,
                     const t_inputrec *ir,
                     const t_forcerec *fr);
        ~QmmmGaussian();

        /*! \brief
         * Fills the MM stuff in in a Qmmm object.
         *
         * The MM atoms are taken from the neighbourlists of the QM atoms -
         * depending on which QM method is used.
         * In a QMMM run this
         * routine should be called at every step, since it updates the MM
         * elements of the object.
         */
        void update(const t_commrec *cr,
                    const t_forcerec *fr,
                    const rvec x[],
                    const t_mdatoms *md,
                    const matrix box,
                    const gmx_localtop_t *top);


        /*! \brief
         * Do the actual QM calculation
         *
         * Computes the QM forces. The binary of the Gaussian package is
         * called by system().
         * Returns the energy.
         */
        real calculate(const t_commrec *cr,
                       const rvec x[],
                       rvec f[],
                       const t_forcerec *fr,
                       const t_mdatoms *md);

  };

} // namespace gmx

#endif
