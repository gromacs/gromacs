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
#ifndef GMX_QMMM_QMMM_H
#define GMX_QMMM_QMMM_H

namespace gmx
{

class QmmmHandle;

/*! \brief
 * Container class for interfacing with QM codes.
 *
 * \inpublicapi
 * \ingroup module_qmmm
 */
class Qmmm : public AbstractQmmm
{
    public:
        //! Creates a qmmm object.
        Qmmm(t_commrec *cr, 
	     matrix box, 
	     gmx_mtop_t *mtop, 
	     t_inputrec *ir,
	     t_forcerec *fr);
        virtual ~Qmmm();

	/*! \brief
	 * Fills the MM stuff in in a Qmmm object.
	 *
	 * The MM atoms are
	 * taken from the neighbourlists of the QM atoms. In a QMMM run this
	 * routine should be called at every step, since it updates the MM
	 * elements of the t_QMMMrec struct.  
	 */
	void update(t_commrec *cr,
		    t_forcerec *fr,
		    rvec x[],
		    t_mdatoms *md,
		    matrix box,
		    gmx_localtop_t *top);


        /*! \brief
         * Do the actual QM calculation
         *
	 * Computes the QM forces. This routine makes either function
	 * calls to gmx QM routines (derived from MOPAC7 (semi-emp.) and MPQC
	 * (ab initio)) or generates input files for an external QM package
	 * (listed in QMMMrec.QMpackage). The binary of the QM package is
	 * called by system().
	 * Returns the energy.
         */
        real calculate(t_commrec *cr,
		       rvec x[], 
		       rvec f[],
		       t_forcerec *fr,
		       t_mdatoms *md);


    private:
        class Impl;

        Impl                *_impl;

        friend class Impl;

        // Copy and assign disallowed by base class.
};

} // namespace gmx

#endif
