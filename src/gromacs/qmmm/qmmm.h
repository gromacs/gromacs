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
 * Declares classes for:
 * Hybrid quantum/classical systems (QMMM and ONIOM).
 * Base class for interface to external QM code.
 *
 * \author David van der Spoel <spoel@xray.bmc.uu.se>
 * \author Lee-Ping Wang <leeping@stanford.edu>
 * \inpublicapi
 * \ingroup module_qmmm
 */
#ifndef GMX_QMMM_QMMM_H
#define GMX_QMMM_QMMM_H

#include "network.h"
#include "typedefs.h"
#include <vector>
//#include "types/qmmmrec.h"

namespace gmx
{

/*! \brief
 * Base class for interface to QM software
 *
 * The QM/MM class will contain one of these, while the
 * ONIOM class can contain multiples of these.
 *
 * \inpublicapi
 * \ingroup module_qmmm
 */
class QMSystem
{
    protected:
        int                               nrQMatoms_;      /* total nr of QM atoms              */
        std::vector<std::vector<real> >   xQM_;            /* shifted to center of box          */
        std::vector<int>                  indexQM_;        /* atom i = atom indexQM[i] in mdrun */
        std::vector<int>                  atomicnumberQM_; /* atomic numbers of QM atoms        */
        std::vector<real>                 QMcharges_;      /* atomic charges of QM atoms(ONIOM) */
        std::vector<int>                  shiftQM_;
        int                               QMcharge_;       /* charge of the QM system           */
        int                               multiplicity_;   /* multipicity (no of unpaired eln)  */
        int                               QMmethod_;       /* see enums.h for all methods       */
        int                               QMbasis_;        /* see enums.h for all bases         */
        int                               nelectrons_;     /* total number of elecs in QM region*/
        gmx_bool                          bTS_;            /* Optimize a TS, only steep, no md  */
        gmx_bool                          bOPT_;           /* Optimize QM subsys, only steep, no md  */
        std::vector<gmx_bool>             frontatoms_;     /* qm atoms on the QM side of a QM-MM bond */

    public:
        QMSystem(int         grpnr,
                 int         nr,
                 int        *atomarray,
                 gmx_mtop_t *mtop,
                 t_inputrec *ir);
        ~QMSystem();
        /*
           virtual real execute(t_commrec *cr,
                 t_forcerec *fr,
                 t_MMrec *mm,
                 rvec f[],
                 rvec fshift[]);
         */
};

/*! \brief
 * Base class for hybrid QM/classical simulation.
 *
 * Contains variables that are common to both QM/classical
 * frameworks (ONIOM and QM/MM).  Does not implement the functions
 * "update" and "calculate" as they are specific to the method.
 *
 * \inpublicapi
 * \ingroup module_qmmm
 */
class HybridQuantumClassical
{
    protected:
        //t_MMrec       *mm_;                    // MM record containing information on MM atoms.
    public:
        HybridQuantumClassical(const t_commrec  *cr,
                               const gmx_mtop_t *mtop,
                               const t_inputrec *ir,
                               const t_forcerec *fr);
        virtual ~HybridQuantumClassical();

        /*! \brief
         * Fills the MM stuff before calling the QM calculation.
         */
        virtual void update(const t_commrec      *cr,
                            const t_forcerec     *fr,
                            const rvec            x[],
                            const t_mdatoms      *md,
                            const matrix          box,
                            const gmx_localtop_t *top);


        /*! \brief
         * Do the actual QM calculation and return the energy.
         */
        virtual real calculate(const t_commrec  *cr,
                               const rvec        x[],
                               rvec              f[],
                               const t_forcerec *fr,
                               const t_mdatoms  *md);

};

/*! \brief
 *
 * Class for performing QM/MM calculations.
 *
 * In QM/MM, there is one QM subsystem and one MM subsystem.
 * They interact via the QM/MM interaction Hamiltonian.
 *
 * The interaction Hamiltonian consists of empirical VdW interactions
 * (these can either be computed by Gromacs or by the QM code) and
 * electrostatic interactions.
 *
 * Implements the methods "update" and "calculate".
 *
 * \inpublicapi
 * \ingroup module_qmmm
 */
class QMMM : public HybridQuantumClassical
{
    public:
        QMMM(const t_commrec  *cr,
             const gmx_mtop_t *mtop,
             const t_inputrec *ir,
             const t_forcerec *fr);
        ~QMMM();

        /*! \brief
         * Fills the MM stuff in a Qmmm object.
         */
        void update(const t_commrec      *cr,
                    const t_forcerec     *fr,
                    const rvec            x[],
                    const t_mdatoms      *md,
                    const matrix          box,
                    const gmx_localtop_t *top);


        /*! \brief
         * Do the actual QM calculation
         */
        real calculate(const t_commrec  *cr,
                       const rvec        x[],
                       rvec              f[],
                       const t_forcerec *fr,
                       const t_mdatoms  *md);
};

/*! \brief
 *
 * Class for performing ONIOM calculations.
 *
 * In ONIOM, there are multiple nexted "layers" of the system, with the smallest
 * layers corresponding to the highest level of theory.  The outermost layer
 * is typically MM.  Each layer includes the atoms of inner layers, such that
 * an outer layer is "AB" while the inner layer is "B".
 *
 * Taken from the Gaussian technical note http://www.gaussian.com/g_whitepap/oniom_technote.htm,
 * the energy for a two-layer system is defined as:
 *
 * EONIOM = Elow(R) + Ehigh(SM) – Elow(SM)
 *
 * where:
 * R stands for Real System (all the atoms),
 * SM stands for Small Model (the high-accuracy region),
 * Elow stands for low level theory (e.g. MM), and
 * Ehigh stands for high level theory (e.g. QM).
 *
 * In three layer ONIOM the energy is approximated as:
 *
 * EONIOM = Elow(R) + Emedium(IM) + Ehigh(SM) – Elow(IM) – Emedium(SM)
 *
 * In contrast to QM/MM there is no explicit interaction between levels of theory.
 * Rather, the low level theory is responsible for all of the interactions.
 *
 * Implements the methods "update" and "calculate".
 *
 * \inpublicapi
 * \ingroup module_qmmm
 */
class ONIOM : public HybridQuantumClassical
{
    private:
        std::vector<QMSystem>        qms;              /* Contains multiple interfaces to QM software. */

    public:
        ONIOM(const t_commrec  *cr,
              const gmx_mtop_t *mtop,
              const t_inputrec *ir,
              const t_forcerec *fr);
        ~ONIOM();

        /*! \brief
         * Fills the MM stuff.
         */
        void update(const t_commrec      *cr,
                    const t_forcerec     *fr,
                    const rvec            x[],
                    const t_mdatoms      *md,
                    const matrix          box,
                    const gmx_localtop_t *top);


        /*! \brief
         * Do the actual QM calculation
         */
        real calculate(const t_commrec  *cr,
                       const rvec        x[],
                       rvec              f[],
                       const t_forcerec *fr,
                       const t_mdatoms  *md);
};

} // namespace gmx

#endif
