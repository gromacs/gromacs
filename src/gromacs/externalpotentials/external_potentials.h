/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2016, by the GROMACS development team, led by
 * Mark Abraham, David van der Spoel, Berk Hess, and Erik Lindahl,
 * and including many others, as listed in the AUTHORS file in the
 * top-level source directory and at http://www.gromacs.org.
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
/*! \brief
 * Declares abstract interface for external potentials in mdrun.
 *
 * In order to maximize encapsulation of functionality and to minimize
 * potential impact on mdrun of errors in external and contributed
 * code, an interface is provided for this kind of code.
 * The idea is that mdrun first queries the code to see what kind of
 * optimizations it can do, and then allocates resources to it that
 * can be optimized on the fly. For this to work the external potential
 * should be able to reorganize the parallellization on the fly as well.
 *
 * \inpublicapi
 */
#ifndef GMX_EXTERNALPOTENTIALS_EXTERNAL_POTENTIALS_H
#define GMX_EXTERNALPOTENTIALS_EXTERNAL_POTENTIALS_H

class ExternalPotential
{
    //! Constructor
    ExternalPotential() {};

    //! Destructor
    virtual ~ExternalPotential();

    //! Return whether the code supports OpenMP parallellization
    virtual bool supportOpenMP();

    //! Return whether the code supports MPI parallellization
    virtual bool supportMPI();

    //! Return whether the code supports GPU optimization
    virtual bool supportGPU();

    /*! \brief Initialization routine
     *
     * \param[in] input String giving information to the external potential
     *                  for instance an input file name.
     * \param[in] mtop  Molecular topology object
     */
    virtual void initialize(const std::string input,
                            const gmx_mtop_t  mtop);

    /*! \brief Set up internal parallellization
     *
     * The code should be able to modify its internal parallellization
     * in response to changes in CPU/core allocation due to load balancing
     *
     * \param[in] cr Communication data structure
     * \return 0 if successful, an error code otherwise.
     */
    virtual int setupParallellization(const t_commrec *cr);

    /*! \brief Do an energy/force evaluation
     *
     * This routine does the actual work. It can not touch the coordinates
     * but should not assume anything about the forces, just add to them.
     * Total energy returned in a scalar variable. Of course stuff can be added
     * here, such as dH/dlambda, or electric fields input or output.
     *
     * \param[in]  x      The atomic coordinates.
     * \param[out] f      The atomic forces.
     * \param[out] energy The scalar energy due to this interaction.
     * \return 0 if successful, or an error code otherwise.
     */
    virtual int computeForceEnergy(const rvec *x,
                                   rvec        f[],
                                   real       *energy);

    /*! \brief Return descriptive error message
     *
     * \param[in] errorCode Code corresponding to a known error state
     *                      for the module
     * \return An error string or NULL if the errorCode is out of range
     */
    virtual const char *errorMessage(const int errorCode);
};

#endif
