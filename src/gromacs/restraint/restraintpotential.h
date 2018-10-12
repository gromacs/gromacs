/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2018, by the GROMACS development team, led by
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

#ifndef GROMACS_RESTRAINT_RESTRAINTPOTENTIAL_H
#define GROMACS_RESTRAINT_RESTRAINTPOTENTIAL_H

/*!
 * \defgroup module_restraint MD restraints
 * \inpublicapi
 * \brief Apply restraints during MD integration.
 *
 * The classes here are available through the public API, but only gmx::RestraintPotential
 * is necessary to implement a restraint plugin.
 */
/*! \file
 * \brief Declare generic interface for restraint implementations.
 *
 * \inpublicapi
 * \ingroup module_restraint
 */

#include <functional>
#include <memory>
#include <ostream>
#include <vector>

#include "gromacs/math/vectypes.h"
#include "gromacs/utility/basedefinitions.h"
#include "gromacs/utility/real.h"

struct gmx_mtop_t;
struct gmx_output_env_t;
struct pull_params_t;
struct pull_t;
struct t_commrec;
struct t_filenm;
struct t_inputrec;
struct t_mdatoms;
struct t_pbc;

namespace gmxapi
{
class Session;
class SessionResources;
}

namespace gmx
{

/*!
 * \brief Provide a vector type name with a more stable interface than RVec and a more stable
 * implementation than vec3<>.
 *
 * \ingroup module_restraint
 *
 * \internal
 * Type alias is used at namespace level
 */
using Vector = ::gmx::RVec;

/*!
 * \brief Structure to hold the results of IRestraintPotential::evaluate().
 *
 * \ingroup module_restraint
 */
class PotentialPointData
{
    public:


        /*!
         * \brief Initialize a new data structure.
         */
        PotentialPointData();

        /*!
         * \brief Initialize from an argument list
         *
         * \param f Force vector.
         * \param e Energy value.
         *
         * Note that if force was calculated as a scalar, it needs to be multiplied by a unit
         * vector in the direction to which it should be applied.
         * If this calculation is in a subclass of gmx::RestraintPotential,
         * you should be able to use the make_force_vec() helper function (not yet implemented).
         */
        PotentialPointData(const Vector &f, real e);

        /*!
         * \brief Force vector calculated for first position.
         */
        Vector force;

        /*!
         * \brief Potential energy calculated for this interaction.
         */
        real energy;
};

/*!
 * \brief Interface for Restraint potentials.
 *
 * Derive from this interface class to implement a restraint potential. The derived class
 * must implement the evaluate() member function that produces a PotentialPointData instance.
 * For various convenience functions and to decouple the internal library
 * interface from implementations of potentials, it is expected that
 * restraints will be programmed by subclassing gmx::RestraintPotential<>
 * rather than gmx::IRestraintPotential.
 *
 * For a set of \f$n\f$ coordinates, generate a force field according to a
 * scalar potential that is a fun. \f$F_i = - \nabla_{q_i} \Phi (q_0, q_1, ... q_n; t)\f$
 *
 * Potentials implemented with these classes may be long ranged and are appropriate
 * for only a small number of particles to avoid substantial performance impact.
 *
 * The indices in the above equation refer to the input and output sites specified
 * before the simulation is started. In the current interface, the evaluate() virtual
 * function allows an implementer to calculate the energy and force acting at the
 * first of a pair of sites. The equal and opposite force is applied at the second site.
 *
 * The potential function is evaluated with a time argument
 * and can be updated during the simulation.
 * For non-time-varying potentials, the time argument may still be useful for
 * internal optimizations, such as managing cached values.
 *
 * In the simplest and most common case, pairs of sites (atoms or groups)
 * are specified by the user and then, during integration, GROMACS provides
 * the current positions of each pair for the restraint potential to be evaluated.
 * In such a case, the potential can be implemented by overriding evaluate().
 * \todo Template headers can help to build compatible calculation methods with different input requirements.
 * For reference, see https://github.com/kassonlab/sample_restraint
 *
 * \ingroup module_restraint
 */
class IRestraintPotential
{
    public:
        virtual ~IRestraintPotential() = default;

        /*!
         * \brief Calculate a force vector according to two input positions at a given time.
         *
         * If not overridden by derived class, returns a zero vector.
         * \param r1 position of first site
         * \param r2 position of second site
         * \param t simulation time in picoseconds
         * \return force vector and potential energy to be applied by calling code.
         *
         * \todo The virtual function call should be replaced by a (table of) function objects retrieved before the run.
         */
        virtual PotentialPointData evaluate(Vector r1,
                                            Vector r2,
                                            double t) = 0;


        /*!
         * \brief Call-back hook for restraint implementations.
         *
         * An update function to be called on the simulation master rank/thread periodically
         * by the Restraint framework.
         * Receives the same input as the evaluate() method, but is only called on the master
         * rank of a simulation to allow implementation code to be thread-safe without knowing
         * anything about the domain decomposition.
         *
         * \param v position of the first site
         * \param v0 position of the second site
         * \param t simulation time
         *
         * \internal
         * We give the definition here because we don't want plugins to have to link against
         * libgromacs right now (complicated header maintenance and no API stability guarantees).
         * But once we've had plugin restraints wrap themselves in a Restraint template,
         * we can set update = 0
         *
         * \todo: Provide gmxapi facility for plugin restraints to wrap themselves
         * with a default implementation to let this class be pure virtual.
         */
        virtual void update(gmx::Vector v,
                            gmx::Vector v0,
                            double      t);


        /*!
         * \brief Find out what sites this restraint is configured to act on.
         * \return
         */
        virtual std::vector<int> sites() const = 0;

        /*!
         * \brief Allow Session-mediated interaction with other resources or workflow elements.
         *
         * \param resources temporary access to the resources provided by the session for additional configuration.
         *
         * A module implements this method to receive a handle to resources configured for this particular workflow
         * element.
         *
         * \internal
         * \todo This should be more general than the RestraintPotential interface.
         */
        virtual void bindSession(gmxapi::SessionResources* resources);
};

}      // end namespace gmx

#endif //GMX_PULLING_PULLPOTENTIAL_H
