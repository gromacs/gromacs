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
//
// Created by Eric Irrgang on 9/23/17.
//

#ifndef GROMACS_RESTRAINT_RESTRAINTPOTENTIAL_H
#define GROMACS_RESTRAINT_RESTRAINTPOTENTIAL_H

/*!
 * \defgroup module_restraint MD restraints
 * \ingroup publicapi
 * \brief Apply restraints during MD integration.
 *
 * More conceptual overview is in the \ref page_pullcodeextension "full documentation".
 *
 * The classes here are available through the public API, but only gmx::RestraintPotential
 * is necessary to implement a restraint plugin.
 */
/*! \file
 * \brief Declare generic interface for restraint implementations.
 *
 * \ingroup module_restraint
 */

#include "vectortype.h"

#include <memory>
#include <vector>
#include <functional>
#include <ostream>

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

class Mdrunner; // defined in src/programs/mdrun/runner.h

/*!
 * \brief Provide a vector type name with a more stable interface than RVec and a more stable
 * implementation than vec3<>.
 *
 * \ingroup module_restraint
 *
 * \internal
 * Type alias is used at namespace level
 */
using Vector = ::gmx::detail::vec3<real>;


/*!
 * \brief Typed unitless time in GROMACS.
 *
 * It may be helpful to explicitly specify the units of time.
 *
 * Example:
 *
 *     Time time;
 *     // float t = time; // error: no implicit conversion because units are ambiguous.
 *     using gmx::time::ps;
 *     float t = time * ps;
 *
 * TODO: If this persists, it will be moved from this header.
 * \ingroup module_restraint
 */
struct Time
{
    using float_type = real;
    float_type t;
};

/*!
 * \brief Structure to hold the results of IRestraintPotential::evaluate().
 *
 * \ingroup module_restraint
 */
class PotentialPointData
{
    public:
        /*!
         * \brief Force vector calculated for first position.
         */
        Vector force;
        /*!
         * \brief Potential energy calculated for this interaction.
         */
        real energy;


        /*!
         * \brief Initialize a new data structure.
         */
        PotentialPointData() : PotentialPointData{Vector(), real(0.0)}
        {};

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
        PotentialPointData(const Vector &f, const real e) :
            force {f},
        energy {e}
        {};
};

/*! \libinternal
 * \brief Library interface for Restraint potentials.
 *
 * Classes implementing IRestraint allow Restraints to be registered with the MD simulation.
 * To implement the interface, a class provides a `getEvaluator()` method that returns a
 * function object for a call signature `std::function<PotentialPointData(Vector, Vector, Time)>`.
 *
 * This interface is a library detail. Code providing Restraint potentials should not implement
 * IRestraintPotential directly, but should specialize the RestraintPotential template with an
 * overload of the calculate() method.
 * \ingroup module_restraint
 */
class IRestraint
{
    public:
        /// \cond
        virtual ~IRestraint() = default;
        /// \endcond

        /// \cond internal
        /*!
         * \brief Provide a function object to evaluate a pairwise restraint.
         *
         * Defines the library-level interface between the restraint management code and restraint classes.
         * Restraint potentials themselves should be implemented with the more
         * stable interface defined by gmx::Restraint.
         *
         * \return Function object to call when evaluating each restraint.
         *
         */
        virtual std::function<PotentialPointData(const Vector &, const Vector &, Time)> getEvaluator() = 0;
        /// \endcond
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
 * The indices in the above equation refer to the input and output sites specified before the simulation is started.
 * In the current interface, the evaluate() virtual function allows an
 * implementer to calculate the energy and force acting at the first of a pair of sites. The equal and opposite force is applied at the second site.
 *
 * The potential function is evaluated with a time argument
 * and can be updated during the simulation.
 * For non-time-varying potentials, the time argument may still be useful for
 * internal optimizations, such as managing cached values.
 *
 * In the simplest and most common case, pairs of sites (atoms or groups)
 * are specified by the user and then, during integration, GROMACS provides
 * the current positions of each pair for the restraint potential to be evaluated.
 * In such a case, the potential can be implemented by overriding...
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
         * But once we've had plugin restraints wrap themselves in a Restraint template, we can set update = 0
         * \todo: Provide gmxapi facility for plugin restraints to wrap themselves with a default implementation to let this class be pure virtual.
         */
        virtual void update(gmx::Vector v,
                            gmx::Vector v0,
                            double      t) { (void)v; (void)v0; (void)t; };


        /*!
         * \brief Find out what sites this restraint is configured to act on.
         * \return
         */
        virtual
        std::vector<unsigned long int> sites() const = 0;

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
        virtual void bindSession(gmxapi::SessionResources* resources)
        {
            (void) resources;
        }
};

/*!
 * \brief Allow stable implementation of restraint potentials.
 *
 * \tparam T implementation type
 *
 * Implementation of a potential is encapsulated in class T. This template
 * then serves to "mix in from below" the functionality to serve as a Restraint
 * in the MD code.
 *
 * To implement a restraint,
 *
 * 1. Write a class with a member function called `calculate` that returns a PotentialPointData object and that uses one of the supported parameter lists.
 * 2. Instantiate the RestraintPotential template specialization for the class
 *    by registering or exporting it for access by a user.
 *
 * The most complete call signature available for a `calculate` function is
 *
 *     PotentialPointData calculate(const Vector& position1, const Vector& position2, Time t);
 *
 * Internally, template matching is used to derive type traits with which to
 * tag the restraint class so that an appropriate function call can be dispatched.
 *
 * Note that the floating point type definition `real` from real.h should be used if your code
 * can flexibly use both double precision and single precision (depending on how GROMACS was compiled).
 *
 * Example:
 * \code
 * #include "gromacs/restraint/restraintpotential.h"
 * #include "gromacs/restraint/vectortype.h"
 *
 * class MyRestraint
 * {
 * public:
 *      gmx::PotentialPointData calculate(real[3] position1, real[3] position2, real t);
 * };
 *
 * gmx::PotentialPointData::calculate(const Vector& r1, const Vector& r2)
 * {
 *      // Apply constant harmonic force between the two sites
 *      real k = 1.0;
 *      real R0 = 1.0;
 *
 *      // Use some vector helper functions
 *      auto diff = norm(r2 - r1) - R0;
 *      auto unitvec = (r2 - r1)/norm(r2 - r1);
 *
 *      auto force = -k * diff * unitvec;
 *      auto energy = 0.5 * k * diff**2;
 *
 *      return {force, energy};
 * }
 *
 * \endcode
 * Todo: explain registration
 *
 * \ingroup module_restraint
 */

//template<class T>
//class RestraintForceProvider : public gmx::IForceProvider
//{
//    public:
//        void calculateForces(const t_commrec *cr,
//                             const t_mdatoms *mdatoms,
//                             const matrix box,
//                             double t,
//                             const rvec *x,
//                             gmx::ArrayRef<gmx::RVec> force)
//        override
//        {
////                T::calculate();
//                force_called++;
//        };
//
//        unsigned int force_called{0};
//};

//template<class T>
//class RestraintMDModule final : public gmx::IMDModule, public T
//{
//    private:
//        class OptionProvider : public gmx::IMdpOptionProvider
//        {};
//        std::shared_ptr<OptionProvider> optionprovider{std::make_shared<OptionProvider>()};
//
//        class OutputProvider : public gmx::IMDOutputProvider
//        {};
//        std::shared_ptr<OutputProvider> outputprovider{std::make_shared<OutputProvider>()};
//
//        std::shared_ptr<RestraintForceProvider<T>> forceprovider{std::make_shared<RestraintForceProvider<T>>()};
//
//        gmx::IForceProvider* getForceProvider()
//        {
//                return forceprovider.get();
//        };
//    public:
//        gmx::IMdpOptionProvider *mdpOptionProvider() override
//        {
//                return optionprovider.get();
//        }
//
//        gmx::IMDOutputProvider *outputProvider() override
//        {
//                return outputprovider.get();
//        }
//
//        void initForceProviders(ForceProviders *forceProviders) override
//        {
//                forceProviders->addForceProvider(getForceProvider());
//        }
//
//        unsigned int force_called() { return forceprovider->force_called; };
//};
template<class T>
class RestraintPotential : public IRestraint, public IRestraintPotential, public T
{
    public:
        /// \cond libapi
        /*!
         * \brief Implement gmx::IRestraint interface
         *
         * \return function object with the appropriate call signature for the library.
         */
        std::function<PotentialPointData(const Vector &,
                                         const Vector &,
                                         Time)> getEvaluator() override;
        /// \endcond

        static std::shared_ptr<RestraintPotential> create();
    private:
        /// Private constructor (use create() method).
        RestraintPotential() = default;

        /// Allow interfaces to keep this alive.
        std::weak_ptr < RestraintPotential < T>> self_;
};
//
//std::shared_ptr<RestraintPotential> RestraintPotential::create()
//{
//        auto newObject = std::make_shared<RestraintPotential>();
//        newObject->self_ = newObject;
//        return newObject;
//}
//
//template<class T>
//std::shared_ptr<::gmx::IMDModule> RestraintPotential<T>::getModuleInterface()
//{
//        class Container : public IMDModule
//        {
//            public:
//                IMdpOptionProvider *mdpOptionProvider() override
//                {
//                        return nullptr;
//                }
//
//                IMDOutputProvider *outputProvider() override
//                {
//                        return nullptr;
//                }
//
//                void initForceProviders(ForceProviders *forceProviders) override
//                {
//
//                }
//        };
//        std::shared_ptr<RestraintPotential<T>> self = self_.lock();
//        auto interface = std::make_shared<Container>();
//        return interface;
//}

/*! \fn PotentialPointData RestraintPotential<T>::calculate(Vector r1, Vector r2, Time t);
 * \brief The most complete call signature available for a RestraintPotential.
 *
 * \param r1 First position to consider for the restraint pair.
 * \param r2 Second position to consider for the restraint pair.
 * \param t Simulation time in picoseconds.
 * \return Energy and force vector, relative to the first member of the pair.
 *
 * \ingroup module_restraint
 */
//        PotentialPointData calculate(Vector r1,
//                                     Vector r2,
//                                     Time t);

/*! \fn PotentialPointData calculate(Vector r1, Vector r2)
 * \brief
 *
 * \param r1 First position to consider for the restraint pair.
 * \param r2 Second position to consider for the restraint pair.
 * \return Energy and force vector, relative to the first member of the pair.
 *
 * \ingroup module_restraint
 */
//        PotentialPointData calculate(Vector r1,
//                                     Vector r2);

/*! \fn PotentialPointData calculate(real distance, Time t);
 * \brief
 *
 * \param distance
 * \param t Simulation time in picoseconds.
 * \return Energy and force vector, relative to the first member of the pair.
 *
 * \ingroup module_restraint
 */
//        PotentialPointData calculate(real distance,
//                                     Time t);

/*! \fn PotentialPointData calculate(real distance)
 * \brief
 *
 * \param distance
 * \return Energy and force vector, relative to the first member of the pair.
 *
 * \ingroup module_restraint
 */
//        PotentialPointData calculate(real distance);

/*!
 * \brief Encapsulate the old pulling schemes.
 *
 * Give the old and new pulling code the same management and wrap up the old set of
 * schemes into a single pulling class. If the input record indicates a pulling protocol,
 * the RestraintPotential may use the associated resources (pull_t, pull_params_t, ...).
 *
 * \ingroup module_restraint
 */
class LegacyPuller
{
    public:
        PotentialPointData calculate(Vector r1,
                                     Vector r2,
                                     Time   t);

        ~LegacyPuller() = default;

        /*!
         * \brief Copy the metadata, but manage the same internals.
         *
         * \param source object to copy
         */
        LegacyPuller(const LegacyPuller &source);

        /*!
         * \brief Copy the metadata, but manage the same internals.
         *
         * \param source object to copy
         */
        LegacyPuller &operator=(const LegacyPuller &source);

        /*!
         * \brief Move ownership of pull code manager.
         *
         * \param old object to move from.
         */
        LegacyPuller(LegacyPuller &&old) noexcept;

        /*!
         * \brief Move ownership of pull code manager.
         *
         * \param old move ownership from rhs of operation.
         */
        LegacyPuller &operator=(LegacyPuller &&old) noexcept;

        /*!
         * \brief Construct manager for legacy pulling code.
         *
         * \param pullWorkPointer pointer to struct created by init_pull() and managed by caller.
         */
        explicit LegacyPuller(pull_t* pullWorkPointer);

        struct pull_t      * getRaw();
        const struct pull_t* getRaw() const;

    private:
        // borrow access to pull_t owned by calling code.
        struct pull_t* pullWorkPointer_;
};

} // end namespace gmx

/*!
 * \brief Hold an arbitrary number of gmx::IRestraintPotential objects.
 *
 * \ingroup module_restraint
 */
class PotentialContainer
{
    public:
        /// Default constructor
        PotentialContainer();
        /// Copy constructor disabled.
        PotentialContainer(const PotentialContainer &) = delete;
        /// Move construction allowed.
        PotentialContainer(PotentialContainer &&) noexcept;
        /// Copy assignment disabled.
        PotentialContainer &operator=(const PotentialContainer &) = delete;
        /// Allow move assignment.
        PotentialContainer &operator=(PotentialContainer &&) noexcept;
        /// non-virtual destructor.
        ~PotentialContainer();

        /*!
         * \brief Add a pulling potential to the managed list.
         *
         * \param puller GROMACS-provided or custom pulling potential
         */
        void addPotential(std::shared_ptr<gmx::IRestraintPotential> puller) noexcept;

        using RestraintIterator = std::vector < std::shared_ptr < gmx::IRestraintPotential>>::iterator;
        RestraintIterator begin();
        RestraintIterator end();

    private:
        /// Private implementation class.
        class Impl;
        /// Implementation object.
        std::unique_ptr<Impl> impl_;
};

#endif //GMX_PULLING_PULLPOTENTIAL_H
