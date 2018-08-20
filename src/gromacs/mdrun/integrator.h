/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2015,2016,2017,2018, by the GROMACS development team, led by
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
/*! \file
 * \brief Declares the integrator interface for mdrun
 *
 * \author David van der Spoel <david.vanderspoel@icm.uu.se>
 * \author Mark Abraham <mark.j.abraham@gmail.com>
 * \ingroup module_mdrun
 */
#ifndef GMX_MDRUN_INTEGRATOR_H
#define GMX_MDRUN_INTEGRATOR_H

#include <cstdio>

#include <memory>

#include "gromacs/compat/make_unique.h"
#include "gromacs/mdrun/context.h"
#include "gromacs/utility/basedefinitions.h"
#include "gromacs/utility/real.h"

class energyhistory_t;
struct gmx_enfrot;
struct gmx_mtop_t;
struct gmx_membed_t;
struct gmx_multisim_t;
struct gmx_output_env_t;
struct gmx_vsite_t;
struct gmx_wallcycle;
struct gmx_walltime_accounting;
struct MdrunOptions;
struct ObservablesHistory;
struct ReplicaExchangeParameters;
struct t_commrec;
struct t_fcdata;
struct t_forcerec;
struct t_filenm;
struct t_inputrec;
struct t_nrnb;
class t_state;

namespace gmx
{

class BoxDeformation;
class Constraints;
class IMDOutputProvider;
class MDLogger;
class MDAtoms;
class IntegratorDispatcherAdapter;

//! Function type for integrator code.
using IntegratorFunctionType = void();

/*!
 * \brief Placeholder named type for integrator methods.
 *
 * Essentially an alias for `unsigned int`, used for enumerated integrator type,
 * until a named enum type is established universally. Has no default constructor,
 * which helps us catch errors with uninitialized members of other objects.
 *
 * \ingroup module_mdrun
 */
class SimulationMethod
{
    public:
        /*! \brief Actual value of enumerated integrator type. */
        unsigned int method_;
        /*!
         * \brief Converting constructor from enumerated type.
         *
         * \param t selected simulation method
         *
         * \see md_enums.h
         */
        explicit SimulationMethod(unsigned int t) : method_ {t}
        {};
        /*!
         * \brief Conversion operator
         *
         * \return value with underlying enumeration type.
         */
        operator unsigned int () { return method_; }
};

/*!
 * \brief Interface for MD integration and other MM methods.
 *
 * A std::unique_ptr<Integrator> serves as a handle to the object produced by
 * IntegratorBuilder::build().
 *
 * \ingroup module_mdrun
 */
class IIntegrator
{
    public:
        /*!
         * \brief Allow client to own an Integrator object.
         */
        virtual ~IIntegrator();

        /*!
         * \brief Call functor.
         *
         * Runs the configured simulator providing this Integrator interface.
         */
        virtual void run() = 0;
};

/*!
 * \brief Mix-in for homeless MD parameters.
 *
 * This is a catch-all parameter container that can be publicly inherited by
 * Integrator classes that need convenient access to the old style parameter pack.
 * It provides "mix-in" data for Integrator classes that
 * have not been migrated to newer data structures. The const reference fields
 * require initialization on construction, but the POD structure allows
 * aggregate initialization. Thus the base class can only be usefully initialized
 * with a template constructor in the derived class that forwards arguments to the
 * Container mix-in.
 *
 * Example:
 *
 *      struct MyIntegrator : public IntegratorParamsContainer
 *      {
 *      template<class... Args> MyIntegrator(Args&&...args) :
 *          IntegratorParamsContainer{std::forward<Args>(args)...} {};
 *      }
 *
 * The IntegratorParamsContainer class is provided to compartmentalize the
 * initialization of the old style parameter pack. As parameters are moved to
 * more modern structures, this class can have members removed without causing
 * source incompatibilities in derived classes.
 *
 * Aggregate initialization is used, for which the chief risk is that
 * if a member is added at the end and not all initializer lists are
 * updated, then the member will be value initialized, which will
 * typically mean initialization to zero.
 */
struct IntegratorParamsContainer
{
    //! Handles logging.
    FILE                            *fplog;
    //! Handles communication.
    t_commrec                       *cr;
    //! Coordinates multi-simulations.
    const gmx_multisim_t            *ms;
    //! Handles logging.
    const MDLogger                  &mdlog;
    //! Count of input file options.
    int                              nfile;
    //! Content of input file options.
    const t_filenm                  *fnm;
    //! Handles writing text output.
    const gmx_output_env_t          *oenv;
    //! Contains command-line options to mdrun.
    const MdrunOptions              &mdrunOptions;
    //! Handles virtual sites.
    gmx_vsite_t                     *vsite;
    //! Handles constraints.
    Constraints                     *constr;
    //! Handles enforced rotation.
    gmx_enfrot                      *enforcedRotation;
    //! Handles box deformation.
    BoxDeformation                  *deform;
    //! Handles writing output files.
    IMDOutputProvider               *outputProvider;
    //! Contains user input mdp options.
    t_inputrec                      *inputrec;
    //! Full system topology.
    gmx_mtop_t                      *top_global;
    //! Helper struct for force calculations.
    t_fcdata                        *fcd;
    //! Full simulation state (only non-nullptr on master rank).
    t_state                         *state_global;
    //! History of simulation observables.
    ObservablesHistory              *observablesHistory;
    //! Atom parameters for this domain.
    MDAtoms                         *mdAtoms;
    //! Manages flop accounting.
    t_nrnb                          *nrnb;
    //! Manages wall cycle accounting.
    gmx_wallcycle                   *wcycle;
    //! Parameters for force calculations.
    t_forcerec                      *fr;
    //! Parameters for replica exchange algorihtms.
    const ReplicaExchangeParameters &replExParams;
    //! Parameters for membrane embedding.
    gmx_membed_t                    *membed;
    //! Manages wall time accounting.
    gmx_walltime_accounting         *walltime_accounting;
    //! We only intend to construct such objects with an initializer list.
#if __GNUC__ > 4 || (__GNUC__ == 4 && __GNUC_MINOR__ >= 9)
    // Aspects of the C++11 spec changed after GCC 4.8.5, and
    // compilation of the initializer list construction in runner.cpp
    // fails in GCC 4.8.5.
    IntegratorParamsContainer() = delete;
#endif
};

/*!
 * \brief Adapter class to help construction of IntegratorParamsContainer.
 *
 * IntegratorParamsContainer cannot be constructed with expressions that use
 * perfect forwarding, such as `make_unique` and `emplace` because it requires
 * aggregate initialization. This is a simple wrapper for the aggregate with a
 * template constructor that can be used as a convenient way to create a
 * IntegratorParamsContainer, or as an example to other classes using
 * IntegratorParamsContainer as a mix-in.
 */
class IntegratorAggregateAdapter : public IntegratorParamsContainer
{
    public:
/*!
 * \brief Constructor templated for forwarding arguments.
 *
 * \tparam Args sequence of types in IntegratorParamsContainer.
 * \param args sequence of values for aggregate initialization of IntegratorParamsContainer.
 *
 * Arguments to this constructor are forwarded to the underlying IntegratorParamsContainer
 * to provide an aggregate initialization adapter. This allows IngtegratorAggregateAdapter
 * to be used in non-forwarding contexts (such as emplace_back or make_unique) to access the
 * aggregate initializer for IntegratorParamsContainer without specifying the full list of
 * parameters anywhere other than the definition of IntegratorParamsContainer.
 */
        template<class ... Args> IntegratorAggregateAdapter(Args && ... args) :
            IntegratorParamsContainer {std::forward<Args>(args) ...}
        {};
};

/*!
 * \brief Wrapper to implement IIntegrator in terms of IntegratorDispatcher.
 *
 * Allows functor objects to provide a run() callable for preconfigured objects
 * implementing different methods. IntegratorDispatcher can be deconstructed
 * over the course of future development work to provide more concise IIntegrator
 * implementations.
 *
 * Also serves as an adapter to the aggregate initialization of IntegratorDispatcher.
 * Since we cannot forward an aggregate initialization list, we provide a
 * template constructor that directly calls the aggregate initialization constructor
 * of IntegratorDispatcher.
 *
 * Having multiple integrators as member functions isn't a good
 * design, and we definitely only intend one to be called, but the
 * goal is to make it easy to change the names and types of members
 * without having to make identical changes in several places in the
 * code. Once many of them have become modules, we should change this
 * approach.
 *
 * Note that the presence of const reference members means that the
 * default constructor would be implicitly deleted. But we only want
 * to make one of these when we know how to initialize these members,
 * so that is perfect. To ensure this remains true even if we would
 * remove those members, we explicitly delete this constructor.
 * Other constructors, copies and moves are OK.
 *
 * In GROMACS 2019, IntegratorDispatcher is both a dispatcher and a container
 * for integrator parameters.
 *
 * \ingroup module_mdrun
 */
struct IntegratorDispatcher : public IIntegrator, public IntegratorParamsContainer
{
    /*!
     * \brief Construct the simulation method dispatcher.
     *
     * This constructor is intended as an internal library interface only. Client
     * code should obtain a IntegratorDispatcher object as a IIntegrator handle
     * from IntegratorBuilder.
     *
     * Depending on the value provided for method, dispatches run() to do MM
     * calculation or MD integration for steepest descent, conjugate
     * gradient energy minimization, normal mode analysis, or test particle insertion.
     *
     * \param method a supported integration method \see md_enums.h
     * \param container old style parameter pack.
     */
    IntegratorDispatcher(SimulationMethod method, const IntegratorParamsContainer &container) :
        IntegratorParamsContainer {container},
    method_ {method}
    {};

    //! Implements steepest descent EM.
    IntegratorFunctionType           do_steep;
    //! Implements conjugate gradient energy minimization
    IntegratorFunctionType           do_cg;
    //! Implements onjugate gradient energy minimization using the L-BFGS algorithm
    IntegratorFunctionType           do_lbfgs;
    //! Implements normal mode analysis
    IntegratorFunctionType           do_nm;
    //! Implements test particle insertion
    IntegratorFunctionType           do_tpi;
    /*!
     * \brief Implement functor callable.
     */
    void run() override;

    /*!
     * \brief Integration method to be dispatched.
     *
     * \return The selected (enumerated) integration method. (\see md_enums.h)
     */
    SimulationMethod getMethod() const;
    private:
        //! simulation method for which we are dispatching.
        SimulationMethod method_;
};


/*!
 * \brief Generic type for client code to configure a new integrator.
 *
 * Client code (such as gmx::Mdrunner) creates a builder for a particular
 * simulation method by calling `auto builder = IntegratorBuilder::create(method)`
 * where `method` is a selection from the enumerated integrator types. Required
 * components for the integrator can be added with the builder member functions.
 * Once configured, a handle to a new integrator is produced with `builder.build()`.
 *
 * While GROMACS evolves, a catch-all member function `setParams` must be called
 * for all builders (unless documented otherwise) before `build` is called.
 */
class IntegratorBuilder final
{
    public:
        /*!
         * \brief Implementation Base class
         *
         * Builders for specific integration methods can include the implementation
         * header and derive from IntegratorBuilder::Base.
         */
        class Base;

        /*!
         * \brief Sentry to help coordinate the validity of externally-owned data.
         *
         * Returned by setParams(). Initial implementation is a trivial and useless
         * object, but provides forward compatibility for an API with which a Builder
         * and client could interact regarding the validity of data provided to setParams().
         */
        class DataSentry final
        {
            public:
                /*!
                 * \brief Expire the monitored data.
                 *
                 * When the DataSentry is destroyed, such as by going out of scope,
                 * the destructor performs notifications or other actions specified by
                 * the API (not yet implemented) to allow interested Integrator instances
                 * to react to the invalidity of pointers and references provided during
                 * setParams().
                 */
                ~DataSentry();

                /*!
                 * \brief Construction protocol not yet specified.
                 *
                 * Default constructor provided for initial implementation. Future API
                 * development will likely require changes to construction methods.
                 */
                DataSentry();

                /*!
                 * \brief DataSentry objects are not copyable.
                 *
                 * DataSentry instances should have one-to-one correspondence to calls
                 * to setParams().
                 *
                 * \{
                 */
                DataSentry(const DataSentry&) = delete;
                DataSentry                 &operator=(const DataSentry &) = delete;
                /// \}

                /*!
                 * \brief DataSentry objects are move-able, but semantics are not specified.
                 *
                 * Client code should transfer ownership of DataSentry objects if necessary
                 * to provide accurate state regarding the validity of data. This feature is
                 * provided for future compatibility, and the initial
                 * implementation does not specify semantics beyond the fact that objects can be moved.
                 *
                 * \{
                 */
                DataSentry(DataSentry &&) noexcept;
                DataSentry &operator=(DataSentry &&) noexcept;
                /*! \} */
        };

        /*!
         * \brief No public default constructor.
         *
         * Clients use createIntegratorBuilder() to get an object of this type.
         */
        IntegratorBuilder() = delete;
        ~IntegratorBuilder();
        /*! \cond */
        IntegratorBuilder(const IntegratorBuilder&)             = delete;
        IntegratorBuilder &operator=(const IntegratorBuilder &) = delete;
        IntegratorBuilder(IntegratorBuilder &&) noexcept;
        IntegratorBuilder &operator=(IntegratorBuilder &&) noexcept;
        /*! \endcond */

        /*!
         * \brief Set GROMACS 2019 style parameter pack required by all integrators.
         *
         * For simulation methods requiring the GROMACS 2019 style integrator parameter pack,
         * setParams be called exactly once before calling `build()`. This method is atypical
         * of a Builder pattern in a number of ways. It does not return a reference to the
         * Builder. It does not cause the Builder to create new data or transfer ownership
         * of anything. setParams() provides a parameter pack that includes non-owning
         * pointers and/or references to data owned by the caller. To help the integrator and
         * caller manage assumptions about data lifetimes, setParams() returns a sentry object.
         *
         * If not assigned to a variable in the calling code, the sentry object persists for
         * the lifetime of the enclosing block, which is assumed to correspond to guarantees
         * on the lifetime of the arguments provided to setParams(). If the caller can or
         * should manage data lifetimes more directly, the client can manage the lifetime of
         * the sentry. Builders that require more interaction with data life times can use the
         * sentry API (not yet implemented) before the sentry object is returned.
         *
         * \tparam ArgsT sequence of types corresponding to the members of IntegratorDispatcher.
         * \param args sequence of values with which to perform an aggregate initialization of a new
         * IntegratorDispatcher.
         * \return a sentry object whose life time should not exceed the life time of any of the args.
         *
         * \note In this implementation, if a client provides a set of arguments that does not correspond
         * to a valid IntegratorDispatcher constructor, the template will be instantiated, but will
         * either not compile or may fail to link, depending on the current implementation. More robust
         * and helpful failure scenarios are pending additional design decisions.
         */
        template<class ... ArgsT>
        DataSentry setParams(ArgsT && ... args);

        /*!
         * \brief Provide a simulation runtime context for the method's implementation.
         *
         * \param context handle to an execution context manager owned by the client code.
         * \return reference to current builder
         */
        IntegratorBuilder &addContext(const md::Context &context);

        /*!
         * \brief Get ownership of a new Integrator object.
         *
         * \return handle to a new object implementing the configured integrator
         */
        std::unique_ptr<IIntegrator> build();

        /*!
         * \brief Create a new builder.
         *
         * Factory function to create a new IntegratorBuilder object. The
         * implementation of the builder returned is dependent on the type of
         * simulation method specified. Required or available builder methods
         * will be affected. Refer to documentation for the various simulation
         * methods or the classes derived from IntegratorBuilder::Base.
         *
         * \param integratorType selection of simulation method
         * \return new builder object
         */
        static IntegratorBuilder create(const SimulationMethod &integratorType);


    private:
        explicit IntegratorBuilder(std::unique_ptr<IntegratorBuilder::Base> impl);

        std::unique_ptr<Base> impl_;
};


/*!
 * \brief Generic interface for IIntegrator builders.
 *
 * Base class from which method-specific builders can be derived. Provides
 * default implementations of methods required to implement IntegratorBuilder.
 * Derived classes must implement the build() method.
 */
class IntegratorBuilder::Base
{
    public:
        virtual ~Base();

        /*!
         * \brief Catch-all setter.
         *
         * This function does not conform to the typical patter of returning
         * a reference to *this (the builder) because it would result in slicing
         * of the derived classes.
         *
         * \todo Return an RAII sentry object to check for parameter validity or changing scope.
         */
        virtual DataSentry setAggregateAdapter(std::unique_ptr<IntegratorAggregateAdapter> container);

        /*!
         * \brief Provide a context manager to the simulation method.
         *
         * Some simulation methods require external resources to be provided
         * related to the execution environment and not already managed by some
         * other module. For such methods (e.g. MD) a md::Context should be
         * provided exactly once before build() is called. Methods that do not
         * use a Context may accept, but ignore calls to addContext().
         *
         * \param context handle to Context owned by the client.
         * \return reference to the current builder implementation.
         */
        virtual Base &addContext(const md::Context &context);

        /*!
         * \brief Build the Integrator product.
         *
         * \return Ownership of a new object providing the chosen simulation method.
         *
         * Implementation depends on the simulation method specified to IntegratorBuilder::create()
         *
         * Builders for different methods will require different builder methods
         * to be called before build() can successfully produce an initialize object.
         * After the call to build(), the builder is in an undefined state and should
         * not be used again.
         */
        virtual std::unique_ptr<IIntegrator> build() = 0;

};

template<class ... ArgsT>
IntegratorBuilder::DataSentry IntegratorBuilder::setParams(ArgsT && ... args)
{
    auto container = gmx::compat::make_unique<IntegratorAggregateAdapter>(std::forward<ArgsT>(args) ...);
    return impl_->setAggregateAdapter(move(container));
}

}      // namespace gmx

#endif // GMX_MDRUN_INTEGRATOR_H
