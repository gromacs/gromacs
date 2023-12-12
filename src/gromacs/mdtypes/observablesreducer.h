/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2021- The GROMACS Authors
 * and the project initiators Erik Lindahl, Berk Hess and David van der Spoel.
 * Consult the AUTHORS/COPYING files and https://www.gromacs.org for details.
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
 * https://www.gnu.org/licenses, or write to the Free Software Foundation,
 * Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA.
 *
 * If you want to redistribute modifications to GROMACS, please
 * consider that scientific software is very special. Version
 * control is crucial - bugs must be traceable. We will be happy to
 * consider code for inclusion in the official distribution, but
 * derived work must not be called official GROMACS. Details are found
 * in the README & COPYING files - if they are missing, get the
 * official version at https://www.gromacs.org.
 *
 * To help us fund GROMACS development, we humbly ask that you cite
 * the research papers on the package. Check out https://www.gromacs.org.
 */
/*! \libinternal \file
 * \brief
 * Declares gmx::ObservablesReducer and builder
 *
 * Periodically modules implementing MD simulations need to
 * communicate with all collaborating ranks to do things like compute
 * observables like total energies, signal conditions, and check
 * internal consistency. This communication synchronizes all
 * participating ranks, which limits scaling and performance, so it is
 * done rarely (typically once per step) and only when required.
 *
 * Modules may provide data of type double to be reduced across
 * all ranks via an MPI all-reduce with MPI_SUM. Double-precision
 * floating-point is chosen so that no meaningful precision is lost
 * e.g. in computing energies, while also permitting integral or
 * boolean messages to be passed as double-precision floating-point
 * values.
 *
 * Different modules typically need to communicate on different MD
 * steps, so in principle one might optimize by filling a
 * std::vector<double> with the values required on the current
 * step. However, that requires that each module produce and then copy
 * to the reduction buffer the data for this step. The typical amount
 * of data required even if all modules need to participate (ie.
 * hundreds of doubles) is smaller than the message headers that are
 * used by the underlying network transport protocol. So optimizing
 * for minimum message size is not particularly effective because it
 * does not meaningfully reduce the total time taken to communicate.
 *
 * Instead, we always reduce a buffer of the size that would be needed
 * if all active modules required communication this step. Then no
 * module needs to copy data merely to achieve reduction. To achieve
 * this, each module needs a stable view of memory into which it can
 * store data for which reduction is desired. It also means that
 * modules not active in the current simulation do not contribute to
 * the workload at run time. Also, modules that are active but don't
 * need communication at any particular MD step can passively opt out
 * and that incurs no overhead.
 *
 * The functionality is separated two main components, one that does
 * work during the simulation, and a builder that is used only during
 * setup time. This separates the responsibilities of
 * - allowing subscription and building the communication buffer, from
 * - orchestrating the minimum activity needed for this MD step.
 *
 * The interaction diagrams for those two workflows are depicted
 * below.
 *
\msc
wordwraparcs=true,
hscale="2";

runner [label="runner"],
builder [label="builder"],
moduleA [label="moduleA"],
moduleB [label="moduleB"],
observablesReducer [label="observablesReducer"];

runner =>> builder [label="makes"];

runner =>> moduleA [label="makes"];
runner =>> moduleA [label="passes builder to"];
moduleA =>> builder [label="subscribes itself to"];

runner =>> moduleB [label="makes"];
runner =>> moduleB [label="passes builder to"];
moduleB =>> builder [label="subscribes itself to"];

runner =>> builder [label="calls build()"];
builder =>> builder [label="makes communication\nbuffer"];
builder =>> moduleA [label="notifies of\ncallback and view"];
builder =>> moduleB [label="notifies of\ncallback and view"];
builder =>> observablesReducer [label="makes"];

\endmsc

Once the \c observablesReducer is built, the builder may be
destructed.

The \c observablesReducer and its modules operate entirely by
passing callbacks.

\msc
wordwraparcs=true,
hscale="2";

runner [label="runner"],
moduleA [label="moduleA"],
moduleB [label="moduleB"],
observablesReducer [label="observablesReducer"],
compute_globals [label="compute_globals()"];

runner =>> moduleA [label="asks for work"];
moduleA =>> moduleA [label="Produces values\nto reduce"];
moduleA =>> observablesReducer [label="requires reduction from"];

runner =>> moduleB [label="asks for work"];
moduleB =>> moduleB [label="Produces values\nto reduce"];
moduleB =>> observablesReducer [label="requires reduction from"];

runner =>> runner [label="Does other things also"];
runner =>> compute_globals [label="asks to do reduction"];
compute_globals =>> compute_globals [label="prepares data to\nreduce in\nlegacy style"];
compute_globals =>> observablesReducer [label="asks for\nbuffer view"];
observablesReducer =>> compute_globals [label="provides\nbuffer view"];
compute_globals =>> compute_globals [label="Does MPI_Allreduce"];
compute_globals =>> observablesReducer [label="notifies after\nreduction"];
observablesReducer =>> moduleA [label="notifies after reduction"];
moduleA =>> moduleA [label="Uses reduced values"];
moduleA =>> observablesReducer [label="returns"];
observablesReducer =>> moduleB [label="notifies after reduction"];
moduleB =>> moduleB [label="Uses reduced values"];
moduleB =>> observablesReducer [label="returns"];
observablesReducer =>> observablesReducer [label="zeroes reduction buffer"];
observablesReducer =>> compute_globals [label="returns"];

runner =>> observablesReducer [label="notifies at end of step"];


\endmsc
 *
 * Three callbacks are produced and called per participating module:
 *
 * 1. One produced by the module and passed to the builder so that
 *    later the ObservablesReducer can call it to notify the module
 *    that reduction is complete.
 * 2. One produced by the builder and returned to the module so the
 *    latter can call it to require reduction when it wishes
 * 3. One produced by the module and passed to the builder so the
 *    latter can call it to notify the former of the buffer view
 *    it should use in the first callback and receive a copy
 *    of the second callback.
 *
 * Modules often request that reduction occur "soon" ie. this step or
 * next step, depending whether reduction has already take place this
 * MD step. However they are also able to request reduction to occur
 * "eventually" ie. only whenever some other module requires it, so
 * the total number of reductions is minimized. Naturally, the
 * callback to such a module happens only after the eventual
 * reduction, which may happen on the same step or a later one. If a
 * module makes more than one "eventually" reduction request before
 * reduction takes place, the callback to that module will be called
 * multiple times when eventually reduction does take place. It is the
 * responsibility of the module to refrain from making those requests
 * if the multiple callbacks would be a problem (e.g. maintain an
 * internal record of whether a reduction request has been made).
 * Modules are not required to set any value for reduction unless they
 * are requesting reduction.
 *
 * An ObservablesReducer object is intended to replace the use of \c
 * compute_globals() by simulations, as
 * https://gitlab.com/gromacs/gromacs/-/issues/3887 progresses. When
 * no modules using the legacy style communication remain, it is
 * anticipated that this class will change to contain an MPI
 * communicator to use to implement the MPI_Allreduce internally.  At
 * that time, communicationBuffer() and reductionComplete() will
 * likely change into a doReduction() method, or similar. The flow of
 * the whole propagator loop will now be less clear inasmuch as the
 * responsibility for requesting reduction now lies with each module,
 * however this is probably still more clear than the large forest of
 * flags that resulted from all modules having to have their control
 * logic in the propagator loop.
 *
 * \inlibraryapi
 * \ingroup module_mdtypes
 */
#ifndef GMX_MDTYPES_OBSERVABLESREDUCER_H
#define GMX_MDTYPES_OBSERVABLESREDUCER_H

#include <cstdint>

#include <functional>
#include <memory>
#include <vector>

namespace gmx
{
template<typename>
class ArrayRef;

class ObservablesReducer;
using Step = int64_t;

/*! \brief Control whether reduction is required soon. */
enum class ReductionRequirement : int
{
    //! Reduce whenever the runner next checks with the ObservablesReducer.
    Soon,
    /*! \brief Reduce whenever the runner next checks with the
     * ObservablesReducer after some module requires reduction Soon */
    Eventually
};

/*! \brief Report whether the reduction has happened this step */
enum class ObservablesReducerStatus : int
{
    //! Reduction has not yet happened this step
    ReadyToReduce,
    //! Reduction has happened this step
    AlreadyReducedThisStep
};

/*! \libinternal \brief
 * Builder for ObservablesReducer
 *
 * Receives subscriptions from MD modules. Caller should call \c
 * build() once all subscriptions have been received, and then not
 * attempt any further subscriptions or builds. At that time, the
 * builder may be destructed.
 *
 * This builder will
 * - receive all subscriptions from MD modules, then
 * - build the communication buffer used by the subscribers,
 * - build the \c ObservablesReducer object that manages the
 *   lifetime of that buffer, and
 * - notify the subscribers via callback of the view of that buffer
 *   that is theirs to use and a further callback to require
 *   reduction of that buffer.
 * See also the interaction diagram in the documentation
 * for observablesreducer.h file.
 *
 * Note that the builder callbacks do not follow the approach of \c
 * MDModulesNotifier because that requires that the same value is
 * passed to all recipients. Here a distinct value goes to each
 * recipient, ie. a different view of the communication buffer.
 *
 * In order to avoid circular build-time dependencies between the
 * ObservablesReducer (and its builder) with the modules that use it,
 * the latter can directly call methods on the former, supplying
 * anonymous callbacks to be used by the former to contact the
 * latter. CallbackAfterReduction and CallbackFromBuilder are of this
 * type.
 *
 * A callback type CallBackToRequireReduction is also used instead of
 * a direct method call on ObservablesReducer to require reduction.
 * This is implemented by calling a method on the Impl object of a
 * ObservablesReducer. This extends the interface of
 * ObservablesReducer in a way that is not directly visible. That
 * complexity provides two benefits:
 * - only registered subscribers can require reduction (which helps
 *   ensure correctness by construction)
 * - the ObservablesReducer::Impl has a stable address from the heap
 *   allocation needed for std::unique_ptr to use in forming the
 *   callback to request reduction.
 * Alternatives exist for the latter, but create requirements on the
 * stability of the address of ObservablesReducer, and/or extra
 * coordination to only pass that address to subscribers once it is
 * stable.
 *
 * It is the subscribers' responsibility to coordinate so that all
 * subscribers on all ranks agree on the need to communicate, e.g. by
 * orchestrating communication based on the current step number or a
 * previous message.
 *
 */
class ObservablesReducerBuilder
{
public:
    //! Constructor
    ObservablesReducerBuilder();
    //! Destructor
    ~ObservablesReducerBuilder();
    //! Move constructor
    ObservablesReducerBuilder(ObservablesReducerBuilder&& other) noexcept;
    //! Move assignment operator
    ObservablesReducerBuilder& operator=(ObservablesReducerBuilder&& other) noexcept;

    /*! \brief Convenience type for the callback subscribers to
     * provide when they require reduction. */
    using CallbackAfterReduction = std::function<void(Step)>;
    /*! \brief Convenience type for the callback subscribers
     * call to require reduction.
     *
     * When called, the status it returns can be used for checking the
     * internal expectations of the subscriber on whether reduction
     * has already occurred this step, or not. */
    using CallbackToRequireReduction = std::function<ObservablesReducerStatus(ReductionRequirement)>;
    /*! \brief Convenience type for the callback from the builder to
     * notify the subscribers of the callback they will own and later
     * use to require reduction and the view of the communication
     * buffer they will later use. */
    using CallbackFromBuilder = std::function<void(CallbackToRequireReduction&&, ArrayRef<double>)>;

    /*! \brief Add a subscriber to the \c ObservablesReducer that will
     * later be built in \c build()
     *
     * Takes ownership of both callbacks supplied by the subscribing
     * module. This approach ensures that time is not spent in the MD
     * loop constructing std::function objects, because constructing
     * one of those requires 1-2 heap allocations (depending on the
     * size of the lambda capture).
     *
     * Must not be called after build() */
    void addSubscriber(int                      sizeRequired,
                       CallbackFromBuilder&&    callbackFromBuilder,
                       CallbackAfterReduction&& callbackAfterReduction);

    /*! \brief Build a \c ObservablesReducer to which any subscribers
     * have been added
     *
     * Must be called only once. Notifies each subscriber (via the
     * CallbackFromBuilder that it supplied) of the view of the
     * reduction buffer that they will use and the
     * CallbackToRequireReduction that they will use. */
    ObservablesReducer build();

private:
    class Impl;
    //! Impl object
    std::unique_ptr<Impl> impl_;
};

/*! \libinternal \brief
 * Manage reduction of observables for registered subscribers
 *
 * Modules can require that the \c ObservablesReducer object to which
 * they have subscribed do communication this step.  After reduction
 * is complete, notifications are made to the callbacks that modules
 * previously supplied to the ObservablesReducerBuilder. Then the
 * reduction buffer is zeroed. Thus the subscribers may not depend on
 * the values in their buffer view after the notification callback
 * returns, so they should do any necessary processing during that
 * callback.
 *
 * Modules are free to request reduction whenever they wish, and have
 * no obligations to do anything at any time. In particular, they
 * do not have to set values for their reduction buffer except when
 * they are requesting reduction.
 *
 * The \c ObservablesReducerBuilder object is responsible for
 * preparing a vector of doubles and notifying the subscribers of the
 * mutually disjoint views of the buffer that they should use for both
 * input and output of the reduction. The ObservablesReducer object
 * that it builds retains no record of the subscribers, because its
 * responsibility is solely to orchestrate the MPI communication and
 * callbacks.
 *
 * Subscribers automatically use the correct \c ObservablesReducer
 * object because the callback they received is bound to the correct
 * one. The only way a module can participate in an \c
 * ObservablesReducer is to have registered with its builder.
 *
 * The owner of an ObservablesReducer must maintain the lifetime of
 * the \c ObservablesReducer object until all subscribers no longer
 * need it. After the destruction of an \c ObservablesReducer, if
 * subscribers access their view of the communication buffer, the
 * behavior is undefined.
 *
 * \inlibraryapi
 * \ingroup module_mdtypes
 */
class ObservablesReducer
{
private:
    class Impl;
    std::unique_ptr<Impl> impl_;

public:
    //! Constructor only usable by ObservablesReducerBuilder
    explicit ObservablesReducer(std::unique_ptr<Impl> impl);
    // Destructor
    ~ObservablesReducer();
    //! Move constructor
    ObservablesReducer(ObservablesReducer&& other) noexcept;
    //! Move assignment operator
    ObservablesReducer& operator=(ObservablesReducer&& other) noexcept;

    /*! \brief Return whether any subscriber has required reduction soon.
     */
    bool isReductionRequired() const;
    /*! \brief Provide view of communication buffer for MPI reduction
     *
     * If no subscriber used ReductionRequirement::Soon since the last
     * call to reductionComplete(), then this method returns an empty
     * buffer. Otherwise it returns a view over the buffer potentially
     * filled by all subscribed modules.
     *
     * For so long as the ObservablesReducer continues to interact
     * with the legacy compute_globals() code, the implementation of
     * ReductionRequirement::Eventually needs to know whether any
     * modules using compute_globals() have already requested a
     * reduction. This value is passed as a parameter to this method.
     *
     * \param[in] reductionRequiredExternally Whether external code
     *                                        has required a reduction
     */
    ArrayRef<double> communicationBuffer(bool reductionRequiredExternally);
    /*! \brief Called by the runner after MPI communication is complete
     *
     * Notifies all subscribers who required reduction since the last
     * call to reductionComplete() and passes the \c step value so
     * they can check internally that the simulation state is
     * consistent.
     *
     * After all notifications, zeroes the communication buffer. It is
     * the responsibility of the subscribers that required reduction
     * to react suitably to the data available during their
     * notification. This ensures that modules cannot get arbitrary
     * but realistic-looking values left behind from previous
     * communication stages. It also ensures that subsequent
     * communication stages will not be able to keep reducing values
     * until they overflow or underflow. This zeroing is most efficient
     * to do centrally in an object of this class.
     *
     * The choice of zero for the sentinel value is not perfect. In
     * principle, a value of zero is potentially significant to any
     * subscriber, so could be provided to a subscriber as the result
     * of an incorrect implementation of ObservablesReducer or
     * inconsistent use by subscribers. However by construction (which
     * is tested), the integration tests never produce a zero result
     * from an reduced value provided by a subscriber. So, if the
     * coverage there is high then there is good reason to expect that
     * when a zero value is used by a subscribers it is the result of
     * a reduction and thus significant, rather than an artefact of
     * the zeroing of the communication buffer after notifications are
     * complete.
     *
     * The choice of zero ensures that the MPI reduction will produce
     * a valid numerical result in all cases except when a module that
     * required reduction set buffer contents that produced a
     * problematic output after reduction.
     */
    void reductionComplete(Step step);
    /*! \brief Notify the ObservablesReducer that it should make
     * ready to receive new values to reduce
     *
     * Any runner using the ObservablesReducer must call this method
     * whenever a step completes, so that subscribed modules can use
     * that information to check whether reduction is happening on the
     * step that they expect.
     *
     * The ObservablesReducer keeps track of whether reduction has
     * already occurred this step, so that when modules request
     * reduction it can notify them of that status. This permits them
     * to check their own requirements, e.g. that
     * ReductionRequirement::Soon will operate this step or next
     * step.
     *
     * For the same reason, it is also necessary to call this method
     * at a suitable point after uses of an ObservablesReducer before
     * the regular steps of a runner. */
    void markAsReadyToReduce();
    //! The builder needs to be able to make the Impl object
    friend class ObservablesReducerBuilder;
};

} // namespace gmx

#endif
