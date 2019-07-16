/*! \file
 * \brief Provide some useful types and templates for GROMACS restraints.
 *
 * \todo This should be part of a template library installed with GROMACS.
 *
 * \author M. Eric Irrgang <ericirrgang@gmail.com>
 */

#ifndef RESTRAINT_SESSIONRESOURCES_H
#define RESTRAINT_SESSIONRESOURCES_H

#include <functional>
#include <memory>
#include <mutex>
#include <vector>

#include "gmxapi/gromacsfwd.h"
#include "gmxapi/session.h"
#include "gmxapi/session/resources.h"
#include "gmxapi/md/mdmodule.h"

#include "gromacs/restraint/restraintpotential.h"
#include "gromacs/utility/real.h"

namespace plugin
{

// Stop-gap for cross-language data exchange pending SharedData implementation and inclusion of Eigen.
// Adapted from pybind docs.
template<class T>
class Matrix
{
    public:
        Matrix(size_t rows,
               size_t cols) :
            rows_(rows),
            cols_(cols),
            data_(rows_ * cols_,
                  0)
        {
        }

        explicit Matrix(std::vector<T>&& captured_data) :
            rows_{1},
            cols_{captured_data.size()},
            data_{std::move(captured_data)}
        {
        }

        std::vector<T>* vector()
        { return &data_; }

        T* data()
        { return data_.data(); };

        size_t rows() const
        { return rows_; }

        size_t cols() const
        { return cols_; }

    private:
        size_t rows_;
        size_t cols_;
        std::vector<T> data_;
};

// Defer implicit instantiation to ensemblepotential.cpp
extern template
class Matrix<double>;

/*!
 * \brief An active handle to ensemble resources provided by the Context.
 *
 * The semantics of holding this handle aren't determined yet, but it should be held as briefly as possible since it
 * may involve locking global resources or preventing the simulation from advancing. Basically, though, it allows the
 * Context implementation flexibility in how or where it provides services.
 *
 * Resources may be incoming input data or functors to trigger output data events.
 *
 * \internal
 * It is not yet clear whether we want to assume that default behavior is for an operation to be called for each edge
 * on every iterative graph execution, leaving less frequent calls as an optimization, or to fully support occasional
 * data events issued by nodes during their execution.
 *
 * In this example, assume the plugin has specified that it provides a `.ostream.stop` port that provides asynchronous
 * boolean events. We can provide a template member function that will handle either execution mode.
 *
 * ResourceHandle::ostream() will return access to a gmxapi::session::OutputStream object, which will provide
 * set("stop", true), to give access to a function pointer from a member vector of function pointers.
 *
 * In the case that we are triggering asynchronous data events, the function will make the appropriate call. In the case
 * that we have output at regular intervals, the function will update internal state for the next time the edge is
 * evaluated.
 *
 * In an alternative implementation, we could maintain a data object that could be queried by subscribers, but a publish
 * and subscribe model seems much more useful, optimizeable, and robust. We can issue the calls to the subscribers and
 * then be done with it.
 *
 * If we need to optimize for reuse of memory locations, we can do one of two things: require that
 * the subscribing object not return until it has done what it needed with the data (including deep copy) or use
 * managed handles for the data, possibly with a custom allocator, that prevents rewriting while there are read handles
 * still open. One way that came up in conversation with Mark to allow some optimization is to allow the recipient of
 * the handle to make either an `open` that gets a potentially blocking read-lock or an `open` that requests ownership.
 * If no other consumers of the data request ownership, the ownership can be transferred without a copy. Otherwise, a
 * copy is made.
 */
class ResourcesHandle
{
    public:
        /*!
         * \brief Ensemble reduce.
         *
         * \param send Matrices to be summed across the ensemble using Context resources.
         * \param receive destination of reduced data instead of updating internal Matrix.
         */
        void reduce(const Matrix<double>& send,
                    Matrix<double>* receive) const;

        /*!
         * \brief Issue a stop condition event.
         *
         * Can be called on any or all ranks. Sets a condition that will cause the current simulation to shut down
         * after the current step.
         */
        void stop();

        // to be abstracted and hidden...
        const std::function<void(const Matrix<double>&,
                                 Matrix<double>*)>* reduce_;

        gmxapi::SessionResources* session_;
};

/*!
 * \brief Reference to workflow-level resources managed by the Context.
 *
 * Provides a connection to the higher-level workflow management with which to access resources and operations. The
 * reference provides no resources directly and we may find that it should not extend the life of a Session or Context.
 * Resources are accessed through Handle objects returned by member functions.
 *
 * gmxapi version 0.1.0 will provide this functionality through SessionResources.
 */
class Resources
{
    public:
        /*!
         * \brief Create a new resources object.
         *
         * This constructor is called by the framework during Session launch to provide the plugin
         * potential with external resources.
         *
         * \note If getHandle() is going to be used, setSession() must be called first.
         *
         * \param reduce ownership of a function object providing ensemble averaging of a 2D matrix.
         */
        explicit Resources(std::function<void(const Matrix<double>&,
                                              Matrix<double>*)>&& reduce) :
            reduce_(reduce),
            session_(nullptr)
        {};

        /*!
         * \brief Grant the caller an active handle for the currently executing block of code.
         *
         * Objects should not keep resource handles open for longer than a single block of code.
         * calculate() and callback() functions get a handle to the resources for the current time step
         * by calling getHandle().
         *
         * \note setSession() must be called before this function can be used.
         * This clumsy protocol requires other infrastructure before it can be
         * cleaned up for gmxapi 0.1
         *
         * \return resource handle
         *
         * In this release, the only facility provided by the resources is a function object for
         * the ensemble averaging function provided by the Context.
         */
        ResourcesHandle getHandle() const;

        /*!
         * \brief Acquires a pointer to a Session managing these resources.
         *
         * \param session non-owning pointer to Session resources.
         */
        void setSession(gmxapi::SessionResources* session);

    private:
        //! bound function object to provide ensemble reduce facility.
        std::function<void(const Matrix<double>&,
                           Matrix<double>*)> reduce_;

        // Raw pointer to the session in which these resources live.
        gmxapi::SessionResources* session_;
};

/*!
 * \brief Template for MDModules from restraints.
 *
 * Allows a GROMACS module to be produced easily from the provided class. Refer to
 * src/pythonmodule/export_plugin.cpp for how this template is used.
 *
 * \tparam R a class implementing the gmx::IRestraintPotential interface.
 *
 * The template type parameter should define a ``input_param_type`` member type.
 *
 * \todo move this to a template header in gmxapi */
template<class R>
class RestraintModule : public gmxapi::MDModule // consider names
{
    public:
        using param_t = typename R::input_param_type;

        /*!
         * \brief Construct a named restraint module.
         *
         * Objects of this type are created during Session launch, so this code really doesn't belong
         * here. The Director / Builder for the restraint uses a generic interface to pass standard
         * parameters for pair restraints: a list of sites, a (custom) parameters structure, and
         * resources provided by the Session.
         *
         * \param name
         * \param sites
         * \param params
         * \param resources
         */
        RestraintModule(std::string name,
                        std::vector<int> sites,
                        const typename R::input_param_type& params,
                        std::shared_ptr<Resources> resources) :
            sites_{std::move(sites)},
            params_{params},
            resources_{std::move(resources)},
            name_{std::move(name)}
        {

        };

        ~RestraintModule() override = default;

        /*!
         * \brief Implement gmxapi::MDModule interface to get module name.
         *
         * name is provided during the building stage.
         * \return
         */
        // \todo make member function const
        const char* name() const override
        {
            return name_.c_str();
        }

        /*!
         * \brief Implement gmxapi::MDModule interface to create a restraint for libgromacs.
         *
         * \return (Possibly shared) Ownership of a restraint instance
         *
         * Creates the restraint instance if it does not already exist. Only creates one restraint
         * instance in the lifetime of the RestraintModule.
         * 
         * Note this interface is not stable but requires other GROMACS and gmxapi infrastructure
         * to mature before it is clear whether we will be creating a new instance or sharing ownership
         * of the object. A future version may use a std::unique_ptr.
         */
        std::shared_ptr<gmx::IRestraintPotential> getRestraint() override
        {
            std::lock_guard<std::mutex> lock(restraintInstantiation_);
            if (!restraint_)
            {
                restraint_ = std::make_shared<R>(sites_,
                                                 params_,
                                                 resources_);
            }
            return restraint_;
        }

    private:
        std::vector<int> sites_;
        param_t params_;

        // Need to figure out if this is copyable or who owns it.
        std::shared_ptr<Resources> resources_;

        const std::string name_;
        std::shared_ptr<R> restraint_{nullptr};
        std::mutex restraintInstantiation_;
};

/*!
 * \brief Filehandle management helper class.
 *
 * Use the RAII pattern to make sure that a (newly) constructed object has an open filehandle and that
 * a the filehandle for a destructed object is closed. Closing a file is not guaranteed to be error-free,
 * so the programmer should explicitly call close() and check for errors (see the C library docs
 * for fclose()).
 *
 * RAIIFile makes sure that fclose() is called exactly once, whether client code issues close()
 * or not.
 */
class RAIIFile
{
    public:

        /*!
         * \brief Open a file in the chosen access mode.
         *
         * \param filename Name of file to be opened.
         * \param mode access mode as described for the fopen C library call.
         */
        RAIIFile(const char* filename,
                 const char* mode) :
            fh_{fopen(filename,
                      mode)}
        {}

        /*!
         * \brief Open a file for writing.
         *
         * \param filename Name of file to be opened.
         *
         * File is opened in mode "w", which truncates data if the file already exists.
         * For other file access modes, use RAIIFile(const char* filename, const char* mode)
         */
        explicit RAIIFile(const char* filename) :
            RAIIFile(filename,
                     "w")
        {}

        /*!
         * \brief Explicitly close the associated filehandle.
         *
         * It is good practice to explicitly close the file at a known point in the client code, though
         * it is not strictly necessary. If the filehandle is still open when the RAIIFile object is
         * destroyed, the fclose will be called then.
         *
         * Calling close() additional times on the same RAIIFile object is fine and has no effect in
         * single-threaded code. However, the destructor and close() routines are not thread-safe, so
         * the client code should make sure that close() is not called at the same time by multiple threads.
         * Standard reference-counting constructs, like std::shared_ptr, can be used to make sure the
         * object destructor is called exactly once if it needs to be shared.
         *
         * Refer to documentation on fclose() on checking for and interpreting `errno`.
         */
        void close()
        {
            if (fh_ != nullptr)
            {
                fclose(fh_);
            }
            fh_ = nullptr;
        }

        /*!
         * \brief RAII destructor.
         *
         * Make sure the filehandle gets closed exactly once.
         */
        ~RAIIFile()
        {
            if (fh_ != nullptr)
            {
                fclose(fh_);
            }
        }

        RAIIFile(const RAIIFile&) = delete;

        RAIIFile& operator=(const RAIIFile&) = delete;

        RAIIFile(RAIIFile&&) = default;

        RAIIFile& operator=(RAIIFile&&) = default;

        /*!
         * \brief Get the managed filehandle.
         *
         * \return raw pointer to the underlying filehandle.
         */
        FILE* fh() const noexcept
        {
            return fh_;
        }

    private:
        /// file handle
        FILE* fh_{nullptr};
};

} // end namespace plugin

#endif //RESTRAINT_SESSIONRESOURCES_H
