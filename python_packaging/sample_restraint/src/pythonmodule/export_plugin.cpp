/*! \file
 * \brief Provide Python bindings and helper functions for setting up restraint potentials.
 *
 * There is currently a lot of boilerplate here that will be generalized and removed in a future
 * version. In the mean time, follow the example for EnsembleRestraint to create the proper helper
 * functions and instantiate the necessary templates.
 *
 * \author M. Eric Irrgang <ericirrgang@gmail.com>
 */

#include "export_plugin.h"

#include <cassert>

#include <memory>

#include "gmxapi/exceptions.h"
#include "gmxapi/md.h"
#include "gmxapi/md/mdmodule.h"
#include "gmxapi/gmxapi.h"

#include "ensemblepotential.h"
#include "nullpotential.h"

// Make a convenient alias to save some typing...
namespace py = pybind11;

////////////////////////////////
// Begin PyRestraint static code
/*!
 * \brief Templated wrapper to use in Python bindings.
 *
 * Boilerplate
 *
 * Mix-in from below. Adds a bind behavior, a getModule() method to get a gmxapi::MDModule adapter,
 * and a create() method that assures a single shared_ptr record for an object that may sometimes
 * be referred to by a raw pointer and/or have shared_from_this called.
 * \tparam T class implementing gmx::IRestraintPotential
 *
 */
template<class T>
class PyRestraint : public T, public std::enable_shared_from_this<PyRestraint<T>>
{
public:
    void bind(py::object object);

    using T::name;

    /*!
     * \brief
     *
     * T must either derive from gmxapi::MDModule or provide a template specialization for
     * PyRestraint<T>::getModule(). If T derives from gmxapi::MDModule, we can keep a weak pointer
     * to ourself and generate a shared_ptr on request, but std::enable_shared_from_this already
     * does that, so we use it when we can.
     * \return
     */
    std::shared_ptr<gmxapi::MDModule> getModule();

    /*!
     * \brief Factory function to get a managed pointer to a new restraint.
     *
     * \tparam ArgsT
     * \param args
     * \return
     */
    template<typename... ArgsT>
    static std::shared_ptr<PyRestraint<T>> create(ArgsT... args)
    {
        auto newRestraint = std::make_shared<PyRestraint<T>>(args...);
        return newRestraint;
    }

    template<typename... ArgsT>
    explicit PyRestraint(ArgsT... args) : T{ args... }
    {
    }
};

/*!
 * \brief Implement the gmxapi binding protocol for restraints.
 *
 * All restraints will use this same code automatically.
 *
 * \tparam T restraint class exported below.
 * \param object Python Capsule object to allow binding with a simple C API.
 */
template<class T>
void PyRestraint<T>::bind(py::object object)
{
    PyObject* capsule = object.ptr();
    if (PyCapsule_IsValid(capsule, gmxapi::MDHolder::api_name))
    {
        auto holder =
                static_cast<gmxapi::MDHolder*>(PyCapsule_GetPointer(capsule, gmxapi::MDHolder::api_name));
        auto workSpec = holder->getSpec();
        std::cout << this->name() << " received " << holder->name();
        std::cout << " containing spec of size ";
        std::cout << workSpec->getModules().size();
        std::cout << std::endl;

        auto module = getModule();
        workSpec->addModule(module);
    }
    else
    {
        throw gmxapi::ProtocolError("bind method requires a python capsule as input");
    }
}
// end PyRestraint static code
//////////////////////////////


/*!
 * \brief Interact with the restraint framework and gmxapi when launching a simulation.
 *
 * This should be generalized and removed from here. Unfortunately, some things need to be
 * standardized first. If a potential follows the example of EnsembleRestraint or HarmonicRestraint,
 * the template specializations below can be mimicked to give GROMACS access to the potential.
 *
 * \tparam T class implementing the gmxapi::MDModule interface.
 * \return shared ownership of a T object via the gmxapi::MDModule interface.
 */
// If T is derived from gmxapi::MDModule, create a default-constructed std::shared_ptr<T>
// \todo Need a better default that can call a shared_from_this()
template<class T>
std::shared_ptr<gmxapi::MDModule> PyRestraint<T>::getModule()
{
    auto module =
            std::make_shared<typename std::enable_if<std::is_base_of<gmxapi::MDModule, T>::value, T>::type>();
    return module;
}


template<>
std::shared_ptr<gmxapi::MDModule> PyRestraint<plugin::RestraintModule<plugin::EnsembleRestraint>>::getModule()
{
    return shared_from_this();
}
//////////////////////////////////////////////////////////////////////////////////////////
// New restraints mimicking EnsembleRestraint should specialize getModule() here as above.
//////////////////////////////////////////////////////////////////////////////////////////

//! Provide the hook used by gmxapi::addSessionRestraint().
template<>
std::shared_ptr<gmxapi::MDModule> PyRestraint<plugin::RestraintModule<plugin::NullRestraint>>::getModule()
{
    return shared_from_this();
}


class EnsembleRestraintBuilder
{
public:
    explicit EnsembleRestraintBuilder(const py::object& element)
    {
        name_ = py::cast<std::string>(element.attr("name"));
        assert(!name_.empty());

        // It looks like we need some boilerplate exceptions for plugins so we have something to
        // raise if the element is invalid.
        assert(py::hasattr(element, "params"));

        // Params attribute should be a Python list
        py::dict parameter_dict = element.attr("params");
        // \todo Check for the presence of these dictionary keys to avoid hard-to-diagnose error.

        // Get positional parameters.
        py::list sites = parameter_dict["sites"];
        for (auto&& site : sites)
        {
            siteIndices_.emplace_back(py::cast<int>(site));
        }

        auto nbins        = py::cast<size_t>(parameter_dict["nbins"]);
        auto binWidth     = py::cast<double>(parameter_dict["binWidth"]);
        auto minDist      = py::cast<double>(parameter_dict["min_dist"]);
        auto maxDist      = pybind11::cast<double>(parameter_dict["max_dist"]);
        auto experimental = pybind11::cast<std::vector<double>>(parameter_dict["experimental"]);
        auto nSamples     = pybind11::cast<unsigned int>(parameter_dict["nsamples"]);
        auto samplePeriod = pybind11::cast<double>(parameter_dict["sample_period"]);
        auto nWindows     = pybind11::cast<unsigned int>(parameter_dict["nwindows"]);
        auto k            = pybind11::cast<double>(parameter_dict["k"]);
        auto sigma        = pybind11::cast<double>(parameter_dict["sigma"]);

        auto params = plugin::makeEnsembleParams(
                nbins, binWidth, minDist, maxDist, experimental, nSamples, samplePeriod, nWindows, k, sigma);
        params_ = std::move(*params);

        // Note that if we want to grab a reference to the Context or its communicator, we can get
        // it here through element.workspec._context. We need a more general API solution, but this
        // code is in the Python bindings code, so we know we are in a Python Context.
        assert(py::hasattr(element, "workspec"));
        auto workspec = element.attr("workspec");
        assert(py::hasattr(workspec, "_context"));
        context_ = workspec.attr("_context");
    }

    /*!
     * \brief Add node(s) to graph for the work element.
     *
     * \param graph networkx.DiGraph object still evolving in gmx.context.
     *
     * \todo This may not follow the latest graph building protocol as described.
     */
    void build(const py::object& graph)
    {
        if (!subscriber_)
        {
            return;
        }
        else
        {
            if (!py::hasattr(subscriber_, "potential"))
                throw gmxapi::ProtocolError("Invalid subscriber");
        }

        // Restraints do not currently add any new nodes to the graph, so we
        // mark this standard 'graph' argument unused.
        (void)graph;

        // Temporarily subvert things to get quick-and-dirty solution for testing.
        // Need to capture Python communicator and pybind syntax in closure so EnsembleResources
        // can just call with matrix arguments.

        // This can be replaced with a subscription and delayed until launch, if necessary.
        if (!py::hasattr(context_, "ensemble_update"))
        {
            throw gmxapi::ProtocolError("context does not have 'ensemble_update'.");
        }
        // make a local copy of the Python object so we can capture it in the lambda
        auto update = context_.attr("ensemble_update");
        // Make a callable with standardizeable signature.
        const std::string name{ name_ };
        auto functor = [update, name](const plugin::Matrix<double>& send, plugin::Matrix<double>* receive)
        { update(send, receive, py::str(name)); };

        // To use a reduce function on the Python side, we need to provide it with a Python buffer-like object,
        // so we will create one here. Note: it looks like the SharedData element will be useful after all.
        auto resources = std::make_shared<plugin::Resources>(std::move(functor));

        auto potential = PyRestraint<plugin::RestraintModule<plugin::EnsembleRestraint>>::create(
                name_, siteIndices_, params_, resources);

        auto     subscriber    = subscriber_;
        py::list potentialList = subscriber.attr("potential");
        potentialList.append(potential);
    };

    /*!
     * \brief Accept subscription of an MD task.
     *
     * \param subscriber Python object with a 'potential' attribute that is a Python list.
     *
     * During build, an object is added to the subscriber's self.potential, which is then bound with
     * system.add_potential(potential) during the subscriber's launch()
     */
    void addSubscriber(const py::object& subscriber)
    {
        assert(py::hasattr(subscriber, "potential"));
        subscriber_ = subscriber;
    };

    py::object       subscriber_;
    py::object       context_;
    std::vector<int> siteIndices_;

    plugin::ensemble_input_param_type params_;

    std::string name_;
};

class NullRestraintBuilder
{
public:
    explicit NullRestraintBuilder(const py::object& element)
    {
        name_ = py::cast<std::string>(element.attr("name"));

        if (name_.empty())
        {
            throw gmxapi::ProtocolError("Restraint must provide a *name*.");
        }
        if (!py::hasattr(element, "params"))
        {
            throw gmxapi::ProtocolError("Invalid WorkflowElement. (missing *params*)");
        }

        // Params attribute should be a Python list
        py::dict parameter_dict = element.attr("params");

        // Get positional parameters.
        py::list sites = parameter_dict["sites"];
        for (auto&& site : sites)
        {
            siteIndices_.emplace_back(py::cast<int>(site));
        }

        params_ = plugin::makeNullParams(std::vector<int>(siteIndices_));

        // Note that if we want to grab a reference to the Context or its communicator, we can get
        // it here through element.workspec._context. We need a more general API solution, but this
        // code is in the Python bindings code, so we know we are in a Python Context.
        assert(py::hasattr(element, "workspec"));
        auto workspec = element.attr("workspec");
        assert(py::hasattr(workspec, "_context"));
        context_ = workspec.attr("_context");
    }

    /*!
     * \brief Add node(s) to graph for the work element.
     *
     * \param graph networkx.DiGraph object still evolving in gmx.context.
     *
     * \todo This may not follow the latest graph building protocol as described.
     */
    void build([[maybe_unused]] const py::object& graph)
    {
        if (!subscriber_)
        {
            return;
        }
        else
        {
            if (!py::hasattr(subscriber_, "potential"))
                throw gmxapi::ProtocolError("Invalid subscriber");
        }

        // This restraint does not need any session resources.
        auto null_resources = std::make_shared<plugin::Resources>(
                [](const plugin::Matrix<double>&, plugin::Matrix<double>*) {});

        auto potential = PyRestraint<plugin::RestraintModule<plugin::NullRestraint>>::create(
                name_, siteIndices_, params_, std::move(null_resources));

        auto     subscriber    = subscriber_;
        py::list potentialList = subscriber.attr("potential");
        potentialList.append(potential);
    };

    /*!
     * \brief Accept subscription of an MD task builder.
     *
     * \param subscriber Python object with a 'potential' attribute that is a Python list.
     *
     * During build, an object is added to the subscriber's self.potential, which is then bound with
     * system.add_potential(potential) during the subscriber's launch()
     */
    void addSubscriber(const py::object& subscriber)
    {
        assert(py::hasattr(subscriber, "potential"));
        subscriber_ = subscriber;
    };

    py::object       subscriber_;
    py::object       context_;
    std::vector<int> siteIndices_;

    plugin::NullRestraint::input_param_type params_;

    std::string name_;
};


namespace
{

/*!
 * \brief Factory function to create a new builder for use during Session launch.
 *
 * \param element WorkElement provided through Context
 * \return ownership of new builder object
 */
std::unique_ptr<EnsembleRestraintBuilder> createEnsembleBuilder(const py::object& element)
{
    using std::make_unique;
    auto builder = make_unique<EnsembleRestraintBuilder>(element);
    return builder;
}

/*!
 * \brief Factory function to create a new builder for use during Session launch.
 *
 * \param element WorkElement provided through Context
 * \return ownership of new builder object
 */
std::unique_ptr<NullRestraintBuilder> createNullRestraintBuilder(const py::object& element)
{
    using std::make_unique;
    auto builder = make_unique<NullRestraintBuilder>(element);
    return builder;
}


} // namespace


////////////////////////////////////////////////////////////////////////////////////////////
// New potentials modeled after EnsembleRestraint should define a Builder class and define a
// factory function here, following the previous two examples. The factory function should be
// exposed to Python following the examples near the end of the PYBIND11_MODULE block.
////////////////////////////////////////////////////////////////////////////////////////////


//////////////////////////////////////////////////////////////////////////////////////////////////
// The PYBIND11_MODULE block uses the pybind11 framework (ref https://github.com/pybind/pybind11 )
// to generate Python bindings to the C++ code elsewhere in this repository. A copy of the pybind11
// source code is included with this repository. Use syntax from the examples below when exposing
// a new potential, along with its builder and parameters structure. In future releases, there will
// be less code to include elsewhere, but more syntax in the block below to define and export the
// interface to a plugin. pybind11 is not required to write a GROMACS extension module or for
// compatibility with the ``gmx`` module provided with gmxapi. It is sufficient to implement the
// various protocols, C API and Python function names, but we do not provide example code
// for other Python bindings frameworks.
//////////////////////////////////////////////////////////////////////////////////////////////////

// The first argument is the name of the module when importing to Python. This should be the same as the name specified
// as the OUTPUT_NAME for the shared object library in the CMakeLists.txt file. The second argument, 'm', can be anything
// but it might as well be short since we use it to refer to aspects of the module we are defining.
PYBIND11_MODULE(myplugin, m)
{
    m.doc() = "sample plugin"; // This will be the text of the module's docstring.

    // Matrix utility class (temporary). Borrowed from http://pybind11.readthedocs.io/en/master/advanced/pycpp/numpy.html#arrays
    py::class_<plugin::Matrix<double>, std::shared_ptr<plugin::Matrix<double>>>(
            m, "Matrix", py::buffer_protocol())
            .def_buffer(
                    [](plugin::Matrix<double>& matrix) -> py::buffer_info
                    {
                        return py::buffer_info(
                                matrix.data(),                           /* Pointer to buffer */
                                sizeof(double),                          /* Size of one scalar */
                                py::format_descriptor<double>::format(), /* Python struct-style format descriptor */
                                2,                                       /* Number of dimensions */
                                { matrix.rows(), matrix.cols() }, /* Buffer dimensions */
                                { sizeof(double) * matrix.cols(), /* Strides (in bytes) for each index */
                                  sizeof(double) });
                    });

    //////////////////////////////////////////////////////////////////////////
    // Begin EnsembleRestraint
    //
    // Define Builder to be returned from ensemble_restraint Python function defined further down.
    pybind11::class_<EnsembleRestraintBuilder> ensembleBuilder(m, "EnsembleBuilder");
    ensembleBuilder.def("add_subscriber", &EnsembleRestraintBuilder::addSubscriber);
    ensembleBuilder.def("build", &EnsembleRestraintBuilder::build);

    // Get more concise name for the template instantiation...
    using PyEnsemble = PyRestraint<plugin::RestraintModule<plugin::EnsembleRestraint>>;

    // Export a Python class for our parameters struct
    py::class_<plugin::EnsembleRestraint::input_param_type> ensembleParams(
            m, "EnsembleRestraintParams");
    m.def("make_ensemble_params", &plugin::makeEnsembleParams);

    // API object to build.
    py::class_<PyEnsemble, std::shared_ptr<PyEnsemble>> ensemble(m, "EnsembleRestraint");
    // EnsembleRestraint can only be created via builder for now.
    ensemble.def("bind", &PyEnsemble::bind, "Implement binding protocol");
    /*
     * To implement gmxapi_workspec_1_0, the module needs a function that a Context can import that
     * produces a builder that translates workspec elements for session launching. The object
     * returned by our function needs to have an add_subscriber(other_builder) method and a
     * build(graph) method. The build() method returns None or a launcher. A launcher has a
     * signature like launch(rank) and returns None or a runner.
     */

    // Generate the name operation that will be used to specify elements of Work in gmxapi workflows.
    // WorkElements will then have namespace: "myplugin" and operation: "ensemble_restraint"
    m.def("ensemble_restraint",
          [](const py::object& element) { return createEnsembleBuilder(element); });
    //
    // End EnsembleRestraint
    ///////////////////////////////////////////////////////////////////////////

    ///////////////////////////////////////////////////////////////////////////
    // Begin NullRestraint
    //
    // Define the Builder Python interface.
    py::class_<NullRestraintBuilder> nullRestraintBuilder(m, "NullBuilder");
    nullRestraintBuilder.def("add_subscriber", &NullRestraintBuilder::addSubscriber);
    nullRestraintBuilder.def("build", &NullRestraintBuilder::build);

    // Get concise name for the template instantiation.
    using PyNullR = PyRestraint<plugin::RestraintModule<plugin::NullRestraint>>;

    // Export a Python class for our data structure.
    py::class_<plugin::NullRestraint::input_param_type> nullParams(m, "NullRestraintParams");
    m.def("make_null_params", &plugin::makeNullParams);

    // Describe the Python API object that is built.
    py::class_<PyNullR, std::shared_ptr<PyNullR>> nullRestraint(m, "NullRestraint");
    nullRestraint.def("bind", &PyNullR::bind, "Implement binding protocol.");
    // We need a protocol for interacting with pluggable extension code and its data.
    // See #3038, #3133, #3145, #4079.
    nullRestraint.def(
            "count",
            [](PyNullR* restraint)
            {
                return plugin::count(
                        dynamic_cast<plugin::NullRestraint*>(restraint->getRestraint().get())->data_);
            });

    // Export the factory method that will resolve for
    // {namespace: "myplugin", operation: "null_restraint"} in a gmxapi WorkElement.
    m.def("null_restraint",
          [](const py::object& element) { return createNullRestraintBuilder(element); });
    //
    // End NullRestraint
    ///////////////////////////////////////////////////////////////////////////
}
