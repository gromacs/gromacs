Extensibility of pull code {#page_pullcodeextension}
==========================

\tableofcontents

External restraint forces can be added to \Gromacs through the <a href="group__module__restraint.xhtml">restraint</a> framework.

Historically, bias potentials have been applied by patching \Gromacs at the pull code.

This approach introduces an immediate difficulty in staying up to date with other
developments in \Gromacs and leaves many problems for the researcher, such as how to
get user input into the additional code and how to transfer data in and out of the
extension code.

\section case_study Case study


2015 \Gromacs fork adding an ensemble biasing scheme after Roux, dx.doi.org/10.1021/jp3110369

The potential applies an additional harmonic force along selected particle-pair separation vectors
to bias the simulation towards experimental sampling distributions.
Until a trajectory is produced that samples the desired conformation space,
pair separation sampling is analyzed across an ensemble of simulations to periodically
refine the biasing potential.

The code (ab)uses the GROMACS "pull" machinery not specifically for free energy calculation,
but because the large length scales and small number of relevant coordinates are a
good match for the computational work managed.

Commit ac9ce76

Summary of changes:

    +11  -0  M docs/user-guide/mdp-options.rst
    +4   -0  M src/gromacs/fileio/tpxio.cpp
    +45  -3  M src/gromacs/gmxana/gmx_mindist.cpp
    +2   -2  M src/gromacs/gmxlib/names.cpp
    +7   -1  M src/gromacs/gmxpreprocess/readpull.cpp
    +3   -2  M src/gromacs/legacyheaders/types/enums.h
    +298 -0  A src/gromacs/pulling/csvparser.cpp
    +48  -0  A src/gromacs/pulling/csvparser.h
    +75  -27 M src/gromacs/pulling/pull.cpp
    +1   -0  M src/gromacs/pulling/pullutil.cpp
    +113 -0  A src/gromacs/pulling/roux.cpp
    +9   -0  A src/gromacs/pulling/roux.h
    +1   -1  M src/programs/mdrun/md.cpp

First, note that over 300 lines are added for processing comma-separated-value data for the extension code.
This sort of dependency certainly should not need to be satisfied in the `libgromacs` code.

Second, note that many of these changes can be separated into a separate unrelated change
for `epullgDISTREF`.

\subsection roux_inputs Roux inputs

    // method parameters
    bin_width = 0.1
    sigma = 0.1
    max_dist = 9.0
    
    // working data (runtime environment)
    path_raw_data = /home/jmh/pmkResearch/exp_DEER

    // system input
    deer_file = experimental.csv
    
    // system parameters
    pair_file = pairs.txt



\subsection globals Globals modified


`gmxlib/names.cpp`

* add "roux" to `const char *epull_names[]`
* add "distance-reference" to `const char *epullg_names[]`

`/legacyheaders/types/enums.h`

* add `epullROUX`
* add `epullgDISTREF`

\subsection functions_modified Functions modified


`fileio/tpxio.cpp`

* `static void do_pull_coord()`
  * call `gmx_fio_do_int(fio, pcrd->group[2])` for `pcrd->eGeom == epullgDIRRELATIVE`

`gmxana/gmx_mindist.cpp`

* `void dist_plot()`
  * add `int rectmat` call parameter
  * handle rectangular distance matrix
* `int gmx_mindist()`
  * additional data member and CLI argument
  * handle new call signature for `dist_plot()`

`gmxpreprocess/readpull.cpp`

* `static void init_pull_coord()`
  * extend scope of calculation to apply for `pcrd->eGeom == epullgDISTREF`
* `char **read_pullparams()`
  * set `ngroup = 3` when `pcrd->eGeom == epullgDISTREF`
* `void make_pull_coords()`
  * Add a line of status output when `pcrd->eGeom == epullgDISTREF`

`pulling/pull.cpp`

* `static void low_get_pull_coord_dr()`
  * add `dvec xref2` to call signature
  * change how `dr[m]` and `dr2` are set when `pcrd->params.eGeom == epullgDISTREF`
  * extend scope of calculation to also apply when `pcrd->params.eGeom == epullgDISTREF`
* `static void get_pull_coord_dr()`
  * add a small data member
  * handle new call signature for `low_get_pull_coord_dr()`
* `static void get_pull_coord_distance()`
  * extend scope of calculation to also apply when `pcrd->params.eGeom == epullgDISTREF`
* `static double get_pull_coord_deviation()`
  * extend scope of calculations to an additional case `epullgDISTREF`
* `static void do_constraint()`
  * add a small data member
  * adapt to new call signature of `low_get_pull_coord_dr()`
  * extend scope of calculations to an additional case `epullgDISTREF`
* `static void calc_pull_coord_force()`
  * add parameter to call signature
  * extend scope of calculations to apply to `pcrd->params.eGeom == epullgDISTREF`
  * add handling for case `epullROUX`
    * set `dev`, `pcrd->f_scal`, `*V`, and `*dVdl`
    * call `getRouxForce(dev, coord_ind, k)` and `getRouxPotential(dev, coord_ind, k)`
* `void set_pull_coord_reference_value()`
  * adapt to new call signature of `calc_pull_coord_force()`
* `static void do_pull_pot_coord()`
  * adapt to new call signature of `calc_pull_coord_force()`
* `struct pull_t * init_pull()`
  * extend scope of calculations to an additional case `epullgDISTREF`
  
\subsection functinos_added Functions added


* `double getRouxForce(double dev, int coord_ind, double K)`
* `double getRouxPotential(double dev, int coord_ind, double K)`

\subsection other Other


* `md.cpp` uses `step_rel` instead of `step` to determine whether this is a neighbor search step.
* Additional input: environment variable provides a filename that is opened and
  parsed by additional code in the extension.

\section old_implementation Old implementation


- do_md()
    - do_force()
      - do_force_cutsVERLET()
        - pull_potential_wrapper()
          - pull_potential()
            - do_pull_pot_coord()
              - calc_pull_coord_force()
                - getRouxForce()
                - getRouxPotential()

\code
static void calc_pull_coord_force(pull_coord_work_t *pcrd, int coord_ind,
                              double dev, real lambda,
                              real *V, tensor vir, real *dVdl)
{   //...
    switch (pcrd->params.eType)
    {   //...
        case epullROUX:
            dev = pcrd->value;
            /* We want to pass the value of the pull coordinate, not the
             * difference between the pull coordinate and pull-init (which
             * would be dev = pcrd->value - pcrd->value_ref)
             */
            pcrd->f_scal = getRouxForce(dev, coord_ind, k);
            *V          += getRouxPotential(dev, coord_ind, k);
            *dVdl       += *V/k * dkdl;
            break;
    }
    //...
}
\endcode

Example mdp file

\code
;	  gmx grompp -f mdp_files/pull_gromacs5.mdp -c opa_config_0/startup_0.gro -p topologies/topol0.top -o atom_setups/pull_0.tpr -maxwarn 1
pull                     = yes
pull-cylinder-r          = 1.5
pull-constr-tol          = 1e-06
pull-print-com1          = no
pull-print-com2          = no
pull-print-ref-value     = no
pull-print-components    = no
pull-nstxout             = 10
pull-nstfout             = 10
pull_ngroups = 3
pull_ncoords = 1
pull-group1-name = Pull_ref
pull-group2-name = first_0
pull-group3-name = second_0
pull-coord1-type = roux
pull-coord1-geometry = distance-reference
pull-coord1-groups = 2 3 1
pull-coord1-dim = Y Y Y
pull-coord1-origin = 0.0 0.0 0.0
pull-coord1-vec = 0.0 0.0 0.0
pull-coord1-start = no
pull-coord1-init = 0.0
pull-coord1-rate = 0.0
pull-coord1-k = 100.000000
pull-coord1-kB = 0
\endcode

\subsection current_gromacs Current pull code implementation


Relevant types:

* t_pull_coord
* pull_t
* pull_coord_work_t

Dependents

* `t_mdebin` tracks pull as an energy data provider

Hooks:

* `register_external_pull_potential()` only checks consistency. It doesn't actually enable pulling implementations.
  That is hard-coded and must be supported by `init_pull()`, `pull_potential()`, `do_potential()`, `apply_external_pull_coord_force()`
  
Code flow:

`Mdrunner::mdrunner()`

conditionally set `inputrec->pull_work = init_pull(...)`

`finish_pull(inputrec->pull_work)` after call to integrator

`do_md()`

call `pull_print_output(ir->pull_work, step, t)` in simulation loop.

- `gmx_mdrun()`
  - `Mdrunner:mdrunner()`
    - `init_pull()` -> `pull_t`

- `do_md()` accesses the `pull_t` member of the input record to call...
  - `do_force()`
    - `do_force_cutsVERLET()`
      - `pull_potential_wrapper()`
        - `pull_potential()`
          - `do_pull_pot_coord()`
            - `calc_pull_coord_scalar_force_and_potential(pcrd, dev, lambda, V, dVdl);`
            - `calc_pull_coord_vector_force(pcrd);`


\section strategy Strategy

Remove `pull_t` from input record. Available restraints are instantiated separately and available through a manager.

Reimplement `init_pull()` and refactor `pull_t`. Allow client code to provide alternatives.

Refactor `do_pull_pot_coord()` to use member functions of `pull_t` argument.

Where control flow is currently guided by enum values, wrap the function to dynamically dispatch to functor or legacy free functions.

A custom `RestraintPotential` class may (re)implement `calc_pull_coord_scalar_force_and_potential()` or `calc_pull_coord_vector_force();` or some subset of their internals.

Note that `pull_potential()` is the lowest public function in this hierarchy and is passed a pull_t from the inputrec.

\section usage Usage

Extension code provides a new class that implements a `calculate()` function as described at  gmx::RestraintPotential

The new class must then be registered and exposed to external access.


Python client code:

\code{.py}
import numpy
import csv
import roux
import gmx

# Load initial data
with open('rawdata.csv', 'r') as datafile:
    histogram = numpy.asarray(csv.reader(datafile), dtype=float)

puller = roux.RestraintPotential()
puller.histogram = histogram

# Set parameters for experimental protocol
puller.sigma = 1.0
puller.k = 30
puller.nbins = 50

# Load system configuration, topology, simulation and run parameters
system = gmx.system.from_tpr('topol.tpr')
# Note: input record does not indicate the custom code
system.md.add_potential(puller)


def update(session):
    """A potentially slow Python-based callback function for the runner hook.
    
    To allow performance optimizations, do not request data that is not needed.
    """
    # get an mpi4py COMM_WORLD handle
    comm = session.communicator
    size = comm.Get_size()
    rank = comm.Get_rank()
    
    # gather and accumulate statistics
    recvbuf = None
    if rank == 0:
        recvbuf_shape = [size] + puller.data.shape
        recvbuf = numpy.empty(recvbuf_shape, dtype=float)
    comm.Gather(puller.data, recvbuf, root=0)

    # perform analysis
    new_histogram = roux.analyze(recvbuf)

    # broadcast updated histogram
    comm.Bcast(new_histogram, root=0)
    puller.histogram = new_histogram
    if rank == 0:
        with open('rawdata.csv', 'w') as datafile:
            csv.writer(datafile).writerows(puller.histogram)

# Alternatively, subclass a functor or generate a closure
# E.g. in roux module, class RestraintPotential has member function:
    def makeCallback(session):
        # get an mpi4py COMM_WORLD handle
        comm = session.communicator
        
        def callable():
            size = comm.Get_size()
            rank = comm.Get_rank()
            # gather and accumulate statistics
            recvbuf = None
            if rank == 0:
                recvbuf_shape = [size] + self.data.shape
                recvbuf = numpy.empty(recvbuf_shape, dtype=float)
            comm.Gather(puller.data, recvbuf, root=0)
        
            # perform analysis
            new_histogram = roux.analyze(recvbuf)
        
            # broadcast updated histogram
            comm.Bcast(new_histogram, root=0)
            puller.histogram = new_histogram
            if rank == 0:
                with open('rawdata.csv', 'w') as datafile:
                    csv.writer(datafile).writerows(self.histogram)
        
        return callable


# In the simple case, let nranks == num_ensemble_members
with gmx.context.MpiEnsemble(system) as session:

    # The C++ API class for gmx.runner.StopCondition provides
    # `unsigned int nextEvaluationStep()` and `bool evaluate(unsigned int currentStep)`
    # for the runner to call. The basic stop condition set by reading a TPR file looks like
    # `unsigned int nextEvaluationStep() { return nsteps; };` and
    # `bool evaluate(unsigned int currentStep) { return currentStep >= nsteps; };`
    # More sophisticated stop conditions can maintain state and subscribe to additional data.
    # The Python side of the interface is Pythonically flexible.
    #
    # assert isinstance(session.stopCondition, gmx.runner.StopCondition)
    # assert session.stopCondition == system.runner.numSteps
    session.stopCondition = puller.convergenceCondition

    session.run(callback = roux.update, callback_period=10000)
    # or
    #session.run(callback = puller.makeCallback(session), callback_period=10000)
    # The most performance and flexibility may be to allow RestraintPotential derivatives
    # to register or bind to one or more periodic updaters.

\endcode

sample C++ client code

\code{.cpp}

RouxPuller::Builder rouxSetup;
rouxSetup.addHistogram(arraydata);
rouxSetup.setDimensions(...);
std::shared_ptr<RouxPuller> myRoux = rouxSetup.build();

gmx::MdrunnerBuilder simulation;
simulation.addPotential(myRoux);
gmx::Mdrunner session = myContext(simulation);

for (auto i=0; i < n_iter; i++)
{
    session->run();
    myRoux->update();
}
\endcode

Implementing a RestraintPotential derived class.

\code{.cpp}

class RouxPuller : public RestraintPotential
{
public:
    /// Something we can use to set parameters from external code.
    class Builder;
    
    /// Provide interface for runner to set parameters
    struct pull_t* init_pull(FILE *fplog,
                             const pull_params_t *pull_params,
                             const t_inputrec *ir,
                             int nfile,
                             const t_filenm fnm[],
                             const gmx_mtop_t *mtop,
                             t_commrec *cr,
                             const gmx_output_env_t *oenv,
                             real lambda,
                             gmx_bool bOutFile,
                             unsigned long Flags) override;
     
    // // Make it easier to make sense of and override this function to operate on a chosen set of indices
    // real pull_potential(struct pull_t *pull, t_mdatoms *md, t_pbc *pbc,
    //                       t_commrec *cr, double t, real lambda,
    //                       rvec *x, rvec *f, tensor vir, real *dvdlambda) override;
    
    
    
    // No need to override
    // void finish_pull() override;
    
    /// Some sort of external interface
    gmxapi::Status update() override;
    
    // custom member functions.
    // ...
    
private:
    /// custom Potential parameters
    double sigma_;
    double k_;
    /// Histogram data for reference
    std::vector< std::vector<double> > histogram_;
};

struct pull_t* RouxPuller::init_pull(FILE *fplog,
                             const pull_params_t *pull_params,
                             const t_inputrec *ir,
                             int nfile,
                             const t_filenm fnm[],
                             const gmx_mtop_t *mtop,
                             t_commrec *cr,
                             const gmx_output_env_t *oenv,
                             real lambda,
                             gmx_bool bOutFile,
                             unsigned long Flags)
{
    assert(!histogram_.empty());
    this->setCommunicator(cr);
    this->setInputRecord(ir);
    return this;
}

void RouxPuller::calc_pull_coord_scalar_force_and_potential(int coord_ind,
                                                              double t,
                                                              real lambda,
                                                              real *V,
                                                              real *dVdl
                                                              )
{
    pull_coord_work_t *pcrd = this->coord[coord_ind];
    // Calculate distances and set V, dVdl
    dev = pcrd->value;
    /* We want to pass the value of the pull coordinate, not the
     * difference between the pull coordinate and pull-init (which
     * would be dev = pcrd->value - pcrd->value_ref)
     */
    pcrd->f_scal = getRouxForce(dev, coord_ind, k);
    *V          += getRouxPotential(dev, coord_ind, k);
    *dVdl       += *V/k * dkdl;
}

PYBIND11_MODULE(roux_, m) {
    m.doc() = somedocstring;

    pybind11::class_< RouxPuller, RestraintPotential, std::shared_ptr<RouxPuller> > roux_plugin(m, "RestraintPotential", "Implements Roux biasing potential.");
    roux_plugin.def_property("data", &RouxPuller::getData);
    // The rest is inherited from RestraintPotential
    
    pybind11::class_< RouxBuilder > roux_builder(m, "Builder", "Set up potential");
    roux_builder.def(pybind11::init());
    roux_builder.def_property("histogram", &RouxBuilder::addHistogram);
    roux_builder.def_property("dimensions", &RouxBuilder::setDimensions);
    roux_builder.def_property("sigma", &RouxBuilder::setSigma);
    roux_builder.def_property("k", &RouxBuilder::setK);
    roux_builder.def_property("nbins", &RouxBuilder::setNBins);
    roux_builder.def("build", &RouxBuilder::build);
}

\endcode

Python wrapper `roux.py`:

\code{.py}
#roux.py
import gmx
import roux_.RestraintPotential
import roux_.Builder

class RestraintPotential(gmx.RestraintPotential):
    #...
    def bind(self, system):
        builder = roux_.Builder()
        builder.histogram = self.histogram
        builder.sigma = self.sigma
        builder.k = self.k
        builder.nbins = self.nbins
        self.api_object = builder.build()
     #...
     @property
     def histogram(self):
        #...
     @property
     def sigma(self):
        #...
     #...

def analyze(data):
    #...
    return histogram
\endcode
