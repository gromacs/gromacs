========================
Using the Python package
========================

After installing GROMACS, sourcing the "GMXRC" (see GROMACS docs), and installing
the gmxapi Python package (see :doc:`install`), import the package in a Python
script or interactive interpreter. This documentation assumes a convenient alias
of ``gmx`` to refer to the ``gmxapi`` Python package.

::

    import gmxapi as gmx

For full documentation of the Python-level interface and API, use the ``pydoc``
command line tool or the :py:func:`help` interactive Python function, or refer to
the :doc:`pythonreference`.

Any Python *exception* raised by gmxapi
should be descended from (and catchable as) :class:`gmxapi.exceptions.Error`.
Additional status messages can be acquired through the :ref:`gmxapi logging`
facility.
Unfortunately, some errors occurring in the GROMACS library are not yet
recoverable at the Python level, and much of the standard GROMACS terminal
output is not yet accessible through Python.
If you find a particularly problematic scenario, please file a GROMACS bug report.

During installation, the *gmxapi* Python package becomes tied to a specific
GROMACS installation.
If you would like to access multiple GROMACS installations
from Python, build and install *gmxapi* in separate
:ref:`virtual environments <gmxapi venv>`.

In some cases *gmxapi* still needs help finding infrastructure from the
GROMACS installation.
For instance, :py:func:`gmxapi.commandline_operation` is not a pure API utility,
but a wrapper for command line tools.
Make sure that the command line tools you intend to use are discoverable in
your :envvar:`PATH`, such as by "source"ing your :file:`GMXRC` before launching
a *gmxapi* script.

.. todo:: Get relevant GROMACS paths in Python environment.

    :py:class:`gmxapi.commandline_operation` relies on the environment :envvar:`PATH`
    to locate executables, including the :command:`gmx` wrapper binary.
    Relates to `#2961 <https://redmine.gromacs.org/issues/2961>`__.

.. _parallelism:

Notes on parallelism and MPI
============================

When launching a *gmxapi* script in an MPI environment,
such as with :command:`mpiexec` or :command:`srun`,
you must help *gmxapi* detect the MPI environment by ensuring that :py:mod:`mpi4py`
is loaded.
Refer to :ref:`mpi_requirements` for more on installing :py:mod:`mpi4py`.

Assuming you use :command:`mpiexec` to launch MPI jobs in your environment,
run a *gmxapi* script on two ranks with something like the following.
Note that it can be helpful to provide :command:`mpiexec` with the full path to
the intended Python interpreter since new process environments are being created.

::

    mpiexec -n 2 `which python` -m mpi4py myscript.py

*gmxapi* 0.1 has limited parallelism, but future versions will include seamless
acceleration as integration improves with the GROMACS library and computing
environment runtime resources.
Currently, *gmxapi* and the GROMACS library do not have an effective way to
share an MPI environment.
Therefore, if you intend to run more than one simulation at a time, in parallel,
in a *gmxapi* script, you should build GROMACS with *thread-MPI* instead of a
standard MPI library.
I.e. configure GROMACS with the CMake flag ``-DGMX_THREAD_MPI=ON``.
Then, launch your *gmxapi* script with one MPI rank per node, and *gmxapi* will
assign each (non-MPI) simulation to its own node, while keeping the full MPI
environment available for use via :py:mod:`mpi4py`.

Running simple simulations
==========================

Once the ``gmxapi`` package is installed, running simulations is easy with
:py:func:`gmxapi.read_tpr`.

::

    import gmxapi as gmx
    simulation_input = gmx.read_tpr(tpr_filename)
    md = gmx.mdrun(simulation_input)

Note that this sets up the work you want to perform, but does not immediately
trigger execution. You can explicitly trigger execution with::

    md.run()

or you can let gmxapi automatically launch work in response to the data you
request.

The :py:func:`gmxapi.mdrun` operation produces a simulation trajectory output.
You can use ``md.output.trajectory`` as input to other operations,
or you can get the output directly by calling ``md.output.trajectory.result()``.
If the simulation has not been run yet when ``result()`` is called,
the simulation will be run before the function returns.

Running ensemble simulations
============================

To run a batch of simulations, just pass an array of inputs.::

    md = gmx.read_tpr([tpr_filename1, tpr_filename2, ...])
    md.run()

Make sure to launch the script in an MPI environment with a sufficient number
of ranks to allow one rank per simulation.

For *gmxapi* 0.1, we recommend configuring the GROMACS build with
``GMX_THREAD_MPI=ON`` and allowing one rank per node in order to allow each
simulation ensemble member to run on a separate node.

.. seealso:: :ref:`parallelism`

.. _commandline:

Accessing command line tools
============================

In *gmxapi* 0.1, most GROMACS tools are not yet exposed as *gmxapi* Python operations.
:class:`gmxapi.commandline_operation` provides a way to convert a :command:`gmx`
(or other) command line tool into an operation that can be used in a *gmxapi*
script.

In order to establish data dependencies, input and output files need to be
indicated with the ``input_files`` and ``output_files`` parameters.
``input_files`` and ``output_files`` key word arguments are dictionaries
consisting of files keyed by command line flags.

For example, you might create a :command:`gmx solvate` operation as::

    solvate = gmx.commandline_operation('gmx',
                                        arguments=['solvate', '-box', '5', '5', '5'],
                                        input_files={'-cs': structurefile},
                                        output_files={'-p': topfile,
                                                      '-o': structurefile,
                                                      }

To check the status or error output of a command line operation, refer to the
``returncode`` and ``erroroutput`` outputs.
To access the results from the output file arguments, use the command line flags
as keys in the ``file`` dictionary output.

Example::

    structurefile = solvate.output.file['-o'].result()
    if solvate.output.returncode.result() != 0:
        print(solvate.output.erroroutput.result())

Preparing simulations
=====================

Continuing the previous example, the output of ``solvate`` may be used as the
input for ``grompp``::

    grompp = gmx.commandline_operation('gmx', 'grompp',
                                       input_files={
                                           '-f': mdpfile,
                                           '-p': solvate.output.file['-p'],
                                           '-c': solvate.output.file['-o'],
                                           '-po': mdout_mdp,
                                       },
                                       output_files={'-o': tprfile})

Then, ``grompp.output.file['-o']`` can be used as the input for :py:func:`gmxapi.read_tpr`.

Simulation input can be modified with the :py:func:`gmxapi.modify_input` operation
before being passed to :py:func:`gmxapi.mdrun`.
For *gmxapi* 0.1, a subset of MDP parameters may be overridden using the
dictionary passed with the ``parameters`` key word argument.

Example::

    simulation_input = gmx.read_tpr(grompp.output.file['-o'])
    modified_input = gmx.modify_input(input=simulation_input, parameters={'nsteps': 1000})
    md = gmx.mdrun(input=modified_input)
    md.run()

Using arbitrary Python functions
================================

Generally, a function in the *gmxapi* package returns an object that references
a node in a work graph,
representing an operation that will be run when the graph executes.
The object has an ``output`` attribute providing access to data Futures that
can be provided as inputs to other operations before computation has actually
been performed.

You can also provide native Python data as input to operations,
or you can operate on native results retrieved from a Future's ``result()``
method.
However, it is trivial to convert most Python functions into *gmxapi* compatible
operations with :py:func:`gmxapi.function_wrapper`.
All function inputs and outputs must have a name and type.
Additionally, functions should be stateless and importable
(e.g. via Python ``from some.module import myfunction``)
for future compatibility.

Simple functions can just use :py:func:`return` to publish their output,
as long as they are defined with a return value type annotation.
Functions with multiple outputs can accept an ``output`` key word argument and
assign values to named attributes on the received argument.

Examples::

    from gmxapi import function_wrapper

    @function_wrapper(output={'data': float})
    def add_float(a: float, b: float) -> float:
        return a + b


    @function_wrapper(output={'data': bool})
    def less_than(lhs: float, rhs: float, output=None):
        output.data = lhs < rhs

.. seealso::

    For more on Python type hinting with function annotations,
    check out :pep:`3107`.

Subgraphs
=========

Basic *gmxapi* work consists of a flow of data from operation outputs to
operation inputs, forming a directed acyclic graph (DAG).
In many cases, it can be useful to repeat execution of a subgraph with updated
inputs.
You may want a data reference that is not tied to the immutable result
of a single node in the work graph, but which instead refers to the most recent
result of a repeated operation.

One or more operations can be staged in a :py:class:`gmxapi.operation.Subgraph`,
a sort of meta-operation factory that can store input binding behavior so that
instances can be created without providing input arguments.

The subgraph *variables* serve as input, output, and mutable internal data
references which can be updated by operations in the subgraph.
Variables also allow state to be propagated between iterations when a subgraph
is used in a *while* loop.

Use :py:func:`gmxapi.subgraph` to create a new empty subgraph.
The ``variables`` argument declares data handles that define the state of the
subgraph when it is run.
To initialize input to the subgraph, give each variable a name and a value.

To populate a subgraph, enter a SubgraphContext by using a :py:func:`with` statement.
Operations created in the *with* block will be captued by the SubgraphContext.
Define the subgraph outputs by assigning operation outputs to subgraph variables
within the *with* block.

After exiting the *with* block, the subgraph may be used to create operation
instances or may be executed repeatedly in a *while* loop.

.. note::

    The object returned by :py:func:`gmxapi.subgraph` is atypical of *gmxapi*
    operations, and has some special behaviors. When used as a Python
    `context manager <https://docs.python.org/3/reference/datamodel.html#context-managers>`__,
    it enters a "builder" state that changes the behavior of its attribute
    variables and of operaton instantiation. After exiting the :py:func:`with`
    block, the subgraph variables are no longer assignable, and operation
    references obtained within the block are no longer valid.

Looping
=======

An operation can be executed an arbitrary number of times with a
:py:func:`gmxapi.while_loop` by providing a factory function as the
*operation* argument.
When the loop operation is run, the *operation* is instantiated and run repeatedly
until *condition* evaluates ``True``.

:py:func:`gmxapi.while_loop` does not provide a direct way to provide *operation*
arguments. Use a *subgraph* to define the data flow for iterative operations.

When a *condition* is a subgraph variable, the variable is evaluated in the
running subgraph instance at the beginning of an iteration.

Example::

    subgraph = gmx.subgraph(variables={'float_with_default': 1.0, 'bool_data': True})
    with subgraph:
        # Define the update for float_with_default to come from an add_float operation.
        subgraph.float_with_default = add_float(subgraph.float_with_default, 1.).output.data
        subgraph.bool_data = less_than(lhs=subgraph.float_with_default, rhs=6.).output.data
    operation_instance = subgraph()
    operation_instance.run()
    assert operation_instance.values['float_with_default'] == 2.

    loop = gmx.while_loop(operation=subgraph, condition=subgraph.bool_data)
    handle = loop()
    assert handle.output.float_with_default.result() == 6

.. _gmxapi logging:

Logging
=======

*gmxapi* uses the Python :py:mod:`logging` module to provide hierarchical
logging, organized by submodule.
You can access the logger at ``gmxapi.logger`` or, after importing *gmxapi*,
through the Python logging framework::

    import gmxapi as gmx
    import logging

    # Get the root gmxapi logger.
    gmx_logger = logging.getLogger('gmxapi')
    # Set a low default logging level
    gmx_logger.setLevel(logging.WARNING)
    # Make some tools very verbose
    #  by descending the hierarchy
    gmx_logger.getChild('commandline').setLevel(logging.DEBUG)
    #  or by direct reference
    logging.getLogger('gmxapi.mdrun').setLevel(logging.DEBUG)

You may prefer to adjust the log format or manipulate the log handlers.
For example, tag the log output with MPI rank::

    try:
        from mpi4py import MPI
        rank_number = MPI.COMM_WORLD.Get_rank()
    except ImportError:
        rank_number = 0
        rank_tag = ''
        MPI = None
    else:
        rank_tag = 'rank{}:'.format(rank_number)

    formatter = logging.Formatter(rank_tag + '%(name)s:%(levelname)s: %(message)s')

    # For additional console logging, create and attach a stream handler.
    ch = logging.StreamHandler()
    ch.setFormatter(formatter)
    logging.getLogger().addHandler(ch)

For more information, refer to the Python `logging documentation <https://docs.python.org/3/library/logging.html>`__.

More
====

Refer to the :doc:`pythonreference` for complete and granular documentation.

For more information on writing or using pluggable simulation extension code,
refer to https://redmine.gromacs.org/issues/3133.
(For gmxapi 0.0.7 and GROMACS 2019, see https://github.com/kassonlab/sample_restraint)

.. todo:: :issue:`3133`: Replace these links as resources for pluggable extension code become available.
