The accelerated weight histogram method (AWH) {#page_awh}
=============================================

Accelerating sampling with AWH
==============================

AWH calculates the free energy along an order parameter of the system.
Free energy barriers are overcome by adaptively tuning a bias potential along
the order parameter such that the biased distribution along the parameter
converges toward a chosen target distribution.
The fundamental equation governing the tuning is: log(target) = bias - free energy, where
the bias and free energy are initially unknown. Typically the target distribution is simply
chosen uniform, such that the bias completely flattens the free energy landscape.


Design of the AWH module
========================

The module implements AWH for the case when the order parameter corresponds to a reaction coordinate,
here referred to as coordinate for short, i.e. a function of the system configuration.
The bias is coupled to the system by a bias potential: either in the form of an harmonic ("umbrella") potential
Monte-Carlo (MC) "jumping" around the current coordinate value, or as a smooth convolution of the umbrellas.

The AWH module is organizes as follows:
The Awh class is the interface between the outside and inside of the module.
The Awh class contains one or more BiasCoupledToSystem objects.
The BiasCoupledToSystem class takes care of the reaction coordinate input
and force output for the single Bias object it containts.
The Bias class is a container and wrapper for a object BiasState + helpers.
All computation takes place in the BiasState object and its sub-classes.
The Bias class also contains a BiasWriter object that takes care of i/o.

Use of AWH in mdrun
===================

The basic use of Awh in mdrun consists of 2 method calls:
Call the constructor Awh() after the pull module has been initialized.
Call applyBiasForcesAndUpdateBias() at every MD step after the pull
potential calculation function has been called.

In grompp the pull potential provider should be registered using
registerAwhWithPull() so grompp can check for unregistered potentials.

The main tasks of AWH are:
- calculate and set the bias force given the current coordinate value.
- after accumulating a number of coordinate samples, update the free energy estimate and the bias.

AWH currently relies on the pull code for the first task. Pull provides AWH with updated coordinate values
and distributes the bias force that AWH calculates to the atoms making up the coordinate. This
also means that there are some order dependencies where pull functions need to be called before AWH
functions (see below).

The implementation is quite general. There can be multiple independent AWH biases coupled to the system
simultaneously. This makes sense if the system is made up of several fairly independent parts,
like monomers in a protein. Each bias acts on exactly one, possibly multidimensional, coordinate.
Each coordinate dimension maps to exactly one pull coordinate. Thus, an n-dimensional
biased coordinate is defined by a set of n pull coordinates. Periodicity is taken care of for coordinate
dimensions that require it (dihedral angles). For increased parallelism, there is the option of
having multiple communicating simulations sharing all samples. All simulations would then share a single
bias and free energy estimate. Alternatively, one may partition the sampling domain into smaller
subdomains with some overlap and have multiple independent simulations sample each subdomain.

Note that internally the AWH module keep tracks of free energies in units
of the thermal energy kT. This is because we mostly deal with free energies
in the form of -log(probability) and using any other unit would be bug prone.
All energy type variables are explicitly documented to be in units of kT.
Also the checkpoint and energy file data is in units of kT. The analysis
tool will by default convert energies to kJ/mol, but there is also
a kT option.

