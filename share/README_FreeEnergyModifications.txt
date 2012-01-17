======== GROMACS 4.5 Free Energy Modifications ===================

New with this version:

- coulombic, vdw, bonded, and restraint transformation lambdas all decoupled
  from each other -- any specified pathway can be taken.
- extended ensemble MC, with a number of different move choices
- Supports either fixed weights, or dynamically adjusted weights
- free energy dependent dihedral, angle, and distance restraints

=== BASIC DOCUMENTATION ===

Important options that may have changed behavior in this version: --
defaults listed first, other options in [ ], then parenthetical
comments.

free-energy              = no [yes,mutate,decouple]

lmc-stats            = no [metropolis-transition,barker-transition,wang-landau,gibbs-wang-landau,minvar]  (the
method used to update the weights)

nst-transition-matrix    = 0 [nonnegativeinteger] frequency at which the transition-matrix is output to the log
lmc-mc-move          = gibbs [metropolis,barker,gibbs,metropolized-gibbs] (the method used to for lambda MC moves)

lmc-seed                 = -1 [for lambda mc transformations, -1 means it's taken from the process number]
mc-temperature          = [positive real] [If omitted, set to the same as the ref_t for the 0th system, if there is one]

lmc-gibbsdelta          = -1 [any integer] (the interval [lambda-1,lambda+1] for moves for gibbs sampling on the lambdas.  -1,
the default, means the entire interval [0,1] will be sampled)

initial-wl-delta        = 1.0 [any positive real] (the initial delta factor for Wang-Landau, in kT)
wl-scale                = 0.8 [real between 0 and 1] (the scaling factor for Wang-Landau)
wl-ratio                = 0.8 [real between 0 and 1] ratio of lowest to highest occupancies for Wang-Landau being seen as flat)    
symmetrized-transition-matrix=no [yes/no]  Whether to compute a symmetrized version of the transition matrix, by averaging the transpose with the matrix.

lmc-forced-nstart       = 0 [any positive integer] (the number of equilibration steps at each lambda when 
			    'forcing through', essentially peforming slow growth, to initialize weights

mininum-var-min        = 100 [any positive integer] must have this many samples in each state before reweighting to 
                       	 get the number of states per state that minimizes variace.  Preferable to make it 
                         something decent so that it the initial weights don't get swallowed up in noise, 
			 as it is using an asymptotic theory.

weight-c-range          = 0 when using minvar, uses the C that is closest to the G*(ln n0/n1).  
			  0 means it defaults to C=0. Otherwise, it computes the free energies 
			  along +/- c-range, with kT spacing.  Should eventually be used for 
                          BAR as well.

lmc-weights-equil       = no [weight-equil-wl-delta,weight-equil-number-all-lambda,weight-equil-number-steps,weight-equil-number-samples,weight-equil-count-ratio]
                        The condition to set when we stop updating the weights and start to equilibrate.

lmc-weights-equil       = no [wl-delta,number-all-lambda,number-steps,number-samples,count-ratio]
weight-equil-wl-delta           =  [positive real] stop when wl-delta drops below this value, if lmc-weights-equil  = wl-delta
weight-equil-number-all-lambda  =  [positive integer] stop when we've got this number of samples at all lambda, if lmc-weights-equil  = number-all-lambda
weight-equil-number-steps       =  [positive integer] stop when we've done this number of steps, if lmc-weights-equil  = number-steps
weight-equil-number-samples     =  [positive integer] stop when we've done this number of samples, if lmc-weights-equil  = number-samples
weight-equil-count-ratio        =  [positive real] stop when the ratio of min to max value is greater than this, if lmc-weights-equil  = count-ratio

fep-lambdas         = 0.0 0.2 0.3 0.4 0.5 0.6 0.7 0.8 1.0 (array of lambdas)
mass-lambdas        = 0.0 0.2 0.3 0.4 0.5 0.6 0.7 0.8 1.0 (array of lambdas)
coul-lambdas        = 0.0 0.2 0.3 0.4 0.5 0.6 0.7 0.8 1.0 (array of lambdas)
vdw-lambdas         = 0.0 0.2 0.3 0.4 0.5 0.6 0.7 0.8 1.0 (array of lambdas)
bonded-lambdas      = 0.0 0.2 0.3 0.4 0.5 0.6 0.7 0.8 1.0 (array of lambdas)
restraint-lambdas   = 0.0 0.2 0.3 0.4 0.5 0.6 0.7 0.8 1.0 (array of lambdas)
init-lambda-weights = 0.0 0.2 0.3 0.4 0.5 0.6 0.7 0.8 1.0 (initial weights for each lambda)

sc-alpha            = 0.5  [any positive real] (soft core alpha - now, only operates on vdw)
sc-coul             = no [yes/no] -  Controls whether coulomb is also softocore 

nstfep              = 20 [any integer multiple of nlist] (the
frequency at which the other energies are determined, and the lambdas are updated.
nstdgdl             = 200 [any integer multiple of nstfep] (the rate at which
the dE terms are output to the dgdl file, per step)
dhdl-print-energy   = Whether to print the energies as well as the energy differences (helps for some analysis codes --
                      for example, required when doing umbrella sampling with different temperatures)
Notes on non-self-explanatory terms:

Weight methods:
*barker-transition -- computes <alpha>_forward/<alpha>_reverse, with
alpha being the Fermi function.  Like Bennett, but without self
consistent iteration.  More efficient than metropolis transition
*metropolis-transition, same as above but with the Metropolis function
(min{1,exp(dE}) instead of the Fermi function
*wang-landau- implements standard wang-landau
*gibbs-wang-landau -- uses p(k)delta for nonlocal updating of weights, instead of delta
*mbar.  Highly experimental, recommend against using.
*minvar - optimizes the number of samples at each states with weights that lower the variance, as described
by Escobedo (ref #???)

weight-equilibration criteria options:

*weight-equil-wl-delta - stop when wl-delta drops below this value set in 'weight-equil-wl-delta'.
*weight-equil-number-all-lambda - stop when all the lambdas are equal or greater than 'weight-equil-number-all-lambda'
*weight-equil-number-steps - stop when we've done the number of steps set in 'weight-equil-number-steps'
*weight-equil-number-samples - stop when we've done the number of samples set in 'weight-equil-number-samples'
*weight-equil-count-ratio - stop when the ratio of top to bottom value is greater than the value set in 'weight-equil-count-ratio'

wl-ratio -- For Wang-Landau, scaling is done when everything is
within mc-ratio of the mean.

Right now, nstfep must be an integer multiple of nslist -- this can
possibly be fixed later.  For now, since we update all the other
energies every nstfep steps, it probably would be too slow if done
more frequently.  Previously, I though if it was too fast, you could
get an unlucky run to an unstable configuration, if it does not
equilibrate the momenta at the new level.  I think I fixed a bug that
was causing this problem.  It might be possible to go down to 10,
though at that point, the cost starts to be large (takes about 20%
longer with nstfep = 10 vs. 20 in my sample system).

This does bring up a point that perhaps for strict detailed balance, the
momenta need to be resampled from a Boltzmann distribution at each
transition.  This would slow down dynamics a lot, though -- is it
necessary for complete detailed balance?  (Answer seems to be that if 
the velocity distribution function is equal in both states, as it is
if the temperature is the same for both, then we're OK).

fep-lambdas is the default lambda array.  All the other four lambda
arrays (except the serial tempering temperatures) use the values in
this array if they are not specified.  If not specificed, it is
assumed to be zero throughout.

For example, if you wanted to first to change coul, then vdw, while
changing bonded at the same time as vdw, but the restraints throughout
the first two thirds of the simulation, then you'd do something like:

init-lambda-state      = 0
coul-lambdas           = 0.0 0.2 0.5 1.0 1.0 1.0 1.0 1.0 1.0 1.0
vdw-lambdas            = 0.0 0.0 0.0 0.0 0.4 0.5 0.6 0.7 0.8 1.0
bonded-lambdas         = 0.0 0.0 0.0 0.0 0.4 0.5 0.6 0.7 0.8 1.0
restraint-lambdas      = 0.0 0.0 0.1 0.2 0.3 0.5 0.7 1.0 1.0 1.0

This is equivalent to:

init-lambda-state      = 0 
fep-lambdas            = 0.0 0.0 0.0 0.0 0.4 0.5 0.6 0.7 0.8 1.0
coul-lambdas           = 0.0 0.2 0.5 1.0 1.0 1.0 1.0 1.0 1.0 1.0
restraint-lambdas      = 0.0 0.0 0.1 0.2 0.3 0.5 0.7 1.0 1.0 1.0

The fep-lambda array, in this case, is being used as the default to
fill in the bonded and vdw lambda arrays.  Usually, it's best to fill
in all array explicitly, just to make sure things are properly
assigned.

If you want to turn on only restraints going from A to B, then it would be:

init-lambda-state      = 1
restraint-lambdas      = 0.0 0.1 0.2 0.4 0.6 1.0

The current weights are written out to the logfile every nstlog steps,
in units of kT.

You can load in weights using init-lambda-weights.  Use 0.00 for the first
(zero because it's the reference) entries, then the rest of the weights,
like:

init-lambda-weights =  0.00  0.54  1.45 3.45 6.78

However, most of the update methods don't do a good job with using
given weights if there is no data added in previously.  So, it's best
to use this for entering the weights for a subsequent equilibrium
simulation.  I'll try to figure out if there is a better way.

Make sure you make the spacing in lambda sufficiently close.  It needs
to be close in order to have good overlap.  I'd keep the free energy
between states less that 5 kT or so, perhaps more like 2-3 kT.

 ============ Coding ====================

The expanded ensemble changes (other than keeping track of global
variables) are almost mostly contained in src/mdlib/sim_util.c, with
additional changes in readir.c and mdebin.c

Things that would be nice to address:

  -nstdgdl not necessarily a multiple of nslist (especially with nstlist = -1)
 
Test to do:

 - throw a number of systems on to find an bugs.

Let me know if there is anything else I should be adding. . .

There's a very large number combinations of different parameters and
approaches that could be adjusted in order to study what works best.

======== Implementation Notes====

I'm still playing around trying to figure out the best combination of
methods to use.  If the weights are poor, equilibration to a flat
histogram in weights can be very slow.  Right now, I think that the
Gibbs Wang-Landau method might be the fastest to equilibrate for less
overlap.  In the limit of large numbers of samples, it probably
doesn't matter very much.  


======Notes on topologies for restraints=======
******This section is old, and may have some issues *****
For distance and dihedral restraints, you need entries like this in
the mdp -- for angle restraints, you don't.  This is a strange gromacs
convention.

Dihedral restraint entries look like this -- note the A->B state, the
A state and the B state (last two commented).

[ dihedral_restraints ]
;  i    j    k    l type label power phiA dphiA  kfacA  phiB  dphiB   kfacB
;mixed state
1410 1393 1391 2610    1     0     2    38    0    0.00     38     0    41.84
1393 1391 2610 2604    1     1     2   111    0    0.00    111     0    41.84
1391 2610 2604 2606    1     2     2   -39    0    0.00    -39     0    41.84
; zero state
;1410 1393 1391 2610    1     0     2   38    0    0.00
;1393 1391 2610 2604    1     1     2  111    0    0.00
;1391 2610 2604 2606    1     2     2  -39    0    0.00
; one state
;1410 1393 1391 2610    1     0     2   38    0   41.84
;1393 1391 2610 2604    1     1     2  111    0   41.84
;1391 2610 2604 2606    1     2     2  -39    0   41.84

Free energy dependent angle restraints look like this:

[ angle_restraints ]
;  i    j    k    l type theta0A fcA mult  theta0B    fcB
;mixed state  -- need to have MULT listed twice!!!
1393 1391 2610 1391    1   62      0.0  1     62    41.84    1
1391 2610 2604 2610    1   79      0.0  1     79    41.84    1
; zero state
;1393 1391 2610 1391   1   62      0.0  1
;1391 2610 2604 2610   1   79      0.0  1
; one state
;1393 1391 2610 1391   1   62    41.84  1
;1391 2610 2604 2610   1   79    41.84  1

Note particularly that even though the multiplicity can't actually be
changed in free energy perturbation, you need to have the multiplicity
listed in both the A and B states.  This is a strange effect of using
proper dihedral codes to do the angle restraints -- I tried to find a
way to fix it, but gave up after a few hours -- it works using this
topology framework.

Free energy dependent harmonic distance restraints are implemented as
follows -- R2 is a ridiculously large number in order to obtain purely
harmonic restraints -- otherwise, it becomes linear outside of R2.

[ distance_restraints ]
;  i    j type label typeprime    r0A   r1A   r2A    fcA  r0B   r1B
r2B      fcB
;mixed state
1391 2610     1     0         1   0.61  0.61   10.0     0.0  0.61  0.61 10.0  4184.0 ;zero state
;1391 2610    1     0         1   0.61  0.61   10.0     0.0 ;one state
;1391 2610    1     0         1   0.61  0.61    0.81 4184.0


