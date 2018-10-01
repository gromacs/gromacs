Simulated Annealing
-------------------

The well known simulated annealing (SA) protocol is supported in
|Gromacs|, and you can even couple multiple groups of atoms separately
with an arbitrary number of reference temperatures that change during
the simulation. The annealing is implemented by simply changing the
current reference temperature for each group in the temperature
coupling, so the actual relaxation and coupling properties depends on
the type of thermostat you use and how hard you are coupling it. Since
we are changing the reference temperature it is important to remember
that the system will NOT instantaneously reach this value - you need to
allow for the inherent relaxation time in the coupling algorithm too. If
you are changing the annealing reference temperature faster than the
temperature relaxation you will probably end up with a crash when the
difference becomes too large.

The annealing protocol is specified as a series of corresponding times
and reference temperatures for each group, and you can also choose
whether you only want a single sequence (after which the temperature
will be coupled to the last reference value), or if the annealing should
be periodic and restart at the first reference point once the sequence
is completed. You can mix and match both types of annealing and
non-annealed groups in your simulation.
