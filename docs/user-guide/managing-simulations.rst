Managing long simulations
=========================

Molecular simulations often extend beyond the lifetime of a single
UNIX command-line process. It is useful to restart the simulation in a
way that would be equivalent to a single run. When :ref:`gmx mdrun` is
halted, it writes a checkpoint file that can restart the simulation
exactly as if there was no interruption. To do this, the checkpoint
retains a full-precision version of the positions and velocities,
along with state information necessary to restart algorithms that
implement coupling to e.g. external thermal reservoirs. A restart can
be attempted using e.g. a :ref:`gro` file but such a restart is not
considered equivalent to one implemented with a checkpoint file.

Such a checkpoint file is also written periodically by :ref:`gmx
mdrun` during the run. The interval is given by the ``-cpt`` flag to
:ref:`gmx mdrun`.


:ref:`gmx mdrun` can be halted in several ways:

* the number of simulation :mdp:`nsteps` can expire
* the user issues a termination signal (e.g. with Ctrl-C on the terminal)
* the job scheduler issues a termination signal when time expires
* when :ref:`gmx mdrun` detects that the length specified with
  ``-maxh`` has elapsed (this option is useful to help cooperate with
  a job scheduler, but can be problematic if jobs can be suspended)
* some kind of catastrophic failure, such as loss of power

To use the checkpoint file for a restart, use a command line such as

::

   gmx mdrun -cpi state

which directs mdrun to use the checkpoint file (which is named
``state.cpt`` by default). You can choose to give the output
checkpoint file a different name with the ``-cpo`` flag, but if so
then you must provide that name as input to ``-cpi`` when you later
use that file.

Appending to output files
-------------------------

By default, :ref:`gmx mdrun` will append to the old output files. If
the previous part ended in a regular way, then the performance data at
the end of the log file will will be removed, some new information
about the run context written, and the simulation will proceed. Otherwise,
mdrun will truncate all the output files back to the time of the last
written checkpoint file, and continue from there, as if the simulation
stopped at that checkpoint in a regular way.

You can choose not to append the output files by using the
``-noappend`` flag, which forces mdrun to write each output to a
separate file, whose name includes a ".partXXXX" string to describe
which simulation part is contained in this file. This numbering starts
from zero and increases monotonically as simulations are restarted,
but does not reflect the number of simulation steps in each part. The
:mdp:`simulation-part` option can be used to set this number manually
in :ref:`gmx grompp`, which can be useful if data has been lost,
e.g. through filesystem failure or user error.

Appending will not work if the files have been modified after mdrun
wrote them, because the checkpoint file maintains a checksum that it
will verify before it writes to them again.

Extending a .tpr file
---------------------

If the simulation described by :ref:`tpr` file has completed and should
be extended, use the :ref:`gmx convert-tpr` tool to extend the run, e.g.

::

   gmx convert-tpr -s previous -extend timetoextendby -o next
   gmx mdrun -s next -cpi state

The time can also be extended using the ``-until`` and ``-nsteps``
options. Note that the original :ref:`mdp` file may have generated
velocities, but that is a one-time operation within :ref:`gmx grompp`
that is never performed again by any other tool.

Changing mdp options for a restart
----------------------------------

If you wish to make changes to your simulations settings other than
length, then you should do so in the :ref:`mdp` file or topology, and
then call

::

   gmx grompp -f possibly-changed -p possibly-changed -cpt old -o new

to instruct :ref:`grompp` to copy the full simulation state in the
checkpoint file into the new :ref:`tpr` file. You should consider your
choices for :mdp:`tinit`, :mdp:`init-step`, :mdp:`nsteps` and
:mdp:`simulation-part`. You should generally not regenerate velocities
with :mdp:`gen-vel`, and generally select :mdp:`continuation` so that
constraints are not re-applied before the first integration step.

Restarts without checkpoint files
---------------------------------

It is possible to perform an exact restart a simulation if you lack a
checkpoint file but have a matching pair of frames in a :ref:`trr` and
:ref:`edr` file written by :ref:`gmx mdrun`. To do this, use

::

   gmx convert-tpr -s old -e matching -t matching -o new

Are continuations exact?
------------------------

If you had a computer with unlimited precision, or if you integrated
the time-discretized equations of motion by hand, exact continuation
would lead to identical results. But since practical computers have
limited precision and MD is chaotic, trajectories will diverge very
rapidly even if one bit is different. Such trajectories will all be
equally valid, but evetunally very different. Continuation using a
checkpoint file, using the same code compiled with the same compiler
and running on the same computer architecture using the same number of
processors would lead to binary identical results. However, by default
the actual work load will be balanced across the hardware according to
the observed execution times. Such trajectories are in principle not
reproducible, and in particular a run that took place in more than one
part will not be identical with an equivalent run in one part - but
neither of them is better in any sense.

Reproducibility
---------------

The following factors affect the reproducibility of a simulation, and thus its output:

* Precision (mixed / double) with double giving "better" reproducibility.
* Number of cores, due to different order in which forces are
  accumulated. For instance (a+b)+c is not necessarily binary
  identical to a+(b+c) in floating-point arithmetic.
* Type of processors. Even within the same processor family there can be slight differences.
* Optimization level when compiling.
* Optimizations at run time: e.g. the FFTW library that is typically
  used for fast Fourier transforms determines at startup which version
  of their algorithms is fastest, and uses that for the remainder of
  the calculations. Since the speed estimate is not deterministic, the
  results may vary from run to run.
* Random numbers used for instance as a seed for generating velocities
  (in GROMACS at the preprocessing stage).
* Uninitialized variables in the code (but there shouldn't be any)
* Dynamic linking to different versions of shared libraries (e.g. for FFTs)
* Dynamic load balancing, since particles are redistributed to
  processors based on elapsed wallclock time, which will lead to
  (a+b)+c != a+(b+c) issues as above
* Number of PME-only ranks (for parallel PME simulations)
* MPI reductions typically do not guarantee the order of the
  operations, and so the absence of associativity for floating-point
  arithmetic means the result of a reduction depends on the order
  actually chosen

The important question is whether it is a problem if simulations are
not completely reproducible. The answer is yes and no. Reproducibility
is a cornerstone of science in general, and hence it is important.
The Central Limit Theorem tells us that in the case of infinitely long
simulation, all observables converge to their equilibrium
values. Molecular simulations in GROMACS adhere to this theorem, and
hence, for instance, the energy of your system will converge to a
finite value, the diffusion constant of your water molecules will
converge to a finite value, and so on. That means all the important
observables, which are the values you would like to get out of your
simulation, are reproducible.

However, it would be useful for debugging if trajectories were
reproducible, too. In order to obtain this it is important to look
over the list above and eliminate the factors that could affect
reproducibility. Further, using

::

   gmx mdrun -reprod

will eliminate all sources of non-reproducibility that it can,
i.e. same executable + same hardware + same shared libraries + same
run input file + same command line parameters will lead to
reproducible results.
