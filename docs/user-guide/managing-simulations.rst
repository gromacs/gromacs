.. _managing long simulations:

Managing long simulations
=========================

Molecular simulations often extend beyond the lifetime of a single
UNIX command-line process. It is useful to be able to stop and
restart the simulation in a
way that is equivalent to a single run. When :ref:`gmx mdrun` is
halted, it writes a checkpoint file that can restart the simulation
exactly as if there was no interruption. To do this, the checkpoint
retains a full-precision version of the positions and velocities,
along with state information necessary to restart algorithms e.g.
that implement coupling to external thermal reservoirs. A restart can
be attempted using e.g. a :ref:`gro` file with velocities, but since
the :ref:`gro` file has significantly less precision, and none of
the coupling algorithms will have their state carried over, such
a restart is less continuous than a normal MD step.

Such a checkpoint file is also written periodically by :ref:`gmx
mdrun` during the run. The interval is given by the ``-cpt`` flag to
:ref:`gmx mdrun`. When :ref:`gmx mdrun` attempts to write each
successive checkpoint file, it first renames the old file with the
suffix ``_prev``, so that even if something goes wrong while writing
the new checkpoint file, only recent progress can be lost.

:ref:`gmx mdrun` can be halted in several ways:

* the number of simulation :mdp:`nsteps` can expire
* the user issues a termination signal (e.g. with Ctrl-C on the terminal)
* the job scheduler issues a termination signal when time expires
* when :ref:`gmx mdrun` detects that the length specified with
  ``-maxh`` has elapsed (this option is useful to help cooperate with
  a job scheduler, but can be problematic if jobs can be suspended)
* some kind of catastrophic failure, such as loss of power, or a
  disk filling up, or a network failing

To use the checkpoint file for a restart, use a command line such as

::

   gmx mdrun -cpi state

which directs mdrun to use the checkpoint file (which is named
``state.cpt`` by default). You can choose to give the output
checkpoint file a different name with the ``-cpo`` flag, but if so
then you must provide that name as input to ``-cpi`` when you later
use that file. You can
query the contents of checkpoint files with :ref:`gmx check` and
:ref:`gmx dump`.

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

Appending will not work if any output files have been modified or
removed after mdrun wrote them, because the checkpoint file maintains
a checksum of each file that it will verify before it writes to them
again. In such cases, you must either restore the file, name them
as the checkpoint file expects, or continue with ``-noappend``. If
your original run used ``-deffnm``, and you want appending, then
your continuations must also use ``-deffnm``.

Backing up your files
---------------------

You should arrange to back up your simulation files frequently. Network
file systems on clusters can be configured in more or less conservative
ways, and this can lead :ref:`gmx mdrun` to be told that a checkpoint
file has been written to disk when actually it is still in memory
somewhere and vulnerable to a power failure or disk that fills or
fails in the meantime. The UNIX tool rsync can be a useful way to
periodically copy your simulation output to a remote storage location,
which works safely even while the simulation is underway. Keeping a copy
of the final checkpoint file from each part of a job submitted to a
cluster can be useful if a file system is unreliable.

Extending a .tpr file
---------------------

If the simulation described by :ref:`tpr` file has completed and should
be extended, use the :ref:`gmx convert-tpr` tool to extend the run, e.g.

::

   gmx convert-tpr -s previous.tpr -extend timetoextendby -o next.tpr
   gmx mdrun -s next.tpr -cpi state.cpt

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

   gmx grompp -f possibly-changed.mdp -p possibly-changed.top -c original.gro -t state.cpt -o new.tpr
   gmx mdrun -s new.tpr -cpi state.cpt

to instruct :ref:`gmx grompp` to copy the full-precision coordinates
and velocities in the checkpoint file into the new :ref:`tpr` file.
You should consider your choices for :mdp:`tinit`, :mdp:`init-step`,
:mdp:`nsteps` and :mdp:`simulation-part`. You should generally not
regenerate velocities with :mdp:`gen-vel`, and generally select
:mdp:`continuation` so that constraints are not re-applied before
the first integration step.

Restarts without checkpoint files
---------------------------------

It used to be possible to continue simulations without the checkpoint
files. As this approach could be unreliable or lead to
unphysical results, only restarts from checkpoints are permitted now.

Are continuations exact?
------------------------

If you had a computer with unlimited precision, or if you integrated
the time-discretized equations of motion by hand, exact continuation
would lead to identical results. But since practical computers have
limited precision and MD is chaotic, trajectories will diverge very
rapidly even if one bit is different. Such trajectories will all be
equally valid, but eventually very different. Continuation using a
checkpoint file, using the same code compiled with the same compiler
and running on the same computer architecture using the same number of
processors without GPUs (see next section) would lead to binary
identical results. However,
by default the actual work load will be balanced across the hardware
according to the observed execution times. Such trajectories are
in principle not reproducible, and in particular a run that took
place in more than one part will not be identical with an equivalent
run in one part - but neither of them is better in any sense.

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
  (in |Gromacs| at the preprocessing stage).
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
* On GPUs, the reduction of e.g. non-bonded forces has a non-deterministic
  summation order, so any fast implementation is non-reproducible by
  design.

The important question is whether it is a problem if simulations are
not completely reproducible. The answer is yes and no. Reproducibility
is a cornerstone of science in general, and hence it is important.
The `Central Limit Theorem <https://en.wikipedia.org/wiki/Central_limit_theorem>`_
tells us that in the case of infinitely long
simulations, all observables converge to their equilibrium
values. Molecular simulations in |Gromacs| adhere to this theorem, and
hence, for instance, the energy of your system will converge to a
finite value, the diffusion constant of your water molecules will
converge to a finite value, and so on. That means all the important
observables, which are the values you would like to get out of your
simulation, are reproducible. Each individual trajectory is not
reproducible, however.

However, there are a few cases where it would be useful if
trajectories were reproducible, too. These include developers doing
debugging, and searching for a rare event in a trajectory when, if
it occurs, you want to have manually saved your checkpoint file so
you can restart the simulation under different conditions, e.g.
writing output much more frequently.

In order to obtain this reproducible trajectory, it is important
to look over the list above and eliminate the factors that could
affect it. Further, using

::

   gmx mdrun -reprod

will eliminate all sources of non-reproducibility that it can,
i.e. same executable + same hardware + same shared libraries + same
run input file + same command line parameters will lead to
reproducible results.
