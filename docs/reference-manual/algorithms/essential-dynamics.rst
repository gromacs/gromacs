Essential Dynamics sampling
---------------------------

The results from Essential Dynamics (see sec. :ref:`covanal`) of a
protein can be used to guide MD simulations. The idea is that from an
initial MD simulation (or from other sources) a definition of the
collective fluctuations with largest amplitude is obtained. The position
along one or more of these collective modes can be constrained in a
(second) MD simulation in a number of ways for several purposes. For
example, the position along a certain mode may be kept fixed to monitor
the average force (free-energy gradient) on that coordinate in that
position. Another application is to enhance sampling efficiency with
respect to usual MD :ref:`65 <refDegroot96a>`, :ref:`66 <refDegroot96b>`.
In this case, the system is
encouraged to sample its available configuration space more
systematically than in a diffusion-like path that proteins usually take.

Another possibility to enhance sampling is flooding. Here a flooding
potential is added to certain (collective) degrees of freedom to expel
the system out of a region of phase space :ref:`67 <refLange2006a>`.

The procedure for essential dynamics sampling or flooding is as follows.
First, the eigenvectors and eigenvalues need to be determined using
covariance analysis (:ref:`gmx covar`) or normal-mode analysis (:ref:`gmx nmeig`).
Then, this information is fed into :ref:`make_edi <gmx make_edi>`, which has many options for
selecting vectors and setting parameters, see ``gmx make_edi -h``. The
generated :ref:`edi` input file is then passed to :ref:`mdrun <gmx mdrun>`.
