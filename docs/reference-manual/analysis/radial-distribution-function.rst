Radial distribution functions
-----------------------------

| :ref:`gmx rdf <gmx rdf>`
| The *radial distribution function* (RDF) or pair correlation function
  :math:`g_{AB}(r)` between particles of type :math:`A` and :math:`B` is
  defined in the following way:

.. math:: \begin{array}{rcl}
          g_{AB}(r)&=&    {\displaystyle \frac{\langle \rho_B(r) \rangle}{\langle\rho_B\rangle_{local}}}         \\
                   &=&    {\displaystyle \frac{1}{\langle\rho_B\rangle_{local}}}{\displaystyle \frac{1}{N_A}}
                          \sum_{i \in A}^{N_A} \sum_{j \in B}^{N_B} 
                          {\displaystyle \frac{\delta( r_{ij} - r )}{4 \pi r^2}}         \\
          \end{array}
          :label: eqnrdfdefine

with :math:`\langle\rho_B(r)\rangle` the particle density of type
:math:`B` at a distance :math:`r` around particles :math:`A`, and
:math:`\langle\rho_B\rangle_{local}` the particle density of type
:math:`B` averaged over all spheres around particles :math:`A` with
radius :math:`r_{max}` (see :numref:`Fig. %s <fig-rdfex>` C).

.. _fig-rdfex:

.. figure:: plots/rdf.*
    :width: 7.00000cm

    Definition of slices in :ref:`gmx rdf <gmx rdf>`: A. :math:`g_{AB}(r)`.
    B. :math:`g_{AB}(r,\theta)`. The slices are colored gray. C.
    Normalization :math:`\langle\rho_B\rangle_{local}`. D. Normalization
    :math:`\langle\rho_B\rangle_{local,\:\theta }`. Normalization volumes
    are colored gray.

Usually the value of :math:`r_{max}` is half of the box length. The
averaging is also performed in time. In practice the analysis program
:ref:`gmx rdf <gmx rdf>` divides the system
into spherical slices (from :math:`r` to :math:`r+dr`, see
:numref:`Fig. %s <fig-rdfex>` A) and makes a histogram in stead of
the :math:`\delta`-function. An example of the RDF of oxygen-oxygen in
SPC water \ :ref::ref:`80 <refBerendsen81>` is given in :numref:`Fig. %s <fig-rdf>`

.. _fig-rdf:

.. figure:: plots/rdfO-O.*
    :width: 8.00000cm

    :math:`g_{OO}(r)` for Oxygen-Oxygen of SPC-water.

With :ref:`gmx rdf <gmx rdf>` it is also possible to calculate an angle
dependent rdf :math:`g_{AB}(r,\theta)`, where the angle :math:`\theta`
is defined with respect to a certain laboratory axis :math:`{\bf e}`,
see :numref:`Fig. %s <fig-rdfex>` B.

.. math:: g_{AB}(r,\theta) = {1 \over \langle\rho_B\rangle_{local,\:\theta }} 
          {1 \over N_A} \sum_{i \in A}^{N_A} \sum_{j \in B}^{N_B} {\delta( r_{ij} - r ) 
          \delta(\theta_{ij} -\theta) \over 2 \pi r^2 sin(\theta)}
          :label: eqnrdfangleaxis1

.. math:: cos(\theta_{ij}) = {{\bf r}_{ij} \cdot {\bf e} \over \|r_{ij}\| \;\| e\| }
          :label: eqnrdfangleaxis2

This :math:`g_{AB}(r,\theta)` is useful for analyzing anisotropic
systems. **Note** that in this case the normalization
:math:`\langle\rho_B\rangle_{local,\:\theta}` is the average density in
all angle slices from :math:`\theta` to :math:`\theta + d\theta` up to
:math:`r_{max}`, so angle dependent, see :numref:`Fig. %s <fig-rdfex>` D.
