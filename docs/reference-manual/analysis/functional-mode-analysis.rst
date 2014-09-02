Functional Mode Analysis (FMA)
------------------------------

| :ref:`gmx fma <gmx fma>`, :ref:`gmx anaeig <gmx anaeig>`


Functional Mode Analysis (FMA)
(see refs. \ :ref:`181 <refHubdeGroot2009FMA>` and
:ref:`182 <refKrivobokova2012PLSFMA>`) is a method to elucidate
the relation of molecular structure and function. To this aim,
FMA identifies a collective motion of a user-defined subset of the
system particles that best correlates with changes in an observable
indicative of the molecule's function. For example, one may want to
identify the gating-motion of a channel protein by supplying the
conductance of the channel as observable. The mathematical method
used to implement FMA is partial least-squares (PLS) regression
(see refs. \ :ref:`183 <refDenham1995PLS>` and :ref:`184 <refNg2013PLS>`).

In mathematical terms, PLS-based FMA tries to construct a linear
model that explains variations in the value of the observable
\ :math:`y` as a function of changes in the particle coordinates
\ :math:`x`. More generally, PLS could also treat a set
of multiple observables, but the current implementation sticks to
a single variable \ :math:`y`. The partial in PLS refers to the
fact that this linear model of the relation between \ :math:`x`
and \ :math:`y` is not constructed directly in terms of these sets
of variables, but in terms of two smaller sets of variables
constructed from linear combinations of the original sets. That is,
the model building is preceded by a step of dimensionality reduction.
This step improves the numerical stability of PLS, especially if
the number of \ :math:`y` variables is small and if there are many
similar observations of \ :math:`x`, and it also allows
PLS to work in cases where the number of samples \ :math:`n` is much
smaller than the number of coordinates \ :math:`m`. Generally, in
case of multiple \ :math:`y`, the observables are expressed in terms
of \ :math:`l` orthogonal, linear combinations of the original
variables \ :math:`y_i` such that the variance in the observable
values described by these new variables is maximized, similarly to
principal component analysis. In contrast, for \ :math:`x` PLS
chooses the reduced set of \ :math:`k` variables such that the
covariance with the reduced variables of the set \ :math:`y`
is maximized. That is, the correlation between the reduced sets of
variables is maximized.


PLS analyzes a set of \ :math:`n` data samples. Sticking to the case
of a single \ :math:`y`, each data sample \ :math:`i` associates an
observable value found in the element \ :math:`i` of the
\ :math:`n \times 1` column vector \ :math:`\mathbf{y}` to a set of
coordinate variable values  \ :math:`x` stored in row
\ :math:`i` of the \ :math:`n \times m` matrix \ :math:`\mathbf{X}`.
The values in \ :math:`\mathbf{y}` and \ :math:`\mathbf{X}` are
centered around their arithmetic mean (the average value over the data
set was substracted from each value). Ultimately, PLS aims to find a
vector \ :math:`\mathbf{u}` that minimizes the vector of residuals
\ :math:`\mathbf{\epsilon}` in

.. math:: \mathbf{y} = \mathbf{X}\mathbf{u} + \mathbf{\epsilon}
          :label: eqnFMA1

The vector \ :math:`\mathbf{u}` represents a collective mode of motion
of the system particles that best captures the variation in \ :math:`y`.

The vector \ :math:`\mathbf{u} = \mathbf{W}\,\mathbf{q}` is a linear
combination of orthogonal vectors found in the rows of \ :math:`\mathbf{W}`
(the PLS regression factors or regressors). Each regressor is chosen such
that it is orthogonal to the previous regressors and maximizes the
correlation to the input observable \ :math:`y`. The \ :math:`k \times 1`
vector \ :math:`\mathbf{q}` contains the weights of the regressors
determined by least squares fitting of

.. math:: \mathbf{y} = \mathbf{X}\,\mathbf{W}\,\mathbf{q}
          :label: eqnFMA2

The optimal number of regressors has to be determined empirically, e.g., by
cross-validation. That is, a subset of the available data is not used in
PLS regression but put aside to compare the predictions of the PLS model
for the observable \ :math:`\mathbf{y} = \mathbf{X}\,\mathbf{u}` as function
of the coordinates in this unused validation set to the actual observable
values.


Functional mode analysis can be performed with :ref:`gmx fma <gmx fma>`.
Structures and trajectories can be projected on the obtained functional mode
with :ref:`gmx anaeig <gmx anaeig>`.
