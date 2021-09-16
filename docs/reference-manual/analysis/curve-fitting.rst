Curve fitting in |Gromacs|
--------------------------

Sum of exponential functions
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Sometimes it is useful to fit a curve to an analytical function, for
example in the case of autocorrelation functions with noisy tails.
|Gromacs| is not a general purpose curve-fitting tool however and
therefore |Gromacs| only supports a limited number of functions.
:numref:`Table %s <table-fitfn>`  lists the available options with the corresponding
command-line options. The underlying routines for fitting use the
Levenberg-Marquardt algorithm as implemented in the lmfit package \ :ref:`162 <reflmfit>`
(a bare-bones version of which is included in |Gromacs| in which an
option for error-weighted fitting was implemented).

.. |exp|  replace:: :math:`e^{-t/{a_0}}`                                                       
.. |aexp| replace:: :math:`a_1e^{-t/{a_0}}`                                                    
.. |exp2| replace:: :math:`a_1e^{-t/{a_0}}+(1-a_1)e^{-t/{a_2}}`                                
.. |exp5| replace:: :math:`a_1e^{-t/{a_0}}+a_3e^{-t/{a_2}}+a_4`                                
.. |exp7| replace:: :math:`a_1e^{-t/{a_0}}+a_3e^{-t/{a_2}}+a_5e^{-t/{a_4}}+a_6`                
.. |exp9| replace:: :math:`a_1e^{-t/{a_0}}+a_3e^{-t/{a_2}}+a_5e^{-t/{a_4}}+a_7e^{-t/{a_6}}+a_8`
.. |nexp2| replace:: :math:`a_2\ge a_0\ge 0`               
.. |nexp5| replace:: :math:`a_2\ge a_0\ge 0`               
.. |nexp7| replace:: :math:`a_4\ge a_2\ge a_0 \ge0`        
.. |nexp9| replace:: :math:`a_6\ge a_4\ge a_2\ge a_0\ge 0` 

.. _table-fitfn:

.. table:: Overview of fitting functions supported in (most) analysis tools 
    that compute autocorrelation functions. The **Note** column describes 
    properties of the output parameters.
    :align: center
    :widths: auto

    +-------------+------------------------------+---------------------+
    | Command     | Functional form :math:`f(t)` | Note                |
    | line option |                              |                     |
    +=============+==============================+=====================+
    | exp         | |exp|                        |                     |
    +-------------+------------------------------+---------------------+
    | aexp        | |aexp|                       |                     |
    +-------------+------------------------------+---------------------+
    | exp_exp     | |exp2|                       | |nexp2|             |
    +-------------+------------------------------+---------------------+
    | exp5        | |exp5|                       | |nexp5|             |
    +-------------+------------------------------+---------------------+
    | exp7        | |exp7|                       | |nexp7|             |
    +-------------+------------------------------+---------------------+
    | exp9        | |exp9|                       | |nexp9|             |
    +-------------+------------------------------+---------------------+


Error estimation
~~~~~~~~~~~~~~~~

Under the hood |Gromacs| implements some more fitting functions, namely a
function to estimate the error in time-correlated data due to Hess \ :ref:`149 <refHess2002a>`:

.. math:: \varepsilon^2(t) =
          2 \alpha\tau_1\left(1+\frac{\tau_1}{t}\left(e^{-t/\tau_1}-1\right)\right)
                + 2 (1-\alpha)\tau_2\left(1+\frac{\tau_2}{t}\left(e^{-t/\tau_2}-1\right)\right)
          :label: eqntimecorrerror

where :math:`\tau_1` and :math:`\tau_2` are time constants (with
:math:`\tau_2 \ge \tau_1`) and :math:`\alpha` usually is close to 1 (in
the fitting procedure it is enforced that :math:`0\leq\alpha\leq 1`).
This is used in :ref:`gmx analyze <gmx analyze>` for error estimation using

.. math:: \lim_{t\rightarrow\infty}\varepsilon(t) = \sigma\sqrt{\frac{2(\alpha\tau_1+(1-\alpha)\tau_2)}{T}}
          :label: eqnanalyzeerrorest

where :math:`\sigma` is the standard deviation of the data set and
:math:`T` is the total simulation time \ :ref:`149 <refHess2002a>`.

Interphase boundary demarcation
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

In order to determine the position and width of an interface,
Steen-Sæthre *et al.* fitted a density profile to the following function

.. math:: f(x) ~=~ \frac{a_0+a_1}{2} - \frac{a_0-a_1}{2}{\rm
          erf}\left(\frac{x-a_2}{a_3^2}\right)
          :label: eqndesprofilefunc

where :math:`a_0` and :math:`a_1` are densities of different phases,
:math:`x` is the coordinate normal to the interface, :math:`a_2` is the
position of the interface and :math:`a_3` is the width of the
interface \ :ref:`163 <refSteen-Saethre2014a>`. This is implemented
in :ref:`gmx densorder <gmx densorder>`.

Transverse current autocorrelation function
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

In order to establish the transverse current autocorrelation function
(useful for computing viscosity  \ :ref:`164 <refPalmer1994a>`) the following function is
fitted:

.. math:: f(x) ~=~ e^{-\nu}\left({\rm cosh}(\omega\nu)+\frac{{\rm
          sinh}(\omega\nu)}{\omega}\right)
          :label: eqntransverseautocorrfunc

with :math:`\nu = x/(2a_0)` and :math:`\omega = \sqrt{1-a_1}`. This is
implemented in :ref:`gmx tcaf <gmx tcaf>`.

Viscosity estimation from pressure autocorrelation function
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The viscosity is a notoriously difficult property to extract from
simulations \ :ref:`149 <refHess2002a>`, :ref:`165 <refWensink2003a>`. It is *in principle*
possible to determine it by integrating the pressure autocorrelation
function \ :ref:`160 <refPSmith93c>`, however this is often hampered by
the noisy tail of the ACF. A workaround to this is fitting the ACF to
the following function \ :ref:`166 <refGuo2002b>`:

.. math:: f(t)/f(0) = (1-C) {\rm cos}(\omega t) e^{-(t/\tau_f)^{\beta_f}} + C
          e^{-(t/\tau_s)^{\beta_s}}
          :label: eqnviscestpressureautocorr

where :math:`\omega` is the frequency of rapid pressure oscillations
(mainly due to bonded forces in molecular simulations), :math:`\tau_f`
and :math:`\beta_f` are the time constant and exponent of fast
relaxation in a stretched-exponential approximation, :math:`\tau_s` and
:math:`\beta_s` are constants for slow relaxation and :math:`C` is the
pre-factor that determines the weight between fast and slow relaxation.
After a fit, the integral of the function :math:`f(t)` is used to
compute the viscosity:

.. math:: \eta = \frac{V}{k_B T}\int_0^{\infty} f(t) dt
          :label: eqncompviscosity

This equation has been applied to computing the bulk and shear
viscosity using different elements from the pressure tensor \ :ref:`167 <refFanourgakis2012a>`.
