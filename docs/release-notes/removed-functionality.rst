Removed functionality
^^^^^^^^^^^^^^^^^^^^^

NVML support removed on NVIDIA GPUs
"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
NVML support (for reporting GPU application clocks  or changing these
for higher throughput) is no longer available. It was only ever supported on
high-end hardware and changing clocks is on recent generations of hardware only
useful when root permissions were available to the user. It may become less useful
as GROMACS evolves, complicated the GROMACS code, and wasn't regularly tested or maintained.
It might return if some of these conditions change.

Support for CUDA compute capability 2.x removed
"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
The Fermi-era GPUs (cira 2010) are no longer in widespread use, are
not tested in Jenkins, complicated the code, and are no longer
supported.

Contrib directory removed
"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
This code had not been maintained in years, so likely didn't work, and
has been removed. The git history retains its memory if anybody needs
it.

BlueGene support removed
"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
As most of these machines are now retired, and the ports have not been actively
maintained since |Gromacs| 5.1, the support for BlueGene and QPX SIMD has been
removed.

Implicit solvent support removed
"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
Since |Gromacs|-4.6, the SIMD and multi-threading support has been
mostly broken. Since nobody wants to fix it, the feature has been
removed. Old force field files with parameters for such simulations can still be
read, but the values are ignored.

Removed ``gmx mdrun -multi``
"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
The implementation of ``gmx mdrun -multidir`` is more reliable and works with more
features. Nobody was willing to maintain the duplicate functionality.
