Removed features
^^^^^^^^^^^^^^^^

NVML support removed on NVIDIA GPUs
"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
NVML support (for reporting GPU application clocks  or changing these
for higher throughput) is no longer available. It was only ever supported on
high-end hardware and changing clocks is on recent generations of hardware only
useful when root permissions were available to the user. It may become less useful
as GROMACS evolves, complicated the GROMACS code, and wasn't regularly tested or maintained.
It might return if some of these conditions change.

Quasi-harmonic analysis removed from gmx anaeig
"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
The quasi-harmonic analysis is removed from gmx anaeig because it is
based on the system under study being in an energy minimum. Since this
generally is not the case when studying condensed phase systems the
results would be erratic. The analysis can still be done based on normal
mode analysis results in gmx nmeig.
