Removed features
^^^^^^^^^^^^^^^^

NVML support removed on NVIDIA GPUs
"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
NVML support (for changing GPU application clocks for higher
throughput) is no longer available. It was only implemented in the
high-end hardware, only useful when root permissions were available to
the user, may become less useful as GROMACS evolves, complicated the
GROMACS code, and wasn't regularly tested or maintained. It might
return if some of these conditions change.
