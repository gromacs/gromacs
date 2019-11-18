Portability
^^^^^^^^^^^

.. Note to developers!
   Please use """"""" to underline the individual entries for fixed issues in the subfolders,
   otherwise the formatting on the webpage is messed up.
   Also, please use the syntax :issue:`number` to reference issues on redmine, without the
   a space between the colon and number!

Added support for Hygon Dhyana CPU architecture
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
Support for hardware detection and related heuristics has been implemented
for the Hygon Dhyana derived from the first-gen AMD Zen which it shares most
of its architectural details with.

Enabled PME offload support with OpenCL on NVIDIA and Intel GPUs
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
Thanks to portability improvements, the previously disabled PME OpenCL offload
is now enabled also on NVIDIA and Intel GPUs.
