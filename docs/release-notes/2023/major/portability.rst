Portability
^^^^^^^^^^^

Full support for RISC-V
"""""""""""""""""""""""
We now provide full support for RISC-V, including hardware
cycle counters for efficient load balancing.

Initial support for Apple silicon GPUs
""""""""""""""""""""""""""""""""""""""
We now recognize Apple-designed GPUs as a supported architecture
in the OpenCL backend.

PME GPU offload on macOS
""""""""""""""""""""""""
Until now, PME calculations could not be offloaded to the GPU on
macOS. They required the clFFT library, which silently crashed Apple's 
OpenCL drivers at runtime. To overcome this incompatibility, we
replaced the clFFT backend with VkFFT on macOS.

Increase of required versions
"""""""""""""""""""""""""""""
* GCC required version is now 9.
* oneMKL required version is now 2021.3.
