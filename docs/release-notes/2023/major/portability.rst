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


VkFFT support for improved portability and performance on GPUs with OpenCL and SYCL
"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

Support for the VkFFT GPU FFT library was added with two goals:
improved portability across GPU platforms and better performance.
VkFFT can be used with OpenCL and SYCL.
For SYCL builds, VkFFT provides a portable backend for AMD and
NVIDIA GPUs, and it is a better-performing alternative recommended 
at least on AMD with runs without PME decomposition (in non-HeFFTe builds).
For OpenCL builds, VkFFT provides an alternative to ClFFT
with much better performance and broader compiler support.
It is the default on macOS and when building with Visual Studio. 
On other platforms, it can be enabled at build-time using
``-DGMX_GPU_FFT_LIBRARY=VKFFT``.

:issue:`4052`

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
