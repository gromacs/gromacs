/* ************************************************************************
 * Copyright 2013-2015 Advanced Micro Devices, Inc.
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 * http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 * ************************************************************************/

/*! @file mainpage.h

This file contains all documentation, no code, in the form of comment text.  It provides
chapter 1 of the documentation that is produced with doxygen. Chapter 1 includes the title page, installation instructions, and prose on the nature of FFT and their use in our library.

@mainpage OpenCL Fast Fourier Transforms (FFTs)

The clFFT library is an OpenCL library implementation of discrete Fast Fourier Transforms. The library:
@li provides a fast and accurate platform for calculating discrete FFTs.
@li works on CPU or GPU backends.
@li supports in-place or out-of-place transforms.
@li supports 1D, 2D, and 3D transforms with a batch size that can be greater than or equal to 1.
@li supports planar (real and complex components are stored in separate arrays) and interleaved (real and complex components are stored as a pair in the same array) formats.
@li supports lengths that are any combination of powers of 2, 3, 5, and 7.
@li supports single and double precision floating point formats.


@section IntroFFT Introduction to clFFT

The FFT is an implementation of the Discrete Fourier Transform (DFT) that makes use of symmetries in the FFT definition to reduce the mathematical intensity required from O(\f$N^2\f$) to O(\f$ N \log N\f$) when the sequence length, *N*, is the product of small prime factors.  Currently, there is no standard API for FFT routines. Hardware vendors usually provide a set of high-performance FFTs optimized for their systems: no two vendors employ the same interfaces for their FFT routines. clFFT provides a set of FFT routines that are optimized for AMD graphics processors, and that are also functional across CPU and other compute devices.

@subsection SupportRadix Supported radices
clFFT supports transform sizes that are powers of 2, 3, 5, and 7. This means that the vector lengths can be  a combination of powers of two, three, five, and seven; examples include \f$2^7, 2^1*3^1, 3^2*5^4, 2^2*3^3*5^5\f$, up to the limit that the device can support.

@subsection SizeLimit Transform size limits
Currently, there is an upper bound on the transform size that the library can support for certain transforms. This limit is \f$2^{24}\f$ for real 1D single precision and \f$2^{22}\f$ for real 1D double precision.

@subsection EnumDim Dimensionality
clFFT currently supports FFTs of (up to) three dimensions, given by the enum @ref clfftDim. This enum
is a required parameter of @ref clfftCreateDefaultPlan() to create an initial plan, where a plan is the collection of (almost) all the parameters needed to specify an FFT computation. For more information about clFFT plans, see the section \ref clFFTPlans.
Depending on the dimensionality that the client requests, clFFT uses the following formulations to compute the DFT:

@li For a 1D complex DFT
\f[
{\tilde{x}}_j = {{1}\over{scale}}\sum_{k=0}^{n-1}x_k\exp\left({\pm i}{{2\pi jk}\over{n}}\right)\hbox{ for } j=0,1,\ldots,n-1
\f]
where, \f$x_k\f$ are the complex data to be transformed, \f$\tilde{x}_j\f$ are the transformed data, and the sign \f$\pm\f$ determines the direction of the transform: \f$-\f$ for forward and \f$+\f$ for backward. Note that you must provide the scaling factor.  By default, the scale is set to 1 for forward transforms, and \f${{1}\over{N}}\f$ for backward transforms, where *N* is the size of transform.

@li For a 2D complex DFT
\f[
{\tilde{x}}_{jk} = {{1}\over{scale}}\sum_{q=0}^{m-1}\sum_{r=0}^{n-1}x_{rq}\exp\left({\pm i} {{2\pi jr}\over{n}}\right)\exp\left({\pm i}{{2\pi kq}\over{m}}\right)
\f]
for \f$j=0,1,\ldots,n-1\hbox{ and } k=0,1,\ldots,m-1\f$, where, \f$x_{rq}\f$ are the complex data to be transformed, \f$\tilde{x}_{jk}\f$ are the transformed data, and the sign \f$\pm\f$ determines the direction of the transform.  By default, the scale is set to 1 for forwards transforms and \f${{1}\over{MN}}\f$ for backwards transforms, where *M* and *N* are the 2D size of the transform.

@li For a 3D complex DFT
\f[
\tilde{x}_{jkl} = {{1}\over{scale}}\sum_{s=0}^{p-1}\sum_{q=0}^{m-1}\sum_{r=0}^{n-1}
x_{rqs}\exp\left({\pm i} {{2\pi jr}\over{n}}\right)\exp\left({\pm i}{{2\pi kq}\over{m}}\right)\exp\left({\pm i}{{2\pi ls}\over{p}}\right)
\f]
for \f$j=0,1,\ldots,n-1\hbox{ and } k=0,1,\ldots,m-1\hbox{ and } l=0,1,\ldots,p-1\f$, where \f$x_{rqs}\f$ are the complex data to be transformed, \f$\tilde{x}_{jkl}\f$ are the transformed data, and the sign \f$\pm\f$ determines the direction of the transform. By default, the scale is set to 1 for forward transforms and \f${{1}\over{MNP}}\f$ for backward transforms, where *M*, *N*, and *P* are the 3D size of the transform.

@subsection InitLibrary Setup and Teardown of clFFT
clFFT is initialized by the API @ref clfftSetup(), which must be called before any other API of
clFFT. This allows the library to create resources needed to manage the plans that you create and
destroy. This API also takes a structure @ref clfftInitSetupData() that is initialized by the
client to control the behavior of the library.

After you use the library, the @ref clfftTeardown() method must be called. This function instructs clFFT to release all resources allocated internally, and resets acquired references to any OpenCL objects.

@subsection ThreadSafety Thread safety
The clFFT API is designed to be thread-safe. It is safe to create plans from multiple threads and to
destroy those plans in separate threads. Multiple threads can call @ref clfftEnqueueTransform() to place work in a command queue at the same time. clFFT does not provide a single-threaded version of the library. The overhead of the synchronization mechanisms inside a clFFT thread-safe is expected to be minor.

Currently, you must manage the multi-device operation. You can create OpenCL contexts that are
associated with multiple devices, but clFFT only uses a single device from that context to transform
the data. You can manage a multi-device operation by creating multiple contexts, in which each
context contains a different device; you are responsible for scheduling and partitioning the work
across multiple devices and contexts.

@subsection MajorFormat Row major formats
clFFT expects all multi-dimensional input passed to it to be in row-major format. This is compatible
with C-based languages. However, clFFT is very flexible in the organization of the input and output data, and it accepts input data by letting you specify a stride for each dimension. This feature can be used to process data in column major arrays and other non-contiguous data formats. See @ref clfftSetPlanInStride() and 
@ref clfftSetPlanOutStride().

@subsection Object OpenCL object creation
Your application must allocate and manage OpenCL objects, such as contexts,  *cl_mem* buffers and command queues. All the clFFT interfaces that interact with OpenCL objects take those objects as references through the API.
Specifically, the plan creation function @ref clfftCreateDefaultPlan() takes an OpenCL context as a parameter reference, increments the reference count on that object, and keeps the object alive until the corresponding plan is destroyed by the call @ref clfftDestroyPlan().

@subsection FlushQueue Flushing of command queues
The clFFT API operates asynchronously; with the exception of thread safety locking with multiple 
threads, all APIs return immediately. Specifically, the @ref clfftEnqueueTransform() API does not
explicitly flush the command queues that are passed by reference to it. It pushes the transform work onto the command queues and returns the modified queues to the client. The client is free to issue its own blocking logic by using OpenCL synchronization mechanisms or push further work onto the queue to continue processing.

@subsection EnvVariables Environment variables
The clFFT library looks for definition of two environment variables: CLFFT_CACHE_PATH and CLFFT_REQUEST_LIB_NOMEMALLOC. 
If the variable CLFFT_CACHE_PATH is defined, the library caches OpenCL binaries. This enables a subsequent run of the application with the same type of transforms to avoid the expensive compilation step. Instead, the stored binaries are loaded and executed. The CLFFT_CACHE_PATH must point to a folder location where the library can store binaries. 
The other variable CLFFT_REQUEST_LIB_NOMEMALLOC when defined, requests the library to do all computations in-place and avoid allocating extra device memory whenever possible. This feature is experimental and currently works only for certain types of transforms  when the library decomposes the input into square matrices or rectangular matrices with dimensions in the ratio 1:2. Currently, it works for 1D complex transforms of size of powers of 2. 

@section clFFTPlans clFFT plans
A plan is the collection of (almost) all the parameters needed to specify an FFT computation.
A clFFT plan includes the following parameters:
<ul>
<li> The OpenCL context that executes the transform
<li> Dimension of the transform (1D, 2D or 3D)
<li> Length or extent of data in each dimension
<li> Number of datasets that are transformed
<li> Precision of the data
<li> Scaling factor to the transformed data
<li> In-place or Out-of-place transform
<li> Format of the input data - interleaved, planar or real
<li> Format of the output data - interleaved, planar or real
</ul>

The clFFT plan does not include the following parameters:
<ul>
<li> The OpenCL handles to the input and output data buffers.
<li> The OpenCL handle to a temporary scratch buffer (if needed).
<li> Direction of execution of the transform (forward or reverse transform).
</ul>
These parameters are specified when the plan is executed.

@subsection Default Default plan values

When a new plan is created by calling @ref clfftCreateDefaultPlan(), its parameters are initialized as
follows:

<ul>
<li> Dimensions: as provided by the caller
<li> Lengths: as provided by the caller
<li> Batch size: 1
<li> Precision: *CLFFT_SINGLE*
<li> Scaling factors:
    <ol>
    <li> for the forward transform, the default is 1.0, or no scale factor is applied
    <li> for the reverse transform, the default is 1.0 / P, where P is the product of the FFT lengths
    </ol>
<li> Location: *CLFFT_INPLACE*
<li> Input layout: *CLFFT_COMPLEX_INTERLEAVED*
<li> Input strides: the strides of a multidimensional array of the lengths specified, where the data is
compactly stored using the row-major convention
<li> Output layout: *CLFFT_COMPLEX_INTERLEAVED*
<li> Output strides: same as input strides
</ul>

Writing client programs that depend on these initial values is <b> not </b> recommended.

@subsection EnumLayout Supported memory layouts
There are two main types of Discrete Fourier Transform (DFT) in clFFT:
<ol>
<li> Transformation of complex data - clFFT supports the following two layouts to store complex numbers:
<ul>
  <li> Planar format - where the real and imaginary components are kept in separate arrays: \n
	   Buffer1: **RRRRR**  \n
	   Buffer2: **IIIII**
  <li> Interleaved format - where the real and imaginary components are stored as contiguous pairs:  \n
	   Buffer1: **RIRIRIRIRIRI**
</ul>
<li> Transformation of real to complex data and vice versa - clFFT provides enums to define these formats.
For transforms involving real data, there are two possibilities:
<ul>
<li> Real data being subject to forward FFT transform that results in complex data.
<li> Complex data being subject to backward FFT transform that results in
real data. See the section \ref RealFFT.
</ul>
</ol>

@subsubsection DistanceStridesandPitches Strides and Distances
For one-dimensional data, if clStrides[0] = strideX = 1, successive elements in the first dimension are stored contiguously in memory. If strideX is an integral value greater than 1, gaps in memory exist between each element of the vectors.
For multi-dimensional data, if clStrides[1] = strideY = LenX for 2 dimensional data and clStrides[2] = strideZ = LenX*LenY for 3 dimensional data, no gaps exist in memory between each element, and all vectors are stored tightly packed in memory. Here, LenX, LenY, and LenZ denote the transform lengths clLengths[0], clLengths[1], and clLengths[2], respectively, which are used to set up the plan.

By specifying non-default strides, it is possible to process either row-major or column-major arrays. Data can be extracted from arrays of structures. Almost any regular data storage pattern can be accommodated.

Distance is the amount of memory that exists between corresponding elements in an FFT primitive in a batch. Distance is measured in units of the FFT primitive; complex data measures in complex units, and real data measures in real units. Stride between tightly packed elements is 1 in either case. Typically, one can measure the distance between any two elements in a batch primitive, be it 1D, 2D, or 3D data. For tightly packed data, the distance between FFT primitives is the size of the FFT primitive, such that dist=LenX for 1D data, dist=LenX*LenY for 2D data, and dist=LenX*LenY*LenZ for 3D data. It is possible to set the distance of a plan to be less than the size of the FFT vector; most often 1 for this case. When computing a batch of 1D FFT vectors, if distance == 1, and strideX == length(vector), a transposed output is produced for a batch of 1D vectors. You must verify that the distance and strides are valid (not intersecting); if not valid, undefined results may occur.
A simple example would be to perform a 1D length 4096 on each row of an array of 1024 rows x 4096 columns of values stored in a column-major array, such as a FORTRAN program might provide. (This would be equivalent to a C or C++ program that has an array of 4096 rows x 1024 columns stored in a row-major manner, on which you want to perform a 1-D length 4096 transform on each column.) In this case, specify the strides as [1024, 1].

A more complex example would be to compute a 2D FFT for each 64 x 64 subtile of the grid that has an input buffer with a raster grid of 1024 x 1024 monochrome pixel values. Specifying strides allows you to treat each horizontal band of 1024 x 64 pixels as an array of 16 64 x 64 matrixes, and process an entire band with a single call @ref clfftEnqueueTransform(). (Specifying strides is not quite flexible enough to transform the entire grid of this example with a single kernel execution.) It is possible to create a Plan to compute arrays of 64 x 64 2D FFTs, then specify three strides: [1, 1024, 64]. The first stride, 1, indicates that the rows of each matrix are stored consecutively; the second stride, 1024, gives the distance between rows, and the third stride, 64, defines the distance between two matrices. Then call @ref clfftEnqueueTransform() 16 times – once for each horizontal band of pixels.

@subsection EnumPrecision Supported precisions in clFFT
Both *CLFFT_SINGLE* and *CLFFT_DOUBLE* precisions are supported by the library for all supported radices. For both these enums the math functions of the host computer are used to produce the sine and cosine tables that are used by the OpenCL kernel.
Both *CLFFT_SINGLE_FAST* and *CLFFT_DOUBLE_FAST* generate faster kernels with reduced accuracy, but are disabled in the current build.
See @ref clfftPrecision, @ref clfftSetPlanPrecision(), and @ref clfftGetPlanPrecision().

@subsection FftDirection clfftDirection
For complex transforms, the direction of the transform is not baked into the plan; the same plan can be used to specify both forward and backward transforms. To specify the direction, @ref clfftDirection is passed as a parameter into @ref clfftEnqueueTransform(). 
For real transforms, the  input and output layouts of the plan determine the direction.

@subsection EnumResultLocation In-place and out-of-place transforms
The clFFT API supports both in-place and out-of-place transforms. With in-place transforms, only the input buffers are provided to the @ref clfftEnqueueTransform() API, and the resulting data is written in the same buffer, overwriting the input data. With out-of-place transforms, distinct output buffers are provided to the @ref clfftEnqueueTransform() API, and the input data is preserved. 
In-place transforms require that the *cl_mem* objects created by the client application, have both read and write permissions. This is given in the nature of the in-place algorithm. Out-of-place transforms require that the destination buffers have read and write permissions, but input buffers can still be created with read-only permissions. This is a clFFT requirement because internally the algorithms may go back and forth between the destination buffers and internally allocated temp buffers. For out-of-place transforms, clFFT never writes back to input buffers.

@subsection clFFTEff Batches
The efficiency of clFFT is improved by utilizing transforms in batches. Sending as much data as possible in a single transform call leverages the parallel compute capabilities of OpenCL devices (and GPU devices in particular), and minimizes the penalty of transfer overhead. It is best to think of an OpenCL device as a high-throughput, high-latency device. Using a networking analogy as an example, this approach is similar to having a massively high-bandwidth pipe with very high ping response times. If the client is ready to send data to the device for compute, it should be sent in as few API calls as possible and this can be done by batching. clFFT plans have a parameter @ref clfftSetPlanBatchSize() to describe the number of transforms being batched, and another parameter @ref clfftSetPlanDistance() to describe how those batches are laid out and spaced in memory. 1D, 2D, or 3D transforms can be batched.

@section Outline  Using clFFT in a client application

To perform FFT calculations using clFFT, the client program must perform the following tasks:
<ul>
	<li> Initialize the library by calling @ref clfftSetup(). </li>
	<li> For each distinct type of FFT needed:
	<ol>
   		 <li> Create an FFT Plan object. 
           			To create an FFT Plan object, do either of the following. 
		<ul>
			<li>Call the factory function @ref clfftCreateDefaultPlan() and specify the value 
			       of the most fundamental parameters, such as plHandle, context, dim, and 
 			       clLengths, while other parameters assume default values.  The 
			       OpenCL context must be provided when the plan is created; it cannot be 
			       changed. </li> <br /> Or
		 	<li>Call @ref clfftCopyPlan(). </li>
		</ul>
            			Note: In both the cases, the function returns an opaque handle to the plan 
			object.
    <li> Complete the specification of all the Plan parameters by calling various parameter-setting 
            functions that have the prefix *clfftSet*. </li>
	  <li> Optionally, "bake" or finalize the plan by calling @ref clfftBakePlan() function. This signals 
	          the library the end of the specification phase, and causes it to generate and compile the 
                         exact OpenCL kernels that perform the specified FFT on the provided OpenCL device.
	          At this point, all performance-enhancing optimizations are applied, possibly including 
                        executing benchmark kernels on the OpenCL device context to maximize the runtime 
                        performance. <br /> <br />
Although the previous step is optional, it is recommended to use it so that you can
have control on when to do this work. Usually, this time consuming step is done when the application is initialized. If you do not call @ref clfftBakePlan(), this work is done during the first call of @ref clfftEnqueueTransform().
		</li>
	</ol>

	<li> Execute the OpenCL FFT kernels as many times as needed. </li>
	<ol>
		<li>  Call @ref clfftEnqueueTransform(). At this point, specify whether you want to 
		         execute a forward or reverse transform; also, provide the OpenCL *cl_mem* 
		         handles for the input buffer(s), output buffer(s)--unless you want the transformed 
		         data to overwrite the input buffers, and (optionally) scratch buffer.
		         @ref clfftEnqueueTransform() performs one or more calls to the OpenCL function 
		          clEnqueueNDRangeKernel. Like clEnqueueNDRangeKernel, @ref 
                                       clfftEnqueueTransform() is a non-blocking call. The commands to execute the FFT 
                                       compute kernel(s) are added to the OpenCL context queue to be executed 
		          asynchronously.
		         An OpenCL event handle is returned to the caller. If multiple NDRangeKernel 
		         operations are queued, the final event handle is returned.
		</li>
		<li>  Add any application OpenCL tasks to the OpenCL context queue. For example, if the    		        next step in the application process is to apply a filter to the transformed data, the 
		        application calls clEnqueueNDRangeKernel, and specifies its output buffer(s) as the  
		        input to the filter kernel, and provides the transform's event handle to ensure 
		        proper synchronization. </li>
		<li>  If the application accessed the transformed data directly, it calls one of the OpenCL 
		        functions for synchronizing the host computer execution with the OpenCL device 
		        (for example: clFinish()). </li>
	</ol>
	<li> Terminate the library by calling @ref clfftTeardown().
</ul>

@section RealFFT  FFTs of real data

When real data is subject to DFT transformation, the resulting complex output data follows a special property. About half of the output is redundant because they are complex conjugates of the other half. This is called the Hermitian redundancy. So, for space and performance considerations, it is only necessary to store the non-redundant part of the data. Most FFT libraries use this property to offer specific storage layouts for FFTs involving real data. clFFT provides three enumerated types to deal with real data FFTs:
<ul>
	<li> *CLFFT_REAL*
	<li> *CLFFT_HERMITIAN_INTERLEAVED*
	<li> *CLFFT_HERMITIAN_PLANAR*
</ul>

The *CLFFT_REAL* enum specifies that the data is purely real. This can be used to feed real input or get back real output. The *CLFFT_HERMITIAN_INTERLEAVED* and *CLFFT_HERMITIAN_PLANAR* enums are similar to the corresponding full complex enums in the way they store real and imaginary components, but store only about half of the complex output. Client applications can do just a forward transform and analyze the output or they can process the output and do a backward transform to get back real data. This is illustrated in the following figure.

@image html realfft_fwdinv.jpg "Forward and Backward Transform Processes"

Let us consider a 1D real FFT of length N. The full output looks as shown in following figure.

@image html realfft_1dlen.jpg "1D Real FFT of Length N"

Here, C* denotes the complex conjugate. Since the values at indices greater than N/2 can be deduced from the first half of the array, clFFT stores data only up to the index N/2. This means that the output contains only 1 + N/2 complex elements, where the division N/2 is rounded down. Examples for even
and odd lengths are given below.

Example for N = 8 is shown in following figure.

@image html realfft_ex_n8.jpg "Example for N = 8"

Example for N = 7 is shown in following figure.

@image html realfft_ex_n7.jpg "Example for N = 7"


For length 8, only (1 + 8/2) = 5 of the output complex numbers are stored, with the index ranging from 0 through 4. Similarly for length 7, only (1 + 7/2) = 4 of the output complex numbers are stored, with the index ranging from 0 through 3.
For 2D and 3D FFTs, the FFT length along the least dimension is used to compute the (1 + N/2) value. This is because the FFT along the least dimension is computed first and is logically a real-to-hermitian transform. The FFTs along other dimensions are computed afterwards; they are simply 'complex-to-complex' transforms. For example, assuming clLengths[2] is used to set up a 2D real FFT, let N1 = clLengths[1], and N0 = clLengths[0]. The output FFT has N1*(1 + N0/2) complex elements. Similarly, for a 3D FFT with clLengths[3] and N2 = clLengths[2], N1 = clLengths[1], and N0 = clLengths[0], the output has N2*N1*(1 + N0/2) complex elements.

@subsection RealModes Supported Modes

Out-of-place transforms:

<ul>
	<li> *CLFFT_REAL* to *CLFFT_HERMITIAN_INTERLEAVED*
	<li> *CLFFT_REAL* to *CLFFT_HERMITIAN_PLANAR*
	<li> *CLFFT_HERMITIAN_INTERLEAVED* to *CLFFT_REAL*
	<li> *CLFFT_ CLFFT_HERMITIAN_PLANAR* to *CLFFT_REAL*
</ul>

In-place transforms:

<ul>
	<li> *CLFFT_REAL* to *CLFFT_HERMITIAN_INTERLEAVED*
	<li> *CLFFT_HERMITIAN_INTERLEAVED* to *CLFFT_REAL*
</ul>

@subsection ExplicitStrides Setting strides

The library currently <b> requires you to explicitly set input and output strides for real transforms.</b> See the following examples to understand what values to use for input and output strides under different scenarios. These examples show typical usages, but you can allocate the buffers and layout data according to your need.

@subsection RealExamples Examples

The following pages provide figures and examples to explain in detail the real FFT features of this library.

@image html realfft_expl_01.jpg "1D FFT - Real to Hermitian"
@image html realfft_expl_02.jpg "1D FFT - Real to Hermitian, Example 1"
@image html realfft_expl_03.jpg "1D FFT - Real to Hermitian, Example 2"
@image html realfft_expl_04.jpg "1D FFT - Real to Hermitian, Example 3"
@image html realfft_expl_05.jpg "1D FFT - Hermitian to Real"
@image html realfft_expl_06.jpg "1D FFT - Hermitian to Real, Example"
@image html realfft_expl_07.jpg "2D FFT - Real to Hermitian In Place"
@image html realfft_expl_08.jpg "2D FFT - Real to Hermitian, Example"

@section Callbacks  clFFT Callbacks

The callback feature of clFFT has the ability to invoke user provided OpenCL™ inline functions to pre-process or post-process data from within the FFT kernel. The inline OpenCL callback function is passed as a string to the library. It is then incorporated into the generated FFT kernel. This eliminates the need for an additional kernel launch to carry out the pre/post processing tasks, thus improving overall performance.
There are 2 types of callback; Pre-callback and Post-callback. Pre-callback invokes user callback function to perform custom  pre-processing of the input data before FFT is executed. Post-callback invokes user callback function to perform custom post-processing of the output data after FFT is executed.

@subsection CallbackWorkflow Callback Workflow

The workflow of FFT execution using callback feature of clFFT is as follows:

<ol>
	<li> Create the clFFT Plan and initialize the standard clFFT parameters.
	<li> Use @ref clfftSetPlanCallback() API to register the callback function with library
		@code
		clfftStatus clfftSetPlanCallback(clfftPlanHandle plHandle,
											const char* funcName,
											const char* funcString,
											int localMemSize,
											clfftCallbackType callbackType,
											void *userdata,
											int numUserdataBuffers)
		@endcode
		The library uses the arguments passed to this API, including callback function string, to stitch the callback code 	into the generated FFT kernel. The arguments for clfftSetPlanCallback are
		<ul>
			<li> clFFT plan handle
			<li> Name of the callback function
			<li> Callback function as character array; the character array can include custom datatype declaration if any custom type is used by callback function
			<li> Size of local memory requested by callback, if any, in bytes
			<li> Type of callback: ‘PRECALLBACK’ or ‘POSTCALLBACK’; this is an enumerator
			<li> Supplementary user data, if any, used by callback function
			<li> Number of user data buffers; the library currently supports only one user data buffer per callback registration
		</ul>
		Multiple callback registration calls to the same type of callback result in overwriting the previously registered callback function
	<li> Invoke Bake Plan step
	Library inserts the callback code into the main FFT kernel during bake plan and compiles it. If there are any compilation errors caused by syntax or incompatible callback function prototype, the library reports failure.
	<li> Enqueue clFFT transform
</ol>

The caller is responsible to provide a callback function that matches the function prototype based on the type of
callback(pre/post), type of transform(real/complex) and whether LDS is used. The bake plan step checks the function prototype.

@subsection CallbackFunctionPrototype Callback Function Prototypes

clFFT expects the callback function to be of a specific prototype depending on the type of callback(pre/post), type of transform(real/complex) and whether LDS is used. These are as follows:

@subsubsection PrecallbackProtyotype Pre-callback Prototypes

 FFT Type                               | Function Prototype
----------------------------------------| -------------
C2C/C2R – Interleaved Single Precision  | Without LDS <br />float2  <precallback_func>  (  __global void *input, uint inoffset, __global void *userdata) <br /> With LDS <br />float2  <precallback_func>  (  __global void *input, uint inoffset, __global void *userdata, __local void *localmem)
C2C/C2R – Interleaved Double Precision  | Without LDS <br />double2  <precallback_func>  (  __global void *input, uint inoffset, __global void *userdata) <br /> With LDS <br />double2  <precallback_func>  (  __global void *input, uint inoffset, __global void *userdata, __local void *localmem)
C2C – Planar Single Precision			| Without LDS <br />float2  <precallback_func>  (  __global void *inputRe, __global void *inputIm, uint inoffset, __global void *userdata)<br /> With LDS <br />float2  <precallback_func>  (  __global void *inputRe, __global void *inputIm, int inoffset, __global void *userdata, __local void *localmem)
C2C – Planar Double Precision			| Without LDS <br />double2  <precallback_func>  (  __global void *inputRe, __global void *inputIm, uint inoffset, __global void *userdata)<br /> With LDS <br />double2  <precallback_func>  (  __global void *inputRe, __global void *inputIm, uint inoffset, __global void *userdata, __local void *localmem)
R2C Single Precision					| Without LDS <br />float  <precallback_func>   (  __global void *input, uint inoffset, __global void *userdata)<br /> With LDS <br />float  <precallback_func>   (  __global void *input, uint inoffset, __global void *userdata, __local void *localmem)
R2C Double Precision					| Without LDS <br />double  <precallback_func>   (  __global void *input, uint inoffset, __global void *userdata)<br /> With LDS <br />double  <precallback_func>   (  __global void *input, uint inoffset, __global void *userdata, __local void *localmem)


Parameters
<ul>
	<li> \c input : The base pointer of the input buffer for R2C and Interleaved C2C/C2R transforms
	<li> \c inputRe : The base pointer of the “Real” input buffer for Planar C2C transforms
	<li> \c inputIm : The base pointer of the “Imaginary” part input buffer for Planar C2C transforms
	<li> \c inoffset : Index of the current element  of the input buffer from the start
	<li> \c userdata : Buffer containing optional caller specified data. The userdata pointer is useful
	for passing any supplementary data to the callback function. For example, buffer having convolution
	filter data or any scalar value. The userdata can be of any custom data type/structure, in which case,
	you have to declare the custom data type and include it along with the callback function string. </li>
	<li> \c localmem : Pointer to local memory. This memory is allocated by library based on the size you specify
	and is subjected to local memory availability. </li>
</ul>

For Planar C2C, the return type of callback is a vector (float2/double2) that contains the   Real
and Imaginary elements as computed in the callback.

@subsubsection PostcallbackProtyotype Post-callback Prototypes

 FFT Type                               | Function Prototype
----------------------------------------| ------------------
C2C/R2C – Interleaved Single Precision  | Without LDS <br />void  <postcallback_func> ( __global void *output, uint outoffset, __global void *userdata, float2 fftoutput) <br /> With LDS <br />void  <postcallback_func> ( __global void *output, uint outoffset, __global void *userdata, float2 fftoutput, __local void *localmem)
C2C/R2C – Interleaved Double Precision  | Without LDS <br />void  <postcallback_func> ( __global void *output, uint outoffset, __global void *userdata, double2 fftoutput) <br /> With LDS <br />void  <postcallback_func> ( __global void *output, uint outoffset, __global void *userdata, double2 fftoutput, __local void *localmem)
C2C/R2C – Planar Single Precision		| Without LDS <br />void  <postcallback_func> ( __global void *outputRe, __global void *outputIm, uint outoffset, __global void *userdata, float fftoutputRe, float fftoutputIm) <br /> With LDS <br />void  <postcallback_func> ( __global void *outputRe, __global void *outputIm, uint outoffset, __global void *userdata, float fftoutputRe, float fftoutputIm, __local void *localmem)
C2C/R2C – Planar Double Precision		| Without LDS <br />void  <postcallback_func> ( __global void *outputRe, __global void *outputIm, uint outoffset, __global void *userdata, double fftoutputRe, double fftoutputIm) <br /> With LDS <br />void  <postcallback_func> ( __global void *outputRe, __global void *outputIm, uint outoffset, __global void *userdata, double fftoutputRe, double fftoutputIm, __local void *localmem)
C2R Single Precision					| Without LDS <br />void  <postcallback_func> ( __global void *output, uint outoffset, __global void *userdata, float fftoutput) <br /> With LDS <br />void  <postcallback_func> ( __global void *output, uint outoffset, __global void *userdata, float fftoutput, __local void *localmem)
C2R Double Precision					| Without LDS <br />void  <postcallback_func> ( __global void *output, uint outoffset, __global void *userdata, double fftoutput) <br /> With LDS <br />void  <postcallback_func> ( __global void *output, uint outoffset, __global void *userdata, double fftoutput, __local void *localmem)


Parameters
<ul>
	<li> \c output  : The base pointer of the output buffer for C2R and Interleaved R2C/C2C transforms
	<li> \c outputRe : The base pointer of the “Real” output buffer for Planar R2C/C2C transforms
	<li> \c outputIm : The base pointer of the “Imaginary” part output buffer for Planar R2C/C2C transforms
	<li> \c outoffset : Index of the current element  of the output buffer from the start
	<li> \c userdata : Buffer containing optional caller specified data. The userdata pointer is useful
	for passing any supplementary data to the callback function. For example, buffer having convolution
	filter data or any scalar value. The userdata can be of any custom data type/structure, in which case,
	you have to declare the custom data type and include it along with the callback function string. </li>
	<li> \c fftoutput : The result computed by clFFT for the element corresponding to outoffset argument
	<li> \c localmem : Pointer to local memory. This memory is allocated by library based on the size you specify
	and is subject to local memory availability. </li>
</ul>

@subsection SampleCallbackCode Sample Callback Code

@code
//**************************************************************************
//* Step 1 : Store the pre and post callback function in a string.
//**************************************************************************
const char* precallbackstr = "float2 pre_mulval(__global void* input, \n
                                  uint inoffset,                      \n
                                  __global void* userdata,            \n
                                  __local void* localmem)             \n
				{                                                            \n
				float scalar = *((__global float*)userdata + inoffset);      \n
				float2 ret = *((__global float2*)input + inoffset) * scalar; \n
				return ret;                                                  \n
				}                                                            \n";

const char* postcallbackstr = "void post_mulval(__global void* output, \n
                                  uint outoffset,                      \n
                                  __global void* userdata,             \n
								  float2 fftoutput,                    \n
                                  __local void* localmem)              \n
				{                                                      \n
				float scalar = *((__global float*)userdata + outoffset);      \n
				*((__global float2*)output + outoffset) = fftoutput * scalar; \n
				}                                                             \n";
				
//**************************************************************************
//* Step 2 : Initialize arguments, if any, required by the callback.
//**************************************************************************
int h_preuserdata[N] = {  };
cl_mem preuserdata = clCreateBuffer(context, CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR, sizeof(float) * N,  (void*)h_preuserdata, NULL);

int h_postuserdata[N] = {  };
cl_mem postuserdata = clCreateBuffer(context, CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR, sizeof(float) * N,  (void*)h_postuserdata, NULL);

//**************************************************************************
//* Step 3 : Register the callback.
//**************************************************************************

status = clfftSetPlanCallback(plan_handle, "pre_mulval", precallbackstr, 0, PRECALLBACK, &preuserdata, 1);

status = clfftSetPlanCallback(plan_handle, "post_mulval", postcallbackstr, 0, POSTCALLBACK, &postuserdata, 1);

//**************************************************************************
//* Step 4 : Bake plan and enqueue transform.
//**************************************************************************
status = clfftBakePlan( plan_handle, 1, &queue, NULL, NULL );

status = clfftEnqueueTransform( plan_handle, dir, 1, &queue, 0, NULL, &outEvent,
			&input_buffers[ 0 ], buffersOut, clMedBuffer );
@endcode

@subsection CallbackNotes Important Notes on Callback

<ol>
	<li> The caller is responsible to provide a callback function in string form that matches the function prototype based on the type of callback, type of transform(real/complex), and whether LDS is used.
	<li> clFFT considers the value returned by the pre-callback function as the new value of the input at the index corresponding to the *inoffset* argument.
	<li> Callback function can request for local memory for its own use. If the requested amount of local memory is available on the device, clFFT passes a pointer to the local memory when it invokes the callback function.
	<li> clFFT may invoke the FFT kernels several times depending on the input parameters. However, the pre-callback function provided by the caller is invoked only once for each point in the input. Similarly, it calls the post-callback function only once for each point in the output.
	<li> If clFFT is implementing a given FFT in multiple phases, it calls the pre-callback function only from the first phase kernel. Similarly, it calls the post-callback function only from the last phase kernel.
</ol>

 */
