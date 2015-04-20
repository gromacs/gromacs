Analysis output data handling {#page_analysisdata}
=============================

The \ref module_analysisdata module provides support for common data analysis
tasks within the \ref page_analysisframework.  The basic approach used in the
module is visualized below:

\dot
  digraph analysisdata_overview {
    rankdir = BT
    dataobject [label="data object\n(subclass of gmx::AbstractAnalysisData)"]
    datamodule1 [label="data module\n(implements gmx::AnalysisDataModuleInterface)"]
    datamodule2 [label="data module\nthat also provides data"]
    datamodule3 [label="data module"]
    datamodule1 -> dataobject
    datamodule2 -> dataobject
    datamodule3 -> datamodule2
  }
\enddot

Typically, an analysis tool provides its raw data output through one or more
gmx::AnalysisData objects (the root _data object_ in the diagram above).
This object provides only storage for the data.

To perform operations on the data, one or more _data modules_ can be attached
to the data object.  Examples of such operations are averaging, histogramming,
and plotting the data into a file.  Some data modules are provided by the \ref
module_analysisdata module.  To implement new ones, it is necessary to create a
class that implements gmx::AnalysisDataModuleInterface.

In many cases, such data modules also provide data that can be processed
further, acting as data objects themselves.  This makes it possible to attach
further data modules to form a processing chain.  In simple cases, such a chain
ends in a module that writes the data into a file, but it is also possible to
access the data in a data object (whether a plain data object or a data module)
programmatically to do further computation or post-processing outside the
framework.  To do this, the data object typically needs to be told in advance
such that it knows to store the data permanently even if attached modules do
not require it.

The modules can do their processing online, i.e., as the data is produced.
If all the attached modules support this, it is not necessary to store all the
raw data in memory.  The module design also supports processing frames in
parallel: in such cases, the data may become available out of order.  In
particular for writing the per-frame data into a file, but also for other types
of post-processing, it is necessary to reorder the data sequentially.  This is
implemented once in the framework, and analysis tools do not need to worry,
other than using the provided API.


Structure of data
=================

At the highest level, data can be structured into separate
gmx::AbstractAnalysisData objects that operate independently.  Each such object
has an independent set of post-processing modules.

Within a gmx::AbstractAnalysisData object, data is structured along three
"dimensions":

 - _frames_: There is one or more frames in each data object.  For raw data
   produced by an analysis tool, these typically correspond to input trajectory
   frames.  For other data set, it can be viewed as an X axis of a graph.
 - _data sets_: There is one or more data sets in each data object.  For most
   purposes, data sets work independently (i.e., the post-processing modules
   operate on each data set separately), but some modules reduce the data sets
   into single columns in the output.  The main purpose for using multiple data
   sets is to share the same post-processing chain for multiple sets of data
   (e.g., multiple RDFs computed by the same tool in one pass), in particular
   for cases where the number of data sets is not known at compile time.
   Note that each data set contains the same number of frames.
 - _columns_: There is one or more columns in each data set.  Different data
   sets can contain a different number of columns.  Each column in a frame can
   contain a single value (see below for supported values).

Programmatically the data within each frame is organized into _point sets_.
Each point set consists of a continuous range of columns from a single data
set.  There are two types of data:

 - _simple_: For each frame, there is exactly one point set for each data set,
   and that point set spans all columns in that data set.
 - _multipoint_: For each frame, there can be any number of point sets, and
   they may span arbitrary columns.  It is allowed that point sets overlap,
   i.e., that multiple point sets specify a value for the same column.

The main purpose of multipoint data is to support cases where it is not known
in advance how many values there will be for each frame, or where that number
is impractically large.  The need to do this is mainly a matter of
performance/implementation complexity tradeoff: with a more complex internal
implementation, it would be possible to support larger data sets without a
performance/memory impact they currently impose.  The current implementation
places the burden of deciding on the appropriate usage pattern on the user
code, allowing for much simpler internal implementation.

An individual value (identified by frame, data set, and column) consists of a
single value of type `real`, an optional error value, and some flags.
The flags identify what parts of the value are really available.  The following
states are possible:
 - _present_: The value is set.
 - _missing_: The value is marked as missing by the data source.  In this
   state, the value can still be accessed, and the returned `real` value has
   some meaning.  Different data modules handle these cases differently.
 - _unset_: The value is not set.  It is not allowed to access the value for
   other than querying the state.  Data modules that ignore missing values
   (by skipping all values not _present_) can also handle unset values.
   Other data modules typically do not allow unset values.


Data provider classes
=====================

The base class for all data objects (including data modules that provide data)
is gmx::AbstractAnalysisData.  This class provides facilities for attaching
data modules to the data, and to query the data.  It does not provide any
methods to alter the data; all logic for managing the actual data is in derived
classes.

The main root (non-module) data object class for use in analysis tools is
gmx::AnalysisData.  This class provides methods to set properties of the data,
and to add frames to it.  The interface is frame-based: you construct one frame
at a time, and after it is finished, you move to the next frame.  The frames
are not constructed directly using gmx::AnalysisData, but instead separate
_data handles_ are used.  This is explained in more detail below under
\ref section_analysisdata_parallelization.

For simple needs and small amounts of data, gmx::AnalysisArrayData is also
provided.  This class allows for all the data to be prepared in memory as a
single big array, and allows random access to the data while setting the
values.  When all the values are set to their final values, it then notifies
the attached data modules by looping over the array.


Parallelization {#section_analysisdata_parallelization}
===============

One major driver for the design of the analysis data module has been to provide
support for transparently processing multiple frames in parallel.  In such
cases, output data for multiple frames may be constructed simultaneously, and
must be ordered correctly for some data modules, such as writing it into a
file.  This ordering is taken care of by the framework, allowing the analysis
tool writer to concentrate on the actual analysis task.

From a user's point of view, the main player in this respect is the
gmx::AnalysisData object.  If there are two threads doing the processing in
parallel, it allows creating a separate gmx::AnalysisDataHandle for each
object.  Each of these handles can be used independently to construct frames
into the output data, and the gmx::AnalysisData object internally takes care of
notifying the modules correctly.  If necessary, it stores finished frames into
a temporary buffer until all preceding frames have also been finished.

For increased efficiency, some data modules are also parallelization-aware:
they have the ability to process the data in any order, allowing
gmx::AnalysisData to notify them as soon as a frame becomes available.
If there are only parallel data modules attached, no frame reordering or
temporary buffers are needed.  If a non-parallel data module is attached to a
parallel data module, then that parallel data module takes the responsibility
of ordering its output frames.  Ideally, such data modules produce
significantly less data than what they take in, making it cheaper to do the
ordering only at this point.

Currently, no parallel runner has been implemented, but it is likely that
applicable tools written to use the framework require minimal or no changes to
take advantage of frame-level parallelism once such a runner materializes.


Provided data processing modules
================================

Data modules provided by the \ref module_analysisdata module are listed below
with a short description.  See the documentation of the individual classes for
more details.
Note that this list is manually maintained, so it may not always be up-to-date.
A comprehensive list can be found by looking at the inheritance graph of
gmx::AnalysisDataModuleInterface, but the list here is more user-friendly.

<dl>
<dt>gmx::AnalysisDataAverageModule</dt>
<dd>
Computes averages and standard deviations for columns in input data.
One output value for each input column.
</dd>
<dt>gmx::AnalysisDataFrameAverageModule</dt>
<dd>
Computes averages for each frame in input data.
One output value for each input data set for each frame.
</dd>
<dt>gmx::AnalysisDataBinAverageModule</dt>
<dd>
Computes averages within bins.  Input is pairs of values, where the first
value defines the bin, and the second value sets the value to accumulate into
the average within the bin.
One output histogram for each input data set.
</dd>
<dt>gmx::AnalysisDataSimpleHistogramModule</dt>
<dd>
Computes histograms.  All values within a data set are added into a histogram.
One output histogram for each input data set.
Provides the histogram for each input frame separately, and also the full
histogram over frames (through an internal submodule).
</dd>
<dt>gmx::AnalysisDataWeightedHistogramModule</dt>
<dd>
Computes histograms.  Input is pairs of values, where the first value defines
the bin, and the second value sets the value to add into that bin.
Output like with gmx::AnalysisDataSimpleHistogramModule.
</dd>
<dt>gmx::AnalysisDataLifetimeModule</dt>
<dd>
Computes lifetime histograms.  For each input column, determines the time
intervals during which a value is continuously present/non-zero, and creates a
histogram from the lengths of these intervals.
One output histogram for each input data set.
</dd>
<dt>gmx::AnalysisDataPlotModule</dt>
<dt>gmx::AnalysisDataVectorPlotModule</dt>
<dd>
Writes data into a file.
</dd>
</dl>
