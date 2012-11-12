/*
 *
 *                This source code is part of
 *
 *                 G   R   O   M   A   C   S
 *
 *          GROningen MAchine for Chemical Simulations
 *
 * Written by David van der Spoel, Erik Lindahl, Berk Hess, and others.
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2009, The GROMACS development team,
 * check out http://www.gromacs.org for more information.

 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 *
 * If you want to redistribute modifications, please consider that
 * scientific software is very special. Version control is crucial -
 * bugs must be traceable. We will be happy to consider code for
 * inclusion in the official distribution, but derived work must not
 * be called official GROMACS. Details are found in the README & COPYING
 * files - if they are missing, get the official version at www.gromacs.org.
 *
 * To help us fund GROMACS development, we humbly ask that you cite
 * the papers on the package - you can find them in the top README file.
 *
 * For more info, check our website at http://www.gromacs.org
 */
/*! \file
 * \brief
 * Declares gmx::AnalysisDataPlotModule for plotting data (into a file).
 *
 * \inpublicapi
 * \ingroup module_analysisdata
 * \author Teemu Murtola <teemu.murtola@cbr.su.se>
 */
#ifndef GMX_ANALYSISDATA_MODULES_PLOT_H
#define GMX_ANALYSISDATA_MODULES_PLOT_H

#include <string>

#include <boost/shared_ptr.hpp>

#include "../datamodule.h"
#include "../../options/timeunitmanager.h"
#include "../../utility/common.h"

namespace gmx
{

class Options;
class SelectionCollection;

/*! \brief
 * Common settings for data plots.
 *
 * \inpublicapi
 * \ingroup module_analysisdata
 */
class AnalysisDataPlotSettings
{
    public:
        //! Constructs default analysis plot settings.
        AnalysisDataPlotSettings();

        //! Returns the selection collection set with setSelectionCollection().
        const SelectionCollection *selectionCollection() const
        {
            return selections_;
        }
        //! Returns the time unit set with setTimeUnit().
        TimeUnit timeUnit() const { return timeUnit_; }
        /*! \brief
         * Returns the plot format.
         *
         * \todo Use a proper enum.
         */
        int plotFormat() const { return plotFormat_; }

        /*! \brief
         * Set selection collection to print as comments into the output.
         *
         * Formatted selection text from all selections in \p selections is
         * printed as comments in the output file.
         * If this method is not called, no selection information is written
         * to the output.
         */
        void setSelectionCollection(const SelectionCollection *selections);
        /*! \brief
         * Sets the time unit for the plot.
         *
         * The value is used only if AbstractPlotModule::setXAxisIsTime() is
         * called, in which case it is used to print the appropriate axis label
         * and to scale the values.
         * If not called, the default time unit is ps.
         */
        void setTimeUnit(TimeUnit timeUnit) { timeUnit_ = timeUnit; }


        /*! \brief
         * Adds common options for setting plot options.
         *
         * \param[in,out] options Options object to which options are added.
         */
        void addOptions(Options *options);

    private:
        const SelectionCollection *selections_;
        TimeUnit timeUnit_;
        int plotFormat_;
};

/*! \brief
 * Abstract data module for writing data into a file.
 *
 * Implements features common to all plotting modules.  Subclasses implement
 * features specific to certain applications (AnalysisDataPlotModule implements
 * straightforward plotting).
 *
 * By default, the data is written into an xvgr file, according to the
 * options read from the Options object given to the constructor.
 * For non-xvgr data, it's possible to skip all headers by calling
 * setPlainOutput().
 *
 * Multipoint data is supported, in which case all the points are written to
 * the output, in the order in which they are added to the data.  A single
 * output line corresponds to a single frame.  In most cases with multipoint
 * data, setPlainOutput() should be called since the output does not make sense
 * as an xvgr file, but this is not enforced.
 *
 * \ingroup module_analysisdata
 */
class AbstractPlotModule : public AnalysisDataModuleInterface
{
    public:
        virtual ~AbstractPlotModule();

        /*! \brief
         * Set common settings for the plotting.
         */
        void setSettings(const AnalysisDataPlotSettings &settings);
        /*! \brief
         * Set the output file name.
         *
         * If no file name is set (or if \p filename is empty), no output occurs.
         */
        void setFileName(const std::string &filename);
        /*! \brief
         * Set plain output.
         *
         * If \p bPlain is true, no xvgr headers are written to the file.
         * In this case, only setOmitX(), setXFormat(), and setYFormat()
         * methods have any effect on the output.
         */
        void setPlainOutput(bool bPlain);
        /*! \brief
         * Omit the X coordinates from the output.
         *
         * This method only makes sense when combined with setPlainOutput().
         */
        void setOmitX(bool bOmitX);
        /*! \brief
         * Set plot title.
         */
        void setTitle(const char *title);
        /*! \brief
         * Set plot subtitle.
         */
        void setSubtitle(const char *subtitle);
        /*! \brief
         * Set X axis label.
         */
        void setXLabel(const char *label);
        /*! \brief
         * Treat X axis as time.
         *
         * Sets the label for the axis accordingly and also scales output to
         * take into account the correct time unit.
         */
        void setXAxisIsTime();
        /*! \brief
         * Set Y axis label.
         */
        void setYLabel(const char *label);
        /*! \brief
         * Add legend from an array of strings.
         *
         * Multiple calls to setLegend() and/or appendLegend() are added
         * together.
         */
        void setLegend(int nsets, const char * const *setname);
        /*! \brief
         * Add a legend string for the next data set.
         *
         * Multiple calls to setLegend() and/or appendLegend() are added
         * together.
         */
        void appendLegend(const char *setname);
        /*! \brief
         * Set field width and precision for X value output.
         */
        void setXFormat(int width, int precision, char format = 'f');
        /*! \brief
         * Set field width and precision for Y value output.
         */
        void setYFormat(int width, int precision, char format = 'f');

        virtual int flags() const;

        virtual void dataStarted(AbstractAnalysisData *data);
        virtual void frameStarted(const AnalysisDataFrameHeader &header);
        virtual void pointsAdded(const AnalysisDataPointSetRef &points) = 0;
        virtual void frameFinished(const AnalysisDataFrameHeader &header);
        virtual void dataFinished();

    protected:
        /*! \cond libapi */
        AbstractPlotModule();
        //! Creates AbstractPlotModule and assign common settings.
        explicit AbstractPlotModule(const AnalysisDataPlotSettings &settings);

        //! Whether an output file has been opened.
        bool isFileOpen() const;
        /*! \brief
         * Appends a single value to the current output line.
         *
         * \param[in] value  Value to append.
         *
         * Should be used from pointsAdded() implementations in derived classes
         * to write out individual y values to the output.
         *
         * Must not be called if isFileOpen() returns false.
         */
        void writeValue(real value) const;
        //! \endcond

    private:
        class Impl;

        PrivateImplPointer<Impl> impl_;
};


/*! \brief
 * Plotting module for straightforward plotting of data.
 *
 * See AbstractPlotModule for common plotting options.
 *
 * \inpublicapi
 * \ingroup module_analysisdata
 */
class AnalysisDataPlotModule : public AbstractPlotModule
{
    public:
        AnalysisDataPlotModule();
        //! Creates AnalysisDataPlotModule and assign common settings.
        explicit AnalysisDataPlotModule(const AnalysisDataPlotSettings &settings);

        virtual void pointsAdded(const AnalysisDataPointSetRef &points);

        // Copy and assign disallowed by base.
};


/*! \brief
 * Plotting module specifically for data consisting of vectors.
 *
 * See AbstractPlotModule for common plotting options.
 *
 * \inpublicapi
 * \ingroup module_analysisdata
 */
class AnalysisDataVectorPlotModule : public AbstractPlotModule
{
    public:
        AnalysisDataVectorPlotModule();
        //! Creates AnalysisDataVectorPlotModule and assign common settings.
        explicit AnalysisDataVectorPlotModule(const AnalysisDataPlotSettings &settings);

        /*! \brief
         * Set whether to write X component.
         */
        void setWriteX(bool bWrite);
        /*! \brief
         * Set whether to write Y component.
         */
        void setWriteY(bool bWrite);
        /*! \brief
         * Set whether to write Z component.
         */
        void setWriteZ(bool bWrite);
        /*! \brief
         * Set whether to write norm of the vector.
         */
        void setWriteNorm(bool bWrite);
        /*! \brief
         * Set mask for what to write.
         */
        void setWriteMask(bool bWrite[4]);

        virtual void pointsAdded(const AnalysisDataPointSetRef &points);

    private:
        bool bWrite_[4];

        // Copy and assign disallowed by base.
};

//! Smart pointer to manage an AnalysisDataPlotModule object.
typedef boost::shared_ptr<AnalysisDataPlotModule>
AnalysisDataPlotModulePointer;
//! Smart pointer to manage an AnalysisDataVectorPlotModule object.
typedef boost::shared_ptr<AnalysisDataVectorPlotModule>
AnalysisDataVectorPlotModulePointer;

} // namespace gmx

#endif
