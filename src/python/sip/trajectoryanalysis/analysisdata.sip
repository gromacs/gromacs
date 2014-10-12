
class AbstractAnalysisData /NoDefaultCtors/ {
%TypeHeaderCode
#include <gromacs/analysisdata/abstractdata.h>
%End
public:
    bool isMultipoint() const;
    int dataSetCount() const;
    int columnCount(int dataSet) const;
    int columnCount() const;
    virtual int frameCount() const = 0;
    AnalysisDataFrameRef tryGetDataFrame(int index) const;
    AnalysisDataFrameRef getDataFrame(int index) const;
    bool requestStorage(int nframes);
};

class AnalysisDataFrameRef /NoDefaultCtors/ {
%TypeHeaderCode
#include <gromacs/analysisdata/dataframe.h>
%End
public:
    bool isValid() const;
    //const AnalysidDataFrameHeader& header() const;
    int frameIndex() const;
    double x() const;
    double dx() const;
    int pointSetCount() const;
    AnalysisDataPointSetRef pointSet(int index) const;
    double y(int i) const;
    double dy(int i) const;
    bool present(int i) const;
    bool allPresent() const;
};

class AnalysisDataFrameHeader /NoDefaultCtors/ {
%TypeHeaderCode
#include <gromacs/analysisdata/dataframe.h>
%End
public:
    bool isValid() const;
    int index() const;
    double x() const;
    double dx() const;
};

class AnalysisDataPointSetRef /NoDefaultCtors/ {
%TypeHeaderCode
#include <gromacs/analysisdata/dataframe.h>
%End
public:
    const AnalysisDataFrameHeader& header() const;
    int frameIndex() const;
    double x() const;
    double dx() const;
    int dataSetIndex() const;
    int firstColumn() const;
    int columnCount() const;
    int lastColumn() const;
    //const std::vector<AnalysisDataValue>& values() const;
    double y(int i) const;
    double dy(int i) const;
    bool present(int i) const;
    bool allPresent() const;
};
