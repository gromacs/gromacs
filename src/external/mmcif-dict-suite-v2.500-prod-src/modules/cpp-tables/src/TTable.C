//$$FILE$$
//$$VERSION$$
//$$DATE$$
//$$LICENSE$$


/*!
** \file TTable.C
**
** \brief Implementation file for TTable class. 
*/


#include <stdexcept>
#include <string>
#include <vector>
#include <iostream>
#include <iomanip>

#include "Exceptions.h"
#include "GenString.h"
#include "TTable.h"


using std::out_of_range;
using std::copy;
using std::string;
using std::vector;
using std::ostream;
using std::setw;
using std::endl;


TTable::TTable() : _numCols(0)
#ifndef TTABLE_COLUMN_DELETE_AS_REMOVE
  , _numDelCols(0)
#endif
{


}


TTable::TTable(const TTable& inTable)
{

    _numCols = inTable._numCols;

    // This statement only copies pointer addresses and not the tuple
    // vectors referenced by those pointers.
    _tuples = inTable._tuples;

    // Here, tuple vectors themselves are copied and previously set
    // pointers are properly initialized.

    for (unsigned int tupleI = 0; tupleI < _tuples.size(); ++tupleI)
    {
        _tuples[tupleI] = new vector<string>(*(inTable._tuples[tupleI]));
    }

#ifndef TTABLE_COLUMN_DELETE_AS_REMOVE
    _delColMap = inTable._delColMap;
    _numDelCols = inTable._numDelCols;
#endif

}


TTable::~TTable()
{

    Clear();

}


TTable& TTable::operator=(const TTable& inTable)
{

    if (this != &inTable)
    {
        Clear();

        _numCols = inTable._numCols;

        // This statement only copies pointer addresses and not the tuple
        // vectors referenced by those pointers.
        _tuples = inTable._tuples;

        // Here, tuple vectors themselves are copied and previously set
        // pointers are properly initialized.

        for (unsigned int tupleI = 0; tupleI < _tuples.size();
          ++tupleI)
        {
            _tuples[tupleI] = new vector<string>(*(inTable._tuples[tupleI]));
        }
    }

#ifndef TTABLE_COLUMN_DELETE_AS_REMOVE
    _delColMap = inTable._delColMap;
    _numDelCols = inTable._numDelCols;
#endif

    return(*this);

}


void TTable::Clear()
{

    if (!_tuples.empty())
    {
        for (unsigned int tupleI = 0; tupleI < _tuples.size(); tupleI++)
        {
            _tuples[tupleI]->clear();
            delete _tuples[tupleI];
        }

        _tuples.clear();
    }

    _numCols = 0;

#ifndef TTABLE_COLUMN_DELETE_AS_REMOVE
    _delColMap.clear();
    _numDelCols = 0;
#endif

}


unsigned int TTable::AddTuple(const vector<string>& tuple)
{

    InsertTuple(_tuples.size(), tuple);

    return(_tuples.size());

}


void TTable::InsertTuple(const unsigned int tupleIndex,
  const vector<string>& tuple)
{

    InsertTuple(tupleIndex, tuple.begin(), tuple.end());

}


void TTable::InsertTuple(const unsigned int tupleIndex,
  vector<string>::const_iterator tupleBeg,
  vector<string>::const_iterator tupleEnd)
{

    if (tupleIndex > _tuples.size())
    {
        // Wrong tuple index.
        throw out_of_range("Invalid tuple index in TTable::InsertTuple");
    }

    if (tupleBeg > tupleEnd)
    {
        // Wrong tuple range.
        throw out_of_range("Invalid tuple range in TTable::InsertTuple");
    }

    unsigned int tupleSize = tupleEnd - tupleBeg;

    if ((GetNumColumns() != 0) && (tupleSize > GetNumColumns()))
    {
        // Wrong tuple index.
        throw out_of_range("Invalid tuple size in TTable::InsertTuple");
    }

    if ((tupleSize != 0) && (tupleIndex != 0) &&
      (_tuples[tupleIndex - 1]->empty()))
    {
        // Previous tuple empty. Cannot insert, since rectangularity will be
        // violated.
        throw out_of_range("Previous tuple empty TTable::InsertTuple");
    }

    if (_numCols == 0)
    {
        // For empty table, new tuple defines the initial number of columns.
        _numCols = tupleSize;
#ifndef TTABLE_COLUMN_DELETE_AS_REMOVE
        EnlargeColMap(_numCols);
#endif
    }

    _tuples.insert(_tuples.begin() + tupleIndex, new vector<string>);

    if (_tuples[tupleIndex]->empty())
    {
        // If table tuple empty, insert empty cells.
        _tuples[tupleIndex]->insert(_tuples[tupleIndex]->begin(), _numCols,
          string());
    }

    copy(tupleBeg, tupleEnd, _tuples[tupleIndex]->begin());

}


void TTable::FillTuple(const unsigned int tupleIndex, const vector<string>& tuple,
  const unsigned int fromColIndex)
{

    if (tupleIndex >= _tuples.size())
    {
        // Wrong tuple index.
        throw out_of_range("Invalid tuple index in TTable::FillTuple");
    }

    if (tuple.empty())
    {
        return;
    }

    if (fromColIndex > GetNumColumns())
    {
        // Wrong from column index.
        throw out_of_range("Invalid from column index in TTable::FillTuple");
    }

    if ((GetNumColumns() != 0) && (tuple.size() > GetNumColumns() - fromColIndex))
    {
        // Wrong tuple index.
        throw out_of_range("Invalid tuple size in TTable::FillTuple");
    }

    copy(tuple.begin(), tuple.end(), _tuples[tupleIndex]->begin() +
      IntColIndex(fromColIndex));

}


void TTable::GetTuple(vector<string>& tuple, const unsigned int tupleIndex,
      const unsigned int fromColIndex, unsigned int toColIndex)
{

    tuple.clear();

    if (tupleIndex >= _tuples.size())
    {
        // Wrong tuple index.
        throw out_of_range("Invalid tuple index in TTable::GetTuple");
    }

    for (unsigned int colI = fromColIndex; colI < toColIndex; ++colI)
    {
        tuple.push_back(operator()(tupleIndex, colI));
    }

}


const vector<string>& TTable::GetTuple(const unsigned int tupleIndex)
{

    if (tupleIndex >= _tuples.size())
    {
        // Wrong tuple index.
        throw out_of_range("Invalid tuple index in TTable::GetTuple");
    }

    return(*(_tuples[tupleIndex]));

}


void TTable::ClearTuple(const unsigned int tupleIndex)
{

    if (tupleIndex >= _tuples.size())
    {
        // Wrong tuple index.
        throw out_of_range("Invalid tuple index in TTable::ClearTuple");
    }

    for (unsigned int colI = 0; colI < _tuples[tupleIndex]->size();
      ++colI)
    {
        (*_tuples[tupleIndex])[colI].clear();
    }

}


void TTable::DeleteTuple(const unsigned int tupleIndex)
{

    if (tupleIndex >= _tuples.size())
    {
        // Wrong tuple index.
        throw out_of_range("Invalid tuple index in TTable::DeleteTuple");
    }

    _tuples[tupleIndex]->clear();

    delete _tuples[tupleIndex];

    _tuples.erase(_tuples.begin() + tupleIndex);

    if (_tuples.empty())
    {
#ifndef TTABLE_COLUMN_DELETE_AS_REMOVE
        ReduceColMap(_numCols);
#endif
        _numCols = 0; 
    }

}


unsigned int TTable::AddColumn(const vector<string>& col)
{

    return(InsertColumn(GetNumColumns(), col));

}


unsigned int TTable::InsertColumn(const unsigned int atColIndex,
  const vector<string>& col)
{

    InsertColumn(atColIndex, col.begin(), col.end());

    return(GetNumColumns());
    
}


void TTable::InsertColumn(const unsigned int atColIndex,
  vector<string>::const_iterator colBeg, vector<string>::const_iterator colEnd)
{

    if (atColIndex > GetNumColumns())
    {
        // Wrong column index.
        throw out_of_range("Invalid column index in TTable::InsertColumn");
    }

    if (colBeg > colEnd)
    {
        // Wrong tuple range.
        throw out_of_range("Invalid column range in TTable::InsertColumn");
    }

    unsigned int colSize = colEnd - colBeg;
    if (colSize == 0)
    {
        throw EmptyValueException("Empty column",
          "TTable::InsertColumn");
    }

    if ((!_tuples.empty()) && (colSize > _tuples.size()))
    {
        throw out_of_range("Invalid column size in TTable::InsertColumn");
    }

    if ((GetNumColumns() == 0) && (colSize > _tuples.size()))
    {
        unsigned int numTuples = colSize - _tuples.size();
        for (unsigned int tupleI = 0; tupleI < numTuples; ++tupleI)
        {
            _tuples.push_back(new vector<string>);
        }
    }

#ifndef TTABLE_COLUMN_DELETE_AS_REMOVE
    if (_numDelCols != 0)
    {
        UnMarkColDeleted(atColIndex);

        if (colSize == 0)
            ClearColumn(atColIndex);
    }
    else
    {
        EnlargeColMap(1);
#endif
        _numCols++;

        unsigned int intRowIndex = IntColIndex(atColIndex);

        for (unsigned int tupleI = 0; tupleI < _tuples.size(); ++tupleI)
        {
            _tuples[tupleI]->insert(_tuples[tupleI]->begin() + intRowIndex,
              1, string());
        }

#ifndef TTABLE_COLUMN_DELETE_AS_REMOVE
    }
#endif

    FillColumn(atColIndex, colBeg, colEnd);

}


void TTable::FillColumn(const unsigned int colIndex, const vector<string>& col,
  const unsigned int fromTupleIndex)
{

    FillColumn(colIndex, col.begin(), col.end(), fromTupleIndex);

}


void TTable::FillColumn(const unsigned int colIndex,
  vector<string>::const_iterator colBeg, vector<string>::const_iterator colEnd,
  const unsigned int fromTupleIndex)
{

    if (colIndex >= GetNumColumns())
    {
        // Wrong column index.
        throw out_of_range("Invalid column index in TTable::FillColumn");
    }

    if (colBeg > colEnd)
    {
        // Wrong tuple range.
        throw out_of_range("Invalid column range in TTable::FillColumn");
    }

    if (fromTupleIndex >= _tuples.size())
    {
        // Wrong tuple range.
        throw out_of_range("Invalid from tuple index in TTable::FillColumn");
    }

    unsigned int colSize = colEnd - colBeg;

    if (colSize == 0)
    {
        return;
    }

    if (colSize > (_tuples.size() - fromTupleIndex))
    {
        throw out_of_range("Invalid column size in TTable::FillColumn");
    }

    unsigned int intRowIndex = IntColIndex(colIndex);

    unsigned int tupleI = 0;
    for (vector<string>::const_iterator currIter = colBeg; currIter < colEnd;
      ++currIter, ++tupleI)
    {
        (*_tuples[fromTupleIndex + tupleI])[intRowIndex] =
          *currIter;
    }

}


void TTable::ClearColumn(const unsigned int colIndex)
{

    if (colIndex >= GetNumColumns())
    {
        // Wrong column index.
        throw out_of_range("Invalid column index in TTable::ClearColumn");
    }

    unsigned int intRowIndex = IntColIndex(colIndex);

    for (unsigned int tupleI = 0; tupleI < _tuples.size(); ++tupleI)
    {
        (*_tuples[tupleI])[intRowIndex].clear();
    }

}


void TTable::DeleteColumn(const unsigned int colIndex)
{

    if (colIndex >= GetNumColumns())
    {
        // Wrong column index.
        throw out_of_range("Invalid column index in TTable::DeleteColumn");
    }

#ifdef TTABLE_COLUMN_DELETE_AS_REMOVE
    for (unsigned int tupleI = 0; tupleI < _tuples.size(); ++tupleI)
    {
        _tuples[tupleI]->erase(_tuples[tupleI]->begin() + colIndex);
    }
#endif

#ifndef TTABLE_COLUMN_DELETE_AS_REMOVE
    MarkColDeleted(colIndex);
#else
    _numCols--;
#endif

}


ostream& operator<<(ostream& out, const TTable& tTable)
{

    for (unsigned int colI = 0; colI < tTable.GetNumColumns(); ++colI)
    {
        for (unsigned int tupleI = 0; tupleI < tTable.GetNumTuples(); ++tupleI)
        {
            out << setw(10) << tTable(tupleI, colI) << " ";
        }

        out << endl;
    }

    return(out);

}


void TTable::GetColumn(vector<string>& col, const unsigned int colIndex,
  const unsigned int fromTupleIndex, unsigned int toTupleIndex)
{

    col.clear();

    if (colIndex >= GetNumColumns())
    {
        // Wrong column index.
        throw out_of_range("Invalid column index in TTable::GetColumn");
    }

    if (fromTupleIndex >= toTupleIndex)
    {
        // To index less than from index.
        throw out_of_range("Invalid tuple index range in TTable::GetColumn");
    }

    if (fromTupleIndex >= _tuples.size())
    {
        // Wrong from tuple index.
        throw out_of_range("Invalid from tuple index in TTable::GetColumn");
    }

    if (toTupleIndex > _tuples.size())
    {
        // Wrong to tuple index.
        throw out_of_range("Invalid to tuple index in TTable::GetColumn");
    }

    // toTupleIndex is inclusive, so one more spot for it needs to be reserved
    col.reserve(toTupleIndex - fromTupleIndex);

    unsigned int intRowIndex = IntColIndex(colIndex);

    for (unsigned int tupleI = fromTupleIndex; tupleI < toTupleIndex;
      tupleI++)
    {
        col.push_back((*_tuples[tupleI])[intRowIndex]);
    }

}


string& TTable::operator()(const unsigned int tupleIndex,
  const unsigned int colIndex)
{

    if (tupleIndex >= _tuples.size())
    {
        // Wrong tuple index.
        throw out_of_range("Invalid tuple index in TTable::UpdateCell");
    }

    if (colIndex >= GetNumColumns())
    {
        // Wrong column index.
        throw out_of_range("Invalid column index in TTable::UpdateCell");
    }

    return((*_tuples[tupleIndex])[IntColIndex(colIndex)]);

}


const string& TTable::operator()(const unsigned int tupleIndex,
  const unsigned int colIndex) const
{

    if (tupleIndex >= _tuples.size())
    {
        // Wrong tuple index.
        throw out_of_range("Invalid tuple index in TTable::GetCell");
    }

    if (colIndex >= GetNumColumns())
    {
        // Wrong column index.
        throw out_of_range("Invalid column index in TTable::GetCell");
    }

    return((*_tuples[tupleIndex])[IntColIndex(colIndex)]);

}


int TTable::Write(Serializer* ser, unsigned int& size)
{

    // Write the number of columns
    // VLAD - THIS MAY NOT BE NECESSARY, SINCE IT CAN BE OBTAINED FROM
    // COLUMN LENGTH
    UInt32 firstIndex = ser->WriteUInt32(GetNumColumns());

    // Write the number of tuples
    UInt32 currIndex = ser->WriteUInt32(_tuples.size());

    // Write the cells
    vector<string> tupleToWrite;
    for (unsigned int tupleI = 0; tupleI < _tuples.size(); ++tupleI)
    {
        for (unsigned int colI = 0; colI < GetNumColumns(); ++colI)
        {
            tupleToWrite.push_back(operator()(tupleI, colI));
        }
        currIndex = ser->WriteStrings(tupleToWrite);
        tupleToWrite.clear();
    }

    size = currIndex - firstIndex + 1;

    return (firstIndex);

}


int TTable::Read(UInt32 index, Serializer* ser)
{

    // unsigned int numCols = ser->ReadUInt32(index);
    // This is not used and is not read, only skipped
    index++;

    unsigned int numTuples = ser->ReadUInt32(index);
    index++;

    for (unsigned int tupleI = 0; tupleI < numTuples; ++tupleI)
    {
        vector<string> tmpTuple;
        ser->ReadStrings(tmpTuple, index);
        index++;
        InsertTuple(tupleI, tmpTuple);
    }

    return(index);

}


#ifndef TTABLE_COLUMN_DELETE_AS_REMOVE
void TTable::EnlargeColMap(const unsigned int numCols)
{

    unsigned int firstIndex = _delColMap.size();

    for (unsigned int colI = firstIndex;
      colI < (firstIndex + numCols); colI++)
    {
        _delColMap.push_back(colI);
    }

}


void TTable::ReduceColMap(const unsigned int numCols)
{

    _delColMap.erase(_delColMap.end() - numCols, _delColMap.end());

}


void TTable::MarkColDeleted(const unsigned int colIndex)
{

    unsigned int saveRowIndex = _delColMap[colIndex];

    unsigned int index;

    for (index = colIndex; index < (_numCols - 1); index++)
    {
        _delColMap[index] = _delColMap[index + 1];
    }

    _delColMap[index] = saveRowIndex;

    ++_numDelCols;

}


void TTable::UnMarkColDeleted(const unsigned int colIndex)
{

    unsigned int availRowIndex = _delColMap[_numCols - 1];

    unsigned int index;

    for (index = _numCols - 1; index > colIndex; --index)
    {
        _delColMap[index] = _delColMap[index - 1];
    }

    _delColMap[index] = availRowIndex;

    --_numDelCols;

}
#endif

