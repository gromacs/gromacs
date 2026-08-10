//$$FILE$$
//$$VERSION$$
//$$DATE$$
//$$LICENSE$$


#include <iostream>
#include <string>
#include <vector> 
#include <algorithm>

#include "GenString.h"
#include "GenCont.h"


using std::string;
using std::vector;
using std::find_if;

using std::ostream;

//  Replacement for the std::bind2nd() call that C++17 removed, as a named
//  functor so the file still compiles as C++98 (see GenString.C).
namespace {

class StringEqualToValue
{
  public:
    StringEqualToValue(const Char::eCompareType compareType, const std::string& value)
      : _compareType(compareType), _value(value) {}

    bool operator()(const std::string& s) const
    {
        return StringEqualTo(_compareType)(s, _value);
    }

  private:
    Char::eCompareType _compareType;
    std::string        _value;
};

}  // namespace



ostream& operator<<(ostream& out, const vector<string>& contVector)
{
    out << "[ ";

    for (vector<string>::const_iterator it = contVector.begin();
      it != contVector.end(); ++it)
    {
        out << "\"" << (*it) << "\" ";
    }

    out << "]";

    return (out);
}


GenCont::GenCont()
{

}


GenCont::~GenCont()
{

}


bool GenCont::IsInVector(const string& element,
  const vector<string>& contVector, const Char::eCompareType compareType)
{
    vector<string>::const_iterator where = find_if(contVector.begin(),
      contVector.end(),
      StringEqualToValue(compareType, element));
 
    return (where != contVector.end());
}


bool GenCont::IsInVectorCi(const string& element,
  const vector<string>& contVector)
{
    return (IsInVector(element, contVector, Char::eCASE_INSENSITIVE));
}

