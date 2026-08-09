//$$FILE$$
//$$VERSION$$
//$$DATE$$
//$$LICENSE$$


#include <string>

#include "Exceptions.h"


using std::string;


RcsbException::RcsbException(const string& message,
  const string& location)
{

    AppendMessage(message, location);

}


RcsbException::~RcsbException() throw()
{

    _message.clear();

}


void RcsbException::AppendMessage(const string& message,
  const string& location)
{

    _message += "Message: \"" + message + "\" generated at: " +
      location + '\n';

}


const char* RcsbException::what() const throw()
{

    return(_message.c_str());

}


EmptyValueException::EmptyValueException(const string& message,
  const string& location) : RcsbException(message, location)
{

}


EmptyValueException::~EmptyValueException() throw()
{

}


NotFoundException::NotFoundException(const string& message,
  const string& location) : RcsbException(message, location)
{

}


NotFoundException::~NotFoundException() throw()
{

}


AlreadyExistsException::AlreadyExistsException(const string& message,
  const string& location) : RcsbException(message, location)
{

}


AlreadyExistsException::~AlreadyExistsException() throw()
{

}


EmptyContainerException::EmptyContainerException(const string& message,
  const string& location) : RcsbException(message, location)
{

}


EmptyContainerException::~EmptyContainerException() throw()
{

}


FileModeException::FileModeException(const string& message,
  const string& location) : RcsbException(message, location)
{

}


FileModeException::~FileModeException() throw()
{

}


InvalidStateException::InvalidStateException(const string& message,
  const string& location) : RcsbException(message, location)
{

}


InvalidStateException::~InvalidStateException() throw()
{

}


FileException::FileException(const string& message,
  const string& location) : RcsbException(message, location)
{

}


FileException::~FileException() throw()
{

}


InvalidOptionsException::InvalidOptionsException(const string& message,
  const string& location)
{
    _message = message;
}


InvalidOptionsException::~InvalidOptionsException() throw()
{

}


VersionMismatchException::VersionMismatchException(const string& message,
  const string& location) : RcsbException(message, location)
{

}


VersionMismatchException::~VersionMismatchException() throw()
{

}

