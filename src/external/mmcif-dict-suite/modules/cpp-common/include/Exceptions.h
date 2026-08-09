//$$FILE$$
//$$VERSION$$
//$$DATE$$
//$$LICENSE$$


/**
** \file Exceptions.h
**
** Declarations of exceptions in RCSB code.
*/


#ifndef EXCEPTIONS_H
#define EXCEPTIONS_H


#include <stdexcept>
#include <string>


/// Base class for all RCSB exceptions
class RcsbException : public std::exception
{
  protected:
    std::string _message;

  public:
    RcsbException(const std::string& message = std::string(),
      const std::string& location = std::string());
    ~RcsbException() throw();

    void AppendMessage(const std::string& message = std::string(),
      const std::string& location = std::string());

    const char* what() const throw();
};


/// Empty value exception (e.g. NULL pointer, empty string)
class EmptyValueException : public RcsbException
{

  public:
    EmptyValueException(const std::string& message = std::string(),
      const std::string& location = std::string());
    ~EmptyValueException() throw();

};


/// Object not found (thrown everywhere except from .find() methods)
class NotFoundException : public RcsbException
{

  public:
    NotFoundException(const std::string& message = std::string(),
      const std::string& location = std::string());
    ~NotFoundException() throw();

};


/// Object already exists
class AlreadyExistsException : public RcsbException
{

  public:
    AlreadyExistsException(const std::string& message = std::string(),
      const std::string& location = std::string());
    ~AlreadyExistsException() throw();

};


/// Empty container
class EmptyContainerException : public RcsbException
{

  public:
    EmptyContainerException(const std::string& message = std::string(),
      const std::string& location = std::string());
    ~EmptyContainerException() throw();

};


/// File mode exception (e.g. attempt to write to read-only file, invalid mode.)
class FileModeException : public RcsbException
{

  public:
    FileModeException(const std::string& message = std::string(),
      const std::string& location = std::string());
    ~FileModeException() throw();

};


/// Invalid state exception (e.g. getting a row reference in a column-wise table/// )
class InvalidStateException : public RcsbException
{

  public:
    InvalidStateException(const std::string& message = std::string(),
      const std::string& location = std::string());
    ~InvalidStateException() throw();

};


/// Generic files related exception (e.g. read error, write errror, etc.)
class FileException : public RcsbException
{

  public:
    FileException(const std::string& message = std::string(),
      const std::string& location = std::string());
    ~FileException() throw();

};

/// Invalid command line options
class InvalidOptionsException : public RcsbException
{
  public:
    InvalidOptionsException(const std::string& message = std::string(),
      const std::string& location = std::string());
    ~InvalidOptionsException() throw();
};

/// Versions do not match
class VersionMismatchException : public RcsbException
{

  public:
    VersionMismatchException(const std::string& message = std::string(),
      const std::string& location = std::string());
    ~VersionMismatchException() throw();

};

#endif
