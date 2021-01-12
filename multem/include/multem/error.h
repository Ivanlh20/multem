/*
 *  error.h
 *
 *  Copyright (C) 2019 Diamond Light Source
 *
 *  Author: James Parkhurst
 *
 *  This code is distributed under the GPLv3 license, a copy of 
 *  which is included in the root directory of this package.
 */
#ifndef MULTEM_ERROR_H
#define MULTEM_ERROR_H

#include <stdexcept>
#include <string>

namespace multem {

  /**
   * An error class that also prints the file and line number
   */
  class Error : public std::runtime_error {
  public:
    
    Error(const std::string& what_arg)
      : std::runtime_error(what_arg) {}
    
    Error(const char* what_arg )
      : std::runtime_error(what_arg) {}
    
    Error(const std::string& file, std::size_t line, const std::string& message)
      : std::runtime_error(file + ":" + std::to_string(line) + " " + message) {}
  };

}

/**
 * Throw an error if the assertion fails
 */
#define MULTEM_ASSERT(assertion) \
  if (!(assertion)) { \
    throw multem::Error(__FILE__, __LINE__, "ASSERT (" #assertion ") failed"); \
  }

#endif
