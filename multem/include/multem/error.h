/*
 * This file is part of MULTEM.
 * Copyright (C) 2019 Diamond Light Source
 *
 *  Author: James Parkhurst
 *
 * MULTEM is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * MULTEM is distributed in the hope that it will be useful, 
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with MULTEM. If not, see <http:// www.gnu.org/licenses/>.
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
