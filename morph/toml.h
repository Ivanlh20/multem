// toml_parser.h

#pragma once

#include <string>
#include <map>
#include <stdexcept>
#include "math_mt.h"

namespace toml {

template <typename T>
T convertToType(const std::string& str);

template <>
bool convertToType<bool>(const std::string& str);

template <>
dt_cfloat32 convertToType<dt_cfloat32>(const std::string& str);

class TomlParser {
public:
    TomlParser(const std::string& filename);

    void parse();

    template <typename T>
    T get(const std::string& section, const std::string& key);

private:
    std::string filename_;
    std::map<std::string, std::map<std::string, std::string>> tomlData_;

    std::string trim(const std::string& str);
};

}  // namespace toml

#include "toml.cpp"
