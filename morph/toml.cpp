#include <iostream>
#include <fstream>
#include <map>
#include <string>
#include <stdexcept>
#include <algorithm>
#include <sstream>
#include "math_mt.h"
#include "toml.h"

namespace toml
{

template <typename T>
T convertToType(const std::string& str) {
    std::istringstream iss(str);
    T value;
    iss >> value;
    return value;
}

template <>
bool convertToType<bool>(const std::string& str) {
    if (str == "true") return true;
    if (str == "false") return false;
    throw std::runtime_error("Failed to convert string to boolean.");
}
    // Helper function to convert a string to a specified type
template <>
dt_cfloat32 convertToType<dt_cfloat32>(const std::string& str) {
    std::string realPartStr, imagPartStr;
    size_t posJ = str.find("j");
    size_t posI = str.find("i");

    // Check if both "j" and "i" appear in the string, which is invalid
    if ((posJ != std::string::npos) && (posI != std::string::npos)) {
        throw std::runtime_error("Error converting value: " + str + " to complex number. Both 'j' and 'i' found.");
    }

    // Determine the position of the imaginary indicator (either "j" or "i")
    size_t posImaginary = (posJ != std::string::npos) ? posJ : posI;

    // If neither "j" nor "i" is found, throw an exception
    if (posImaginary == std::string::npos) {
        throw std::runtime_error("Error converting value: " + str + " to complex number. No 'j' or 'i' found.");
    }

    // Find the '+' or '-' symbol before the imaginary indicator to split the real and imaginary parts
    size_t posSign = str.rfind('+', posImaginary);
    if (posSign == std::string::npos) {
        posSign = str.rfind('-', posImaginary);
    }

    if (posSign != std::string::npos) {
        realPartStr = str.substr(0, posSign);
        imagPartStr = str.substr(posSign, posImaginary - posSign);
    } else {
        realPartStr = "0";
        imagPartStr = str.substr(0, posImaginary);
    }

    // Convert the string parts to floats
    float realPart, imagPart;
    std::istringstream(realPartStr) >> realPart;
    std::istringstream(imagPartStr) >> imagPart;

    // Return the complex number
    return std::complex<float>(realPart, imagPart);
}



    TomlParser::TomlParser(const std::string& filename) : filename_(filename) {}

    void TomlParser::parse() {
        std::ifstream file(filename_);

        if (!file) {
            std::cerr << "Error opening file: " << filename_ << std::endl;
        }

        std::string line;
        std::string currentSection;
        while (std::getline(file, line)) {
            line = trim(line);  // Trim the line first
            if (!line.empty() && line[0] == '[' && line[line.length() - 1] == ']') {
                // This line is a section header
                currentSection = line.substr(1, line.length() - 2);
            } else {
                // This line is a key-value pair
                size_t equalsPos = line.find('=');
                if (equalsPos != std::string::npos) {
                    std::string key = line.substr(0, equalsPos);
                    std::string value = line.substr(equalsPos + 1);
                    key = trim(key); 
                    value = trim(value); 
                    // Remove quotes around string values
                    if (value.length() >= 2 && value[0] == '"' && value[value.length() - 1] == '"') {
                        value = value.substr(1, value.length() - 2);
                    }
                    tomlData_[currentSection][key] = value;
                }
            }
        }
    }

    template <typename T>
    T TomlParser::get(const std::string& section, const std::string& key) {
        std::string valueStr = tomlData_[section][key];
        try {
            return convertToType<T>(valueStr);
        } catch (const std::exception& e) {
            throw std::runtime_error("Error converting value: " + valueStr + " to type " + typeid(T).name());
        }
    }

    std::string TomlParser::trim(const std::string& str) {
        size_t first = str.find_first_not_of(' ');
        if (std::string::npos == first) {
            return str;
        }
        size_t last = str.find_last_not_of(' ');
        return str.substr(first, (last - first + 1));
    }

    
}  // namespace toml

// int main() {
//     toml::TomlParser config("input.toml");

//     config.parse();
//     auto nkr = config.get<int>("Parameters", "nkr");
//     auto floatValue = config.get<float>("Parameters", "floatValue");
//     auto operation = config.get<std::string>("Parameters", "operation");
//     auto complex = config.get<dt_cfloat32>("Parameters", "complex");
//     std::cout << "nkr = " << nkr << std::endl;
//     std::cout << "floatValue = " << floatValue << std::endl;
//     std::cout << "operation = " << operation << std::endl;
//     std::cout << "complex = " << complex << std::endl;

//     return 0;
// }