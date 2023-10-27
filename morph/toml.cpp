#include <iostream>
#include <fstream>
#include <map>
#include <string>
#include <stdexcept>
#include <algorithm>
#include <sstream>
#include "math_mt.h"

namespace toml
{
    // Function to trim whitespace from both ends of a string
    std::string trim(const std::string& str) {
        std::string trimmed = str;
        trimmed.erase(trimmed.begin(), std::find_if(trimmed.begin(), trimmed.end(), [](int ch) {
            return !std::isspace(ch);
        }));
        trimmed.erase(std::find_if(trimmed.rbegin(), trimmed.rend(), [](int ch) {
            return !std::isspace(ch);
        }).base(), trimmed.end());
        return trimmed;
    }

    template <typename T>
    T convertToType(const std::string& str) {
        std::istringstream iss(str);
        T value;
        iss >> value;
        return value;
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


    class TomlParser {
    public:
        TomlParser(const std::string& filename) : filename_(filename) {}

        bool parse() {
            std::ifstream file(filename_);

            if (!file) {
                std::cerr << "Error opening file: " << filename_ << std::endl;
                return false;
            }

            std::string line;
            std::string currentSection;
            while (std::getline(file, line)) {
                if (!line.empty() && line[0] == '[' && line[line.length() - 1] == ']') {
                    // This line is a section header
                    currentSection = line.substr(1, line.length() - 2);
                } else {
                    // This line is a key-value pair
                    size_t equalsPos = line.find('=');
                    if (equalsPos != std::string::npos) {
                        std::string key = line.substr(0, equalsPos);
                        std::string value = line.substr(equalsPos + 1);
                        key = trim(key); // Remove leading and trailing spaces
                        tomlData_[currentSection][key] = value;
                    }
                }
            }

            return true;
        }

        // Function to retrieve values as a specified type
        template <typename T>
        T get(const std::string& section, const std::string& key) {
            std::string valueStr = tomlData_[section][key];
            try {
                return convertToType<T>(valueStr);
            } catch (const std::exception& e) {
                throw std::runtime_error("Error converting value: " + valueStr + " to type " + typeid(T).name());
            }
        }

    private:
        std::string filename_;
        std::map<std::string, std::map<std::string, std::string>> tomlData_;

        // Helper function to convert a string to a specified type
    };
}

int main() {
    toml::TomlParser parser("input.toml");

    if (parser.parse()) {
        auto nkr = parser.get<int>("Parameters", "nkr");
        auto floatValue = parser.get<float>("Parameters", "floatValue");
        auto operation = parser.get<std::string>("Parameters", "operation");
        auto complex = parser.get<dt_cfloat32>("Parameters", "complex");
        std::cout << "nkr = " << nkr << std::endl;
        std::cout << "floatValue = " << floatValue << std::endl;
        std::cout << "operation = " << operation << std::endl;
        std::cout << "complex = " << complex << std::endl;
    }

    return 0;
}