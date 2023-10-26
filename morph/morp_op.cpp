#include "tomlplusplus/toml.hpp"
#include <iostream>
#include <map>
#include "math_mt.h"
#include "vctr_cpu.h"
#include "fcns_image_cpu.h"
#include "H5Cpp.h"

class MorphOperation {
private:
    std::string input_file;
    std::string input_dataset;
    std::string output_file;
    std::string output_dataset;
    std::string operation;
    H5::DataType data_type;
    int nkr;

    H5::DataType read_data_type() const;

    template <class T>
    mt::Vctr_cpu<T> read_h5();

    template <class T>
    void write_h5(mt::Vctr_cpu<T> &mx_o);

    template <class T>
    void run();

public:
    MorphOperation(const std::string &config_file);

    void run_by_data_type();
};

H5::DataType MorphOperation::read_data_type() const {
    H5::H5File file(input_file, H5F_ACC_RDONLY);
    H5::DataSet dataset = file.openDataSet("image");
    return dataset.getDataType();
}

template <class T>
mt::Vctr_cpu<T> MorphOperation::read_h5() {
    H5::H5File file(input_file, H5F_ACC_RDONLY);
    H5::DataSet dataset = file.openDataSet(input_dataset);
    H5::DataSpace dataspace = dataset.getSpace();
    H5::DataType datatype = dataset.getDataType();

    hsize_t dims_out[2];
    dataspace.getSimpleExtentDims(dims_out, NULL);

    mt::Vctr_cpu<T> mx_i(dt_shape({dims_out[1], dims_out[0]}));
    dataset.read(mx_i.m_data, datatype);
    file.close();

    return mx_i;
}

template <class T>
void MorphOperation::write_h5(mt::Vctr_cpu<T> &mx_o) {
    H5::H5File file(output_file, H5F_ACC_TRUNC);
    hsize_t dimsf[2];
    dimsf[1] = mx_o.m_s0;
    dimsf[0] = mx_o.m_s1;
    H5::DataSpace dataspace(2, dimsf);
    H5::DataSet dataset = file.createDataSet(output_dataset, data_type, dataspace);

    dataset.write(mx_o.m_data, data_type);

    file.close();
}

template <class T>
void MorphOperation::run() {
    auto mx_i = read_h5<T>();
    mt::Vctr_cpu<T> mx_o;

    // I want to choose the function to call based on the operation string with a map
    if (this->operation == "dilate") {
        mx_o = mt::fcns_image_cpu::fcn_morp_op_dilate(mx_i, nkr);
    } else if (this->operation == "erode") {
        mx_o = mt::fcns_image_cpu::fcn_morp_op_erode(mx_i, nkr);
    } else if (this->operation == "open")
    {
        mx_o = mt::fcns_image_cpu::fcn_morp_op_open(mx_i, nkr);
    }
    else if (this->operation == "close")
    {
        mx_o = mt::fcns_image_cpu::fcn_morp_op_close(mx_i, nkr);
    }
    else if (this->operation == "tophat")
    {
        mx_o = mt::fcns_image_cpu::fcn_morp_op_tophat(mx_i, nkr);
    }
    else {
        throw std::runtime_error("Unsupported morphological operation");
    }
    
    write_h5<T>(mx_o);
}

MorphOperation::MorphOperation(const std::string &config_file) {
    auto config = toml::parse_file(config_file);
    input_file = config["Input"]["file"].as_string()->get();
    input_dataset = config["Input"]["dataset"].as_string()->get();
    output_file = config["Output"]["file"].as_string()->get();
    output_dataset = config["Output"]["dataset"].as_string()->get();
    nkr = config["Parameters"]["nkr"].as_integer()->get();
    operation = config["Parameters"]["operation"].as_string()->get();
    data_type = read_data_type();
}

void MorphOperation::run_by_data_type() 
{
    if (data_type == H5::PredType::NATIVE_FLOAT) {
        this->run<dt_float32>();
    } else if (data_type == H5::PredType::NATIVE_DOUBLE) {
        this->run<dt_float64>();
    } else if (data_type == H5::PredType::NATIVE_INT8) {
        this->run<int>();
    } else if (data_type == H5::PredType::NATIVE_UINT8) {
        this->run<dt_uint8>();
    } else if (data_type == H5::PredType::NATIVE_SHORT) {
        this->run<short>();
    } else if (data_type == H5::PredType::NATIVE_USHORT) {
        this->run<unsigned short>();
    } else if (data_type == H5::PredType::NATIVE_LLONG) {
        this->run<long long>();
    } else if (data_type == H5::PredType::NATIVE_ULLONG) {
        this->run<unsigned long long>();
    } else {
        throw std::runtime_error("Unsupported datatype in HDF5 file");
    }
}

int main(int argc, char** argv) {
    if (argc < 2) {
        std::cerr << "Error: Configuration file not provided!" << std::endl;
        return 1;
    }

    MorphOperation morphOp(argv[1]);
    morphOp.run_by_data_type();

    return 0;
}
