#include <gtest/gtest.h>
#include "vctr_cpu.h"  // Ensure this includes the path to your vctr_cpu.h file

class VctrTest : public ::testing::Test {
protected:
    void SetUp() override {
        // Common setup code for all tests, if any.
    }

    void TearDown() override {
        // Common cleanup code for all tests, if any.
    }
};

// Test for the default constructor
TEST_F(VctrTest, DefaultConstructor) {
    mt::Vctr<double, mt::edev_cpu> vctr;  // Replace 'double' and 'edev_cpu' with appropriate types and device as per your implementation

    // Example Assertions
    // Check if the vector is empty after default construction
    EXPECT_TRUE(vctr.empty()); // Assuming 'empty' is a method that checks if the vector is empty

    // Check if the size of the vector is 0 after default construction
    EXPECT_EQ(0, vctr.size()); // Assuming 'size' is a method that returns the size of the vector
}

// Add additional tests for other constructors and methods...

int main(int argc, char **argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
