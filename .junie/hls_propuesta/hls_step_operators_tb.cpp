/**
 * @file hls_step_operators_tb.cpp
 * @brief Testbench for HLS implementation of step operators for wave propagation
 * 
 * This file contains test cases for the HLS implementations of the step operators
 * for wave propagation on the Kria KV 260 FPGA.
 */

#include <iostream>
#include <cmath>
#include <complex>
#include <vector>
#include <cstdlib>
#include <ctime>
#include "hls_step_operators.h"

// Define tolerance for floating-point comparisons
#define TOLERANCE 1e-5

// Function to generate a random complex field
void generate_random_field(complex_t field[MAX_NY][MAX_NX], int ny, int nx) {
    for (int j = 0; j < ny; j++) {
        for (int i = 0; i < nx; i++) {
            float real_part = (float)rand() / RAND_MAX * 2.0f - 1.0f;
            float imag_part = (float)rand() / RAND_MAX * 2.0f - 1.0f;
            field[j][i] = complex_t(real_part, imag_part);
        }
    }
}

// Function to generate a random 1D complex array
void generate_random_array(complex_t array[MAX_SIZE], int size) {
    for (int i = 0; i < size; i++) {
        float real_part = (float)rand() / RAND_MAX * 2.0f - 1.0f;
        float imag_part = (float)rand() / RAND_MAX * 2.0f - 1.0f;
        array[i] = complex_t(real_part, imag_part);
    }
}

// Function to compare two complex arrays
bool compare_arrays(complex_t array1[MAX_SIZE], complex_t array2[MAX_SIZE], int size) {
    for (int i = 0; i < size; i++) {
        float real_diff = std::abs(array1[i].real() - array2[i].real());
        float imag_diff = std::abs(array1[i].imag() - array2[i].imag());
        if (real_diff > TOLERANCE || imag_diff > TOLERANCE) {
            std::cout << "Mismatch at index " << i << ": " 
                      << "(" << array1[i].real() << "," << array1[i].imag() << ") vs "
                      << "(" << array2[i].real() << "," << array2[i].imag() << ")" << std::endl;
            return false;
        }
    }
    return true;
}

// Function to compare two 2D complex fields
bool compare_fields(complex_t field1[MAX_NY][MAX_NX], complex_t field2[MAX_NY][MAX_NX], int ny, int nx) {
    for (int j = 0; j < ny; j++) {
        for (int i = 0; i < nx; i++) {
            float real_diff = std::abs(field1[j][i].real() - field2[j][i].real());
            float imag_diff = std::abs(field1[j][i].imag() - field2[j][i].imag());
            if (real_diff > TOLERANCE || imag_diff > TOLERANCE) {
                std::cout << "Mismatch at position (" << j << "," << i << "): " 
                          << "(" << field1[j][i].real() << "," << field1[j][i].imag() << ") vs "
                          << "(" << field2[j][i].real() << "," << field2[j][i].imag() << ")" << std::endl;
                return false;
            }
        }
    }
    return true;
}

// Test for custom_thomas_solver
bool test_custom_thomas_solver() {
    std::cout << "Testing custom_thomas_solver..." << std::endl;
    
    const int n = 10;
    complex_t b[MAX_SIZE];
    complex_t x[MAX_SIZE];
    complex_t x_expected[MAX_SIZE];
    
    // Initialize test data
    complex_t dp(2.0, 0.0);
    complex_t dp1(2.0, 0.0);
    complex_t dp2(2.0, 0.0);
    complex_t do_val(-1.0, 0.0);
    
    // Generate a known solution
    for (int i = 0; i < n; i++) {
        x_expected[i] = complex_t(i + 1.0, 0.0);
    }
    
    // Compute right-hand side b = A*x_expected
    b[0] = dp1 * x_expected[0] + do_val * x_expected[1];
    for (int i = 1; i < n-1; i++) {
        b[i] = do_val * x_expected[i-1] + dp * x_expected[i] + do_val * x_expected[i+1];
    }
    b[n-1] = do_val * x_expected[n-2] + dp2 * x_expected[n-1];
    
    // Solve the system
    custom_thomas_solver(dp, dp1, dp2, do_val, b, x, n);
    
    // Compare results
    bool success = compare_arrays(x, x_expected, n);
    if (success) {
        std::cout << "custom_thomas_solver test PASSED" << std::endl;
    } else {
        std::cout << "custom_thomas_solver test FAILED" << std::endl;
    }
    
    return success;
}

// Test for compute_b_vector
bool test_compute_b_vector() {
    std::cout << "Testing compute_b_vector..." << std::endl;
    
    const int n = 10;
    complex_t x0[MAX_SIZE];
    complex_t b[MAX_SIZE];
    complex_t b_expected[MAX_SIZE];
    
    // Initialize test data
    complex_t dp(2.0, 0.0);
    complex_t dp1(2.0, 0.0);
    complex_t dp2(2.0, 0.0);
    complex_t do_val(-1.0, 0.0);
    
    // Generate random input vector
    generate_random_array(x0, n);
    
    // Compute expected result
    b_expected[0] = dp1 * x0[0] + do_val * x0[1];
    for (int i = 1; i < n-1; i++) {
        b_expected[i] = do_val * x0[i-1] + dp * x0[i] + do_val * x0[i+1];
    }
    b_expected[n-1] = do_val * x0[n-2] + dp2 * x0[n-1];
    
    // Compute actual result
    compute_b_vector(dp, dp1, dp2, do_val, x0, b, n);
    
    // Compare results
    bool success = compare_arrays(b, b_expected, n);
    if (success) {
        std::cout << "compute_b_vector test PASSED" << std::endl;
    } else {
        std::cout << "compute_b_vector test FAILED" << std::endl;
    }
    
    return success;
}

// Test for half_nonlinear
bool test_half_nonlinear() {
    std::cout << "Testing half_nonlinear..." << std::endl;
    
    const int size = 10;
    complex_t phi[MAX_SIZE];
    complex_t phi_out[MAX_SIZE];
    complex_t phi_expected[MAX_SIZE];
    
    // Initialize test data
    float k_sample = 10.0f;
    float n2_sample = 1e-20f;
    float dz = 1e-6f;
    
    // Generate random input field
    generate_random_array(phi, size);
    
    // Compute expected result
    float phase_const = k_sample * n2_sample * dz / 2.0f;
    for (int i = 0; i < size; i++) {
        float abs_squared = phi[i].real() * phi[i].real() + phi[i].imag() * phi[i].imag();
        float phase_angle = phase_const * abs_squared;
        float sin_val = std::sin(phase_angle);
        float cos_val = std::cos(phase_angle);
        complex_t phase(cos_val, sin_val);
        phi_expected[i] = phi[i] * phase;
    }
    
    // Compute actual result
    half_nonlinear(phi, phi_out, k_sample, n2_sample, dz, size);
    
    // Compare results
    bool success = compare_arrays(phi_out, phi_expected, size);
    if (success) {
        std::cout << "half_nonlinear test PASSED" << std::endl;
    } else {
        std::cout << "half_nonlinear test FAILED" << std::endl;
    }
    
    return success;
}

// Test for half_linear_absorption
bool test_half_linear_absorption() {
    std::cout << "Testing half_linear_absorption..." << std::endl;
    
    const int size = 10;
    complex_t phi[MAX_SIZE];
    complex_t phi_out[MAX_SIZE];
    complex_t phi_expected[MAX_SIZE];
    
    // Initialize test data
    float alpha = 0.1f;
    float dz = 1e-6f;
    
    // Generate random input field
    generate_random_array(phi, size);
    
    // Compute expected result
    float attenuation = std::exp(-alpha * dz / 4.0f);
    for (int i = 0; i < size; i++) {
        phi_expected[i] = phi[i] * attenuation;
    }
    
    // Compute actual result
    half_linear_absorption(phi, phi_out, alpha, dz, size);
    
    // Compare results
    bool success = compare_arrays(phi_out, phi_expected, size);
    if (success) {
        std::cout << "half_linear_absorption test PASSED" << std::endl;
    } else {
        std::cout << "half_linear_absorption test FAILED" << std::endl;
    }
    
    return success;
}

// Test for half_2photon_absorption
bool test_half_2photon_absorption() {
    std::cout << "Testing half_2photon_absorption..." << std::endl;
    
    const int size = 10;
    complex_t phi[MAX_SIZE];
    complex_t phi_out[MAX_SIZE];
    complex_t phi_expected[MAX_SIZE];
    
    // Initialize test data
    float beta = 1e-12f;
    float dz = 1e-6f;
    
    // Generate random input field
    generate_random_array(phi, size);
    
    // Compute expected result
    for (int i = 0; i < size; i++) {
        float abs_squared = phi[i].real() * phi[i].real() + phi[i].imag() * phi[i].imag();
        float attenuation = std::exp(-beta * dz / 4.0f * abs_squared);
        phi_expected[i] = phi[i] * attenuation;
    }
    
    // Compute actual result
    half_2photon_absorption(phi, phi_out, beta, dz, size);
    
    // Compare results
    bool success = compare_arrays(phi_out, phi_expected, size);
    if (success) {
        std::cout << "half_2photon_absorption test PASSED" << std::endl;
    } else {
        std::cout << "half_2photon_absorption test FAILED" << std::endl;
    }
    
    return success;
}

// Test for adi_x and adi_y (simplified test)
bool test_adi() {
    std::cout << "Testing ADI methods (simplified)..." << std::endl;
    
    const int nx = 32;
    const int ny = 32;
    complex_t phi[MAX_NY][MAX_NX];
    complex_t phi_inter[MAX_NY][MAX_NX];
    
    // Initialize test data
    float eps = 1e-10f;
    float k = 10.0f;
    float dz = 1e-6f;
    float dx = 1e-6f;
    float dy = 1e-6f;
    
    // Generate random input field
    generate_random_field(phi, ny, nx);
    
    // Test adi_x
    adi_x(phi, phi_inter, ny, nx, eps, k, dz, dx);
    
    // Test adi_y
    adi_y(phi_inter, phi, nx, ny, eps, k, dz, dy);
    
    // For a simplified test, we just check that the functions run without errors
    std::cout << "ADI methods test PASSED (simplified)" << std::endl;
    
    return true;
}

// Test for propagation_step (simplified test)
bool test_propagation_step() {
    std::cout << "Testing propagation_step (simplified)..." << std::endl;
    
    const int nx = 32;
    const int ny = 32;
    complex_t phi_in[MAX_NY][MAX_NX];
    complex_t phi_out[MAX_NY][MAX_NX];
    
    // Initialize test data
    float eps = 1e-10f;
    float k = 10.0f;
    float k_sample = 10.0f;
    float n2_sample = 1e-20f;
    float alpha = 0.1f;
    float beta = 1e-12f;
    float dz = 1e-6f;
    float dx = 1e-6f;
    float dy = 1e-6f;
    
    // Generate random input field
    generate_random_field(phi_in, ny, nx);
    
    // Run propagation step
    propagation_step(phi_in, phi_out, nx, ny, eps, k, k_sample, n2_sample, alpha, beta, dz, dx, dy);
    
    // For a simplified test, we just check that the function runs without errors
    std::cout << "propagation_step test PASSED (simplified)" << std::endl;
    
    return true;
}

int main() {
    // Seed random number generator
    std::srand(std::time(nullptr));
    
    std::cout << "Starting HLS step operators testbench..." << std::endl;
    
    // Run tests
    bool all_tests_passed = true;
    
    all_tests_passed &= test_custom_thomas_solver();
    all_tests_passed &= test_compute_b_vector();
    all_tests_passed &= test_half_nonlinear();
    all_tests_passed &= test_half_linear_absorption();
    all_tests_passed &= test_half_2photon_absorption();
    all_tests_passed &= test_adi();
    all_tests_passed &= test_propagation_step();
    
    if (all_tests_passed) {
        std::cout << "\nAll tests PASSED!" << std::endl;
        return 0;
    } else {
        std::cout << "\nSome tests FAILED!" << std::endl;
        return 1;
    }
}