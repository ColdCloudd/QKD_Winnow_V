#pragma once
#include <random>
#include <vector>
#include <numeric>
#include <cstdlib>

#include <BS_thread_pool.hpp>
#include <indicators/progress_bar.hpp>
#include <indicators/cursor_control.hpp>

#include "config.hpp"
#include "bit_array_operations.hpp"
#include "winnow_algorithm.hpp"

struct test_combination
{
    int test_number{};
    std::vector<size_t> trial_combination{};
    float QBER{};
};

struct test_result
{
    int test_number{};
    std::vector<size_t> trial_combination{};
    double initial_qber{};      // Exact QBER in the key.

    double frame_error_rate{};
    
    double final_qber_mean{};
    double final_qber_std_dev{};
    double final_qber_min{};
    double final_qber_max{};

    double final_fraction_mean{};
    double final_fraction_std_dev{};
    double final_fraction_min{};
    double final_fraction_max{};
};

std::vector<std::vector<size_t>> cartesian_product(std::vector<std::vector<size_t>> trial_elements);
std::vector<test_combination> prepare_combinations(const std::vector<std::vector<size_t>>& trial_combinations,
                                                   std::vector<double> bit_error_rates);
void run_trial(std::vector<int> &alice_bit_array,
                 std::vector<int> &bob_bit_array,
                 const std::vector<size_t> &trial_combination,
                 size_t seed);
test_result run_test(const test_combination combination,
                     size_t seed);
std::vector<test_result> run_simulation(const std::vector<test_combination> &combinations);
std::string get_trial_combination_string(const std::vector<size_t> &combination);
std::string get_num_pass_with_block_size_sequence_string(const std::vector<size_t> &combination);
std::string get_header_block_size_string(const std::vector<size_t> &combination);
void write_file(const std::vector<test_result> &data, 
                fs::path directory);
