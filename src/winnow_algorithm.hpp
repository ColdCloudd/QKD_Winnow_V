#pragma once
#include <cmath>
#include <vector>
#include <numeric>
#include <algorithm>

#include <fmt/core.h>
#include <fmt/color.h>
#include <fmt/ranges.h>

#include "config.hpp"
#include "utils.hpp"

void discard_bits_for_parity_check(const std::vector<int> &source_bit_array,
                                   std::vector<int> &destination_bit_array, 
                                   const size_t &syndrome_length);
void discard_bits_for_syndrome(std::vector<int>::const_iterator source_bit_block, 
                               std::vector<int>::iterator destination_bit_block,
                               const std::vector<int> &discarded_bit_positions);
bool calculate_block_parity(std::vector<int>::const_iterator bit_block, 
                            const size_t &block_length);
std::vector<std::vector<int>> construct_Hamming_hash_matrix(const size_t &syndrome_length);
void calculate_syndrome(std::vector<int>::const_iterator bit_block,
                        const size_t &block_length,
                        const std::vector<std::vector<int>> &hash_matrix, 
                        std::vector<int> &output_syndrome);                               
void correct_error(std::vector<int>::iterator bit_block, 
                   const std::vector<int> &first_syndrome, 
                   const std::vector<int> &second_syndrome);
bool winnow(std::vector<int> &alice_bit_array,
            std::vector<int> &bob_bit_array, 
            const size_t &syndrome_length,
            const std::vector<std::vector<int>> &hash_mat);                

