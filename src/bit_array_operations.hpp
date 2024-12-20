#pragma once 
#include <random>
#include <algorithm>

#include <XoshiroCpp.hpp>

void fill_random_bits(XoshiroCpp::Xoshiro256PlusPlus &prng,
                      std::vector<int> &bit_array);
void introduce_errors(XoshiroCpp::Xoshiro256PlusPlus &prng,
                      const std::vector<int> &bit_array,
                      float QBER,
                      std::vector<int> &output_bit_array_with_errors);
void shuffle_array_bits(std::vector<int> &alice_bit_array,
                        std::vector<int> &bob_bit_array, 
                        size_t seed);
void calculate_error_positions(const std::vector<int> &alice_bit_array, 
                               const std::vector<int> &bob_bit_array,
                               std::vector<int> &output_error_positions);                      