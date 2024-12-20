#include "bit_array_operations.hpp"


// Generates Alice's key
void fill_random_bits(XoshiroCpp::Xoshiro256PlusPlus &prng,
                      std::vector<int> &bit_array)
{
    std::uniform_int_distribution<int> distribution(0, 1);

    // Generate random bits and fill the vector
    for (int i = 0; i < bit_array.size(); ++i)
    {
        bit_array[i] = distribution(prng);
    }
}

// Generates Bob's key by making errors in Alice's key. Generates the exact number of errors in the key and returns the exact QBER.
void introduce_errors(XoshiroCpp::Xoshiro256PlusPlus &prng,
                      const std::vector<int> &bit_array,
                      float QBER,
                      std::vector<int> &output_bit_array_with_errors)
{
    size_t array_length = bit_array.size();
    size_t num_errors = static_cast<size_t>(array_length * QBER);
    
    output_bit_array_with_errors = bit_array;

    if (num_errors > 0)
    {
        std::vector<size_t> error_positions(array_length);
        for (size_t i = 0; i < array_length; ++i)
        {
            error_positions[i] = i;
        }

        std::shuffle(error_positions.begin(), error_positions.end(), prng);

        for (size_t i = 0; i < num_errors; ++i)
        {
            output_bit_array_with_errors[error_positions[i]] ^= 1;
        }
    }
}

// Shuffles Alice's and Bob's bit arrays by seed
void shuffle_array_bits(std::vector<int> &alice_bit_array,
                        std::vector<int> &bob_bit_array, 
                        size_t seed)
{
    XoshiroCpp::Xoshiro256PlusPlus rng1(seed);
    std::shuffle(alice_bit_array.begin(), alice_bit_array.end(), rng1);
    XoshiroCpp::Xoshiro256PlusPlus rng2(seed);
    std::shuffle(bob_bit_array.begin(), bob_bit_array.end(), rng2);
}

// Calculates an array that consists of 0 and 1, where 1 denotes that Alice's and Bob's bits are different
void calculate_error_positions(const std::vector<int> &alice_bit_array, 
                               const std::vector<int> &bob_bit_array,
                               std::vector<int> &output_error_positions)
{
    for (size_t i = 0; i < alice_bit_array.size(); ++i)
    {
        output_error_positions[i] = alice_bit_array[i] ^ bob_bit_array[i];
    }
}