#include <iostream>

#include "config.hpp"
#include "utils.hpp"
#include "winnow_algorithm.hpp"

config_data CFG {};

int main()
{
    CFG.TRACE_WINNOW = true;

    const size_t ARR_LEN = 32;
    std::vector<int> alice_arr = 
    {
        0, 1, 0, 1, 1, 1, 0, 0,
        1, 1, 0, 1, 0, 1, 1, 1,
        0, 1, 0, 1, 1, 1, 0, 0,
        0, 1, 1, 1, 0, 1, 0, 1
    };
    std::vector<int> bob_arr = alice_arr;

    bob_arr[10] = bob_arr[10] ^ 1; // block #2 (one error)

    bob_arr[17] = bob_arr[17] ^ 1; // block #3 (two errors)
    bob_arr[23] = bob_arr[23] ^ 1;

    bob_arr[24] = bob_arr[24] ^ 1; // block #4 (three errors)
    bob_arr[25] = bob_arr[25] ^ 1;
    bob_arr[26] = bob_arr[26] ^ 1;

    size_t syndrome_length = 3;
    std::vector<std::vector<int>> hash_mat = construct_Hamming_hash_matrix(syndrome_length);

    winnow(alice_arr, bob_arr, syndrome_length, hash_mat);

}