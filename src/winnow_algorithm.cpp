#include "winnow_algorithm.hpp"

// Discards bits at the first position in the block for privacy amplification
void discard_bits_for_parity_check(const std::vector<int> &source_bit_array,
                                   std::vector<int> &destination_bit_array, 
                                   const size_t &syndrome_length)
{
    size_t source_block_size = static_cast<size_t>(pow(2, syndrome_length));
    size_t destination_block_size = source_block_size - 1;
    for (size_t i = 0, j = 0; i < source_bit_array.size(); i += source_block_size, j += destination_block_size)
    {
        std::copy(source_bit_array.begin() + i + 1, source_bit_array.begin() + i + source_block_size, destination_bit_array.begin() + j);
    }
}

// Discards bits at positions 2^n, where n=(0,1, ... syndrome_length - 1) in the block for privacy amplification
void discard_bits_for_syndrome(std::vector<int>::const_iterator source_bit_block, 
                               std::vector<int>::iterator destination_bit_block,
                               const std::vector<int> &discarded_bit_positions)
{
    size_t destination_current_start = 0;
    // Copying bits that are between the positions described in discarded_bit_positions
    for (size_t i = 1; i < discarded_bit_positions.size() - 1; ++i)
    {
        destination_current_start += discarded_bit_positions[i] - discarded_bit_positions[i - 1] - 1;
        std::copy(source_bit_block + discarded_bit_positions[i] + 1,
                  source_bit_block + discarded_bit_positions[i + 1],
                  destination_bit_block + destination_current_start);
    }
}

// Calculates the parity of a bit block
bool calculate_parity(std::vector<int>::const_iterator bit_block, 
                      const size_t &block_length)
{
    return std::accumulate(bit_block, bit_block + block_length, 0) % 2 == 0;
}

// Computes a matrix based on the Hamming hash function
std::vector<std::vector<int>> construct_Hamming_hash_matrix(const size_t &syndrome_length)
{
    size_t block_length = static_cast<size_t>(pow(2, syndrome_length)) - 1;

    std::vector<std::vector<int>> hash_matrix(syndrome_length, std::vector<int>(block_length));

    for (size_t i = 0; i < syndrome_length; ++i)
    {
        for (size_t j = 0; j < block_length; ++j)
        {
            hash_matrix[i][j] = static_cast<int>(floor((j + 1) / pow(2, i))) % 2;
        }
    }
    return hash_matrix;
}

// Calculates the syndrome by multiplying the Hamming matrix by the bit block column vector
void calculate_syndrome(std::vector<int>::const_iterator bit_block,
                        const size_t &block_length,
                        const std::vector<std::vector<int>> &hash_matrix, 
                        std::vector<int> &output_syndrome)
{
    int xor_sum = 0;
    size_t syndrome_length = hash_matrix.size();

    for (size_t i = 0; i < syndrome_length; ++i)
    {
        xor_sum = 0;
        for (size_t j = 0; j < block_length; ++j)
        {
            if (*(bit_block + j))
            {
                xor_sum ^= hash_matrix[i][j];
            }
        }
        output_syndrome[i] = xor_sum;
    }
}

// Calculates the error position in a block based on Alice's and Bob's syndromes and inverts this bit in Alice's block
void correct_error(std::vector<int>::iterator bit_block, 
                   const std::vector<int> &first_syndrome, 
                   const std::vector<int> &second_syndrome)
{
    int error_bit_position = -1;
    for (size_t i = 0; i < first_syndrome.size(); ++i)
    {
        error_bit_position += (first_syndrome[i] ^ second_syndrome[i]) * static_cast<int>(pow(2, i));
    }
    if (error_bit_position >= 0)
    {
        *(bit_block + error_bit_position) = !*(bit_block + error_bit_position);
    }
}

bool winnow_trace(std::vector<int> &alice_bit_array,
                  std::vector<int> &bob_bit_array, 
                  const size_t &syndrome_length,
                  const std::vector<std::vector<int>> &hash_mat)
{
    size_t block_len = static_cast<size_t>(pow(2, syndrome_length));
    size_t array_length = alice_bit_array.size();
    size_t blocks_cnt = array_length / block_len;

    std::vector<int> diff_par_blocks; // Contains the numbers of blocks whose parity bits did not match for Alice and Bob
    diff_par_blocks.reserve(blocks_cnt);
    for (size_t i = 0; i < array_length; i += block_len)
    {
        if (calculate_parity(alice_bit_array.begin() + i, block_len) != calculate_parity(bob_bit_array.begin() + i, block_len))
        {
            diff_par_blocks.push_back(static_cast<int>(i / block_len));
        }
    }

    fmt::print(fg(fmt::color::blue), "______________WINNOW_TRACE______________\n");
    fmt::print(fg(fmt::color::blue), "{}      - Numbers of blocks with detected errors\n", fmt::join(diff_par_blocks, " "));
    
    size_t priv_amp_arr_len = array_length - static_cast<size_t>(array_length / block_len);
    std::vector<int> alice_priv_amp(priv_amp_arr_len);
    std::vector<int> bob_priv_amp(priv_amp_arr_len);

    // Privacy amplification by discarding the first bit in each block.
    discard_bits_for_parity_check(alice_bit_array, alice_priv_amp, syndrome_length);
    discard_bits_for_parity_check(bob_bit_array, bob_priv_amp, syndrome_length);

    print_array(alice_bit_array, block_len);
    fmt::print(fg(fmt::color::blue), "      - Alice's key\n");
    print_array(bob_bit_array, block_len);
    fmt::print(fg(fmt::color::blue), "      - Bob's key\n");

    print_array(alice_priv_amp, block_len - 1);
    fmt::print(fg(fmt::color::blue), "      - Alice's key after first privacy maintenance\n");
    print_array(bob_priv_amp, block_len - 1);
    fmt::print(fg(fmt::color::blue), "      - Bob's key after first privacy maintenance\n");     

    block_len -= 1;
    // Alice and Bob syndromes
    std::vector<int> alice_syn(syndrome_length);
    std::vector<int> bob_syn(syndrome_length);
    for (size_t i = 0; i < diff_par_blocks.size(); i++)
    {
        // Calculation of syndromes for blocks with non-matching parity bits, followed by error correction
        calculate_syndrome(alice_priv_amp.begin() + (diff_par_blocks[i] * block_len), block_len, hash_mat, alice_syn);
        calculate_syndrome(bob_priv_amp.begin() + (diff_par_blocks[i] * block_len), block_len, hash_mat, bob_syn);
        correct_error(alice_priv_amp.begin() + (diff_par_blocks[i] * block_len), alice_syn, bob_syn);
    }

    print_array(alice_priv_amp, block_len);
    fmt::print(fg(fmt::color::blue), "      - Alice's key after error correction\n");
    print_array(bob_priv_amp, block_len);
    fmt::print(fg(fmt::color::blue), "      - Bob's key after error correction\n");

    size_t remain_bits_cnt = block_len - syndrome_length; // Number of remaining bits in blocks for which syndromes were calculated
    size_t out_arr_len = priv_amp_arr_len - diff_par_blocks.size() * syndrome_length;
    alice_bit_array.resize(out_arr_len);
    bob_bit_array.resize(out_arr_len);

    // Contains bounds that specify valid bits [from 2^0, 2^1, ... , 2^(syndrome_length-1)],
    // and the last element (2^syndrome_length) which is the right boundary
    std::vector<int> disc_bit_pos(syndrome_length + 1);
    for (size_t i = 0; i < disc_bit_pos.size(); i++)
    {
        disc_bit_pos[i] = static_cast<int>(pow(2, i) - 1);
    }

    size_t j = 0;            // Used to move through diff_par_blocks
    size_t copy_delta = 0;   // Specifies the number of bits that can be copied
    size_t dest_cur_pos = 0; // The current position in the destination array from which new bits can be inserted
    size_t diff_par_len = diff_par_blocks.size();
    for (size_t i = 0; i < priv_amp_arr_len;)
    {
        if (j < diff_par_len && (i / block_len == diff_par_blocks[j])) // In a block whose number is in diff_par_blocks, bits are discarded
        {
            discard_bits_for_syndrome(alice_priv_amp.begin() + i, alice_bit_array.begin() + dest_cur_pos, disc_bit_pos);
            discard_bits_for_syndrome(bob_priv_amp.begin() + i, bob_bit_array.begin() + dest_cur_pos, disc_bit_pos);
            dest_cur_pos += remain_bits_cnt;
            i += block_len;
            j++;
        }
        else if (diff_par_len == 0) // No errors in arrays -> copy the array completely
        {
            std::copy(alice_priv_amp.begin(), alice_priv_amp.end(), alice_bit_array.begin());
            std::copy(bob_priv_amp.begin(), bob_priv_amp.end(), bob_bit_array.begin());
            i += priv_amp_arr_len;
        }
        else if (j >= diff_par_len) // Remaining blocks with no errors are copied
        {
            copy_delta = priv_amp_arr_len - i;
            std::copy(alice_priv_amp.begin() + i, alice_priv_amp.begin() + i + copy_delta, alice_bit_array.begin() + dest_cur_pos);
            std::copy(bob_priv_amp.begin() + i, bob_priv_amp.begin() + i + copy_delta, bob_bit_array.begin() + dest_cur_pos);
            i += copy_delta;
        }
        else
        {
            copy_delta = diff_par_blocks[j] * block_len - i; // Blocks between two erroneous blocks are copied
            std::copy(alice_priv_amp.begin() + i, alice_priv_amp.begin() + i + copy_delta, alice_bit_array.begin() + dest_cur_pos);
            std::copy(bob_priv_amp.begin() + i, bob_priv_amp.begin() + i + copy_delta, bob_bit_array.begin() + dest_cur_pos);
            dest_cur_pos += copy_delta;
            i += copy_delta;
        }
    }

    // Used to account for the number of bits to be removed for padding
    bool last_block_with_error = false; 
    if(diff_par_len > 0 && (diff_par_blocks.back() == (blocks_cnt - 1)))
    {
        last_block_with_error = true;
    }

    print_array(alice_bit_array, block_len);
    fmt::print(fg(fmt::color::blue), "      - Alice's key after second privacy maintenance\n");
    print_array(bob_bit_array, block_len);
    fmt::print(fg(fmt::color::blue), "      - Bob's key after second privacy maintenance\n");
    fmt::print(fg(fmt::color::blue), "______________WINNOW_TRACE______________\n");

    return last_block_with_error;
}

bool winnow(std::vector<int> &alice_bit_array,
            std::vector<int> &bob_bit_array, 
            const size_t &syndrome_length,
            const std::vector<std::vector<int>> &hash_mat)
{
    size_t init_block_len = static_cast<size_t>(pow(2, static_cast<double>(syndrome_length)));
    size_t block_len = init_block_len - 1;              // Block length after removing one bit to maintain privacy
    size_t init_arr_len = alice_bit_array.size();
    
    std::vector<int> alice_bit_arr = alice_bit_array;
    std::vector<int> bob_bit_arr = bob_bit_array;

    // Contains bounds that specify valid bits [from 2^0, 2^1, ... , 2^(syndrome_length-1)],
    // and the last index (2^syndrome_length) which is the right boundary
    std::vector<int> disc_bit_pos(syndrome_length + 1);
    for (size_t i = 0; i < disc_bit_pos.size(); ++i)
    {
        disc_bit_pos[i] = static_cast<int>(pow(2, i) - 1);
    }

    size_t dest_curr_index = 0;        // Is used to keep track of the current starting position for writing to the destination array
    size_t synd_msg_cnt = 0;         // Number of calculated syndromes based on the results of the pass 
    size_t remain_bits_cnt = block_len - syndrome_length;       // The number of bits remaining in the block after deletion to maintain privacy (as a result of the syndrome calculation)

    std::vector<int> alice_syn(syndrome_length);
    std::vector<int> bob_syn(syndrome_length);

    size_t curr_block_begin = 0;
    bool last_block_with_error = false;    // Used to account for the number of bits to be removed for padding
    size_t last_block_index = init_arr_len - init_block_len;    // Index of first element of the last block

    std::vector<int>::iterator alice_curr_block_begin {};
    std::vector<int>::iterator bob_curr_block_begin {};
    for (size_t i = 0; i < init_arr_len; i += init_block_len)
    {
        curr_block_begin = i + 1;   // !!! Starting with the second bit in the block, since the first one is removed as a result of maintaining privacy !!!
        alice_curr_block_begin = (alice_bit_arr.begin() + curr_block_begin);     
        bob_curr_block_begin = (bob_bit_arr.begin() + curr_block_begin);
        if (calculate_parity(alice_bit_arr.begin() + i, init_block_len) == calculate_parity(bob_bit_arr.begin() + i, init_block_len))   // If the parity matches, remove the first bit from each block
        {
            std::copy(alice_curr_block_begin, alice_bit_arr.begin() + (i + init_block_len), alice_bit_array.begin() + dest_curr_index);
            std::copy(bob_curr_block_begin, bob_bit_arr.begin() + (i + init_block_len), bob_bit_array.begin() + dest_curr_index);
            dest_curr_index += block_len;
        }
        else    // Remove the first bit from each block + `syndrome_length`-bits 
        {
            calculate_syndrome(alice_curr_block_begin, block_len, hash_mat, alice_syn);
            calculate_syndrome(bob_curr_block_begin, block_len, hash_mat, bob_syn);
            correct_error(alice_curr_block_begin, alice_syn, bob_syn);
            discard_bits_for_syndrome(alice_curr_block_begin, alice_bit_array.begin() + dest_curr_index, disc_bit_pos);
            discard_bits_for_syndrome(bob_curr_block_begin, bob_bit_array.begin() + dest_curr_index, disc_bit_pos);
            dest_curr_index += remain_bits_cnt;
            ++synd_msg_cnt;
            if (i == last_block_index)
            {
                last_block_with_error=true;
            }
        }
    }

    size_t blocks_cnt = init_arr_len / init_block_len;
    size_t arr_len = init_arr_len - blocks_cnt - synd_msg_cnt * syndrome_length;
    alice_bit_array.resize(arr_len);
    bob_bit_array.resize(arr_len);

    return last_block_with_error;
}