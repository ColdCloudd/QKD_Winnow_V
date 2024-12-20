#include "utils.hpp"

void print_array(const std::vector<int> &bit_array, 
                 const size_t &block_length)
{
    for (size_t i = 0; i < bit_array.size(); i++)
    {
        if (i % block_length == 0 && i != 0)
        {
            fmt::print(" ");
        }
        fmt::print(fg(fmt::color::blue), "{}", bit_array[i]);
    }
}

