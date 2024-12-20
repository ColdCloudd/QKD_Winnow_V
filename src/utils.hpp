#pragma once
#include <vector>
#include <string>
#include <fstream>
#include <filesystem>

#include <fmt/core.h>
#include <fmt/color.h>
#include <fmt/ranges.h>

namespace fs = std::filesystem;

void print_array(const std::vector<int> &bit_array,
                 const size_t &block_length);

