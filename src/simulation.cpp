#include "simulation.hpp"

// Computes combinations of runs, based on given elements as the Cartesian product of vectors
std::vector<std::vector<size_t>> cartesian_product(std::vector<std::vector<size_t>> trial_elements)
{
    auto product = [](long long a, std::vector<size_t> &b)
    { return a * b.size(); };
    const long long combination_number = accumulate(trial_elements.begin(), trial_elements.end(), 1LL, product);
    std::vector<std::vector<size_t>> result(combination_number, std::vector<size_t>(trial_elements.size()));
    for (long long n = 0; n < combination_number; ++n)
    {
        std::lldiv_t q{n, 0};
        for (long long i = trial_elements.size() - 1; 0 <= i; --i)
        {
            q = std::div(q.quot, trial_elements[i].size());
            result[n][i] = trial_elements[i][q.rem];
        }
    }
    return result;
}

// Generates combinations of runs consisting of the run number, combinations that determine
// the number of Winnow runs with a given block size, and the QBER probability
std::vector<test_combination> prepare_combinations(const std::vector<std::vector<size_t>>& trial_combinations,
                                                   std::vector<double> bit_error_rates)
{
    std::vector<test_combination> combinations(trial_combinations.size() * bit_error_rates.size());
    size_t test_number = 0;
    for (size_t i = 0; i < trial_combinations.size(); i++)
    {
        for (size_t j = 0; j < bit_error_rates.size(); j++)
        {
            combinations[test_number].test_number = test_number;
            combinations[test_number].trial_combination = trial_combinations[i];
            combinations[test_number].QBER = bit_error_rates[j];
            test_number++;
        }
    }
    return combinations;
}

// Runs the Winnow algorithm sequentially several times with different block sizes
void run_trial(std::vector<int> &alice_bit_array,
               std::vector<int> &bob_bit_array,
               const std::vector<size_t> &trial_combination,
               size_t seed)
{
    bool shuffle_bits = CFG.SHUFFLE_MODE;
    size_t syndrome_length = CFG.INITIAL_SYNDROME_LENGTH;

    bool padding = false;
    bool last_block_with_error = false;
    size_t block_length = 0;
    size_t padding_length = 0;
    size_t last_incompl_block_length = 0;
    size_t initial_array_length = alice_bit_array.size();
    size_t curr_array_length = alice_bit_array.size();
   
    for (size_t i = 0; i < trial_combination.size(); ++i)
    {
        block_length = static_cast<size_t>(pow(2, syndrome_length));
        std::vector<std::vector<int>> hash_mat = construct_Hamming_hash_matrix(syndrome_length);
        for (size_t j = 0; j < trial_combination[i]; ++j)
        {
            last_incompl_block_length = curr_array_length % block_length; // Before each Winnow run, сalculates the length of the last block
            
            if(last_incompl_block_length > 0)   // Padding 
            {
                padding_length = block_length - last_incompl_block_length;
                if(curr_array_length + padding_length > initial_array_length)    // If padding extends beyond the initial array, the last block is cut off
                {
                    curr_array_length -= last_incompl_block_length;
                    alice_bit_array.resize(curr_array_length);
                    bob_bit_array.resize(curr_array_length);
                }
                else
                {
                    alice_bit_array.resize(curr_array_length + padding_length);
                    bob_bit_array.resize(curr_array_length + padding_length);
                    std::fill(alice_bit_array.begin() + curr_array_length, alice_bit_array.end(), 0);
                    std::fill(bob_bit_array.begin() + curr_array_length, bob_bit_array.end(), 0);
                    padding = true;
                }
            }

            if (CFG.TRACE_WINNOW)
                last_block_with_error = winnow_trace(alice_bit_array, bob_bit_array, syndrome_length, hash_mat);
            else
                last_block_with_error = winnow(alice_bit_array, bob_bit_array, syndrome_length, hash_mat);
            
            curr_array_length = alice_bit_array.size();
            if(padding)
            {
                if(last_block_with_error)
                {
                    for (size_t k = 0; k < syndrome_length; k++)
                    {
                        if (static_cast<size_t>(pow(2, k)) >= last_incompl_block_length)
                        {
                            padding_length--;       // Number of padding bits that were removed by the winnow process 
                        }
                    } 
                }
                curr_array_length -= padding_length;
                alice_bit_array.resize(curr_array_length);
                bob_bit_array.resize(curr_array_length);
                padding = false;
            }
            
            if (shuffle_bits)
            {
                shuffle_array_bits(alice_bit_array, bob_bit_array, seed);
                seed++;
            }

        }
        syndrome_length++;
    }
}

// Runs the experiment with the given TRIAL_NUMBER- times combination and calculates
// the final average key error rate and the final average key fraction
test_result run_test(const test_combination combination,
                     size_t seed)
{
    size_t errors_number = 0;
    std::vector<int> alice_bit_array(CFG.SIFTED_KEY_LENGTH);
    std::vector<int> bob_bit_array(CFG.SIFTED_KEY_LENGTH);
    std::vector<int> error_positions_array(CFG.SIFTED_KEY_LENGTH);

    double frame_error_rate = 0;
    std::vector<double> final_qbers(CFG.TRIALS_NUMBER); 
    std::vector<double> final_fractions(CFG.TRIALS_NUMBER); 

    // Pseudo-random number generator
    XoshiroCpp::Xoshiro256PlusPlus prng(seed);
    std::uniform_int_distribution<size_t> distribution(0, std::numeric_limits<size_t>::max());

    std::vector<std::chrono::microseconds> runtimes;
    if (CFG.ENABLE_THROUGHPUT_MEASUREMENT)
        runtimes.resize(CFG.TRIALS_NUMBER);

    for (size_t i = 0; i < CFG.TRIALS_NUMBER; ++i)
    {
        alice_bit_array.resize(CFG.SIFTED_KEY_LENGTH);
        bob_bit_array.resize(CFG.SIFTED_KEY_LENGTH);

        fill_random_bits(prng, alice_bit_array);
        introduce_errors(prng, alice_bit_array, combination.QBER, bob_bit_array);
        
        if (CFG.ENABLE_THROUGHPUT_MEASUREMENT)
        {
            auto start = std::chrono::high_resolution_clock::now();
            run_trial(alice_bit_array, bob_bit_array, combination.trial_combination, distribution(prng));
            auto end = std::chrono::high_resolution_clock::now();
            runtimes[i] = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
        }
        else
        {
            run_trial(alice_bit_array, bob_bit_array, combination.trial_combination, distribution(prng));
        }

        error_positions_array.resize(alice_bit_array.size());
        calculate_error_positions(alice_bit_array, bob_bit_array, error_positions_array);
        errors_number = std::accumulate(error_positions_array.begin(), error_positions_array.end(), 0);
        if(errors_number != 0)
            frame_error_rate += 1.;
        final_qbers[i] = static_cast<double>(errors_number) / static_cast<double>(alice_bit_array.size());
        final_fractions[i] = static_cast<double>(alice_bit_array.size()) / static_cast<double>(CFG.SIFTED_KEY_LENGTH);
    }
    test_result result;

    frame_error_rate /= static_cast<double>(CFG.TRIALS_NUMBER);

    double final_qber_mean = std::accumulate(final_qbers.begin(), final_qbers.end(), 0.) / static_cast<double>(CFG.TRIALS_NUMBER);
    auto final_qber_minmax = std::minmax_element(final_qbers.begin(), final_qbers.end());
    double final_qber_min = *final_qber_minmax.first;
    double final_qber_max = *final_qber_minmax.second;
    double final_qber_variance = std::accumulate(final_qbers.begin(), final_qbers.end(), 0., 
        [final_qber_mean](double acc, double val) {
            return acc + std::pow(val - final_qber_mean, 2);
        });
    double final_qber_std_dev = sqrt(final_qber_variance / static_cast<double>(CFG.TRIALS_NUMBER));

    double final_fraction_mean = std::accumulate(final_fractions.begin(), final_fractions.end(), 0.) / static_cast<double>(CFG.TRIALS_NUMBER);
    auto final_fraction_minmax = std::minmax_element(final_fractions.begin(), final_fractions.end());
    double final_fraction_min = *final_fraction_minmax.first;
    double final_fraction_max = *final_fraction_minmax.second;
    double final_fraction_variance = std::accumulate(final_fractions.begin(), final_fractions.end(), 0., 
        [final_fraction_mean](double acc, double val) {
            return acc + std::pow(val - final_fraction_mean, 2);
        });
    double final_fraction_std_dev = sqrt(final_fraction_variance / static_cast<double>(CFG.TRIALS_NUMBER));

    if (CFG.ENABLE_THROUGHPUT_MEASUREMENT)
    {
        double out_key_length{};
        const double MICROSECONDS_IN_SECOND = 1000000.;
        const double MICROSECONDS_IN_MILLISECOND = 1000.;
        double curr_throughput{};
        double throughput_max = 0;
        double throughput_min = std::numeric_limits<double>::max();
        double throughput_mean = 0;
        double throughput_std_dev = 0;
        double message_count = 2. * static_cast<double>(std::accumulate(
            combination.trial_combination.begin(), combination.trial_combination.end(), 0));    // Number of network interactions between Alice and Bob 

        for (size_t i = 0; i < runtimes.size(); ++i)
        {
            out_key_length = final_fractions[i] * static_cast<double>(CFG.SIFTED_KEY_LENGTH);
            if (CFG.CONSIDER_RTT)
            {
                curr_throughput = out_key_length * MICROSECONDS_IN_SECOND / 
                (static_cast<double>(runtimes[i].count()) + static_cast<double>(CFG.RTT) * MICROSECONDS_IN_MILLISECOND * message_count);    // bits/s
            }
            else
            {
                curr_throughput = out_key_length * MICROSECONDS_IN_SECOND / static_cast<double>(runtimes[i].count());   // bits/s
            }

            throughput_mean += curr_throughput;
            if (curr_throughput > throughput_max)
            {
                throughput_max = curr_throughput;
            }
            if (curr_throughput < throughput_min)
            {
                throughput_min = curr_throughput;
            }
        }
        throughput_mean /= static_cast<double>(CFG.TRIALS_NUMBER);
        
        for (size_t i = 0; i < runtimes.size(); ++i)
        {
            out_key_length = final_fractions[i] * static_cast<double>(CFG.SIFTED_KEY_LENGTH);
            if (CFG.CONSIDER_RTT)
            {
                curr_throughput = out_key_length * MICROSECONDS_IN_SECOND / 
                (static_cast<double>(runtimes[i].count()) + static_cast<double>(CFG.RTT) * MICROSECONDS_IN_MILLISECOND * message_count);    // bits/s
            }
            else
            {
                curr_throughput = out_key_length * MICROSECONDS_IN_SECOND / static_cast<double>(runtimes[i].count());   // bits/s
            }
            throughput_std_dev += pow((curr_throughput - throughput_mean), 2);
        }
        throughput_std_dev /= static_cast<double>(CFG.TRIALS_NUMBER);
        throughput_std_dev = sqrt(throughput_std_dev);
        
        result.throughput_max = static_cast<size_t>(throughput_max);
        result.throughput_min = static_cast<size_t>(throughput_min);
        result.throughput_mean = static_cast<size_t>(throughput_mean);
        result.throughput_std_dev = static_cast<size_t>(throughput_std_dev);
    }
    
    result.test_number = combination.test_number;
    result.trial_combination = combination.trial_combination;
    result.initial_qber = static_cast<double>(static_cast<size_t>(static_cast<float>(CFG.SIFTED_KEY_LENGTH) 
    * combination.QBER)) / static_cast<double>(CFG.SIFTED_KEY_LENGTH); // Exact QBER in the key.

    result.frame_error_rate = frame_error_rate;

    result.final_qber_mean = final_qber_mean;
    result.final_qber_std_dev = final_qber_std_dev;
    result.final_qber_min = final_qber_min;
    result.final_qber_max = final_qber_max;

    result.final_fraction_mean = final_fraction_mean;
    result.final_fraction_std_dev = final_fraction_std_dev;
    result.final_fraction_min = final_fraction_min;
    result.final_fraction_max = final_fraction_max;

    return result;
}

// Distributes all combinations of the experiment evenly across the CPU threads and runs it
std::vector<test_result> run_simulation(const std::vector<test_combination> &combinations)
{
    using namespace indicators;

    std::vector<test_result> results(combinations.size());
    BS::thread_pool pool(CFG.THREADS_NUMBER);

    indicators::show_console_cursor(false);
    indicators::ProgressBar bar{
        option::BarWidth{50},
        option::Start{" ["},
        option::Fill{"="},
        option::Lead{">"},
        option::Remainder{"-"},
        option::End{"]"},
        option::PrefixText{"PROGRESS"},
        option::ForegroundColor{Color::green},
        option::ShowElapsedTime{true},
        option::ShowRemainingTime{true},
        option::FontStyles{std::vector<FontStyle>{FontStyle::bold}},
        option::MaxProgress{combinations.size()}};

    size_t iteration = 0;
    std::mt19937 prng(CFG.SIMULATION_SEED);
    std::uniform_int_distribution<size_t> distribution(0, std::numeric_limits<size_t>::max());
    std::vector<size_t> seeds(combinations.size());
    for (size_t i = 0; i < seeds.size(); i++)
    {
        seeds[i] = distribution(prng);
    }
    
    pool.detach_loop<size_t>(0, combinations.size(),
                             [&combinations, &results, &prng, &seeds, &bar, &iteration](size_t i)
                             {
                                 bar.set_option(option::PostfixText{
                                     std::to_string(iteration) + "/" + std::to_string(combinations.size())});
                                 bar.tick();
                                 iteration++;

                                 results[i] = run_test(combinations[i], seeds[i]);
                             });
    pool.wait();
    indicators::show_console_cursor(true);

    return results;
}

// Returns the combination as a python tuple in string format
std::string get_trial_combination_string(const std::vector<size_t> &combination)
{
    std::string comb_str = "(";
    for (size_t i = 0; i < combination.size(); i++)
    {
        comb_str += std::to_string(combination[i]);
        if (i < combination.size() - 1)
        {
            comb_str += ", ";
        }
    }
    comb_str += ")";
    return comb_str;
}

// Returns the run represented by the sequence (example: 0;2;1 - 0 passes with block size 8, then 2 passes with block size 16 and 1 pass with block size 32).
std::string get_num_pass_with_block_size_sequence_string(const std::vector<size_t> &combination)
{
    std::string seq_str = "";
    for (size_t i = 0; i < combination.size(); i++)
    {
        seq_str += std::to_string(combination[i]);
        if (i < combination.size() - 1)
        {
            seq_str += ";";
        }
    }
    return seq_str;
}

// Returns a block size header (example: N=8;N=16;N=32 ...).
std::string get_header_block_size_string(const std::vector<size_t> &combination)
{
    std::string header_str = "";
    size_t syndrome_length = CFG.INITIAL_SYNDROME_LENGTH;
    for (size_t i = 0; i < combination.size(); i++)
    {
        header_str += "N=";
        header_str += std::to_string(static_cast<size_t>(pow(2, syndrome_length)));
        if (i < combination.size() - 1)
        {
            header_str += ";";
        }
        syndrome_length++;
    }
    return header_str;
}

// Records the results of the simulation in a ".csv" format file
void write_file(const std::vector<test_result> &data,
                fs::path directory)
{
    try
    {
        if (!fs::exists(directory))
        {
            fs::create_directories(directory);
        }
        std::string base_filename = "winnow(trial_num=" + std::to_string(CFG.TRIALS_NUMBER) + ",shuff_mode=" + ((CFG.SHUFFLE_MODE)?"on":"off") 
        +  ",key_length=" + std::to_string(CFG.SIFTED_KEY_LENGTH) + ((CFG.ENABLE_THROUGHPUT_MEASUREMENT && CFG.CONSIDER_RTT) ? 
        (",RTT=" + std::to_string(CFG.RTT)):"") + ",seed=" + std::to_string(CFG.SIMULATION_SEED) + ")";
        std::string extension = ".csv";
        fs::path result_file_path = directory / (base_filename + extension);

        size_t file_count = 1;
        while (fs::exists(result_file_path))
        {
            result_file_path = directory / (base_filename + "_" + std::to_string(file_count) + extension);
            file_count++;
        }

        std::fstream fout;
        fout.open(result_file_path, std::ios::out | std::ios::trunc);
        fout << "№;TRIAL_COMBINATION;"<< get_header_block_size_string(data[0].trial_combination) << ";INITIAL_QBER;FINAL_QBER_MEAN;FINAL_QBER_STD_DEV;FINAL_QBER_MIN;FINAL_QBER_MAX;" 
             << "FINAL_FRACTION_MEAN;FINAL_FRACTION_STD_DEV;FINAL_FRACTION_MIN;FINAL_FRACTION_MAX;FER" 
             << ((CFG.ENABLE_THROUGHPUT_MEASUREMENT) ? ";THROUGHPUT_MEAN;THROUGHPUT_STD_DEV;THROUGHPUT_MIN;THROUGHPUT_MAX" : "") <<"\n";
        for (size_t i = 0; i < data.size(); i++)
        {
            fout << data[i].test_number << ";" << get_trial_combination_string(data[i].trial_combination) << ";"<< get_num_pass_with_block_size_sequence_string(data[i].trial_combination) << ";"
                 << data[i].initial_qber << ";" << data[i].final_qber_mean << ";" << data[i].final_qber_std_dev << ";" << data[i].final_qber_min << ";" << data[i].final_qber_max << ";"
                 << data[i].final_fraction_mean << ";" << data[i].final_fraction_std_dev << ";" << data[i].final_fraction_min << ";" << data[i].final_fraction_max << ";" << data[i].frame_error_rate 
                 << ((CFG.ENABLE_THROUGHPUT_MEASUREMENT) ? (";" + std::to_string(data[i].throughput_mean) + ";" + std::to_string(data[i].throughput_std_dev) + ";" + std::to_string(data[i].throughput_min) + ";"
                 + std::to_string(data[i].throughput_max)):"") << "\n";
        }
        fout.close();
    }
    catch (const std::exception &e)
    {
        fmt::print(stderr, fg(fmt::color::red),"An error occurred while writing to the file.\n");
        throw;
    }
}