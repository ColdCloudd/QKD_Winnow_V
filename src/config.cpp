#include "config.hpp"

config_data get_config_data(fs::path config_path)
{
    if (!fs::exists(config_path))
    {
        throw std::runtime_error("Configuration file not found: " + config_path.string());
    }

    std::ifstream config_file(config_path);
    if (!config_file.is_open())
    {
        throw std::runtime_error("Failed to open configuration file: " + config_path.string());
    }

    json config = json::parse(config_file);
    config_file.close();
    if (config.empty())
    {
        throw std::runtime_error("Configuration file is empty: " + config_path.string());
    }
    
    try
    {
        config_data cfg{};
        cfg.THREADS_NUMBER = config["threads_number"].template get<size_t>();
        if (cfg.THREADS_NUMBER < 1)
        {
            throw std::runtime_error("Number of threads must be >= 1!!");
        }

        cfg.TRIALS_NUMBER = config["trials_number"].template get<size_t>();
        if (cfg.TRIALS_NUMBER < 1)
        {
            throw std::runtime_error("Number of trials must be >= 1!!");
        }

        if (config["use_config_simulation_seed"].template get<bool>())
        {
            cfg.SIMULATION_SEED = config["simulation_seed"].template get<size_t>();
        }
        else
        {
            cfg.SIMULATION_SEED = time(nullptr);
        }
        
        cfg.SHUFFLE_MODE = config["shuffle_mode"].template get<bool>();
        cfg.ENABLE_THROUGHPUT_MEASUREMENT = config["enable_throughput_measurement"].template get<bool>();
        if (cfg.ENABLE_THROUGHPUT_MEASUREMENT)
        {
            fmt::print(fg(fmt::color::purple), "WARNING: Throughput measurement is enabled. It is recommended to perform experiments in single-threaded mode.\n");

            cfg.CONSIDER_RTT = config["consider_RTT"].template get<bool>();
            if (cfg.CONSIDER_RTT)
                cfg.RTT = config["RTT"].template get<size_t>();
        }

        cfg.SIFTED_KEY_LENGTH = config["sifted_key_length"].template get<size_t>();
        if (cfg.SIFTED_KEY_LENGTH < 8)
        {
            throw std::runtime_error("Minimum sifted key length is 8!");
        }

        cfg.INITIAL_SYNDROME_LENGTH = config["initial_syndrome_length"].template get<size_t>();
        if (cfg.INITIAL_SYNDROME_LENGTH < 3)
        {
            throw std::runtime_error("Minimum initial syndrome length is 3!");
        }
        
        cfg.QBER = config["qber"].template get<std::vector<double>>();
        for (size_t i = 0; i < cfg.QBER.size(); i++)
        {
            if (cfg.QBER[i] <= 0. || cfg.QBER[i] >= 1.)
            {
                throw std::runtime_error("QBER must be: 0 < QBER < 1!");
            }
        }

        cfg.USE_SPECIFIED_COMBINATIONS = config["use_specified_combinations"].template get<bool>();
        if (cfg.USE_SPECIFIED_COMBINATIONS)
        {
            auto combinations = config["specified_combinations"];
            for (size_t i = 0; i < combinations.size(); i++)
            {
                if (combinations[i].size() != 2)
                {
                    throw std::runtime_error("Format error of a combination located in the 'specified_combinations' vector!");
                }
                std::vector<size_t> schedule = combinations[i][0].get<std::vector<size_t>>();
                std::vector<double> qbers = combinations[i][1].get<std::vector<double>>();
                if (schedule.empty() || qbers.empty())
                {
                    throw std::runtime_error("There should be no empty vectors in the 'specified_combinations' vector element!");
                }
                for (size_t i = 0; i < qbers.size(); i++)
                {
                    if (qbers[i] <= 0. || qbers[i] >= 1.)
                    {
                        throw std::runtime_error("QBER must be: 0 < QBER < 1!");
                    }
                }
                specified_combination spec_comb{};
                spec_comb.schedule = schedule;
                spec_comb.qbers = qbers;
                cfg.SPECIFIED_COMBINATIONS.push_back(spec_comb);
            }
        }
        else
        {
            cfg.SCHEDULE_ELEMENTS = config["schedule_elements"].template get<std::vector<std::vector<size_t>>>();
            if (cfg.SCHEDULE_ELEMENTS.empty())
            {
                throw std::runtime_error("Vector 'schedule_elements' cannot be empty!");
            }
            for (size_t i = 0; i < cfg.SCHEDULE_ELEMENTS.size(); i++)
            {
                if (cfg.SCHEDULE_ELEMENTS[i].empty())
                {
                    throw std::runtime_error("Vector 'schedule_elements' cannot be empty!");
                }
            }
        }
        
        cfg.TRACE_WINNOW = config["trace_winnow"].template get<bool>();
        return cfg;
    }
    catch(const std::exception& e)
    {
        fmt::print(stderr, fg(fmt::color::red),"An error occurred while reading a configuration parameter.\n");
        throw;
    }
}