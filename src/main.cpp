#include "config.hpp"
#include "utils.hpp"
#include "simulation.hpp"

config_data CFG {};

#ifdef USE_CURRENT_DIR
    const fs::path CONFIG_PATH = std::filesystem::current_path() / "config.json";
    const fs::path RESULTS_DIRECTORY_PATH = std::filesystem::current_path() / "results";
#else
    const fs::path CONFIG_PATH = fs::path(SOURCE_DIR) / "config.json";
    const fs::path RESULTS_DIRECTORY_PATH = fs::path(SOURCE_DIR) / "results";
#endif

int main()
{
    try
    {
        CFG = get_config_data(CONFIG_PATH);
        
        std::vector<test_combination> combinations{};
        if (CFG.USE_SPECIFIED_COMBINATIONS)
        {
            combinations = prepare_combinations(CFG.SPECIFIED_COMBINATIONS);
        }
        else
        {
            std::vector<std::vector<size_t>> schedules = cartesian_product(CFG.SCHEDULE_ELEMENTS);
            combinations = prepare_combinations(schedules, CFG.QBER);
        }

        std::vector<test_result> result = run_simulation(combinations);

        fmt::print(fg(fmt::color::green),"All tests were completed successfully! \nThe results will be written to the directory: {}\n", RESULTS_DIRECTORY_PATH.string());
        write_file(result, RESULTS_DIRECTORY_PATH);
    }
    catch(const std::exception& e)
    {
        fmt::print(stderr, fg(fmt::color::red), "ERROR: {}\n", e.what());
        return EXIT_FAILURE;
    }

    fmt::print(fg(fmt::color::green), "Press Enter to exit...");
    std::cin.get();

    return EXIT_SUCCESS;
}
